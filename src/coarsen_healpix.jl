const _COARSEN_DIMS = ("cell", "time")

"""
    coarsen_healpix(input, output; levels=nothing, nside=nothing, legacy_std=:error)

Write a coarser NESTED HEALPix NetCDF from a SatelliteGridding archive.
By default, coarsen one level. Specify the number of levels or the target
`nside`. Combine weighted means and population standard deviations in
Float64 and write them unrounded. For equivalence with direct gridding,
use a Float64 accumulator and unrounded Float64 output with a NaN fill
value. Use the same geometry precision and fixed `n_oversample` in both runs.

Tagged archives use their explicit statistic roles and links. For untagged
archives with `<name>` and `<name>_std` pairs, the default `legacy_std=:error`
requires an explicit choice. Use `:by_name` for complete, disjoint population
standard deviation pairs or `:none` to treat all science fields as independent
means. Untagged archives without pairs are coarsened as means.

Sum optional `sum_w2` and `n_soundings` across children. These sums
underestimate squared weights and overcount soundings when a footprint
touches multiple descendants of the same parent.
Do not use the summed statistics to estimate standard errors.
"""
function coarsen_healpix(input::AbstractString, output::AbstractString;
                         levels=nothing, nside=nothing, legacy_std=:error, command=nothing)
    (abspath(input) == abspath(output) || Base.samefile(input, output)) &&
        throw(ArgumentError("HEALPix coarsening requires distinct input and output files"))
    Dataset(input) do src
        level = _coarsen_source_level(src)
        levels = _coarsen_levels(level, levels, nside)
        previous_levels = get(src.attrib, "coarsen_levels", 0)
        previous_levels isa Integer && 0 <= previous_levels <= trailing_zeros(_HEALPIX_MAX_NSIDE) - level ||
            throw(ArgumentError("Invalid coarsen_levels in HEALPix archive"))
        spec = HEALPixGridSpec(nside=1 << (level - levels), T=Float64,
                               cell_bounds=haskey(src, "lon_bounds"))
        weights, fields, stds, statistics = _coarsen_fields(src, legacy_std)
        warning = _coarsen_precision_warning(src, weights, !isempty(stds))
        if !isempty(warning)
            @warn "HEALPix coarsening exactness is limited by the source archive: $warning. See the exact-archive recipe in docs/src/healpix.md."
        end
        Dataset(output, "c") do dst
            for (key, value) in pairs(src.attrib)
                key == "coarsen_command" || (dst.attrib[key] = value)
            end
            write_coordinates!(dst, spec)
            defDim(dst, "time", Inf)
            dst.attrib["coarsen_source"] = get(src.attrib, "coarsen_source", abspath(input))
            dst.attrib["coarsen_levels"] = Int32(previous_levels + levels)
            dst.attrib["coarsen_source_nside"] = 1 << (level + previous_levels)
            dst.attrib["coarsen_target_nside"] = spec.nside
            command === nothing || (dst.attrib["coarsen_command"] = command)
            isempty(warning) || (dst.attrib["coarsen_precision_warning"] = warning)
            dst.attrib["significant_digits_mean"] = "none"
            isempty(stds) || (dst.attrib["significant_digits_std"] = "none")

            for name in keys(src)
                var = src[name]
                dimnames(var) == ("time",) || continue
                dest = defVar(dst, name, eltype(var.var), ("time",); attrib=Dict(pairs(var.attrib)))
                isempty(var) || (dest.var[1:length(var)] = var.var[:])
            end
            options = NetCDFOptions(; precision="float64")
            storage = _storage_kwargs(options, spec)
            payload = Dict{String,Any}()
            mean_for_std = Dict(std => mean for (mean, std) in stds)
            for name in keys(src)
                var = src[name]
                dimnames(var) == _COARSEN_DIMS || continue
                attrs = Dict{String,Any}(pairs(var.attrib))
                pop!(attrs, "_FillValue", nothing)
                role = get(statistics, name, "")
                T = role == "n_soundings" ? Int32 : Float64
                # A valid parent mean can equal a finite child fill marker.
                fillvalue = name == weights || !isempty(role) ? T(-1) : T(NaN)
                if name == weights
                    attrs[_STATISTIC_ROLE] = "weight_sum"
                elseif haskey(mean_for_std, name)
                    attrs[_STATISTIC_ROLE] = "population_std"
                    attrs[_STATISTIC_MEAN] = mean_for_std[name]
                elseif !isempty(role)
                    note = role == "sum_w2" ?
                        "Summed over descendant cells. Underestimates the squared overlap when a footprint straddles descendants." :
                        "Summed over descendant cells. Overcounts footprints that straddle descendants."
                    comment = get(attrs, "comment", "")
                    occursin(note, comment) || (attrs["comment"] = string(comment, " ", note))
                    attrs["coarsening_approximation"] = "$note Do not use for sample variance or standard error."
                else
                    attrs[_STATISTIC_ROLE] = "mean"
                    haskey(stds, name) && (attrs[_STATISTIC_STD] = stds[name])
                end
                payload[name] = defVar(dst, name, T, _COARSEN_DIMS; storage..., fillvalue, attrib=attrs)
            end
            _coarsen_payload!(src, payload, weights, fields, stds, keys(statistics),
                              nside2npix(spec.nside), 1 << (2levels))
        end
    end
    nothing
end

function _coarsen_source_level(src)
    haskey(src, "healpix") && haskey(src, "cell_ids") &&
        haskey(src, "time") && haskey(src.dim, "cell") ||
        throw(ArgumentError("Expected a NESTED HEALPix NetCDF archive from SatelliteGridding.jl"))
    mapping = src["healpix"].attrib
    level = get(mapping, "refinement_level", -1)
    ids = src["cell_ids"]
    get(mapping, "grid_mapping_name", "") == "healpix" &&
        get(mapping, "indexing_scheme", "") == "nested" &&
        level isa Integer && 0 <= level <= trailing_zeros(_HEALPIX_MAX_NSIDE) &&
        src.dim["cell"] == nside2npix(1 << level) &&
        dimnames(ids) == ("cell",) ||
        throw(ArgumentError("Invalid NESTED HEALPix grid mapping or cell axis"))
    for cells in Iterators.partition(1:length(ids), _HEALPIX_COORD_CHUNK)
        ids.var[cells] == (first(cells) - 1):(last(cells) - 1) ||
            throw(ArgumentError("HEALPix coarsening requires dense, zero-based NESTED cell IDs in order"))
    end
    Int(level)
end

function _coarsen_levels(level, levels, nside)
    levels === nothing || nside === nothing ||
        throw(ArgumentError("Specify either levels or nside for HEALPix coarsening"))
    requested = nside === nothing ? something(levels, 1) : level - trailing_zeros(_healpix_nside(nside))
    requested isa Integer && !(requested isa Bool) && 1 <= requested <= level ||
        throw(ArgumentError("Target nside must be smaller than source nside $(1 << level) and at least 1"))
    Int(requested)
end

function _coarsen_fields(src, legacy_std)
    payload = String[]
    tagged = false
    for name in keys(src)
        var = src[name]
        has_links = any(attr -> haskey(var.attrib, attr), (_STATISTIC_ROLE, _STATISTIC_MEAN, _STATISTIC_STD))
        if dimnames(var) == _COARSEN_DIMS
            push!(payload, name)
            tagged |= has_links
        elseif has_links
            throw(ArgumentError("Invalid dimensions on tagged HEALPix variable '$name'"))
        end
    end
    if !tagged
        "n" in payload || throw(ArgumentError("HEALPix archive is missing its weight variable 'n' on (cell, time)"))
        return _coarsen_untagged_fields(src, filter(!=("n"), payload), legacy_std)
    end

    roles = Dict(name => get(src[name].attrib, _STATISTIC_ROLE, nothing) for name in payload)
    weights = [name for name in payload if roles[name] == "weight_sum"]
    length(weights) == 1 || throw(ArgumentError("HEALPix archive must have one weight variable"))
    fields = [name for name in payload if roles[name] == "mean"]
    stats = Dict(name => roles[name] for name in payload if roles[name] in _OVERLAP_STATISTICS)
    stds = Dict{String,String}()
    for name in payload
        role = roles[name]
        attrs = src[name].attrib
        if role == "mean"
            haskey(attrs, _STATISTIC_STD) || continue
            std = attrs[_STATISTIC_STD]
            get(roles, std, nothing) == "population_std" && get(src[std].attrib, _STATISTIC_MEAN, nothing) == name ||
                throw(ArgumentError("Invalid standard deviation reference on HEALPix mean '$name'"))
            stds[name] = std
        elseif role == "population_std"
            mean = get(attrs, _STATISTIC_MEAN, nothing)
            get(roles, mean, nothing) == "mean" && get(src[mean].attrib, _STATISTIC_STD, nothing) == name ||
                throw(ArgumentError("Invalid mean reference on HEALPix standard deviation '$name'"))
        elseif role != "weight_sum" && role ∉ _OVERLAP_STATISTICS
            throw(ArgumentError("Missing or unknown HEALPix statistic role '$role' on '$name'"))
        end
    end
    only(weights), fields, stds, stats
end

function _coarsen_untagged_fields(src, payload, legacy_std)
    legacy_std isa Symbol && legacy_std in (:error, :by_name, :none) ||
        throw(ArgumentError("legacy_std must be :error, :by_name, or :none"))
    stds = Dict{String,String}()
    if legacy_std !== :none
        for name in payload
            std = name * _STD_SUFFIX
            haskey(src, std) && (stds[name] = std)
        end
        if !isempty(stds)
            legacy_std === :error && throw(ArgumentError(
                "Ambiguous untagged HEALPix archive: a '<name>_std' field can be a population standard deviation or an independent mean. " *
                "Use legacy_std=:by_name (--legacy-std by_name) for standard deviations or legacy_std=:none (--legacy-std none) for independent means."))
            # The legacy writer generated a deviation for every mean when enabled.
            2 * length(stds) == length(payload) &&
                all(std -> dimnames(src[std]) == _COARSEN_DIMS && !haskey(stds, std), values(stds)) ||
                throw(ArgumentError("Ambiguous untagged HEALPix layout: expected complete, disjoint <mean>_std pairs on (cell, time)"))
        end
    end
    fields = filter(name -> name ∉ values(stds), payload)
    "n", fields, stds, Dict{String,String}()
end

function _coarsen_precision_warning(src, weights, has_std)
    reasons = String.(split(get(src.attrib, "coarsen_precision_warning", ""), ", "; keepempty=false))
    eltype(src[weights].var) === Float64 || push!(reasons, "Float32 output")
    precision = get(src.attrib, "accumulator_precision", "unknown")
    precision == "float32" && push!(reasons, "Float32 accumulation")
    precision == "unknown" && push!(reasons, "accumulator precision metadata is unavailable")
    names = has_std ? ("significant_digits_mean", "significant_digits_std") : ("significant_digits_mean",)
    digits = [get(src.attrib, name, "unknown") for name in names]
    "unknown" in digits && push!(reasons, "rounding settings are unavailable")
    any(d -> d ∉ ("none", "unknown"), digits) && push!(reasons, "decimal rounding")
    join(unique(reasons), ", ")
end

# Parent-aligned batches amortize I/O. Split a single parent's descendants
# into smaller reads when many levels are removed at once.
function _coarsen_payload!(src, dst, weight_name, fields, stds, stats, npix, descendants)
    weight_var = src[weight_name]
    # Older files can declare 0 as fill even though the writer stored 0 for empty cells.
    weight_data = get(weight_var.attrib, "_FillValue", nothing) == 0 ? weight_var.var : weight_var
    parent_batch = max(1, _HEALPIX_COORD_CHUNK ÷ descendants)
    for t in 1:src.dim["time"], parents in Iterators.partition(1:npix, parent_batch)
        n = zeros(Float64, length(parents))
        mean = zeros(Float64, length(parents), length(fields))
        reference = zeros(Float64, size(mean))
        m2 = zeros(Float64, size(mean))
        sums = [_coarsen_statistic_buffer(eltype(dst[name].var), length(parents)) for name in stats]
        first_child = (first(parents) - 1) * descendants + 1
        last_child = last(parents) * descendants
        for children in Iterators.partition(first_child:last_child, _HEALPIX_COORD_CHUNK)
            weights = weight_data[children, t]
            offset = first(children) - first_child
            fraction, previous = _coarsen_weights!(n, weights, offset, descendants)
            for (v, name) in enumerate(fields)
                values = src[name][children, t]
                deviations = haskey(stds, name) ? src[stds[name]][children, t] : nothing
                _coarsen_moments!(view(reference, :, v), view(mean, :, v), view(m2, :, v),
                                  values, deviations, weights, fraction, previous, offset, descendants)
            end
            for (sum, name) in zip(sums, stats)
                _coarsen_statistic!(sum, src[name][children, t], weights, offset, descendants)
            end
        end
        for (sum, name) in zip(sums, stats)
            _write_coarsened_statistic!(dst[name], sum, parents, t)
        end
        dst[weight_name].var[parents, t] = n
        for (v, name) in enumerate(fields)
            values = view(reference, :, v) .+ view(mean, :, v)
            _write_coarsened_field!(dst[name], values, n, parents, t)
            if haskey(stds, name)
                deviation = sqrt.(view(m2, :, v) ./ n)
                _write_coarsened_field!(dst[stds[name]], deviation, n, parents, t)
            end
        end
    end
    nothing
end

function _coarsen_weights!(n, weights, offset, descendants)
    fraction = zeros(Float64, length(weights))
    previous = similar(fraction)
    for j in eachindex(weights)
        weight = weights[j]
        ismissing(weight) && throw(ArgumentError("HEALPix archive contains unwritten weights"))
        w = Float64(weight)
        isfinite(w) && w >= 0 || throw(ArgumentError("Invalid HEALPix weight $w"))
        iszero(w) && continue
        p = cld(offset + j, descendants)
        previous[j] = n[p]
        n[p] += w
        fraction[j] = w / n[p]
    end
    fraction, previous
end

function _coarsen_moments!(reference, mean, m2, values, deviations, weights, fraction, previous, offset, descendants)
    for j in eachindex(weights)
        w = Float64(weights[j])
        w > 0 || continue
        p = cld(offset + j, descendants)
        child_mean = Float64(coalesce(values[j], NaN))
        if iszero(previous[j])
            reference[p] = child_mean
            delta = 0.0
        else
            # Keep the parent mean relative to its first observed child until writing.
            delta = (child_mean - reference[p]) - mean[p]
            mean[p] += fraction[j] * delta
        end
        if deviations !== nothing
            std = Float64(coalesce(deviations[j], NaN))
            cross = previous[j] * fraction[j]
            m2[p] += isfinite(child_mean) ? (w * std) * std + (cross * delta) * delta : NaN
        end
    end
    nothing
end

function _coarsen_statistic_buffer(::Type{T}, n) where {T}
    fill!(Vector{Union{Missing,T}}(undef, n), zero(T))
end

_coarsen_sum(::Type{Float64}, a, b) = a + b
_coarsen_sum(::Type{Int32}, a, b) = ismissing(a) || ismissing(b) ? missing : _add_sounding_counts(a, Int32(b))

function _coarsen_statistic!(totals::Vector{Union{Missing,T}}, values, weights, offset, descendants) where {T}
    for j in eachindex(weights)
        weights[j] > 0 || continue
        p = cld(offset + j, descendants)
        totals[p] = _coarsen_sum(T, totals[p], values[j])
    end
    nothing
end

function _write_coarsened_statistic!(out, values, cells, t)
    fill = out.attrib["_FillValue"]
    eltype(out.var) === Int32 && _check_sounding_counts(skipmissing(values))
    out.var[cells, t] = coalesce.(values, fill)
    nothing
end

function _write_coarsened_field!(out, values, weights, cells, t)
    fill = out.attrib["_FillValue"]
    out.var[cells, t] = [w > 0 && isfinite(x) ? x : fill for (x, w) in zip(values, weights)]
    nothing
end
