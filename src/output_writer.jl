_convert_output(::Type{T}, data::Array{T}) where {T} = data
_convert_output(::Type{T}, data::Array) where {T} = T.(data)
_convert_output(::Type, ::Nothing) = nothing

# Time-slice output shared by rectangular, cubed-sphere, HEALPix, and center gridding.
_round_output(value, ::Nothing) = value
function _round_output(value, sigdigits::Int)
    rounded = round(value; sigdigits)
    # Julia's decimal rounding can return signed zero on overflow. A nonzero
    # input rounded to a positive number of significant digits cannot round to zero.
    iszero(rounded) && !iszero(value) ? oftype(value, NaN) : rounded
end

_output_digits(::Type{T}, digits) where {T} =
    digits !== nothing && digits <= (T === Float32 ? 9 : 17) ? digits : nothing

function _masked_output_slice(data, weights, fill_value::T, sigdigits, default_sigdigits;
                              min_weight=1e-10) where {T}
    result = similar(data, T)
    digits = _output_digits(T, sigdigits)
    # Legacy Float32 output rounded in grid precision before narrowing.
    legacy_rounding = T === Float32 && eltype(data) !== T && digits === default_sigdigits
    @inbounds for i in eachindex(result, data, weights)
        value = T(data[i])
        if weights[i] < min_weight || !isfinite(value)
            result[i] = fill_value
            continue
        end
        rounded = legacy_rounding ? T(_round_output(data[i], digits)) : _round_output(value, digits)
        result[i] = isfinite(rounded) ? rounded : fill_value
    end
    result
end

function _write_output_slice!(nc_var, data, ct, shape)
    spatial = ntuple(_ -> Colon(), length(shape))
    # Unfold cubed-sphere accumulators or flatten HEALPix cells in NESTED
    # order only at the I/O boundary. Rectangular arrays have the same
    # shape in memory and on disk.
    nc_var[spatial..., ct] = reshape(data, shape)
    nothing
end

function _write_time_coordinates!(ds, ct, date; window_start=nothing, window_end=nothing)
    ds["time"][ct] = date
    window_start === nothing || (ds["window_start"][ct] = window_start)
    window_end === nothing || (ds["window_end"][ct] = window_end)
    nothing
end

function _write_field!(var, data, weights, ct, shape, options, digits, default_digits;
                       min_weight=1e-10)
    values = _masked_output_slice(data, weights, options.fill_value, digits, default_digits; min_weight)
    _write_output_slice!(var, values, ct, shape)
end

function _write_time_slice!(ds_out, nc_vars, config, spec::AbstractGridSpec, output_options,
                            grid_data, grid_std, grid_weights, compute_std, ct, date;
                            window_start=nothing, window_end=nothing, overlap=nothing)
    shape = _output_spatial_shape(spec)
    _check_overlap_slice(overlap)
    _write_output_slice!(nc_vars["n"], grid_weights, ct, shape)
    _write_overlap_slice!(nc_vars, overlap, ct, shape)
    _write_time_coordinates!(ds_out, ct, date; window_start=window_start, window_end=window_end)
    maximum(grid_weights) > 0 || return nothing
    for (column, key) in enumerate(keys(config.grid_vars))
        _write_field!(nc_vars[key], @view(grid_data[:, :, column]), grid_weights, ct, shape,
                      output_options, output_options.significant_digits_mean, _DEFAULT_MEAN_DIGITS)
        if compute_std
            _write_field!(nc_vars[key * "_std"], @view(grid_std[:, :, column]), grid_weights, ct, shape,
                          output_options, output_options.significant_digits_std, _DEFAULT_STD_DIGITS)
        end
    end
    nothing
end

_check_overlap_slice(::Nothing) = nothing
function _check_overlap_slice(overlap::NamedTuple)
    haskey(overlap, :n_soundings) && _check_sounding_counts(overlap.n_soundings)
    nothing
end

_write_overlap_slice!(nc_vars, ::Nothing, ct, shape) = nothing
function _write_overlap_slice!(nc_vars, overlap, ct, shape)
    for (name, data) in pairs(overlap)
        _write_output_slice!(nc_vars[string(name)], data, ct, shape)
    end
    nothing
end

_output_arrays(::NetCDFOptions{Float64}, T, data, std, weights) = (data, std, weights)

# Preserve the existing rounding inputs when the grid and KA accumulator differ.
_output_arrays(::NetCDFOptions{Float32}, T, data, std, weights) =
    (_convert_output(T, data), _convert_output(T, std), _convert_output(T, weights))

function _emit_l2_snapshot!(ds, nc_vars, zg, config, spec::AbstractGridSpec{T}, output_options,
                            data, std, weights,
                            compute_std, ct, date; window_start, window_end, zarr_epoch,
                            overlap=nothing) where {T}
    data, std, weights = _output_arrays(output_options, T, data, std, weights)
    ds === nothing || _write_time_slice!(ds, nc_vars, config, spec, output_options, data, std,
        weights, compute_std, ct, date; window_start, window_end, overlap)
    zg === nothing || write_zarr_day!(zg, date, data, std, weights, config.grid_vars,
                                      compute_std; epoch=zarr_epoch)
    nothing
end

function _write_center_time_slice!(ds, nc_vars, vars, veg_names, options, data, weights, ct, d, min_count)
    shape = size(weights)
    _write_time_coordinates!(ds, ct, d)
    _write_output_slice!(nc_vars["n"], weights, ct, shape)

    for (column, key) in enumerate(Iterators.flatten((keys(vars), veg_names)))
        # Center outputs always write a masked slice, including for empty snapshots.
        values = @view(data[:, :, column]) ./ max.(weights, eps(eltype(data)))
        _write_field!(nc_vars[key], values, weights, ct, shape, options,
                      options.significant_digits_center, _DEFAULT_CENTER_DIGITS; min_weight=min_count)
    end
    nothing
end
