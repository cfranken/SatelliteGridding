# Backend-specific accumulation of prepared L2 observations.

"""Sequential weighted running means, with an optional Welford M2 array."""
struct WelfordAccumulator{T,M}
    data::Array{T,3}
    weights::Matrix{T}
    m2::M
end

"""Rectangular KA weighted sums and weights, stored on the selected backend."""
struct RectKAAccumulator{B,D,W}
    backend::B
    data::D
    weights::W
end

"""
Spherical KA shifted sums and weights, with a shared per-variable reference.
When STD is enabled, M2 and residual store bounds during the first pass and
deviations during replay.
"""
struct SphericalKAAccumulator{B,O,D,W,M,R,H}
    backend::B
    options::O
    data::D
    weights::W
    m2::M
    residual::M
    reference::R
    host_batch::H
end

function new_accumulator(spec::AbstractGridSpec{T}, ::Nothing, ::Nothing, n_vars, compute_std) where {T}
    shape = accumulator_shape(spec)
    data = zeros(T, shape..., n_vars)
    WelfordAccumulator(data, zeros(T, shape), compute_std ? zeros(T, size(data)) : nothing)
end

function new_accumulator(spec::AbstractGridSpec{T}, backend, ::Nothing, n_vars, compute_std) where {T}
    shape = accumulator_shape(spec)
    data = KernelAbstractions.zeros(backend, T, shape..., n_vars)
    weights = KernelAbstractions.zeros(backend, T, shape...)
    RectKAAccumulator(backend, data, weights)
end

function new_accumulator(spec, backend, options::_SphericalKAOptions, n_vars, compute_std)
    A = options.accumulator_precision
    shape = accumulator_shape(spec)
    data = KernelAbstractions.zeros(backend, A, shape..., n_vars)
    weights = KernelAbstractions.zeros(backend, A, shape...)
    m2 = compute_std ? KernelAbstractions.allocate(backend, A, size(data)) : nothing
    residual = compute_std ? KernelAbstractions.allocate(backend, A, size(data)) : nothing
    SphericalKAAccumulator(backend, options, data, weights, m2, residual,
                           _SphericalKAReference(backend, A, n_vars),
                           _SphericalKAHostBatch(options.geometry_precision, A))
end

_zero!(a) = fill!(a, zero(eltype(a)))
_zero!(::Nothing) = nothing

function reset!(acc::Union{RectKAAccumulator,SphericalKAAccumulator})
    _zero!(acc.data)
    _zero!(acc.weights)
    nothing
end

function reset!(acc::WelfordAccumulator)
    _zero!(acc.data)
    _zero!(acc.m2)
    _zero!(acc.weights)
    nothing
end

# A tiled snapshot chooses its own reference. Rolling days share one reference.
_reset_reference!(acc) = nothing
_reset_reference!(acc::SphericalKAAccumulator) = (acc.reference.initialized = false)

_start_pass!(acc) = nothing
_start_pass!(acc::SphericalKAAccumulator) =
    acc.m2 === nothing ? nothing : _spherical_mean_bounds!(acc.m2, acc.residual)

_l2_value_type(acc, ::AbstractGridSpec{T}) where {T} = T
# Screen native Float32/Float64 science values before narrowing for a
# Float32 accumulator. Finite values outside its range must reach the
# kernel so their cells get fill for the non-finite moments.
_l2_value_type(::SphericalKAAccumulator, ::AbstractGridSpec) = Float64

function accumulate_file!(acc, ctx, filepath, progress, pass)
    (; config, spec, run) = ctx
    replay = pass isa Tuple
    _with_l2_observations(filepath, config, spec, run.method, ctx.nc_vars,
                          !replay && ctx.fill_attrib[], progress;
                          metadata=replay ? nothing : ctx.metadata,
                          value_type=_l2_value_type(acc, spec)) do observations
        _accumulate!(acc, spec, observations, run.method, pass)
    end
end

_accumulate!(acc, spec, obs, method, ::Nothing) = _accumulate!(acc, spec, obs, method)

function _accumulate!(acc::Union{WelfordAccumulator,RectKAAccumulator},
                      spec::RectangularGridSpec, obs, method)
    ilat, ilon, centers = _rectangular_indices(obs, spec, method)
    n = _l2_oversampling(ilat, ilon, method.n_oversample)
    _accumulate_rectangular!(acc, ilat, ilon, centers, obs.values, n)
end

function _accumulate_rectangular!(acc::WelfordAccumulator{T}, ilat, ilon,
                                  ::Nothing, values, n) where {T}
    n_pixels, n_vars = size(values)
    points = zeros(T, n, n, 2)
    ix, iy = zeros(Int32, n^2), zeros(Int32, n^2)
    lats0, lons0, lats1, lons1 = zeros(n), zeros(n), zeros(n), zeros(n)
    accumulate_footprint!(acc.data, acc.m2, acc.weights, acc.m2 !== nothing,
                          ilat, ilon, values, n_pixels, n_vars, n,
                          points, ix, iy, lats0, lons0, lats1, lons1)
end

function _accumulate_rectangular!(acc::WelfordAccumulator, ilat, ilon,
                                  centers::Tuple, values, n)
    n_pixels, n_vars = size(values)
    accumulate_circular_footprint!(acc.data, acc.m2, acc.weights, acc.m2 !== nothing,
                                   centers..., ilat, ilon, values, n_pixels, n_vars, n)
end

function _accumulate_rectangular!(acc::RectKAAccumulator, ilat, ilon, centers, values, n)
    # L2 corners are already CCW-sorted in degrees. Sorting the fractional
    # indices again can rotate the starting corner due to Float32 rounding.
    _accumulate_rectangular_batch!(acc.backend, acc.data, acc.weights, centers,
                                   ilat, ilon, values, n, size(values, 2), false)
end

# Spherical grids locate each subpoint independently. Their subdivision
# trait selects the geometry and sampling buffers.
function _accumulate!(acc::WelfordAccumulator, spec::AbstractSphericalGridSpec, obs, method::SubpixelGridding)
    _accumulate_l2_spherical!(_spherical_subdivision(spec), acc.data, acc.m2,
                              acc.weights, acc.m2 !== nothing, spec, obs, method.n_oversample)
end

function _accumulate_l2_spherical!(::Val{:lonlat}, data, m2, weights, compute_std,
                                   spec, observations, n_override)
    n = n_override === nothing ? 4 : n_override
    points = zeros(_spherical_subpixel_type(spec), n, n, 2)
    lats0, lons0, lats1, lons1 = zeros(n), zeros(n), zeros(n), zeros(n)
    values = observations.values
    n_pixels, n_vars = size(values)
    _accumulate_spherical_footprint!(data, m2, weights, compute_std, spec,
                                     observations.lat_bounds, observations.lon_bounds,
                                     values, n_pixels, n_vars, n,
                                     points, lats0, lons0, lats1, lons1)
end

function _accumulate_l2_spherical!(::Val{:unit_vectors}, data, m2, weights, compute_std,
                                   spec, observations, n_override)
    capacity = n_override === nothing ? _HEALPIX_MAX_OVERSAMPLE : n_override
    cells = Vector{Tuple{Int,Int}}(undef, capacity^2)
    _accumulate_unit_vector_footprints!(data, m2, weights, compute_std, spec,
                                       observations.lat_bounds, observations.lon_bounds,
                                       observations.values, n_override, cells)
end

function _accumulate!(acc::SphericalKAAccumulator, spec::AbstractSphericalGridSpec, obs, method::SubpixelGridding)
    _accumulate_spherical_pass!(acc, obs, method.n_oversample, acc.data, nothing)
end

function _accumulate!(acc::SphericalKAAccumulator, spec::AbstractSphericalGridSpec,
                     obs, method::SubpixelGridding, bounds::NamedTuple)
    _accumulate_spherical_pass!(acc, obs, method.n_oversample, acc.data, bounds)
end

function _accumulate!(acc::SphericalKAAccumulator, spec::AbstractSphericalGridSpec,
                     obs, method::SubpixelGridding, replay::Tuple)
    _accumulate_spherical_pass!(acc, obs, method.n_oversample, acc.m2, replay)
end

function _accumulate_spherical_pass!(acc, obs, n_override, data, pass)
    _accumulate_spherical_batch!(acc.backend, data, acc.weights, pass, acc.options.locator,
        acc.options.geometry_precision, obs.lat_bounds, obs.lon_bounds, obs.values,
        n_override, acc.reference, acc.host_batch)
end

finish_pass!(acc, ctx, files, ::Nothing, keep_shifted) = nothing

"""
Replay the first-pass files around each cell mean, then correct the mean
and central M2 using the accumulated linear and squared deviations. A replay
error invalidates the window because that file already contributed to the mean.
With `keep_shifted=true`, store the reference-shifted mean and its low-order
correction for a later rolling-window merge.
"""
function finish_pass!(acc::SphericalKAAccumulator, ctx, files, bounds::NamedTuple, keep_shifted)
    _finish_spherical_pass!(acc, bounds, keep_shifted) do replay
        p_files = _progress_bar(length(files), "  STD files: ", ctx.run.progress)
        for filepath in files
            try
                accumulate_file!(acc, ctx, filepath, p_files, replay)
            catch e
                error("Error replaying file $filepath for standard deviation: $(sprint(showerror, e))")
            end
        end
    end
end

function _finish_spherical_pass!(replay!, acc::SphericalKAAccumulator, bounds, keep_shifted)
    (; backend, data, weights, m2, residual, reference) = acc
    _finalize_spherical_mean!(backend, data, weights, reference, bounds)
    _zero!(m2)
    _zero!(residual)
    replay!((data, residual))
    _correct_spherical_moments_ka!(backend, _SPHERICAL_WORKGROUP_SIZE)(
        data, m2, residual, weights, reference.device, Val(keep_shifted); ndrange=length(weights))
    KernelAbstractions.synchronize(backend)
    nothing
end

function finish!(acc::WelfordAccumulator)
    acc.m2 === nothing || finalize_std!(acc.m2, acc.weights)
    acc.data, acc.m2, acc.weights
end

function finish!(acc::RectKAAccumulator)
    data, weights = Array(acc.data), Array(acc.weights)
    finalize_mean!(data, weights)
    data, nothing, weights
end

function finish!(acc::SphericalKAAccumulator)
    if acc.m2 === nothing
        _finalize_spherical_mean!(acc.backend, acc.data, acc.weights, acc.reference)
    end
    data, weights = Array(acc.data), Array(acc.weights)
    std = if acc.m2 !== nothing
        m2 = Array(acc.m2)
        _finalize_spherical_std!(m2, weights)
        m2
    end
    data, std, weights
end
