# Backend-specific accumulation of prepared L2 observations.

"""Sequential weighted running means, with an optional Welford M2 array."""
struct WelfordAccumulator{T,M,S}
    data::Array{T,3}
    weights::Matrix{T}
    m2::M
    overlap::S
end

"""Buffers shared by the first-pass bounds and the refinement replay."""
struct KAReplay{S,D,B}
    squared::S
    residual::D
    bounds::B
end

_ka_m2(::Nothing) = nothing
_ka_m2(replay::KAReplay) = replay.squared

"""Shifted KA sums and weights with a shared reference and optional replay."""
struct KAMoments{B,D,W,P,R,S}
    backend::B
    data::D
    weights::W
    replay::P
    reference::R
    overlap::S
end

abstract type AbstractKAAccumulator end

struct RectKAAccumulator{M} <: AbstractKAAccumulator
    moments::M
end

struct SphericalKAAccumulator{M,O,H} <: AbstractKAAccumulator
    moments::M
    options::O
    host_batch::H
end

_l2_accumulator_precision(acc) = eltype(acc.weights)
_l2_accumulator_precision(acc::AbstractKAAccumulator) = _l2_accumulator_precision(acc.moments)
_l2_geometry_precision(::WelfordAccumulator) = Float64
_l2_geometry_precision(acc::SphericalKAAccumulator) = acc.options.geometry_precision

function new_accumulator(spec::AbstractGridSpec{T}, ::Nothing, ::Nothing, n_vars, compute_std,
                         step=nothing, statistics=()) where {T}
    # A Float64 merge cannot recover variation already lost in a day's mean
    # or in the input values. Keep tiled Welford accumulation in grid precision.
    A = step === nothing ? T : Float64
    shape = accumulator_shape(spec)
    data = zeros(A, shape..., n_vars)
    overlap = _new_overlap(T -> zeros(T, shape), A, statistics)
    WelfordAccumulator(data, zeros(A, shape), compute_std ? zeros(A, size(data)) : nothing, overlap)
end

function _ka_replay(backend, data::AbstractArray{A}, ::Val{Std}) where {A,Std}
    (Std || A === Float32) || return nothing
    squared = Std ? KernelAbstractions.allocate(backend, A, size(data)) : nothing
    residual = KernelAbstractions.allocate(backend, A, size(data))
    KAReplay(squared, residual, _ka_replay_bounds(backend, size(data), A, squared, residual))
end

# Float32 rounds enough during replay to require the bounds a second time.
function _ka_replay_bounds(backend, shape, ::Type{Float32}, squared, residual)
    lo = KernelAbstractions.allocate(backend, UInt32, shape)
    hi = KernelAbstractions.allocate(backend, UInt32, shape)
    return (; lo, hi)
end

# Float64 consumes its bounds before replay and reuses that storage for q and r.
_ka_replay_bounds(backend, shape, ::Type{Float64}, squared, residual) =
    (lo=reinterpret(UInt64, squared), hi=reinterpret(UInt64, residual))

_ka_refinement_bounds(replay, ::Type{Float32}) = replay.bounds
_ka_refinement_bounds(replay, ::Type{Float64}) = nothing

function KAMoments(backend, ::Type{A}, shape, n_vars, compute_std, statistics) where {A}
    data = KernelAbstractions.zeros(backend, A, shape..., n_vars)
    weights = KernelAbstractions.zeros(backend, A, shape...)
    replay = _ka_replay(backend, data, Val(compute_std))
    overlap = _new_overlap(T -> KernelAbstractions.zeros(backend, T, shape...), A, statistics)
    KAMoments(backend, data, weights, replay, _KAReference(backend, A, n_vars), overlap)
end

function new_accumulator(spec::RectangularGridSpec, backend, ::Type{A}, n_vars, compute_std,
                         step=nothing, statistics=()) where {A}
    RectKAAccumulator(KAMoments(backend, A, accumulator_shape(spec), n_vars, compute_std, statistics))
end

function new_accumulator(spec, backend, options::_SphericalKAOptions, n_vars, compute_std,
                         step=nothing, statistics=())
    moments = KAMoments(backend, options.accumulator_precision, accumulator_shape(spec),
                        n_vars, compute_std, statistics)
    SphericalKAAccumulator(moments, options, _SphericalKAHostBatch(options.geometry_precision))
end

_zero!(a) = fill!(a, zero(eltype(a)))
_zero!(::Nothing) = nothing

reset!(acc::AbstractKAAccumulator) = reset!(acc.moments)

function reset!(acc::KAMoments)
    _zero!(acc.data)
    _zero!(acc.weights)
    _map_overlap(_zero!, acc.overlap)
    nothing
end

function reset!(acc::WelfordAccumulator)
    _zero!(acc.data)
    _zero!(acc.m2)
    _zero!(acc.weights)
    _map_overlap(_zero!, acc.overlap)
    nothing
end

# A tiled snapshot chooses its own reference. Rolling days share one reference.
_reset_reference!(acc) = nothing
_reset_reference!(acc::AbstractKAAccumulator) = (acc.moments.reference.initialized = false)

_start_pass!(acc) = nothing
_start_pass!(acc::AbstractKAAccumulator) = _start_pass!(acc.moments)
_start_pass!(acc::KAMoments) = _start_pass!(acc.replay)
_start_pass!(replay::KAReplay) = _reset_ka_bounds!(replay.bounds)

_l2_value_type(acc::WelfordAccumulator{T}) where {T} = T
# Screen native science values before narrowing for a Float32 KA accumulator.
# Finite values outside its range must reach the kernel so affected cells get fill.
_l2_value_type(::AbstractKAAccumulator) = Float64

function accumulate_file!(acc, ctx, filepath, progress, pass)
    (; config, spec, run) = ctx
    replay = pass isa Tuple
    _with_l2_observations(filepath, config, spec, run.method, ctx.nc_vars,
                          !replay && ctx.fill_attrib[], progress;
                          metadata=replay ? nothing : ctx.metadata,
                          value_type=_l2_value_type(acc)) do observations
        _accumulate!(acc, spec, observations, run.method, pass)
    end
end

_accumulate!(acc, spec, obs, method, ::Nothing) = _accumulate!(acc, spec, obs, method)

function _accumulate!(acc::WelfordAccumulator, spec::RectangularGridSpec{T}, obs, method) where {T}
    ilat, ilon, centers = _rectangular_indices(obs, spec, method)
    n = _l2_oversampling(ilat, ilon, method.n_oversample)
    _accumulate_rectangular!(acc, T, ilat, ilon, centers, obs.values, n)
end

function _accumulate_rectangular!(acc::WelfordAccumulator, ::Type{T}, ilat, ilon,
                                  ::Nothing, values, n) where {T}
    n_pixels, n_vars = size(values)
    points = zeros(T, n, n, 2)
    ix, iy = zeros(Int32, n^2), zeros(Int32, n^2)
    lats0, lons0, lats1, lons1 = zeros(n), zeros(n), zeros(n), zeros(n)
    accumulate_footprint!(acc.data, acc.m2, acc.weights, acc.m2 !== nothing,
                          ilat, ilon, values, n_pixels, n_vars, n,
                          points, ix, iy, lats0, lons0, lats1, lons1, acc.overlap)
end

function _accumulate_rectangular!(acc::WelfordAccumulator, ::Type, ilat, ilon,
                                  centers::Tuple, values, n)
    n_pixels, n_vars = size(values)
    accumulate_circular_footprint!(acc.data, acc.m2, acc.weights, acc.m2 !== nothing,
                                   centers..., ilat, ilon, values, n_pixels, n_vars, n, acc.overlap)
end

function _accumulate!(acc::WelfordAccumulator, spec::AbstractSphericalGridSpec, obs, method::SubpixelGridding)
    _accumulate_l2_spherical!(acc.data, acc.m2, acc.weights, acc.m2 !== nothing,
                              spec, obs, method.n_oversample, acc.overlap)
end

function _accumulate_l2_spherical!(data, m2, weights, compute_std,
                                   spec, observations, n_override, overlap=nothing)
    capacity = n_override === nothing ? _SPHERICAL_MAX_OVERSAMPLE : n_override
    cells = Vector{Int}(undef, capacity^2)
    _accumulate_unit_vector_footprints!(data, m2, weights, compute_std, spec,
                                       observations.lat_bounds, observations.lon_bounds,
                                       observations.values, n_override, cells, overlap)
end

_accumulate!(acc::AbstractKAAccumulator, spec, obs, method) =
    _accumulate_ka_pass!(acc, spec, obs, method, acc.moments.data, nothing)

_accumulate!(acc::AbstractKAAccumulator, spec, obs, method, bounds::NamedTuple) =
    _accumulate_ka_pass!(acc, spec, obs, method, acc.moments.data, bounds)

_accumulate!(acc::AbstractKAAccumulator, spec, obs, method, replay::Tuple) =
    _accumulate_ka_pass!(acc, spec, obs, method, _ka_m2(acc.moments.replay), replay)

_ka_scatter_pass(pass, ::Nothing) = pass
_ka_scatter_pass(pass, overlap::NamedTuple) = (bounds=pass, overlap=overlap)
_ka_scatter_pass(pass::Tuple, ::NamedTuple) = pass

# Opted-in quads use one update per footprint and cell in both passes.
_rect_sample_counts(counts, nsq, overlap) = counts
_rect_sample_counts(::Nothing, nsq, ::NamedTuple) = nsq

function _accumulate_ka_pass!(acc::RectKAAccumulator, spec::RectangularGridSpec{T}, obs, method,
                              data, pass) where {T}
    (; backend, weights, reference, overlap) = acc.moments
    ilat, ilon, centers = _rectangular_indices(obs, spec, method)
    n = _l2_oversampling(ilat, ilon, method.n_oversample)
    # The points must round in grid precision to mirror sequential geometry.
    # Preserve the corner order established in degrees by input preparation.
    (; ix, iy, counts, skip) = _rect_batch_geometry(backend, size(weights), T, centers,
                                                    ilat, ilon, n, false)
    values = _ka_batch_values!(backend, reference, obs.values)
    _scatter_rectangular!(backend, data, weights, _ka_scatter_pass(pass, overlap),
                          ix, iy, _rect_sample_counts(counts, n^2, overlap), skip,
                          values, reference.device, n^2, size(values, 2))
end

function _accumulate_ka_pass!(acc::SphericalKAAccumulator, spec, obs, method, data, pass)
    _accumulate_spherical_pass!(acc, spec, obs, method.n_oversample, data, pass)
end

function _accumulate_spherical_pass!(acc, spec, obs, n_override, data, pass)
    (; backend, weights, reference, overlap) = acc.moments
    _accumulate_spherical_batch!(backend, data, weights, _ka_scatter_pass(pass, overlap),
        spec, acc.options.locator,
        acc.options.geometry_precision, obs.lat_bounds, obs.lon_bounds, obs.values,
        n_override, reference, acc.host_batch)
end

finish_pass!(acc, ctx, files, ::Nothing, keep_mean_parts) = nothing

"""
Replay the first-pass files around each cell mean, then correct the mean
and central M2 using the accumulated linear and squared deviations. A replay
error invalidates the window because that file already contributed to the mean.
With `keep_mean_parts=true`, store the mean and its low-order correction
for a later rolling-window merge.
"""
function finish_pass!(acc::AbstractKAAccumulator, ctx, files, bounds::NamedTuple, keep_mean_parts)
    _finish_ka_pass!(acc.moments, bounds, keep_mean_parts) do replay
        p_files = _progress_bar(length(files), "  Replay files: ", ctx.run.progress)
        for filepath in files
            try
                accumulate_file!(acc, ctx, filepath, p_files, replay)
            catch e
                error("Error replaying file $filepath for moment refinement: $(sprint(showerror, e))")
            end
        end
    end
end

function _finish_ka_pass!(replay!, acc::KAMoments, bounds, keep_mean_parts)
    (; backend, data, weights, replay, reference) = acc
    m2 = _ka_m2(replay)
    residual = replay.residual
    _finalize_ka_mean!(backend, data, weights, reference, bounds)
    _zero!(m2)
    _zero!(residual)
    replay!((data, residual))
    _correct_ka_moments_ka!(backend, _KA_WORKGROUP_SIZE)(
        data, m2, residual, weights, reference.device, _ka_refinement_bounds(replay, eltype(data)),
        Val(keep_mean_parts); ndrange=length(weights))
    KernelAbstractions.synchronize(backend)
    nothing
end

function finish!(acc::WelfordAccumulator)
    acc.m2 === nothing || finalize_std!(acc.m2, acc.weights)
    acc.data, acc.m2, acc.weights, acc.overlap
end

finish!(acc::AbstractKAAccumulator) = finish!(acc.moments)

function finish!(acc::KAMoments{B,D,W,Nothing}) where {B,D,W}
    _finalize_ka_mean!(acc.backend, acc.data, acc.weights, acc.reference)
    Array(acc.data), nothing, Array(acc.weights), _map_overlap(Array, acc.overlap)
end

function finish!(acc::KAMoments)
    data, weights = Array(acc.data), Array(acc.weights)
    data, _ka_std(_ka_m2(acc.replay), weights), weights, _map_overlap(Array, acc.overlap)
end

_ka_std(::Nothing, weights) = nothing
function _ka_std(m2, weights)
    std = Array(m2)
    _finalize_centered_std!(std, weights)
    std
end
