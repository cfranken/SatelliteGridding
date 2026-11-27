"""
    supported_backend_names() -> Tuple

Return the backend names accepted by [`resolve_backend`](@ref) and the CLI.
`sequential` uses the single-threaded Welford path; the others use the
KernelAbstractions accumulation path.
"""
supported_backend_names() = ("sequential", "cpu", "cuda", "metal")

"""
    backend_help_lines() -> Vector{String}

Short descriptions for CLI help and documentation.
"""
backend_help_lines() = [
    "sequential  Sequential Welford path (default, supports --compSTD)",
    "cpu         KernelAbstractions CPU backend",
    "cuda        KernelAbstractions CUDA backend (requires CUDA.jl)",
    "metal       KernelAbstractions Metal backend (requires Metal.jl on macOS/Apple GPU)",
]

"""
    resolve_backend(name; grid_spec=nothing, geometry_precision=nothing, accumulator_precision=nothing) -> Union{Nothing, KernelAbstractions.Backend}

Resolve a user-facing backend name to a KernelAbstractions backend object.
The optional GPU packages are weak dependencies and are loaded only when the
corresponding backend is requested. Pass `grid_spec` to check that it supports
the backend before loading any optional packages. Pass explicit KA
precisions here to validate them before loading a backend.
"""
function resolve_backend(name::AbstractString; grid_spec=nothing,
                         geometry_precision=nothing, accumulator_precision=nothing)
    backend = lowercase(strip(String(name)))
    backend = replace(backend, "_" => "-", "ka-" => "")

    if backend in ("", "none", "sequential")
        _require_spherical_ka_geometry(geometry_precision)
        _require_ka_accumulator(accumulator_precision)
        return nothing
    end
    backend in supported_backend_names() ||
        error("Unknown backend '$name'. Use one of: $(join(supported_backend_names(), ", ")).")
    locator = _validate_ka_grid(grid_spec)
    if grid_spec !== nothing && locator === nothing
        _require_spherical_ka_geometry(geometry_precision)
        grid_spec isa RectangularGridSpec || _require_ka_accumulator(accumulator_precision)
    end
    _validate_precision_type(geometry_precision, :geometry_precision)
    _validate_precision_type(accumulator_precision, :accumulator_precision)
    backend == "metal" && _validate_metal_precision(geometry_precision, accumulator_precision)
    if backend == "cpu"
        return CPU()
    elseif backend == "cuda"
        return _load_optional_backend(:CUDA, :SatelliteGriddingCUDAExt,
                                      :cuda_backend, "CUDA.jl", "CUDA")
    else
        return _load_optional_backend(:Metal, :SatelliteGriddingMetalExt,
                                      :metal_backend, "Metal.jl", "Metal")
    end
end

_validate_ka_grid(spec) = nothing

function _validate_ka_grid(spec::AbstractSphericalGridSpec)
    locator = _device_locator(spec)
    locator === nothing &&
        error("$(typeof(spec)) supports only the sequential CPU path. Omit `backend` or use --backend sequential.")
    _spherical_subdivision(spec) === Val(:unit_vectors) ||
        throw(ArgumentError("The device locator for $(typeof(spec)) requires _spherical_subdivision(spec) = Val(:unit_vectors)"))
    isbits(locator) || throw(ArgumentError("The device locator for $(typeof(spec)) must be isbits"))
    return locator
end

_default_ka_precision(backend) = Float64
_validate_ka_precision(backend, F, A) = nothing

function _validate_precision_type(precision, name)
    (precision === nothing || precision === Float32 || precision === Float64) ||
        throw(ArgumentError("$name must be Float32 or Float64"))
end

function _require_spherical_ka_geometry(geometry_precision)
    geometry_precision === nothing ||
        throw(ArgumentError("geometry_precision requires a spherical grid and a KA backend"))
end

function _require_ka_accumulator(accumulator_precision)
    accumulator_precision === nothing ||
        throw(ArgumentError("accumulator_precision requires a KA backend on a supported grid"))
end

function _validate_metal_precision(F, A)
    F === Float64 && throw(ArgumentError("Metal requires Float32 geometry_precision"))
    A === Float64 && throw(ArgumentError("Metal requires Float32 accumulator_precision"))
    return nothing
end

"""Validated locator and geometry and accumulator precisions for spherical KA."""
struct _SphericalKAOptions{L,F,A}
    locator::L
    geometry_precision::Type{F}
    accumulator_precision::Type{A}
end

function _spherical_ka_options(spec, backend, geometry_precision, accumulator_precision)
    locator = backend === nothing ? nothing : _validate_ka_grid(spec)
    if locator === nothing
        _require_spherical_ka_geometry(geometry_precision)
        _require_ka_accumulator(accumulator_precision)
        return nothing
    end
    default = _default_ka_precision(backend)
    _validate_precision_type(geometry_precision, :geometry_precision)
    _validate_precision_type(accumulator_precision, :accumulator_precision)
    F = something(geometry_precision, default)
    A = something(accumulator_precision, default)
    _validate_ka_precision(backend, F, A)
    return _SphericalKAOptions(locator, F, A)
end

_ka_options(spec, backend, F, A) = _spherical_ka_options(spec, backend, F, A)

function _ka_options(spec::RectangularGridSpec, backend, geometry_precision, accumulator_precision)
    _require_spherical_ka_geometry(geometry_precision)
    if backend === nothing
        _require_ka_accumulator(accumulator_precision)
        return nothing
    end
    _validate_precision_type(accumulator_precision, :accumulator_precision)
    A = something(accumulator_precision, _default_ka_precision(backend))
    _validate_ka_precision(backend, nothing, A)
    A
end

function _load_optional_backend(pkg::Symbol, extname::Symbol, fn::Symbol,
                                package_label::String, install_name::String)
    try
        Base.eval(Main, :(using $(pkg)))
    catch e
        error("$package_label backend requested but $package_label could not be loaded. Install it with `import Pkg; Pkg.add(\"$install_name\")` and make sure the backend is supported on this platform. Original error: $(sprint(showerror, e))")
    end

    ext = Base.get_extension(@__MODULE__, extname)
    ext === nothing &&
        error("$package_label loaded, but the SatelliteGridding extension $extname did not load.")
    Base.invokelatest(getproperty(ext, fn))
end

"""
    cuda_backend()

Return a CUDA KernelAbstractions backend after `CUDA.jl` has been loaded.
"""
function cuda_backend()
    ext = Base.get_extension(@__MODULE__, :SatelliteGriddingCUDAExt)
    ext === nothing &&
        error("CUDA backend requested but CUDA.jl is not loaded. Run `using CUDA` or use `resolve_backend(\"cuda\")`.")
    ext.cuda_backend()
end

"""
    metal_backend()

Return a Metal KernelAbstractions backend after `Metal.jl` has been loaded.
"""
function metal_backend()
    ext = Base.get_extension(@__MODULE__, :SatelliteGriddingMetalExt)
    ext === nothing &&
        error("Metal backend requested but Metal.jl is not loaded. Run `using Metal` or use `resolve_backend(\"metal\")`.")
    ext.metal_backend()
end

"""
    _ka_float64_type(backend) -> Type{<:AbstractFloat}

Return `Float64` when `backend` supports double precision (KA CPU, CUDA) and
`Float32` otherwise (Metal). Used only where the sequential path performs
`Float64` arithmetic internally (circular footprint geometry and intermediate
quadrilateral baseline interpolation).
"""
_ka_float64_type(backend) =
    KernelAbstractions.supports_float64(backend) ? Float64 : Float32

"""
    _ka_supported_float(backend, T) -> Type{<:AbstractFloat}

Return `T` unchanged when `backend` supports it, falling back to `Float32` when
`T === Float64` on a backend without double-precision support (Metal).
"""
_ka_supported_float(backend, ::Type{T}) where {T<:AbstractFloat} =
    (T === Float64 && !KernelAbstractions.supports_float64(backend)) ? Float32 : T
_ka_supported_float(backend, ::Type{<:Real}) = _ka_float64_type(backend)

_backend_input(::CPU, host) = host

function _backend_input(backend, host)
    device = KernelAbstractions.allocate(backend, eltype(host), size(host))
    copyto!(device, host)
    device
end

# Rectangular kernels preserve input precision where the backend supports it.
_rect_backend_input(::CPU, host) = host

function _rect_backend_input(backend, host)
    T = _ka_supported_float(backend, eltype(host))
    _backend_input(backend, eltype(host) === T ? host : T.(host))
end
