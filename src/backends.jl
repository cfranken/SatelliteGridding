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
the backend before loading any optional packages. Pass explicit spherical
precisions here to validate them before loading a backend.
"""
function resolve_backend(name::AbstractString; grid_spec=nothing,
                         geometry_precision=nothing, accumulator_precision=nothing)
    backend = lowercase(strip(String(name)))
    backend = replace(backend, "_" => "-", "ka-" => "")

    if backend in ("", "none", "sequential")
        _require_spherical_ka(geometry_precision, accumulator_precision)
        return nothing
    end
    backend in supported_backend_names() ||
        error("Unknown backend '$name'. Use one of: $(join(supported_backend_names(), ", ")).")
    locator = _validate_ka_grid(grid_spec)
    if grid_spec !== nothing && locator === nothing
        _require_spherical_ka(geometry_precision, accumulator_precision)
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

function _validate_ka_grid(spec)
    spec isa AbstractSphericalGridSpec || return nothing
    locator = _device_locator(spec)
    locator === nothing &&
        error("$(typeof(spec)) supports only the sequential CPU path. Omit `backend` or use --backend sequential.")
    _spherical_subdivision(spec) === Val(:unit_vectors) ||
        throw(ArgumentError("The device locator for $(typeof(spec)) requires _spherical_subdivision(spec) = Val(:unit_vectors)"))
    isbits(locator) || throw(ArgumentError("The device locator for $(typeof(spec)) must be isbits"))
    return locator
end

_default_spherical_ka_precision(backend) = Float64
_validate_spherical_ka_precision(backend, F, A) = nothing

function _validate_precision_type(precision, name)
    (precision === nothing || precision === Float32 || precision === Float64) ||
        throw(ArgumentError("$name must be Float32 or Float64"))
end

function _require_spherical_ka(geometry_precision, accumulator_precision)
    geometry_precision === nothing ||
        throw(ArgumentError("geometry_precision requires a spherical grid and a KA backend"))
    accumulator_precision === nothing ||
        throw(ArgumentError("accumulator_precision requires a spherical grid and a KA backend"))
end

function _validate_metal_precision(F, A)
    F === Float64 && throw(ArgumentError("Metal requires Float32 geometry_precision"))
    A === Float64 && throw(ArgumentError("Metal requires Float32 accumulator_precision"))
    return nothing
end

function _spherical_ka_options(spec, backend, geometry_precision, accumulator_precision)
    locator = backend === nothing ? nothing : _validate_ka_grid(spec)
    if locator === nothing
        _require_spherical_ka(geometry_precision, accumulator_precision)
        return nothing
    end
    default = _default_spherical_ka_precision(backend)
    _validate_precision_type(geometry_precision, :geometry_precision)
    _validate_precision_type(accumulator_precision, :accumulator_precision)
    F = something(geometry_precision, default)
    A = something(accumulator_precision, default)
    _validate_spherical_ka_precision(backend, F, A)
    return (; locator, geometry_precision=F, accumulator_precision=A)
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
