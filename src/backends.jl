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
    resolve_backend(name) -> Union{Nothing, KernelAbstractions.Backend}

Resolve a user-facing backend name to a KernelAbstractions backend object.
The optional GPU packages are weak dependencies and are loaded only when the
corresponding backend is requested.
"""
function resolve_backend(name::AbstractString)
    backend = lowercase(strip(String(name)))
    backend = replace(backend, "_" => "-", "ka-" => "")

    if backend in ("", "none", "sequential")
        return nothing
    elseif backend == "cpu"
        return CPU()
    elseif backend == "cuda"
        return _load_optional_backend(:CUDA, :SatelliteGriddingCUDAExt,
                                      :cuda_backend, "CUDA.jl", "CUDA")
    elseif backend == "metal"
        return _load_optional_backend(:Metal, :SatelliteGriddingMetalExt,
                                      :metal_backend, "Metal.jl", "Metal")
    end

    error("Unknown backend '$name'. Use one of: $(join(supported_backend_names(), ", ")).")
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
    getproperty(ext, fn)()
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

function _to_backend_array(backend, a::AbstractArray)
    backend isa CPU && return a
    E = _ka_supported_float(backend, eltype(a))
    host = eltype(a) === E ? a : E.(a)
    dev = KernelAbstractions.allocate(backend, E, size(host))
    copyto!(dev, host)
    return dev
end
