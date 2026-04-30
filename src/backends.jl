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
