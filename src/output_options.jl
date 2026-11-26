"""
    NetCDFOptions(; compression="deflate", deflate_level=4, shuffle=false,
                    fill_value=-999f0, chunksizes=nothing, checksum="none")

Lossless Float32 NetCDF storage settings. `compression` is `"deflate"` or
`"none"`; levels are 0–9. `fill_value` may be a finite number or NaN.
Chunk sizes follow Julia dimension order, including time (last). By default,
chunks contain at most 1024 × 1024 spatial cells (2^20 cells on a flat
cell axis), one cube face, and one time.
`checksum` is `"none"` or `"fletcher32"`.
"""
struct NetCDFOptions
    compression::String
    deflate_level::Int
    shuffle::Bool
    fill_value::Float32
    chunksizes::Union{Nothing,Tuple{Vararg{Int}}}
    checksum::Symbol

    function NetCDFOptions(; compression="deflate", deflate_level=4, shuffle=false,
                            fill_value=-999f0, chunksizes=nothing, checksum="none")
        compression in ("deflate", "none") ||
            throw(ArgumentError("output.compression must be deflate or none"))
        deflate_level isa Integer && !(deflate_level isa Bool) && 0 <= deflate_level <= 9 ||
            throw(ArgumentError("output.deflate_level must be an integer from 0 to 9"))
        shuffle isa Bool || throw(ArgumentError("output.shuffle must be a boolean"))
        fill_value isa Real && !(fill_value isa Bool) ||
            throw(ArgumentError("output.fill_value must be a number or nan"))
        fv = Float32(fill_value)
        isinf(fv) && throw(ArgumentError("output.fill_value must be finite Float32 or nan"))
        checksum in ("none", "fletcher32") ||
            throw(ArgumentError("output.checksum must be none or fletcher32"))
        if chunksizes !== nothing
            (chunksizes isa AbstractVector || chunksizes isa Tuple) &&
                !isempty(chunksizes) && all(x -> x isa Integer && !(x isa Bool) && x > 0, chunksizes) ||
                throw(ArgumentError("output.chunksizes must contain positive integers"))
            chunksizes = Tuple(Int.(chunksizes))
        end
        new(compression, Int(deflate_level), shuffle, fv, chunksizes,
            checksum == "none" ? :nochecksum : :fletcher32)
    end
end

function _output_options(options::AbstractDict)
    settings = get(options, "output", Dict())
    settings isa AbstractDict || throw(ArgumentError("output must be a table"))
    known = ("compression", "deflate_level", "shuffle", "fill_value", "chunksizes", "checksum")
    for key in keys(settings)
        key in known || throw(ArgumentError("Unknown [output] option: $key"))
    end
    NetCDFOptions(; (Symbol(k) => v for (k, v) in settings)...)
end

# At most 2^20 spatial values per time chunk, for flat and 2D spatial layouts.
const _NETCDF_CHUNK_EDGE = 1024

_default_netcdf_chunks(spatial::Tuple{Int}) = (min(spatial[1], _NETCDF_CHUNK_EDGE^2), 1)
_default_netcdf_chunks(spatial::Tuple) =
    (min(spatial[1], _NETCDF_CHUNK_EDGE), min(spatial[2], _NETCDF_CHUNK_EDGE),
     ntuple(_ -> 1, length(spatial) - 2)..., 1)

function _storage_kwargs(options::NetCDFOptions, spec::AbstractGridSpec)
    spatial = _output_spatial_shape(spec)
    chunks = options.chunksizes === nothing ? _default_netcdf_chunks(spatial) : options.chunksizes
    length(chunks) == length(spatial) + 1 ||
        throw(ArgumentError("output.chunksizes needs $(length(spatial) + 1) entries, including time"))
    all(chunks[i] <= spatial[i] for i in eachindex(spatial)) ||
        throw(ArgumentError("output.chunksizes exceeds a fixed spatial dimension"))
    prod(big.(chunks)) * sizeof(Float32) < big(2)^32 ||
        throw(ArgumentError("NetCDF chunks must be smaller than 4 GiB"))
    level = options.compression == "none" ? 0 : options.deflate_level
    # NCDatasets 0.12 passes chunk sizes directly to libnetcdf as an array.
    storage = (; fillvalue=options.fill_value, chunksizes=collect(chunks), checksum=options.checksum)
    # On 0.12, supplying a level enables deflate even when the level is zero.
    level == 0 && return storage
    return (; storage..., deflatelevel=level, shuffle=options.shuffle)
end

_output_spatial_shape(spec::RectangularGridSpec) = grid_shape(spec)
_output_spatial_shape(spec::CubedSphereGridSpec) = (spec.Nc, spec.Nc, 6)
_output_spatial_shape(spec::HEALPixGridSpec) = (nside2npix(spec.nside),)

_output_fill(nc_var) = Float32(nc_var.attrib["_FillValue"])
