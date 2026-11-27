const _DEFAULT_MEAN_DIGITS = 8
const _DEFAULT_STD_DIGITS = 6
const _DEFAULT_CENTER_DIGITS = 5
const _OVERLAP_STATISTICS = ("sum_w2", "n_soundings")

"""
    NetCDFOptions(; compression="deflate", deflate_level=4, shuffle=false,
                    fill_value=-999f0, chunksizes=nothing, checksum="none",
                    precision="float32", significant_digits_mean=8,
                    significant_digits_std=6, significant_digits_center=5,
                    statistics=[])

NetCDF storage and rounding settings. `precision` selects `"float32"` or
`"float64"` payloads. `compression` is `"deflate"` or `"none"`; levels are
0–9. `fill_value` may be a number finite in the selected precision or NaN.
Each significant-digit setting accepts a positive integer or `"none"` to
disable rounding. Requests above 9 digits for Float32 or 17 for Float64
leave the output value unrounded. To preserve Float64 results exactly,
use Float64 accumulators with `precision="float64"` and set the applicable
significant-digit settings to `"none"`. See [`grid_l2`](@ref) and
[`grid_center`](@ref) for accumulator precision.
Chunk sizes follow Julia dimension order, including time (last). By default,
chunks contain at most 1024 × 1024 spatial cells (2^20 cells on a flat
cell axis), one cube face, and one time.
`checksum` is `"none"` or `"fletcher32"`.
`statistics` selects optional L2 NetCDF fields `"sum_w2"` and `"n_soundings"`.
"""
struct NetCDFOptions{T<:Union{Float32,Float64}}
    compression::String
    deflate_level::Int
    shuffle::Bool
    fill_value::T
    chunksizes::Union{Nothing,Tuple{Vararg{Int}}}
    checksum::Symbol
    significant_digits_mean::Union{Nothing,Int}
    significant_digits_std::Union{Nothing,Int}
    significant_digits_center::Union{Nothing,Int}
    statistics::Tuple{Vararg{Symbol}}

    function NetCDFOptions(; compression="deflate", deflate_level=4, shuffle=false,
                            fill_value=-999f0, chunksizes=nothing, checksum="none",
                            precision="float32", significant_digits_mean=_DEFAULT_MEAN_DIGITS,
                            significant_digits_std=_DEFAULT_STD_DIGITS,
                            significant_digits_center=_DEFAULT_CENTER_DIGITS, statistics=[])
        compression in ("deflate", "none") ||
            throw(ArgumentError("output.compression must be deflate or none"))
        deflate_level isa Integer && !(deflate_level isa Bool) && 0 <= deflate_level <= 9 ||
            throw(ArgumentError("output.deflate_level must be an integer from 0 to 9"))
        shuffle isa Bool || throw(ArgumentError("output.shuffle must be a boolean"))
        fill_value isa Real && !(fill_value isa Bool) ||
            throw(ArgumentError("output.fill_value must be a number or nan"))
        precision isa AbstractString && precision in ("float32", "float64") ||
            throw(ArgumentError("output.precision must be float32 or float64"))
        T = precision == "float32" ? Float32 : Float64
        fv = T(fill_value)
        isinf(fv) && throw(ArgumentError("output.fill_value must be finite in $precision or nan"))
        checksum in ("none", "fletcher32") ||
            throw(ArgumentError("output.checksum must be none or fletcher32"))
        if chunksizes !== nothing
            (chunksizes isa AbstractVector || chunksizes isa Tuple) &&
                !isempty(chunksizes) && all(x -> x isa Integer && !(x isa Bool) && x > 0, chunksizes) ||
                throw(ArgumentError("output.chunksizes must contain positive integers"))
            chunksizes = Tuple(Int.(chunksizes))
        end
        new{T}(compression, Int(deflate_level), shuffle, fv, chunksizes,
            checksum == "none" ? :nochecksum : :fletcher32,
            _significant_digits(significant_digits_mean, "significant_digits_mean"),
            _significant_digits(significant_digits_std, "significant_digits_std"),
            _significant_digits(significant_digits_center, "significant_digits_center"),
            _overlap_statistics(statistics))
    end
end

function _overlap_statistics(names)
    (names isa AbstractVector || names isa Tuple) &&
        all(name -> name isa AbstractString && name in _OVERLAP_STATISTICS, names) ||
        throw(ArgumentError("output.statistics must be an array of names chosen from sum_w2 and n_soundings"))
    Tuple(Symbol(name) for name in _OVERLAP_STATISTICS if name in names)
end

function _significant_digits(value, name)
    value isa AbstractString && value == "none" && return nothing
    value isa Integer && !(value isa Bool) && 0 < value <= typemax(Int) ||
        throw(ArgumentError("output.$name must be a positive integer or \"none\""))
    Int(value)
end

function _output_options(options::AbstractDict)
    settings = get(options, "output", Dict())
    settings isa AbstractDict || throw(ArgumentError("output must be a table"))
    known = ("compression", "deflate_level", "shuffle", "fill_value", "chunksizes", "checksum",
             "precision", "significant_digits_mean", "significant_digits_std",
             "significant_digits_center", "statistics")
    for key in keys(settings)
        key in known || throw(ArgumentError("Unknown [output] option: $key"))
    end
    NetCDFOptions(; (Symbol(k) => v for (k, v) in settings)...)
end

# At most 2^20 spatial values per time chunk, for flat and 2D spatial layouts.
const _NETCDF_CHUNK_EDGE = 1024
# Accumulated overlap weights are nonnegative. Gridders write zero for empty cells.
const _NETCDF_WEIGHT_FILL = -1f0

_default_netcdf_chunks(spatial::Tuple{Int}) = (min(spatial[1], _NETCDF_CHUNK_EDGE^2), 1)
_default_netcdf_chunks(spatial::Tuple) =
    (min(spatial[1], _NETCDF_CHUNK_EDGE), min(spatial[2], _NETCDF_CHUNK_EDGE),
     ntuple(_ -> 1, length(spatial) - 2)..., 1)

function _storage_kwargs(options::NetCDFOptions{T}, spec::AbstractGridSpec) where {T}
    spatial = _output_spatial_shape(spec)
    chunks = options.chunksizes === nothing ? _default_netcdf_chunks(spatial) : options.chunksizes
    length(chunks) == length(spatial) + 1 ||
        throw(ArgumentError("output.chunksizes needs $(length(spatial) + 1) entries, including time"))
    all(chunks[i] <= spatial[i] for i in eachindex(spatial)) ||
        throw(ArgumentError("output.chunksizes exceeds a fixed spatial dimension"))
    prod(big.(chunks)) * sizeof(T) < big(2)^32 ||
        throw(ArgumentError("NetCDF chunks must be smaller than 4 GiB"))
    level = options.compression == "none" ? 0 : options.deflate_level
    # NCDatasets 0.12 passes chunk sizes directly to libnetcdf as an array.
    storage = (; fillvalue=options.fill_value, chunksizes=collect(chunks), checksum=options.checksum)
    # On 0.12, supplying a level enables deflate even when the level is zero.
    level == 0 && return storage
    return (; storage..., deflatelevel=level, shuffle=options.shuffle)
end
