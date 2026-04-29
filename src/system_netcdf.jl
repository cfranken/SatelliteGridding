const _SYSTEM_NETCDF_LIB = Ref{Union{Nothing,String}}(nothing)
const _SYSTEM_NETCDF_HANDLE = Ref{Ptr{Cvoid}}(C_NULL)
const _SYSTEM_NETCDF_CHECKED = Ref(false)

const _NC_NOWRITE = Cint(0)
const _NC_BYTE = Cint(1)
const _NC_SHORT = Cint(3)
const _NC_INT = Cint(4)
const _NC_FLOAT = Cint(5)
const _NC_DOUBLE = Cint(6)
const _NC_UBYTE = Cint(7)
const _NC_USHORT = Cint(8)
const _NC_UINT = Cint(9)

function _system_netcdf_lib()
    if !_SYSTEM_NETCDF_CHECKED[]
        lib = Libdl.find_library(["libnetcdf"],
                                 ["/lib64", "/usr/lib64",
                                  "/usr/lib/x86_64-linux-gnu", "/usr/lib"])
        _SYSTEM_NETCDF_LIB[] = isempty(lib) ? nothing : lib
        _SYSTEM_NETCDF_CHECKED[] = true
    end
    _SYSTEM_NETCDF_LIB[]
end

function _system_netcdf_handle()
    lib = _system_netcdf_lib()
    lib === nothing && error("No system libnetcdf found for HDF4 fallback")
    if _SYSTEM_NETCDF_HANDLE[] == C_NULL
        _SYSTEM_NETCDF_HANDLE[] = Libdl.dlopen(lib)
    end
    _SYSTEM_NETCDF_HANDLE[]
end

_nc_sym(name::Symbol) = Libdl.dlsym(_system_netcdf_handle(), name)

function _nc_error(lib::String, code::Cint)
    f = _nc_sym(:nc_strerror)
    msg = unsafe_string(ccall(f, Cstring, (Cint,), code))
    ErrorException("NetCDF C error $code: $msg")
end

function _nc_check(lib::String, code::Cint)
    code == 0 || throw(_nc_error(lib, code))
    nothing
end

"""
    read_system_netcdf_variable(path, varname)

Read a flat variable through the system `libnetcdf`. This is primarily a fallback
for HDF4/HDF-EOS2 files when the Julia NetCDF_jll build lacks HDF4 support.
"""
function read_system_netcdf_variable(path::String, varname::String)
    occursin("/", varname) && error("System NetCDF fallback only supports flat variable names: $varname")

    lib = _system_netcdf_lib()
    lib === nothing && error("No system libnetcdf found for HDF4 fallback")

    ncid = Ref{Cint}()
    f_open = _nc_sym(:nc_open)
    _nc_check(lib, ccall(f_open, Cint,
                         (Cstring, Cint, Ref{Cint}), path, _NC_NOWRITE, ncid))
    try
        varid = Ref{Cint}()
        f_inq_varid = _nc_sym(:nc_inq_varid)
        _nc_check(lib, ccall(f_inq_varid, Cint,
                             (Cint, Cstring, Ref{Cint}), ncid[], varname, varid))

        xtype = Ref{Cint}()
        ndims = Ref{Cint}()
        f_inq_var = _nc_sym(:nc_inq_var)
        _nc_check(lib, ccall(f_inq_var, Cint,
                             (Cint, Cint, Ptr{UInt8}, Ref{Cint}, Ref{Cint},
                              Ptr{Cint}, Ptr{Cint}),
                             ncid[], varid[], C_NULL, xtype, ndims, C_NULL, C_NULL))

        dimids = Vector{Cint}(undef, ndims[])
        f_inq_vardimid = _nc_sym(:nc_inq_vardimid)
        _nc_check(lib, ccall(f_inq_vardimid, Cint,
                             (Cint, Cint, Ptr{Cint}), ncid[], varid[], dimids))

        dims = Vector{Int}(undef, ndims[])
        f_inq_dimlen = _nc_sym(:nc_inq_dimlen)
        for i in eachindex(dimids)
            len = Ref{Csize_t}()
            _nc_check(lib, ccall(f_inq_dimlen, Cint,
                                 (Cint, Cint, Ref{Csize_t}), ncid[], dimids[i], len))
            dims[i] = Int(len[])
        end

        T, getter = if xtype[] == _NC_BYTE
            Int8, :nc_get_var_schar
        elseif xtype[] == _NC_UBYTE
            UInt8, :nc_get_var_uchar
        elseif xtype[] == _NC_SHORT
            Int16, :nc_get_var_short
        elseif xtype[] == _NC_USHORT
            UInt16, :nc_get_var_ushort
        elseif xtype[] == _NC_INT
            Int32, :nc_get_var_int
        elseif xtype[] == _NC_UINT
            UInt32, :nc_get_var_uint
        elseif xtype[] == _NC_FLOAT
            Float32, :nc_get_var_float
        elseif xtype[] == _NC_DOUBLE
            Float64, :nc_get_var_double
        else
            error("Unsupported NetCDF variable type $(xtype[]) for $varname")
        end

        data = Vector{T}(undef, prod(dims))
        f_getter = _nc_sym(getter)
        _nc_check(lib, ccall(f_getter, Cint,
                             (Cint, Cint, Ptr{Cvoid}), ncid[], varid[], pointer(data)))
        return _reshape_c_order(data, dims)
    finally
        f_close = _nc_sym(:nc_close)
        ccall(f_close, Cint, (Cint,), ncid[])
    end
end

function _reshape_c_order(data::AbstractVector, dims::Vector{Int})
    if length(dims) <= 1
        return data
    end
    permutedims(reshape(data, Tuple(reverse(dims))))
end
