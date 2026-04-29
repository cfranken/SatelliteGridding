"""
    read_nc_variable(dataset, path::String; bounds=false)

Read a variable from a NetCDF dataset, handling nested groups via `/` separators.

If `bounds=true`, reshape the result so that footprint vertices are in the last dimension,
returning an N×4 matrix (N soundings, 4 corners).
"""
function read_nc_variable(dataset, path::String; bounds::Bool=false)
    loc = split(path, "/")

    if length(loc) == 1
        return _read_flat_variable(dataset, path, bounds)
    end

    # Navigate nested groups
    gr = dataset.group[loc[1]]
    for i in 2:length(loc)-1
        gr = gr.group[loc[i]]
    end

    varname = loc[end]
    si = size(gr[varname])

    if bounds
        if si[1] == 4
            return reshape(gr[varname].var[:], 4, prod(si[2:end]))'
        elseif si[end] == 4
            return reshape(gr[varname].var[:], prod(si[1:end-1]), 4)
        end
    end

    return reshape(gr[varname].var[:], prod(si))
end

function _read_flat_variable(dataset, varname::String, bounds::Bool)
    if !bounds
        return dataset[varname].var[:]
    end
    si = size(dataset[varname])
    if si[1] == 4
        return reshape(dataset[varname].var[:], 4, prod(si[2:end]))'
    elseif si[end] == 4
        return reshape(dataset[varname].var[:], prod(si[1:end-1]), 4)
    end
    return dataset[varname].var[:]
end

"""
    read_variable_from_file(path, variable; bounds=false)

Read a variable from a NetCDF-like file. Falls back to the system `libnetcdf`
for flat HDF4/HDF-EOS variables when the Julia NetCDF build cannot open them.
"""
function read_variable_from_file(path::String, variable::String; bounds::Bool=false)
    try
        ds = Dataset(path)
        try
            return read_nc_variable(ds, variable; bounds=bounds)
        finally
            close(ds)
        end
    catch e
        if bounds
            rethrow(e)
        end
        try
            return read_system_netcdf_variable(path, variable)
        catch fallback_error
            error("Could not read variable '$variable' from '$path' with NCDatasets ($(sprint(showerror, e))) or system libnetcdf fallback ($(sprint(showerror, fallback_error)))")
        end
    end
end

"""
    read_array_from_file(path, variable)

Read a variable while preserving its native dimensionality.
"""
function read_array_from_file(path::String, variable::String)
    try
        ds = Dataset(path)
        try
            return _read_nc_array(ds, variable)
        finally
            close(ds)
        end
    catch e
        try
            return read_system_netcdf_variable(path, variable)
        catch fallback_error
            error("Could not read variable '$variable' from '$path' with NCDatasets ($(sprint(showerror, e))) or system libnetcdf fallback ($(sprint(showerror, fallback_error)))")
        end
    end
end

function _read_nc_array(dataset, path::String)
    loc = split(path, "/")
    if length(loc) == 1
        return _read_full_array(dataset[path])
    end

    gr = dataset.group[loc[1]]
    for i in 2:length(loc)-1
        gr = gr.group[loc[i]]
    end
    _read_full_array(gr[loc[end]])
end

function _read_full_array(var)
    selectors = ntuple(_ -> Colon(), ndims(var))
    var[selectors...]
end

"""
    read_nc_attribute(dataset, path::String, attr::String)

Read an attribute from a possibly group-nested variable.
"""
function read_nc_attribute(dataset, path::String, attr::String)
    loc = split(path, "/")

    if length(loc) == 1
        return dataset[path].attrib[attr]
    end

    gr = dataset.group[loc[1]]
    for i in 2:length(loc)-1
        gr = gr.group[loc[i]]
    end
    return gr[loc[end]].attrib[attr]
end

"""
    create_output_dataset(outfile, grid_spec, time_spec, grid_vars, compute_std) -> Dataset, Dict

Create the output NetCDF4 file with dimensions, coordinate variables, and data variable
definitions. Returns the dataset and a Dict mapping variable names to NCDatasets variables.
"""
function create_output_dataset(outfile::String, grid_spec::GridSpec,
                               n_times::Int, grid_vars::OrderedDict{String,String},
                               compute_std::Bool)
    ds = Dataset(outfile, "c")

    defDim(ds, "lon", length(grid_spec.lon))
    defDim(ds, "lat", length(grid_spec.lat))
    defDim(ds, "time", n_times)

    ds_lat = defVar(ds, "lat", Float32, ("lat",),
                    attrib=["units" => "degrees_north", "long_name" => "Latitude"])
    ds_lon = defVar(ds, "lon", Float32, ("lon",),
                    attrib=["units" => "degrees_east", "long_name" => "Longitude"])
    ds_time = defVar(ds, "time", Float32, ("time",),
                     attrib=["units" => "days since 1970-01-01",
                              "long_name" => "Time (UTC), start of interval"])

    ds_lat[:] = grid_spec.lat
    ds_lon[:] = grid_spec.lon

    ds.attrib["title"] = "Gridded satellite data"
    ds.attrib["created_with"] = "SatelliteGridding.jl"

    # Define output data variables
    nc_vars = Dict{String,Any}()
    for (key, _) in grid_vars
        nc_vars[key] = defVar(ds, key, Float32, ("time", "lon", "lat"),
                              deflatelevel=4, fillvalue=-999.0f0)
        if compute_std
            key_std = key * "_std"
            nc_vars[key_std] = defVar(ds, key_std, Float32, ("time", "lon", "lat"),
                                      deflatelevel=4, fillvalue=-999.0f0)
        end
    end

    nc_vars["n"] = defVar(ds, "n", Float32, ("time", "lon", "lat"),
                          deflatelevel=4, fillvalue=-999.0f0,
                          attrib=["units" => "", "long_name" => "Number of pixels in average"])

    return ds, nc_vars
end

"""
    copy_variable_attributes!(nc_vars, dataset, grid_vars)

Copy attributes (units, long_name, etc.) from the first input file to the output variables.
"""
function copy_variable_attributes!(nc_vars::Dict, dataset, grid_vars::OrderedDict{String,String})
    attrib_names = ["units", "long_name", "valid_range", "description", "unit", "longname"]
    for (key, value) in grid_vars
        for at in attrib_names
            try
                nc_vars[key].attrib[at] = read_nc_attribute(dataset, value, at)
            catch
            end
        end
    end
end
