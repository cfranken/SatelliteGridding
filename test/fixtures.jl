# Numerical expectations below come from rectangle areas and elementary moments,
# independently of the package's geometry / accumulation helpers.
function _accuracy_input(file, rectangles, values; qa=ones(Float32, length(values)), version="B1",
                         vertex_first=false, bounds_group=nothing)
    Dataset(file, "c") do ds
        defDim(ds, "sounding", length(values))
        defDim(ds, "vertex", 4)
        lat = defVar(ds, "lat", Float32, ("sounding",))
        lon = defVar(ds, "lon", Float32, ("sounding",))
        group = bounds_group === nothing ? ds : defGroup(ds, bounds_group)
        dims = vertex_first ? ("vertex", "sounding") : ("sounding", "vertex")
        lb = defVar(group, "lat_bounds", Float32, dims)
        ob = defVar(group, "lon_bounds", Float32, dims)
        a = defVar(ds, "signal", Float32, ("sounding",); fillvalue=-9999f0,
                   attrib=["units" => "ppm", "long_name" => "Synthetic signal"])
        # Physical packed values = 2*signal + 7; raw integers = 4*signal.
        b = defVar(ds, "packed", Int16, ("sounding",); fillvalue=Int16(-32768),
                   attrib=["scale_factor" => 0.5f0, "add_offset" => 7f0])
        quality = defVar(ds, "qa", Float32, ("sounding",))
        for (i, (x0, x1, y0, y1)) in enumerate(rectangles)
            lat[i] = (y0 + y1) / 2
            lon[i] = (x0 + x1) / 2
            lat_bounds = Float32[y0, y0, y1, y1]
            lon_bounds = Float32[x0, x1, x1, x0]
            if vertex_first
                lb[:, i] = lat_bounds
                ob[:, i] = lon_bounds
            else
                lb[i, :] = lat_bounds
                ob[i, :] = lon_bounds
            end
        end
        a[:] = values
        b.var[:] = [isfinite(v) && v != -9999 ? Int16(4v) : Int16(-32768) for v in values]
        quality[:] = qa
        ds.attrib["product_version"] = version
        ds.attrib["institution"] = "Synthetic lab"
        group = defGroup(ds, "Metadata")
        defVar(group, "BuildId", String, ())[] = version * ".build"
    end
end

