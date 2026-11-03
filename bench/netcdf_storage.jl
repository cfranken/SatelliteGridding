# Compare lossless storage for a representative 3D rectangular gridded variable.
# Usage: julia --project=. bench/netcdf_storage.jl file.nc variable [start=1] [count=20]
using NCDatasets
using Printf

function write_sample(file, data, level, shuffle, width)
    Dataset(file, "c") do ds
        defDim(ds, "lon", size(data, 1))
        defDim(ds, "lat", size(data, 2))
        defDim(ds, "time", Inf)
        chunks = width == 0 ? nothing : (min(width, size(data, 1)), min(width, size(data, 2)), 1)
        v = defVar(ds, "sample", Float32, ("lon", "lat", "time");
                   deflatelevel=level, shuffle=shuffle, chunksizes=chunks, fillvalue=-999f0)
        for t in axes(data, 3)
            v[:, :, t] = data[:, :, t]
        end
    end
end

function main(args)
    length(args) in 2:4 || error("Usage: netcdf_storage.jl file.nc variable [start=1] [count=20]")
    start = length(args) >= 3 ? parse(Int, args[3]) : 1
    count = length(args) >= 4 ? parse(Int, args[4]) : 20
    count > 0 || error("count must be positive")
    data = Dataset(args[1]) do ds
        v = ds[args[2]]
        ndims(v) == 3 || error("Choose a rectangular gridded variable with dimensions lon, lat, time")
        1 <= start <= size(v, 3) || error("start is outside the time axis")
        last = min(start + count - 1, size(v, 3))
        values = v[:, :, start:last]
        println("Shape: ", size(values), "; populated cells: ", Base.count(!ismissing, values))
        Float32.(coalesce.(values, -999f0))
    end
    cases = [("uncompressed", 0, false, 1024),
             ("old automatic", 4, false, 0),
             ("default", 4, false, 1024),
             ("default + shuffle", 4, true, 1024),
             ("fast", 1, false, 1024),
             ("small", 9, false, 1024),
             ("small + shuffle", 9, true, 1024)]
    mktempdir() do dir
        file = joinpath(dir, "sample.nc")
        write_sample(file, data, 4, false, 1024) # warm compilation before timing
        @printf("%-20s %12s %10s\n", "Settings", "KiB", "Seconds")
        for (name, level, shuffle, width) in cases
            seconds = @elapsed write_sample(file, data, level, shuffle, width)
            Dataset(file) do ds
                isequal(ds["sample"].var[:, :, :], data) || error("Round-trip mismatch for $name")
            end
            @printf("%-20s %12.1f %10.3f\n", name, filesize(file) / 1024, seconds)
        end
    end
end

main(ARGS)
