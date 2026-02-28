#!/usr/bin/env julia
"""
GPU vs CPU benchmark for the SatelliteGridding oversampling pipeline.
Uses native CUDA.jl kernels for sort + subpixel computation on GPU,
then downloads indices and does cache-friendly Welford accumulation on CPU.
"""

println("Loading packages...")
flush(stdout)
t_load = @elapsed begin
    using SatelliteGridding, Dates, NCDatasets, OrderedCollections, Statistics
    using CUDA
end
println("  Loaded in $(round(t_load, digits=1))s")
println("GPU: ", CUDA.name(CUDA.device()))
flush(stdout)

# ─── Native CUDA Kernels ─────────────────────────────────────────────────────

function gpu_sort_ccw_kernel!(lat_mat, lon_mat, n)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if i > n; return; end
    @inbounds begin
        clat = (lat_mat[i,1] + lat_mat[i,2] + lat_mat[i,3] + lat_mat[i,4]) / 4
        clon = (lon_mat[i,1] + lon_mat[i,2] + lon_mat[i,3] + lon_mat[i,4]) / 4
        a1 = atan(lat_mat[i,1] - clat, lon_mat[i,1] - clon)
        a2 = atan(lat_mat[i,2] - clat, lon_mat[i,2] - clon)
        a3 = atan(lat_mat[i,3] - clat, lon_mat[i,3] - clon)
        a4 = atan(lat_mat[i,4] - clat, lon_mat[i,4] - clon)
        la1 = lat_mat[i,1]; la2 = lat_mat[i,2]; la3 = lat_mat[i,3]; la4 = lat_mat[i,4]
        lo1 = lon_mat[i,1]; lo2 = lon_mat[i,2]; lo3 = lon_mat[i,3]; lo4 = lon_mat[i,4]
        if a1 > a2; a1,a2 = a2,a1; la1,la2 = la2,la1; lo1,lo2 = lo2,lo1; end
        if a3 > a4; a3,a4 = a4,a3; la3,la4 = la4,la3; lo3,lo4 = lo4,lo3; end
        if a1 > a3; a1,a3 = a3,a1; la1,la3 = la3,la1; lo1,lo3 = lo3,lo1; end
        if a2 > a4; a2,a4 = a4,a2; la2,la4 = la4,la2; lo2,lo4 = lo4,lo2; end
        if a2 > a3; a2,a3 = a3,a2; la2,la3 = la3,la2; lo2,lo3 = lo3,lo2; end
        lat_mat[i,1] = la1; lat_mat[i,2] = la2; lat_mat[i,3] = la3; lat_mat[i,4] = la4
        lon_mat[i,1] = lo1; lon_mat[i,2] = lo2; lon_mat[i,3] = lo3; lon_mat[i,4] = lo4
    end
    return nothing
end

function gpu_subpixel_floor_kernel!(ix_out, iy_out, lat_mat, lon_mat,
                                    lat_min, lat_range, lon_min, lon_range,
                                    n_lat, n_lon, N, n_fp)
    fp = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if fp > n_fp; return; end
    @inbounds begin
        il1 = (lat_mat[fp,1] - lat_min) / lat_range * n_lat + 1.0f0
        il2 = (lat_mat[fp,2] - lat_min) / lat_range * n_lat + 1.0f0
        il3 = (lat_mat[fp,3] - lat_min) / lat_range * n_lat + 1.0f0
        il4 = (lat_mat[fp,4] - lat_min) / lat_range * n_lat + 1.0f0
        io1 = (lon_mat[fp,1] - lon_min) / lon_range * n_lon + 1.0f0
        io2 = (lon_mat[fp,2] - lon_min) / lon_range * n_lon + 1.0f0
        io3 = (lon_mat[fp,3] - lon_min) / lon_range * n_lon + 1.0f0
        io4 = (lon_mat[fp,4] - lon_min) / lon_range * n_lon + 1.0f0

        dlat_02 = (il2 - il1) / (2.0f0 * N)
        dlon_02 = (io2 - io1) / (2.0f0 * N)
        dlat_43 = (il3 - il4) / (2.0f0 * N)
        dlon_43 = (io3 - io4) / (2.0f0 * N)

        for i in Int32(1):N
            lat_0 = il1 + dlat_02 * (2.0f0 * i - 1.0f0)
            lon_0 = io1 + dlon_02 * (2.0f0 * i - 1.0f0)
            lat_1 = il4 + dlat_43 * (2.0f0 * i - 1.0f0)
            lon_1 = io4 + dlon_43 * (2.0f0 * i - 1.0f0)

            dlat_c = (lat_1 - lat_0) / (2.0f0 * N)
            dlon_c = (lon_1 - lon_0) / (2.0f0 * N)
            for j in Int32(1):N
                sub_lat = lat_0 + dlat_c * (2.0f0 * j - 1.0f0)
                sub_lon = lon_0 + dlon_c * (2.0f0 * j - 1.0f0)
                k = (fp - Int32(1)) * N * N + (i - Int32(1)) * N + j
                ix_out[k] = floor(Int32, sub_lat)
                iy_out[k] = floor(Int32, sub_lon)
            end
        end
    end
    return nothing
end

# ─── Warmup ──────────────────────────────────────────────────────────────────

function warmup_gpu(N::Int32)
    println("Warming up GPU kernels...")
    flush(stdout)
    t = @elapsed begin
        d_lat = CUDA.rand(Float32, 1000, 4)
        d_lon = CUDA.rand(Float32, 1000, 4)
        threads = 256; blocks = cld(1000, threads)
        @cuda threads=threads blocks=blocks gpu_sort_ccw_kernel!(d_lat, d_lon, Int32(1000))
        CUDA.synchronize()

        d_ix = CUDA.zeros(Int32, 1000 * N * N)
        d_iy = CUDA.zeros(Int32, 1000 * N * N)
        @cuda threads=threads blocks=blocks gpu_subpixel_floor_kernel!(
            d_ix, d_iy, d_lat, d_lon,
            -90.0f0, 180.0f0, -180.0f0, 360.0f0,
            360.0f0, 720.0f0, N, Int32(1000))
        CUDA.synchronize()
    end
    println("  Warmup done in $(round(t, digits=1))s\n")
    flush(stdout)
end

# ─── Benchmark ───────────────────────────────────────────────────────────────

function benchmark()
    config = load_config("examples/tropomi_sif.toml")
    grid_spec = GridSpec(lat_min=-90f0, lat_max=90f0, lon_min=-180f0, lon_max=180f0,
                         dlat=0.5f0, dlon=0.5f0)

    d = DateTime("2019-07-01")
    end_date = d + Dates.Day(13)
    files = SatelliteGridding.find_files_range(config.file_pattern, config.folder, d, end_date)

    n_lon = length(grid_spec.lon)
    n_lat = length(grid_spec.lat)
    n_vars = length(config.grid_vars)
    T = Float32
    N = Int32(10)
    Ni = Int(N)

    warmup_gpu(N)

    println("=== Benchmark: $(length(files)) files, n_oversample=$N ===\n")
    flush(stdout)

    total_gpu_compute = 0.0
    total_accum_gpu_path = 0.0
    total_cpu_full = 0.0
    total_io = 0.0
    total_soundings = 0

    # Grid accumulators
    grid_data_cpu = zeros(T, n_lon, n_lat, n_vars)
    grid_std_cpu = zeros(T, n_lon, n_lat, n_vars)
    grid_w_cpu = zeros(T, n_lon, n_lat)
    grid_data_gpu = zeros(T, n_lon, n_lat, n_vars)
    grid_std_gpu = zeros(T, n_lon, n_lat, n_vars)
    grid_w_gpu = zeros(T, n_lon, n_lat)

    # CPU buffers
    points_buf = zeros(T, Ni, Ni, 2)
    ix_buf = zeros(Int32, Ni^2)
    iy_buf = zeros(Int32, Ni^2)
    lats0 = zeros(T, Ni); lons0 = zeros(T, Ni)
    lats1 = zeros(T, Ni); lons1 = zeros(T, Ni)

    # Second set of buffers for GPU path (so CPU and GPU don't share mutable state)
    points_buf2 = zeros(T, Ni, Ni, 2)
    ix_buf2 = zeros(Int32, Ni^2)
    iy_buf2 = zeros(Int32, Ni^2)
    lats0_2 = zeros(T, Ni); lons0_2 = zeros(T, Ni)
    lats1_2 = zeros(T, Ni); lons1_2 = zeros(T, Ni)

    lat_range = T(grid_spec.lat_max - grid_spec.lat_min)
    lon_range = T(grid_spec.lon_max - grid_spec.lon_min)
    threads = 256

    for (fi, filepath) in enumerate(files)
        # --- IO ---
        t_io = @elapsed begin
            fin = Dataset(filepath)
            lat_bnd = read_nc_variable(fin, config.basic["lat_bnd"]; bounds=true)
            lon_bnd = read_nc_variable(fin, config.basic["lon_bnd"]; bounds=true)
            ns = size(lat_bnd, 1)
            mat_in = zeros(T, ns, n_vars)
            for (co, (_, value)) in enumerate(config.grid_vars)
                mat_in[:, co] = read_nc_variable(fin, value)
            end
            idx = apply_filters(fin, config, lat_bnd, lon_bnd, grid_spec)
            close(fin)
        end
        total_io += t_io

        np = length(idx)
        if np == 0; continue; end
        total_soundings += np

        lat_f = T.(lat_bnd[idx, :])
        lon_f = T.(lon_bnd[idx, :])
        vals_f = mat_in[idx, :]

        # === CPU path (reference) ===
        t_cpu = @elapsed begin
            lat_cpu = copy(lat_f)
            lon_cpu = copy(lon_f)
            sort_corners_ccw!(lat_cpu, lon_cpu)
            ilat_cpu = ((lat_cpu .- grid_spec.lat_min) / lat_range * T(n_lat)) .+ T(1)
            ilon_cpu = ((lon_cpu .- grid_spec.lon_min) / lon_range * T(n_lon)) .+ T(1)
            accumulate_footprint!(grid_data_cpu, grid_std_cpu, grid_w_cpu, true,
                                  ilat_cpu, ilon_cpu, vals_f,
                                  np, n_vars, Ni,
                                  points_buf, ix_buf, iy_buf,
                                  lats0, lons0, lats1, lons1)
        end
        total_cpu_full += t_cpu

        # === GPU path: sort on GPU, then use same CPU accumulation logic ===
        np32 = Int32(np)
        blocks = cld(np, threads)

        # GPU sort + download sorted corners
        t_gpu_sort = @elapsed begin
            d_lat = CuArray(lat_f)
            d_lon = CuArray(lon_f)
            @cuda threads=threads blocks=blocks gpu_sort_ccw_kernel!(d_lat, d_lon, np32)
            CUDA.synchronize()
            lat_gpu = Array(d_lat)
            lon_gpu = Array(d_lon)
            CUDA.unsafe_free!(d_lat); CUDA.unsafe_free!(d_lon)
        end
        total_gpu_compute += t_gpu_sort

        # Use same CPU accumulation as reference (index conversion + subpixels + Welford)
        t_accum = @elapsed begin
            ilat_gpu = ((lat_gpu .- grid_spec.lat_min) / lat_range * T(n_lat)) .+ T(1)
            ilon_gpu = ((lon_gpu .- grid_spec.lon_min) / lon_range * T(n_lon)) .+ T(1)
            accumulate_footprint!(grid_data_gpu, grid_std_gpu, grid_w_gpu, true,
                                  ilat_gpu, ilon_gpu, vals_f,
                                  np, n_vars, Ni,
                                  points_buf2, ix_buf2, iy_buf2,
                                  lats0_2, lons0_2, lats1_2, lons1_2)
        end
        total_accum_gpu_path += t_accum

        # Detailed GPU kernel timing on file 2
        if fi == 2
            d_lat2 = CuArray(lat_f)
            d_lon2 = CuArray(lon_f)
            CUDA.synchronize()
            t_sort_only = @elapsed begin
                @cuda threads=threads blocks=blocks gpu_sort_ccw_kernel!(d_lat2, d_lon2, np32)
                CUDA.synchronize()
            end
            d_ix2 = CUDA.zeros(Int32, np * Ni * Ni)
            d_iy2 = CUDA.zeros(Int32, np * Ni * Ni)
            t_sub_only = @elapsed begin
                @cuda threads=threads blocks=blocks gpu_subpixel_floor_kernel!(
                    d_ix2, d_iy2, d_lat2, d_lon2,
                    T(grid_spec.lat_min), lat_range,
                    T(grid_spec.lon_min), lon_range,
                    T(n_lat), T(n_lon), N, np32)
                CUDA.synchronize()
            end
            CUDA.unsafe_free!(d_lat2); CUDA.unsafe_free!(d_lon2)
            CUDA.unsafe_free!(d_ix2); CUDA.unsafe_free!(d_iy2)
            println("  GPU kernel timing (file 2, $(np) soundings):")
            println("    Sort CCW:       $(round(t_sort_only*1000, digits=1))ms")
            println("    Subpixel+floor: $(round(t_sub_only*1000, digits=1))ms")
            flush(stdout)
        end

        println("File $fi/$(length(files)): $(np) soundings  " *
                "CPU=$(round(t_cpu*1000, digits=0))ms  " *
                "GPU_sort+xfer=$(round(t_gpu_sort*1000, digits=0))ms  " *
                "accum=$(round(t_accum*1000, digits=0))ms")
        flush(stdout)
    end

    gpu_total = total_gpu_compute + total_accum_gpu_path
    cpu_pipeline = total_io + total_cpu_full
    gpu_pipeline = total_io + gpu_total

    println("\n" * "="^60)
    println("Results: $(length(files)) files, $(round(total_soundings/1e6, digits=1))M soundings")
    println("="^60)
    println("  IO (read+filter):        $(round(total_io, digits=2))s")
    println()
    println("  CPU path (full):         $(round(total_cpu_full, digits=2))s")
    println()
    println("  GPU path breakdown:")
    println("    GPU compute + xfer:    $(round(total_gpu_compute, digits=2))s")
    println("    Welford accum (CPU):   $(round(total_accum_gpu_path, digits=2))s")
    println("    GPU path total:        $(round(gpu_total, digits=2))s")
    println()
    println("  Full CPU pipeline:       $(round(cpu_pipeline, digits=2))s")
    println("  Full GPU pipeline:       $(round(gpu_pipeline, digits=2))s")
    println("  Speedup:                 $(round(cpu_pipeline / gpu_pipeline, digits=2))x")

    # Verify results match
    mask = grid_w_cpu .> 0
    ncells = sum(mask)
    if ncells > 0
        max_w_diff = maximum(abs.(grid_w_cpu[mask] - grid_w_gpu[mask]) ./ max.(grid_w_cpu[mask], T(1e-10)))
        max_d_diff = T(0)
        for z in 1:n_vars
            dc = grid_data_cpu[:,:,z]
            dg = grid_data_gpu[:,:,z]
            valid = mask .& (abs.(dc) .> T(1e-10))
            if sum(valid) > 0
                md = maximum(abs.(dc[valid] - dg[valid]) ./ abs.(dc[valid]))
                max_d_diff = max(max_d_diff, md)
            end
        end
        println()
        println("  Verification ($ncells grid cells with data):")
        println("  Max relative weight diff: $(round(max_w_diff * 100, sigdigits=3))%")
        println("  Max relative data diff:   $(round(max_d_diff * 100, sigdigits=3))%")
        match = max_w_diff < 0.01 && max_d_diff < 0.01
        println("  Results match: $(match ? "YES" : "NO")")
    end
    flush(stdout)
end

benchmark()
