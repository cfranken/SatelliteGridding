#!/usr/bin/env julia
#
# Benchmark: Sequential vs KA CPU vs KA CUDA vs Old Script
#
# Compares timing and correctness of all backends on TROPOMI SIF data.
#
# Usage:
#   julia --project=SatelliteGridding SatelliteGridding/bench/cuda_vs_old.jl
#
# Optional env vars:
#   SKIP_OLD=1     — skip the old gridL2_Dates.jl run
#   SKIP_CUDA=1    — skip the CUDA backend (if no GPU available)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using SatelliteGridding
using NCDatasets
using Statistics
using Dates
using Printf
using KernelAbstractions

const HAS_CUDA = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end

# --- Configuration ---
const DLAT      = 0.5f0
const DLON      = 0.5f0
const START     = "2019-07-01"
const STOP      = "2019-07-31"
const N_DAYS    = 31
const N_OVER    = 10
const TOML_CFG  = joinpath(@__DIR__, "..", "examples", "tropomi_sif.toml")
const JSON_CFG  = joinpath(@__DIR__, "..", "legacy", "jsonFiles", "tropomi_all.json")
const OLD_SCRIPT = joinpath(@__DIR__, "..", "legacy", "gridL2_Dates.jl")

const OUT_OLD   = "/tmp/bench_old.nc"
const OUT_SEQ   = "/tmp/bench_sequential.nc"
const OUT_CPU   = "/tmp/bench_cpu.nc"
const OUT_CUDA  = "/tmp/bench_cuda.nc"

const SKIP_OLD  = get(ENV, "SKIP_OLD", "0") == "1"
const SKIP_CUDA = get(ENV, "SKIP_CUDA", "0") == "1" || !HAS_CUDA

# --- Runners ---

function run_old_script()
    println("=" ^ 60)
    println("Backend: OLD gridL2_Dates.jl")
    println("=" ^ 60)
    flush(stdout)

    t = @elapsed begin
        cmd = `julia --project=$(joinpath(@__DIR__, "..")) $OLD_SCRIPT
               --Dict $JSON_CFG
               --startDate $START --stopDate $STOP
               --dDays $N_DAYS
               --dLat $DLAT --dLon $DLON
               -o $OUT_OLD`
        try
            run(cmd)
        catch e
            println("  WARNING: Old script failed (likely NCDatasets API change)")
            println("  ", e)
            flush(stdout)
            return NaN
        end
    end
    @printf("  Completed in %.1fs\n\n", t)
    flush(stdout)
    return t
end

function run_sequential()
    println("=" ^ 60)
    println("Backend: Sequential (Welford)")
    println("=" ^ 60)
    flush(stdout)

    config = load_config(TOML_CFG)
    grid_spec = GridSpec(lat_min=-90f0, lat_max=90f0,
                         lon_min=-180f0, lon_max=180f0,
                         dlat=DLAT, dlon=DLON)
    time_spec = TimeSpec(DateTime(START), DateTime(STOP), Dates.Day(N_DAYS))

    t = @elapsed grid_l2(config, grid_spec, time_spec;
                         n_oversample=N_OVER, outfile=OUT_SEQ)
    @printf("  Completed in %.1fs\n\n", t)
    flush(stdout)
    return t
end

function run_ka_cpu()
    println("=" ^ 60)
    println("Backend: KA CPU")
    println("=" ^ 60)
    flush(stdout)

    config = load_config(TOML_CFG)
    grid_spec = GridSpec(lat_min=-90f0, lat_max=90f0,
                         lon_min=-180f0, lon_max=180f0,
                         dlat=DLAT, dlon=DLON)
    time_spec = TimeSpec(DateTime(START), DateTime(STOP), Dates.Day(N_DAYS))

    t = @elapsed grid_l2(config, grid_spec, time_spec;
                         n_oversample=N_OVER, backend=CPU(),
                         outfile=OUT_CPU)
    @printf("  Completed in %.1fs\n\n", t)
    flush(stdout)
    return t
end

function run_ka_cuda()
    println("=" ^ 60)
    println("Backend: KA CUDA")
    println("=" ^ 60)
    println("  GPU: ", CUDA.name(CUDA.device()))
    flush(stdout)

    config = load_config(TOML_CFG)
    grid_spec = GridSpec(lat_min=-90f0, lat_max=90f0,
                         lon_min=-180f0, lon_max=180f0,
                         dlat=DLAT, dlon=DLON)
    time_spec = TimeSpec(DateTime(START), DateTime(STOP), Dates.Day(N_DAYS))

    t = @elapsed grid_l2(config, grid_spec, time_spec;
                         n_oversample=N_OVER, backend=CUDABackend(),
                         outfile=OUT_CUDA)
    @printf("  Completed in %.1fs\n\n", t)
    flush(stdout)
    return t
end

# --- Comparison ---

function compare_two(label_a, file_a, label_b, file_b)
    !isfile(file_a) && return nothing
    !isfile(file_b) && return nothing

    println("\n--- $label_a vs $label_b ---")
    ds_a = Dataset(file_a)
    ds_b = Dataset(file_b)

    # Weights (coalesce Missing → 0)
    n_a = Float64.(coalesce.(ds_a["n"][:, :, :], 0f0))
    n_b = Float64.(coalesce.(ds_b["n"][:, :, :], 0f0))

    mask = (n_a .> 0) .| (n_b .> 0)
    n_active = sum(mask)

    if n_active > 0
        w_diff = abs.(n_a .- n_b)
        max_w = maximum(w_diff[mask])
        mean_w = mean(w_diff[mask])
        rel_w = maximum(w_diff[mask] ./ max.(abs.(n_a[mask]), 1e-10))
        @printf("  Weights:  %d active cells, max_abs=%.4e, mean_abs=%.4e, max_rel=%.4e\n",
                n_active, max_w, mean_w, rel_w)
    else
        println("  No data in either file!")
        close(ds_a); close(ds_b)
        return false
    end

    all_pass = rel_w < 1e-3

    # Grid variables
    for varname in keys(ds_a)
        varname in ["lat", "lon", "time", "n"] && continue
        haskey(ds_b, varname) || continue

        v_a = Float64.(coalesce.(ds_a[varname][:, :, :], -999f0))
        v_b = Float64.(coalesce.(ds_b[varname][:, :, :], -999f0))

        if size(v_a) != size(v_b)
            @printf("  %-20s SKIP: shapes differ %s vs %s\n", varname, size(v_a), size(v_b))
            continue
        end

        valid = (v_a .> -900) .& (v_b .> -900)
        nv = sum(valid)
        if nv > 0
            d = abs.(v_a[valid] .- v_b[valid])
            max_d = maximum(d)
            mean_d = mean(d)
            denom = max.(abs.(v_a[valid]), 1e-10)
            rel_d = maximum(d ./ denom)
            pass = rel_d < 1e-3
            @printf("  %-20s %d valid, max_abs=%.4e, max_rel=%.4e  %s\n",
                    varname, nv, max_d, rel_d, pass ? "PASS" : "FAIL")
            all_pass &= pass
        else
            @printf("  %-20s no overlapping valid cells\n", varname)
        end
    end

    close(ds_a)
    close(ds_b)
    return all_pass
end

# --- Main ---

function main()
    println("=" ^ 60)
    println("CUDA vs Old Benchmark")
    println("=" ^ 60)
    println("  Config:     $TOML_CFG")
    println("  Date range: $START to $STOP ($N_DAYS days)")
    println("  Resolution: $(DLAT)deg x $(DLON)deg")
    println("  Oversample: $(N_OVER)x")
    println("  GPU:        $(HAS_CUDA ? "available" : "not available")")
    println("  Skip old:   $SKIP_OLD")
    println("  Skip CUDA:  $SKIP_CUDA")
    println()
    flush(stdout)

    timings = Dict{String,Float64}()

    # Run backends
    if !SKIP_OLD && isfile(OLD_SCRIPT) && isfile(JSON_CFG)
        timings["old"] = run_old_script()
    end

    timings["sequential"] = run_sequential()
    timings["ka_cpu"] = run_ka_cpu()

    if !SKIP_CUDA
        timings["ka_cuda"] = run_ka_cuda()
    end

    # Compare outputs
    println("\n" * "=" ^ 60)
    println("RESULT COMPARISON")
    println("=" ^ 60)

    if haskey(timings, "old") && !isnan(timings["old"])
        compare_two("Old", OUT_OLD, "Sequential", OUT_SEQ)
    end
    compare_two("Sequential", OUT_SEQ, "KA CPU", OUT_CPU)
    if !SKIP_CUDA
        compare_two("Sequential", OUT_SEQ, "KA CUDA", OUT_CUDA)
        compare_two("KA CPU", OUT_CPU, "KA CUDA", OUT_CUDA)
    end

    # Timing summary
    println("\n" * "=" ^ 60)
    println("TIMING SUMMARY")
    println("=" ^ 60)

    ref = timings["sequential"]
    for (name, t) in sort(collect(timings), by=x -> isnan(x[2]) ? Inf : x[2])
        if isnan(t)
            @printf("  %-15s   FAILED\n", name)
        else
            @printf("  %-15s %6.1fs  (%.2fx vs sequential)\n", name, t, ref / t)
        end
    end

    println("=" ^ 60)
    flush(stdout)

    # Cleanup
    for f in [OUT_OLD, OUT_SEQ, OUT_CPU, OUT_CUDA]
        rm(f, force=true)
    end
end

main()
