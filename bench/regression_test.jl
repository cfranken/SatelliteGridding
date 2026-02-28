#!/usr/bin/env julia
#
# Regression Test: Old gridL2_Dates.jl vs New SatelliteGridding.grid_l2()
#
# Runs both codepaths on the same TROPOMI SIF data and compares outputs.
# Both use the SAME legacy JSON config to ensure identical data sources.
#
# Usage:
#   julia --project=. bench/regression_test.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using SatelliteGridding
using NCDatasets
using Statistics
using Dates
using Printf

# --- Configuration ---
const DLAT      = 0.5f0
const START     = "2019-07-01"
const STOP      = "2019-07-14"
const N_DAYS    = 14
const N_OVER    = 10

# Paths — both scripts use the SAME JSON config (same data source)
const OLD_SCRIPT = joinpath(@__DIR__, "..", "legacy", "gridL2_Dates.jl")
const JSON_CONFIG = joinpath(@__DIR__, "..", "legacy", "jsonFiles", "tropomi_all.json")
const OUT_OLD = "/tmp/regression_old_output.nc"
const OUT_NEW = "/tmp/regression_new_output.nc"

function run_old_script()
    println("=" ^ 60)
    println("Running OLD gridL2_Dates.jl")
    println("=" ^ 60)

    t_start = time()
    cmd = `julia --project=$(joinpath(@__DIR__, "..")) $OLD_SCRIPT
           --Dict $JSON_CONFIG
           --startDate $START --stopDate $STOP
           --dDays $N_DAYS
           --dLat $DLAT --dLon $DLAT
           -o $OUT_OLD`
    run(cmd)
    t_old = time() - t_start
    println("\nOld script completed in $(@sprintf("%.1f", t_old)) seconds")
    return t_old
end

function run_new_package()
    println("\n" * "=" ^ 60)
    println("Running NEW SatelliteGridding.grid_l2()")
    println("=" ^ 60)

    # Load the SAME JSON config as the old script
    config = load_config(JSON_CONFIG)
    grid_spec = GridSpec(lat_min=-90f0, lat_max=90f0,
                         lon_min=-180f0, lon_max=180f0,
                         dlat=DLAT, dlon=DLAT)
    time_spec = TimeSpec(DateTime(START), DateTime(STOP), Dates.Day(N_DAYS))

    t_start = time()
    grid_l2(config, grid_spec, time_spec;
            n_oversample=N_OVER,
            outfile=OUT_NEW)
    t_new = time() - t_start
    println("\nNew package completed in $(@sprintf("%.1f", t_new)) seconds")
    return t_new
end

function compare_outputs()
    println("\n" * "=" ^ 60)
    println("Comparing outputs")
    println("=" ^ 60)

    ds_old = Dataset(OUT_OLD)
    ds_new = Dataset(OUT_NEW)

    all_pass = true

    println("\nOld dims: ", [(string(d), s) for (d, s) in zip(keys(ds_old.dim), values(ds_old.dim))])
    println("New dims: ", [(string(d), s) for (d, s) in zip(keys(ds_new.dim), values(ds_new.dim))])

    # Compare weights (n)
    n_old = Float64.(ds_old["n"][:, :, :])
    n_new = Float64.(ds_new["n"][:, :, :])

    println("\n--- Weight comparison (n) ---")
    println("  Shape old: $(size(n_old)), new: $(size(n_new))")
    mask = (n_old .> 0) .| (n_new .> 0)
    n_active = sum(mask)

    if n_active > 0
        w_diff = abs.(n_old .- n_new)
        max_w = maximum(w_diff[mask])
        mean_w = mean(w_diff[mask])
        rel_w = maximum(w_diff[mask] ./ max.(abs.(n_old[mask]), 1e-10))
        println("  Active cells:       $n_active")
        println("  Total old weight:   $(@sprintf("%.1f", sum(n_old)))")
        println("  Total new weight:   $(@sprintf("%.1f", sum(n_new)))")
        println("  Max absolute diff:  $(@sprintf("%.6e", max_w))")
        println("  Mean absolute diff: $(@sprintf("%.6e", mean_w))")
        println("  Max relative diff:  $(@sprintf("%.6e", rel_w))")

        if rel_w > 1e-4
            println("  FAIL")
            all_pass = false
        else
            println("  PASS")
        end
    else
        println("  No data in either file!")
        all_pass = false
    end

    # Compare grid variables present in both files
    for varname in keys(ds_old)
        varname in ["lat", "lon", "time", "n"] && continue
        haskey(ds_new, varname) || continue

        println("\n--- Variable: $varname ---")
        v_old = Float64.(ds_old[varname][:, :, :])
        v_new = Float64.(ds_new[varname][:, :, :])

        if size(v_old) != size(v_new)
            println("  SKIP: shapes differ old=$(size(v_old)) new=$(size(v_new))")
            continue
        end

        valid = (v_old .> -900) .& (v_new .> -900)
        n_valid = sum(valid)
        n_old_only = sum((v_old .> -900) .& (v_new .<= -900))
        n_new_only = sum((v_old .<= -900) .& (v_new .> -900))

        println("  Valid in both: $n_valid, old only: $n_old_only, new only: $n_new_only")

        if n_valid > 0
            d = abs.(v_old[valid] .- v_new[valid])
            max_d = maximum(d)
            mean_d = mean(d)
            denom = max.(abs.(v_old[valid]), 1e-10)
            rel_d = maximum(d ./ denom)
            println("  Max absolute diff:  $(@sprintf("%.6e", max_d))")
            println("  Mean absolute diff: $(@sprintf("%.6e", mean_d))")
            println("  Max relative diff:  $(@sprintf("%.6e", rel_d))")

            if rel_d > 1e-4
                println("  FAIL")
                all_pass = false
            else
                println("  PASS")
            end
        end
    end

    close(ds_old)
    close(ds_new)
    return all_pass
end

# --- Main ---
println("Regression Test: Old vs New TROPOMI SIF Gridding")
println("  Config:       $JSON_CONFIG")
println("  Date range:   $START to $STOP ($N_DAYS days)")
println("  Resolution:   $(DLAT)°")
println("  Oversample:   $(N_OVER)x")
println()

t_old = run_old_script()
t_new = run_new_package()
all_pass = compare_outputs()

println("\n" * "=" ^ 60)
println("TIMING SUMMARY")
println("=" ^ 60)
println("  Old script:  $(@sprintf("%.1f", t_old))s")
println("  New package: $(@sprintf("%.1f", t_new))s")
println("  Speedup:     $(@sprintf("%.2f", t_old / t_new))x")
println()
if all_pass
    println("RESULT: PASS")
else
    println("RESULT: FAIL")
end
println("=" ^ 60)

# Cleanup
rm(OUT_OLD, force=true)
rm(OUT_NEW, force=true)
