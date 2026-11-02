# In-suite ordering / layout test for the Zarr daily emit. The full cross-language
# gate (xarray + zarrita) lives in test/zarr_probe/run.sh; this guards the Julia
# transform, the reversed-dim on-disk layout, idempotency, and SEM derivation.
# Requires Zarr (in the test target); `using Zarr` at the top of runtests.jl loads
# the SatelliteGriddingZarrExt extension.

@testset "Zarr daily emit — ordering, layout, SEM" begin
    dir = mktempdir()
    path = joinpath(dir, "probe.zarr")
    gs = GridSpec()                       # 1° global: 180 lat × 360 lon
    gv = OrderedDict("probe" => "probe")

    # In-memory [lon, lat] south-up ramp engineered so the logical store cell
    # (iy, ix) (0-based, iy=0 NORTH, ix=0 lon=−180) holds iy*1000 + ix.
    G_mem = Float32[(180 - iy1) * 1000 + (ix1 - 1) for ix1 in 1:360, iy1 in 1:180]
    grid_data = reshape(G_mem, 360, 180, 1)
    grid_w = ones(Float32, 360, 180)
    grid_std = zeros(Float32, 360, 180, 1)

    g = init_zarr_store(path, gs, gv; compute_std = false)
    ti = write_zarr_day!(g, Date(2021, 1, 1), grid_data, grid_std, grid_w, gv, false)
    consolidate_zarr!(path)
    @test ti == 1                          # epoch day → absolute index 1

    # Reopen. Julia layout is reversed (NX,NY,NT); logical A[ti,iy,ix] == julia z[ix+1,iy+1,ti+1].
    z = zopen(path, "r")
    P = z["probe"]
    @test size(P) == (360, 180, 3653)
    for (iy, ix) in [(0, 0), (44, 59), (89, 180), (179, 359)]
        @test P[ix + 1, iy + 1, 1] == iy * 1000 + ix
    end
    @test z["lat"][1] ≈ 89.5           # north-up, descending
    @test z["lat"][end] ≈ -89.5
    @test z["lon"][1] ≈ -179.5         # ascending east
    @test z["lon"][end] ≈ 179.5

    # On-disk metadata must read as logical (time,lat,lon) for xarray/zarrita.
    zmeta = SatelliteGridding.JSON.parsefile(joinpath(path, "probe", ".zarray"))
    @test zmeta["shape"] == [3653, 180, 360]
    @test zmeta["chunks"] == [1, 180, 360]
    @test zmeta["order"] == "C"
    zattrs = SatelliteGridding.JSON.parsefile(joinpath(path, "probe", ".zattrs"))
    @test zattrs["_ARRAY_DIMENSIONS"] == ["time", "lat", "lon"]

    # _n present (Int16); _sem absent when compute_std=false.
    @test isdir(joinpath(path, "probe_n"))
    @test !isdir(joinpath(path, "probe_sem"))

    # .zmetadata consolidated, with _ARRAY_DIMENSIONS preserved.
    cons = SatelliteGridding.JSON.parsefile(joinpath(path, ".zmetadata"))
    @test cons["zarr_consolidated_format"] == 1
    @test haskey(cons["metadata"], "probe/.zarray")

    # Idempotent re-write of the same date → same absolute chunk index.
    @test write_zarr_day!(g, Date(2021, 1, 1), grid_data, grid_std, grid_w, gv, false) == 1

    @testset "SEM derivation + Int16 counts" begin
        path2 = joinpath(dir, "probe2.zarr")
        grid_std2 = fill(2.0f0, 360, 180, 1)      # std = 2 everywhere
        grid_w2 = fill(4.0f0, 360, 180)           # n = 4 → sem = 2/√4 = 1
        g2 = init_zarr_store(path2, gs, gv; compute_std = true)
        write_zarr_day!(g2, Date(2021, 1, 1), grid_data, grid_std2, grid_w2, gv, true)
        z2 = zopen(path2, "r")
        @test z2["probe_sem"][1, 1, 1] ≈ 1.0f0
        @test z2["probe_n"][1, 1, 1] == Int16(4)
    end

    @testset "empty cells → NaN / −1 fill" begin
        path3 = joinpath(dir, "probe3.zarr")
        w_sparse = zeros(Float32, 360, 180)
        w_sparse[10, 20] = 3.0f0                  # one populated cell
        gd = fill(7.0f0, 360, 180, 1)
        g3 = init_zarr_store(path3, gs, gv; compute_std = false)
        write_zarr_day!(g3, Date(2021, 1, 1), gd, zeros(Float32, 360, 180, 1), w_sparse, gv, false)
        z3 = zopen(path3, "r")
        @test z3["probe"][10, 180 - 20 + 1, 1] == 7.0f0     # lat-flip: in-mem row 20 → north-up row 161
        @test z3["probe_n"][10, 180 - 20 + 1, 1] == Int16(3)
        @test isnan(z3["probe"][1, 1, 1])                    # unpopulated → NaN
        @test z3["probe_n"][1, 1, 1] == Int16(-1)            # unpopulated → −1
    end

    @testset "non-1° rectangular grid (store dims, no hard-coded 360×180)" begin
        p = joinpath(dir, "deg2.zarr")
        gs2 = GridSpec(dlat = 2.0f0, dlon = 2.0f0)        # 90 lat × 180 lon
        gv2 = OrderedDict("v" => "v")
        nx, ny = length(gs2.lon), length(gs2.lat)
        @test (nx, ny) == (180, 90)
        gd = Float32[(ny - iy1) * 1000 + (ix1 - 1) for ix1 in 1:nx, iy1 in 1:ny]
        g4 = init_zarr_store(p, gs2, gv2; compute_std = false)
        write_zarr_day!(g4, Date(2021, 1, 1), reshape(gd, nx, ny, 1),
                        zeros(Float32, nx, ny, 1), ones(Float32, nx, ny), gv2, false)
        z4 = zopen(p, "r")
        @test size(z4["v"]) == (nx, ny, 3653)
        @test z4["lat"][1] ≈ 89.0       # 2° centers: 89, 87, …, −89 (descending)
        @test z4["lat"][end] ≈ -89.0
        for (iy, ix) in [(0, 0), (10, 59), (ny - 1, nx - 1)]
            @test z4["v"][ix + 1, iy + 1, 1] == iy * 1000 + ix   # logical A[ti,iy,ix]
        end
    end

    @testset "custom epoch bounds against the store's own time axis" begin
        p = joinpath(dir, "ep.zarr")
        gv3 = OrderedDict("v" => "v")
        ep = Date(2020, 1, 1)
        hz = Date(2022, 1, 1)
        g5 = init_zarr_store(p, GridSpec(), gv3; epoch = ep, horizon = hz, compute_std = false)
        NTexp = Dates.value(hz - ep) + 1
        @test size(g5["time"])[1] == NTexp
        z = zeros(Float32, 360, 180, 1); w = ones(Float32, 360, 180)
        ti = write_zarr_day!(g5, Date(2021, 12, 31), z, z, w, gv3, false; epoch = ep)
        @test ti == Dates.value(Date(2021, 12, 31) - ep) + 1
        # Beyond the custom horizon → clear error (previously slipped past the ZARR_NT constant).
        @test_throws ErrorException write_zarr_day!(g5, Date(2025, 1, 1), z, z, w, gv3, false; epoch = ep)
    end

    @testset "errors without Zarr loaded would surface helpfully" begin
        # Sanity: the parent delegator resolves to the extension (Zarr is loaded here).
        @test init_zarr_store isa Function
        @test write_zarr_day! isa Function
    end
end
