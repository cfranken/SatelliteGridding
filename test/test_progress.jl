# Progress bars are drawn on request, or by default only when stderr is a terminal: in a
# file or a pipe (batch jobs, pipelines) every redraw would become new lines.
@testset "Progress display" begin
    # Every bar is made by _progress_bar, so this is where the option takes effect.
    @test !SatelliteGridding._progress_bar(3, "Time steps: ", false).enabled
    @test SatelliteGridding._progress_bar(3, "Time steps: ", true).enabled

    # Both entry points take the option, and default to drawing only in a terminal.
    for f in (grid_l2, grid_center)
        keywords = Base.kwarg_decl(only(methods(f)))
        @test :progress in keywords
        @test :verbose in keywords
    end
    mktemp() do _, io
        redirect_stderr(io) do
            @test SatelliteGridding._default_progress() == false  # stderr is a file here
        end
    end
    @test SatelliteGridding._default_progress() == (stderr isa Base.TTY)
end

module ProgressTestCLI
include(joinpath(@__DIR__, "..", "bin", "grid.jl"))
end

using Logging

@testset "CLI quiet and noProgress" begin
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        for day in ("01", "02", "03")
            _accuracy_input(joinpath(dir, "input_2020-01-$day.nc"), [(0.1, 0.9, 0.1, 0.9)], [2f0])
        end
        config = joinpath(dir, "input.toml")
        write(config, """
            filePattern = "input_YYYY-MM-DD.nc"
            folder = "$dir"
            [basic]
            lat = "lat"
            lon = "lon"
            lat_bnd = "lat_bounds"
            lon_bnd = "lon_bounds"
            [grid]
            signal = "signal"
            [center]
            min_count = 1
            """)
        outfile = joinpath(dir, "out.nc")
        args = ["--config", config, "--outFile", outfile,
                "--startDate", "2020-01-01", "--stopDate", "2020-01-03", "--dDays", "1",
                "--latMin", "0", "--latMax", "2", "--lonMin", "0", "--lonMax", "2"]
        for (command, parse) in (("l2", parse_l2_args), ("center", parse_center_args))
            @test parse([args; "--quiet"])["noProgress"]
            @test !parse([args; "--noProgress"])["quiet"]
            mktemp(dir) do _, io
                redirect_stdout(io) do
                    redirect_stderr(io) do
                        # Simulate a TTY default even though output is captured.
                        @test_logs min_level=Logging.Info ProgressTestCLI.main(
                            [command; args; "--quiet"]; default_progress=true)
                    end
                end
                @test isempty(read(seekstart(io), String))

                redirect_stdout(io) do
                    redirect_stderr(io) do
                        @test_logs (:info, r"Output (accumulator|file) dimension") (:info, r"Output written to") ProgressTestCLI.main(
                            [command; args; "--noProgress"]; default_progress=true)
                    end
                end
                @test occursin(r"^Total time: .* seconds\n$", read(seekstart(io), String))
            end
            Dataset(outfile) do ds
                @test maximum(ds["n"][:]) > 0
            end
        end
    end
end
