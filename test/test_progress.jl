# Progress bars are drawn on request, or by default only when stderr is a terminal: in a
# file or a pipe (batch jobs, pipelines) every redraw would become new lines.
@testset "Progress display" begin
    # Every bar is made by _progress_bar, so this is where the option takes effect.
    @test !SatelliteGridding._progress_bar(3, "Time steps: ", false).enabled
    @test SatelliteGridding._progress_bar(3, "Time steps: ", true).enabled

    # Both entry points take the option, and default to drawing only in a terminal.
    for f in (grid_l2, grid_center)
        @test :progress in Base.kwarg_decl(only(methods(f)))
    end
    mktemp() do _, io
        redirect_stderr(io) do
            @test SatelliteGridding._default_progress() == false  # stderr is a file here
        end
    end
    @test SatelliteGridding._default_progress() == (stderr isa Base.TTY)
end
