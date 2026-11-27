@testset "Backend Resolution" begin
    @test supported_backend_names() == ("sequential", "cpu", "cuda", "metal")
    @test resolve_backend("sequential") === nothing
    @test resolve_backend("none") === nothing
    @test resolve_backend("") === nothing
    @test resolve_backend("cpu") isa KernelAbstractions.CPU
    @test resolve_backend("ka-cpu") isa KernelAbstractions.CPU

    @test_throws ErrorException resolve_backend("not-a-backend")

    spec = GridSpec()
    @test resolve_backend("cpu"; grid_spec=spec, accumulator_precision=Float32) isa CPU
    @test_throws ArgumentError("geometry_precision requires a spherical grid and a KA backend") resolve_backend(
        "cpu"; grid_spec=spec, geometry_precision=Float32)
    @test_throws ArgumentError("Metal requires Float32 accumulator_precision") resolve_backend(
        "metal"; grid_spec=spec, accumulator_precision=Float64)

    help = join(backend_help_lines(), "\n")
    @test occursin("cuda", help)
    @test occursin("metal", help)
end
