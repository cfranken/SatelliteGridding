@testset "Backend Resolution" begin
    @test supported_backend_names() == ("sequential", "cpu", "cuda", "metal")
    @test resolve_backend("sequential") === nothing
    @test resolve_backend("none") === nothing
    @test resolve_backend("") === nothing
    @test resolve_backend("cpu") isa KernelAbstractions.CPU
    @test resolve_backend("ka-cpu") isa KernelAbstractions.CPU

    @test_throws ErrorException resolve_backend("not-a-backend")

    help = join(backend_help_lines(), "\n")
    @test occursin("cuda", help)
    @test occursin("metal", help)
end
