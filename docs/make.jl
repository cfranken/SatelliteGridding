using Documenter
using Literate

# Generate tutorial pages from Literate.jl scripts
LITERATE_SRC = joinpath(@__DIR__, "src")
LITERATE_OUT = joinpath(@__DIR__, "src", "generated")
mkpath(LITERATE_OUT)

for file in ["tutorial.jl"]
    Literate.markdown(joinpath(LITERATE_SRC, file), LITERATE_OUT;
                      documenter=true, credit=false)
end

# Load the package (assumes docs/Project.toml has dev'd the parent package)
using SatelliteGridding

makedocs(
    sitename = "SatelliteGridding.jl",
    modules  = [SatelliteGridding],
    pages = [
        "Home" => "index.md",
        "Tutorial" => "generated/tutorial.md",
        "Algorithm" => "algorithm.md",
        "Configuration" => "config.md",
        "CLI Reference" => "cli.md",
        "GPU Acceleration" => "gpu.md",
        "API Reference" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://cfranken.github.io/SatelliteGridding/",
    ),
    checkdocs = :exports,
    warnonly  = get(ENV, "DOCUMENTER_STRICT", "false") != "true",
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(repo = "github.com/cfranken/SatelliteGridding.git", devbranch = "master")
end
