using Documenter, Literate, OWENSPreComp

# Build documentation
makedocs(;
    modules = [OWENSPreComp],
    pages = [
        "Home" => "index.md",
        "Guide" => joinpath("guide", "precomp.md"),
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "OWENSPreComp.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
)

deploydocs(
    repo = "github.com/sandialabs/OWENSPreComp.jl.git",
)