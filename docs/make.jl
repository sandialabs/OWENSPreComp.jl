using Documenter, Literate, OWENSPreComp

# Build documentation
makedocs(;
    modules = [OWENSPreComp],
    pages = [
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Inputs and Outputs" => "inputs_outputs.md",
        "Theory, Frames, and Units" => joinpath("theory", "frames_units.md"),
        "Validation and Testing" => "validation.md",
        "Guide" => joinpath("guide", "precomp.md"),
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "OWENSPreComp.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
    remotes = nothing,
)

deploydocs(
    repo = "github.com/sandialabs/OWENSPreComp.jl.git",
)
