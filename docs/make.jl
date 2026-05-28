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
        "Developer Guide" => "developer.md",
        "Legacy PreComp Guide" => joinpath("guide", "precomp.md"),
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "OWENSPreComp.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
    remotes = nothing,
    format = Documenter.HTML(
        repolink = "https://github.com/sandialabs/OWENSPreComp.jl",
        edit_link = "master",
    ),
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/sandialabs/OWENSPreComp.jl.git",
        devbranch = "master",
    )
end
