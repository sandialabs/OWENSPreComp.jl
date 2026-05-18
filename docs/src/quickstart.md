# Quickstart

The preferred Julia workflow constructs an `OWENSPreComp.Input`, calls `properties`, and reads the typed fields of the resulting `OWENSPreComp.Output`.

## Install

```julia
using Pkg
Pkg.add(PackageSpec(url = "https://github.com/sandialabs/OWENSPreComp.jl.git"))
```

For toolkit development from sibling checkouts:

```julia
using Pkg
Pkg.develop(path = "../OWENSPreComp.jl")
```

The docs project uses local package sources and should be built with Julia 1.11
or newer.

## Run the Tested File Workflow

The regression tests demonstrate the full legacy input path:

```julia
using OWENSPreComp

fixture = joinpath(pkgdir(OWENSPreComp), "test", "test01_composite_blade.pci")
main_input = OWENSPreComp.readmain(fixture)
```

The tests then read the airfoil profile, material file, and composite-section schedule, construct an `OWENSPreComp.Input` for each station, and call:

```julia
output = OWENSPreComp.properties(input)
```

Use `test/runtests.jl` as the canonical executable example until a compact public fixture is added to the docs.

## First Checks

After installation, run:

```julia
using Pkg
Pkg.test("OWENSPreComp")
```

The current suite pins twist-rate behavior, input validation errors, warning behavior, file parsing, and numerical section-property values against PreComp-generated fixture output.
