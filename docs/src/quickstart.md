# Quickstart

The preferred Julia workflow constructs an `OWENSPreComp.Input`, calls `properties`, and reads the typed fields of the resulting `OWENSPreComp.Output`.

## Install

Install the package from the repository in the active Julia environment:

```julia
using Pkg
Pkg.add(PackageSpec(url = "https://github.com/sandialabs/OWENSPreComp.jl.git"))
```

For toolkit development from a local checkout, develop the package path instead:

```julia
using Pkg
Pkg.develop(path = "/path/to/OWENSPreComp.jl")
```

The docs project points at the local package source through `docs/Project.toml`.

## Run The Tests

From the package directory:

```julia
using Pkg
Pkg.test()
```

The current suite exercises legacy file parsing, twist-rate calculation, input
validation, warning behavior, and section-property values against
PreComp-generated fixture output.

## Build One Section From Legacy Files

The regression fixture demonstrates the complete legacy input path. The main
file supplies span stations, leading-edge locations, chord, twist, airfoil-file
names, structural-layup-file names, and web end locations:

```julia
using OWENSPreComp

fixture = joinpath(pkgdir(OWENSPreComp), "test", "test01_composite_blade.pci")
(
    sloc,
    le_loc,
    chord,
    tw_aero_d,
    af_shape_file,
    int_str_file,
    ib_sp_stn,
    ob_sp_stn,
    ib_ch_loc,
    ob_ch_loc,
) = OWENSPreComp.readmain(fixture)
```

Use `tw_rate` to compute the twist derivative required by the kernel:

```julia
naf = length(sloc)
tw_prime_d = OWENSPreComp.tw_rate(naf, sloc, tw_aero_d)
```

For each station, read the airfoil profile, material table, and composite layup,
then build an `Input`:

```julia
i = 1
testdir = joinpath(pkgdir(OWENSPreComp), "test")

xu, yu, xl, yl = OWENSPreComp.readprecompprofile(joinpath(testdir, af_shape_file[i]))
e1, e2, g12, nu12, rho, material_name = OWENSPreComp.readmaterials(joinpath(testdir, "materials.inp"))

locW = Float64[]

(
    locU,
    n_laminaU,
    n_pliesU,
    tU,
    thetaU,
    mat_idxU,
    locL,
    n_laminaL,
    n_pliesL,
    tL,
    thetaL,
    mat_idxL,
    locW,
    n_laminaW,
    n_pliesW,
    tW,
    thetaW,
    mat_idxW,
) = OWENSPreComp.readcompositesection(joinpath(testdir, int_str_file[i]), locW)

input = OWENSPreComp.Input(
    chord[i],
    tw_aero_d[i],
    tw_prime_d[i],
    le_loc[i],
    vcat(xu, reverse(xl)[2:(end - 1)]),
    vcat(yu, reverse(yl)[2:(end - 1)]),
    e1,
    e2,
    g12,
    nu12,
    rho,
    locU,
    n_laminaU,
    float.(n_pliesU),
    tU,
    thetaU,
    mat_idxU,
    locL,
    n_laminaL,
    float.(n_pliesL),
    tL,
    thetaL,
    mat_idxL,
    locW,
    n_laminaW,
    float.(n_pliesW),
    tW,
    thetaW,
    mat_idxW,
)
```

Then evaluate the section:

```julia
output = OWENSPreComp.properties(input)
output.ei_flap
output.gj
output.mass
```

The full fixture includes webs from station 7 through station 16. See
`test/runtests.jl` for the web-location interpolation used when converting the
main input file into per-section `loc_web` values.

## Direct Julia Input

For coupled OWENS workflows, construct `OWENSPreComp.Input` directly from model
data instead of round-tripping through text files. The direct path avoids parser
format assumptions and is better suited to automatic differentiation because the
numeric arrays can be supplied by the caller.
