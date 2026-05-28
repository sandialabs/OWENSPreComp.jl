# Validation and Testing

OWENSPreComp has a package-level regression suite built around the bundled
legacy PreComp fixture in `test/`. The tests are the authoritative executable
example for the current file-reader-to-typed-input workflow.

## Current Test Evidence

| Area | Evidence |
| --- | --- |
| Legacy file input | `test/runtests.jl` reads `test01_composite_blade.pci`, the referenced airfoil files, `materials.inp`, and `int01.inp`. |
| Input-to-output workflow | The tests construct one `OWENSPreComp.Input` per station and call `OWENSPreComp.properties(input)`. |
| Twist rate | `tw_rate` is pinned for linear twist and used for the 16-station fixture. |
| Airfoil parsing | `TEtoTEdata`, `LEtoLEdata`, and `readprofile` are checked for upper/lower surface ordering. |
| Material parsing | `readmaterials` is checked for material-id sorting, numeric types, and optional material names. |
| Input validation | Mismatched vector lengths, invalid material properties, and invalid airfoil/sector geometry throw pinned messages. |
| Warning behavior | A leading-edge-aft-of-reference-axis warning is checked explicitly. |
| Numerical properties | Output fields are compared against PreComp-generated `.out_gen` fixture values with concrete tolerances. |

## Pinned Fixture Values

The first station of `test01_composite_blade.out_gen` is pinned as a simple
reference point:

| Quantity | Expected value |
| --- | --- |
| normalized span | `0.0` |
| chord | `0.664 m` |
| aerodynamic twist | `0.0 deg` |
| `ei_flap` | `0.5481e8 N m^2` |
| `ei_lag` | `0.2739e8 N m^2` |
| `gj` | `0.3003e8 N m^2` |
| `ea` | `0.8840e9 N` |
| `mass` | `0.1080e3 kg/m` |
| `flap_iner` | `0.6343e1 kg m` |
| `lag_iner` | `0.4066e1 kg m` |

For all 16 stations, the tests compare three-decimal fields with an absolute
tolerance of `5e-4` and scientific-notation fields with a relative tolerance of
`5e-4`.

## What The Fixture Covers

The fixture includes:

- 16 blade stations over a 20.15 m blade;
- root stations with no webs and outboard stations with two webs;
- changing chord, leading-edge location, aerodynamic twist, and airfoil files;
- five materials: unidirectional FRP, double-bias FRP, gelcoat, nexus, and
  balsa;
- three upper and three lower laminate sectors with flattened lamina schedules;
- web laminate stacks for the stations inside the web span.

This provides useful regression coverage for parser compatibility and broad
section-property behavior, but it is still one blade model. New physics changes
should add focused fixtures that isolate the intended behavior.

## Acceptance Rules For New Tests

- Pin named `Output` fields with numeric values and tolerances.
- Include at least one anisotropic or coupled-stiffness fixture for coupling-sensitive changes.
- Check offsets and signs in addition to stiffness magnitudes.
- Add AD checks for any new branch intended for optimization or gradient-based design.
- Keep legacy file-reader tests separate from typed `Input` tests so parser regressions and kernel regressions are easy to diagnose.
- Exercise both no-web and webbed stations when changing web-location handling.
- For unit or sign-convention changes, include a comment in the fixture that
  states the source convention and expected downstream convention.

## Known Documentation/API Gaps

`properties` and `tw_rate` are exported, but the tested typed workflow also uses module-qualified `Input`, `Output`, and reader helpers. A future API cleanup should either export the supported user-facing constructors/readers or explicitly mark them internal and provide a smaller public builder.
