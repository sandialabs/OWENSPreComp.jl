# Developer Guide

This page documents the local conventions for changing OWENSPreComp without
breaking the tested PreComp compatibility workflow.

## Source Layout

| Path | Role |
| --- | --- |
| `src/OWENSPreComp.jl` | Module definition and exports. Currently exports `properties` and `tw_rate`. |
| `src/main.jl` | Section-property kernel, `Input` and `Output` types, geometry helpers, twist-rate helper, and material stiffness transforms. |
| `src/io.jl` | Legacy PreComp text-file readers. |
| `src/deprecated.jl` | Deprecated compatibility surface. |
| `test/runtests.jl` | Executable reference for file parsing, `Input` construction, validation behavior, and pinned numerical outputs. |

`Input`, `Output`, and the file readers are used in the documented workflow but
are not exported. Use module-qualified names in tests and examples unless the
public API is deliberately changed.

## Change Checklist

When changing the kernel or readers, check these points before trusting a
numerical result:

- units: chord and offsets are dimensional, airfoil coordinates and sector
  locations are normalized by chord, twist inputs are degrees, and the kernel
  converts to radians internally;
- signs: coupled stiffnesses, `y` offsets, and principal-axis angles must match
  the legacy wind-turbine output convention;
- geometry: upper and lower surfaces must remain single-valued in chordwise
  `x`, with a leading-edge first point at `(0, 0)`;
- material indexing: material ids in layups are one-based indices into the
  material-property arrays returned by `readmaterials`;
- web state: an empty web-location vector is the no-web case, and web locations
  for the legacy workflow are interpolated before reading the composite section;
- physical consistency: material Poisson ratios and transformed stiffness terms
  must remain physically realizable;
- automatic differentiation: new branches intended for gradient-based design
  should be covered by AD checks.

## Testing Strategy

Keep parser and kernel tests separable:

- reader tests should verify text-file ordering, type stability of returned
  arrays where relevant, material-id sorting, optional names, and airfoil
  upper/lower surface ordering;
- typed-input tests should construct `OWENSPreComp.Input` directly so a kernel
  failure is not hidden behind a parser change;
- fixture regression tests should compare named `Output` fields with explicit
  tolerances and include sign-sensitive quantities such as offsets and coupled
  stiffnesses;
- validation tests should pin error messages only when callers depend on them
  or when the message documents a clear input contract.

The current fixture compares all 16 stations from
`test01_composite_blade.out_gen`. It uses `5e-4` absolute tolerance for
three-decimal output fields and `5e-4` relative tolerance for scientific-format
stiffness, coupling, mass, and inertia fields.

## Documentation Updates

When behavior changes, update the docs in the same patch:

- `quickstart.md` for user workflow changes;
- `inputs_outputs.md` for field names, units, parser returns, or constructor
  expectations;
- `theory/frames_units.md` for assumptions, sign conventions, frames, or units;
- `validation.md` for new fixtures, changed tolerances, or known gaps;
- `reference/reference.md` if the public or documented helper surface changes.

Build the docs from the package docs environment when local dependencies are
available:

```julia
julia --project=OWENSPreComp.jl/docs OWENSPreComp.jl/docs/make.jl
```
