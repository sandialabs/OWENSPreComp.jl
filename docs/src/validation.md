# Validation and Testing

OWENSPreComp has one of the stronger package-level regression suites in the toolkit. The main remaining work is to make the public docs and exported API match that test evidence.

## Current Test Evidence

| Area | Evidence |
| --- | --- |
| Legacy file input | `test/runtests.jl` reads PreComp-style main, airfoil, material, and composite-section files. |
| Twist rate | `tw_rate` is pinned for linear and fixture-based spanwise twist behavior. |
| Input validation | Mismatched vector lengths, invalid material properties, and invalid airfoil/sector geometry throw pinned messages. |
| Warning behavior | Leading-edge warning behavior is checked explicitly. |
| Numerical properties | Output fields are compared against PreComp-generated fixture files with concrete tolerances. |

## Acceptance Rules for New Tests

- Pin named `Output` fields with numeric values and tolerances.
- Include at least one anisotropic or coupled-stiffness fixture for coupling-sensitive changes.
- Check offsets and signs in addition to stiffness magnitudes.
- Add AD checks for any new branch intended for optimization or gradient-based design.
- Keep legacy file-reader tests separate from typed `Input` tests so parser regressions and kernel regressions are easy to diagnose.

## Known Documentation/API Gaps

`properties` and `tw_rate` are exported, but the tested typed workflow also uses module-qualified `Input`, `Output`, and reader helpers. A future API cleanup should either export the supported user-facing constructors/readers or explicitly mark them internal and provide a smaller public builder.
