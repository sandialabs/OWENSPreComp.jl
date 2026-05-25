# Inputs and Outputs

OWENSPreComp has both a typed Julia interface and legacy file readers. New coupled code should prefer the typed interface because it is clearer, easier to test, and friendlier to automatic differentiation.

## Required Input Groups

| Group | Representative fields | Notes |
| --- | --- | --- |
| Span and airfoil | `chord`, `tw_aero_d`, `tw_prime_d`, `le_loc`, `xnode`, `ynode` | Airfoil coordinates describe the section perimeter. Twist values with `_d` are degrees. |
| Materials | `e1`, `e2`, `g12`, `anu12`, `rho`, `ath1`, `ath2` | Lengths must match the material count and material consistency is checked. |
| Upper/lower laminates | `xsec_nodeU`, `xsec_nodeL`, ply thickness, material id, fiber angle | Sector boundaries must be valid chordwise locations. |
| Web laminates | `loc_web`, `n_laminasW`, web ply data | Web fields are active only when webs are present. |
| Offsets and options | `precurve`, reference-axis location, output flags | These fields control section coordinate references and output format. |

The validation tests intentionally mutate these groups to confirm that mismatched lengths, invalid material properties, and invalid section geometry throw explicit errors.

## Output Fields

`OWENSPreComp.Output` stores the computed section properties. The most important quantities for OWENS coupling are:

| Field category | Meaning |
| --- | --- |
| axial, flap, lag, torsional stiffness | Beam stiffness terms such as `ea`, `eiflap`, `eilag`, and `gj`. |
| coupled stiffness | Cross-coupling terms for anisotropic laminate behavior. |
| inertia and mass | Section mass, mass moments, and principal inertias. |
| offsets | Shear-center, tension-center, and center-of-mass locations relative to the reference axes. |
| orientation | Principal inertia-axis orientation for downstream beam models. |

When adding tests, pin the actual typed fields needed by the consumer instead of checking only that the output is real-valued.

## File Readers

The file-based readers remain useful for compatibility:

- `readmain` reads the PreComp main input file and returns station and chord-location data.
- `readprecompprofile` reads upper and lower airfoil coordinates.
- `readmaterials` reads material-property tables.
- `readcompositesection` reads upper, lower, and web laminate schedules.

These readers are part of the documented workflow even though some are not exported. Use the module-qualified names in user code.
