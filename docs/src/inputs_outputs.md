# Inputs and Outputs

OWENSPreComp has both a typed Julia interface and legacy file readers. New coupled code should prefer the typed interface because it is clearer, easier to test, and friendlier to automatic differentiation.

## Section Workflow

`properties(input::OWENSPreComp.Input)` evaluates one blade section. A spanwise
blade model therefore creates one `Input` per station and stores one
`Output` per station:

```julia
tw_prime_d = OWENSPreComp.tw_rate(naf, sloc, tw_aero_d)
input[i] = OWENSPreComp.Input(section_fields...)
output[i] = OWENSPreComp.properties(input[i])
```

The legacy reader workflow used by the tests is:

1. `readmain` reads station data and web-end definitions from the `.pci` file.
2. `tw_rate` differentiates the spanwise twist distribution.
3. `readprecompprofile` reads each airfoil file and returns upper and lower
   surface arrays.
4. `readmaterials` reads and material-id-sorts the material table.
5. `readcompositesection` reads upper, lower, and web laminate schedules for
   the station.
6. The caller combines those arrays into `OWENSPreComp.Input`.
7. `properties` returns the typed `OWENSPreComp.Output`.

## `Input` Field Groups

| Group | Representative fields | Notes |
| --- | --- | --- |
| Span and airfoil | `chord`, `tw_aero_d`, `tw_prime_d`, `le_loc`, `xnode`, `ynode` | Airfoil coordinates describe the normalized perimeter. Twist values ending in `_d` are degrees. |
| Materials | `e1`, `e2`, `g12`, `anu12`, `density` | All material-property vectors must have the same length. Material ids in the layup arrays are one-based indices into these vectors. |
| Upper laminates | `xsec_nodeU`, `n_laminaU`, `n_pliesU`, `t_lamU`, `tht_lamU`, `mat_lamU` | `xsec_nodeU` has one more entry than `n_laminaU`; laminate property vectors are flattened over sectors. |
| Lower laminates | `xsec_nodeL`, `n_laminaL`, `n_pliesL`, `t_lamL`, `tht_lamL`, `mat_lamL` | Same structure as the upper surface. |
| Web laminates | `loc_web`, `n_laminaW`, `n_pliesW`, `t_lamW`, `tht_lamW`, `mat_lamW` | Empty arrays represent no webs at that station. Each web has one laminate stack. |

The validation tests intentionally mutate these groups to confirm that mismatched lengths, invalid material properties, and invalid section geometry throw explicit errors.

## Geometry Contract

Airfoil coordinates are normalized by chord before being passed to
`properties`. The combined `xnode` and `ynode` arrays start at the leading edge
`(0, 0)`, traverse the upper surface to the trailing edge, then return along the
lower surface. The implementation checks that:

- `xnode` and `ynode` have the same length;
- at least three nodes are present;
- the first point is the leading edge at `(0, 0)`;
- the trailing edge does not exceed normalized chord location `1`;
- upper and lower surfaces are single-valued in chordwise `x`;
- sector boundary coordinates are positive, ascending, and within the allowed
  surface bounds.

`readprecompprofile` returns `(xu, yu, xl, yl)`. The tests combine those into
the perimeter form with:

```julia
xnode = vcat(xu, reverse(xl)[2:(end - 1)])
ynode = vcat(yu, reverse(yl)[2:(end - 1)])
```

## Materials And Layups

The material table contains `E1`, `E2`, `G12`, `Nu12`, and density. In the
legacy file reader, rows are sorted by material id before returning property
vectors, so a layup `Mat_id` of `1` indexes the first returned material
property. The kernel checks the orthotropic consistency condition used by
PreComp before computing transformed laminate stiffnesses.

For each upper or lower surface:

- `xsec_node*` gives normalized chordwise sector boundaries;
- `n_lamina*` gives the number of laminate rows in each sector;
- `n_plies*`, `t_lam*`, `tht_lam*`, and `mat_lam*` are flattened over all
  sectors in file order;
- `t_lam*` is ply thickness in meters, and the physical lamina thickness is
  `n_plies * t_lam`;
- `tht_lam*` is the ply material-axis angle in degrees.

For webs, `loc_web` is a normalized chordwise web location at the evaluated
station. The `.pci` file stores web locations at the inboard and outboard web
stations; the test workflow linearly interpolates those locations in physical
space and converts them back to station-local normalized chord coordinates
before calling `readcompositesection`.

## File Reader Returns

| Reader | Return values |
| --- | --- |
| `readmain(path)` | `sloc`, `le_loc`, `chord`, `tw_aero_d`, `af_shape_file`, `int_str_file`, `ib_sp_stn`, `ob_sp_stn`, `ib_ch_loc`, `ob_ch_loc`. `sloc` is dimensionalized by blade length. |
| `readprecompprofile(path)` | `xu`, `yu`, `xl`, `yl` upper and lower normalized airfoil coordinates. |
| `readprofile(path, nheader, LEtoLE)` | Same airfoil-coordinate return as `readprecompprofile`, with selectable header count and coordinate ordering. |
| `readmaterials(path)` | `e1`, `e2`, `g12`, `nu12`, `rho`, `name`, sorted by material id. |
| `readcompositesection(path, locW)` | Upper, lower, and web sector boundaries and flattened laminate arrays. The input `locW` is returned as the web-location vector. |

## `Output` Fields

`OWENSPreComp.Output` stores the computed section properties. The most important quantities for OWENS coupling are:

| Field | Meaning | Unit |
| --- | --- |
| `ei_flap`, `ei_lag` | Flapwise and lag/edgewise bending stiffness | N m^2 |
| `gj` | Torsional stiffness | N m^2 |
| `ea` | Axial stiffness | N |
| `s_fl` | Flap-lag coupled stiffness | N m^2 |
| `s_af`, `s_al` | Axial-flap and axial-lag coupled stiffness | N m |
| `s_ft`, `s_lt` | Flap-torsion and lag-torsion coupled stiffness | N m^2 |
| `s_at` | Axial-torsion coupled stiffness | N m |
| `x_sc`, `y_sc` | Shear-center offset from the reference axes | m |
| `x_tc`, `y_tc` | Tension-center offset from the reference axes | m |
| `mass` | Section mass per unit length | kg/m |
| `flap_iner`, `lag_iner` | Principal flap and lag mass moments per unit length | kg m |
| `tw_iner_d` | Principal inertia-axis orientation | deg |
| `x_cm`, `y_cm` | Center-of-mass offset from the reference axes | m |

When adding tests, pin the actual typed fields needed by the consumer instead of checking only that the output is real-valued.

## Output Shape

`properties(input)` returns only the first 20 section-property values in the
typed `Output`. The lower-level positional `properties(chord, tw_aero_d, ...)`
method returns those 20 values followed by counts used internally and by legacy
callers: airfoil nodes, materials, upper and lower sector counts, web count, and
flattened laminate counts.

The readers and `Input`/`Output` constructors are not currently exported. Use
module-qualified names in user code.
