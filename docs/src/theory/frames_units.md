# Theory, Frames, and Units

OWENSPreComp follows the PreComp section-property approach: a modified classical laminate theory is combined with a shear-flow section solve to compute equivalent beam properties for composite blades.

## Modeling Assumptions

- Each blade station is represented as a thin-walled closed multicell section
  made from airfoil perimeter segments and optional straight webs.
- Upper and lower surfaces are divided into chordwise laminate sectors. Each
  sector is a stack of laminae, and each lamina is represented by number of
  plies, ply thickness, material id, and fiber angle.
- Webs are normal to the chord at a station and each web has one laminate stack.
- The section is assumed not to distort in its own plane, transverse shear is
  neglected, and warping is treated with the PreComp free-warping assumption.
- The output is a set of equivalent beam properties, not a three-dimensional
  stress or strain field.
- Shear-center calculations follow the PreComp approximation and should be
  validated for unusual anisotropic layouts, open-section-like geometry, or
  root-region sections where constrained warping may matter.

## Frames

| Frame or point | Meaning |
| --- | --- |
| airfoil coordinates | Normalized section coordinates beginning at the leading edge. The chordwise coordinate runs from leading edge to trailing edge. |
| reference axes, `XR-YR` | Section reference axes used for returned offsets. In the legacy files, `le_loc` locates the leading edge relative to this reference and is normalized by chord. |
| shear center | Elastic reference for bending and torsion properties. Returned as `x_sc`, `y_sc` relative to the reference axes. |
| tension center | Reference point for axial loading. Returned as `x_tc`, `y_tc` relative to the reference axes. |
| center of mass | Inertial reference point for mass and rotary inertia. Returned as `x_cm`, `y_cm` relative to the reference axes. |
| principal inertia axes | Axes used for principal section inertias. `tw_iner_d` is their orientation in degrees after the PreComp wind-turbine sign conversion. |

Offsets returned by `properties` should always be passed downstream with their associated frame. Mixing airfoil, reference, shear-center, and mass-center frames is a common source of coupled-model errors.

## Sign Conventions

Input airfoil coordinates and sector locations are normalized by chord, while
returned offsets are dimensional. In the wind-turbine output form hardwired in
the current implementation, selected `y`-direction offsets and coupling terms
are sign-converted to match the legacy PreComp output convention. The pinned
fixture compares:

- `x_sc`, `y_sc`, `x_tc`, `y_tc`, `x_cm`, `y_cm`;
- `s_fl`, `s_af`, `s_al`, `s_ft`, `s_lt`, `s_at`;
- `tw_iner_d`.

Downstream coupling code should pin signs for offsets and coupled stiffnesses,
not only their magnitudes.

## Units

The Julia interface uses the same conventions as the PreComp input files:

| Quantity | Unit |
| --- | --- |
| chord, offsets, coordinates, ply thickness | meters |
| material stiffness and shear moduli | pascal |
| material density | kg/m^3 |
| aerodynamic twist fields ending in `_d` | degrees |
| twist rate | degrees per meter |
| bending and torsion stiffness outputs | N m^2 |
| axial stiffness output | N |
| axial-bending and axial-torsion coupling outputs | N m |
| bending-torsion and flap-lag coupling outputs | N m^2 |
| inertia outputs | kg/m and kg m as appropriate for section properties |

Tests for new consumers should include a units note in fixture comments and pin both magnitudes and signs for offsets.

## Twist Rate

`tw_rate(naf, sloc, tw_aero_d)` converts aerodynamic twist from degrees to
radians, computes finite-difference derivatives with respect to dimensional
span location, and returns the twist rate in degrees per meter. Interior points
use a nonuniform-grid centered derivative; the first and last stations use
one-sided differences. The derivative is a required `Input` field because the
section-property solve includes twist-gradient terms.
