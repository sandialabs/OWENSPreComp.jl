# Theory, Frames, and Units

OWENSPreComp follows the PreComp section-property approach: a modified classical laminate theory is combined with a shear-flow section solve to compute equivalent beam properties for composite blades.

## Modeling Assumptions

- The blade section is represented by airfoil perimeter segments, laminate sectors, and optional webs.
- Laminates are modeled with ply material properties, ply thickness, and fiber orientation.
- The output is a set of equivalent section properties suitable for beam models, not a three-dimensional stress field.
- Shear-center calculations follow the PreComp approximation and should be validated for unusual laminate layouts.

## Frames

| Frame or point | Meaning |
| --- | --- |
| reference axis | User-specified section reference, often related to the pitch axis. |
| airfoil coordinates | Section perimeter coordinates used for upper and lower surface geometry. |
| shear center | Elastic reference for bending and torsion properties. |
| tension center | Reference point for axial loading. |
| mass center | Inertial reference point for mass and rotary inertia. |
| principal inertia axes | Axes used for principal section inertias and their orientation. |

Offsets returned by `properties` should always be passed downstream with their associated frame. Mixing airfoil, reference, shear-center, and mass-center frames is a common source of coupled-model errors.

## Units

The Julia interface uses the same conventions as the PreComp input files:

| Quantity | Unit |
| --- | --- |
| chord, offsets, coordinates, ply thickness | meters |
| material stiffness and shear moduli | pascal |
| material density | kg/m^3 |
| aerodynamic twist fields ending in `_d` | degrees |
| twist rate | degrees per meter |
| stiffness outputs | SI beam-property units, for example N, N m^2, or N m^2/rad depending on term |
| inertia outputs | kg/m and kg m as appropriate for section properties |

Tests for new consumers should include a units note in fixture comments and pin both magnitudes and signs for offsets.
