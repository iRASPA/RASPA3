# Potentials
\page potentials Potentials

The energy of a system is a sum of terms, each a function of the positions of
a few atoms. This page gives an overview of the kinds of terms `RASPA3` knows;
each links to a page with the physical background, the functional forms that
are available, their parameters and units, and how they are written in the
input files.

\f[
U = \underbrace{U_\text{dispersion} + U_\text{electrostatic} + U_\text{polarization}}_{\text{non-bonded}}
+ \underbrace{U_\text{bond} + U_\text{bend} + U_\text{torsion} + U_\text{cross} + \dots}_{\text{bonded}}
+ U_\text{external}
\f]

## Where the parameters are given

| File | Contains |
|---|---|
| `force_field.json` | the pseudo-atoms (mass, charge), the van der Waals interactions between them, the mixing rule, the cutoffs and the treatment of the tail |
| the molecule file of a component | the atoms of the molecule, its connectivity, and the bonded terms |
| `simulation.json` | the charge method of a system, and whether polarization and an external field are used |

An atom of a molecule or framework refers to a *pseudo-atom* of the force
field by name. A pseudo-atom is an interaction site: an atom, or a group of
atoms treated as one (a united-atom `CH3`, for instance).

## Non-bonded interactions

These act between atoms of different molecules, between molecules and the
framework, and between atoms of one flexible molecule that are far enough
apart along its bonds.

| Interaction | What it describes | Page |
|---|---|---|
| Dispersion and repulsion | The attraction between fluctuating electron clouds and the repulsion of overlapping ones: Lennard-Jones and twenty other pair potentials, with the mixing rules, the cutoff and the tail correction | \subpage potential_dispersion |
| Electrostatics | The interaction between partial charges, summed over the periodic images with the Ewald method or a cutoff-based alternative | \subpage potential_electrostatics |
| Polarization | The energy of dipoles induced on the atoms by the electric field of the charges around them | \subpage potential_polarization |
| Non-bonded interactions within a molecule | Which pairs of atoms of one flexible molecule interact, which are excluded, and how the 1-4 pairs are scaled | \subpage potential_intramolecular_nonbonded |

## Bonded interactions

These act between atoms of one flexible molecule (or flexible framework) that
are connected by bonds. Each term has a list in the molecule file. Bonds,
bends and torsions are given per combination of pseudo-atom *types* and
apply to every matching set of bonded atoms; the other terms are listed
explicitly, with the *indices* of the atoms they act on.

| Interaction | Atoms | What it describes | Page |
|---|---|---|---|
| Bonds | 2 | The stretching of a chemical bond | \subpage potential_bonds |
| Urey-Bradley | 2 (1-3) | A distance term between the two outer atoms of an angle | \subpage potential_urey_bradley |
| Bends | 3 | The bending of the angle between two bonds | \subpage potential_bends |
| Inversion bends | 4 | The out-of-plane motion of an atom bonded to three others | \subpage potential_inversion_bends |
| Torsions | 4 | The rotation about a bond, including improper torsions | \subpage potential_torsions |
| Cross terms | 3-4 | Couplings between bonds, bends and torsions (class-II force fields) | \subpage potential_cross_terms |
| Correction maps (CMAP) | 5 | A tabulated energy of two consecutive torsion angles | \subpage potential_cmap |

## External fields

| Interaction | What it describes | Page |
|---|---|---|
| External field | A potential that depends on the position of an atom only: a confining channel, or a tabulated field | \subpage potential_external_fields |

## Rigid, flexible and constrained

How many of these terms a molecule has depends on how it is modelled.

- A **rigid** molecule has no bonded terms and no interactions within itself.
  Its geometry is that of the molecule file, and it moves by translation and
  rotation.
- A **flexible** molecule has bonded terms, and non-bonded interactions
  between its atoms that are separated by more than three bonds.
- A **semi-flexible** molecule consists of rigid bodies (a phenyl ring, for
  instance) connected by flexible parts. Atoms within one rigid body do not
  interact with each other.
- A bond, a Urey-Bradley distance, a bend or a torsion can be *fixed*: the
  type `"FIXED"` is not a potential but a constraint, which holds the
  coordinate at its value. In Monte Carlo the molecule is grown and moved
  with that coordinate fixed; in molecular dynamics it is held with
  SHAKE and RATTLE.

## Units

Parameters are given in the units of the input: energies as \f$E/k_B\f$ in
kelvin, lengths in ångström, angles in degrees, charges in units of the
elementary charge. Force constants follow from these: K/Å² for a harmonic
bond, K/rad² for a harmonic bend. A few types keep the units of the force
field they come from (the MM3 terms use mdyne and kcal/mol); the pages say so
per type.

## Scaling with the coupling parameter

A *fractional* molecule, used in the continuous-fractional-component moves
and in thermodynamic integration, interacts with a strength set by its
coupling parameter \f$\lambda\f$. Its van der Waals and electrostatic
interactions with the rest of the system are scaled; the bonded terms of the
molecule itself are not. The pages on
[dispersion](\ref potential_dispersion) and
[electrostatics](\ref potential_electrostatics) describe the scaling.
