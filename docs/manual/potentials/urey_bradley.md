# Urey-Bradley
\page potential_urey_bradley Urey-Bradley

A Urey-Bradley term is a distance potential between the two *outer* atoms of
an angle A-B-C: atoms A and C, which are two bonds apart (a 1-3 pair).

\f[
U_\text{UB} = \sum_{\text{1-3 pairs}} u(r_{AC}) , \qquad
u(r) = \frac{1}{2}\,k\,(r - r_0)^2 .
\f]

## The physics

The 1-3 distance is determined by the two bond lengths and the angle between
them,

\f[
r_{AC}^2 = r_{AB}^2 + r_{BC}^2 - 2\,r_{AB}\,r_{BC}\cos\theta ,
\f]

so a potential in \f$r_{AC}\f$ acts on all three at once. It describes the
repulsion between the two outer atoms, which is what couples the bending of
an angle to the stretching of its bonds: when the angle closes, the outer
atoms approach, and the bonds lengthen to relieve it. A separate
[bend](\ref potential_bends) and [bond](\ref potential_bonds) term cannot
produce this coupling. The Urey-Bradley term is the older and simpler way of
including it; the [cross terms](\ref potential_cross_terms) of class-II force
fields are the more systematic one.

In practice the term is used *together with* a bend term on the same angle,
and improves the vibrational frequencies, in particular the splitting between
symmetric and antisymmetric stretching modes. It is part of the CHARMM force
fields.

## In the molecule file

Each entry of `"UreyBradleys"` gives the two outer atoms, the type of
potential and its parameters. The atoms are given by their *indices* in the
list of atoms of the molecule (counting from zero), not by their types:

```json
"UreyBradleys" : [
  [[0, 2], "HARMONIC", [15000.0, 1.80]]
]
```

## Available forms

The term is a function of one distance, so it takes the same forms as a
[bond](\ref potential_bonds), with the same parameters and units:

| Type | \f$u(r)\f$ |
|---|---|
| `"HARMONIC"` | \f$\frac{1}{2}p_0\,(r - p_1)^2\f$ |
| `"MORSE"` | \f$p_0\left\{\left[1 - e^{-p_1(r - p_2)}\right]^2 - 1\right\}\f$ |
| `"QUARTIC"` | \f$\frac{1}{2}p_0(r-p_1)^2 + \frac{1}{3}p_2(r-p_1)^3 + \frac{1}{4}p_3(r-p_1)^4\f$ |
| `"CFF_QUARTIC"` | \f$p_0(r-p_1)^2 + p_2(r-p_1)^3 + p_3(r-p_1)^4\f$ |
| `"MM3"` | \f$p_0(r-p_1)^2\left[1 - 2.55(r-p_1) + \frac{7}{12}\,2.55^2\,(r-p_1)^2\right]\f$ |
| `"RESTRAINED_HARMONIC"` | harmonic within \f$p_2\f$ of \f$p_1\f$, linear beyond |
| `"CORE_SHELL_SPRING"` | \f$\frac{1}{2}p_0\,r^2\f$ |
| `"LENNARD_JONES"` | \f$4p_0\left[(p_1/r)^{12} - (p_1/r)^{6}\right]\f$ |
| `"LJ_12_6"` | \f$p_0/r^{12} - p_1/r^{6}\f$ |
| `"BUCKINGHAM"` | \f$p_0\,e^{-p_1 r} - p_2/r^{6}\f$ |
| `"FIXED"` | a constraint on the 1-3 distance |

## Fixed 1-3 distances

A `"FIXED"` Urey-Bradley distance, together with two fixed bonds, fixes the
angle between them. This is how a rigid angle is imposed with distance
constraints only, as in rigid water models held with SHAKE.

## Relation to the non-bonded interaction

A 1-3 pair has no van der Waals or Coulomb interaction (see
[non-bonded interactions within a molecule](\ref potential_intramolecular_nonbonded)).
The Urey-Bradley term is the explicit, fitted replacement for it.

## References

- H.C. Urey and C.A. Bradley, *The Vibrations of Pentatonic Tetrahedral
  Molecules*, Phys. Rev. **38**, 1969 (1931).
- A.D. MacKerell Jr. *et al.*, *All-Atom Empirical Potential for Molecular
  Modeling and Dynamics Studies of Proteins*,
  J. Phys. Chem. B **102**, 3586 (1998).
