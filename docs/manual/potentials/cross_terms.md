# Cross terms
\page potential_cross_terms Cross terms

In the simplest force fields the bonded energy is a sum of independent terms:
each bond, bend and torsion has its own potential, as if the others were not
there. In a real molecule the internal coordinates are coupled. *Cross terms*
are potentials that depend on two internal coordinates at once and describe
that coupling. They are what distinguishes the *class-II* force fields (CFF,
PCFF, COMPASS, MM3) from the class-I force fields (AMBER, CHARMM, OPLS,
TraPPE).

## The physics

Expand the energy of a molecule around its equilibrium geometry in its
internal coordinates \f$q_i\f$ (bond lengths, angles, torsions):

\f[
U = U_0 + \frac{1}{2}\sum_i k_{ii}\,(q_i - q_i^0)^2
+ \sum_{i<j} k_{ij}\,(q_i - q_i^0)(q_j - q_j^0) + \dots
\f]

The diagonal terms are the ordinary bond and bend potentials. The
off-diagonal terms \f$k_{ij}\f$ are the cross terms. They have a clear
origin in the electronic structure. When an angle closes, the two outer atoms
come closer and the bonds lengthen to relieve the repulsion: a bond-bend
coupling. When one bond of an angle is stretched, the electrons redistribute
and the other bond shortens: a bond-bond coupling. When a torsion is rotated
towards an eclipsed conformation, the angles open: a bend-torsion coupling.

Cross terms hardly change equilibrium structures. They matter for:

- **vibrational frequencies**, where they determine the splitting between
  modes that involve the same atoms (symmetric and antisymmetric stretches);
  without them frequencies are typically off by tens of wavenumbers;
- **strained molecules**, such as small rings, where the coordinates are far
  from their reference values;
- **flexible frameworks**, where the collective soft modes (breathing,
  negative thermal expansion) involve concerted changes of many bonds and
  angles.

The cross terms of a force field are fitted together with its diagonal terms,
usually to quantum-chemical second derivatives. They cannot be added to a
force field that was fitted without them.

## In the molecule file

Each kind has its own list. An entry gives the atoms, the type and the
parameters; the atoms are given by their *indices* in the list of atoms of
the molecule (counting from zero), not by their types:

```json
"BondBonds" : [
  [[0, 1, 2], "CFF", [5000.0, 1.53, 1.53]]
]
```

Lengths are in Å, angles in degrees, and the force constants in K per unit of
the coordinates involved.

## Bond-bond

`"BondBonds"`: three atoms A-B-C. The coupling between the two bonds of an
angle.

| Type | \f$u\f$ | Parameters |
|---|---|---|
| `"CFF"` | \f$p_0\,(r_{AB} - p_1)(r_{BC} - p_2)\f$ | [K/Å²], [Å], [Å] |
| `"CVFF"` | \f$p_0\,(r_{AB} - p_1)(r_{BC} - p_2)\f$ | [K/Å²], [Å], [Å] |

A positive constant means that stretching one bond makes it favourable to
shorten the other.

## Bond-bend

`"BondBends"`: three atoms A-B-C. The coupling between an angle and its two
bonds.

| Type | \f$u\f$ |
|---|---|
| `"CFF"`, `"CVFF"` | \f$(\theta - p_0)\left[p_1(r_{AB} - p_2) + p_3(r_{BC} - p_4)\right]\f$ |
| `"MM3"` | \f$p_0\left[(r_{AB} - p_1) + (r_{BC} - p_2)\right](\theta - p_3)\f$ |
| `"TRUNCATED_HARMONIC"` | \f$\frac{1}{2}p_0(\theta - p_1)^2\,\exp\!\left[-(r_{AB}^8 + r_{BC}^8)/p_2^8\right]\f$ |
| `"SCREENED_HARMONIC"` | \f$\frac{1}{2}p_0(\theta - p_1)^2\,\exp\!\left[-(r_{AB}/p_2 + r_{BC}/p_3)\right]\f$ |
| `"SCREENED_VESSAL"` | \f$\dfrac{p_0}{8(\theta - \pi)^2}\left[(p_1 - \pi)^2 - (\theta - \pi)^2\right]^2 \exp\!\left[-(r_{AB}/p_2 + r_{BC}/p_3)\right]\f$ |
| `"TRUNCATED_VESSAL"` | a Vessal bend with a power-law prefactor in the angle |

The first two are the stretch-bend terms of the class-II and MM3 force
fields. The *truncated* and *screened* forms are of a different kind: a bend
potential whose strength decays when the bonds are stretched, so that the
angle term switches off as a bond breaks. They come from force fields for
inorganic solids (the forms of DL_POLY), where bonds are not permanent.

## Bend-bend

`"BendBends"`: four atoms, two angles that share a central atom and a bond.

| Type | \f$u\f$ |
|---|---|
| `"CFF"`, `"CVFF"` | \f$p_0\,(\theta_1 - p_1)(\theta_2 - p_2)\f$ |
| `"MM3"` | \f$-p_0\,(\theta_1 - p_1)(\theta_2 - p_2)\f$ |

It describes that the angles around one atom are not independent: opening one
closes its neighbours.

## Bond-torsion

`"BondTorsions"`: four atoms A-B-C-D. The coupling between a torsion and the
length of its central bond.

| Type | \f$u\f$ |
|---|---|
| `"MM3"` | \f$\frac{1}{2}(r - p_3)\left[p_0(1 + \cos\phi) + p_1(1 + \cos 2\phi) + p_2(1 + \cos 3\phi)\right]\f$ |

The central bond is longer in the eclipsed conformation than in the staggered
one; this term reproduces that.

## Bend-torsion

`"BendTorsions"`: four atoms A-B-C-D. The coupling between a torsion and its
two angles \f$\theta_1\f$ (A-B-C) and \f$\theta_2\f$ (B-C-D).

| Type | \f$u\f$ |
|---|---|
| `"CFF"`, `"CVFF"` | \f$p_0\,(\theta_1 - p_1)(\theta_2 - p_2)\cos\phi\f$ |
| `"SMOOTHED_CFF3"` | the same, optionally smoothed |
| `"SMOOTHED"` | \f$S(\theta_1)\,p_0\left[1 + \cos(p_1\phi - p_2)\right]S(\theta_2)\f$ |
| `"SMOOTHED_THREE_COSINE"`, `"NICHOLAS"` | \f$S(\theta_1)\left[\frac{1}{2}p_0(1 + \cos\phi) + \frac{1}{2}p_1(1 - \cos 2\phi) + \frac{1}{2}p_2(1 + \cos 3\phi)\right]S(\theta_2)\f$ |
| `"SMOOTHED_CFF"` | \f$S(\theta_1)\left[p_0(1 - \cos\phi) + p_1(1 - \cos 2\phi) + p_2(1 - \cos 3\phi)\right]S(\theta_2)\f$ |
| `"SMOOTHED_CFF2"` | a smoothed three-term cosine series |

**The smoothing function.** A torsion angle is undefined when three of its
four atoms are collinear: as an angle \f$\theta\f$ approaches 180°, a small
displacement changes \f$\phi\f$ by a large amount, and the forces of a plain
torsion potential diverge. This happens in frameworks with nearly linear
linkages, such as the Si-O-Si angles of zeolites. The smoothed types multiply
the torsion by a function \f$S(\theta)\f$ of each angle that is one below
170° and goes smoothly to zero at 180°,

\f[
S(\theta) = \frac{(\theta_\text{off} - \theta)^2\,(\theta_\text{off} + 2\theta - 3\theta_\text{on})}
{(\theta_\text{off} - \theta_\text{on})^3}
\quad\text{for } \theta_\text{on} < \theta < \theta_\text{off},
\f]

with \f$\theta_\text{on} = 170^\circ\f$ and \f$\theta_\text{off} = 180^\circ\f$,
so that the energy and the force stay finite. These are therefore torsion
potentials made safe for flexible frameworks rather than couplings in the
sense of the expansion above.

## Relation to Urey-Bradley

A [Urey-Bradley](\ref potential_urey_bradley) term, a distance potential
between the outer atoms of an angle, also couples the two bonds and the bend.
It is the simpler, older way of including the bond-bond and bond-bend
couplings, with a single force constant instead of separate ones.

## References

- J.R. Maple, M.-J. Hwang, T.P. Stockfisch, U. Dinur, M. Waldman, C.S. Ewig
  and A.T. Hagler, *Derivation of class II force fields. I. Methodology and
  quantum force field for the alkyl functional group and alkane molecules*,
  J. Comput. Chem. **15**, 162 (1994).
- N.L. Allinger, Y.H. Yuh and J.-H. Lii, *Molecular mechanics. The MM3 force
  field for hydrocarbons. 1*, J. Am. Chem. Soc. **111**, 8551 (1989).
- H. Sun, *COMPASS: An ab Initio Force-Field Optimized for Condensed-Phase
  Applications: Overview with Details on Alkane and Benzene Compounds*,
  J. Phys. Chem. B **102**, 7338 (1998).
- J.B. Nicholas, A.J. Hopfinger, F.R. Trouw and L.E. Iton, *Molecular modeling
  of zeolite structure. 2. Structure and dynamics of silica sodalite and
  silicate force field*, J. Am. Chem. Soc. **113**, 4792 (1991).
- E.B. Wilson Jr., J.C. Decius and P.C. Cross, *Molecular Vibrations*,
  McGraw-Hill (1955).
