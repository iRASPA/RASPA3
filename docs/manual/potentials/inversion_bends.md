# Inversion bends
\page potential_inversion_bends Inversion bends

An inversion bend describes the out-of-plane motion of an atom that is bonded
to three others. It is the term that keeps a planar centre planar (the carbon
of a carbonyl group, an aromatic carbon with its substituent) and that sets
the height of the pyramid of a pyramidal centre (the nitrogen of an amine).

## The coordinate

Take a central atom B bonded to three atoms A, C and D. The inversion angle
\f$\chi\f$ is the angle between the bond B-A and a plane through three of
the atoms:

\f[
\sin\chi = \frac{\mathbf{r}_{BA}\cdot(\mathbf{r}_{BC}\times\mathbf{r}_{BD})}
{r_{BA}\,\left|\mathbf{r}_{BC}\times\mathbf{r}_{BD}\right|} .
\f]

For a planar centre \f$\chi = 0\f$. When the central atom moves out of the
plane of its neighbours, \f$\chi\f$ grows. An *inversion* is the passage of
the central atom through the plane to the other side, as in the umbrella
motion of ammonia.

Two choices of plane are available, which is the difference between the
types with and without the suffix `2`:

- the plane through the central atom and the two other neighbours, B-C-D;
- the plane through the three neighbours, A-C-D.

For small displacements the two angles are proportional, so the choice
follows the definition used by the force field the parameters come from.

## Why the bends are not enough

A centre with three neighbours has three bend angles. For a planar centre
they add up to 360°, and for a small out-of-plane displacement \f$h\f$ their
sum decreases only as \f$h^2\f$: each individual angle changes to second
order. The harmonic bends therefore give almost no restoring force against
out-of-plane motion, and the potential they produce for it is quartic rather
than quadratic. A separate term in the out-of-plane coordinate is needed to
give planar groups their stiffness and the out-of-plane vibrations their
frequencies.

## In the molecule file

Each entry of `"InversionBends"` gives the four atoms A, B, C, D, with the
central atom second, the type of potential and its parameters. The atoms are
given by their *indices* in the list of atoms of the molecule (counting from
zero), not by their types:

```json
"InversionBends" : [
  [[3, 1, 0, 2], "HARMONIC", [20000.0, 0.0]]
]
```

For a centre with three different neighbours there are three such angles,
one for each choice of the atom A; force fields either apply the term to all
three and average, or to one. The entries of the list determine which.

## Available forms

| Type | Plane | \f$u(\chi)\f$ | Parameters |
|---|---|---|---|
| `"HARMONIC"` | B-C-D | \f$\frac{1}{2}p_0\,(\chi - p_1)^2\f$ | \f$k\f$, \f$\chi_0\f$ [°] |
| `"HARMONIC_COSINE"` | B-C-D | \f$\frac{1}{2}p_0\,(\cos\chi - \cos p_1)^2\f$ | \f$k\f$ [K], \f$\chi_0\f$ [°] |
| `"PLANAR"` | B-C-D | \f$p_0\,(1 - \cos\chi)\f$ | \f$k\f$ [K] |
| `"HARMONIC2"` | A-C-D | \f$\frac{1}{2}p_0\,(\chi - p_1)^2\f$ | \f$k\f$, \f$\chi_0\f$ [°] |
| `"HARMONIC_COSINE2"` | A-C-D | \f$\frac{1}{2}p_0\,(\cos\chi - \cos p_1)^2\f$ | \f$k\f$ [K], \f$\chi_0\f$ [°] |
| `"PLANAR2"` | A-C-D | \f$p_0\,(1 - \cos\chi)\f$ | \f$k\f$ [K] |
| `"MM3"` | A-C-D | \f$p_0\,\Delta^2\left[1 - 0.014\Delta + 5.6\times10^{-5}\Delta^2 - 7\times10^{-7}\Delta^3 + 2.2\times10^{-8}\Delta^4\right]\f$, \f$\Delta = \chi - p_1\f$ in degrees | [mdyne Å/rad²], \f$\chi_0\f$ [°] |

Notes on the forms:

- **Harmonic** with \f$\chi_0 = 0\f$ holds a centre planar. With
  \f$\chi_0 \ne 0\f$ the minimum is pyramidal, but the harmonic form has one
  minimum only, on one side of the plane.
- **Harmonic cosine** is even in \f$\chi\f$: it has two equivalent minima at
  \f$\pm\chi_0\f$ with a barrier at the planar geometry of height
  \f$\frac{1}{2}k(1 - \cos\chi_0)^2\f$. This is the form for centres that can
  invert, such as amines.
- **Planar** is the form of UFF and DREIDING for planar centres; it is
  harmonic for small \f$\chi\f$ with force constant \f$p_0\f$.

## Inversion bend or improper torsion

Force fields use one of two coordinates for the same physics. The inversion
angle of this page is the choice of DREIDING, UFF, MM3 and the class-II force
fields. AMBER, CHARMM and OPLS use an *improper torsion* instead: the angle
between the planes A-B-C and B-C-D of four atoms that are not bonded in a
chain. In `RASPA3` improper torsions are entered in their own list with the
forms of a [torsion](\ref potential_torsions). The two are not
interchangeable: parameters fitted with one coordinate cannot be used with
the other.

## Chirality

An inversion term with one minimum fixes the side of the plane on which the
central atom sits, and so preserves the handedness of a chiral centre in a
united-atom model, where the fourth substituent (a hydrogen) is not present
explicitly. Without it such a centre would invert freely.

## References

- S.L. Mayo, B.D. Olafson and W.A. Goddard III, *DREIDING: a generic force
  field for molecular simulations*, J. Phys. Chem. **94**, 8897 (1990).
- A.K. Rappé, C.J. Casewit, K.S. Colwell, W.A. Goddard III and W.M. Skiff,
  *UFF, a full periodic table force field for molecular mechanics and
  molecular dynamics simulations*, J. Am. Chem. Soc. **114**, 10024 (1992).
- N.L. Allinger, Y.H. Yuh and J.-H. Lii, *Molecular mechanics. The MM3 force
  field for hydrocarbons. 1*, J. Am. Chem. Soc. **111**, 8551 (1989).
- E.B. Wilson Jr., J.C. Decius and P.C. Cross, *Molecular Vibrations*,
  McGraw-Hill (1955).
