# Torsions
\page potential_torsions Torsions

A torsion term describes the energy of rotating one part of a molecule
relative to another about a bond. It is the softest of the bonded terms, and
the one that determines the *conformation* of a flexible molecule: whether a
chain is extended or folded.

## The coordinate

For four atoms bonded in a chain, A-B-C-D, the torsion (dihedral) angle
\f$\phi\f$ is the angle between the plane through A, B, C and the plane
through B, C, D: the rotation about the central bond B-C that carries the
bond B-A onto the direction of C-D, seen along B-C. With the normals of the
two planes,

\f[
\mathbf{m} = \mathbf{r}_{AB}\times\mathbf{r}_{BC}, \qquad
\mathbf{n} = \mathbf{r}_{BC}\times\mathbf{r}_{CD}, \qquad
\cos\phi = \frac{\mathbf{m}\cdot\mathbf{n}}{|\mathbf{m}|\,|\mathbf{n}|} ,
\f]

and the sign of \f$\phi\f$ follows from the sense of the rotation. `RASPA3`
uses the IUPAC convention: the *trans* (anti) conformation, with A and D on
opposite sides, is \f$\phi = 180^\circ\f$, and the *cis* conformation is
\f$\phi = 0\f$.

*A note on conventions.* The polymer literature often uses the opposite
choice, trans at \f$\phi = 0\f$. A potential written for one convention is
converted to the other by \f$\phi \to \phi - 180^\circ\f$, which changes the
sign of the odd cosine terms. Parameters taken from the literature must be
checked for this.

## The physics

The energy of a rotation about a single bond has its origin in the
interaction between the bonds and substituents at the two ends: the
repulsion between electron pairs when they eclipse each other, and
hyperconjugation. For butane the result is a profile with three minima: the
*trans* conformation, lowest in energy, and two equivalent *gauche*
conformations near \f$\pm 60^\circ\f$, about 0.9 kcal/mol higher, separated
by barriers of 3-5 kcal/mol. These numbers are of the order of the thermal
energy, so all three are populated at room temperature and the molecule
changes between them on a time scale of tens of picoseconds.

Because the profile is periodic in \f$\phi\f$, it is naturally written as a
Fourier series in the angle. The symmetry of the groups sets which terms
appear: a rotation between two \f$sp^3\f$ centres has a threefold term, a
rotation about a bond with double-bond character a twofold term.

## In the molecule file

Each entry of `"Torsions"` names the four pseudo-atom types in the order of
the chain, the type of potential and its parameters:

```json
"Torsions" : [
  [["CH3", "CH2", "CH2", "CH3"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
]
```

Energies are in K and angles in degrees.

## Available forms

| Type | \f$u(\phi)\f$ |
|---|---|
| `"TRAPPE"` | \f$p_0 + p_1(1 + \cos\phi) + p_2(1 - \cos 2\phi) + p_3(1 + \cos 3\phi)\f$ |
| `"TRAPPE_EXTENDED"` | \f$p_0 + p_1\cos\phi + p_2\cos 2\phi + p_3\cos 3\phi + p_4\cos 4\phi\f$ |
| `"MOD_TRAPPE"` | the TraPPE form with the angle shifted: \f$\phi \to \phi - p_4\f$ |
| `"OPLS"` | \f$\frac{1}{2}p_0 + \frac{1}{2}p_1(1 + \cos\phi) + \frac{1}{2}p_2(1 - \cos 2\phi) + \frac{1}{2}p_3(1 + \cos 3\phi)\f$ |
| `"THREE_COSINE"` | \f$\frac{1}{2}p_0(1 + \cos\phi) + \frac{1}{2}p_1(1 - \cos 2\phi) + \frac{1}{2}p_2(1 + \cos 3\phi)\f$ |
| `"MM3"` | as `"THREE_COSINE"`, with the parameters in kcal/mol |
| `"FOURIER_SERIES"` | \f$\frac{1}{2}\left[p_0(1 + \cos\phi) + p_1(1 - \cos 2\phi) + p_2(1 + \cos 3\phi) + p_3(1 - \cos 4\phi) + p_4(1 + \cos 5\phi) + p_5(1 + \cos 6\phi)\right]\f$ |
| `"FOURIER_SERIES2"` | as `"FOURIER_SERIES"`, with \f$p_3(1 + \cos 4\phi)\f$ |
| `"CFF"` | \f$p_0(1 - \cos\phi) + p_1(1 - \cos 2\phi) + p_2(1 - \cos 3\phi)\f$ |
| `"CFF2"` | \f$p_0(1 + \cos\phi) + p_1(1 + \cos 2\phi) + p_2(1 + \cos 3\phi)\f$ |
| `"CVFF"` | \f$p_0\left[1 + \cos(p_1\phi - p_2)\right]\f$ |
| `"RYCKAERT_BELLEMANS"` | \f$\sum_{i=0}^{5} p_i\cos^i\phi\f$ |
| `"POLYNOMIAL"` (`"NHARMONIC"`) | \f$\sum_{i=0}^{5} p_i\cos^i\phi\f$ |
| `"SIX_COSINE"` | a six-term series in powers of \f$\cos\phi\f$ |
| `"HARMONIC"` | \f$\frac{1}{2}p_0\,(\phi - p_1)^2\f$ |
| `"HARMONIC_COSINE"` | \f$\frac{1}{2}p_0\,(\cos\phi - \cos p_1)^2\f$ |
| `"CVFF_BLOCKED"` | no energy; used to detect a blocked region of the angle |
| `"Fixed"` | a constraint: the angle is held at its value |

The parameters are energies in K, except the multiplicity \f$p_1\f$ and the
phase \f$p_2\f$ [°] of `"CVFF"`, the reference angle [°] of the harmonic
forms, with a force constant in K/rad² for `"HARMONIC"`, and the shift
\f$p_4\f$ [°] of `"MOD_TRAPPE"`.

Notes on the forms:

- **Cosine series** (TraPPE, OPLS, three-cosine, Fourier). These are the same
  kind of function, a sum of cosines of multiples of the angle, with different
  conventions for the constant and for the factor \f$\frac{1}{2}\f$. Each term
  has a meaning: the \f$\cos\phi\f$ term distinguishes cis from trans, the
  \f$\cos 2\phi\f$ term favours planar or perpendicular arrangements, and the
  \f$\cos 3\phi\f$ term is the threefold barrier of \f$sp^3\f$ centres.
- **CVFF.** One cosine with a multiplicity \f$n\f$ and a phase, the form of
  AMBER and CHARMM. A profile that needs several multiplicities is entered as
  several terms on the same four atoms.
- **Powers of the cosine** (Ryckaert-Bellemans, polynomial). A polynomial in
  \f$\cos\phi\f$ is equivalent to a cosine series, since \f$\cos n\phi\f$ is a
  polynomial of degree \f$n\f$ in \f$\cos\phi\f$. It is cheap to evaluate,
  because \f$\cos\phi\f$ follows directly from the positions, and its
  derivative has no singularity at 0° or 180°. The polynomial form is the
  `nharmonic` dihedral of LAMMPS. The Ryckaert-Bellemans potential was
  published in the polymer convention (trans at 0°); check the sign of the odd
  coefficients against the convention used here.
- **Harmonic.** For a torsion that is held near one value, as an
  [improper torsion](#potential-torsions-improper). The plain harmonic form is
  not periodic and has a discontinuous derivative where the angle wraps; the
  harmonic-cosine form avoids that.

## Improper torsions {#potential-torsions-improper}

An *improper* torsion uses the same angle between two planes, for four atoms
that are **not** bonded in a chain: a central atom and its three neighbours.
It serves the purpose of an
[inversion bend](\ref potential_inversion_bends), keeping a planar centre
planar or a chiral centre in its configuration, in the force fields that use
this coordinate (AMBER, CHARMM, OPLS). The entries are given in the list
`"ImproperTorsions"`, with the forms of the table above. Because the four
atoms are not a chain, they are given by their *indices* in the list of atoms
of the molecule (counting from zero), not by their types:

```json
"ImproperTorsions" : [
  [[1, 0, 2, 3], "CVFF", [5280.0, 2.0, 180.0]]
]
```

The order of the four atoms matters and differs between force fields; it
must be taken over together with the parameters.

## The 1-4 interaction

The two outer atoms of a torsion are a 1-4 pair. Their non-bonded interaction
varies with the angle and so contributes to the rotational profile. Whether
it is included, and with what scaling, is part of the force field; see
[non-bonded interactions within a molecule](\ref potential_intramolecular_nonbonded).
Torsion parameters are only valid with the 1-4 scaling they were fitted
with.

## Sampling torsions

Torsions are the coordinates that configurational-bias Monte Carlo was
designed for. When a chain is grown atom by atom, the bond length and the
bend angle of a new atom are drawn from their own distributions, which leaves
the torsion angle: the position of the atom on the circle around the previous
bond. Several trial positions are generated on that circle with a probability
proportional to \f$e^{-\beta u(\phi)}\f$, and one is chosen according to its
non-bonded energy. The growth thereby favours conformations that are both low
in torsional energy and fit in the surroundings, which is what makes the
insertion of chain molecules in a dense phase or a narrow pore possible.

In molecular dynamics the torsional barriers set the time scale on which
conformations change. Barriers of several \f$k_BT\f$ make this slow, and for
chains it is the reason to use replica exchange.

## Fixed torsions

A torsion of type `"Fixed"` is held at its value. A molecule with all bonds,
bends and torsions fixed is a rigid body, and is better declared as one.

## References

- J.-P. Ryckaert and A. Bellemans, *Molecular dynamics of liquid alkanes*,
  Faraday Discuss. Chem. Soc. **66**, 95 (1978).
- M.G. Martin and J.I. Siepmann, *Transferable Potentials for Phase
  Equilibria. 1. United-Atom Description of n-Alkanes*,
  J. Phys. Chem. B **102**, 2569 (1998).
- W.L. Jorgensen, D.S. Maxwell and J. Tirado-Rives, *Development and Testing
  of the OPLS All-Atom Force Field on Conformational Energetics and
  Properties of Organic Liquids*, J. Am. Chem. Soc. **118**, 11225 (1996).
- P. Dauber-Osguthorpe, V.A. Roberts, D.J. Osguthorpe, J. Wolff, M. Genest
  and A.T. Hagler, *Structure and energetics of ligand binding to proteins:
  Escherichia coli dihydrofolate reductase-trimethoprim, a drug-receptor
  system*, Proteins **4**, 31 (1988).
- J.I. Siepmann and D. Frenkel, *Configurational bias Monte Carlo: a new
  sampling scheme for flexible chains*, Mol. Phys. **75**, 59 (1992).
- IUPAC-IUB Commission on Biochemical Nomenclature, *Abbreviations and
  symbols for the description of the conformation of polypeptide chains*,
  Biochemistry **9**, 3471 (1970).
