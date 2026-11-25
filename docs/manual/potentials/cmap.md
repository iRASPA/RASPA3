# Correction maps (CMAP)
\page potential_cmap Correction maps (CMAP)

A correction map is a *tabulated* energy as a function of two consecutive
torsion angles. It is used in protein force fields (CHARMM, and AMBER
ff19SB), where it corrects the energy surface of the two backbone torsions of
each amino-acid residue.

## The coordinate

Five atoms bonded in a chain, A-B-C-D-E, define two torsion angles that share
the bond in the middle of the chain:

- \f$\phi\f$ over the atoms A-B-C-D,
- \f$\psi\f$ over the atoms B-C-D-E.

Both follow the IUPAC sign convention (trans is 180°). The CMAP term is a
function of the pair,

\f[
U_\text{CMAP} = \sum_{\text{maps}} E(\phi, \psi) .
\f]

For the backbone of a protein the five atoms are the carbonyl carbon of the
previous residue, and the nitrogen, the \f$\alpha\f$-carbon and the carbonyl
carbon of the residue, and the nitrogen of the next one; \f$\phi\f$ and
\f$\psi\f$ are the angles of the Ramachandran plot.

## Why a map

A sum of one-dimensional [torsion](\ref potential_torsions) potentials,
\f$u_1(\phi) + u_2(\psi)\f$, can only produce an energy surface that is a sum
of a function of \f$\phi\f$ and a function of \f$\psi\f$. The true surface is
not of that form: the energy of rotating about one bond depends on the
conformation about the next, because groups on both sides come close to each
other in some combinations of the angles and not in others. The 1-4 and
longer-range non-bonded interactions include part of this coupling, but not
accurately enough for the balance between secondary structures (helices and
sheets), which depends on energy differences of a few tenths of a kcal/mol.

Fitting an analytical two-dimensional function would need many terms. A
correction map sidesteps the problem. The difference between a
quantum-chemical energy surface of a model compound and the surface the rest
of the force field gives is computed on a grid of \f$(\phi, \psi)\f$ values
and simply stored. The term is then, by construction, exactly the correction
that makes the force field reproduce the reference surface.

## The map

A map is an \f$n \times n\f$ grid of energies, with both angles running from
\f$-180^\circ\f$ in steps of \f$360^\circ/n\f$ (15° for the usual
\f$n = 24\f$). It is periodic in both directions.

Between the grid points the energy is obtained with a **periodic bicubic
spline**. Within each cell of the grid the energy is a polynomial of third
degree in each angle,

\f[
E(\phi, \psi) = \sum_{i=0}^{3}\sum_{j=0}^{3} c_{ij}\,t^i\,u^j ,
\f]

with \f$t\f$ and \f$u\f$ the fractional positions in the cell. The sixteen
coefficients of a cell follow from the energy, its two first derivatives and
its cross derivative at the four corners. The derivatives at the grid points
are not given but derived from the energies: the first derivatives from
periodic cubic splines along each angle, the cross derivative from a spline
through \f$\partial E/\partial\phi\f$ along \f$\psi\f$. The result is a
surface with continuous first derivatives everywhere, so the forces are
continuous, which molecular dynamics needs. This is the interpolation used by
CHARMM and OpenMM, so a map gives the same energies as in those programs.

## In the input

The maps are given in `force_field.json` under `"CMAPs"`, each with a name,
the resolution \f$n\f$ and the \f$n \times n\f$ energies in K, in the order of
the AMBER and CHARMM parameter files (the first angle \f$\phi\f$ is the slow
index, \f$\psi\f$ the fast one):

```json
"CMAPs" : [
  {"Name" : "CMAP_2", "Resolution" : 24, "Energies" : [ ... ]}
]
```

A molecule file can define maps of its own in the same way; a map of the
molecule file takes precedence over a force-field map of the same name.

The molecule file lists where the maps apply, with `"CMAPTorsions"`. Each
entry gives the five atoms by their indices in the molecule and the name of
the map:

```json
"CMAPTorsions" : [
  [[4, 6, 8, 14, 16], "CMAP_2"]
]
```

## Where the term enters

The CMAP energy is part of the intramolecular energy and is listed separately
in the output. It is included

- in Monte Carlo, where configurational-bias growth takes a term into account
  as soon as all five of its atoms have been placed;
- in molecular dynamics, with its forces and its contribution to the pressure
  tensor;
- in the analytic second derivatives used by the minimization and the
  normal-mode analysis.

## Things to know

- A map belongs to the force field it was made for. It is a *correction* to
  that force field's torsions, charges and 1-4 scaling, and has no meaning
  with another.
- The ordinary torsion terms on \f$\phi\f$ and \f$\psi\f$ remain; the map is
  added to them.
- Different residues use different maps (glycine and proline differ from the
  others), which is why maps are referred to by name.

## References

- A.D. MacKerell Jr., M. Feig and C.L. Brooks III, *Extending the treatment
  of backbone energetics in protein force fields: Limitations of gas-phase
  quantum mechanics in reproducing protein conformational distributions in
  molecular dynamics simulations*, J. Comput. Chem. **25**, 1400 (2004).
- C. Tian, K. Kasavajhala, K.A.A. Belfon, L. Raguette, H. Huang, A.N. Migues,
  J. Bickel, Y. Wang, J. Pincay, Q. Wu and C. Simmerling, *ff19SB:
  Amino-Acid-Specific Protein Backbone Parameters Trained against Quantum
  Mechanics Energy Surfaces in Solution*,
  J. Chem. Theory Comput. **16**, 528 (2020).
- G.N. Ramachandran, C. Ramakrishnan and V. Sasisekharan, *Stereochemistry of
  polypeptide chain configurations*, J. Mol. Biol. **7**, 95 (1963).
- W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery, *Numerical
  Recipes*, 3rd edition, Cambridge University Press (2007), section on
  bicubic interpolation.
