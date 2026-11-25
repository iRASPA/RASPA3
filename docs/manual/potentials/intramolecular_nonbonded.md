# Non-bonded interactions within a molecule
\page potential_intramolecular_nonbonded Non-bonded interactions within a molecule

Two atoms of the same flexible molecule that are far apart along the chain
interact with each other like atoms of different molecules: with the
[dispersion](\ref potential_dispersion) and the
[electrostatic](\ref potential_electrostatics) potential. Atoms that are close
along the chain do not, because their interaction is already contained in the
bonded terms. This page describes which pairs interact.

## Counting bonds

Pairs of atoms within a molecule are named by the number of bonds between
them:

| Pair | Separation | Covered by | Non-bonded interaction |
|---|---|---|---|
| 1-2 | one bond | the [bond](\ref potential_bonds) | excluded |
| 1-3 | two bonds | the [bend](\ref potential_bends) (and [Urey-Bradley](\ref potential_urey_bradley)) | excluded |
| 1-4 | three bonds | the [torsion](\ref potential_torsions), partly | scaled |
| 1-5 and beyond | four or more bonds | nothing | full |

The 1-2 and 1-3 distances are fixed, or nearly so, by the bond and the bend.
Their non-bonded interaction would be a large, nearly constant energy at a
distance far inside the repulsive wall, and force fields leave it out.

The 1-4 distance varies with the torsion angle, and the non-bonded
interaction between the two outer atoms is part of what gives the torsion its
energy profile. Force fields differ in how the two are combined, which is why
the 1-4 interaction has scaling factors.

## The 1-4 scaling

`"Intra14VanDerWaalsScalingValue"` and `"Intra14ChargeChargeScalingValue"` in
the molecule file multiply the van der Waals and the Coulomb interaction of
the 1-4 pairs:

| Convention | van der Waals | Coulomb |
|---|---|---|
| TraPPE and united-atom alkane models (the default) | 0 | 0 |
| AMBER | 0.5 | 0.8333 |
| OPLS | 0.5 | 0.5 |
| CHARMM | 1, with separate 1-4 parameters | 1 |

With a factor of zero the torsion potential carries the whole barrier. With
a non-zero factor the torsion parameters were fitted *together with* the
scaled 1-4 interaction, so the factors belong to the force field as much as
the torsion parameters do.

**Separate 1-4 parameters.** CHARMM uses different Lennard-Jones parameters
for 1-4 pairs. They are given as `"parameters14"` in the interactions of the
force-field file. When any interaction defines them, a second table is built,
with the same mixing rule, and used for the 1-4 pairs of every component
(still multiplied by the scaling value). A site without 1-4 parameters uses
its regular ones.

## Rigid bodies

Atoms within one rigid body have fixed distances, so their interaction is a
constant and is left out. This holds for a rigid molecule as a whole, which
has no intramolecular interactions at all, and for each rigid part of a
semi-flexible molecule. Atoms of two *different* rigid bodies of one molecule
interact according to the number of bonds between them.

## How it is evaluated

The interacting pairs of a molecule use the van der Waals parameters of the
force field, truncated at the molecule-molecule cutoff, and the partial
charges of the sites with the same charge method as the intermolecular
interactions.

For the electrostatics the excluded and scaled pairs need a correction. The
reciprocal-space sum of the Ewald method includes *all* pairs of charges,
whatever molecule they belong to. The pairs that should not interact, or only
partly, are corrected by subtracting their reciprocal-space contribution,

\f[
U_\text{excl} = \frac{1}{4\pi\epsilon_0}\sum_{\text{pairs}} (1 - s_{ij})\,q_i q_j\,\frac{\mathrm{erf}(\alpha r_{ij})}{r_{ij}} ,
\f]

with \f$s_{ij}\f$ the scaling of the pair (zero for excluded pairs). The
cutoff-based charge methods have the corresponding correction.

## Specifying the pairs

By default the lists of interacting pairs are derived from the connectivity
of the molecule (`"VanDerWaals" : "auto"`), following the rules of the table
above.

## Why it matters for sampling

For a long flexible molecule the intramolecular non-bonded interaction is
what keeps the chain from overlapping with itself. In configurational-bias
growth it is part of the energy by which each trial position is weighted, so
that a chain is grown avoiding its own earlier segments.

## References

- M.G. Martin and J.I. Siepmann, *Transferable Potentials for Phase
  Equilibria. 1. United-Atom Description of n-Alkanes*,
  J. Phys. Chem. B **102**, 2569 (1998).
- W.D. Cornell, P. Cieplak, C.I. Bayly, I.R. Gould, K.M. Merz, D.M. Ferguson,
  D.C. Spellmeyer, T. Fox, J.W. Caldwell and P.A. Kollman, *A Second
  Generation Force Field for the Simulation of Proteins, Nucleic Acids, and
  Organic Molecules*, J. Am. Chem. Soc. **117**, 5179 (1995).
- W.L. Jorgensen, D.S. Maxwell and J. Tirado-Rives, *Development and Testing
  of the OPLS All-Atom Force Field on Conformational Energetics and
  Properties of Organic Liquids*, J. Am. Chem. Soc. **118**, 11225 (1996).
- A.D. MacKerell Jr. *et al.*, *All-Atom Empirical Potential for Molecular
  Modeling and Dynamics Studies of Proteins*,
  J. Phys. Chem. B **102**, 3586 (1998).
