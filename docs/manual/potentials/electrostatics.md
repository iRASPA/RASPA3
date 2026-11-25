# Electrostatics
\page potential_electrostatics Electrostatics

The uneven distribution of the electrons in a molecule or a framework is
represented by *partial charges* on the interaction sites. Two charges
interact with the Coulomb potential,

\f[
u(r) = \frac{q_i\,q_j}{4\pi\epsilon_0\,r} ,
\f]

and the electrostatic energy is the sum over all pairs of charges. In polar
and ionic systems (water, CO<sub>2</sub> in a zeolite with cations, ionic
liquids) this term is as important as the dispersion.

## Where it is defined

The charge of a site is the `"charge"` of its pseudo-atom in
`force_field.json`, in units of the elementary charge. An atom of a molecule
file can override it with its own charge, so that atoms of one type carry
different charges, as in protein force fields. The charges of a framework can
also be taken from its CIF file or computed by charge equilibration
(`"UseChargesFrom"`).

The method is chosen per system with `"ChargeMethod"`. `"None"` skips the
electrostatics entirely, which is correct only when no site carries a charge.

## Why a cutoff does not work

The dispersion falls off as \f$r^{-6}\f$, and the part beyond a cutoff is a
small correction. The Coulomb potential falls off as \f$r^{-1}\f$, while the
number of charges at distance \f$r\f$ grows as \f$r^2\f$. The sum over a
periodic system is *conditionally convergent*: its value depends on the order
in which the terms are added, and truncating it at a cutoff gives an error
that does not vanish as the cutoff grows, together with artificial structure
at the cutoff distance.

## The Ewald summation

The Ewald method gives the energy of the infinite periodic system exactly, by
splitting each point charge into two parts. Around every charge a Gaussian
cloud of the *opposite* charge is added, of width \f$1/\alpha\f$. The charge
plus its screening cloud has no net charge, so its potential decays within a
short distance and can be summed in real space. The compensating clouds, with
the *same* sign as the charges, form a smooth periodic charge density, whose
energy is summed in reciprocal space, where a smooth function has a rapidly
converging Fourier series. The energy is

\f[
U = U_\text{real} + U_\text{rec} - U_\text{self} - U_\text{excl} .
\f]

**Real space.** The screened interaction, within the cutoff \f$r_c\f$:

\f[
U_\text{real} = \frac{1}{4\pi\epsilon_0}\sum_{i<j} q_i q_j\,\frac{\mathrm{erfc}(\alpha r_{ij})}{r_{ij}} .
\f]

**Reciprocal space.** The energy of the Gaussian clouds, as a sum over the
wave vectors \f$\mathbf{k}\f$ of the periodic box:

\f[
U_\text{rec} = \frac{1}{2V\epsilon_0}\sum_{\mathbf{k}\ne 0}
\frac{e^{-k^2/4\alpha^2}}{k^2}\,\left|S(\mathbf{k})\right|^2, \qquad
S(\mathbf{k}) = \sum_j q_j\,e^{i\mathbf{k}\cdot\mathbf{r}_j} .
\f]

The *structure factor* \f$S(\mathbf{k})\f$ is a sum over atoms, not over
pairs. When one molecule moves, only its contribution changes, so the stored
structure factors are updated rather than recomputed. This is what makes the
method affordable in Monte Carlo.

**Self term.** The reciprocal sum includes the interaction of each Gaussian
cloud with itself, which must be removed:

\f[
U_\text{self} = \frac{\alpha}{\sqrt{\pi}}\,\frac{1}{4\pi\epsilon_0}\sum_i q_i^2 .
\f]

**Exclusions.** The reciprocal sum also includes the interaction between the
clouds of atoms of the *same* molecule that should not interact (bonded
neighbours, atoms of a rigid molecule). It is removed pair by pair:

\f[
U_\text{excl} = \frac{1}{4\pi\epsilon_0}\sum_{\substack{i<j\\ \text{excluded}}} q_i q_j\,\frac{\mathrm{erf}(\alpha r_{ij})}{r_{ij}} .
\f]

**The framework.** For a rigid framework the structure factor of the
framework atoms is computed once and stored; only the interaction of the
moving molecules with it, and with each other, is evaluated during the run.

## Accuracy

The splitting parameter \f$\alpha\f$ moves work between the two sums: a
large \f$\alpha\f$ makes the real-space term decay faster but needs more wave
vectors. The result does not depend on it once both sums are converged.
`"EwaldPrecision"` sets the relative error \f$\epsilon\f$, from which
\f$\alpha\f$ and the wave vectors follow, given the Coulomb cutoff
\f$r_c\f$:

\f[
t = \sqrt{|\ln(\epsilon\,r_c)|},\qquad
\alpha = \frac{\sqrt{|\ln(\epsilon\,r_c\,t)|}}{r_c},\qquad
k_\text{max} = 2\alpha\sqrt{-\ln\!\left[\epsilon\,r_c\,(2t\alpha)^2\right]} .
\f]

The reciprocal sum runs over all wave vectors with \f$|\mathbf{k}| <
k_\text{max}\f$, for any shape of the cell. The parameters are derived again
when the box changes. They can also be given explicitly with
`"EwaldParameters"`.

The cost of the reciprocal sum grows with the number of wave vectors, which
grows with the volume. For large boxes a Coulomb cutoff of about half the
shortest box length keeps that number down, and for very large boxes the
[spatial-decomposition driver](\ref driver_molecular_dynamics_spatial_decomposition)
replaces the sum by a particle-mesh method.

## Net charge

The Ewald sum assumes the box is neutral. For a box with a net charge the
\f$\mathbf{k} = 0\f$ term diverges; leaving it out is equivalent to adding a
uniform neutralizing background charge. The energy then contains an
additional term that depends on the net charge and the volume, which the
moves that change the net charge take into account. Inserting ions as neutral
pairs or groups avoids the issue; see the pair and group swaps on the
[Monte Carlo moves](\ref mc_moves) page.

## Cutoff-based alternatives

In a dense, disordered system the charges screen each other, and the
effective interaction is short-ranged. Several methods use this to avoid the
reciprocal sum. They are cheaper and scale linearly, but they are
approximations whose accuracy must be checked against the Ewald sum for the
system at hand, and they need a large cutoff.

| `"ChargeMethod"` | Idea |
|---|---|
| `"Ewald"` | The exact sum described above. The default. |
| `"Coulomb"` | The bare Coulomb potential, truncated at the cutoff. For isolated molecules and testing. |
| `"Wolf"` | Wolf *et al.*: the charge within the cutoff sphere is neutralized by image charges on its surface, which shifts the damped potential \f$\mathrm{erfc}(\alpha r)/r\f$ to zero at the cutoff. |
| `"DampedShiftedForce"` | Fennell and Gezelter: the Wolf potential with the *force* shifted to zero at the cutoff as well, which makes it suitable for molecular dynamics. |
| `"ModifiedShiftedForce"` | Waibel, Feinler and Gross: a shifted-force form with a modified short-range part, with the parameter `"ModifiedShiftedForceBeta"`. |
| `"ZeroDipole"` | Fukuda: the cutoff sphere is neutralized in its dipole as well as its charge. |

These methods share the bookkeeping of the Ewald method: a self term and
exclusion corrections for the intramolecular pairs.

## Scaling for a fractional molecule

The charges of a fractional molecule are scaled with the coupling parameter.
As for the dispersion, scaling the potential alone leaves a singularity at
\f$r = 0\f$, where a partially coupled molecule may overlap with another. The
distance in the real-space term is therefore offset,

\f[
u_\lambda(r) = \lambda\,\frac{q_i q_j}{4\pi\epsilon_0}\,\frac{\mathrm{erfc}(\alpha s)}{s},
\qquad s = r + \delta\,(1 - \lambda) ,
\f]

which is finite at \f$r = 0\f$ for \f$\lambda < 1\f$ and becomes the ordinary
term at \f$\lambda = 1\f$. The reciprocal sum uses the scaled charges.

## Things to check

- The system should be neutral, or its net charge intended.
- With electrostatics the truncation of the *dispersion* and the Ewald
  precision are independent settings; a low Ewald precision shows up as an
  energy drift in molecular dynamics.
- The Coulomb cutoff must not exceed half the smallest perpendicular width of
  the box.

## References

- P.P. Ewald, *Die Berechnung optischer und elektrostatischer
  Gitterpotentiale*, Ann. Phys. **369**, 253 (1921).
- S.W. de Leeuw, J.W. Perram and E.R. Smith, *Simulation of electrostatic
  systems in periodic boundary conditions. I. Lattice sums and dielectric
  constants*, Proc. R. Soc. Lond. A **373**, 27 (1980).
- D. Frenkel and B. Smit, *Understanding Molecular Simulation*, 3rd edition,
  Academic Press (2023), chapter on long-range interactions.
- D. Wolf, P. Keblinski, S.R. Phillpot and J. Eggebrecht, *Exact method for
  the simulation of Coulombic systems by spherically truncated, pairwise
  \f$r^{-1}\f$ summation*, J. Chem. Phys. **110**, 8254 (1999).
- C.J. Fennell and J.D. Gezelter, *Is the Ewald summation still necessary?
  Pairwise alternatives to the accepted standard for long-range
  electrostatics*, J. Chem. Phys. **124**, 234104 (2006).
- I. Fukuda, Y. Yonezawa and H. Nakamura, *Molecular dynamics scheme for
  precise estimation of electrostatic interaction via zero-dipole summation
  principle*, J. Chem. Phys. **134**, 164107 (2011).
- T.J.H. Vlugt, E. García-Pérez, D. Dubbeldam, S. Ban and S. Calero,
  *Computing the Heat of Adsorption using Molecular Simulations: The Effect of
  Strong Coulombic Interactions*, J. Chem. Theory Comput. **4**, 1107 (2008).
