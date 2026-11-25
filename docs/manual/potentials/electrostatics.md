# Electrostatics
\page potential_electrostatics Electrostatics

The uneven distribution of the electrons in a molecule or a framework is
represented by *partial charges* on the interaction sites. Two charges
interact with the Coulomb potential,

\f[
u(r) = C\,\frac{q_i\,q_j}{r}, \qquad C = \frac{1}{4\pi\epsilon_0} ,
\f]

and the electrostatic energy is the sum over all pairs of charges. In polar
and ionic systems (water, CO<sub>2</sub> in a zeolite with cations, ionic
liquids) this term is as important as the dispersion. With charges in units
of the elementary charge, distances in ångström and energies as
\f$E/k_B\f$ in kelvin, \f$C \approx 167101\f$ K Å.

This page first explains why the sum is difficult, then derives the Ewald
summation and describes how each of its terms is evaluated, and ends with
the other methods that are implemented.

- [Where the charges and the method are defined](#electrostatics-input)
- [The problem with a cutoff](#electrostatics-problem)
- [The Ewald summation](#electrostatics-ewald)
- [How the terms are evaluated](#electrostatics-terms)
- [Accuracy and cost](#electrostatics-accuracy)
- [Monte Carlo: energy differences](#electrostatics-monte-carlo)
- [Forces and pressure](#electrostatics-forces)
- [A box with a net charge](#electrostatics-net-charge)
- [Fractional molecules](#electrostatics-fractional)
- [The particle-mesh method](#electrostatics-pme)
- [The other methods](#electrostatics-other-methods)

## Where the charges and the method are defined {#electrostatics-input}

The charge of a site is the `"charge"` of its pseudo-atom in
`force_field.json`. An atom of a molecule file can override it with its own
charge, so that atoms of one type carry different charges, as in protein
force fields. The charges of a framework can also be taken from its CIF file
or computed by charge equilibration (`"UseChargesFrom"`).

| Key | Where | Meaning |
|---|---|---|
| `"ChargeMethod"` | system | The method: `"Ewald"` (default), `"None"`, or one of the [other methods](#electrostatics-other-methods) |
| `"CutOffCoulomb"` | force field | The cutoff \f$r_c\f$ of the real-space sum [Å] |
| `"EwaldPrecision"` | force field | The relative precision: a number, or `"High"` (\f$10^{-8}\f$), `"Medium"` (\f$10^{-6}\f$, the default), `"Low"` (\f$10^{-4}\f$) |
| `"EwaldParameters"` | force field | \f$[\alpha, k_x, k_y, k_z]\f$: the splitting parameter and the largest wave-vector indices, given by hand |
| `"OmitEwaldFourier"` | force field | Leaves out the reciprocal sum. For testing. |

`"None"` skips the electrostatics entirely, which is correct only when no
site carries a charge.

## The problem with a cutoff {#electrostatics-problem}

A simulation box is surrounded by its periodic images, and a charge interacts
with all charges in all images,

\f[
U = \frac{C}{2}\sum_{\mathbf{n}}{}^{'}\sum_{i,j} \frac{q_i\,q_j}{|\mathbf{r}_{ij} + \mathbf{n}|} ,
\f]

where \f$\mathbf{n}\f$ runs over the lattice vectors of the box and the prime
excludes \f$i = j\f$ for \f$\mathbf{n} = 0\f$.

The dispersion falls off as \f$r^{-6}\f$, and the part beyond a cutoff is a
small correction. The Coulomb potential falls off as \f$r^{-1}\f$, while the
number of charges at distance \f$r\f$ grows as \f$r^2\f$: the contribution of
a shell does not decrease with its radius. For a neutral system the positive
and negative contributions cancel, but only *conditionally*: the value of the
sum depends on the order in which the terms are added.

Truncating the sum at a cutoff sphere therefore fails in a specific way. The
sphere around a charge contains, at any instant, an unbalanced amount of
charge, because the cutoff cuts through molecules and through the screening
cloud of every ion. That net charge is of the order of the square root of the
number of charges at the surface, and its energy does not vanish as the
cutoff grows. The result is a large error in the energy, artificial ordering
of the charges at the cutoff distance, and in molecular dynamics a force that
jumps whenever a charge crosses the sphere.

## The Ewald summation {#electrostatics-ewald}

The Ewald method computes the sum of the infinite periodic system without
truncating the interaction, by rewriting one slowly converging sum as two
rapidly converging ones.

### The idea

Around every point charge \f$q_i\f$ a Gaussian cloud of the *opposite*
total charge is added,

\f[
\rho_i^\text{screen}(\mathbf{r}) = -q_i \left(\frac{\alpha^2}{\pi}\right)^{3/2}
e^{-\alpha^2 |\mathbf{r} - \mathbf{r}_i|^2} .
\f]

Seen from a distance, a charge with its cloud is neutral, so its potential
decays within a length of the order of \f$1/\alpha\f$ and can be summed over
neighbours within a cutoff. To leave the system unchanged, clouds of the
*same* sign as the charges are added as well. Those form a smooth, periodic
charge density, and a smooth periodic function has a Fourier series that
converges quickly.

\f[
\underbrace{\text{point charges}}_{\text{long-ranged, singular}}
= \underbrace{\text{point charges} + \text{screening clouds}}_{\text{short-ranged: real space}}
+ \underbrace{\text{compensating clouds}}_{\text{smooth: reciprocal space}}
\f]

### The real-space part

The potential of a unit point charge with its screening cloud follows from
Poisson's equation. The cloud contributes \f$-\mathrm{erf}(\alpha r)/r\f$,
so the screened potential is

\f[
\frac{1}{r} - \frac{\mathrm{erf}(\alpha r)}{r} = \frac{\mathrm{erfc}(\alpha r)}{r} ,
\f]

with \f$\mathrm{erfc}\f$ the complementary error function, which decays as
\f$e^{-\alpha^2 r^2}\f$. The energy of the screened charges is

\f[
U_\text{real} = C \sum_{i<j} q_i\,q_j\,\frac{\mathrm{erfc}(\alpha r_{ij})}{r_{ij}} ,
\f]

summed over the pairs within \f$r_c\f$, with the nearest image.

### The reciprocal part

The compensating clouds have the charge density

\f[
\rho(\mathbf{r}) = \sum_{\mathbf{n}}\sum_j q_j \left(\frac{\alpha^2}{\pi}\right)^{3/2}
e^{-\alpha^2 |\mathbf{r} - \mathbf{r}_j - \mathbf{n}|^2} ,
\f]

which is periodic and can be written as a Fourier series over the wave
vectors \f$\mathbf{k}\f$ of the reciprocal lattice, with the coefficients

\f[
\tilde\rho(\mathbf{k}) = \frac{1}{V}\,e^{-k^2/4\alpha^2}\,S(\mathbf{k}), \qquad
S(\mathbf{k}) = \sum_j q_j\,e^{i\mathbf{k}\cdot\mathbf{r}_j} .
\f]

\f$S(\mathbf{k})\f$ is the *structure factor* of the charges. The Gaussian
factor is what makes the series converge: the Fourier transform of a
Gaussian of width \f$1/\alpha\f$ is a Gaussian of width \f$2\alpha\f$.

In Fourier space Poisson's equation \f$-\nabla^2\phi = \rho/\epsilon_0\f$ is
algebraic, \f$k^2\,\tilde\phi(\mathbf{k}) = \tilde\rho(\mathbf{k})/\epsilon_0\f$.
The energy of the point charges in the potential of the clouds,
\f$\frac{1}{2}\sum_i q_i\,\phi(\mathbf{r}_i)\f$, then becomes

\f[
U_\text{rec} = \frac{1}{2V\epsilon_0}\sum_{\mathbf{k}\ne 0}
\frac{e^{-k^2/4\alpha^2}}{k^2}\,\lvert S(\mathbf{k})\rvert^2 .
\f]

The term \f$\mathbf{k} = 0\f$ is left out. For a neutral box
\f$S(0) = \sum_j q_j = 0\f$, and what remains of that term depends on the
shape and the surroundings of the infinite crystal of boxes. Leaving it out
corresponds to a crystal embedded in a conductor (*tin-foil* boundary
conditions), which is the usual choice and the one made here: there is no
surface-dipole term.

### The self term

The potential of the clouds at the position of charge \f$i\f$ includes the
potential of the cloud of \f$i\f$ itself. That interaction is not physical.
The potential at the centre of a Gaussian cloud of unit charge is
\f$\lim_{r\to0}\mathrm{erf}(\alpha r)/r = 2\alpha/\sqrt{\pi}\f$, so

\f[
U_\text{self} = -C\,\frac{\alpha}{\sqrt{\pi}}\sum_i q_i^2 .
\f]

It is a constant as long as the charges and the number of atoms do not
change. It matters when molecules are inserted or deleted, when a fractional
molecule changes its coupling, and when the box changes if \f$\alpha\f$
follows the box.

### The exclusion term

The reciprocal sum is a sum over all pairs of clouds, including pairs of
atoms of the *same* molecule that must not interact: bonded neighbours, and
all atoms of a rigid molecule. In real space such a pair is simply skipped.
From the reciprocal sum its contribution, the interaction of a point charge
with a cloud, is removed pair by pair:

\f[
U_\text{excl} = -C \sum_{\substack{i<j\\ \text{excluded}}} q_i\,q_j\,\frac{\mathrm{erf}(\alpha r_{ij})}{r_{ij}} .
\f]

This term has no cutoff, and it depends on the conformation of a flexible
molecule. The self term is its limit for \f$i = j\f$.

### The total

\f[
U_\text{Coulomb} = U_\text{real} + U_\text{rec} + U_\text{self} + U_\text{excl}
\f]

The splitting parameter \f$\alpha\f$ moves work between the two sums. A
large \f$\alpha\f$ gives narrow clouds: the real-space term decays quickly,
but the clouds are sharp and need many wave vectors. Once both sums are
converged, the result does not depend on \f$\alpha\f$. That is the best test
of an implementation and of a choice of parameters.

## How the terms are evaluated {#electrostatics-terms}

| Term | Evaluated | Reported as |
|---|---|---|
| Real space | With the pair interactions, within `"CutOffCoulomb"` | Framework-molecule and molecule-molecule Coulomb energy |
| Reciprocal | A sum over wave vectors of the squared structure factor | Ewald Fourier |
| Self | A sum over atoms | Ewald self |
| Exclusion | A sum over the excluded pairs of each molecule | Ewald exclusion |

**The real-space term** is computed in the same loops as the dispersion. For
a fully coupled pair \f$\mathrm{erfc}(\alpha r)/r\f$ is read from a table: a
cubic Hermite spline on a grid in \f$r^2\f$, so that a pair needs no error
function, exponential, square root or division. Its relative error is below
\f$10^{-9}\f$, and the force is the exact derivative of the interpolated
energy, so the table does not affect energy conservation. Pairs closer than
1 Å, scaled pairs, and the second derivatives use the closed form.

**The wave vectors.** For a box with lattice vectors \f$\mathbf{a}_i\f$ and
reciprocal vectors \f$\mathbf{b}_i\f$,

\f[
\mathbf{k} = 2\pi\,(h_1\mathbf{b}_1 + h_2\mathbf{b}_2 + h_3\mathbf{b}_3),
\qquad h_i \text{ integer} ,
\f]

for any shape of the cell. Since \f$S(-\mathbf{k})\f$ is the complex
conjugate of \f$S(\mathbf{k})\f$, half of the vectors is enough. The factors
\f$e^{i\mathbf{k}\cdot\mathbf{r}_j}\f$ are products of powers of three
complex numbers per atom, so they are built by multiplication, without a
trigonometric function per wave vector.

**The framework.** The structure factor of a rigid framework is computed
once and stored. During the run the total is the stored framework part plus
the part of the molecules, and the energy of the framework with itself, a
constant, is left out. The atoms of a flexible framework take part as the
atoms of molecules do. When
the framework-molecule interaction is read from an interpolation grid, the
grid holds the real-space term; the reciprocal part is computed as described
here.

**Excluded and included pairs within a molecule.** The excluded pairs are the
1-2 and 1-3 pairs of a flexible molecule, the pairs within a rigid part, and
all pairs of a rigid molecule; see
[non-bonded interactions within a molecule](\ref potential_intramolecular_nonbonded).
The other pairs of a flexible molecule interact like any pair: a real-space
term (scaled for a 1-4 pair) plus their share of the reciprocal sum.

**Polarization.** The electric field at the atoms, needed for the
[induced dipoles](\ref potential_polarization), is obtained from the same
real-space and reciprocal sums.

**Second derivatives.** All four terms have analytic second derivatives,
which the minimization and the normal-mode analysis use.

## Accuracy and cost {#electrostatics-accuracy}

`"EwaldPrecision"` sets the relative precision \f$\epsilon\f$. Together with
the Coulomb cutoff \f$r_c\f$ it determines \f$\alpha\f$ and the wave
vectors:

\f[
t = \sqrt{\lvert\ln(\epsilon\,r_c)\rvert},\qquad
\alpha = \frac{\sqrt{\lvert\ln(\epsilon\,r_c\,t)\rvert}}{r_c},\qquad
t_1 = \sqrt{-\ln\!\left[\epsilon\,r_c\,(2t\alpha)^2\right]},\qquad
k_\text{max} = 2\,\alpha\,t_1 .
\f]

The first two choose \f$\alpha\f$ so that \f$\mathrm{erfc}(\alpha r_c)/r_c\f$
is of the order of \f$\epsilon\f$ at the cutoff. The last choose
\f$k_\text{max}\f$ so that \f$e^{-k^2/4\alpha^2}\f$ is as small at the edge
of the reciprocal sum. For \f$r_c = 12\f$ Å and \f$\epsilon = 10^{-6}\f$ this
gives \f$\alpha = 0.265\f$ Å⁻¹.

The reciprocal sum runs over every wave vector inside the sphere
\f$\lvert\mathbf{k}\rvert < k_\text{max}\f$. The indices that are visited run
up to

\f[
\lvert h_i \rvert \le \left\lfloor \frac{\lvert\mathbf{a}_i\rvert\,\alpha\,t_1}{\pi} \right\rfloor ,
\f]

with \f$\lvert\mathbf{a}_i\rvert\f$ the length of lattice vector \f$i\f$.
These ranges are the "Ewald k-vectors" printed in the output. The error is
then the same fraction of the energy in a cubic, an elongated and a skewed
cell. The parameters are derived again whenever the box changes.

With `"EwaldParameters"` \f$\alpha\f$ and the index ranges are given by
hand. The sum then runs over all index triplets in the range, without the
spherical cut, and the values are kept when the box changes.

**Cost.** The real-space sum costs a number of operations proportional to
\f$N r_c^3\f$. The number of wave vectors is proportional to
\f$V k_\text{max}^3\f$, so it grows with the volume of the box and with
\f$\alpha^3\f$, and each needs a sum over all \f$N\f$ atoms. A larger cutoff
allows a smaller \f$\alpha\f$ and fewer wave vectors. For a large box a
Coulomb cutoff of about half the shortest box length keeps the number of wave
vectors down. With the best balance between the two sums the cost of a full
evaluation grows as \f$N^{3/2}\f$; beyond some tens of thousands of atoms the
[particle-mesh method](#electrostatics-pme) is faster.

## Monte Carlo: energy differences {#electrostatics-monte-carlo}

A Monte Carlo move changes a few atoms and needs an energy difference. The
real-space, self and exclusion terms are sums over atoms or pairs, and their
difference involves only the moved atoms. The reciprocal energy is a sum over
wave vectors of \f$\lvert S(\mathbf{k})\rvert^2\f$, which is not a sum over
pairs; but the structure factor itself is a sum over atoms. So the structure
factor of the current configuration is *stored*, and for a trial move

\f[
S_\text{new}(\mathbf{k}) = S_\text{old}(\mathbf{k})
+ \sum_{j\,\in\,\text{moved}} \left[ q_j^\text{new}\,e^{i\mathbf{k}\cdot\mathbf{r}_j^\text{new}}
- q_j^\text{old}\,e^{i\mathbf{k}\cdot\mathbf{r}_j^\text{old}} \right] ,
\f]

\f[
\Delta U_\text{rec} = \frac{1}{2V\epsilon_0}\sum_{\mathbf{k}\ne 0}
\frac{e^{-k^2/4\alpha^2}}{k^2}
\left[\lvert S_\text{new}(\mathbf{k})\rvert^2 - \lvert S_\text{old}(\mathbf{k})\rvert^2\right] .
\f]

The cost is the number of wave vectors times the number of moved atoms,
independent of the size of the system. The trial structure factor is kept
next to the stored one and replaces it when the move is accepted.

- An **insertion** has no old atoms and a **deletion** no new ones; the self
  and exclusion terms of the molecule are part of the difference.
- A move of **part of a molecule** (a bead displacement or a crankshaft of a
  long chain) evaluates the difference for the moved atoms only, and the
  exclusion term for the pairs that contain a moved atom.
- A **volume change** moves every atom and changes the wave vectors, so the
  structure factor is computed anew.
- In configurational-bias growth the reciprocal energy is not part of the
  selection of trial positions. It is added as a correction in the
  acceptance rule, which keeps the growth cheap and the sampling exact.

## Forces and pressure {#electrostatics-forces}

The force on atom \f$i\f$ from the real-space term is that of a pair
potential. From the reciprocal term it is

\f[
\mathbf{F}_i^\text{rec} = \frac{q_i}{V\epsilon_0}\sum_{\mathbf{k}\ne 0}
\frac{e^{-k^2/4\alpha^2}}{k^2}\;\mathbf{k}\;
\mathrm{Im}\!\left[S^*(\mathbf{k})\,e^{i\mathbf{k}\cdot\mathbf{r}_i}\right] ,
\f]

which depends on the structure factor of the whole system but only on the
position of atom \f$i\f$. The force on one molecule can therefore be obtained
from a stored structure factor at the cost of that molecule alone, which the
force-biased Monte Carlo moves use. The exclusion term gives a pair force
within a molecule; the self term gives none.

The contribution of the reciprocal term to the pressure tensor follows from
the derivative of the energy with respect to a strain
\f$\varepsilon_{\alpha\beta}\f$ of the box, which changes both the wave
vectors and the volume:

\f[
-\frac{\partial U_\text{rec}}{\partial \varepsilon_{\alpha\beta}} = \frac{1}{2V\epsilon_0}\sum_{\mathbf{k}\ne 0}
\frac{e^{-k^2/4\alpha^2}}{k^2}\,\lvert S(\mathbf{k})\rvert^2
\left[\delta_{\alpha\beta} - 2\left(\frac{1}{k^2} + \frac{1}{4\alpha^2}\right) k_\alpha k_\beta\right] .
\f]

The real-space and exclusion terms contribute as pair interactions. For rigid
molecules the atomic expression is corrected to the molecular
(centre-of-mass) pressure.

## A box with a net charge {#electrostatics-net-charge}

The derivation assumes a neutral box. With a net charge
\f$Q = \sum_j q_j \ne 0\f$ the energy of the infinite system diverges.
Leaving out \f$\mathbf{k} = 0\f$ is then equivalent to adding a uniform
background charge that neutralizes the box, and the energy becomes finite
but contains an artefact: the interaction of the net charge with its own
periodic images and with the background. That term depends on \f$Q\f$ and on
the volume, and it distorts the pressure and the free energy of inserting a
charged molecule.

It is removed with the correction of Bogusz, Cheatham and Brooks,

\f[
U_\text{net} = -Q^2 \left[\frac{1}{2V\epsilon_0}\sum_{\mathbf{k}\ne 0}
\frac{e^{-k^2/4\alpha^2}}{k^2} \;-\; C\,\frac{\alpha}{\sqrt{\pi}}\right] ,
\f]

where the term in brackets is the Ewald energy of a single unit charge in the
box. \f$Q\f$ is the net charge of the framework plus the molecules, with the
charges of fractional molecules scaled. The correction is zero for a neutral
box, is included in the reported Ewald Fourier energy, and its change is part
of the energy difference of every move that changes the net charge.

A framework with a net charge that is balanced by mobile cations is a neutral
box. To insert ions at constant neutrality, use the pair and group swap moves
of the [Monte Carlo moves](\ref mc_moves) page.

## Fractional molecules {#electrostatics-fractional}

A *fractional* molecule interacts with a strength set by a coupling
parameter. The van der Waals interaction is switched on over the first half
of the range of \f$\lambda\f$ and the electrostatics over the second half, so
that a charge is never coupled before the repulsion that keeps other atoms
away from it. Each atom carries a Coulomb scaling \f$s_i\f$ between 0 and 1,
and a pair is scaled with \f$\lambda_t = s_i\,s_j\f$.

Scaling the potential alone is not enough: a partly coupled molecule can
overlap with another atom, and \f$\lambda_t\,\mathrm{erfc}(\alpha r)/r\f$
still diverges at \f$r = 0\f$. The distance is therefore offset,

\f[
u_\lambda(r) = \lambda_t\,C\,q_i\,q_j\,\frac{\mathrm{erfc}(\alpha s)}{s},
\qquad s = r + \delta\,(1 - \lambda_t), \qquad \delta = 0.01\ \text{Å} ,
\f]

which is finite at \f$r = 0\f$ for \f$\lambda_t < 1\f$ and is the ordinary
term at \f$\lambda_t = 1\f$. The exclusion term has the same form, with
\f$\mathrm{erf}\f$. The reciprocal and self terms use the scaled charges
\f$s_i q_i\f$.

The derivative \f$\partial U/\partial\lambda\f$, needed for thermodynamic
integration, is accumulated for all four terms. For the reciprocal term this
needs, per wave vector, the sum \f$\sum_j q_j e^{i\mathbf{k}\cdot\mathbf{r}_j}\f$
over the atoms of the fractional molecule alone, which is stored next to the
structure factor.

## The particle-mesh method {#electrostatics-pme}

The direct reciprocal sum visits every atom for every wave vector. The
[spatial-decomposition driver](\ref driver_molecular_dynamics_spatial_decomposition),
made for large systems, replaces it by the *smooth particle-mesh Ewald*
method of Essmann *et al.*:

1. The charges are assigned to the points of a regular mesh, each charge
   spread over the nearest points with a cardinal B-spline of order \f$p\f$.
2. A fast Fourier transform of the mesh gives an approximation of the
   structure factor at all wave vectors at once.
3. The transform is multiplied by the *influence function*: the factor
   \f$e^{-k^2/4\alpha^2}/k^2\f$ of the Ewald sum, corrected for the smoothing
   by the B-splines.
4. A transform back gives the potential on the mesh, from which the forces on
   the atoms are interpolated with the derivatives of the same B-splines.

The cost grows as \f$N \log N\f$ instead of \f$N^{3/2}\f$. The same
\f$\alpha\f$ is used, so the real-space, self, exclusion and net-charge terms
are unchanged, and the result converges to the Ewald sum as the mesh is
refined.

| Key | Default | Meaning |
|---|---|---|
| `"PPPMMeshSpacing"` | 1.0 | The target spacing of the mesh [Å] |
| `"PPPMInterpolationOrder"` | 5 | The order of the B-spline (3 to 7) |

With the default order the relative error of the reciprocal energy and
forces is about \f$10^{-4}\f$ at a spacing of 1 Å and \f$10^{-6}\f$ at
0.5 Å, for \f$\alpha \approx 0.3\f$ Å⁻¹. Triclinic cells are supported, and
the influence function is recomputed when the cell changes.

## The other methods {#electrostatics-other-methods}

In a dense, disordered system the charges screen each other, and the
effective interaction of two charges is short-ranged even though the bare
potential is not. Several methods use this to do without the reciprocal sum.
The pair potential is modified so that the cutoff sphere is neutral, and
everything is a sum over pairs within the cutoff. These methods are cheaper,
scale linearly with the number of atoms and are simple to run in parallel.
They are approximations, and their result depends on the damping and the
cutoff.

| `"ChargeMethod"` | Also accepted | Method |
|---|---|---|
| `"Ewald"` | | The Ewald summation. The default. |
| `"Coulomb"` | | The bare potential, truncated |
| `"Wolf"` | | Damped, shifted potential |
| `"DampedShiftedForce"` | `"DSF"`, `"Damped-Shifted-Force"` | Damped, shifted potential and force |
| `"ModifiedShiftedForce"` | `"MSF"`, `"mDSF"`, `"Modified-Shifted-Force"` | Shifted force with an exponential switch |
| `"ZeroDipole"` | `"ZD"`, `"Zero-Dipole"` | Neutralized charge and dipole |
| `"None"` | | No electrostatics |

In the formulas below \f$B(r) = \mathrm{erfc}(\alpha r)/r\f$ is the damped
potential, \f$B'\f$ its derivative, and the pair energy is
\f$C\,q_i\,q_j\,V(r)\f$ for \f$r < r_c\f$ and zero beyond.

### Truncated Coulomb

\f[
V(r) = \frac{1}{r}
\f]

The bare potential, cut at \f$r_c\f$, without any correction. It has the
problems described [above](#electrostatics-problem) and is meant for isolated
molecules, for tests, and for comparison.

### Wolf

Wolf *et al.* traced the error of a truncated sum to the net charge inside
the cutoff sphere, and removed it by placing, for every charge inside the
sphere, an opposite charge on its surface. For the pair potential this is a
shift to zero at the cutoff. Combined with the damping of the Ewald
real-space term:

\f[
V(r) = B(r) - B(r_c)
\f]

The energy converges quickly with the cutoff. The potential is continuous at
\f$r_c\f$, but the force is not, so the method suits Monte Carlo better than
molecular dynamics.

### Damped shifted force

Fennell and Gezelter shifted the force to zero at the cutoff as well:

\f[
V(r) = B(r) - B(r_c) - B'(r_c)\,(r - r_c)
\f]

Both the energy and the force are continuous, so the energy is conserved in
molecular dynamics. It is the most widely used of these methods. The linear
term changes the potential at all distances, which is the price for the
continuous force.

### Modified shifted force

Waibel, Feinler and Gross replaced the linear term by an exponential one,

\f[
V(r) = B(r) - B(r_c) - B'(r_c)\,\frac{e^{\beta (r - r_c)} - 1}{\beta} ,
\f]

which also makes the energy and the force zero at the cutoff, but the
correction decays towards short distances over a length \f$1/\beta\f$, so the
potential near contact is closer to the unshifted one. \f$\beta\f$ is set
with `"ModifiedShiftedForceBeta"` (default 0.3 Å⁻¹). For
\f$\beta \to 0\f$ the method becomes the damped shifted force.

### Zero dipole

Fukuda *et al.* required that the cutoff sphere is neutral in its dipole
moment as well as in its charge:

\f[
V(r) = B(r) - B(r_c) - B'(r_c)\,\frac{r^2 - r_c^2}{2\,r_c}
\f]

The energy and the force are zero at the cutoff. The correction is quadratic
in \f$r\f$, the potential of a uniform field gradient, and the method is
accurate for dipolar liquids such as water at moderate cutoffs.

### What these methods share

**The damping.** \f$\alpha\f$ is the same parameter as in the Ewald sum, and
is set in the same way: from `"EwaldPrecision"` and the cutoff, or with the
first value of `"EwaldParameters"`. Here it is not a numerical parameter
that drops out. A small \f$\alpha\f$ approaches the undamped shifted
potential, which needs a large cutoff; a large \f$\alpha\f$ converges at a
short cutoff but to the wrong energy, because the part of the interaction
that the reciprocal sum would restore is missing. Values around
0.2 Å⁻¹ with a cutoff of 12 Å or more are common in the literature.

**The self term.** Neutralizing the sphere includes the charge at its
centre. Each atom therefore has a constant energy, half the interaction of
the charge with its own neutralizing charge,

\f[
u_{\text{self},i} = -\frac{C\,q_i^2}{2}\,\lim_{r\to 0}\left[\frac{1}{r} - V(r)\right] .
\f]

| Method | \f$u_{\text{self},i}/(C\,q_i^2)\f$ |
|---|---|
| Wolf | \f$-\dfrac{\alpha}{\sqrt{\pi}} - \dfrac{\mathrm{erfc}(\alpha r_c)}{2 r_c}\f$ |
| Damped shifted force | \f$-\dfrac{\alpha}{\sqrt{\pi}}\left(1 + e^{-\alpha^2 r_c^2}\right) - \dfrac{\mathrm{erfc}(\alpha r_c)}{r_c}\f$ |
| Zero dipole | \f$-\dfrac{\alpha}{\sqrt{\pi}}\left(1 + \tfrac{1}{2}e^{-\alpha^2 r_c^2}\right) - \dfrac{3\,\mathrm{erfc}(\alpha r_c)}{4 r_c}\f$ |

The modified shifted force has a corresponding term that contains
\f$\beta\f$. The truncated Coulomb potential has none.

**Excluded pairs.** The neutralization counts every charge in the sphere,
including the other atoms of the same molecule. An excluded pair within the
cutoff therefore contributes

\f[
u_\text{excl}(r) = C\,q_i\,q_j\left[V(r) - \frac{1}{r}\right] ,
\f]

the pair term of the method without the bare interaction that the pair must
not have. It plays the role of the exclusion term of the Ewald sum, and is
reported in the same place, as is the self term.

**The isolated molecule.** With the Ewald sum, the self and exclusion terms
of a molecule are cancelled by the reciprocal sum when the molecule is alone
in a large box. Without a reciprocal sum they remain: a molecule has this
energy even in an empty box. It is part of the ideal-gas reference of the
insertion and deletion moves, so that a fugacity imposes the same density
with these methods as with the Ewald sum.

**Fractional molecules** have their charges scaled. Forces and the pressure
are available for all methods.

### Which to use

| | Ewald | Particle mesh | Shifted methods |
|---|---|---|---|
| Result | The periodic sum, to a set precision | The same, to the precision of the mesh | An approximation that depends on \f$\alpha\f$ and \f$r_c\f$ |
| Cost | \f$N^{3/2}\f$ | \f$N\log N\f$ | \f$N\f$ |
| Monte Carlo | Efficient, through the stored structure factor | Not for moves of single molecules | Efficient |
| Suited to | Any system; the reference | Large systems in molecular dynamics | Dense, homogeneous, disordered systems |
| Take care with | Many wave vectors in a large or elongated box | The mesh spacing | Interfaces, slabs, low density, ionic crystals, a framework with few strong charges |

The shifted methods rely on screening. Where a system is not screened
uniformly, as at an interface or in a sparsely filled pore, the assumption
does not hold and the error can be large. Validate a shifted method against
the Ewald sum for the system and the property of interest before relying on
it.

## Things to check

- The system should be neutral, or its net charge intended.
- The Coulomb cutoff must not exceed half the smallest perpendicular width of
  the box.
- The truncation of the dispersion and the Ewald precision are independent
  settings. A precision that is too low shows as an energy drift in molecular
  dynamics.
- The reported number of wave vectors: a very large number means that the
  cutoff is small for the box.

## References

- P.P. Ewald, *Die Berechnung optischer und elektrostatischer
  Gitterpotentiale*, Ann. Phys. **369**, 253 (1921).
- S.W. de Leeuw, J.W. Perram and E.R. Smith, *Simulation of electrostatic
  systems in periodic boundary conditions. I. Lattice sums and dielectric
  constants*, Proc. R. Soc. Lond. A **373**, 27 (1980).
- D. Frenkel and B. Smit, *Understanding Molecular Simulation*, 3rd edition,
  Academic Press (2023), chapter on long-range interactions.
- U. Essmann, L. Perera, M.L. Berkowitz, T. Darden, H. Lee and L.G. Pedersen,
  *A smooth particle mesh Ewald method*, J. Chem. Phys. **103**, 8577 (1995).
- S. Bogusz, T.E. Cheatham III and B.R. Brooks, *Removal of pressure and free
  energy artifacts in charged periodic systems via net charge corrections to
  the Ewald potential*, J. Chem. Phys. **108**, 7070 (1998).
- S. Stenberg and B. Stenqvist, *An Exact Ewald Summation Method in Theory
  and Practice*, J. Phys. Chem. A **124**, 3943 (2020).
- R. Hens, A. Rahbari, S. Caro-Ortiz, N. Dawass, M. Erdős, A. Poursaeidesfahani,
  H.S. Salehi, A.T. Celebi, M. Ramdin, O.A. Moultos, D. Dubbeldam and
  T.J.H. Vlugt, *Brick-CFCMC: Open Source Software for Monte Carlo
  Simulations of Phase and Reaction Equilibria Using the Continuous
  Fractional Component Method*, J. Chem. Inf. Model. **60**, 2678 (2020).
- D. Wolf, P. Keblinski, S.R. Phillpot and J. Eggebrecht, *Exact method for
  the simulation of Coulombic systems by spherically truncated, pairwise
  \f$r^{-1}\f$ summation*, J. Chem. Phys. **110**, 8254 (1999).
- C.J. Fennell and J.D. Gezelter, *Is the Ewald summation still necessary?
  Pairwise alternatives to the accepted standard for long-range
  electrostatics*, J. Chem. Phys. **124**, 234104 (2006).
- C. Waibel, M.S. Feinler and J. Gross, *A Modified Shifted Force Approach to
  the Wolf Summation*, J. Chem. Theory Comput. **15**, 572 (2019).
- I. Fukuda, Y. Yonezawa and H. Nakamura, *Molecular dynamics scheme for
  precise estimation of electrostatic interaction via zero-dipole summation
  principle*, J. Chem. Phys. **134**, 164107 (2011).
- T.J.H. Vlugt, E. García-Pérez, D. Dubbeldam, S. Ban and S. Calero,
  *Computing the Heat of Adsorption using Molecular Simulations: The Effect of
  Strong Coulombic Interactions*, J. Chem. Theory Comput. **4**, 1107 (2008).
