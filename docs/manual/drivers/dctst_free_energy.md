# dcTST: the free-energy profile
\page driver_dctst_free_energy dcTST: the free-energy profile

`"SimulationType" : "DynamicallyCorrectedTSTFreeEnergy"`

Dynamically corrected transition state theory (dcTST) computes diffusion
coefficients of molecules that hop between cages, in the regime where a hop
is too rare to be observed in molecular dynamics. The method has two halves.
This driver computes the first: the free-energy profile along the hopping
coordinate and, from it, the transition-state rate. The
[transmission driver](\ref driver_dctst) computes the correction for
recrossings.

## Hopping as a rare event

In a cage-type material (LTA, CHA, ERI and many others) a molecule spends a
long time in a cage and occasionally passes a narrow window to the next one.
On time scales much longer than a hop the motion is a random walk on the
lattice of cages, and the self-diffusion coefficient in one direction is

\f[
D_S = k_{A\to B}\;\lambda^2 ,
\f]

with \f$k_{A\to B}\f$ the rate of hopping to one particular neighbouring cage
and \f$\lambda\f$ the distance between the cage centres. The problem is the
rate. When the barrier is many \f$k_BT\f$ high, a simulation of nanoseconds
contains few hops or none.

## The Bennett-Chandler approach

The rate is written as a product of two factors that can each be computed
efficiently,

\f[
k_{A\to B} = \kappa \; k^\text{TST} .
\f]

Define a reaction coordinate \f$q\f$, here the position of one bead of the
molecule projected on the normal of the window plane, with the dividing
surface at \f$q^*\f$. Transition state theory assumes that every molecule
that reaches the dividing surface moving towards \f$B\f$ ends up in \f$B\f$.
The rate is then the equilibrium flux through the surface, per molecule in
\f$A\f$:

\f[
k^\text{TST} = \frac{1}{2}\,\langle|\dot q|\rangle\; P(q^*) , \qquad
P(q^*) = \frac{e^{-\beta F(q^*)}}{\displaystyle\int_\text{cage A} e^{-\beta F(q)}\,dq} .
\f]

\f$P(q^*)\f$ is the probability density of finding the molecule on top of the
barrier, relative to finding it anywhere in cage \f$A\f$, and
\f$\frac12\langle|\dot q|\rangle\f$ is the mean velocity towards \f$B\f$ of
the molecules that are there. For a bead that moves with the molecular mass
\f$m\f$ the thermal average is

\f[
\langle|\dot q|\rangle = \sqrt{\frac{2 k_B T}{\pi m}} .
\f]

For an off-centre bead of a rigid molecule the rotation contributes to the
velocity of the bead, and the average is taken over the sampled orientations.

The factor \f$\kappa\f$, between 0 and 1, corrects the assumption: it is the
fraction of crossings that actually commit to \f$B\f$. It needs dynamics and
is the subject of the other driver. The product \f$\kappa\,k^\text{TST}\f$ is
*exact* and does not depend on the choice of \f$q\f$ or \f$q^*\f$: a poor
choice gives a larger \f$k^\text{TST}\f$ and a smaller \f$\kappa\f$.

## The free-energy profile

The profile is the logarithm of the probability density of the reaction
coordinate,

\f[
\beta F(q) = -\ln P(q) .
\f]

At finite loading this is the profile of one molecule in the presence of all
the others, so the effect of the loading on the barrier is included. Because
the rate is proportional to \f$e^{-\beta F(q^*)}\f$, an error \f$\delta\f$ in
\f$\beta F(q^*)\f$ is a relative error \f$\delta\f$ in the diffusion
coefficient. The barrier top, which an unbiased simulation almost never
visits, must therefore be sampled deliberately.

## Umbrella sampling

One molecule is *tagged* and a bias \f$w_k(q)\f$ is added to its energy that
holds its reaction coordinate near a chosen value \f$q_k\f$,

\f[
w_k(q) = \frac{1}{2}\,\kappa_w\,(q - q_k)^2 + w_0(q) .
\f]

The harmonic term makes a *window*; a series of windows with centres from the
cage to the barrier top covers the whole coordinate. The optional term
\f$w_0\f$ is a flattening bias, an estimate of \f$-F(q)\f$ from the unbiased
stage, which lets each window sample its whole width instead of sliding down
the slope of the barrier. In window \f$k\f$ the tagged molecule samples

\f[
P_k(q) \propto P(q)\; e^{-\beta w_k(q)} ,
\f]

so every window gives the unbiased profile over its own range up to an
unknown constant: \f$\beta F(q) = -\ln P_k(q) - \beta w_k(q) + c_k\f$.

Two conditions make a window useful. Its histogram must overlap with those of
its neighbours, which sets the spacing. And it must *hold* the barrier top:
the bias must be stiffer than the barrier is curved,
\f$\beta\kappa_w > |\beta F''(q^*)|\f$, or the molecule slides off to one
side. The default stiffness, \f$\beta\kappa_w = 4/\text{spacing}^2\f$, makes
the biased distribution half a spacing wide, so that neighbouring windows
overlap at about one standard deviation each.

## Combining the windows

The constants \f$c_k\f$ are found by WHAM, with the reaction coordinate as
the only variable. With \f$H_k(q)\f$ the histogram of window \f$k\f$ and
\f$n_k\f$ its number of samples,

\f[
P(q) = \frac{\sum_k H_k(q)}{\sum_k n_k\; e^{\,g_k - \beta w_k(q)}} , \qquad
e^{-g_k} = \int P(q)\; e^{-\beta w_k(q)}\,dq ,
\f]

solved self-consistently (see the
[page on histogram reweighting](\ref driver_reweighted_histogram) for the
derivation). The bias is tabulated as a constant per fine bin, so the
reweighting is exact rather than approximate within a bin; for a bias that
depends on \f$q\f$ only, binned WHAM is then equivalent to MBAR.

The unbiased stage enters the same equations as one more state, with
\f$w = 0\f$ and the histogram of *all* molecules of the component. It samples
the bottom of the cage far better than any window does, and it fixes the
normalization \f$\int_\text{cage} P\,dq\f$ that turns the height of the
barrier into the probability \f$P(q^*)\f$. Without it, windows on the barrier
alone would miss the wells, and \f$P(q^*)\f$ would be far too high.

## The value at the barrier

\f$P(q^*)\f$ is taken as the value of a parabola fitted to \f$\ln P\f$ over a
short range around the dividing plane, not as the content of one bin. The fit
uses the neighbouring bins to reduce the noise, and gives the curvature
\f$\beta F''(q^*)\f$ as a by-product, which is negative at a barrier and
tells whether the window stiffness was sufficient.

## The run

After the initialization the run has two stages:

1. the **pilot**: unbiased Monte Carlo, in which the histogram of \f$q\f$ of
   every molecule of the tagged component is accumulated, together with the
   distribution of the occupancies of neighbouring cages;
2. the **windows**: biased Monte Carlo per window, the windows distributed
   over the replicas (threads). The tagged molecule is moved by biased
   translations and biased reinsertions within the support of its window, and
   a fixed fraction of the steps is spent on it, so that its histogram does
   not starve when there are many other molecules.

Mirror symmetry of the cage about the window can be used to sample only one
side of the barrier.

## Occupancy-resolved rates

The hopping rate of a molecule depends on how full its own cage and the
receiving cage are. The pilot records the joint distribution
\f$P(n_A, n_B)\f$ of the two occupancies, and the window stage records the
same distribution conditioned on the tagged bead being at the dividing plane.
The bias acts on \f$q\f$ only, so the second is unbiased in the occupancies,
and their ratio is the factor by which the transition-state rate of a
particular pair of occupancies differs from the average,

\f[
\frac{k^\text{TST}(n_A \to n_B)}{k^\text{TST}} = \frac{P(n_A, n_B\,|\,q^*)}{P(n_A, n_B)} .
\f]

These rates are the input of a lattice model of the diffusion; see the
[transmission driver](\ref driver_dctst).

## Things to check

The summary lists, per window, the mean and width of the sampled coordinate,
the overlap with the next window and the reweighting free energy. To look
for:

- a window whose mean lies far from its centre: it does not hold its
  position, and the stiffness or the spacing must change;
- a small overlap between two windows: the profile has an unconstrained step
  there;
- the deviation between the pilot histogram and the reweighted profile where
  the pilot is well sampled: a consistency check of the windows and the bias
  table;
- the bin estimate of \f$P(q^*)\f$ next to the fitted one;
- a warning that part of the cage holds no samples.

## Limits

- The framework must be rigid, and the tagged component rigid or fully
  flexible.
- Moves that change the number of molecules or the connectivity are refused.
- Cages that Monte Carlo can reach but dynamics cannot (the sodalite cages of
  LTA) must be blocked.
- Binary restart files are not supported.

## References

- C.H. Bennett, *Molecular dynamics and transition state theory: the
  simulation of infrequent events*, in *Algorithms for Chemical
  Computations*, ACS Symposium Series **46**, 63 (1977).
- D. Chandler, *Statistical mechanics of isomerization dynamics in liquids
  and the transition state approximation*, J. Chem. Phys. **68**, 2959 (1978).
- G.M. Torrie and J.P. Valleau, *Nonphysical sampling distributions in Monte
  Carlo free-energy estimation: Umbrella sampling*,
  J. Comput. Phys. **23**, 187 (1977).
- S. Kumar, J.M. Rosenberg, D. Bouzida, R.H. Swendsen and P.A. Kollman, *The
  weighted histogram analysis method for free-energy calculations on
  biomolecules. I. The method*, J. Comput. Chem. **13**, 1011 (1992).
- E. Beerdsen, B. Smit and D. Dubbeldam, *Molecular Simulation of Loading
  Dependent Slow Diffusion in Confined Systems*,
  Phys. Rev. Lett. **93**, 248301 (2004).
- D. Dubbeldam, E. Beerdsen, T.J.H. Vlugt and B. Smit, *Molecular simulation
  of loading-dependent diffusion in nanoporous materials using extended
  dynamically corrected transition state theory*,
  J. Chem. Phys. **122**, 224712 (2005).
