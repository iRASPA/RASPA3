# dcTST: the transmission coefficient
\page driver_dctst dcTST: the transmission coefficient

`"SimulationType" : "DynamicallyCorrectedTST"`

This driver computes the second half of dynamically corrected transition
state theory: the transmission coefficient \f$\kappa\f$, the fraction of
barrier crossings that lead to a hop. Combined with the transition-state rate
of the [free-energy driver](\ref driver_dctst_free_energy) it gives the
hopping rate and the self and collective diffusion coefficients,

\f[
k_{A\to B} = \kappa_S\,k^\text{TST}, \qquad
D_S = \kappa_S\,k^\text{TST}\,\lambda^2, \qquad
D_C = \kappa_C\,k^\text{TST}\,\lambda^2 .
\f]

The background (the lattice picture of hopping, the reaction coordinate
\f$q\f$, the dividing plane \f$q^*\f$ and \f$k^\text{TST}\f$) is on the page
of the free-energy driver and is assumed here.

## Why transition state theory overestimates

Transition state theory counts every crossing of the dividing plane towards
\f$B\f$ as a hop. In reality a molecule that crosses may turn back: it
collides with the window, with a molecule in the receiving cage, or it simply
has not lost the energy that carried it over. Each real hop is then counted
several times, once per crossing. The exact rate is the transition-state rate
multiplied by the fraction of crossings that are *final*.

## The reactive flux

Bennett and Chandler expressed that fraction as a correlation function. Start
trajectories from the equilibrium distribution *constrained to the dividing
plane*, with velocities from the Maxwell-Boltzmann distribution, and ask
where they are a time \f$t\f$ later:

\f[
\kappa(t) = \frac{\left\langle \dot q(0)\;\theta\!\left[q(t) - q^*\right]\right\rangle_{q^*}}
{\left\langle \tfrac{1}{2}\,|\dot q(0)|\right\rangle_{q^*}} ,
\f]

with \f$\theta\f$ the step function. The numerator weights each trajectory
with its initial velocity: a trajectory that starts towards \f$B\f$ and stays
there contributes positively, one that starts towards \f$A\f$ but ends in
\f$B\f$ contributes negatively, and recrossings cancel.

At \f$t = 0^+\f$ every trajectory is still on the side it started towards,
and \f$\kappa = 1\f$: transition state theory. As the recrossings happen
\f$\kappa(t)\f$ decays, and once every molecule has committed to a cage it
reaches a *plateau*. The plateau value is the transmission coefficient. It
holds for times longer than the recrossing dynamics but much shorter than the
time between hops, and the existence of the plateau is the sign that the
separation of time scales, on which the whole rate picture rests, is real.

Each starting configuration is used for a pair of trajectories with opposite
velocities. The time-reversed partner costs one more trajectory and removes
most of the noise of the velocity weighting.

## The constrained ensemble

The starting points must be distributed as the equilibrium ensemble with the
reaction bead on the plane. The driver samples that ensemble with Monte
Carlo: the tagged molecule has its reaction bead held on the barrier plane
and is moved by in-plane translations and in-plane reinsertions, while all
other molecules are sampled with their ordinary moves. The starting
configuration is built with the tagged molecule first, on the plane in the
empty framework, and the others grown around it; placing it into a full
system would fail at loadings where every window is occupied.

At regular intervals of this Monte Carlo run a number of trajectory pairs is
launched from the current configuration and integrated for a short time (of
the order of picoseconds).

## Self and collective transmission

**Self.** \f$\kappa_S\f$ follows the tagged molecule: the cage it is in at
time \f$t\f$, counted along the barrier normal. It gives the hopping rate of
a molecule and \f$D_S\f$.

**Collective.** The collective (corrected, Maxwell-Stefan) diffusion
coefficient describes the motion of the centre of mass of all molecules, the
quantity that determines transport under a gradient. \f$\kappa_C\f$ is
accumulated from the *net cage displacement of all molecules* during the same
trajectories: when the hop of the tagged molecule pushes another molecule
onward, the collective displacement is larger than the self displacement;
when a molecule moves back into the cage that was left, it is smaller. The
contributions are resolved per shell of cages around the hop and per
component.

In a cage system part of the correlation between hops acts on the time scale
of the hops themselves: a hop leaves a vacancy that is refilled much later
than a picosecond trajectory lasts. The plateau of \f$\kappa_C(t)\f$ captures
the fast, concerted part only. The modes below address the rest.

## Estimators of \f$\kappa_S\f$

Several estimators are formed from the same trajectories. They agree when
molecules thermalize in the cage, and their disagreement is diagnostic.

- **Net cage count**: the cage index of the tagged molecule, as above.
- **Two-state (Bennett-Chandler)**: only the side of the dividing plane
  counts.
- **Committed**: the side is frozen once the molecule has moved a given
  distance from the plane (at most to the cage centre).
- **Committor-smoothed** (Ruiz-Montero, Frenkel and Brey): the step function
  is replaced by the committor of the free-energy profile,

\f[
p_B(q) = \frac{\displaystyle\int_{-\lambda/2}^{q} e^{\beta F(q')}\,dq'}
{\displaystyle\int_{-\lambda/2}^{\lambda/2} e^{\beta F(q')}\,dq'} ,
\f]

  the probability that a molecule at \f$q\f$ ends in \f$B\f$ in the diffusive
  limit. It is flat inside the cages and switches across the window. A
  trajectory that lingers near the plane then contributes a value between 0
  and 1 instead of jumping between them, which removes the variance of the
  decision. This is the recommended value when a profile is available.

## Thermalization: the low-loading problem

In a rigid framework at low loading nothing takes the energy from a molecule
that has crossed a window. It runs on to a second window, or bounces back
through the first, and \f$\kappa(t)\f$ never reaches a plateau. This is
Kramers' energy-diffusion regime, in which the rate is limited by the
exchange of energy, and it is an artefact of the rigid lattice: in a real
material the lattice vibrations take the energy.

A stochastic bath on all molecules therefore stands in for the lattice during
the trajectories:

- **Andersen** (default): every molecule gets a fresh Maxwell-Boltzmann
  velocity with a given rate; between collisions the dynamics is exactly
  Hamiltonian.
- **Langevin**: friction and noise on the translational and angular
  velocities.
- **None**: plain \f$NVE\f$, appropriate when the other molecules thermalize
  a molecule within the cage (high loading), and as a diagnostic.

The bath rate must lie between two time scales: slower than the recrossing
dynamics in the window, so that the barrier crossing is unperturbed, and
faster than the residence time in a cage. Between them, a range of several
decades, \f$\kappa\f$ does not depend on the rate; this should be verified
with two or three values. The *long-jump fraction* (trajectories in which the
tagged molecule reached a second window) reports the problem directly and
should be close to zero.

## Occupancy-resolved rates

A hop depends on the occupancies \f$n_A\f$ and \f$n_B\f$ of the cage a
molecule leaves and the cage it enters. Every trajectory is classified by the
occupancies at launch, which gives \f$\kappa(n_A, n_B)\f$ per class. With the
transition-state rates per class from the free-energy driver,

\f[
k(n_A \to n_B) = \kappa(n_A, n_B)\; k^\text{TST}\;\frac{P(n_A, n_B\,|\,q^*)}{P(n_A, n_B)} ,
\f]

and the consistency check is that the weighted sum over the classes returns
the average rate, \f$\sum P(n_A,n_B)\,k(n_A\to n_B) = \kappa\,k^\text{TST}\f$.
The set of rates is the input of a lattice model (a master equation or
kinetic Monte Carlo) that yields \f$D_S\f$ and \f$D_C\f$ *including* the slow,
vacancy-mediated correlations. A run can also be conditioned on one pair of
occupancies, held by rejecting the moves that would change them.

## Event-terminated trajectories

Instead of integrating every trajectory for a fixed time, a trajectory can be
stopped at the \f$n\f$-th *local event* after the tagged molecule has
committed: a hop of any molecule through one of the three windows of the
column of cages \f$A', A, B, B'\f$ around the tagged window. The values at
termination define two further estimators, \f$\kappa_S^\text{event}\f$ and
\f$\kappa_C^\text{event}\f$, and two numbers that summarize the correlations:

- the **back-correlation factor**
  \f$f = \kappa_S^\text{event}/\kappa_S^\text{committed}\f$: the recrossings of
  the tagged molecule after commitment but before any other local hop;
- the **collective excess**
  \f$\epsilon = \kappa_C^\text{event}/\kappa_S^\text{committed} - 1\f$: the net
  correlated flux of the other molecules up to that event.

Stopping at the first event is a first-order truncation of the collective
correlations; raising \f$n\f$ includes the next orders, and \f$f\f$ and
\f$\epsilon\f$ should not change when converged. A run at order \f$n\f$
contains all lower orders, so one run gives the convergence table.

## The hybrid estimator

The first-event estimator samples a race between hop channels. Every
trajectory ends with a sign, so \f$\epsilon\f$ has a variance of order one
per trajectory, and needs of the order of \f$1/\epsilon^2\f$ trajectories,
each as long as the local exchange time. The hybrid estimator avoids sampling
the race. At a handover time \f$\tau\f$ after commitment it splits the
collective excess in two,

\f[
\epsilon_\text{hyb}(\tau) = \epsilon_\text{conc}(\tau) + \epsilon_\text{Markov}(\tau) .
\f]

The *concerted* part is measured: the displacement up to \f$\tau\f$ of the
tagged molecule and of the others in the column, with the committor in place
of the step function. The *Markov* part is computed: from the occupancies of
the column at \f$\tau\f$ and the occupancy-resolved rates, the expected sign
of the next event,

\f[
\epsilon_\text{Markov}(\tau) = \sum_{n} P(n)\;\frac{\sum_c s_c\,k_c(n)}{\sum_c k_c(n)} ,
\f]

over the hop channels \f$c\f$ with their signs \f$s_c\f$. If the hopping
beyond \f$\tau\f$ is Markovian the sum does not depend on \f$\tau\f$; a drift
with \f$\tau\f$ measures the non-Markovian part. The variance is now set by
the molecules that actually cross a window before \f$\tau\f$, which is of
order \f$\epsilon\f$ rather than one, and the cost per trajectory is the
recrossing time plus \f$\tau\f$. That makes the estimator affordable at low
loading, where the exchange time is long.

## Things to check

- **The plateau.** \f$\kappa(t)\f$ must be flat over the window used for the
  average.
- **The estimators** of \f$\kappa_S\f$ should agree; the long-jump fraction
  should be near zero.
- **The mean speed** at the plane, compared with
  \f$\sqrt{2k_BT/\pi m}\f$, and the energy drift of the trajectories.
- **The cage centre** printed in the summary, which anchors the cage lattice.

## Recovering and pooling runs

Binary restart files are not supported. Instead, a termination signal stops
the replicas at the next cycle boundary and the analysis is written over the
completed trajectory pairs; one line per completed pair is written to a file,
from which every estimator can be recomputed, so that runs can be
concatenated and pooled exactly; and the accumulators are saved at regular
intervals.

## References

- C.H. Bennett, *Molecular dynamics and transition state theory: the
  simulation of infrequent events*, in *Algorithms for Chemical
  Computations*, ACS Symposium Series **46**, 63 (1977).
- D. Chandler, *Statistical mechanics of isomerization dynamics in liquids
  and the transition state approximation*, J. Chem. Phys. **68**, 2959 (1978).
- M.J. Ruiz-Montero, D. Frenkel and J.J. Brey, *Efficient schemes to compute
  diffusive barrier crossing rates*, Mol. Phys. **90**, 925 (1997).
- H.A. Kramers, *Brownian motion in a field of force and the diffusion model
  of chemical reactions*, Physica **7**, 284 (1940).
- E. Beerdsen, B. Smit and D. Dubbeldam, *Molecular Simulation of Loading
  Dependent Slow Diffusion in Confined Systems*,
  Phys. Rev. Lett. **93**, 248301 (2004).
- D. Dubbeldam, E. Beerdsen, T.J.H. Vlugt and B. Smit, *Molecular simulation
  of loading-dependent diffusion in nanoporous materials using extended
  dynamically corrected transition state theory*,
  J. Chem. Phys. **122**, 224712 (2005).
- E. Beerdsen, D. Dubbeldam and B. Smit, *Understanding Diffusion in
  Nanoporous Materials*, Phys. Rev. Lett. **96**, 044501 (2006).
