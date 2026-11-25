# Parallel tempering
\page driver_parallel_tempering Parallel tempering

`"SimulationType" : "ParallelTempering"`

Parallel tempering (replica exchange) runs several copies of one system at
different temperatures and lets them exchange configurations. It is used when
a single Monte Carlo run gets trapped: a configuration that cannot cross a
barrier at a low temperature can travel to a high temperature, cross there,
and come back.

## The extended ensemble

Take \f$M\f$ replicas at inverse temperatures
\f$\beta_1 > \beta_2 > \dots > \beta_M\f$. They do not interact, so the
probability of the set of configurations is the product of the individual
weights,

\f[
\pi(\mathbf{x}_1,\dots,\mathbf{x}_M) \propto \prod_{k=1}^{M} e^{-\beta_k U(\mathbf{x}_k)} .
\f]

Two kinds of moves sample this distribution. Within a replica the ordinary
[Monte Carlo](\ref driver_monte_carlo) moves act. Between two replicas a
*swap* exchanges their configurations. A swap changes the weight by

\f[
\frac{\pi(\text{new})}{\pi(\text{old})}
= \frac{e^{-\beta_A U_B}\, e^{-\beta_B U_A}}{e^{-\beta_A U_A}\, e^{-\beta_B U_B}}
= e^{(\beta_B - \beta_A)(U_B - U_A)} ,
\f]

and the proposal is symmetric, so the acceptance rule is

\f[
\text{acc} = \min\left(1,\; e^{(\beta_B - \beta_A)(U_B - U_A)}\right) .
\f]

No energy has to be computed for a swap: the two energies are already known.
Because every replica keeps sampling its own Boltzmann distribution, the
averages of each replica are ordinary canonical averages at its temperature;
the swaps only help the replicas reach equilibrium.

When the replicas are open systems the weights contain the fugacity, and the
swap acceptance gets the factor of Yan and de Pablo,

\f[
\text{acc} = \min\left(1,\; e^{(\beta_B - \beta_A)(U_B - U_A)}
\prod_i \left[\frac{\beta_A f_{A,i}}{\beta_B f_{B,i}}\right]^{N_{B,i} - N_{A,i}}\right) ,
\f]

with \f$N_{X,i}\f$ the number of molecules of component \f$i\f$ in the
configuration that is at state point \f$X\f$ before the swap. A
\f$\beta\,p\,\Delta V\f$ term enters in the same way when the boxes travel
with the configurations.

## Why neighbouring temperatures

The acceptance is governed by the overlap of the energy distributions of the
two replicas. The energy of a replica has a mean that grows with the
temperature and a width \f$\sigma_U = \sqrt{k_B T^2 C_V}\f$. Two temperatures
that are far apart have distributions that do not overlap, and a swap between
them is essentially never accepted. Swaps are therefore attempted between
neighbours of the ladder only.

Since \f$C_V\f$ is proportional to the number of degrees of freedom
\f$N_f\f$, the relative width of the energy distribution shrinks as
\f$1/\sqrt{N_f}\f$: the number of replicas needed to span a temperature range
grows as \f$\sqrt{N_f}\f$. For a constant heat capacity, an equal acceptance
between all neighbours is obtained with a *geometric* ladder,
\f$T_{k+1}/T_k = \text{constant}\f$. Near a phase transition the heat capacity
peaks and the ladder must be denser there.

## The swap schedule

The replicas run in their own threads. Every `"ParallelTemperingSwapEvery"`
cycles they wait for each other, and swaps are attempted for a set of
disjoint neighbour pairs: alternately the pairs \f$(1,2),(3,4),\dots\f$ and
the pairs \f$(2,3),(4,5),\dots\f$. The alternation is what allows a
configuration to move through the whole ladder. The replicas keep their
temperatures; the configurations migrate.

## Measuring whether it works

The acceptance of the swaps is a local measure: it says whether two
neighbours exchange, not whether configurations travel. The quantity that
matters is the *round trip*: a configuration that goes from the coldest
replica to the hottest and back. The driver reports

- the number of round trips and the mean round-trip time,
- the round trips of every individual configuration,
- the *up-moving fraction* \f$f(T)\f$: of the configurations that visit a
  temperature, the fraction that visited the coldest replica more recently
  than the hottest.

For a ladder in which the configurations diffuse freely, \f$f\f$ falls
linearly from 1 at the coldest replica to 0 at the hottest. A plateau in
\f$f(T)\f$ marks a bottleneck: configurations arrive from both sides and turn
back. The remedy is more replicas at that place, which is not necessarily
where the swap acceptance is lowest.

## Limits

- A swap exchanges whole configurations between replicas with the same
  Hamiltonian and the same topology. Systems with fractional molecules,
  flexible components or reactions reject all swap attempts in this driver.
- Each temperature is a thread. The thread pool used inside one energy
  evaluation should stay serial, or the machine is oversubscribed.
- Restart files are not supported by this driver.

## References

- R.H. Swendsen and J.-S. Wang, *Replica Monte Carlo Simulation of Spin
  Glasses*, Phys. Rev. Lett. **57**, 2607 (1986).
- C.J. Geyer, *Markov chain Monte Carlo maximum likelihood*, in *Computing
  Science and Statistics: Proceedings of the 23rd Symposium on the Interface*,
  156 (1991).
- K. Hukushima and K. Nemoto, *Exchange Monte Carlo Method and Application to
  Spin Glass Simulations*, J. Phys. Soc. Jpn. **65**, 1604 (1996).
- Q. Yan and J.J. de Pablo, *Hyper-parallel tempering Monte Carlo: Application
  to the Lennard-Jones fluid and the restricted primitive model*,
  J. Chem. Phys. **111**, 9509 (1999).
- D.J. Earl and M.W. Deem, *Parallel tempering: Theory, applications, and new
  perspectives*, Phys. Chem. Chem. Phys. **7**, 3910 (2005).
- H.G. Katzgraber, S. Trebst, D.A. Huse and M. Troyer, *Feedback-optimized
  parallel tempering Monte Carlo*, J. Stat. Mech. P03018 (2006).
