# Histogram reweighting (WHAM)
\page driver_reweighted_histogram Histogram reweighting (WHAM)

`"SimulationType" : "ReweightedHistogram"`

This driver runs the replica-exchange simulation of
[hyper-parallel tempering](\ref driver_hyper_parallel_tempering) over a grid
of temperatures and pressures, and then combines the samples of all replicas
with the *weighted histogram analysis method* (WHAM). The result is not a set
of isotherm points but a continuous surface: the loading, and other averages,
at any temperature and pressure inside (and somewhat beyond) the simulated
grid. For a fluid in a box the same analysis gives the vapour-liquid
coexistence curve.

## The idea: one simulation, many state points

A grand-canonical simulation at inverse temperature \f$\beta\f$ and fugacity
\f$f\f$ visits a configuration with \f$N\f$ molecules and energy \f$U\f$ with
probability

\f[
P(N, U;\,\beta, f) = \frac{\Omega(N, U)\; e^{-\beta U}\,(\beta f)^{N}}{\Xi(\beta, f)} ,
\f]

where \f$\Omega(N,U)\f$ is the *density of states*: the number of
configurations with that \f$N\f$ and \f$U\f$ (the factors \f$V^N/N!\f$ are
absorbed in it). The density of states does not depend on the temperature or
the fugacity. It is a property of the system alone; the state point only
enters through the known factor \f$e^{-\beta U}(\beta f)^N\f$.

So a histogram \f$H(N,U)\f$ measured at one state point gives an estimate of
\f$\Omega\f$, up to a constant,

\f[
\Omega(N, U) \propto H(N, U)\; e^{+\beta U}\,(\beta f)^{-N} ,
\f]

and with \f$\Omega\f$ the distribution at *another* state point
\f$(\beta', f')\f$ follows by putting the new factor back,

\f[
P(N, U;\,\beta', f') \propto H(N, U)\; e^{-(\beta' - \beta) U}\,
\left(\frac{\beta' f'}{\beta f}\right)^{N} .
\f]

This is single-histogram reweighting (Ferrenberg and Swendsen, 1988). It
works only close to the simulated state point. Far from it the new
distribution has its weight where the histogram has few or no counts: the
estimate is then dominated by the noise in the tail of the histogram.

## Combining several simulations

Simulations at several state points cover different regions of \f$(N,U)\f$.
Each gives its own estimate of \f$\Omega\f$, good where its histogram is well
filled. The question is how to combine them, when each estimate has an
unknown normalization \f$\Xi_i\f$.

Let run \f$i\f$ (at \f$\beta_i, f_i\f$) have \f$n_i\f$ samples and histogram
\f$H_i(N,U)\f$. Write \f$g_i = -\ln \Xi_i\f$ for its dimensionless free
energy. The expected count in a bin is

\f[
\langle H_i(N,U)\rangle = n_i\;\Omega(N,U)\; e^{\,g_i - \beta_i U + N\ln(\beta_i f_i)} .
\f]

The counts are (nearly) Poisson distributed, so the relative error of the
estimate of \f$\Omega\f$ from run \f$i\f$ is \f$1/\sqrt{\langle H_i\rangle}\f$.
The combination of the estimates with the smallest variance weights each run
in proportion to its expected count. Carrying that out gives the WHAM
equations,

\f[
\Omega(N,U) = \frac{\displaystyle\sum_i H_i(N,U)}
{\displaystyle\sum_i n_i\; e^{\,g_i - \beta_i U + N\ln(\beta_i f_i)}} ,
\f]

\f[
e^{-g_i} = \sum_{N,U} \Omega(N,U)\; e^{-\beta_i U + N\ln(\beta_i f_i)} .
\f]

The first equation has a transparent reading: the numerator is the total
number of samples that fell in the bin, from all runs together, and the
denominator is the number that *should* have fallen there, per unit of
\f$\Omega\f$, given where the runs were. The second equation is just the
normalization of each run. The two depend on each other and are solved
together.

The same equations follow from maximizing the likelihood of the observed
histograms, which shows that WHAM is the statistically optimal estimate
given binned data. In the limit of infinitely narrow bins it becomes the
multistate Bennett acceptance ratio (MBAR).

## Solving the equations

The equations are solved by iteration: start from a guess for the \f$g_i\f$, compute
\f$\Omega\f$ from the first equation, compute new \f$g_i\f$ from the second,
and repeat until the largest change of any \f$g_i\f$ is below a tolerance.
Everything is done with logarithms,

\f[
\ln\Omega(N,U) = \ln M(N,U) - \ln\sum_i n_i\; e^{\,g_i - \beta_i U + N\ln(\beta_i f_i)} ,
\f]

with \f$M = \sum_i H_i\f$, because the terms span hundreds of orders of
magnitude; the sums are log-sum-exp operations. Only differences of the
\f$g_i\f$ are determined, since \f$\Omega\f$ has an arbitrary overall
constant. The plain iteration converges slowly when the runs overlap poorly,
because information then travels only one grid point per iteration; the
tolerance and the maximum number of iterations are input keys.

## What is sampled

During production every replica records, every `"SampleReweightingEvery"`
cycles, one pair \f$(N, U)\f$: the number of molecules of the adsorbate and
the potential energy. At the end the samples of all replicas are pooled and
binned. The replica-exchange swaps of the simulation do not disturb this:
each replica samples its own distribution at all times, which is all WHAM
assumes.

A cycle of this driver is a fixed number of Monte Carlo moves, equal to the
maximum number of molecules, not the current number. Empty and full replicas
then do the same work per cycle and contribute samples at the same rate.

## Reweighting to other state points

With \f$\Omega\f$ known, any average that depends on \f$N\f$ and \f$U\f$ is
available at any temperature and fugacity,

\f[
\langle A\rangle_{\beta, f} =
\frac{\sum_{N,U} A(N,U)\;\Omega(N,U)\; e^{-\beta U}(\beta f)^N}
{\sum_{N,U} \Omega(N,U)\; e^{-\beta U}(\beta f)^N} .
\f]

Taking \f$A = N\f$ on a fine grid of pressures gives a smooth isotherm at
each requested temperature (each pressure is converted to a fugacity with the
Peng-Robinson equation of state). The fluctuations give the isosteric heat by
the same formula as in a single simulation, now as a continuous function of
loading.

**How far can one reweight?** The answer is in the samples. The reweighted
average is a weighted sum over the pooled samples, and when a few samples
carry most of the weight the result is unreliable. The *effective sample
size*,

\f[
n_\text{eff} = \frac{\left(\sum_s w_s\right)^2}{\sum_s w_s^2} ,
\f]

with \f$w_s\f$ the weight of sample \f$s\f$ at the target state point, is
reported with every reweighted point. It is of the order of the number of
samples inside the simulated grid and collapses outside the sampled region of
\f$(N,U)\f$. A second check is the *self-consistency table*: the loading
reweighted to each simulated grid point next to the loading measured directly
there. They must agree within the errors.

## Errors

The WHAM equations are solved again for each block of the production run
separately, and the spread of the results over the blocks gives the error of
every reweighted quantity. This captures both the statistical noise of the
histograms and its propagation through the \f$g_i\f$, which a formula based on
the counts alone would not.

## Vapour-liquid equilibrium

For a fluid in a box (no framework) the distribution of the number of
molecules,

\f[
P(N;\,\beta, f) \propto (\beta f)^N \sum_U \Omega(N,U)\, e^{-\beta U} ,
\f]

is bimodal below the critical temperature near the saturation pressure: one
peak at the vapour density, one at the liquid density. At coexistence the two
phases have the same pressure, and since \f$\beta p V = \ln\Xi\f$, each phase
contributes the same amount to the partition function: the two peaks have
*equal weight* (equal area, not equal height). The driver finds the
coexistence fugacity by bisection on the difference of the two areas. The
positions of the peaks give the saturated densities, and the saturation
pressure follows from \f$\beta p V = \ln\Xi\f$ with \f$\Xi\f$ normalized by
the empty box, \f$P(N{=}0) = 1/\Xi\f$.

Two practical points follow from the theory:

- The histograms must connect the two phases. That happens near the critical
  point, where the vapour and liquid peaks merge. The top of the temperature
  ladder should therefore be close to the critical temperature; the lower
  temperatures are then reached by reweighting along the ladder.
- Grand-canonical sampling does not nucleate the liquid from the vapour
  within a simulation. The run should start from a liquid-like configuration,
  so that both phases are visited.

## Requirements and limits

- **Overlap.** The \f$(N,U)\f$ histograms of neighbouring grid points must
  overlap. This is the same condition as a healthy swap acceptance, so the
  per-pair acceptance table of the replica exchange doubles as an overlap
  diagnostic.
- **One adsorbate.** The histograms are collected over the number of
  molecules of one component. Cations of a fixed number may be present.
- **Temperature-independent potentials.** The method reuses each sampled
  energy at all temperatures, so potentials that depend on the temperature
  (Feynman-Hibbs) are refused.
- **Bin width.** Wide bins in \f$U\f$ bias the reweighting in temperature,
  because the energy is treated as constant within a bin; narrow bins have
  few counts but this does not hurt WHAM, which sums the counts.
- Transition-matrix biasing, reactions and identity changes are refused.

## Comparison with transition-matrix Monte Carlo

[Transition-matrix Monte Carlo](\ref driver_parallel_tmmc) also yields the
distribution \f$P(N)\f$ and reweights it in fugacity. It obtains the
distribution from acceptance probabilities instead of visit counts and
forces the simulation to visit all \f$N\f$, which makes it the stronger
method across a first-order transition at one temperature. WHAM keeps the
energy as a second variable and therefore reweights in *temperature* as well,
from one run.

## References

- A.M. Ferrenberg and R.H. Swendsen, *New Monte Carlo technique for studying
  phase transitions*, Phys. Rev. Lett. **61**, 2635 (1988).
- A.M. Ferrenberg and R.H. Swendsen, *Optimized Monte Carlo data analysis*,
  Phys. Rev. Lett. **63**, 1195 (1989).
- S. Kumar, J.M. Rosenberg, D. Bouzida, R.H. Swendsen and P.A. Kollman, *The
  weighted histogram analysis method for free-energy calculations on
  biomolecules. I. The method*, J. Comput. Chem. **13**, 1011 (1992).
- M.R. Shirts and J.D. Chodera, *Statistically optimal analysis of samples
  from multiple equilibrium states*, J. Chem. Phys. **129**, 124105 (2008).
- N.B. Wilding, *Critical-point and coexistence-curve properties of the
  Lennard-Jones fluid: A finite-size scaling study*,
  Phys. Rev. E **52**, 602 (1995).
- A.Z. Panagiotopoulos, *Monte Carlo methods for phase equilibria of fluids*,
  J. Phys.: Condens. Matter **12**, R25 (2000).
- Q. Yan and J.J. de Pablo, *Hyper-parallel tempering Monte Carlo: Application
  to the Lennard-Jones fluid and the restricted primitive model*,
  J. Chem. Phys. **111**, 9509 (1999).
