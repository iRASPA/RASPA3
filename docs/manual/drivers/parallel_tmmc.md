# Parallel transition-matrix Monte Carlo
\page driver_parallel_tmmc Parallel transition-matrix Monte Carlo

`"SimulationType" : "ParallelTMMC"`

This driver computes the macrostate distribution \f$\Pi(N)\f$ of
[transition-matrix Monte Carlo](\ref driver_transition_matrix) with many
independent walkers, each confined to its own *window* of molecule numbers.
The theory of the collection matrix, the bias and the reweighting is that of
the serial method and is not repeated here; this page explains the windows
and what is done with the result.

## Why windows

A flat-histogram walk over a range of \f$\Delta N\f$ macrostates is a random
walk in \f$N\f$. The time to cross the range once grows as
\f$\Delta N^2\f$, and the estimate of \f$\Pi\f$ needs many crossings. Dividing
the range into \f$W\f$ windows of width \f$\Delta N/W\f$ reduces the crossing
time of each walker by \f$W^2\f$, and the walkers run at the same time. Even
on one processor the total work to cross all windows is a factor \f$W\f$
smaller.

The second gain is in the start-up. Each walker only has to flatten its own
window, which is a much smaller free-energy range than the whole.

## Stitching without overlap fitting

In histogram methods the windows are joined by matching the histograms in
their overlap, which adds an error per junction. Here that is not needed. A
walker records the unbiased acceptance probability of every attempted
insertion and deletion, *also of those that are rejected because they would
leave the window*. The collection matrix of a window therefore holds correct
transition statistics for all of its macrostates, including the two at its
edges. Neighbouring windows share their end-point macrostate, and since the
entries are sums over attempts, the collection matrices of all windows of
one temperature simply add,

\f[
C(N \to M) = \sum_{w} C_w(N \to M) .
\f]

The distribution over the full range then follows from the one recursion

\f[
\ln\Pi(N+1) = \ln\Pi(N) + \ln P(N \to N+1) - \ln P(N+1 \to N) .
\f]

There is no junction error, because there is no junction: every link of the
chain is estimated from the pooled statistics of whichever walkers sampled
it.

## The walkers

The macrostate range from the minimum to the maximum number of molecules is
split into `"NumberOfWindows"` contiguous windows, and one walker is made for
every (temperature, window) pair. Each walker

- starts inside its window: molecules are grown up to the lower window bound,
- runs grand-canonical Monte Carlo in its own thread, with moves that would
  leave the window rejected (after their acceptance probability has been
  recorded),
- is flattened by its own transition-matrix bias, re-derived from its
  collection matrix at regular intervals, with a Wang-Landau bias driving the
  first passes,
- does a fixed number of Monte Carlo moves per cycle, equal to the global
  maximum number of molecules, so that all walkers do the same work.

The walkers never exchange anything, so the threads only meet at the
boundaries of the stages. The collection matrix is not reset between the
equilibration and production stages, since its entries stay valid when the
bias changes; the final distribution uses all of it.

## Errors

The error bars of \f$\ln\Pi\f$, and of everything derived from it, come from
the blocks of the production stage: the increments of the collection matrices
per block give one estimate of \f$\ln\Pi\f$ each, and their spread is the
error.

## Isotherms

Each temperature is solved from its own walkers. The distribution is exact at
the reference fugacity (from the reference pressure through the
Peng-Robinson equation of state) and is reweighted to any other,

\f[
\ln\Pi(N;\,f) = \ln\Pi(N;\,f_\text{ref}) + N\ln\frac{f}{f_\text{ref}} .
\f]

On a logarithmic grid of pressures three isotherms are written:

- the **equilibrium isotherm**, \f$\langle N\rangle = \sum_N N\,\Pi(N;f)\f$
  over the whole range;
- the **adsorption branch**, the average over the low-density peak only,
  which is the metastable state followed when the pressure is raised;
- the **desorption branch**, the average over the high-density peak only,
  followed when the pressure is lowered.

Where \f$\Pi(N;f)\f$ has one peak the three coincide. Where it has two, the
branches separate and trace the hysteresis loop around the equilibrium step.
The loop of a real experiment is narrower than the two branches drawn this
way to the ends of their existence, because a real system leaves a metastable
state by nucleation before it becomes unstable.

Because temperatures are solved separately, potentials that depend on the
temperature are allowed, unlike in the
[reweighted-histogram driver](\ref driver_reweighted_histogram).

## Vapour-liquid equilibrium

For a fluid in a box the fugacity at which the vapour and liquid peaks of
\f$\Pi(N)\f$ have equal area is found by bisection. The peak positions give
the saturated densities. The saturation pressure follows from the
normalization: the empty box has
\f$\Pi(0) = 1/\Xi\f$ and \f$\beta p V = \ln\Xi\f$, so

\f[
\beta\,p\,V = -\ln\Pi(0;\,f) .
\f]

This is exact when the macrostate range starts at zero, which is the reason
to set the minimum to zero for such a calculation. The flat-histogram walk
crosses the gap between vapour and liquid by construction, so no special
starting configuration is needed, in contrast to grand-canonical sampling at
one state point.

## Things to check

- **Coverage.** The combined output lists the macrostates visited by each
  walker. Every state of every window must have been visited for the stitched
  distribution to be reliable.
- **The upper bound.** The maximum number of molecules must comfortably
  exceed the liquid peak (or the saturation loading); otherwise the isotherms
  saturate artificially.
- **The number of windows.** More windows shorten the equilibration, but each
  window must still be crossed many times, and each is a thread.

## Limits

- Exactly one system with exactly one component.
- The number of threads is the number of temperatures times the number of
  windows.
- Restart files are not supported by this driver.

## References

- J.R. Errington, *Direct calculation of liquid-vapor phase equilibria from
  transition matrix Monte Carlo simulation*, J. Chem. Phys. **118**, 9915
  (2003).
- V.K. Shen and J.R. Errington, *Determination of fluid-phase behavior using
  transition-matrix Monte Carlo: Binary Lennard-Jones mixtures*,
  J. Chem. Phys. **122**, 064508 (2005).
- D.W. Siderius and V.K. Shen, *Use of the Grand Canonical Transition-Matrix
  Monte Carlo Method to Model Gas Adsorption in Porous Materials*,
  J. Phys. Chem. C **117**, 5861 (2013).
- N.B. Wilding, *Critical-point and coexistence-curve properties of the
  Lennard-Jones fluid: A finite-size scaling study*,
  Phys. Rev. E **52**, 602 (1995).
