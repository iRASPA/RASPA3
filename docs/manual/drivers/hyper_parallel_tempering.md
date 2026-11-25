# Hyper-parallel tempering
\page driver_hyper_parallel_tempering Hyper-parallel tempering

`"SimulationType" : "HyperParallelTempering"`

Hyper-parallel tempering extends [parallel tempering](\ref driver_parallel_tempering)
from a ladder of temperatures to a grid of temperatures *and* chemical
potentials. Every replica is a grand-canonical simulation at its own state
point, so one run gives a set of adsorption isotherms, one point per replica,
and the exchanges help the replicas at difficult state points (low
temperature, high loading) to equilibrate.

## The grand-canonical weight

A replica at inverse temperature \f$\beta\f$ and fugacities \f$f_i\f$ samples
configurations with \f$N_i\f$ molecules of component \f$i\f$ with the weight

\f[
\pi(\mathbf{x}) \propto e^{-\beta U(\mathbf{x})}\,
\prod_i \frac{(\beta f_i V)^{N_i}}{N_i!} .
\f]

The input gives pressures, not fugacities. The driver converts every
(temperature, pressure) grid point with the Peng-Robinson equation of state,
\f$f_i = \phi_i\,y_i\,p\f$, with the fugacity coefficient \f$\phi_i\f$
computed at that state point. A fugacity coefficient given in the input is
ignored, since one value cannot hold over the whole grid.

## The swap

Exchanging the configurations of two replicas \f$A\f$ and \f$B\f$ changes the
product of their weights by

\f[
\frac{\pi(\text{new})}{\pi(\text{old})}
= e^{(\beta_B - \beta_A)(U_B - U_A)}
\prod_i \left[\frac{\beta_A f_{A,i}}{\beta_B f_{B,i}}\right]^{N_{B,i} - N_{A,i}} ,
\f]

where \f$U_X\f$ and \f$N_{X,i}\f$ belong to the configuration that is at
state point \f$X\f$ before the swap. The volume and the factorials are the
same in both replicas and cancel. The acceptance is the minimum of one and
this ratio.

The two factors have a simple meaning. The first is the temperature factor of
ordinary parallel tempering. The second compares the numbers of molecules: a
configuration with more molecules is favoured at the state point with the
higher activity \f$\beta f\f$. Two special cases follow directly:

- **Same pressure, different temperatures.** Both factors contribute, because
  \f$\beta f\f$ depends on the temperature.
- **Same temperature, different pressures.** The energy factor is one, and
  the acceptance is \f$(f_A/f_B)^{N_B - N_A}\f$: the swap is accepted easily
  when the two replicas hold similar numbers of molecules.

## Overlap

A swap is likely when the joint distributions of \f$(N, U)\f$ of the two
replicas overlap. Along an isotherm the loading changes fastest where the
isotherm is steep, on a logarithmic pressure scale: at the onset of
adsorption, at a pore-filling step, at capillary condensation. Those are the
places where the pressure ladder must be dense. A pressure ladder that is
uniform in \f$\ln p\f$ is the natural starting point, because the loading is
closer to linear in \f$\ln p\f$ than in \f$p\f$.

A first-order transition is the hard case. On one side of the step the
replicas are empty, on the other they are full, and their distributions do
not overlap at the same temperature. The grid offers a way around: at a
higher temperature the step is smeared out or gone, so a configuration can
go up in temperature, across in pressure, and down again. That detour is the
reason for tempering in two dimensions.

## The schedule

Each grid point runs in its own thread. Every
`"ParallelTemperingSwapEvery"` cycles the threads wait for each other. The
sweeps alternate between the temperature direction (neighbouring temperatures
at the same pressure) and the pressure direction (neighbouring pressures at
the same temperature), and within a direction between the two sets of
disjoint pairs, so that a configuration can reach every point of the grid.
The replicas keep their state points; the configurations move.

## Results

Every replica reports ordinary grand-canonical averages at its state point:
they are unbiased, because every replica samples its own distribution. The
loadings are assembled into one isotherm file per temperature, with the
fugacity, the absolute loading with its error in several units, and the
pressure. The files are rewritten during the run, so the convergence of the
isotherms can be followed.

The combined output holds the acceptance per pair, separately for the two
directions. A pair with a low acceptance is a gap in the grid.

## From here to WHAM

The replicas of this driver already visit the states of each other. The
[reweighted-histogram driver](\ref driver_reweighted_histogram) runs the same
simulation and, in addition, combines the \f$(N, U)\f$ samples of all
replicas into one density of states, which gives the isotherms *between* the
grid points as well.

## Limits

- The number of threads is the number of temperatures times the number of
  pressures.
- The swap requires rigid, whole-molecule replicas: fractional molecules,
  flexible components and reactions reject all swap attempts.
- Restart files are not supported by this driver.

## References

- Q. Yan and J.J. de Pablo, *Hyper-parallel tempering Monte Carlo: Application
  to the Lennard-Jones fluid and the restricted primitive model*,
  J. Chem. Phys. **111**, 9509 (1999).
- D.-Y. Peng and D.B. Robinson, *A New Two-Constant Equation of State*,
  Ind. Eng. Chem. Fundam. **15**, 59 (1976).
- D.J. Earl and M.W. Deem, *Parallel tempering: Theory, applications, and new
  perspectives*, Phys. Chem. Chem. Phys. **7**, 3910 (2005).
