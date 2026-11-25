# Drivers
\page drivers Drivers

A *driver* is the top-level routine that runs a simulation: it owns the
systems, decides what a cycle is, runs the stages of the simulation
(initialization, equilibration, production), and writes the results. The
driver is chosen with the top-level input key `"SimulationType"`; the keys
that each driver reads are documented on the [Commands](\ref commands) page.

The pages below explain the *theory* behind each driver: what it samples or
solves, the equations it is based on, how the estimators and their errors are
formed, and what to watch for when using it. They are meant to be read next
to the input reference, not instead of it.

## Sampling at one state point

| Driver | `"SimulationType"` | What it does |
|---|---|---|
| \subpage driver_monte_carlo | `"MonteCarlo"` | Markov-chain Monte Carlo in the ensemble defined by the moves that are switched on |
| \subpage driver_molecular_dynamics | `"MolecularDynamics"` | Integration of the equations of motion with thermostats, barostats and constraints |
| \subpage driver_molecular_dynamics_spatial_decomposition | `"MolecularDynamicsSpatialDecomposition"` | Molecular dynamics of large boxes with domain decomposition and a particle-mesh Ewald sum |
| \subpage driver_minimization | `"Minimization"` | Energy minimization by eigenvector following, with elastic constants, normal modes and phonons |

## Sampling many state points at once

| Driver | `"SimulationType"` | What it does |
|---|---|---|
| \subpage driver_parallel_tempering | `"ParallelTempering"` | Replica-exchange Monte Carlo over a ladder of temperatures |
| \subpage driver_parallel_tempering_molecular_dynamics | `"ParallelTemperingMolecularDynamics"` | Replica-exchange molecular dynamics over a ladder of temperatures |
| \subpage driver_hyper_parallel_tempering | `"HyperParallelTempering"` | Replica-exchange Monte Carlo over a grid of temperatures and pressures |
| \subpage driver_reweighted_histogram | `"ReweightedHistogram"` | Hyper-parallel tempering combined with multiple-histogram reweighting (WHAM) |

## Flat-histogram sampling

| Driver | `"SimulationType"` | What it does |
|---|---|---|
| \subpage driver_transition_matrix | `"MonteCarloTransitionMatrix"` | Transition-matrix Monte Carlo over the number of molecules |
| \subpage driver_parallel_tmmc | `"ParallelTMMC"` | Transition-matrix Monte Carlo with one walker per window of molecule numbers |

## Free energies

| Driver | `"SimulationType"` | What it does |
|---|---|---|
| \subpage driver_thermodynamic_integration | `"ThermodynamicIntegration"` | The average of \f$\partial U/\partial\lambda\f$ at one fixed coupling parameter |
| \subpage driver_parallel_thermodynamic_integration | `"ParallelThermodynamicIntegration"` | The whole \f$\langle\partial U/\partial\lambda\rangle\f$ curve in one run, with exchanges of \f$\lambda\f$ between replicas |

## Rare events and diffusion

| Driver | `"SimulationType"` | What it does |
|---|---|---|
| \subpage driver_dctst_free_energy | `"DynamicallyCorrectedTSTFreeEnergy"` | The free-energy profile over a barrier by umbrella sampling, and the transition-state rate |
| \subpage driver_dctst | `"DynamicallyCorrectedTST"` | The transmission coefficient from trajectories started on the barrier, and the diffusion coefficients |

## Other methods

| Driver | `"SimulationType"` | What it does |
|---|---|---|
| \subpage driver_nldft | `"NLDFT"` | An adsorption isotherm from classical density functional theory on the framework's own energy grid |
| \subpage driver_polymerization | `"Polymerization"` | Reactive molecular dynamics in which bonds form between monomers according to reaction rules |

## How the drivers relate

Several drivers are built on others, and reading the pages in this order
avoids repetition:

- [Monte Carlo](\ref driver_monte_carlo) and
  [molecular dynamics](\ref driver_molecular_dynamics) are the two basic
  engines. Every other sampling driver runs one of them inside each replica,
  walker or window.
- [Parallel tempering](\ref driver_parallel_tempering) adds exchanges between
  replicas at different temperatures.
  [Hyper-parallel tempering](\ref driver_hyper_parallel_tempering) extends the
  ladder to a grid of temperatures and pressures, and the
  [reweighted-histogram driver](\ref driver_reweighted_histogram) analyses the
  samples of that grid with WHAM.
- [Transition-matrix Monte Carlo](\ref driver_transition_matrix) obtains the
  same kind of result as WHAM, the probability of every number of molecules,
  from the acceptance probabilities instead of from histograms;
  [parallel TMMC](\ref driver_parallel_tmmc) divides the range of molecule
  numbers over threads.
- [Thermodynamic integration](\ref driver_thermodynamic_integration) and its
  [parallel version](\ref driver_parallel_thermodynamic_integration) use the
  fractional molecule of the CFCMC moves, held at a fixed coupling.
- The two dcTST drivers form one method: the
  [free-energy driver](\ref driver_dctst_free_energy) gives the probability of
  being at the barrier and the
  [transmission driver](\ref driver_dctst) the fraction of barrier crossings
  that succeed.
- [Polymerization](\ref driver_polymerization) is built on
  [spatial-decomposition molecular dynamics](\ref driver_molecular_dynamics_spatial_decomposition).

## Notation

The pages share the notation of the [Monte Carlo moves](\ref mc_moves) page:
\f$\beta = 1/(k_B T)\f$, \f$U\f$ the potential energy, \f$N\f$ the number of
molecules, \f$V\f$ the volume, \f$p\f$ the pressure, \f$f\f$ the fugacity,
\f$\mu\f$ the chemical potential and \f$\lambda\f$ the coupling parameter of a
fractional molecule. Angular brackets \f$\langle\cdot\rangle\f$ are ensemble
averages.
