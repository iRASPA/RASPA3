# Simulation input options
\page commands Commands

`RASPA` is driven by a single JSON input file, `simulation.json`. This file
selects the type of simulation, sets the global run length, and describes one
or more *systems* and the *components* (molecules) that live in them. This page
documents the available keywords. Keyword names are matched case-insensitively.

## Table of Contents

<!-- TOC -->
* [Input sections](#input-sections)
* [RASPA stages](#raspa-stages)
* [General options](#general-options)
  * [Simulation types](#simulation-types)
  * [Simulation duration](#simulation-duration)
  * [Restart and crash-recovery](#restart-and-crash-recovery)
  * [Printing options](#printing-options)
  * [Parameter tuning](#parameter-tuning)
  * [Threading and reproducibility](#threading-and-reproducibility)
  * [Systems & Components](#systems-components)
* [System options](#system-options)
  * [Operating conditions and thermostat/barostat-parameters](#operating-conditions-and-thermostatbarostat-parameters)
  * [Box/Framework options](#boxframework-options)
  * [Force field definition](#force-field-definition)
  * [System `MC`-moves](#system-mc-moves)
  * [Cross-links between molecules](#cross-links)
  * [Molecular dynamics parameters](#molecular-dynamics-parameters)
  * [Options to measure properties](#options-to-measure-properties)
    * [Output pdb-movies](#output-pdb-movies)
    * [Histogram of the energy](#histogram-of-the-energy)
    * [Histogram of the number of molecules](#histogram-of-the-number-of-molecules)
    * [Radial Distribution Function (RDF) force-based](#radial-distribution-function-rdf-force-based)
    * [Radial Distribution Function (RDF) conventional](#radial-distribution-function-rdf-conventional)
    * [Mean-Squared Displacement (MSD) order-N](#mean-squared-displacement-msd-order-n)
    * [Density grids](#density-grids)
* [Force field options](#force-field-options)
  * [Configurational-bias and recoil-growth options](#cbmc-options)
* [Component options](#component-options)
  * [Component properties](#component-properties)
  * [Component analyses](#component-analyses)
    * [Histograms of the intra-molecular geometry](#histograms-of-the-intra-molecular-geometry)
    * [Molecule shape: the gyration-tensor family](#molecule-shape)
    * [Molecule backbone: chain statistics](#molecule-backbone)
    * [End-to-end vector autocorrelation function and relaxation time](#end-to-end-autocorrelation-function)
  * [Component `MC`-moves](#component-mc-moves)
<!-- TOC -->

----------------------------------------------------------------------------------

## Input sections <a name="input-sections"></a>

A minimal input file has three parts: a set of top-level (general) options, a
`"Systems"` list, and a `"Components"` list. The example below runs a molecular
dynamics simulation of CO<sub>2</sub> in the framework Cu-BTC.

```json
{
  "SimulationType" : "MolecularDynamics",
  "NumberOfProductionCycles" : 100000,
  "NumberOfInitializationCycles" : 1000,
  "NumberOfEquilibrationCycles" : 10000,
  "PrintEvery" : 1000,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "Cu-BTC",
      "NumberOfUnitCells" : [1, 1, 1],
      "ChargeMethod" : "Ewald",
      "ExternalTemperature" : 323.0,
      "ExternalPressure" : 1.0e4,
      "OutputPDBMovie" : false,
      "SampleMovieEvery" : 10
    }
  ],

  "Components" : [
    {
      "Name" : "CO2",
      "FugacityCoefficient" : 1.0,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "SwapProbability" : 0.0,
      "WidomProbability" : 0.0,
      "CreateNumberOfMolecules" : 20
    }
  ]
}
```

Unknown keywords are rejected: if a general, system, component, or reaction key
is not recognized, `RASPA` stops with an error naming the offending key. This
helps catch typos early.

----------------------------------------------------------------------------------

## RASPA stages <a name="raspa-stages"></a>

A Monte Carlo simulation in `RASPA` is executed as a sequence of four
consecutive stages. Each stage runs a number of cycles (set through the
corresponding `NumberOf...Cycles` command) and writes its own status report to
the output file. The stages always run in the order below, and a stage is simply
skipped when its number of cycles is `0`.

1.  **Pre-initialization**
    (`"NumberOfPreInitializationCycles"`)\
    An optional relaxation stage that runs *before* the regular initialization.
    It uses a restricted set of moves only: translation, rotation, reinsertion,
    and partial-reinsertion. Because none of these moves changes the number of
    molecules, this stage keeps the composition of the system fixed while
    relaxing the initial configuration. It is mainly used to remove close
    contacts and overlaps that can appear right after molecules are created or
    read from a restart file, so that the subsequent stages start from a
    reasonable configuration. Statistics gathered here are not used for the
    final averages.

2.  **Initialization**
    (`"NumberOfInitializationCycles"`)\
    The first stage in which the full set of configured Monte Carlo moves is
    used, including moves that insert and delete molecules (e.g. swap moves in
    the grand-canonical ensemble). This stage brings the system towards its
    equilibrium state, for example towards the equilibrium loading in an
    adsorption simulation. Statistics gathered here are not used for the final
    averages.

3.  **Equilibration**
    (`"NumberOfEquilibrationCycles"`)\
    Continues to equilibrate the system with the full set of moves. For
    Continuous Fractional Component Monte Carlo (`CFCMC`) this stage is also used
    to measure the λ biasing factors using Wang-Landau estimation, so that the
    fractional molecule samples all λ values uniformly during production.
    Statistics gathered here are not used for the final averages.

4.  **Production**
    (`"NumberOfProductionCycles"`)\
    The main stage during which the thermodynamic properties of interest
    (loadings, energies, pressures, enthalpies of adsorption, radial distribution
    functions, etc.) are sampled and averaged. Block averages are computed over
    this stage to provide error estimates. Any biasing factors determined during
    equilibration are kept fixed.

In all cases a Monte Carlo *cycle* consists of $N$ steps, where $N$ is the
number of molecules (with a minimum of 20). During each cycle, on average, one
Monte Carlo move is attempted per molecule. The CPU time spent in each stage is
reported separately at the end of the simulation.

----------------------------------------------------------------------------------

## General options <a name="general-options"></a>

### Simulation types <a name="simulation-types"></a>

-   `"SimulationType" : "MonteCarlo"`\
    Runs the Monte Carlo engine. The ensemble is not stated explicitly but is
    deduced from the Monte Carlo moves that are switched on. A hybrid MC/MD
    scheme can be obtained by enabling the hybrid `MD`-move.

-   `"SimulationType" : "MolecularDynamics"`\
    Runs the Molecular Dynamics engine. The ensemble must be specified
    explicitly through the `"Ensemble"` key.

-   `"SimulationType" : "MolecularDynamicsSpatialDecomposition"`\
    (aliases `"SpatialDecompositionMolecularDynamics"`, `"MolecularDynamicsSD"`)\
    Molecular dynamics of molecules in a box with a multithreaded
    spatial-decomposition force engine, for large systems (tens of thousands
    of atoms and more). The pre-initialization and initialization stages are
    the ordinary serial Monte Carlo cycles of `"MolecularDynamics"`; the
    equilibration and production stages integrate the equations of motion
    with forces from a domain-decomposed engine that replaces the
    all-pairs and direct k-space code:

    -   the box is divided into `"NumberOfThreads"` sub-domains (a
        \f$p_x \times p_y \times p_z\f$ grid of staggered slabs with the
        smallest surface, or the grid given by `"DomainGrid"`) whose cut
        planes are placed on the atom positions at every neighbour-list
        rebuild so that all threads own the same number of atoms; each
        thread works on a compact private copy of its atoms and of the
        *ghost images* it interacts with (atoms of other sub-domains, or
        periodic images, with the periodic shift resolved when the Verlet
        list is built, so the pair kernel needs no minimum-image
        operation). Every pair of the system is evaluated once with
        Newton's third law; the forces on the ghost images are collected
        by the owning threads after the pair phase. The neighbour search
        uses cells of a third of the cutoff plus skin and tests the
        candidates per cell pair on wrapped positions with the box
        translation of the cell pair applied once, and the lists are
        rebuilt only when an atom has moved more than half the skin;
    -   for a force field of plain 12-6 Lennard-Jones pairs (truncated or
        shifted) with Ewald or no electrostatics and fully coupled atoms, a
        specialised pair kernel is used with the Ewald real-space term
        \f$\mathrm{erfc}(\alpha r)/r\f$ tabulated as a cubic Hermite spline
        in \f$r^2\f$ (no square root or division per pair; relative error
        below \f$10^{-9}\f$, force consistent with the interpolated energy);
        other pair potentials, charge methods or scaled (fractional) atoms
        use the generic kernels of the rest of the code (the status output
        names the kernel);
    -   the Ewald reciprocal sum is replaced by a particle-mesh (SPME/PPPM)
        solver: B-spline charge assignment of order
        `"PPPMInterpolationOrder"` on a mesh of spacing `"PPPMMeshSpacing"`,
        FFTs by FFTW, with the same \f$\alpha\f$ as the Ewald summation of
        the force field, so the real-space, self and exclusion terms are
        unchanged and the result converges to the exact Ewald energy and
        forces as the mesh is refined;
    -   bonded terms, the Ewald self and intra-molecular exclusion
        corrections and the molecular pressure tensor (pair, reciprocal,
        exclusion and tail contributions, corrected to the molecular
        center-of-mass virial) are computed by the same threads, so the
        NPT and NPT-PR barostats work without the \f$O(N^2)\f$ pressure
        evaluation.

    The threads form a persistent team (one worker per sub-domain, the
    calling thread included) synchronized by barriers; `"NumberOfThreads" :
    1` is a valid serial cell-list/particle-mesh run through the same code
    path. Every sub-domain must be at least one cell wide along every
    axis, so the number of threads must fit the box: with
    \f$c_i = \lfloor 4 L_i^{\perp} / (r_c + \text{skin}) \rfloor\f$
    cells along axis \f$i\f$ at the finest subdivision, at most
    \f$c_x c_y c_z\f$ sub-domains are possible and each sub-domain axis
    needs \f$p_i \le c_i\f$; the input is rejected with a message when
    they do not. The minimum-image neighbour lists require
    \f$r_c + \text{skin} \le L^{\perp}_{\min} / 2\f$: since the automatic
    Coulomb cutoff (`"CutOffCoulomb" : "auto"`) is exactly half the box,
    set an explicit `"CutOffCoulomb"` (and `"CutOffVDW"`) in the force
    field.

    At the start of the equilibration and of the production stage the
    engine's energy and forces are compared against the exact code once and
    the differences are reported in the output (`Spatial-decomposition
    force engine check`), together with the cell grid, the sub-domain
    layout, the number of atoms and pairs per sub-domain and the mesh size.
    Timings per phase (neighbour-list rebuilds, pairs, mesh, bonded terms)
    are printed at the end.

    Scope: molecules in a box only (rigid, semi-flexible and flexible
    components, Lennard-Jones plus Coulomb with the Ewald, damped
    shifted-force, Wolf, modified shifted-force or zero-dipole methods,
    bonded potentials), ensembles `"NVE"`, `"NVT"`, `"NPT"` and `"NPTPR"`.
    Frameworks, external fields, polarization, cross-link bonds,
    CFCMC/fractional molecules, `"OmitInterInteractions"`,
    `"UseDualCutOff"`, the stress-fluctuation elastic constants and the
    particle-exchange ensembles (`"MuVT"`, `"MuPT"`, `"MuPTPR"`) are rejected
    at input time with a message that points to `"SimulationType" :
    "MolecularDynamics"`. The per-component energy decomposition
    (`Energy averages and statistics`) is sampled with the exact
    \f$O(N^2)\f$ code every `"PrintEvery"` cycles only; the total energies,
    the conserved-energy drift, the pressure tensor and all property samplers
    use the engine every step. Binary restarts continue with the engine
    settings of the (new) input file, not of the restart file.

-   `"SimulationType" : "MonteCarloTransitionMatrix"`\
    Runs Monte Carlo with transition-matrix (TMMC) biasing enabled for every
    system. See the macro-state keywords in the system options.

-   `"SimulationType" : "ParallelTempering"`\
    Runs a multithreaded parallel-tempering (replica-exchange) Monte Carlo
    simulation. Exactly one system is declared in the input, with a temperature
    ladder given by the system key `"ExternalTemperatures"` (a sorted list of
    at least two temperatures); the system is replicated internally into one
    replica per temperature, and every replica runs in its own thread with its
    own random-number stream.

    Every `"ParallelTemperingSwapEvery"` cycles (default `10`, `0` disables)
    the threads synchronize on a barrier and configuration swaps between
    replicas at neighboring temperatures are attempted with acceptance rule
    min(1, exp[(β_B − β_A)(U_B − U_A)]) (extended with the Yan & de Pablo
    fugacity factor Π_i [(β_A f_A,i)/(β_B f_B,i)]^(N_B,i − N_A,i) when the
    molecule counts differ, and with the PV work term when the boxes travel
    with the configurations). The replicas keep their temperatures;
    only the configurations migrate through the ladder. The pairing offset
    alternates between sweeps, so a configuration can traverse the whole
    ladder. This is the only synchronization point between the threads besides
    the start and end of each stage.

    Every replica writes its own output file
    `output/output_{T_k}_{P}.parallel_tempering.r{k}.txt` (and `.json`) with
    the standard status reports and final averages; the combined file
    `output/output.parallel_tempering.txt` (and `.json`) holds the swap
    statistics, including a per-pair acceptance table (low acceptance for a
    particular pair marks a bottleneck in the ladder; use a denser ladder
    there) and the replica round-trip statistics: the number of round trips
    (a configuration returning to the coldest replica after having visited the
    hottest one), the mean round-trip time, the round trips of every
    individual configuration, and the up-moving fraction \f$f(T)\f$ per
    temperature (Katzgraber et al. 2006; linear from 1 to 0 for an optimal
    ladder, a plateau marks a bottleneck). The round-trip count, not the swap
    acceptance, is the measure of how well the configurations actually mix
    through the ladder. The optional property files (RDFs, density grids, energy and
    molecule-count histograms, molecule properties, movies, and the
    number-of-molecules/volume evolution files) are written per replica, keyed
    by the replica index (`.s{k}`). Restart files (JSON and binary) are not
    supported by this driver.

    The driver spawns one worker thread per temperature (plain C++ threads);
    leave `"NumberOfThreads"` at its default of `1` so the
    per-energy-evaluation thread pool stays serial and the machine is not
    oversubscribed. Note that the swap move requires rigid, whole-molecule
    replicas: systems with fractional (CFCMC) molecules, flexible components,
    or reactions reject all swap attempts.

-   `"SimulationType" : "ParallelTemperingMolecularDynamics"`\
    (aliases `"ReplicaExchangeMolecularDynamics"`, `"ParallelTemperingMD"`)\
    Runs a multithreaded replica-exchange molecular-dynamics (REMD) simulation
    (Sugita & Okamoto, Chem. Phys. Lett. 314, 141-151, 1999). Exactly one
    system is declared in the input, with a temperature ladder given by the
    system key `"ExternalTemperatures"` (a sorted list of at least two
    temperatures) and an `"Ensemble"` of `"NVE"` or `"NVT"` (the swap
    exchanges configurations and momenta only; the cell and the number of
    molecules of every replica are fixed, so the barostat and
    particle-exchange ensembles are rejected). The system is replicated
    internally into one replica per temperature, and every replica is
    integrated in its own thread with its own random-number stream and its
    own Nosé–Hoover chain at its own temperature.

    As in `"MolecularDynamics"`, the pre-initialization and initialization
    stages are Monte Carlo cycles that relax the initial configuration; the
    equilibration and production stages integrate the equations of motion,
    one time step per cycle (Maxwell–Boltzmann velocities at the replica
    temperature are drawn at the start of the equilibration stage). Every
    `"ParallelTemperingSwapEvery"` cycles (default `10`, `0` disables; for
    MD a cycle is one time step, so a larger value such as `100`–`1000` is
    usually appropriate) the threads synchronize on a barrier and
    configuration swaps between replicas at neighboring temperatures are
    attempted with acceptance rule min(1, exp[(β_B − β_A)(U_B − U_A)]) on the
    potential energies. After an accepted swap the momenta that travelled
    with the configuration are rescaled by sqrt(T_new/T_old), so the kinetic
    parts of the Boltzmann factors cancel and the kinetic temperature of each
    replica is unchanged by the exchange; the thermostat chain is a property
    of the heat bath and stays with the replica; the conserved-energy
    reference used for the drift bookkeeping is reset (the extended-system
    energy jumps at a swap by construction). The pairing offset alternates
    between sweeps, so a configuration can traverse the whole ladder.

    Every replica writes its own output file
    `output/output_{T_k}_{P}.parallel_tempering_md.r{k}.txt` (and `.json`)
    with the MD status reports (kinetic temperatures, conserved-energy drift)
    and final averages; the combined file
    `output/output.parallel_tempering_md.txt` (and `.json`) holds the swap
    statistics, including a per-pair acceptance table (low acceptance for a
    particular pair marks a bottleneck in the ladder; use a denser ladder
    there), the replica round-trip statistics (round trips, mean round-trip
    time, up-moving fraction \f$f(T)\f$ per temperature; see
    `ParallelTempering`), and the potential- and conserved-energy drift of
    every replica.
    The optional property files (RDFs, density grids, MSD, VACF, energy and
    molecule-count histograms, molecule properties, and the
    number-of-molecules/volume evolution files) are written per replica,
    keyed by the replica index (`.s{k}`). Binary restart files are written
    at the barrier synchronization points (`"WriteBinaryRestartEvery"`, on a
    shutdown signal) and resumed with `"RestartFromBinaryFile"`.

    The driver spawns one worker thread per temperature (plain C++ threads);
    leave `"NumberOfThreads"` at its default of `1` so the
    per-energy-evaluation thread pool stays serial and the machine is not
    oversubscribed. The swap has the same compatibility requirements as the
    Monte Carlo variant (same Hamiltonian and topology in all replicas; no
    reactions or pair/group/Gibbs fractional molecules).

-   `"SimulationType" : "HyperParallelTempering"`\
    Runs a multithreaded hyper-parallel-tempering (replica-exchange) Monte
    Carlo simulation over a two-dimensional grid of state points (Yan & de
    Pablo, JCP 111(21), 9509-9516, 1999). Exactly one system is declared in the
    input, with a temperature ladder `"ExternalTemperatures"` and a pressure
    ladder `"ExternalPressures"` (both sorted lists); the system is replicated
    internally into one replica per (temperature, pressure) grid point, and
    every replica runs in its own thread with its own random-number stream.
    The pressures are converted to per-component fugacities internally: the
    fugacity coefficients are recomputed with the Peng-Robinson equation of
    state at every grid point (an explicitly given `"FugacityCoefficient"` is
    ignored, since a single value cannot be valid at all state points).

    Every `"ParallelTemperingSwapEvery"` cycles (default `10`, `0` disables)
    the threads synchronize on a barrier and configuration swaps between
    replicas at neighboring grid points are attempted with the Yan & de Pablo
    acceptance rule

    min(1, exp[(β_B − β_A)(U_B − U_A)] × Π_i [(β_A f_A,i)/(β_B f_B,i)]^(N_B,i − N_A,i))

    with per-component fugacities f_X,i (the second factor accounts for the
    different numbers of adsorbed molecules in the two configurations). The
    sweeps alternate between the temperature direction (neighboring
    temperatures at the same pressure) and the pressure direction (neighboring
    pressures at the same temperature), each with an alternating pairing
    offset, so a configuration can traverse the whole grid. The replicas keep
    their (temperature, pressure) state points; only the configurations
    migrate. This is the only synchronization point between the threads
    besides the start and end of each stage.

    Every replica writes its own output file
    `output/output_{T}_{P}.hyper_parallel_tempering.r{k}.txt` (and `.json`)
    with the standard status reports and final averages (one adsorption
    isotherm/isobar point per replica); the combined file
    `output/output.hyper_parallel_tempering.txt` (and `.json`) holds the swap
    statistics, with separate per-pair acceptance tables for the temperature
    and the pressure direction (low acceptance for a particular pair marks a
    bottleneck in the grid; use a denser ladder there). The optional property
    files (RDFs, density grids, energy and molecule-count histograms, molecule
    properties, movies, and the number-of-molecules/volume evolution files)
    are written per replica, keyed by the replica index (`.s{k}`). Restart
    files (JSON and binary) are not supported by this driver.

    The assembled adsorption isotherms are additionally written to one
    gnuplot-friendly file per temperature,
    `output/isotherm_{T}.hyper_parallel_tempering.txt`, with one block per
    component and rows ordered by increasing pressure (columns: fugacity [Pa],
    absolute loading with its confidence-interval error in molecules/cell,
    molecules/unit-cell, mol/kg-framework and mg/g-framework, and the pressure
    [Pa]). The files are overwritten every `"PrintEvery"` cycles at a swap
    synchronization point during production, and once more at the end of the
    run, so the convergence of the isotherms can be monitored while the
    simulation is running.

    The driver spawns one worker thread per grid point (plain C++ threads) — with N_T temperatures and N_P pressures that is N_T × N_P
    threads, so size the grid to the machine. Leave `"NumberOfThreads"` at its
    default of `1` so the per-energy-evaluation thread pool stays serial. The
    swap move requires rigid, whole-molecule replicas: systems with fractional
    (CFCMC) molecules, flexible components, or reactions reject all swap
    attempts.

-   `"SimulationType" : "ReweightedHistogram"`\
    Runs the same multithreaded replica-exchange simulation over a
    (temperature, pressure) grid as `"HyperParallelTempering"` (one system with
    the ladders `"ExternalTemperatures"` and `"ExternalPressures"`, one thread
    per grid point, Yan & de Pablo configuration swaps, per-replica and
    combined output files, directly-measured isotherm files) and additionally
    combines the results of all threads with multiple-histogram reweighting
    (WHAM; Ferrenberg & Swendsen, PRL 63, 1195, 1989; Kumar et al., J. Comput.
    Chem. 13, 1011, 1992) into a continuous isotherm surface.

    A cycle is a fixed number of Monte Carlo moves equal to
    `"MacroStateMaximumNumberOfMolecules"` (the same filling ceiling TMMC uses),
    not the current occupancy, so empty and filled replicas do the same amount of
    work per cycle. With `"ComputeBET"` the key may be `"auto"` or omitted: an
    unbiased GCMC scout at nitrogen P0 places N_max just above the plateau
    loading.

    Every `"SampleReweightingEvery"` production cycles (default `5`) each
    replica records a raw (N, U) sample: the molecule count of the adsorbate
    and the potential energy. At the end of the run the pooled samples are
    binned over (N, U) and the WHAM self-consistent equations for the
    grand-canonical density of states Ω(N, U) are solved in log space,

    ln Ω(N,U) = ln M(N,U) − ln Σ_i n_i exp[g_i − β_i U + N ln(β_i f_i)],
    g_i = −ln Σ_{N,U} Ω(N,U) exp[−β_i U + N ln(β_i f_i)]

    (M is the total count of a bin, n_i the sample count of grid point i, f_i
    its fugacity). The density of states is then reweighted to arbitrary
    (temperature, fugacity) points, giving smooth isotherms in between (and
    somewhat beyond) the simulated grid points from a single run. The error
    bars are obtained by re-solving the WHAM equations per block.

    Outputs, in addition to the hyper-parallel-tempering ones (which use the
    tag `reweighted_histogram` in the filenames): one reweighted isotherm per
    requested temperature, `output/reweighted_isotherm_{T}.reweighted_histogram.txt`,
    evaluated on a fine log-spaced pressure grid (fugacity coefficients from
    the Peng-Robinson equation of state at every point; columns as in the
    directly-measured isotherm files, plus the effective sample size that
    diagnoses the overlap — small values flag extrapolation beyond the sampled
    (N, U) region); the per-state free energies g_i together with a
    self-consistency table (reweighted vs. directly-measured loading at every
    simulated grid point) in
    `output/reweighted_free_energies.reweighted_histogram.txt`; and the WHAM
    convergence report in the combined output file. The analysis controls are
    the optional top-level keys `"ReweightingTemperatures"` (default: the
    temperature ladder), `"ReweightingPressureRange"` (default: the span of
    the pressure ladder) and `"ReweightingNumberOfPressures"` (default `100`).

    For bulk boxes (no framework) the analysis additionally computes the
    vapor-liquid equilibrium at every requested subcritical temperature with
    the equal-weight criterion (Wilding): the fugacity is bisected until the
    vapor and liquid peaks of the bimodal reweighted molecule-number
    distribution P(N) carry equal probability weight. The coexistence table
    `output/vle_coexistence.reweighted_histogram.txt` holds, per temperature,
    the coexistence fugacity, the saturation pressure (from β p V = ln Ξ,
    normalized exactly by the empty-box state when it was sampled, otherwise
    approximately by an ideal-gas reference at the dilute end of the pressure
    range), and the saturated vapor and liquid densities in kg/m³, all with
    per-block error bars; the distribution at coexistence is written to
    `output/vle_distribution_{T}.reweighted_histogram.txt` (for inspection and
    finite-size scaling). Temperatures where no bimodal distribution is found
    in the scanned pressure range (supercritical, or the sampling does not
    connect the phases) are flagged. Practical setup: place the top of the
    temperature ladder close to the critical point (there the vapor and
    liquid histograms overlap and configurations can cross between the
    phases), span the coexistence pressures with the pressure ladder, and
    start from a liquid-density configuration (`"CreateNumberOfMolecules"`) so
    both phases are visited — the liquid does not nucleate spontaneously from
    the vapor in grand-canonical sampling.

    The reweighting reuses each sampled energy U(x) at all temperatures, so
    temperature-dependent potentials (Feynman-Hibbs) are rejected, and exactly
    one component is required (the histograms are collected over its molecule
    count). Reliable reweighting requires the (N, U) histograms of neighboring
    grid points to overlap — the same requirement as a healthy swap acceptance,
    so the per-pair acceptance tables double as an overlap diagnostic.

-   `"SimulationType" : "ParallelTMMC"`\
    Runs a multithreaded transition-matrix Monte Carlo simulation with
    windowed macrostate walkers (Errington, JCP 118, 9915, 2003; Shen &
    Errington, JCP 122, 064508, 2005). Exactly one system with exactly one
    component is declared. The macrostate range
    `"MacroStateMinimumNumberOfMolecules"` to
    `"MacroStateMaximumNumberOfMolecules"` is split into `"NumberOfWindows"`
    contiguous windows that share their endpoint macrostates, and the system
    is replicated into one walker per (temperature, window) pair of the
    ladder `"ExternalTemperatures"` × windows (a single
    `"ExternalTemperature"` gives a one-temperature run). A cycle is a fixed
    number of Monte Carlo moves equal to the global
    `"MacroStateMaximumNumberOfMolecules"` (not the current occupancy or the
    window width), so every walker does the same amount of work per cycle.
    Every walker runs
    grand-canonical Monte Carlo in its own thread with its own random-number
    stream, confined to its window and flattened by the transition-matrix
    bias, which is re-derived from the collection matrix every
    `"TMMCUpdateEvery"` steps (default `100000`). Each walker starts inside
    its window: molecules are grown with CBMC up to the lower window boundary
    (`"CreateNumberOfMolecules"` must not exceed the macrostate minimum). The
    walkers are fully independent — there are no swaps and no barriers, so
    the threads only join at the stage boundaries.

    The collection matrix records the unbiased acceptance probabilities of
    all attempted insertions and deletions — also those rejected at the
    window bounds — so the collection matrices of the windows of one
    temperature simply add, and the macrostate probability distribution over
    the full range follows from detailed balance,
    ln Π(N+1) = ln Π(N) + ln P(N→N+1) − ln P(N+1→N). The distribution is
    exact at the reference fugacity (from `"ExternalPressure"` through the
    Peng-Robinson equation of state at every temperature) and reweights
    exactly to any other fugacity, ln Π(N; f) = ln Π(N; f_ref) + N ln(f/f_ref).
    Every temperature is solved from its own walkers only, so
    temperature-dependent potentials (Feynman-Hibbs) are allowed. The error
    bars come from re-deriving ln Π from the production-only per-block
    increments of the collection matrices. The collection matrix itself is
    never reset (its entries are valid across bias updates), so the final
    ln Π uses the statistics of the equilibration and production stages
    combined.

    Outputs: per temperature, the macrostate distribution
    `output/lnpi_{T}.parallel_tmmc.txt` (ln Π(N) with per-block errors and
    the visit histogram) and the reweighted isotherms on a log-spaced
    pressure grid spanning `"ReweightingPressureRange"` (default: four
    decades around the reference pressure) with
    `"ReweightingNumberOfPressures"` points (fugacity coefficients from the
    Peng-Robinson equation of state at every point; note the isotherms
    saturate artificially when ⟨N⟩ approaches the upper macrostate bound):
    the equilibrium isotherm
    `output/reweighted_isotherm_{T}.parallel_tmmc.txt` (averaged over both
    basins of Π(N)) and the adsorption and desorption branches
    `output/adsorption_isotherm_{T}.parallel_tmmc.txt` and
    `output/desorption_isotherm_{T}.parallel_tmmc.txt` — where Π(N; f) is
    bimodal (a first-order transition: capillary condensation in a pore,
    vapor-liquid in a box) the adsorption branch is the conditional average
    over the low-density basin (the metastable states followed on the way
    up), the desorption branch the conditional average over the high-density
    basin, and together they trace the hysteresis loop around the
    equilibrium step; where Π(N; f) is unimodal all three coincide (the last
    column flags the bimodal points).
    Every walker writes its own output file
    `output/output_{T}_w{w}.parallel_tmmc.r{k}.txt` (and `.json`; the direct
    averages there are biased flat-histogram averages, diagnostics only) and
    its transition-matrix statistics to
    `tmmc/tmmc_statistics_{T}_w{w}.parallel_tmmc.txt`. The combined file
    `output/output.parallel_tmmc.txt` (and `.json`) holds the macrostate
    coverage per walker (every state of a window must be visited for the
    stitched ln Π to be reliable) and the analysis report. Restart files
    (JSON and binary) are not supported by this driver.

    For bulk boxes (no framework) the analysis additionally computes the
    vapor-liquid equilibrium at every simulated subcritical temperature with
    the equal-weight criterion (Wilding): the fugacity is bisected until the
    vapor and liquid peaks of the bimodal Π(N) carry equal probability
    weight. The coexistence table `output/vle_coexistence.parallel_tmmc.txt`
    holds, per temperature, the coexistence fugacity, the saturation pressure
    (from β p V = ln Ξ, normalized exactly by the empty-box state — set the
    macrostate minimum to `0` for this), and the saturated vapor and liquid
    densities in kg/m³, all with per-block error bars; the distribution at
    coexistence is written to `output/vle_distribution_{T}.parallel_tmmc.txt`.
    Practical setup: the macrostate maximum must comfortably exceed the
    liquid peak (ρ_liq V), and the scanned pressure range must bracket the
    saturation pressures of all the temperatures. In contrast to
    grand-canonical sampling at a single state point, the flat-histogram walk
    crosses the vapor-liquid gap by construction — no starting configuration
    tricks are needed.

    The driver spawns one worker thread per walker (plain C++ threads) — with N_T temperatures and N_W windows that is N_T × N_W
    threads, so size the grid to the machine. Leave `"NumberOfThreads"` at
    its default of `1` so the per-energy-evaluation thread pool stays serial.
    More windows shorten the equilibration (each walker only needs to flatten
    its own window) but every window must still be crossed many times for a
    reliable stitched distribution.

-   `"SimulationType" : "Minimization"`\
    Performs an energy minimization of the initial configuration.

    `"ComputeElasticConstants" : boolean` optionally computes the static,
    relaxed elastic tensor after convergence. The calculation uses the full
    six-component symmetric strain basis even when the minimization used a
    fixed or restricted cell. It reports the Born, internal-relaxation, and
    hydrostatic-pressure terms, the stiffness and compliance matrices, Born
    stability eigenvalues, and Voigt/Reuss/Hill moduli. Physical-unit output is
    in GPa (compliance in GPa^-1); shear entries use engineering-strain Voigt
    order `xx, yy, zz, yz, xz, xy`.

    `"ComputeNormalModes" : boolean` optionally performs a Gamma-point
    normal-mode analysis after convergence. The generalized Hessian is
    mass-weighted (atomic masses for Cartesian degrees of freedom; molecular
    mass and the space-frame inertia tensor for rigid-molecule center-of-mass
    and orientation degrees of freedom) and diagonalized. Frequencies are
    reported per mode as omega^2, THz, cm^-1, and meV (reduced units in the
    reduced unit system); imaginary modes appear as negative frequencies.
    Orientation directions with vanishing moment of inertia (single-bead or
    linear rigid molecules) are excluded from the mass metric and show up as
    zero modes.

    `"NormalModeMovies" : boolean` optionally writes an animated PDB movie of
    every normal mode into a `normal_modes` directory (implying the normal-mode
    analysis). Each file `mode_XXXX.s{system}.pdb` animates the atoms
    oscillating along the mode's displacement pattern. `"NormalModeMoviePeriods"
    : integer` sets the number of full oscillation periods shown per movie
    (default `1`), `"NormalModeMoviePointsPerPeriod" : integer` sets the number
    of frames sampled within one period (default `16`), and
    `"NormalModeMovieAmplitude" : number` sets the maximum atomic displacement in
    Angstrom used to scale each mode (default `0.5`).

    `"ComputePhononDispersion" : boolean` optionally computes the phonon band
    structure after convergence. The image-resolved force constants are Fourier
    transformed to the k-dependent dynamical matrix `D(k)` (including the
    reciprocal-space Ewald term for charged systems), mass-weighted, and
    diagonalized along a high-symmetry path. Rigid molecules are handled in
    generalized center-of-mass/orientation coordinates, so at the Gamma point the
    result matches `ComputeNormalModes`. `"PhononDispersionPath" : array` defines
    the path as a list of nodes in fractional reciprocal-lattice coordinates,
    each node being either a bare `[kx, ky, kz]` array or an object
    `{"Label": "X", "kPoint": [kx, ky, kz]}`; consecutive nodes form the segments.
    When omitted, a default connected `G-X-G-Y-G-Z` star along the reciprocal axes
    is used. `"PhononDispersionPointsPerSegment" : integer` sets the number of
    sampled k-points per segment (default `20`). Results are written to
    `output/minimization.s{system}.json` and to a gnuplot-friendly band file
    `output/phonon_dispersion.s{system}.txt` (frequencies in cm^-1, negative for
    imaginary/unstable modes).

    `"ElasticEigenvalueTolerance" : real` controls the relative spectral
    threshold used to remove translational and rotational zero modes from the
    internal Hessian (default `1.0e-8`). A significant negative internal mode
    is treated as an unstable minimized structure.

-   `"SimulationType" : "ThermodynamicIntegration"`\
    Runs a Monte Carlo simulation at a fixed value of the CFCMC coupling
    parameter λ. Every component with a `"LambdaBinIndex"` gets fractional
    molecule(s) pinned at λ = binIndex / (`NumberOfLambdaBins` − 1); no
    λ-changing moves and no Wang-Landau biasing are involved. The pinned λ
    coordinate follows the component definition: with `"GroupComponents"` the
    group-swap λ is pinned (fractional molecules for the central component and
    all satellites), with `"PairComponent"` the ion-pair λ (fractional
    molecules for both components of the pair), otherwise the grand-canonical
    λ (a single fractional molecule). During production the ensemble average
    ⟨∂U/∂λ⟩ is accumulated at that λ and reported (with a block-error
    estimate) at the end of the run. Running one simulation per λ-bin and
    integrating ⟨∂U/∂λ⟩ from λ=0 to λ=1 yields the excess chemical potential
    via thermodynamic integration.

-   `"SimulationType" : "ParallelThermodynamicIntegration"`\
    Computes the full ⟨∂U/∂λ⟩(λ) curve and the excess chemical potential in a
    single multithreaded run. Exactly one system is declared in the input; it
    is replicated internally into one replica per λ-bin (replica *k* starts
    pinned at λ-bin *k*), and every replica runs in its own thread with its own
    random-number stream. Exactly one component is marked as the
    thermodynamic-integration component with `"ThermodynamicIntegration" :
    true` (a `"LambdaBinIndex"` marking also works; its value is ignored). The
    pinned λ coordinate is inferred from the component definition exactly as
    for `"ThermodynamicIntegration"` (group-swap, ion-pair, or grand-canonical
    λ).

    Every `"LambdaExchangeEvery"` cycles (default `10`, `0` disables) the
    threads synchronize on a barrier and Hamiltonian parallel-tempering
    exchanges of the λ values between replicas at neighboring λ-bins are
    attempted, with acceptance rule
    min(1, exp[−β(ΔU_A + ΔU_B)]) where ΔU_A = U_A(λ_B) − U_A(λ_A) and
    ΔU_B = U_B(λ_A) − U_B(λ_B). This is the only synchronization point between
    the threads besides the start and end of each stage. Because the exchanges
    permute the λ values over the replicas, every λ-bin is occupied by exactly
    one replica at all times and the whole curve is sampled every cycle.

    At the end the per-bin ⟨∂U/∂λ⟩ book-keeping of all replicas is stitched
    together and integrated from λ=0 to λ=1 with three quadrature rules, each
    with a confidence interval from its block-wise integrals: a natural cubic
    spline (the recommended value), composite Simpson (1/3 rule with a 3/8
    tail for odd interval counts), and the trapezoidal rule as a baseline.
    Because the data live on fixed equidistant λ-bins, Gaussian quadrature is
    not applicable (it would require samples at the non-equidistant Gauss
    nodes); the Newton–Cotes rules and the spline are the applicable choices.
    A quadrature spread far below the sampling error confirms the λ-grid
    resolves the curvature of the curve; a large spread signals that more
    λ-bins are needed. The combined
    output file `output/output_{T}_{P}.parallel_ti.txt` (and `.json`) holds
    the stitched results and the λ-exchange statistics. In addition every
    replica writes its own output file
    `output/output_{T}_{P}.parallel_ti.r{k}.txt` (and `.json`) with the
    standard status reports every `"PrintEvery"` cycles (including the λ-bin
    the replica currently occupies), and at the end its energy-drift check,
    Monte-Carlo move statistics, averages, and its own per-bin ⟨∂U/∂λ⟩
    book-keeping (the bins it visited through accepted λ-exchanges).
    During production the current stitched curve and its running spline
    integral are additionally written to the gnuplot-friendly snapshot file
    `output/dudlambda_{T}_{P}.parallel_ti.txt` (columns: λ, ⟨∂U/∂λ⟩ [K],
    error [K]; overwritten every `"PrintEvery"` cycles at a λ-exchange
    synchronization point, and once more at the end of the run), so the
    convergence of the curve can be monitored while the simulation is
    running. Binary restart files are not supported by this driver.

    The driver spawns `NumberOfLambdaBins` worker threads (plain C++ threads);
    leave `"NumberOfThreads"` at its default of `1` so the
    per-energy-evaluation thread pool stays serial and the machine is not
    oversubscribed.

### Simulation duration <a name="simulation-duration"></a>

-   `"NumberOfProductionCycles" : integer`\
    The number of cycles in the production run. For Monte Carlo a cycle consists
    of $N$ steps, where $N$ is the number of molecules with a minimum of 20. On
    average, one Monte Carlo move is therefore attempted per molecule per cycle
    (whether accepted or rejected). For Molecular Dynamics the number of cycles
    is simply the number of integration steps.

-   `"NumberOfBlocks" : integer`\
    Number of contiguous production blocks used for confidence-interval
    estimates. At least three blocks are required; default: `5`.

-   `"NumberOfInitializationCycles" : integer`\
    The number of Monte Carlo cycles used to bring the system towards
    equilibrium. This applies to both Monte Carlo and Molecular Dynamics runs and
    is useful to relax the initial atomic positions before production.

-   `"NumberOfEquilibrationCycles" : integer`\
    For Molecular Dynamics, the number of steps used to equilibrate the system
    velocities before production starts. For Monte Carlo, and in particular
    `CFCMC`, the equilibration phase is used to measure the biasing factors using
    Wang-Landau estimation.

-   `"NumberOfPreInitializationCycles" : integer`\
    The number of cycles for the optional pre-initialization stage described in
    [RASPA stages](#raspa-stages), which relaxes the configuration using only
    moves that keep the number of molecules fixed. Default: `0` (stage skipped).

### Restart and crash-recovery <a name="restart-and-crash-recovery"></a>

-   `"RestartFromBinaryFile" : boolean`\
    When `true`, `RASPA` resumes from a binary restart file that contains the
    complete state of the program, continuing from the point at which that file
    was written. The file can be large (up to several hundred megabytes) and is
    written every `"WriteBinaryRestartEvery"` cycles during a run. Default:
    `false`.

-   `"BinaryRestartFileName" : string`\
    The name of the binary restart file used by `"RestartFromBinaryFile"`.
    Default: `restart_data.bin`.

-   `"WriteBinaryRestartEvery" : integer`\
    How often (in cycles) the binary crash-recovery file is written. Default:
    `5000`.

-   `"RestartFileName" : string` *(per system)*\
    Reads the atomic positions of each component, and the simulation box, from a
    JSON restart file for that system. Any molecules requested with
    `"CreateNumberOfMolecules"` are created *in addition* to, and after, the
    positions read from this file. This is convenient for loading a fixed set of
    positions (for example cations) and then creating adsorbates on top of them.

### Printing options <a name="printing-options"></a>

-   `"PrintEvery" : integer`\
    Prints the loadings (when a framework is present) and energies every `int`
    cycles. For Molecular Dynamics, quantities such as energy conservation and
    the stress are also reported. Every status report starts with a progress
    line (wall time since the previous report, cycles per second, and the
    estimated time to the end of the current stage). Default: `5000`.

-   `"PrintMoveStatistics" : boolean`\
    Whether the periodic Monte Carlo status reports include a compact table of
    the move statistics: one line per move and component with the number of
    attempts and the acceptance since the previous report and cumulative, the
    current maximum change of the adaptive moves, and the share of the CPU
    time (the remainder, property sampling and energy/pressure computation, is
    listed as the last line). Widom insertions have no acceptance and show `-`.
    A move with at least 1000 attempts and no acceptance is flagged
    (`<-- never accepted`); a CBMC-type move whose trial growth fails for a
    large fraction of the attempts shows the constructed fraction. The
    acceptance since the previous report reveals a change of regime (a
    collapsing chain, a filling pore) that the cumulative acceptance averages
    out; the CPU share exposes moves that cost much and contribute nothing.
    Default: `true`.

### Parameter tuning <a name="parameter-tuning"></a>

-   `"RescaleWangLandauEvery" : integer`\
    How often (in cycles) the λ biasing factor is rescaled during the
    equilibration phase, for example in Continuous Fractional Component Monte
    Carlo. Default: `5000`.

-   `"OptimizeMCMovesEvery" : integer`\
    How often (in cycles) the maximum change of each Monte Carlo move is adjusted
    towards an optimal acceptance ratio (target: 0.5). The translation move tunes
    its maximum displacement, the rotation move its maximum angle, the hybrid MC
    move its time step, and so on. Default: `5000`.

-   `"ParallelTemperingSwapEvery" : integer`\
    For `"SimulationType" : "ParallelTempering"`,
    `"ParallelTemperingMolecularDynamics"`, `"HyperParallelTempering"`
    and `"ReweightedHistogram"`: how often (in cycles) a sweep of
    configuration swaps between replicas at neighboring state points is
    attempted (`0` disables the swaps). For the molecular-dynamics variant a
    cycle is one time step. Default: `10`.

-   `"SampleReweightingEvery" : integer`\
    For `"SimulationType" : "ReweightedHistogram"`: every this many production
    cycles each replica records a raw (N, U) sample for the reweighting
    analysis (controls the memory use and the sample correlation).
    Default: `5`.

-   `"ReweightingTemperatures" : [T_0, T_1, ...]`\
    For `"SimulationType" : "ReweightedHistogram"`: the temperatures (in
    Kelvin) the reweighted isotherms are evaluated and written at; they may
    lie in between the simulated temperatures. Default: the temperature
    ladder `"ExternalTemperatures"`.

-   `"ReweightingPressureRange" : [P_min, P_max]`\
    For `"SimulationType" : "ReweightedHistogram"` and `"ParallelTMMC"`: the
    pressure range (in Pascal) of the reweighted isotherms. Default: the span
    of the pressure ladder `"ExternalPressures"` (WHAM) or four decades around
    the reference pressure (TMMC). Written as the string `"auto"` — only with
    `"ComputeBET" : true` — the range is placed from a Widom Henry coefficient
    of the empty framework to nitrogen P0 (101325 Pa), the same rule
    the simulation BET path uses. With `"ComputeBET"` the key
    may also be omitted, which is the same as `"auto"`.

-   `"ReweightingNumberOfPressures" : integer`\
    For `"SimulationType" : "ReweightedHistogram"` and `"ParallelTMMC"`: the
    number of log-spaced pressures the reweighted isotherms are evaluated at.
    Default: `100`. When `"ComputeBET"` places the range automatically and this
    key is omitted, WHAM uses 400 points and TMMC uses 12 points per decade
    (at least 100).

-   `"WHAMTolerance" : floating-point-number`\
    For `"SimulationType" : "ReweightedHistogram"`: stop the WHAM free-energy
    iteration when the maximum absolute change in any state free energy \(g_i\)
    between successive iterates falls below this value (and likewise for the
    per-block solves used for error bars). Default: `1e-6`.

-   `"WHAMIterations" : integer`\
    For `"SimulationType" : "ReweightedHistogram"`: maximum number of WHAM
    free-energy iterations for the pooled solve and for each production-block
    solve. Default: `1000000`.

-   `"BETScoutMaximumCycles" : integer`\
    Cap on the unbiased GCMC occupancy scouts used when `"ComputeBET"` places
    the auto filling ceiling (`"MacroStateMaximumNumberOfMolecules": "auto"`)
    and when WHAM scouts Langmuir-coverage probes for the pressure ladder.
    Scouts still stop early once the loading plateaus (three successive
    1000-cycle blocks). Default: `15000`.

-   `"ComputeBET" : boolean`\
    For `"SimulationType" : "ReweightedHistogram"` and `"ParallelTMMC"`: after
    the reweighted isotherm is written, extract a BET surface area
    (Rouquerol consistency) using the adsorbate component's
    `"SaturationPressure"`, `"CrossSection"`, and `"LiquidVolume"`, and write
    it to the combined output. The same isotherm is also fit to the
    finite-layer (n-layer / BDDT) BET equation (\(n_m\), \(C\), integer \(n\)),
    reported next to the classical line in the text and JSON (`finiteLayer`).
    `"ComputeBTE"` is accepted as an alias. With this
    flag, `"ExternalPressures"` (WHAM) and `"ReweightingPressureRange"` may be
    `"auto"` or omitted: a Widom Henry coefficient of the empty framework
    places the bottom of the span at min(1 Pa, 10⁻³ n_Gurvich / K_H) and the
    top at P0. `"NumberOfThreads"` then sizes the WHAM pressure ladder
    (8–32 rungs), matching `raspa3 --threads`. WHAM then scouts occupancy
    at Langmuir θ = 0.1, 0.5, 0.9, fits Langmuir vs Langmuir–Freundlich vs Toth,
    and places rungs at equal Fisher overlap on a log-spaced skeleton (a point
    at least every two equal-log steps, at most two extras per interval). For WHAM,
    each converged production-block density of states also rebuilds the isotherm and
    refits slope/intercept inside the full-data Rouquerol window; the spread of
    those block areas is reported as the confidence-interval error on the BET
    area (and on \(n_m\) and \(C\)) when at least three blocks succeed. The
    finite-layer fit is jackknifed the same way with frozen \(n\) and window.
    TMMC does the
    same from the per-block collection-matrix increments (rebuilt `ln Π(N)`). TMMC still samples at
    `"ExternalPressure"` (typically P0); only the reweighting grid is placed
    automatically. `"MacroStateMaximumNumberOfMolecules"` may also be `"auto"`
    or omitted: an unbiased GCMC scout at P0 places the filling ceiling used
    as TMMC N_max and as the WHAM cycle length (length capped by
    `"BETScoutMaximumCycles"`). Default: `false`.

-   `"NumberOfWindows" : integer`\
    For `"SimulationType" : "ParallelTMMC"`: the number of contiguous
    macrostate windows the range `"MacroStateMinimumNumberOfMolecules"` to
    `"MacroStateMaximumNumberOfMolecules"` is split into (the windows share
    their endpoint macrostates). One walker (thread) is run per (temperature,
    window) pair. Default: `1`.

-   `"TMMCUpdateEvery" : integer`\
    For `"SimulationType" : "ParallelTMMC"`: the number of Monte Carlo steps
    between updates of the flattening transition-matrix bias (re-derived from
    the collection matrix). Default: `100000`.

### Threading and reproducibility <a name="threading-and-reproducibility"></a>

-   `"NumberOfThreads" : integer`\
    The number of worker threads. A value greater than 1 selects the thread-pool
    backend; otherwise the simulation runs serially. Default: `1`.

-   `"ThreadingType" : string`\
    Selects the threading backend explicitly. Either `"Serial"` or
    `"ThreadPool"`.

    For `"SimulationType" : "MolecularDynamicsSpatialDecomposition"`,
    `"NumberOfThreads"` is the number of sub-domains of the force engine
    (one persistent thread each); the thread pool is not used. The following
    keys tune that engine:

-   `"VerletSkin" : number`\
    For `"MolecularDynamicsSpatialDecomposition"`: the Verlet skin in
    &Aring; added to the largest cutoff when the neighbour lists are built.
    Lists are rebuilt when an atom has moved more than half the skin. A
    larger skin means fewer rebuilds but more pairs per force evaluation;
    values of 1–3 &Aring; are typical for a 1 fs time step. Default: `2.0`.

-   `"PPPMMeshSpacing" : number`\
    For `"MolecularDynamicsSpatialDecomposition"` with `"ChargeMethod" :
    "Ewald"`: the target spacing in &Aring; of the particle-mesh grid; the
    number of mesh points per cell vector is the smallest FFT-friendly
    integer (factors 2, 3, 5, 7) of at least the cell length divided by the
    spacing and at least twice the interpolation order. With the default
    order the relative error of the reciprocal energy and forces is about
    \f$10^{-4}\f$ at 1 &Aring; and \f$10^{-6}\f$ at 0.5 &Aring; (for
    \f$\alpha \approx 0.3\f$ &Aring;\f$^{-1}\f$). Default: `1.0`.

-   `"PPPMInterpolationOrder" : integer`\
    For `"MolecularDynamicsSpatialDecomposition"`: the order of the cardinal
    B-spline used to assign charges to the mesh and to interpolate the
    potential back (3–7; higher is more accurate per mesh point and costs
    order\f$^3\f$ mesh operations per atom). Default: `5`.

-   `"DomainGrid" : [integer, integer, integer]`\
    For `"MolecularDynamicsSpatialDecomposition"`: the sub-domain grid
    \f$p_x \times p_y \times p_z\f$ used by the threads; the product must
    equal `"NumberOfThreads"` and each \f$p_i\f$ must not exceed the number
    of neighbour-search cells along that axis (three per cutoff plus skin).
    The cut positions along each axis are balanced on the atom count at
    every list rebuild. Default: chosen automatically (the factorization
    with the smallest sub-domain surface).

-   `"RandomSeed" : integer`\
    Seeds the random-number generator for reproducible runs. When omitted a
    non-deterministic seed is used.

-   `"Units" : string`\
    Set to `"Reduced"` to run in reduced (dimensionless) Lennard-Jones units
    instead of the default physical units.

### Systems & Components <a name="systems-components"></a>

-   `"Systems" : list`\
    A list of system definitions, each a dictionary of the key-value pairs
    described in [System options](#system-options). A single process can own
    multiple systems; Gibbs-ensemble and parallel-tempering moves act on a pair
    of systems.

-   `"Components" : list`\
    A list of component (molecule) definitions, each a dictionary of the
    key-value pairs described in [Component options](#component-options).

----------------------------------------------------------------------------------

## System options <a name="system-options"></a>

### Operating conditions and thermostat/barostat-parameters <a name="operating-conditions-and-thermostatbarostat-parameters"></a>

-   `"ExternalTemperature" : floating-point-number`\
    The external temperature of the system in Kelvin. The inverse temperature
    β is derived from it and enters all Boltzmann statistics. This key is
    required for every system. Default: `300`.

-   `"ExternalTemperatures" : [T_0, T_1, ...]`\
    The temperature ladder for `"SimulationType" : "ParallelTempering"` and
    `"ParallelTemperingMolecularDynamics"` (a sorted list of at least two
    temperatures in Kelvin),
    `"HyperParallelTempering"`, `"ReweightedHistogram"` or `"ParallelTMMC"`
    (at least one). The single declared system is replicated into one replica
    per temperature (per (temperature, pressure) grid point for the
    grid-based types, per (temperature, window) pair for `"ParallelTMMC"`).
    Replaces `"ExternalTemperature"` for those simulation types.

-   `"ExternalPressure" : floating-point-number`\
    The external pressure of the system in Pascal.

-   `"ExternalPressures" : [P_0, P_1, ...]`\
    The pressure ladder for `"SimulationType" : "HyperParallelTempering"` or
    `"ReweightedHistogram"`: a sorted list of pressures in Pascal, converted
    to per-component fugacities internally with the Peng-Robinson equation of
    state at each grid point. Together with `"ExternalTemperatures"` it spans
    the (temperature, pressure) replica grid. Replaces `"ExternalPressure"`
    for those simulation types. Written as the string `"auto"` — only with
    `"SimulationType" : "ReweightedHistogram"` and `"ComputeBET" : true` — a
    log-spaced ladder is placed from the nitrogen Henry limit to P0 after the
    system is built; `"NumberOfThreads"` sets the rung count (8–32). Occupancy
    is then scouted at Langmuir θ = 0.1, 0.5, 0.9; Langmuir vs Langmuir–Freundlich
    vs Toth is selected, and the same count is re-placed at equal Fisher overlap
    on a log-spaced skeleton (gap at most twice equal-log spacing, at most two
    extra rungs per interval), keeping the Henry and P0 endpoints.
    With `"ComputeBET"` the key may also be omitted, which is the same as `"auto"`.

-   `"ExternalPressureX" / "ExternalPressureY" / "ExternalPressureZ" : floating-point-number`\
    Override individual diagonal components of the pressure tensor, for
    anisotropic (directional) pressure control. Each defaults to
    `"ExternalPressure"` when not given.

-   `"ChemicalPotential" : floating-point-number`\
    Sets the imposed chemical potential (in internal units); the corresponding
    fugacity is derived from it and the temperature.

-   `"MacroStateMinimumNumberOfMolecules" / "MacroStateMaximumNumberOfMolecules" : integer or "auto"`\
    For `"SimulationType" : "MonteCarloTransitionMatrix"` and
    `"ParallelTMMC"`: the macrostate range (the total molecule count) the
    transition-matrix walk is confined to. `"ReweightedHistogram"` uses the
    maximum as the number of Monte Carlo moves per cycle. For `"ParallelTMMC"`
    the range is split into `"NumberOfWindows"` windows, and a minimum of `0`
    enables the exact normalization of the saturation pressure by the
    empty-box state. With `"ComputeBET"`, `"MacroStateMaximumNumberOfMolecules"`
    may be `"auto"` or omitted: occupancy is scouted at nitrogen P0 and N_max
    is placed just above the plateau (capped by Gurvich packing; scout length
    by `"BETScoutMaximumCycles"`). Defaults: `0` and `100`.

-   `"MacroStateUseBias" : boolean`\
    For `"SimulationType" : "MonteCarloTransitionMatrix"`: whether the
    flattening transition-matrix bias is applied to the insertion/deletion
    acceptance (the collection-matrix statistics are unbiased either way).
    Default: `true`.

-   `"ThermostatChainLength" : integer`\
    The length of the Nosé-Hoover chain used to thermostat the system. Default:
    `5`.

-   `"NumberOfRespaSteps" : integer`\
    The number of RESPA substeps used by the thermostat and barostat chains.
    Default: `5`.

-   `"NumberOfYoshidaSuzukiSteps" : integer`\
    The number of Yoshida/Suzuki multiple-timestep integration steps. Default:
    `5`.

-   `"TimeScaleParameterThermostat" : floating-point-number`\
    The time scale on which the thermostat evolves. Default: `0.15`.

-   `"BarostatChainLength" : integer`\
    The length of the Nosé-Hoover chain coupled to the isotropic or cell
    barostat. Default: `5`.

-   `"TimeScaleParameterBarostat" : floating-point-number`\
    The pressure-coupling time scale $\tau_b$ in picoseconds. Default: `1.0`.
    The barostat mass follows Martyna, Tobias and Klein,
    $W = (N_f + 3)\,k_B T\,\tau_b^2$ for the strain $\epsilon = \ln(V)/3$
    (the cell-matrix mass of `NPTPR` is $W/3$), and the barostat chain masses
    are $k_B T\,\tau_b^2$. The integrator is the measure-preserving scheme of
    Tuckerman et al., J. Phys. A 39, 5629 (2006); the reported conserved
    quantity includes $W\dot\epsilon^2/2 + PV$ and the two chain energies.

### Box/Framework options <a name="boxframework-options"></a>

-   `"Type" : string`\
    Sets the system type:

    -   `"Box"`
        A simulation cell whose lengths and angles are specified directly.

    -   `"Framework"`
        A framework read from a `CIF`-file; the cell lengths and angles follow
        from that file.

-   `"BoxLengths" : [floating-point-number, floating-point-number, floating-point-number]`\
    The cell edge lengths of a `"Box"` system, in Ångström. Default:
    `[25, 25, 25]`.

-   `"BoxAngles" : [floating-point-number, floating-point-number, floating-point-number]`\
    The cell angles of a `"Box"` system, in degrees. Default: `[90, 90, 90]`.

-   `"Name" : string`\
    For `"Type" : "Framework"`, the name of the framework; it is used in output
    filenames and, when `"FileName"` is not given, the framework is loaded from
    the file `string.cif` (looked up in the working directory or `RASPA_DIR`).
    The name may not contain directory separators; use `"FileName"` to load
    from a path.

-   `"FileName" : string`\
    For `"Type" : "Framework"`, loads the framework from this `CIF`-file path
    (the `.cif` extension may be omitted). The framework name is derived from
    the file stem, unless `"Name"` overrides it.

-   `"NumberOfUnitCells" : [integer, integer, integer] or "auto"`\
    The number of unit cells in the `x`, `y`, and `z` directions. The super-cell
    contains these unit cells, and periodic boundary conditions are applied at
    the super-cell level (*not* at the unit-cell level). Default: `[1, 1, 1]`.
    With `"auto"`, the smallest integer replication in each direction is chosen
    so the super-cell satisfies the minimum-image convention for the force-field
    cut-offs (`ceil(2 * cutOff / perpendicularWidth)` per axis, at least 1).

-   `"HeliumVoidFraction" : floating-point-number`\
    The void fraction obtained by probing the structure with helium at room
    temperature. This value comes from a separate simulation and is required to
    compute the *excess* adsorption.

-   `"UseChargesFrom" : string`\
    Selects where framework charges are taken from:

    -   `"PseudoAtoms"`
        Uses the charges from the force-field definition file.

    -   `"CIF_File"`
        Uses the charges listed in the `CIF`-file via the `_atom_site_charge`
        tag. This allows an individual charge per framework atom, even for atoms
        of the same type.

    -   `"ChargeEquilibration"`
        Computes the framework charges with the charge-equilibration scheme of
        Wilmer and Snurr. The charges are symmetrized over symmetry-equivalent
        atoms.

### Force field definition <a name="force-field-definition"></a>

-   `"ForceField" : string`\
    Reads the force field from `string/force_field.json`. If a local
    `force_field.json` is present in the working directory it is used instead;
    otherwise the file is looked up under:

        ${RASPA_DIR}/simulations/share/raspa3/forcefield/string/force_field.json

### System `MC`-moves <a name="system-mc-moves"></a>

-   `"VolumeMoveProbability" : floating-point-number`\
    The probability per cycle of attempting a volume change. Rigid molecules are
    scaled by their center of mass, while flexible molecules and the framework
    are scaled atom by atom.

-   `"AnisotropicVolumeMoveProbability" : floating-point-number`\
    The probability per cycle of attempting an anisotropic volume change, in
    which the box edges are scaled independently.

-   `"GibbsVolumeMoveProbability" : floating-point-number`\
    The probability per cycle of attempting a Gibbs volume-change move in a Gibbs
    ensemble simulation. The total volume of the two boxes (typically a gas and a
    liquid phase) is kept constant while the individual box volumes change; the
    change is drawn randomly in $\ln(V_\mathrm{I}/V_\mathrm{II})$.

-   `"HybridMCProbability" : floating-point-number`\
    The probability per cycle of attempting a hybrid MC move. This move
    propagates the Hamiltonian through a short Molecular Dynamics trajectory and
    accepts or rejects the new state based on the energy drift. Rigid and
    flexible adsorbates are supported in a rigid or flexible framework. Use
    `"HybridMCMoveNumberOfSteps"` to set the number of MD steps.

-   `"ParallelTemperingSwapProbability" : floating-point-number`\
    The probability per cycle of attempting a parallel-tempering swap between two
    systems. Ignored with `"SimulationType" : "ParallelTempering"`,
    `"ParallelTemperingMolecularDynamics"`, `"HyperParallelTempering"` and
    `"ReweightedHistogram"`, where the swaps are performed by the driver at
    the barrier synchronization points (see `"ParallelTemperingSwapEvery"`).

-   `"TranslationSmartMCAllProbability" : floating-point-number`\
    The probability per cycle of attempting a translation smart-MC move that
    displaces all molecules simultaneously along the forces acting on them.
    Alias: `"ForceBiasTranslationAllProbability"`.

-   `"RotationSmartMCAllProbability" : floating-point-number`\
    The probability per cycle of attempting a rotation smart-MC move that
    rotates all rigid multi-atomic molecules simultaneously along the torques
    acting on them (quaternion update).

-   `"CrossLinkSwapProbability" : floating-point-number`\
    The probability per cycle of attempting a cross-link swap move (see
    [Cross-links between molecules](#cross-links)). An existing cross-link is
    picked at random together with one of its two ends as the pivot; the other
    end is detached and reattached to a free reactive site of a compatible type
    on a different molecule within the `"CaptureRadius"` of the pivot. The
    number of cross-links is conserved, so the move samples the topology of a
    network at fixed connectivity. Accepted with the Metropolis rule on the
    energy difference of the two bonded terms (the candidate set is state
    independent once the link is removed, so no extra bias factor is needed).

-   `"CrossLinkFormationProbability" : floating-point-number`\
    The probability per cycle of attempting a cross-link formation or scission
    move (each chosen with 50% probability). Formation picks a free reactive
    site and a free partner site of a compatible type within the
    `"CaptureRadius"`, and adds the link; scission removes a randomly chosen
    link. The acceptance rules contain the ratio of the number of free sites,
    the number of links and the number of partner candidates of both ends, so
    that the two moves are each other's reverse and sample the Boltzmann
    distribution over topologies, including the constant `"FormationEnergy"`
    of a link. A link is never formed when it would exceed the valence of one
    of its sites.

-   `"CrossLinkExchangeProbability" : floating-point-number`\
    The probability per cycle of attempting a cross-link exchange move, in
    which two links trade partners: $(a,b) + (c,d) \rightarrow (a,c) + (b,d)$.
    A link and one of its ends $a$ (the pivot) are picked at random; the
    exchange partner $c$ is drawn uniformly among the *linked* reactive sites of
    a compatible type on other molecules within the `"CaptureRadius"` of $a$
    that are not linked to $a$, and the link $(c,d)$ that $c$ gives up uniformly
    among the links of $c$. Neither the number of links nor the number of
    links of any site changes, so the valences are automatically respected and
    the move rewires a network whose sites are all saturated (where the swap
    move finds no free partner): the metathesis / transesterification-type
    exchange of vitrimers. Accepted with the Metropolis rule on the energy
    difference of the four bonded terms times the proposal factor
    $n_\mathrm{links}(c)/n_\mathrm{links}(b)$ (equal to one for valence-1
    sites). Rejected when $b$ lies outside the capture radius of $a$ or when
    the new link $(b,d)$ is not admissible (same molecule, no bond type, or
    already present).

### Cross-links between molecules <a name="cross-links"></a>

Cross-links are bonds *between* molecules. They are the building blocks for
reversibly cross-linked networks, associating polymers, vitrimers and gels:
the molecules keep their own (immutable) intra-molecular force field, and the
system owns a table of inter-molecular bonds that the three topology moves
above (swap, formation/scission, exchange) create, remove and rewire. Each molecule that can take part declares its
reactive atoms with `"ReactiveSites"` in its molecule definition file (see
[Component properties](#component-properties)); the bond potentials between
site types are declared per system with `"CrossLinkBonds"`.

The energy of a cross-linked state is the ordinary inter-molecular energy of
all molecules, corrected for the bonded pairs, plus the bonded terms:

$$U = U_\mathrm{inter} - \sum_{\mathrm{links}} u_\mathrm{pair}(r_{ij})
    + \sum_{\mathrm{links}} \left[ u_\mathrm{bond}(r_{ij}) + \epsilon_\mathrm{form}
    + \sum u_\mathrm{junction}(\theta) \right]$$

The van der Waals and real-space Coulomb interaction of the two bonded atoms is
subtracted (their exclusion is booked in the molecule-molecule VDW and Coulomb
slots and, for Ewald summation, in the Ewald exclusion slot, exactly as an
intra-molecular 1-2 exclusion would be), so that two linked monomers have the
same energy as one molecule with that bond. The bonded terms are reported in a
separate `cross-link` energy slot. The optional junction bends run over every
angle *neighbour–site–partner* formed by the link with the intra-molecular
neighbours of each site. Cross-links contribute to the gradients and the
virial, so they can be used with molecular dynamics, hybrid MC, volume moves
and Gibbs volume moves.

Linked molecules may be moved by all displacement moves (translation,
rotation, bead displacement, bead flip, crankshaft, pivot, concerted rotation,
smart MC, hybrid MC); the cross-link energy difference is included in the
acceptance. Reinsertion and partial reinsertion regrow a linked molecule with
its linked site atoms kept in place (a fixed-endpoint regrowth: the parts of
the molecule between two fixed sites are closed with the bridge-closure steps
of the CBMC engine); the link's bond, junction bends and exclusion corrections
enter the Rosenbluth weights as tethers to the frozen partner molecule, so the
regrowth samples the junction geometry exactly. Both the CBMC and the
recoil-growth chain scheme support these tethers. A rigid-body molecule with
two linked sites in one rigid fragment can not be regrown and is rejected when
the input is read; the regrowth plans of a reactive component with a
reinsertion move are built when the input is read (for up to twelve reactive
sites per molecule), so an impossible plan is reported before the simulation
starts. Moves that remove a molecule as a whole (deletion, tethered proton
hop) are skipped for molecules that currently carry a link; a molecule can
only leave the system after its links have been broken. Reactive components
can not use CFCMC-type moves,
identity changes, pair or group swaps, reptation, double bridging, the Gibbs
swap moves or reactions. Cross-links are written to and read from the restart
and crash-recovery files, and are swapped along with the configurations in
parallel tempering.

-   `"CrossLinkBonds" : list of objects`\
    The bond types that can be formed between reactive sites. Each object has

    -   `"Sites" : [string, string]`, the two reactive-site type names the bond
        connects (they may be equal, e.g. `["X", "X"]`);
    -   `"Bond" : ["POTENTIAL", [parameters...]]`, the bond potential and its
        parameters in the same form and units as the `"Bonds"` of a molecule
        definition file, e.g. `["HARMONIC", [100000.0, 1.54]]` (force constant
        in K/Å², equilibrium distance in Å). `"FIXED"` and `"NONE"` are not
        allowed;
    -   `"JunctionBend" : ["POTENTIAL", [parameters...]]` (optional), a bend
        potential applied to every angle formed by an intra-molecular neighbour
        of a site, the site, and its partner across the link, in the same form
        and units as the `"Bends"` of a molecule definition file, e.g.
        `["HARMONIC", [62500.0, 114.0]]` (K/rad², degrees);
    -   `"CaptureRadius" : floating-point-number` (default `2.0` Å), the
        maximum site–site distance at which the topology moves propose a new
        link. It only affects the proposal distribution (and therefore the
        efficiency), not the sampled distribution, but the bond potential should
        be able to bring pairs from the capture radius to the equilibrium
        distance within reasonable energies;
    -   `"FormationEnergy" : floating-point-number` (default `0.0`, in K), a
        constant energy added per link. A negative value favours bonded states
        and controls the degree of cross-linking (association constant) at
        equilibrium.

    Example:

        "CrossLinkBonds" : [
          {
            "Sites" : ["X", "X"],
            "Bond" : ["HARMONIC", [100.0, 3.8]],
            "JunctionBend" : ["HARMONIC", [100.0, 120.0]],
            "CaptureRadius" : 6.5,
            "FormationEnergy" : -2500.0
          }
        ]

-   `"InitialCrossLinks" : list of [[c, m, a], [c, m, a]]`\
    Links present at the start of the simulation, each given as two
    `[component, molecule, atom]` triples. The molecule indices refer to the
    molecules present at the start (those read from `"RestartFileName"`
    followed by those created with `"CreateNumberOfMolecules"`), the atoms must be
    reactive sites of a type for which a bond type exists, the two sites must
    belong to different molecules, and the valence of every site is respected.
    Without this key the simulation starts without cross-links and the
    formation/scission move builds them up.

### Molecular dynamics parameters <a name="molecular-dynamics-parameters"></a>

-   `"TimeStep" : floating-point-number`\
    The integration time step in picoseconds for `MD`. Default: `0.0005`.

-   `"HybridMCMoveNumberOfSteps" : integer`\
    The number of Molecular Dynamics steps used per hybrid MC move.

-   `"Ensemble" : string`\
    Sets the Molecular Dynamics ensemble:

    -   `"NVE"`\
        The micro-canonical ensemble: the number of particles $N$, the volume
        $V$, and the energy $E$ are constant.

    -   `"NVT"`\
        The canonical ensemble: the number of particles $N$, the volume $V$, and
        the average temperature $\left\langle T\right\rangle$ are constant, while
        the instantaneous temperature fluctuates. A Nosé-Hoover thermostat is
        attached.

    -   `"NPT"`\
        Isothermal-isobaric molecular dynamics with isotropic log-volume
        coupling. The cell shape is fixed and all three lengths scale together.
        The barostat couples to the centres of mass of rigid molecules and
        rigid groups and to every atom of a flexible molecule individually, and
        is driven by the virial and kinetic energy of exactly those points (for
        a flexible fluid: the atomic virial, bonded forces included, with
        $N_{\text{atoms}} k_B T$). The pressure reported in the output remains
        the molecular (centre-of-mass) pressure that the Monte Carlo volume
        move uses; the two are different estimators of the same average.

    -   `"NPTPR"`\
        Martyna-Parrinello-Rahman isothermal-isobaric dynamics with a flexible
        cell. `"CellType"` selects the constrained cell space:
        `"Regular"` (6 degrees of freedom), `"Monoclinic"` (4),
        `"Isotropic"` (1), `"Anisotropic"` (3),
        `"RegularUpperTriangle"`/`"REGULAR_UPPER_TRIANGLE"` (6), or
        `"MonoclinicUpperTriangle"`/`"MONOCLINIC_UPPER_TRIANGLE"` (4).
        `"MonoclinicAngleType"` selects `"Alpha"`, `"Beta"` (default), or
        `"Gamma"` for the single shear degree of freedom. Upper-triangular modes
        preserve forbidden lower-triangle entries exactly.

    -   `"MuVT"`\
        Grand-canonical molecular dynamics at fixed volume. RASPA attempts one
        configurational-bias insertion or deletion every three MD steps and
        otherwise uses the NVT integrator.

    -   `"MuPT"`\
        Osmotic molecular dynamics with grand-canonical particle exchange and
        the isotropic NPT pressure controller.

    -   `"MuPTPR"`\
        Osmotic molecular dynamics with grand-canonical particle exchange and
        the flexible-cell NPTPR pressure controller.

    NPT, NPTPR, MuPT, and MuPTPR require `"ExternalPressure"` and use hydrostatic pressure
    coupling. Their reported conserved quantity is the extended enthalpy:
    physical energy plus thermostat energy, pressure-volume work, cell/barostat
    kinetic energy, and barostat-chain energy. Complete Nosé-Hoover
    thermobarostat trajectories are tested for time reversibility and bounded
    extended-enthalpy drift; canonical symplecticity applies only to isolated
    Hamiltonian submaps.

    MuVT, MuPT, and MuPTPR also require `"ExternalPressure"` as the reservoir
    pressure and `"SwapProbability"` greater than zero for at least one
    component. The component fugacity is computed from reservoir pressure,
    mole fraction, and `"FugacityCoefficient"`. Accepted insertions receive
    Maxwell-Boltzmann velocities, and thermostat/barostat masses are refreshed
    for the new number of degrees of freedom. Because particle exchange is a
    stochastic Monte Carlo step, conserved-energy drift is meaningful only
    between accepted exchanges.

-   `"ComputeElasticConstantsFromFluctuations" : boolean`\
    Enables an isothermal elastic-tensor calculation during fixed-cell NVT
    molecular dynamics. Each observation accumulates the instantaneous affine
    Born tensor and configurational stress covariance. The molecular kinetic
    contribution is added consistently with the active translational
    constraints. Output contains the separate Born, kinetic, covariance, and
    hydrostatic prestress terms, the Helmholtz and tangent stiffness matrices,
    block confidence intervals, stability eigenvalues, and derived moduli.
    Voigt order is `xx, yy, zz, yz, xz, xy`, with engineering shear strains.
    External fields and polarization are not currently supported.

-   `"ElasticConstantsSampleEvery" : integer`\
    Number of MD steps between elastic observations (default: `100`). Born
    Hessians are substantially more expensive than ordinary force evaluations.
    Choose the interval at least as long as the short-time stress correlation,
    use production blocks much longer than the integrated stress
    autocorrelation time, and check convergence against trajectory length,
    block length, sampling interval, and system size. The NVT reference cell
    should first be equilibrated at the desired temperature and mean pressure.

### Options to measure properties <a name="options-to-measure-properties"></a>

#### Output pdb-movies <a name="output-pdb-movies"></a>

`"OutputPDBMovie" : boolean`

Whether to write simulation snapshots as PDB movies. Output is written to the
directory `movies`.

-   `"SampleMovieEvery" : integer`\
    Write a snapshot every `int` cycles. Default: `1`.

-   `"RestrictMoviePositionsToBox" : boolean`\
    Whether to wrap the written positions back into the simulation box. Default:
    `true`.

#### Histogram of the energy <a name="histogram-of-the-energy"></a>

`"ComputeEnergyHistogram" : boolean`

Whether to accumulate a histogram of the energy for the system. During
adsorption, for example, it tracks the total energy together with the Van der
Waals, Coulombic, and polarization contributions. Output is written to the
directory `energy_histogram`.

-   `"SampleEnergyHistogramEvery" : integer`\
    Sample the energy histogram every `int` cycles. Default: `1`.

-   `"WriteEnergyHistogramEvery" : integer`\
    Write the energy histogram every `int` cycles. Default: `5000`.

-   `"NumberOfBinsEnergyHistogram" : integer`\
    The number of bins in the histogram. Default: `128`.

-   `"LowerLimitEnergyHistogram" : floating-point-number`\
    The lower bound of the histogram. Default: `-5000`.

-   `"UpperLimitEnergyHistogram" : floating-point-number`\
    The upper bound of the histogram. Default: `1000`.

#### Histogram of the number of molecules <a name="histogram-of-the-number-of-molecules"></a>

`"ComputeNumberOfMoleculesHistogram" : boolean`

Whether to accumulate histograms of the number of molecules for the system. In
open ensembles the number of molecules fluctuates. Output is written to the
directory `number_of_molecules_histogram`.

-   `"SampleNumberOfMoleculesHistogramEvery" : integer`\
    Sample the histogram every `int` cycles. Default: `1`.

-   `"WriteNumberOfMoleculesHistogramEvery" : integer`\
    Write the histogram every `int` cycles. Default: `5000`.

-   `"LowerLimitNumberOfMoleculesHistogram" : integer`\
    The lower bound of the histograms. Default: `0`.

-   `"UpperLimitNumberOfMoleculesHistogram" : integer`\
    The upper bound of the histograms. Default: `200`.

#### Radial Distribution Function (RDF) force-based <a name="radial-distribution-function-rdf-force-based"></a>

`"ComputeRDF" : boolean`

Whether to compute the force-based (Borgis) radial distribution function using
site gradients. Output is written to the directory `rdf`.

In molecular dynamics the integrator forces are reused. In Monte Carlo a full
site-gradient evaluation is performed when sampling (framework + intermolecular +
Ewald + intramolecular), so flexible molecules are handled correctly. This is
independent of the molecular-pressure gradient path, which omits intramolecular
forces for the atomic-to-molecular virial correction.

-   `"SampleRDFEvery" : integer`\
    Sample the RDF every `int` cycles. Default: `10`.

-   `"WriteRDFEvery" : integer`\
    Write the RDF every `int` cycles. Default: `5000`.

-   `"NumberOfBinsRDF" : integer`\
    The number of bins in the RDF. Default: `128`.

-   `"UpperLimitRDF" : floating-point-number`\
    The upper distance limit of the RDF, in Ångström. Default: `15.0`.

#### Radial Distribution Function (RDF) conventional <a name="radial-distribution-function-rdf-conventional"></a>

`"ComputeConventionalRDF" : boolean`

Whether to compute the conventional (histogram-based) radial distribution
function. Output is written to the directory `conventional_rdf`.

-   `"SampleConventionalRDFEvery" : integer`\
    Sample the RDF every `int` cycles. Default: `10`.

-   `"WriteConventionalRDFEvery" : integer`\
    Write the RDF every `int` cycles. Default: `5000`.

-   `"NumberOfBinsConventionalRDF" : integer`\
    The number of bins in the RDF. Default: `128`.

-   `"RangeConventionalRDF" : floating-point-number`\
    The upper distance limit of the RDF, in Ångström. Default: `15.0`.

#### Mean-Squared Displacement (MSD) order-N <a name="mean-squared-displacement-msd-order-n"></a>

`"ComputeMSD" : boolean`

Whether to compute the mean-squared displacement (MSD) using the order-N
algorithm, from which self-diffusion coefficients can be obtained. Output is
written to the directory `msd`. Besides the self-MSD per component, the
collective (Onsager) MSDs are written per component pair, normalized by the
total number of molecules \(N\) so that
\(\text{MSD}_{ij} = \langle \Delta \mathbf{R}_i \cdot \Delta \mathbf{R}_j \rangle / N\)
is symmetric in the components. Computing the MSD requires a fixed number of
molecules; do not combine it with insertion/deletion moves.

-   `"SampleMSDEvery" : integer`\
    Sample the MSD every `int` cycles. Default: `10`.

-   `"WriteMSDEvery" : integer`\
    Write the MSD every `int` cycles. Default: `5000`.

-   `"NumberOfBlockElementsMSD" : integer`\
    The number of elements per block in the order-N scheme. Default: `25`.

#### Velocity Auto-Correlation Function (VACF) <a name="velocity-auto-correlation-function-vacf"></a>

`"ComputeVACF" : boolean`

Whether to compute the velocity auto-correlation function (VACF) using
multiple staggered buffers, from which self-diffusion coefficients can be
obtained via the Green-Kubo relation. Output is written to the directory
`vacf`. Besides the self-VACF per component, the collective (Onsager) VACFs
are written per component pair, normalized by the total number of molecules
\(N\) so that
\(\text{VACF}_{ij} = \langle \mathbf{V}_i(t) \cdot \mathbf{V}_j(0) \rangle / N\)
is symmetric in the components. Computing the VACF requires a fixed number of
molecules; do not combine it with insertion/deletion moves.

-   `"SampleVACFEvery" : integer`\
    Sample the VACF every `int` cycles. Default: `10`.

-   `"WriteVACFEvery" : integer`\
    Write the VACF every `int` cycles. Default: `5000`.

-   `"NumberOfBuffersVACF" : integer`\
    The number of staggered buffers (time origins in use at any moment).
    Default: `20`.

-   `"BufferLengthVACF" : integer`\
    The length of each buffer, i.e. the number of correlation times.
    Default: `1000`.

#### Density grids <a name="density-grids"></a>

`"ComputeDensityGrid" : boolean`

Whether to compute three-dimensional density grids. Output is written to the
directory `density_grids`.

-   `"SampleDensityGridEvery" : integer`\
    Sample the density grids every `int` cycles. Default: `10`.

-   `"WriteDensityGridEvery" : integer`\
    Write the density grids every `int` cycles. Default: `5000`.

-   `"DensityGridSize" : [integer, integer, integer]`\
    The number of voxels along each axis. Default: `[128, 128, 128]`.

-   `"DensityGridNormalization" : string`\
    How the grid values are normalized: `"Max"` (default, scaled to the maximum)
    or `"NumberDensity"`.

-   `"DensityGridBinning" : string`\
    The strategy used to accumulate the density grids:

    -   `"Standard"`\
        Conventional histogram binning: each particle contributes fully to the
        voxel it resides in. Default: `"Standard"`.

    -   `"Equitable"`\
        Each particle contributes fractionally to neighboring voxels based on its
        position. This produces smoother grids and reduces discretization
        artifacts, especially for fine grids.

-   `"DensityGridPseudoAtomsList" : [string, string, ...]`\
    Restricts the density grid to a subset of pseudo-atoms of a component. When
    given, a separate grid is produced for each listed pseudo-atom type instead
    of a single combined grid. When omitted, all pseudo-atoms of the component
    are accumulated into one grid. This is useful for resolving atom-specific
    adsorption within a molecule, for example separating the carbon and oxygen
    sites of CO<sub>2</sub>.

----------------------------------------------------------------------------------

## Force field options <a name="force-field-options"></a>

The following keywords control the force field. `"MixingRule"`,
`"TruncationMethod"`, `"TailCorrections"`, `"PseudoAtoms"`, `"SelfInteractions"`,
and `"BinaryInteractions"` are read from the force field file
(`force_field.json`), while the cutoffs and `"ChargeMethod"` are set per system.

-   `"MixingRule" : string`
    -   `"Lorentz-Berthelot"`
        The geometric mean for the strength parameter and the arithmetic mean for
        the size parameter. For Lennard-Jones:
        \begin{equation}
        \varepsilon_{ij}=\sqrt{\varepsilon_i \varepsilon_j}
        \end{equation}
        \begin{equation}
        \sigma_{ij}=\frac{\sigma_i+\sigma_j}{2}
        \end{equation}

    -   `"Jorgensen"`
        The geometric mean for both parameters. For Lennard-Jones:
        \begin{equation}
        \varepsilon_{ij}=\sqrt{\varepsilon_i \varepsilon_j}
        \end{equation}
        \begin{equation}
        \sigma_{ij}=\sqrt{\sigma_i \sigma_j}
        \end{equation}

-   `"TruncationMethod" : string`
    -   `"truncated"`
        Truncates the potential at the cutoff.
    -   `"shifted"`
        Truncates the potential at the cutoff and shifts it so that the potential
        energy is zero at the cutoff radius.

-   `"TailCorrections" : boolean`\
    Whether to apply analytic tail corrections for the truncated Van der Waals
    potential.

-   `"CutOffVDW" : floating-point-number`\
    The cutoff of the Van der Waals potential (both framework-molecule and
    molecule-molecule interactions). Interactions beyond this distance are
    omitted from the energy and force evaluation.

-   `"CutOffCoulomb" : floating-point-number`\
    The cutoff of the charge-charge potential, which is truncated at the cutoff.
    Tail corrections are not applied; the long-range part is instead recovered
    with the Ewald summation (`"ChargeMethod" : "Ewald"`). Together with the
    Ewald precision, this cutoff also determines the number of wave vectors and
    the Ewald parameter α. For large unit cells a Coulomb cutoff of about half
    the shortest box length avoids an excessive number of wave vectors. For
    non-Ewald calculations the cutoff should be as large as possible (greater
    than about 30 Å).

-   `"CutOff" : floating-point-number`\
    A convenience key that sets both `"CutOffVDW"` cutoffs (framework-molecule
    and molecule-molecule) at once.

-   `"OmitEwaldFourier" : boolean`\
    Skips the Fourier (reciprocal-space) part of the Ewald summation. Intended
    for testing only.

-   `"ComputePolarization" : boolean`\
    Whether to include polarization (induced-dipole) energy in the interactions.

-   `"ChargeMethod" : string`\
    Sets the method used for the electrostatics:

    -   `"None"`
        Skips the entire charge calculation. Use only when none of the species
        carry a charge.

    -   `"Ewald"`
        Uses the Ewald summation for the charge calculation.

-   `"PseudoAtoms" : list` <br>
    A list of pseudo-atoms, each with
    - `"name" : string`
    - `"framework" : boolean`
    - `"print_to_output" : boolean`
    - `"element" : string`
    - `"print_as" : string`
    - `"mass" : floating-point-number`
    - `"charge" : floating-point-number`
    - `"source" : string`

-   `"SelfInteractions" : list` <br>
    A list of self-interactions, each with
    - `"name" : string`
    - `"type" : string`
    - `"parameters" : [floating-point-number]`
    - `"source" : string`

-   `"BinaryInteractions" : []` <br>
    A list of binary interactions, each with
    - `"names" : [string, string]`
    - `"type" : string`
    - `"parameters" : [floating-point-number]`
    - `"source" : string`

### Configurational-bias and recoil-growth options <a name="cbmc-options"></a>

The following keywords, read from the force field file (`force_field.json`),
control how molecules are grown in the CBMC-based moves (insertion, deletion,
reinsertion, partial reinsertion, reptation, identity change, Gibbs and CFCMC
swaps, Widom insertion). A molecule is grown along a deterministic *growth
plan* over its fragment graph: the first bead is placed with the
multiple-first-bead scheme, flexible beads are attached one branch point at a
time with an exact bond/bend sampler and a Rosenbluth-selected torsion spin,
rigid bodies (`"RigidBodies"` in the molecule file) are hinged as one unit,
and rings are closed as one cluster with an internal Monte-Carlo. All options
affect sampling efficiency only; every value yields the same Boltzmann
distribution.

The torsion-spin selection of a flexible attach step is *guided*: the trial
spins are weighted not only by the bonded and short-range intramolecular terms
that are known at that step, but also by a tabulated *lookahead* factor
`g(φ)`, the Boltzmann average over the ideal-chain conformations of the beads
of the next two bonds of the interactions that couple those future beads to
the already placed ones (the 1-5 van der Waals and Coulomb pairs across the
junction). This steers the spin away from torsion states that the bare torsion
potential favours but that a bead placed two steps later cannot accommodate
(e.g. the cis state of an alkyl-ester C-C-O-C(=O) dihedral, whose carbonyl
oxygen would then clash with the chain). The guide is divided out of the
Rosenbluth weight of the selected spin, so the sampled distribution is exact
for any `g`; only the variance of the growth weights (and hence the acceptance
of the CBMC moves and the quality of `"CreateNumberOfMolecules"` growths)
improves. The tables are built once per distinct junction environment and
temperature at start-up (typically a few seconds for a molecule with tens of
beads) and require no input. The guide is local: it does not know about
interactions between distant parts of the same molecule, so the acceptance of
regrowing a molecule that is collapsed on itself (a long chain in vacuum) is
not improved by it.

-   `"NumberOfFirstBeadPositions" : integer`\
    Number of trial positions of the first bead, drawn uniformly in the box
    (default: 10).

-   `"NumberOfTrialDirections" : integer`\
    Number of trial directions `k` per growth step of the configurational-bias
    scheme; the Rosenbluth weight of a step is the sum of their Boltzmann
    factors divided by `k` (default: 10).

-   `"NumberOfTorsionTrialDirections" : integer`\
    Number of trial spins about the junction bond among which the torsion
    orientation of a step is Rosenbluth-selected, per trial direction
    (default: 100). Every trial direction of a step costs this many bonded
    energy evaluations; 10 to 20 is usually sufficient for a single torsion.

-   `"NumberOfTrialMovesPerOpenBead" : integer`\
    Number of internal Metropolis moves per placed bead used to relax the
    orientation of a hinged rigid body and the conformation of a ring cluster
    (default: 150). These moves carry no Rosenbluth weight.

-   `"CBMCRingCrankshaftProbability" : floating-point-number`\
    Probability, per internal move of a ring cluster, of attempting a large-angle
    crankshaft rotation of one ring atom about its two neighbours (the move that
    hops between ring conformers, e.g. chair and twist-boat) instead of a local
    displacement (default: 0.2).

-   `"CBMCRingTiltProbability" : floating-point-number`\
    Probability, per internal move of a ring cluster with a junction, of tilting
    the whole ring about its anchor instead of moving a single unit
    (default: 0.25).

-   `"UseDualCutOff" : boolean`\
    Grow and retrace with a short inner cut-off (`"DualCutOff"`) for all
    framework-molecule, molecule-molecule, and Coulomb interactions, and correct
    the Rosenbluth weights and energies of the grown and retraced configurations
    to the full cut-offs afterwards (Vlugt et al.). The acceptance rule is exact;
    the gain is that every trial direction is evaluated with the cheap short
    cut-off and only the selected configuration with the full one
    (default: `false`).

-   `"DualCutOff" : floating-point-number`\
    The inner cut-off of the dual cut-off scheme, in Å. Must be smaller than
    every full cut-off (default: 6.0).

-   `"UseRecoilGrowth" : boolean`\
    Grow and retrace the chain beyond the first bead with recoil growth
    (Consta, Vlugt, Wichers Hoeth, Smit, and Frenkel, *Mol. Phys.* **97**, 1243
    (1999)) instead of configurational bias (default: `false`). At every step
    `k` trial directions are generated; a direction is *open* with probability
    `min(1, exp(-β(u - u_ref)))`, with `u_ref` a fixed per-step reference energy:
    the maximum intramolecular strain of that step over 50 equilibrated
    ideal-gas conformations of the molecule, grown once at setup with a fixed
    seed (so a molecule with intrinsic non-bonded strain is not penalised, and
    the reference is a constant of the run). A direction is *available* when a
    *feeler* of `l - 1` further steps can be grown from it. The growth backtracks
    (recoils) over at most `l` steps when it dead-ends. Recoil growth is more
    efficient than configurational bias for long chains in dense or strongly
    confining environments, where a configurational-bias grow commits to a
    direction that has no future. It applies to every CBMC-based move. The
    component statistics report, per component, how many recoil grows
    completed, dead-ended (recoiled all the way back; more directions or a
    longer recoil helps), or were discarded (dead-ended past a committed
    direction; a longer recoil helps), and the mean fraction `m_i/k` of
    available directions per step (close to 1: an open environment, close to
    `1/k`: the weight is dominated by single available directions). Two
    caveats:
    - The recoil-growth weight is a valid factor of a Metropolis acceptance
      ratio but it is not the Rosenbluth weight whose average is the Widom
      estimator of the excess chemical potential; Widom insertion therefore
      always uses configurational bias, whatever this option is.
    - The retrace divides by the openness probability of the existing
      configuration; in crowded, repulsive environments this gives the weight a
      larger variance than configurational bias.

-   `"RecoilGrowthNumberOfTrialDirections" : integer`\
    The number of trial directions `k` per step of recoil growth (default: 5).

-   `"RecoilGrowthMaximumRecoilLength" : integer`\
    The recoil length `l` (default: 2). The feelers are exhaustive searches, so
    the cost per step scales as `k^l`; `l = 2` is usually sufficient and a
    warning is printed for `l ≥ 3`.

----------------------------------------------------------------------------------

## Component options <a name="component-options"></a>

### Component properties <a name="component-properties"></a>

-   `"Name" : string`\
    The name of the component; it is used in output filenames and must be
    unique among the components. When `"FileName"` is not given, it is also the
    base name of the definition file `Name.json` (looked up in the working
    directory or `RASPA_DIR`). The name may not contain directory separators;
    use `"FileName"` to load from a path.

-   `"FileName" : string`\
    Loads the component definition from this file path (the `.json` extension
    may be omitted). The component name is derived from the file stem, unless
    `"Name"` overrides it — useful when two components share one definition
    file.

-   `"Type" : string`\
    The component type: `"Adsorbate"` (default) or `"Cation"`.

-   `"MolFraction" : floating-point-number`\
    The mole fraction of this component in the mixture. Values may be given
    relative to the other components, as the fractions are normalized afterwards.
    Per-component partial pressures follow from the total pressure and the mole
    fractions.

-   `"FugacityCoefficient" : floating-point-number`\
    The fugacity coefficient of the component. When set to 0 (or omitted), the
    fugacity coefficient is computed automatically from the Peng-Robinson
    equation of state; this requires the critical pressure, critical temperature,
    and acentric factor to be present in the molecule file.

-   `"IdealGasRosenbluthWeight" : floating-point-number`\
    The ideal-gas Rosenbluth weight, i.e. the `CBMC` growth factor of a single
    chain in an empty box. It depends only on temperature and therefore needs to
    be computed once. Supplying it in advance is convenient for adsorption,
    because the applied pressure then needs no correction afterwards (the
    Rosenbluth weight shifts the chemical-potential reference, and the chemical
    potential follows directly from the fugacity). For equimolar mixtures this is
    essential.

-   `"CrossSection" : floating-point-number`\
    Probe cross-section σ for BET surface-area conversion [Å²]. Required when
    `"ComputeBET"` is true (nitrogen is usually 16.2 Å²).

-   `"LiquidVolume" : floating-point-number`\
    Liquid molecular volume v_L for Gurvich packing and the t-plot [Å³ per
    molecule]. Required when `"ComputeBET"` is true (nitrogen is usually
    57.7 Å³).

-   `"SaturationPressure" : floating-point-number`\
    Saturation pressure P0 used for relative pressure x = P/P0 in the BET fit
    [Pa]. Required when `"ComputeBET"` is true (nitrogen at 77 K is usually
    101325 Pa).

-   `"CreateNumberOfMolecules" : integer`\
    The number of molecules to create for this component at start-up. These
    molecules are created *in addition* to anything read from a restart file, so
    when restarting this value is usually set back to zero. Setting it
    unreasonably high can cause an infinite loop: the routine only accepts
    molecules whose growth causes no overlap (energy below the overlap
    criterion). A flexible molecule is not taken from a single growth: the first
    valid growth starts a short chain of further growths, each replacing the
    current one with the reinsertion acceptance min(1, W_new/W_old), which
    removes the conformations a single growth over-represents (a torsion chosen
    before its 1-5 partners exist, such as a twisted or E ester at the end of a
    diacrylate that molecular dynamics could not undo). In a dense fluid the
    weights are dominated by the fit into the surroundings, so there the
    initialization cycles (partial reinsertion with fixed endpoints) remain
    necessary to finish the job. The starting configurations are far from
    optimal, so substantial equilibration is needed to relax the energy; the
    `CBMC` growth can, however, reach very high densities.

-   `"StartingBead" : integer`\
    The index of the bead from which `CBMC` growth starts. Must be smaller than
    the number of atoms in the molecule.

-   `"BlockingPockets" : [[3 x floating-point-number, floating-point-number]] or string`\
    Blocks certain pockets of the simulation volume so molecules cannot grow
    into them. A typical example is the sodalite cages in FAU- and LTA-type
    zeolites, which are inaccessible to methane and larger molecules.

    Written as a list, each pocket is four numbers: the fractional positions
    $s_x$, $s_y$, $s_z$ and a radius in Ångström. For example, the blocking
    pockets of ITQ-29 for small molecules are:

        "BlockingPockets" : [
                   [0.0,       0.0,        0.0,       4.0],
                   [0.5,       0.0,        0.0,       0.5],
                   [0.0,       0.5,        0.0,       0.5],
                   [0.0,       0.0,        0.5,       0.5]
                 ]

    Written as the string `"auto"`, the pockets are computed from the framework
    read from the CIF-file, before the first molecule is placed:

        "BlockingPockets" : "auto"

    The void is split into what a helium probe can reach and what it cannot,
    and each unreachable cavity is covered by a sphere at its centre, of the
    lesser of the radius that holds the cavity and the radius past which the
    sphere would reach a channel. The spheres are a property of the framework
    and of the probe, so every component asking for them gets the same ones,
    and the framework is measured once. Silicalite comes back with none, being
    all channel; KFI comes back with eight, its two *lta* cages and its six
    *pau* cages. This needs the force field to give the framework atoms a van
    der Waals size of their own: a force field that leaves them without
    self-interactions and instead names every framework-guest pair outright
    describes a framework of points, and the run stops and says so rather than
    reporting that nothing is blocked.

    Written as any other string, it is the name of a `.block` file to read the
    pockets from, in the format the structural analysis writes: the number of
    spheres on the first line, then one line of $s_x$, $s_y$, $s_z$ and a radius
    per sphere, with no comments. The `.block` extension is added when absent
    and the file is looked for in the working directory and then in `RASPA_DIR`:

        "BlockingPockets" : "ITQ-29"

    All three forms may be given in the molecule definition file as well, and
    what the two files say is added together.

-   `"LambdaBiasFileName" : string`\
    Points to a JSON file of preset λ values, allowing optimized CFCMC
    simulations to run without re-estimating the biasing weights with
    Wang-Landau.

-   `"ThermodynamicIntegration" : boolean or string`\
    Enables thermodynamic integration of dU/dλ for the fractional molecule. As a
    boolean, `true` integrates the default (grand-canonical) λ. As a string it
    selects which λ coordinate to follow: `"CFCMC"` (default), `"CFCMC_PairSwap"`,
    or `"CFCMC_CBMC_PairSwap"`. With
    `"SimulationType" : "ParallelThermodynamicIntegration"`, `true` marks the
    component whose λ is pinned per replica (the λ coordinate is inferred from
    the component definition as for `"LambdaBinIndex"`).

-   `"LambdaBinIndex" : integer`\
    Used with `"SimulationType" : "ThermodynamicIntegration"`. Creates
    fractional molecule(s) pinned at the fixed lambda-bin
    λ = binIndex / (`NumberOfLambdaBins` − 1). The value must be smaller than
    `NumberOfLambdaBins`, and cannot be combined with λ-changing CFCMC moves.
    When the component defines `"GroupComponents"` the group-swap λ is pinned
    and the whole group (central component plus satellites) becomes fractional;
    when it defines `"PairComponent"` the ion-pair λ is pinned and both
    components of the pair become fractional (set `"LambdaBinIndex"` on the
    lowest-index component of the pair); otherwise the grand-canonical λ is
    pinned with a single fractional molecule. During production ⟨∂U/∂λ⟩ is
    sampled at this λ only; the average and its block-error estimate are
    written to the text and JSON output
    (`"properties" > "thermodynamicIntegration"`), giving one point of the
    ⟨∂U/∂λ⟩(λ) curve.

-   `"ReactiveSites" : list` (molecule definition file)\
    Declares which atoms of the molecule can form cross-links with atoms of
    other molecules (see [Cross-links between molecules](#cross-links)). Each
    entry is `[atom, "type"]` or `[atom, "type", valence]`, or equivalently an
    object `{"Atom" : atom, "Type" : "type", "Valence" : valence}`. The atom
    index refers to the `"pseudoAtoms"` list of the molecule, the type is a
    free name matched against the `"Sites"` of the system's
    `"CrossLinkBonds"`, and the valence (default `1`) is the maximum number of
    simultaneous cross-links the site can carry. An atom may be listed once.
    Example for a telechelic chain whose two end beads can each form one bond:

        "ReactiveSites" : [
          [0, "X"],
          [3, "X"]
        ]

    A component with reactive sites can not have a fractional molecule and
    can not use identity-change, pair-swap, group-swap, reptation, double
    bridging or Gibbs swap moves.

-   `"LnPartitionFunction" : number or string`\
    The natural logarithm of the (reduced) partition function used for reactions.
    Give a number to set it directly, or a species name (or `"auto"`, which uses
    the component name) to look it up in the embedded thermochemical database.
    The lookup is evaluated at each system's `"ExternalTemperature"` and uses the
    database selected by `"ThermochemicalDatabase"`.

### Component analyses <a name="component-analyses"></a>

The intra-molecular analyses below are switched on per component in its entry of
the `"Components"` list. Each has its own sampling and writing schedule and is
accumulated in every system that contains the component. The accumulated data of
these analyses is included in the binary restart files, so a continued run picks
up where it left off.

#### Histograms of the intra-molecular geometry <a name="histograms-of-the-intra-molecular-geometry"></a>

`"ComputeMoleculeProperties" : boolean`

Whether to accumulate probability histograms of the intra-molecular geometry
(bond lengths, bend angles, and torsion/dihedral angles) of this
component. This mirrors the "molecule properties" analysis from RASPA2. Output
is written to the directory `molecule_properties`.

-   `"SampleMoleculePropertiesEvery" : integer`\
    Sample the histograms every `int` cycles. Default: `10`.

-   `"WriteMoleculePropertiesEvery" : integer`\
    Write the histograms every `int` cycles. Default: `5000`.

-   `"NumberOfBinsMoleculeProperties" : integer`\
    The number of bins of the bond, bend and torsion histograms (fixed ranges).
    Default: `128`.

-   `"BondRangeMoleculeProperties" : floating-point-number`\
    The upper bound of the bond-length histogram, in Ångström. The bend range is
    fixed to `[0, 180]` degrees and the torsion range to `[-180, 180]` degrees.
    Default: `4.0`.

-   `"EndToEndRangeMoleculeProperties" : floating-point-number`\
    The upper bound of the end-to-end distance histogram, in Ångström (sampled
    when the component has end-to-end atoms, see `EndToEndAtoms`). By default
    1.05 times the contour length between the two ends, so every reachable
    conformation fits.

-   `"BinWidthEndToEndMoleculeProperties" : floating-point-number`\
    The bin width of the end-to-end distance histogram, in Ångström. Because the
    range grows with the chain length, this histogram is sized by its bin width
    rather than a bin count: the number of bins is the range divided by the width,
    rounded up (the last bin closes at or just beyond the range). Default: `0.25`.

#### Molecule shape: the gyration-tensor family <a name="molecule-shape"></a>

`"ComputeMoleculeShape" : boolean`

Whether to sample the gyration tensor
\f$S_{\alpha\beta} = \sum_i w_i (r_{i\alpha}-r_{\mathrm{cm},\alpha})(r_{i\beta}-r_{\mathrm{cm},\beta})\f$
of every molecule of this component (which needs at least two atoms), and the shape
descriptors that follow from its eigenvalues \f$\lambda_1 \ge \lambda_2 \ge \lambda_3\f$:

-   the radius of gyration \f$R_g^2 = \lambda_1+\lambda_2+\lambda_3\f$,
-   the asphericity \f$b = \lambda_1 - \tfrac{1}{2}(\lambda_2+\lambda_3)\f$ and the
    acylindricity \f$c = \lambda_2-\lambda_3\f$,
-   the relative shape anisotropy \f$\kappa^2 = (b^2 + \tfrac{3}{4}c^2)/R_g^4\f$
    (0 for a sphere, 1 for a rod),
-   the prolateness \f$S = 27\prod_i(\lambda_i-\bar\lambda)/R_g^6\f$ with
    \f$\bar\lambda = R_g^2/3\f$ (\f$-1/4\f$ for an oblate disk, 0 for a sphere, 2 for a
    prolate rod),
-   the hydrodynamic radius in the Kirkwood approximation
    \f$R_h^{-1} = \langle N^{-2}\sum_{i\ne j} 1/r_{ij}\rangle\f$ and the ratio
    \f$\sqrt{\langle R_g^2\rangle}/R_h\f$ (1.5 for a Gaussian chain, 0.78 for a hard sphere).

Both the per-molecule averages (\f$\langle\kappa^2\rangle\f$, \f$\langle S\rangle\f$) and the
ensemble-ratio forms of Theodorou and Suter (\f$\langle\lambda_1\rangle:\langle\lambda_2\rangle:\langle\lambda_3\rangle\f$,
\f$\langle b\rangle/\langle R_g^2\rangle\f$, \f$\kappa^2 = 1 - 3\langle\lambda_1\lambda_2+\lambda_2\lambda_3+\lambda_3\lambda_1\rangle/\langle R_g^4\rangle\f$)
are reported. The total (ensemble-averaged, lab-frame) tensor
\f$\langle S_{\alpha\beta}\rangle\f$ is written as well, with its eigenvalues, their
fractions of the trace (1/3 each when isotropic), the lab-frame anisotropy
\f$(\lambda_1-\lambda_3)/\mathrm{tr}\f$ and the principal axes; in a framework these show
along which box direction the chains are aligned. When the component has end-to-end
atoms (see `EndToEndAtoms`) the end-to-end distance is sampled as well so the ratio
\f$\langle R^2\rangle/\langle R_g^2\rangle\f$ (6 for an ideal chain) appears in the same summary.
All averages carry 95% confidence intervals from block averaging. Output is written
to the directory `molecule_shape`: a summary `molecule_shape_<component>.s<system>.txt`
and probability-density histograms of \f$R_g\f$, \f$\kappa^2\f$ and \f$S\f$.

For components that declare `RepeatUnits` the descriptors are also sampled per
monomer, from the gyration tensor of the atoms of each repeat unit (same weighting
convention, renormalized within the unit). The table
`monomer_shape_<component>.s<system>.txt` lists, per unit index along the chain and
pooled over all units, \f$\langle R_g\rangle\f$, \f$\langle R_g^2\rangle\f$, the eigenvalues,
\f$\langle b\rangle/\langle R_g^2\rangle\f$, \f$\langle c\rangle/\langle R_g^2\rangle\f$,
\f$\langle\kappa^2\rangle\f$ and \f$\langle S\rangle\f$ with their errors; it shows, for
instance, whether end monomers are more extended than interior ones.

-   `"SampleMoleculeShapeEvery" : integer`\
    Sample the shape descriptors every `int` cycles. Default: `10`.

-   `"WriteMoleculeShapeEvery" : integer`\
    Write the output every `int` cycles. Default: `5000`.

-   `"NumberOfBinsMoleculeShape" : integer`\
    The number of bins of the \f$\kappa^2\f$ and \f$S\f$ histograms (fixed ranges).
    Default: `128`.

-   `"MassWeightedMoleculeShape" : boolean`\
    Weight the atoms by their pseudo-atom masses (\f$w_i = m_i/M\f$) instead of
    uniformly (\f$w_i = 1/N\f$, the polymer-physics convention). Default: `false`.

-   `"RadiusOfGyrationRangeMoleculeShape" : floating-point-number`\
    The upper bound of the radius-of-gyration histogram, in Ångström. By default
    0.6 times the contour length of the bond-graph diameter of the component, which
    bounds every reachable conformation; values beyond the range are dropped from
    the histogram but not from the averages.

-   `"BinWidthRadiusOfGyrationMoleculeShape" : floating-point-number`\
    The bin width of the radius-of-gyration histogram, in Ångström. The range
    grows with the chain length, so this histogram is sized by its bin width: the
    number of bins is the range divided by the width, rounded up. Default: `0.1`.

#### Molecule backbone: chain statistics <a name="molecule-backbone"></a>

`"ComputeMoleculeBackbone" : boolean`

Whether to sample chain statistics along the backbone of this component: the shortest topological path between the end-to-end atoms (see
`EndToEndAtoms`; inferred from `RepeatUnits` or the bond-graph diameter otherwise).
The backbone needs at least three beads. With \f$N_b\f$
backbone beads and bond vectors \f$\mathbf b_i = \mathbf r_{i+1}-\mathbf r_i\f$:

-   mean squared internal distances \f$\langle r^2(k)\rangle\f$ for \f$k = 1\ldots N_b-1\f$
    bonds apart, written with \f$\langle r^2(k)\rangle/(k\langle l\rangle^2)\f$, which is flat
    for an ideal chain;
-   the bond-vector correlation \f$C(k) = \langle\hat{\mathbf b}_i\cdot\hat{\mathbf b}_{i+k}\rangle\f$;
-   the single-chain form factor \f$P(q) = \langle N^{-2}\sum_{ij}\sin(qr_{ij})/(qr_{ij})\rangle\f$
    over all atoms, on a logarithmic \f$q\f$ grid, next to the Debye function of a Gaussian
    chain with the same \f$\langle R_g^2\rangle\f$. It is accumulated as a histogram of the
    intramolecular pair distances (bin width 0.005 Å) and transformed when the output is
    written, so its cost is independent of the number of wave vectors.

The summary `molecule_backbone_<component>.s<system>.txt` reports the contour length
\f$R_{\max}\f$, \f$\langle l\rangle\f$, \f$\langle R^2\rangle\f$, \f$\langle R_g^2\rangle\f$,
the characteristic ratio \f$C_N = \langle R^2\rangle/((N_b-1)\langle l\rangle^2)\f$, the Kuhn
length \f$b_K = \langle R^2\rangle/R_{\max}\f$ and number of Kuhn segments, the persistence
length both from the projection \f$\langle\sum_j\hat{\mathbf b}_{\mathrm{end}}\cdot\mathbf b_j\rangle\f$
and from a fit of \f$\ln C(k) = -k\langle l\rangle/l_p\f$ over the initial decay, and the
Flory exponent from the log-log slope of \f$\langle r^2(k)\rangle\f$ over
\f$N_b/8 \le k \le N_b/2\f$ (zero when the chain is too short to fit). All carry 95%
confidence intervals from block averaging. Output is written to the directory
`molecule_backbone`.

-   `"SampleMoleculeBackboneEvery" : integer`\
    Sample every `int` cycles. The internal distances and the pair-distance histogram
    of the form factor are \f$O(N^2)\f$ per molecule. Default: `10`.

-   `"WriteMoleculeBackboneEvery" : integer`\
    Write the output every `int` cycles. Default: `5000`.

-   `"NumberOfWaveVectorsMoleculeBackbone" : integer`\
    The number of logarithmically spaced wave vectors of the form factor. Default: `64`.

-   `"LowerLimitWaveVectorMoleculeBackbone" : floating-point-number`\
    The smallest wave vector, in 1/Ångström. Default: `0.01`.

-   `"UpperLimitWaveVectorMoleculeBackbone" : floating-point-number`\
    The largest wave vector, in 1/Ångström. Default: `5.0`.

#### End-to-end vector autocorrelation function and relaxation time <a name="end-to-end-autocorrelation-function"></a>

`"ComputeEndToEndACF" : boolean`

Whether to compute the autocorrelation function of the end-to-end vector
\(\mathbf{R}\) of the molecules of this
component, \(C(t) = \langle \mathbf{R}(0) \cdot \mathbf{R}(t) \rangle\). The
component needs end-to-end atoms (`"EndToEndAtoms"` in the molecule definition,
or the two ends of a linear chain). The function is accumulated with the order-N
blocking scheme of the MSD, so a single run covers lags from one sampling
interval up to \(\text{sample interval} \times n^{\text{blocks}}\) with a
logarithmic density of points. Output is written to the directory
`end_to_end_acf`, one file per component with the lag, \(C(t)\), the normalized
function \(C(t)/\langle R^2 \rangle\) and the number of samples per lag. The
header of the file reports \(\langle R^2 \rangle\) and three estimates of the
end-to-end relaxation time \(\tau_R\), the spacing of statistically independent
samples of the end-to-end distance: the integral of \(C(t)/C(0)\) up to its
first zero crossing (flagged as a lower bound when the function has not yet
crossed zero), the time at which \(C(t)/C(0)\) falls to \(1/e\), and the decay
time of a single exponential fitted to \(0.05 < C(t)/C(0) \le 0.5\) (the range
dominated by the slowest Rouse mode). A production run of length \(T\) yields
roughly \(T / (2 \tau_R)\) independent samples of \(R\) per chain.

In molecular dynamics the lag is in picoseconds; in Monte Carlo (no time step)
the lag is in cycles and measures how fast the Monte Carlo moves decorrelate the
chain conformations. Computing the function requires a fixed number of
molecules; do not combine it with insertion/deletion moves.

-   `"SampleEndToEndACFEvery" : integer`\
    Sample the end-to-end vectors every `int` cycles. Default: `10`.

-   `"WriteEndToEndACFEvery" : integer`\
    Write the autocorrelation function every `int` cycles. Default: `5000`.

-   `"NumberOfBlockElementsEndToEndACF" : integer`\
    The number of elements \(n\) per block in the order-N scheme. Default: `25`.

### Component `MC`-moves <a name="component-mc-moves"></a>

-   `"TranslationProbability" : floating-point-number`\
    The relative probability of a translation move. A random displacement is
    drawn along the allowed directions; the internal configuration of the
    molecule is unchanged. The maximum displacement is tuned during the run
    towards a 50% acceptance ratio.

-   `"RandomTranslationProbability" : floating-point-number`\
    The relative probability of a random translation move, in which the
    displacement can reach any position in the box. It is therefore similar to
    reinsertion, except that reinsertion also changes the internal conformation
    and uses biasing.

-   `"TranslationSmartMCProbability" : floating-point-number`\
    The relative probability of a translation smart-MC (force-biased) move of a
    single molecule. Alias: `"ForceBiasTranslationProbability"`.

-   `"RotationProbability" : floating-point-number`\
    The relative probability of a rotation move about the starting bead. A random
    vector on the unit sphere is generated and the molecule is rotated by a
    random angle around it.

-   `"RotationSmartMCProbability" : floating-point-number`\
    The relative probability of a rotation smart-MC (torque-biased) move of a
    single molecule. The trial orientation is updated with a quaternion.

-   `"TranslationRotationSmartMCProbability" : floating-point-number`\
    The relative probability of a combined translation-rotation smart-MC move of
    a single molecule: the displacement is biased along the force and the
    rotation along the torque in a single trial move, using one shared gradient
    evaluation. The two step sizes (Angstrom and radians) are optimized
    independently.

-   `"ReinsertionProbability" : floating-point-number`\
    The relative probability of a full `CBMC` reinsertion move. Several first
    beads are trial-placed and one is chosen from its Boltzmann weight; the rest
    of the molecule is then grown with biasing. This move is very useful, and
    often necessary, to change the internal configuration of flexible molecules.

-   `"PartialReinsertionProbability" : floating-point-number`\
    The relative probability of a partial `CBMC` reinsertion move, which regrows
    only part of the molecule. The parts are declared in the molecule file as
    `"Partial-reinsertion" : [[fixed atoms], [fixed atoms], ...]`, a list of
    sets of atom indices; each move picks one set at random, keeps those atoms
    where they are, and regrows all other atoms with `CBMC`.
    A fixed set may be disconnected in the bond graph. Fixed atoms on both sides
    of the regrown part make the move a fixed-endpoint (bridging) regrowth: the
    interior segment is grown from one fixed side and its last bead is closed
    onto the other, e.g. `[0, 1, 4, 5]` for a six-bead chain regrows beads 2
    and 3 between them, `[0, 5]` regrows the whole interior. The closure is
    exact (the closing bead is drawn from its two bond-length distributions in
    bipolar coordinates, `FIXED` bonds are satisfied exactly, and the beads
    before it are steered towards closable geometries by a bias that cancels
    from the acceptance rule), so the acceptance rule is the usual
    `W_new / W_old`. Limits: the regrown part between fixed atoms must be
    flexible single-atom beads (a rigid body or ring bonded to fixed atoms on
    more than one side is rejected when the molecule is read), a bead may not be
    bonded to more than two fixed atoms, and declared chiral centres may not
    include the closing bead. Long or stiff bridged segments close with low
    acceptance; keep bridged segments short.

-   `"PivotProbability" : floating-point-number`\
    The relative probability of a pivot move for flexible molecules. A bond
    that is not part of a ring or interior to a rigid fragment is chosen at
    random and the smaller of the two chain parts hanging off it is rotated
    rigidly about the bond axis by a random angle. All bond lengths and bend
    angles are preserved; only the torsions through the pivot bond and the
    non-bonded energies change. A fraction `"PivotRandomizationFraction"`
    (default 0.2) of the attempts draws the angle uniformly from
    $[-\pi, \pi]$; the remainder uses an adaptive window tuned towards a 50%
    acceptance ratio. Not available with polarization.

-   `"CrankshaftProbability" : floating-point-number`\
    The relative probability of a crankshaft move for flexible molecules. A
    small connected segment (at most `"CrankshaftMaxSegmentSize"` atoms,
    default 4) attached to exactly two anchor atoms is rotated rigidly about
    the axis through the anchors. All bond lengths are preserved; the bend
    angles and torsions at the two junctions change. The move relaxes the
    chain interior locally without moving the chain ends and can also rotate
    segments of flexible rings. A fraction `"CrankshaftRandomizationFraction"`
    (default 0.2) of the attempts draws the angle uniformly from
    $[-\pi, \pi]$. Not available with polarization.

-   `"ConcertedRotationProbability" : floating-point-number`\
    The relative probability of a concerted-rotation (ConRot) move (Dodd,
    Boone and Theodorou, Mol. Phys. 78, 961 (1993)) for flexible chain
    molecules. Alias: `"ConRotProbability"`. A window of eight consecutive
    backbone atoms $a_0 \ldots a_7$ is chosen at random; atom $a_2$ is
    rotated about the $a_0 a_1$ bond by a random driver angle and the
    trimer $a_3, a_4, a_5$ is rebridged onto the fixed atoms $a_6, a_7$
    such that all bond lengths and bend angles of the window are conserved
    exactly (the trimer positions are the discrete solutions of a loop-closure
    problem, found by a one-dimensional scan). Only the seven torsions of the
    window, the bends and torsions to side groups at $a_1$ and $a_6$, and
    the non-bonded energies change, and everything outside the window stays in
    place. The move is therefore strictly local, which makes it the move of
    choice for long chains and for dense melts, where the pivot move fails.
    The acceptance rule contains, besides the Metropolis factor, the ratio of
    the numbers of closure solutions and the ratio of the closure Jacobians of
    the old and new window; it is exact for `FIXED` and for harmonic bonds
    and bends. Side groups hanging off $a_2 \ldots a_5$ are carried rigidly
    with their backbone atom. A window is only valid if no ring passes through
    it, rigid fragments do not straddle its moving groups, and there is no
    `FIXED` or `RIGID` bend at $a_1$ or $a_6$ that would be violated. The
    first and last two backbone atoms of a linear chain are never moved, so
    the move must be combined with a move that samples the chain ends (pivot,
    reptation or (partial) reinsertion). A fraction
    `"ConcertedRotationRandomizationFraction"` (default 0.2) of the attempts
    draws the driver angle uniformly from $[-\pi, \pi]$; the remainder uses
    an adaptive window. Not available with polarization.

-   `"BeadDisplacementProbability" : floating-point-number`\
    The relative probability of a single-bead displacement for flexible
    molecules. One bead is chosen at random and displaced along one random
    Cartesian direction by a random amount within an adaptive maximum
    displacement (one per direction, tuned towards a 50% acceptance ratio,
    initial value 0.3 Å); all other atoms stay in place. Every bond, bend and
    torsion the bead takes part in changes, so the move relaxes flexible bond
    lengths and bend angles locally, which none of the rotation moves (pivot,
    crankshaft, concerted rotation, bead flip) can do. Beads inside a rigid
    fragment of more than one atom and beads that take part in a `FIXED` bond
    or a `FIXED`/`RIGID` bend or torsion are never chosen; a molecule without a
    displaceable bead rejects the move. Not available with polarization.

-   `"BeadFlipProbability" : floating-point-number`\
    The relative probability of a bead flip for flexible molecules: a
    bond-length-preserving rotation of a single bead. An interior bead with
    two neighbours is rotated about the axis through its neighbours (kink
    jump), which keeps both bonds and the bend centred on the bead and changes
    the bends and torsions at the neighbours. A terminal bead is rotated about
    a random axis through its neighbour (end rotation), so its bond direction
    moves over the sphere and the bend and torsions at the neighbour change.
    Beads with three or more neighbours, beads inside a rigid fragment, and
    beads whose rotation would violate a `FIXED`/`RIGID` bend or torsion are
    never chosen; `FIXED` bonds are allowed. The move is the smallest local
    move that respects fixed bond lengths, moves the chain ends (which the
    pivot cannot) and is cheap; the interior flip coincides with a crankshaft
    of a one-atom segment. A fraction `"BeadFlipRandomizationFraction"`
    (default 0.2) of the attempts randomizes the rotation completely (angle
    uniform on $[-\pi, \pi]$, or a uniformly random direction on the sphere
    for a terminal bead); the remainder uses an adaptive angle window. Not
    available with polarization.

-   `"DoubleBridgingProbability" : floating-point-number`\
    The relative probability of a double-bridging (DB) move (Karayiannis,
    Mavrantzas and Theodorou, Phys. Rev. Lett. 88, 105503 (2002); J. Chem.
    Phys. 117, 5465 (2002)) for acyclic chain molecules of the same
    component: a connectivity-altering move that exchanges the tails of two
    chains. A backbone site $s$ is chosen at random; on both chains the trimer
    of backbone units $s+1, s+2, s+3$ is excised and the head of each chain
    (units up to $s$) is bridged onto the tail of the other (units from $s+4$)
    with a new trimer that keeps all bond lengths and bend angles (the same
    loop-closure problem as the concerted rotation). Using the same site on
    both chains keeps every chain length unchanged, so the move samples the
    same monodisperse ensemble as the other moves; the large-scale
    conformations of two chains change at once while every atom outside the
    two trimers stays in place, which decorrelates the end-to-end vectors
    orders of magnitude faster than local moves. The partner chain is chosen
    uniformly among the chains whose anchor atoms are within bridging reach
    (the sum of the four backbone bond lengths); the acceptance rule contains
    the Metropolis factor, the solution-count and closure-Jacobian ratios of
    both bridges, and the ratio of the partner-selection probabilities. Side
    groups travel with their backbone atom (those of the re-bridged trimers
    rigidly with its local frame). The backbone is the path between the
    `"EndToEndAtoms"` (by default the graph diameter); a site is valid when no
    rigid fragment straddles the moving groups and no `FIXED`/`RIGID` bend at
    the anchors would be violated. Requires at least two whole (non-fractional)
    molecules of the component and a backbone of at least seven units. Not
    available with polarization.

-   `"IntramolecularDoubleRebridgingProbability" : floating-point-number`\
    The relative probability of an intramolecular double rebridging (IDR)
    move (Karayiannis et al., J. Chem. Phys. 117, 5465 (2002)), the
    single-chain analogue of double bridging. Alias:
    `"DoubleRebridgingProbability"`. Two sites $a < b$ ($b \ge a + 5$) are
    chosen on the same chain; the trimers at both sites are excised and the
    segment of units $a+4 \ldots b$ between them is reversed by bridging the
    chain head onto its far end and its near end onto the chain tail, again
    with bond lengths and bend angles preserved. Chain length and
    connectivity are unchanged; the whole segment changes conformation while
    every atom outside the two trimers stays in place. The segment must be
    congruent under reversal (unit $k$ and unit $a+4+b-k$ have the same
    pseudo-atom types, charges and side-group structure; a homopolymer
    backbone qualifies). The acceptance rule contains the Metropolis factor,
    the closure factors of both bridges and the ratio of the chart volume
    elements of the reversed segment (unity for uniform fixed bonds and
    bends). Requires a backbone of at least twelve units. Not available with
    polarization.

-   `"ReptationProbability" : floating-point-number`\
    The relative probability of a reptation (slithering-snake) move for chain
    molecules that declare their repeat units in the molecule file as
    `"RepeatUnits" : [[atoms of unit 1], [atoms of unit 2], ...]`, an ordered
    partition of the atoms into monomer blocks that is checked at read time to
    be shift-periodic (identical types, charges, connectivity, potential terms
    and rigid fragments under the one-unit shift). The move removes the repeat
    unit at one chain end (chosen with 50% probability) and grows a new unit
    at the opposite end with `CBMC`, translating the chain by one monomer along
    its own contour; the acceptance rule is the usual ratio of the grow and
    retrace Rosenbluth weights.

-   `"SwapConventionalProbability" : floating-point-number`\
    The relative probability of a conventional (non-CBMC) insertion or deletion
    move, each chosen with 50% probability. The swap move imposes chemical
    equilibrium between the system and an imaginary particle reservoir.

-   `"SwapProbability" : floating-point-number`\
    The relative probability of a `CBMC` insertion or deletion move (insertion or
    deletion chosen with 50% probability each). Like the conventional swap it
    imposes chemical equilibrium with a reservoir, but it grows the molecule from
    multiple first beads using biasing.

-   `"CFCMC_SwapProbability" : floating-point-number`\
    The relative probability of an insertion or deletion move via the `CFCMC`
    scheme.

-   `"CFCMC_CBMC_SwapProbability" : floating-point-number`\
    The relative probability of an insertion or deletion move via the combined
    `CB/CFCMC` scheme.

-   `"GibbsSwapCBMCProbability" : floating-point-number`\
    The relative probability of a Gibbs swap move, which transfers a randomly
    selected molecule from one box to the other (50% from box `I` to `II`, 50%
    the other way).

-   `"GibbsSwapCFCMCProbability" : floating-point-number`\
    The relative probability of a Gibbs swap move using the `CFCMC` scheme.

-   `"WidomProbability" : floating-point-number`\
    The relative probability of a Widom particle-insertion move, which measures
    the chemical potential and relates directly to the Henry coefficient and the
    heat of adsorption.
