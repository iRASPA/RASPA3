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
  * [Molecular dynamics parameters](#molecular-dynamics-parameters)
  * [Options to measure properties](#options-to-measure-properties)
    * [Output pdb-movies](#output-pdb-movies)
    * [Histogram of the energy](#histogram-of-the-energy)
    * [Histogram of the number of molecules](#histogram-of-the-number-of-molecules)
    * [Histograms of the intra-molecular geometry](#histograms-of-the-intra-molecular-geometry)
    * [Radial Distribution Function (RDF) force-based](#radial-distribution-function-rdf-force-based)
    * [Radial Distribution Function (RDF) conventional](#radial-distribution-function-rdf-conventional)
    * [Mean-Squared Displacement (MSD) order-N](#mean-squared-displacement-msd-order-n)
    * [Density grids](#density-grids)
* [Force field options](#force-field-options)
* [Component options](#component-options)
  * [Component properties](#component-properties)
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
    there). The optional property files (RDFs, density grids, energy and
    molecule-count histograms, molecule properties, movies, and the
    number-of-molecules/volume evolution files) are written per replica, keyed
    by the replica index (`.s{k}`). Restart files (JSON and binary) are not
    supported by this driver.

    The driver spawns one worker thread per temperature (plain C++ threads,
    no OpenMP); leave `"NumberOfThreads"` at its default of `1` so the
    per-energy-evaluation thread pool stays serial and the machine is not
    oversubscribed. Note that the swap move requires rigid, whole-molecule
    replicas: systems with fractional (CFCMC) molecules, flexible components,
    or reactions reject all swap attempts.

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

    The driver spawns one worker thread per grid point (plain C++ threads, no
    OpenMP) — with N_T temperatures and N_P pressures that is N_T × N_P
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
    `"ExternalTemperature"` gives a one-temperature run). Every walker runs
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

    The driver spawns one worker thread per walker (plain C++ threads, no
    OpenMP) — with N_T temperatures and N_W windows that is N_T × N_W
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

    The driver spawns `NumberOfLambdaBins` worker threads (plain C++ threads,
    no OpenMP); leave `"NumberOfThreads"` at its default of `1` so the
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
    the stress are also reported. Default: `5000`.

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
    For `"SimulationType" : "ParallelTempering"`, `"HyperParallelTempering"`
    and `"ReweightedHistogram"`: how often (in cycles) a sweep of
    configuration swaps between replicas at neighboring state points is
    attempted (`0` disables the swaps). Default: `10`.

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
    For `"SimulationType" : "ReweightedHistogram"`: the pressure range (in
    Pascal) of the reweighted isotherms. Default: the span of the pressure
    ladder `"ExternalPressures"`.

-   `"ReweightingNumberOfPressures" : integer`\
    For `"SimulationType" : "ReweightedHistogram"` and `"ParallelTMMC"`: the
    number of log-spaced pressures the reweighted isotherms are evaluated at.
    Default: `100`.

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
    Selects the threading backend explicitly. One of `"Serial"`, `"ThreadPool"`,
    `"OpenMP"`, or `"GPU-Offload"`.

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
    The temperature ladder for `"SimulationType" : "ParallelTempering"` (a
    sorted list of at least two temperatures in Kelvin),
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
    for those simulation types.

-   `"ExternalPressureX" / "ExternalPressureY" / "ExternalPressureZ" : floating-point-number`\
    Override individual diagonal components of the pressure tensor, for
    anisotropic (directional) pressure control. Each defaults to
    `"ExternalPressure"` when not given.

-   `"ChemicalPotential" : floating-point-number`\
    Sets the imposed chemical potential (in internal units); the corresponding
    fugacity is derived from it and the temperature.

-   `"MacroStateMinimumNumberOfMolecules" / "MacroStateMaximumNumberOfMolecules" : integer`\
    For `"SimulationType" : "MonteCarloTransitionMatrix"` and
    `"ParallelTMMC"`: the macrostate range (the total molecule count) the
    transition-matrix walk is confined to. For `"ParallelTMMC"` the range is
    split into `"NumberOfWindows"` windows, and a minimum of `0` enables the
    exact normalization of the saturation pressure by the empty-box state.
    Defaults: `0` and `100`.

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
    The pressure-coupling time scale in picoseconds. Default: `1.0`.

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
    For `"Type" : "Framework"`, loads the framework from the file `string.cif`.

-   `"NumberOfUnitCells" : [integer, integer, integer]`\
    The number of unit cells in the `x`, `y`, and `z` directions. The super-cell
    contains these unit cells, and periodic boundary conditions are applied at
    the super-cell level (*not* at the unit-cell level). Default: `[1, 1, 1]`.

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
    `"HyperParallelTempering"` and `"ReweightedHistogram"`, where the swaps
    are performed by the driver at the barrier synchronization points (see
    `"ParallelTemperingSwapEvery"`).

-   `"TranslationSmartMCAllProbability" : floating-point-number`\
    The probability per cycle of attempting a translation smart-MC move that
    displaces all molecules simultaneously along the forces acting on them.
    Alias: `"ForceBiasTranslationAllProbability"`.

-   `"RotationSmartMCAllProbability" : floating-point-number`\
    The probability per cycle of attempting a rotation smart-MC move that
    rotates all rigid multi-atomic molecules simultaneously along the torques
    acting on them (quaternion update).

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

#### Histograms of the intra-molecular geometry <a name="histograms-of-the-intra-molecular-geometry"></a>

`"ComputeMoleculeProperties" : boolean`

Whether to accumulate probability histograms of the intra-molecular geometry
(bond lengths, bend angles, and torsion/dihedral angles) for every flexible
component. This mirrors the "molecule properties" analysis from RASPA2. Output
is written to the directory `molecule_properties`.

-   `"SampleMoleculePropertiesEvery" : integer`\
    Sample the histograms every `int` cycles. Default: `10`.

-   `"WriteMoleculePropertiesEvery" : integer`\
    Write the histograms every `int` cycles. Default: `5000`.

-   `"NumberOfBinsMoleculeProperties" : integer`\
    The number of bins in each histogram. Default: `128`.

-   `"BondRangeMoleculeProperties" : floating-point-number`\
    The upper bound of the bond-length histogram, in Ångström. The bend range is
    fixed to `[0, 180]` degrees and the torsion range to `[-180, 180]` degrees.
    Default: `4.0`.

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

----------------------------------------------------------------------------------

## Component options <a name="component-options"></a>

### Component properties <a name="component-properties"></a>

-   `"Name" : string`\
    The name of the component. For a rigid or flexible molecule this is also the
    base name of its definition file.

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

-   `"CreateNumberOfMolecules" : integer`\
    The number of molecules to create for this component at start-up. These
    molecules are created *in addition* to anything read from a restart file, so
    when restarting this value is usually set back to zero. Setting it
    unreasonably high can cause an infinite loop: the routine only accepts
    molecules whose growth causes no overlap (energy below the overlap
    criterion). The starting configurations are far from optimal, so substantial
    equilibration is needed to relax the energy; the `CBMC` growth can, however,
    reach very high densities.

-   `"StartingBead" : integer`\
    The index of the bead from which `CBMC` growth starts. Must be smaller than
    the number of atoms in the molecule.

-   `"BlockingPockets" : [[3 x floating-point-number, floating-point-number]]`\
    Blocks certain pockets of the simulation volume so molecules cannot grow
    into them. A typical example is the sodalite cages in FAU- and LTA-type
    zeolites, which are inaccessible to methane and larger molecules. Each pocket
    is a list of four numbers: the fractional positions $s_x$, $s_y$, $s_z$ and a
    radius in Ångström. For example, the blocking pockets of ITQ-29 for small
    molecules are:

        "BlockingPockets" : [
                   [0.0,       0.0,        0.0,       4.0],
                   [0.5,       0.0,        0.0,       0.5],
                   [0.0,       0.5,        0.0,       0.5],
                   [0.0,       0.0,        0.5,       0.5]
                 ]

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

-   `"LnPartitionFunction" : number or string`\
    The natural logarithm of the (reduced) partition function used for reactions.
    Give a number to set it directly, or a species name (or `"auto"`, which uses
    the component name) to look it up in the embedded thermochemical database.
    The lookup is evaluated at each system's `"ExternalTemperature"` and uses the
    database selected by `"ThermochemicalDatabase"`.

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
    only part of the molecule.

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
