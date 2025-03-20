# Simulation input options
\page commands Commands

## Table of Contents

<!-- TOC -->
* [Input sections](#input-sections)
* [General options](#general-options)
  * [Simulation types](#simulation-types)
  * [Simulation duration](#simulation-duration)
  * [Restart and crash-recovery](#restart-and-crash-recovery)
  * [Printing options](#printing-options)
* [System options](#system-options)
  * [Operating conditions and thermostat/barostat-parameters](#operating-conditions-and-thermostatbarostat-parameters)
  * [Box/Framework options](#boxframework-options)
  * [Force field definitions](#force-field-definitions)
  * [System `MC`-moves](#system-mc-moves)
  * [Molecular dynamics parameters](#molecular-dynamics-parameters)
  * [Options to measure properties](#options-to-measure-properties)
    * [Output pdb-movies](#output-pdb-movies)
    * [Histogram of the energy](#histogram-of-the-energy)
    * [Histogram of the number of molecules](#histogram-of-the-number-of-molecules)
    * [Radial Distribution Function (RDF) force-based](#radial-distribution-function-rdf-force-based)
    * [Radial Distribution Function (RDF) conventional](#radial-distribution-function-rdf-conventional)
    * [Mean-Squared Displacement (MSD) order-N](#mean-squared-displacement-msd-order-n)
    * [Density grids](#density-grids)
* [Component options](#component-options)
  * [Component properties](#component-properties)
  * [Component **`MC`**-moves](#component-mc-moves)
<!-- TOC -->

----------------------------------------------------------------------------------

## Input sections

```json
{
  "SimulationType" : "MolecularDynamics",
  "NumberOfCycles" : 100000,
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

----------------------------------------------------------------------------------

## General options

### Simulation types

-   `"SimulationType" : "MonteCarlo"`
    Starts the Monte Carlo part of `RASPA`. The particular ensemble is
    not specified but implicitly deduced from the specified Monte Carlo
    moves. Note that a `MD`-move can be used for hybrid `MC`/`MD`.

-   `"SimulationType" : "MolecularDynamics"`
    Starts the Molecular Dynamics part of `RASPA`. The ensemble must be
    explicitly specified.

### Simulation duration

-   `"NumberOfCycles" : integer`
    The number of cycles for the production run. For Monte Carlo a cycle
    consists of $N$ steps, where $N$ is the amount of molecules with a
    minimum of 20 steps. This means that on average during each cycle on
    each molecule a Monte Carlo move has been attempted (either
    successful or unsuccessful). For MD the number of cycles is simply
    the amount of integration steps.

-   `"NumberOfInitializationCycles": integer`
    The number of cycles used to initialize the system using Monte
    Carlo. This can be used for both Monte Carlo as well as Molecular
    Dynamics to quickly minimize the positions of the atoms in the
    system.

-   `"NumberOfEquilibrationCycles" : integer`
    For Molecular Dynamics it is the number of MD steps to equilibrate
    the velocities in the systems. After this equilibration the
    production run is started. For Monte Carlo, in particular `CFCMC`,
    the equilibration-phase is used to measure the biasing factors, using
    Wang-Landau estimation.

### Restart and crash-recovery

-   `"RestartFile" : boolean`
    Reads the positions, velocities, and force from the directory
    'RestartInitial'. Any creation of molecules in the
    'simulation.input' file will be in addition and after this first
    read from file. This is useful to load initial positions of cations
    for example, and after that create adsorbates. The restart file is
    written at 'PrintEvery' intervals.

-   `"ContinueAfterCrash" : boolean`
    Write a binary file containing the complete status of the program.
    The file name is 'restart_data.bin'. With this option to `true` the 
    presence of this file will result in continuation from the point where
    the program was at the moment of outputting this file. The file can be
    quite big (several hundreds of megabytes) and will be outputted
    every 'WriteBinaryRestartFileEvery' cycles.

-   `"WriteBinaryRestartFileEvery" : integer`
    The output frequency (i.e. every `int` cycles) of writing the
    crash-recovery file.

### Printing options

-   `"PrintEvery" : integer`
    Prints the loadings (when a framework is present) and energies every
    `int` cycles. For MD information like energy conservation and
    stress are printed.

### Parameter tuning

-   `"RescaleWangLandauEvery" : integer`
    Determines the frequency of updates for optimizing the λ parameter in
    for example Continuous Fractional Component Monte Carlo during the
    equilibration phase of the simulation.

-   `"OptimizeMCMovesEvery" : integer`
    Determines the frequency of updating the maximum change in the Monte
    Carlo moves towards an optimal acceptance ratio (default: 0.5). For
    Translation, the maximum displacement is optimized, for rotation the
    maximum angle, for hybrid MC the maximum timestep etc.

### Systems & Components

-   `"Systems" : list`
    List of system settings, each containing a dictionary with possible
    key-value pairs described below. Multiple systems can be owned by one
    process and Monte Carlo moves starting with Gibbs and parallel tempering
    act on two systems.

-   `"Components" : list`
    List of component settings, each containing a dictionary with possible
    key-value pairs described below.

----------------------------------------------------------------------------------

## System options

### Operating conditions and thermostat/barostat-parameters

-   `"ExternalTemperature" : floating-point-number`
    The external temperature in Kelvin for the system. From this, the system
    beta is calculated, which is central in all statistics. Default: `298`

-   `"ExternalPressure" : floating-point-number`
    The external pressure in Pascal for the system. Default: `0`

-   `"ThermostatChainLength" : integer`
    The length of the chain to thermostat the system. Default: `5`

-   `"NumberOfYoshidaSuzukiSteps" : integer`
    The number of Yoshida/Suzuki multiple timesteps. Default: `5`

-   `"TimeScaleParameterThermostat" : floating-point-number`
    The time scale on which the system thermostat evolves. Default:
    `0.15`

### Box/Framework options

-   `"Type" : string`
    Sets the system type. The type string can be:

    -   `"Box"`
        Sets the system to a simulation cell where the lengths and
        angles of the cell can be explicitely be specified.

    -   `"Framework"`
        Set the system to type 'Framework'. The cell lengths and cell
        angles follows from the specified framework file.

-   `"BoxLengths" : [floating-point-number, floating-point-number, floating-point-number]`
    The cell dimensions of rectangular box of system in Angstroms.

-   `"BoxAngles" : [floating-point-number, floating-point-number, floating-point-number]`
    The cell angles of rectangular box of system in Degrees.

-   `"Name" : string`
    For `"Type" : "Framework"`, loads the framework with filename
    `string.cif`.

-   `"NumberOfUnitCells" : [integer, integer, integer]`
    The number of unit cells in `x`, `y`, and `z` direction for the
    system. The super-cell will contain the unit cells, and periodic
    boundary conditions will be applied on the super-cell level (*not*
    on a unit cell level).

-   `"HeliumVoidFraction" : floating-point-number`
    Sets the void fraction as measured by probing the structure with
    helium a room temperature. This quantity has to be obtained from a
    separate simulation and is essential to compute the
    *excess*-adsorption during the simulation.

-   `"UseChargesFrom" : string`
    Specifies the way to define the framework charges. The string can
    be:

    -   `"PseudoAtoms"`
        Takes the charges from the force-field definition file.

    -   `"CIF_File"`
        Takes the charges from the definitions in the `CIF`-file. The
        charges for the atoms list in the file needs to be given using
        the `_atom_site_charge` tag. Using this option allows for
        individual charges for each framework atom, even when having the
        same atom-type.

    -   `"ChargeEquilibration"`
        Computes the framework charges uses the charge equilibration
        scheme of Wilmer and Snurr. The charges are symmetrized over the
        asymmetric atoms.

### Force field definitions

-   `"ForceField" : string`
    Reads in the force field file `string.json`, Note that if this file
    is in the working directory then this will be read and used instead
    of:

        ${RASPA_DIR}/simulations/share/raspa3/forcefield/string/force_field.json

-   `"CutOffFrameworkVDW" : floating-point-number`
    The cutoff of the Van der Waals potentials for framework-molecule
    interactions. Interactions longer then this distance are omitted
    from the energy and force computations.

-   `"CutOffMoleculeVDW" : floating-point-number`
    The cutoff of the Van der Waals potentials for molecule-molecule
    interactions. Interactions longer then this distance are omitted
    from the energy and force computations.

-   `"CutOffVDW" : floating-point-number`
    The cutoff of the Van der Waals potentials. Interactions longer then
    this distance are omitted from the energy and force computations.
    The option `CutOffVDW` implies setting both `CutOffFrameworkVDW` and
    `CutOffMoleculeVDW`.

-   `"CutOffCoulombic" : floating-point-number`
    The cutoff of the charge-charge potential. The potential is
    truncated at the cutoff. No tail-corrections are (or can be)
    applied. The only way to include the long-range part is to use
    'ChargeMethod Ewald'. The parameter is also used in combination with
    the Ewald precision to compute the number of wave vectors and Ewald
    parameter $\alpha$. For the Ewald summation using rather large unit
    cells, a charge-charge cutoff of about half the smallest box-length
    would be advisable in order to avoid the use of an excessive amount
    of wave-vectors in Fourier space. For non-Ewald methods the cutoff
    should be as large as possible (greater than about 30 Å).

-   `"CutOff" : floating-point-number`
    Implies setting `CutOffFrameworkVDW`, `CutOffMoleculeVDW`, and
    `CutOffCoulomb`

-   `"ChargeMethod" : string` Sets the method to compute charges. The
    string can be:

    -   `"None"`
        Skips the entire charge calculation and should only be used when
        all adsorbates do not contain any charges.

    -   `"Ewald"`
        Switches on the Ewald summation for the charge calculation.

### System `MC`-moves

-   `"VolumeChangeProbability" : floating-point-number`
    The probability per cycle to attempt a volume-change. Rigid
    molecules are scaled by center-of-mass, while flexible molecules and
    the framework is atomically scaled.

-   `"GibbVolumeChangeProbability" : floating-point-number`
    The probability per cycle to attempt a Gibbs volume-change `MC` move
    during a Gibbs ensemble simulation. The total volume of the two
    boxes (usually one for the gas phase, one for the liquid phase)
    remains constant, but the individual volume of the boxes are
    changed. The volumes are changed by a random change in
    $\ln(V_I/V_{II})$.

-   `"GibbVolumeChangeProbability" : floating-point-number`
    The probability per cycle to attempt a Hybrid MC move. A hybrid MC
    move propagates the Hamiltonian via a short Molecular Dynamics
    simulation and accepts the new state based on the drift.

### Molecular dynamics parameters

-   `"TimeStep" : floating-point-number`
    The time step in picoseconds for `MD` integration. Default value:
    `0.0005`

-   `"Ensemble" : string`
    Sets the ensemble. The ensemble string can be:

    -   `"NVE"`\
        The micro canonical ensemble, the number of particle $N$, the
        volume $V$, and the energy $E$ are constant.

    -   `"NVT"`\
        The canonical ensemble, the number of particle $N$, the volume
        $V$, and the average temperature $\left\langle T\right\rangle$
        are constant. Instantaneous values for the temperature are
        fluctuating.

### Options to measure properties

#### Output pdb-movies

`"OutputPDBMovie" : boolean`

Sets whether or not to output simulation snapshots to pdb movies.
Output is written to the directory `movies`.

-   `"SampleMovieEvery" : integer`
    Sample the movie every `int` cycles. Default: `1`

#### Histogram of the energy


`"ComputeEnergyHistogram" : boolean`


Sets whether or not to compute a histogram of the energy for the current
system. For example, during adsorption it keeps track of the total
energy, the VDW energy, the Coulombic energy, and the polarization
energy.\
Output is written to the directory `energy_histogram`.

-   `"SampleEnergyHistogramEvery" : integer`
    Sample the energy histogram of the system every `int` cycles.
    Default: `1`

-   `"WriteEnergyHistogramEvery" : integer`
    Writes the energy histogram of the system every `int` cycles.
    Default: `5000`

-   `"NumberOfBinsEnergyHistogram" : integer`
    Sets the number of elements of the histogram. Default: `128`

-   `"LowerLimitEnergyHistogram" : floating-point-number`
    The lower limit of the histogram. Default: `-5000`

-   `"UpperLimitEnergyHistogram" : floating-point-number`
    The upper limit of the histogram. Default: `1000`

#### Histogram of the number of molecules


`"ComputeNumberOfMoleculesHistogram" : boolean`


Sets whether or not to compute the histograms of the number of molecules
for the current system. In open ensembles the number of molecules
fluctuates.\
Output is written to the directory `number_of_molecules_histogram`.

-   `"SampleNumberOfMoleculesHistogramEvery" : integer`
    Sample the histogram every `int` cycles. Default: `1`

-   `"WriteNumberOfMoleculesHistogramEvery" : integer`
    Output the histogram every `int` cycles. Default: `5000`

-   `"LowerLimitNumberOfMoleculesHistogram" : floating-point-number`
    The lower limit of the histograms. Default: `0`

-   `"UpperLimitNumberOfMoleculesHistogram" : floating-point-number`
    The upper limit of the histograms. Default: `200`

#### Radial Distribution Function (RDF) force-based


`"ComputeRDF" : boolean`


Sets whether or not to compute the radial distribution function (RDF).\
Output is written to the directory `rdf`.

-   `"SampleRDFEvery" : integer`
    Sample the rdf every `int` cycles. Default: `10`

-   `"WriteRDFEvery" : integer`
    Output the rdf every `int` cycles. Default: `5000`

-   `"NumberOfBinsRDF" : integer`
    Sets the number of elements of the rdf. Default: `128`

-   `"UpperLimitRDF" : floating-point-number`
    The upper limit of the rdf. Default: `15.0`

#### Radial Distribution Function (RDF) conventional


`"ComputeConventionalRDF" : boolean`


Sets whether or not to compute the radial distribution function (RDF).
Output is written to the directory `conventional_rdf`.

-   `"SampleConventionalRDFEvery" : integer`
    Sample the rdf every `int` cycles. Default: `10`

-   `"WriteConventionalRDFEvery" : integer`
    Output the rdf every `int` cycles. Default: `5000`

-   `"NumberOfBinsConventionalRDF" : integer`
    Sets the number of elements of the rdf. Default: `128`

-   `"UpperLimitConventionalRDF" : floating-point-number`
    The upper limit of the rdf. Default: `15.0`

#### Mean-Squared Displacement (MSD) order-N


`"ComputeMSD" : boolean`


Sets whether or not to compute the mean-squared displacement (MSD).
Output is written to the directory `msd`.

-   `"SampleMSDEvery" : integer`
    Sample the msd every `int` cycles. Default: `10`

-   `"WriteMSDEvery" : integer`
    Output the msd every `int` cycles. Default: `5000`

-   `"NumberOfBlockElementsMSD" : integer`
    The number of elements per block of the msd. Default: `25.0`

#### Density grids


`"ComputeDensityGrid" : boolean`


Sets whether or not to compute the density grids.
Output is written to the directory `density_grids`.

-   `"SampleDensityGridEvery" : integer`
    Sample the density grids every `int` cycles. Default: `10`

-   `"WriteDensityGridEvery" : integer`
    Output the density grids every `int` cycles. Default: `5000`

-   `"DensityGridSize" : [integer, integer, integer]`
    Sets the size of the density grids. Default: `[128, 128, 128]`

----------------------------------------------------------------------------------

## Component options

### Component properties

-   `"Name" : string`
    The descriptive name of the component.

-   `"MolFraction" : floating-point-number`
    The mol fraction of this component in the mixture. The values can be
    specified relative to other components, as the fractions are
    normalized afterwards. The partial pressures for each component are
    computed from the total pressure and the mol fraction per component.

-   `"FugacityCoefficient" : floating-point-number`
    The fugacity coefficient for the current component. For values 0 (or
    by not specifying this line), the fugacity coefficients are
    automatically computed using the Peng-Robinson equation of state.
    Note the critical pressure, critical temperature, and acentric
    factor need to be specified in the molecule file.

-   `"IdealGasRosenbluthWeight" : floating-point-number`
    The ideal Rosenbluth weight is the growth factor of the `CBMC`
    algorithm for a single chain in an empty box. The value only depends
    on temperature and therefore needs to be computed only once. For
    adsorption, specifying the value in advance is convenient because
    the applied pressure does not need to be corrected afterwards (the
    Rosenbluth weight corresponds to a shift in the chemical potential
    reference value, and the chemical potential is directly obtained
    from the fugacity). For equimolar mixtures this is essential.

-   `"CreateNumberOfMolecules" : integer`
    The number of molecule to create for the current component. Note
    these molecules are *in addition* to anything read in by using a
    restart-file. Usually, when the restart-file is used the amount here
    should be put back to zero. A warning, putting this value
    unreasonably high results in an infinite loop. The routine accepts
    molecules that are grown causing no overlap (energy smaller than
    'EnergyOverlapCriteria'). Also the initial starting configurations
    are far from optimal and substantial equilibration is needed to
    reduce the energy. However, the `CBMC` growth is able to reach very
    high densities.

-   `"BlockingPockets" : [[3 x floating-point-number, floating-point-number]]`
    Block certain pockets in the simulation volume. The growth of a
    molecule is not allowed in a blocked pocket. A typical example is
    the sodalite cages in FAU and LTA-type zeolites, these are not
    accessible to molecules like methane and bigger. The pockets are
    specified as a list of 4 floating points numbers: the $s_x$, $s_y$,
    $s_z$ fractional positions and a radius in Angstrom.

    For example, blocking pockets for ITQ-29 for small molecules are
    specified as:

        "BlockingPockets" : [
                   [0.0,       0.0,        0.0,       4.0],
                   [0.5,       0.0,        0.0,       0.5],
                   [0.0,       0.5,        0.0,       0.5],
                   [0.0,       0.0,        0.5,       0.5]
                 ]

-   `"LambdaBiasFileName" : string`
    Pointing this parameter to a json file containing preset λ values
    allows optimized CFCMC simulations to be run without Wang-Landau
    estimation of the biasing weights.

-   `"ThermodynamicIntegration" : boolean`
    Boolean switch to determine whether to do a thermodynamic integration
    over dU/dλ for the fractional component.

### Component **`MC`**-moves

-   `"TranslationProbability" : floating-point-number`
    The relative probability to attempt a translation move for the
    current component. A random displacement is chosen in the allowed
    directions (see 'TranslationDirection'). Note that the internal
    configuration of the molecule is unchanged by this move. The maximum
    displacement is scaled during the simulation to achieve an
    acceptance ratio of 50%.

-   `"RandomTranslationProbability" : floating-point-number`
    The relative probability to attempt a random translation move for
    the current component. The displacement is chosen such that any
    position in the box can reached. It is therefore similar as
    reinsertion, but 'reinsertion' changes the internal conformation of
    a molecule and uses biasing.

-   `"RotationProbability" : floating-point-number`
    The relative probability to attempt a random rotation move for the
    current component. The rotation is around the starting bead. A
    random vector on a sphere is generated, and the rotation is random
    around this vector.

-   `"ReinsertionProbability" : floating-point-number`
    The relative probability to attempt a full reinsertion move for the
    current component. Multiple first beads are chosen, and one of these
    is selected according to its Boltzmann weight. The remaining part of
    the molecule is grown using biasing. This move is very useful, and
    often necessary, to change the internal configuration of flexible
    molecules.

-   `"SwapConventionalProbability" : floating-point-number`
    The relative probability to attempt a insertion or deletion move.
    Whether to insert or delete is decided randomly with a probability
    of 50% for each. The swap move imposes a chemical equilibrium
    between the system and an imaginary particle reservoir for the
    current component.

-   `"SwapProbability" : floating-point-number`
    The relative probability to attempt a insertion or deletion move using
    CBMC. Whether to insert or delete is decided randomly with a probability
    of 50% for each. The swap move imposes a chemical equilibrium
    between the system and an imaginary particle reservoir for the
    current component. The move starts with multiple first bead, and
    grows the remainder of the molecule using biasing.

-   `"CFCMC_SwapProbability" : floating-point-number`
    The relative probability to attempt a insertion or deletion move via
    `CFCMC` scheme.

-   `"CFCMC_CBMC_SwapProbability" : floating-point-number`
    The relative probability to attempt a insertion or deletion move via
    `CB/CFCMC` scheme.

-   `GibbsSwapProbability" : floating-point-number`
    The relative probability to attempt a Gibbs swap MC move for the
    current component. The 'GibbsSwapMove' transfers a randomly selected
    particle from one box to the other (50% probability to transfer a
    particle from box `I` to `II`, an 50% visa versa).

-   `Gibbs_CFCMC_SwapProbability" : floating-point-number`
    The relative probability to attempt a Gibbs swap MC move for the
    current component using the `CFCMC` scheme.

-   `"WidomProbability" : floating-point-number`
    The relative probability to attempt a Widom particle insertion move
    for the current component. The Widom particle insertion moves
    measure the chemical potential and can be directly related to Henry
    coefficients and heats of adsorption.
