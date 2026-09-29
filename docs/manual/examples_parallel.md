# Examples Parallel
\page examples_parallel Examples Parallel

These examples use the multithreaded drivers. The replica-exchange drivers
(examples 1–10) declare a single system that RASPA replicates internally onto
a temperature, pressure, or λ ladder and exchange configurations between
neighbouring replicas; the spatial-decomposition MD driver (example 11) runs
one large system with the threads owning sub-domains of the box.

## Table of Contents
1. [Parallel tempering: methane in a box](#Example_parallel_1)
2. [Hyper-parallel tempering: methane in MFI](#Example_parallel_2)
3. [Parallel thermodynamic integration: NaCl in water](#Example_parallel_3)
4. [Parallel TMMC: methane vapour–liquid equilibrium](#Example_parallel_4)
5. [WHAM-CBMC: nitrogen BET in MFI](#Example_parallel_5)
6. [WHAM-CBCFCMC: nitrogen BET in MFI](#Example_parallel_6)
7. [TMMC-CBMC: nitrogen BET in MFI](#Example_parallel_7)
8. [TMMC-CBCFCMC: nitrogen BET in MFI](#Example_parallel_8)
9. [Replica-exchange molecular dynamics: pHDDA-20 end-to-end distance](#Example_parallel_9)
10. [Replica-exchange molecular dynamics: end-to-end distance of liquid HDDA](#Example_parallel_10)
11. [Spatial-decomposition molecular dynamics: 2000 HDDA monomers](#Example_parallel_11)

----------------------------------------------------------------------------------

#### Parallel tempering: methane in a box <a name="Example_parallel_1"></a>

A parallel-tempering Monte Carlo simulation of 100 methane molecules in a
\f$30 \times 30 \times 30\f$ &Aring; box. The system is replicated onto eight
temperatures from 240 K to 380 K. Neighbouring replicas exchange configurations
every 10 cycles.

Run from `examples/parallel/1_parallel_tempering_methane_in_box`:

```json
{
  "SimulationType" : "ParallelTempering",
  "ParallelTemperingSwapEvery" : 10,
  "NumberOfThreads" : 1,
  "NumberOfInitializationCycles" : 1000,
  "NumberOfEquilibrationCycles" : 0,
  "NumberOfProductionCycles" : 10000,
  "PrintEvery" : 1000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperatures" : [240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0],
      "ChargeMethod" : "None"
    }
  ],

  "Components" :
  [
    {
      "Name" : "methane",
      "MoleculeDefinition" : "ExampleDefinitions",
      "TranslationProbability" : 1.0,
      "CreateNumberOfMolecules" : 100
    }
  ]
}
```

Leave `"NumberOfThreads"` at 1 so the per-energy-evaluation thread pool stays
serial; the driver already spawns one worker thread per temperature.

#### Hyper-parallel tempering: methane in MFI <a name="Example_parallel_2"></a>

Hyper-parallel tempering of methane adsorption in siliceous MFI over a
\f$3 \times 3\f$ grid of temperatures (300, 350, 400 K) and pressures
(\f$10^4\f$, \f$10^5\f$, \f$10^6\f$ Pa). Swaps alternate between the temperature
and pressure directions. Directly measured isotherms are written per temperature
to `output/isotherm_{T}.hyper_parallel_tempering.txt`.

Run from `examples/parallel/2_hyper_parallel_tempering_methane_in_mfi`:

```json
{
  "SimulationType" : "HyperParallelTempering",
  "NumberOfInitializationCycles" : 10000,
  "NumberOfProductionCycles" : 100000,
  "PrintEvery" : 5000,
  "ParallelTemperingSwapEvery" : 10,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "MFI_SI",
      "NumberOfUnitCells" : [2, 2, 2],
      "ExternalTemperatures" : [300.0, 350.0, 400.0],
      "ExternalPressures" : [1.0e4, 1.0e5, 1.0e6],
      "ChargeMethod" : "None"
    }
  ],

  "Components" : [
    {
      "Name" : "methane",
      "IdealGasRosenbluthWeight" : 1.0,
      "TranslationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "SwapProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

#### Parallel thermodynamic integration: NaCl in water <a name="Example_parallel_3"></a>

The ion-pair coupling parameter λ of NaCl in 300 water molecules is sampled on
16 λ-bins in one multithreaded run. Neighbouring replicas exchange λ values
every 10 cycles. The stitched \f$\langle \partial U/\partial\lambda \rangle\f$
curve is integrated from λ = 0 to 1 to give the excess chemical potential.

Run from `examples/parallel/3_parallel_ti_nacl_in_water`:

```json
{
  "SimulationType" : "ParallelThermodynamicIntegration",
  "NumberOfLambdaBins" : 16,
  "LambdaExchangeEvery" : 10,
  "NumberOfInitializationCycles" : 1000,
  "NumberOfEquilibrationCycles" : 100000,
  "NumberOfProductionCycles" : 100000,
  "PrintEvery" : 5000,

  "Systems" : [
    {
      "Type" : "Box",
      "BoxLengths" : [20.8, 20.8, 20.8],
      "ExternalTemperature" : 298.0,
      "ExternalPressure" : 1.0e5,
      "ChargeMethod" : "Ewald",
      "CutOff" : 10.0,
      "VolumeMoveProbability" : 0.01,
      "HybridMCProbability" : 0.001,
      "HybridMCMoveNumberOfSteps" : 10
    }
  ],

  "Components" : [
    {
      "Name" : "water",
      "FugacityCoefficient" : 1.0,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "CreateNumberOfMolecules" : 300
    },
    {
      "Name" : "sodium",
      "FugacityCoefficient" : 1.0,
      "PairComponent" : 2,
      "ThermodynamicIntegration" : true,
      "TranslationProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "chloride",
      "FugacityCoefficient" : 1.0,
      "PairComponent" : 1,
      "TranslationProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

#### Parallel TMMC: methane vapour–liquid equilibrium <a name="Example_parallel_4"></a>

Windowed transition-matrix Monte Carlo of methane in a
\f$30 \times 30 \times 30\f$ &Aring; box at 160, 170 and 180 K. The macrostate
range \f$N = 0\ldots 420\f$ is split into 6 overlapping windows. The driver
reweights the collected ln Π(N) over a pressure range to obtain the
vapour–liquid coexistence.

Run from `examples/parallel/4_parallel_tmmc_vle_methane`:

```json
{
  "SimulationType" : "ParallelTMMC",
  "NumberOfInitializationCycles" : 2000,
  "NumberOfEquilibrationCycles" : 20000,
  "NumberOfProductionCycles" : 100000,
  "PrintEvery" : 5000,
  "NumberOfWindows" : 6,
  "TMMCUpdateEvery" : 100000,
  "ReweightingPressureRange" : [2.0e5, 8.0e6],

  "Systems" : [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperatures" : [160.0, 170.0, 180.0],
      "ExternalPressure" : 2.0e6,
      "ChargeMethod" : "None",
      "MacroStateMinimumNumberOfMolecules" : 0,
      "MacroStateMaximumNumberOfMolecules" : 420
    }
  ],

  "Components" : [
    {
      "Name" : "methane",
      "IdealGasRosenbluthWeight" : 1.0,
      "TranslationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "SwapProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

#### WHAM-CBMC: nitrogen BET in MFI <a name="Example_parallel_5"></a>

Nitrogen at 77.355 K in siliceous MFI, sampled with CBMC swaps on a WHAM
replica grid. `"ComputeBET" : true` extracts a Rouquerol BET area from the
reweighted isotherm and places the pressure ladder and reweighting grid
automatically from a Widom Henry coefficient up to P0 (101325 Pa). Short GCMC
probes at Langmuir θ = 0.1, 0.5, 0.9 select Langmuir, Langmuir–Freundlich, or
Toth; rungs are then re-placed at equal Fisher overlap on a 2× equal-log
skeleton plus at most two extra rungs per interval.
`"MacroStateMaximumNumberOfMolecules" : "auto"` scouts occupancy at P0 for
the filling ceiling (WHAM cycle length / TMMC N_max). `"NumberOfThreads"`
sizes the ladder (8–32 rungs). Analysis files go to `wham/`.

Run from `examples/parallel/5_parallel_wham_cbmc_n2_bet_in_mfi`:

```json
{
  "SimulationType" : "ReweightedHistogram",
  "NumberOfInitializationCycles" : 5000,
  "NumberOfEquilibrationCycles" : 5000,
  "NumberOfProductionCycles" : 40000,
  "PrintEvery" : 5000,
  "ParallelTemperingSwapEvery" : 10,
  "SampleReweightingEvery" : 5,
  "ReweightingTemperatures" : [77.355],
  "ReweightingPressureRange" : "auto",
  "NumberOfThreads" : 16,
  "ComputeBET" : true,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "MFI_SI",
      "NumberOfUnitCells" : [2, 2, 2],
      "ExternalTemperatures" : [77.355],
      "ExternalPressures" : "auto",
      "ChargeMethod" : "Ewald",
      "MacroStateMinimumNumberOfMolecules" : 0,
      "MacroStateMaximumNumberOfMolecules" : "auto"
    }
  ],

  "Components" : [
    {
      "Name" : "N2",
      "IdealGasRosenbluthWeight" : 1.0,
      "CrossSection" : 16.2,
      "LiquidVolume" : 57.7,
      "SaturationPressure" : 101325.0,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "SwapProbability" : 1.0,
      "BlockingPockets" : "auto",
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

#### WHAM-CBCFCMC: nitrogen BET in MFI <a name="Example_parallel_6"></a>

The same nitrogen BET measurement as example 5, with CB/CFCMC swaps and a
longer equilibration so the Wang–Landau λ-bias can flatten before production.

Run from `examples/parallel/6_parallel_wham_cbcfcmc_n2_bet_in_mfi`:

```json
{
  "SimulationType" : "ReweightedHistogram",
  "NumberOfInitializationCycles" : 5000,
  "NumberOfEquilibrationCycles" : 20000,
  "NumberOfProductionCycles" : 40000,
  "PrintEvery" : 5000,
  "ParallelTemperingSwapEvery" : 10,
  "SampleReweightingEvery" : 5,
  "ReweightingTemperatures" : [77.355],
  "ReweightingPressureRange" : "auto",
  "NumberOfThreads" : 16,
  "ComputeBET" : true,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "MFI_SI",
      "NumberOfUnitCells" : [2, 2, 2],
      "ExternalTemperatures" : [77.355],
      "ExternalPressures" : "auto",
      "ChargeMethod" : "Ewald",
      "MacroStateMinimumNumberOfMolecules" : 0,
      "MacroStateMaximumNumberOfMolecules" : "auto"
    }
  ],

  "Components" : [
    {
      "Name" : "N2",
      "IdealGasRosenbluthWeight" : 1.0,
      "CrossSection" : 16.2,
      "LiquidVolume" : 57.7,
      "SaturationPressure" : 101325.0,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "CFCMC_CBMC_SwapProbability" : 1.0,
      "BlockingPockets" : "auto",
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

#### TMMC-CBMC: nitrogen BET in MFI <a name="Example_parallel_7"></a>

Windowed TMMC of the same N2/MFI system. The walk samples at P0; `"ComputeBET"`
places only the reweighting grid from the Henry coefficient to P0. The driver
writes equilibrium, adsorption, and desorption isotherms (and BET-plot tables)
under `tmmc/`.

Run from `examples/parallel/7_parallel_tmmc_cbmc_n2_bet_in_mfi`:

```json
{
  "SimulationType" : "ParallelTMMC",
  "NumberOfInitializationCycles" : 5000,
  "NumberOfEquilibrationCycles" : 20000,
  "NumberOfProductionCycles" : 40000,
  "PrintEvery" : 5000,
  "NumberOfWindows" : 16,
  "NumberOfThreads" : 16,
  "TMMCUpdateEvery" : 10000,
  "ReweightingPressureRange" : "auto",
  "ComputeBET" : true,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "MFI_SI",
      "NumberOfUnitCells" : [2, 2, 2],
      "ExternalTemperature" : 77.355,
      "ExternalPressure" : 1.01325e5,
      "ChargeMethod" : "Ewald",
      "MacroStateMinimumNumberOfMolecules" : 0,
      "MacroStateMaximumNumberOfMolecules" : "auto"
    }
  ],

  "Components" : [
    {
      "Name" : "N2",
      "IdealGasRosenbluthWeight" : 1.0,
      "CrossSection" : 16.2,
      "LiquidVolume" : 57.7,
      "SaturationPressure" : 101325.0,
      "TranslationProbability" : 1.0,
      "RotationProbability" : 1.0,
      "ReinsertionProbability" : 2.0,
      "SwapProbability" : 2.0,
      "BlockingPockets" : "auto",
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

#### TMMC-CBCFCMC: nitrogen BET in MFI <a name="Example_parallel_8"></a>

The same TMMC nitrogen BET measurement as example 7, with CB/CFCMC swaps so
the collection matrix is the flattened (N, λ) chain.

Run from `examples/parallel/8_parallel_tmmc_cbcfcmc_n2_bet_in_mfi`:

```json
{
  "SimulationType" : "ParallelTMMC",
  "NumberOfInitializationCycles" : 5000,
  "NumberOfEquilibrationCycles" : 20000,
  "NumberOfProductionCycles" : 40000,
  "PrintEvery" : 5000,
  "NumberOfWindows" : 16,
  "NumberOfThreads" : 16,
  "TMMCUpdateEvery" : 10000,
  "ReweightingPressureRange" : "auto",
  "ComputeBET" : true,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "MFI_SI",
      "NumberOfUnitCells" : [2, 2, 2],
      "ExternalTemperature" : 77.355,
      "ExternalPressure" : 1.01325e5,
      "ChargeMethod" : "Ewald",
      "MacroStateMinimumNumberOfMolecules" : 0,
      "MacroStateMaximumNumberOfMolecules" : "auto"
    }
  ],

  "Components" : [
    {
      "Name" : "N2",
      "IdealGasRosenbluthWeight" : 1.0,
      "CrossSection" : 16.2,
      "LiquidVolume" : 57.7,
      "SaturationPressure" : 101325.0,
      "TranslationProbability" : 1.0,
      "RotationProbability" : 1.0,
      "ReinsertionProbability" : 2.0,
      "CFCMC_CBMC_SwapProbability" : 2.0,
      "BlockingPockets" : "auto",
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

#### Replica-exchange molecular dynamics: pHDDA-20 end-to-end distance <a name="Example_parallel_9"></a>

Replica-exchange molecular dynamics (REMD, Sugita & Okamoto 1999) of a single
flexible pHDDA-20 chain in a \f$60 \times 60 \times 60\f$ &Aring; box, the MD
counterpart of the parallel-tempering Monte Carlo run in
`examples/polymers/4_mc_parallel_tempering_end_to_end_distance_phdda_20_in_box`
and of the plain NVT MD run in
`examples/polymers/3_md_end_to_end_distance_phdda_20_in_box`. The chain is
replicated onto a geometric ladder of 16 temperatures from 300 K to 500 K; each
replica is integrated with its own Nosé–Hoover chain, and every 500 time steps
neighbouring replicas attempt to exchange their configurations. After an
accepted exchange the momenta are rescaled by \f$\sqrt{T_\text{new}/T_\text{old}}\f$,
so every replica stays canonical at its own temperature while the chain
conformations diffuse through the ladder and the low-temperature replicas
escape the collapsed states that trap a single MD trajectory. The MC moves of
the component are used only in the initialization stage to relax the chain
before the dynamics starts; the molecule-property histograms (end-to-end
distance, radius of gyration) are written per replica.

Run from `examples/polymers/5_md_parallel_tempering_end_to_end_distance_phdda_20_in_box`:

```json
{
  "SimulationType" : "ParallelTemperingMolecularDynamics",
  "ParallelTemperingSwapEvery" : 500,
  "NumberOfInitializationCycles" : 2000,
  "NumberOfEquilibrationCycles" : 20000,
  "NumberOfProductionCycles" : 1000000,
  "NumberOfThreads" : 1,
  "PrintEvery" : 10000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [60.0, 60.0, 60.0],
      "ExternalTemperatures" : [300.0, 310.4, 321.1, 332.3, 343.8, 355.7, 368.0, 380.8,
                                393.9, 407.6, 421.7, 436.3, 451.4, 467.1, 483.3, 500.0],
      "Ensemble" : "NVT",
      "TimeStep" : 0.001,
      "ChargeMethod" : "Ewald",
      "ComputeMoleculeProperties" : true,
      "SampleMoleculePropertiesEvery" : 10,
      "WriteMoleculePropertiesEvery" : 100000,
      "NumberOfBinsMoleculeProperties" : 128
    }
  ],

  "Components" :
  [
    {
      "Name" : "pHDDA-20",
      "PivotProbability" : 1.0,
      "CrankshaftProbability" : 1.0,
      "BeadFlipProbability" : 1.0,
      "BeadDisplacementProbability" : 1.0,
      "CreateNumberOfMolecules" : 1
    }
  ]
}
```

For MD a cycle is one time step, so `"ParallelTemperingSwapEvery"` is set
well above the default of 10. The combined output file
`output/output.parallel_tempering_md.txt` reports the exchange acceptance per
neighbouring pair (a pair with a low acceptance marks a gap in the ladder) and
the potential- and conserved-energy drift of every replica; the per-replica
files `output/output_{T}_0.parallel_tempering_md.r{k}.txt` carry the usual MD
status reports with the kinetic temperatures, which should average to the
ladder temperature of the replica.

#### Replica-exchange molecular dynamics: end-to-end distance of liquid HDDA <a name="Example_parallel_10"></a>

Reproduces Fig. S5 of Torres-Knoop, Kryven, Schamboeck and Iedema, *Soft
Matter* **14**, 3404 (2018): the histogram of the end-to-end distance of
1,6-hexanediol diacrylate (HDDA) monomers in the liquid at 300 K, the
observation behind the small-cycle mechanism of that paper. In the liquid the
monomer is mostly in a semi-coiled state, but it coils (vinyl groups about
4 &Aring; apart) and uncoils (about 15 &Aring;) across a free-energy barrier that
the paper puts at 5–10 kcal/mol by metadynamics. At 300 K a single MD
trajectory crosses such a barrier only every few to hundreds of nanoseconds
per molecule (depending on where in that range the barrier lies); the
replica-exchange ladder lets the conformations equilibrate at up to 450 K,
where the crossing is one to two orders of magnitude faster, and diffuse back
to 300 K.

The monomer `HDDA.json` is the TraPPE-UA acrylate model of Maerzke et al. with
harmonic bonds (Tables S1–S4 of the paper's ESI): 16 united atoms,
CH2=CH–C(=O)–O–CH2–(CH2)4–CH2–O–C(=O)–CH=CH2, with the ester charges,
1-4 Coulomb scaling of 0.5 and the same force field file as the pHDDA
examples. The end-to-end distance is measured between the two terminal vinyl
carbons (`"EndToEndAtoms" : [0, 13]`), and the histogram is written on the
0–20 &Aring; axis of the figure (`"EndToEndRangeMoleculeProperties" : 20.0`, 80
bins of 0.25 &Aring;). The box holds 100 monomers at the experimental density of
1.025 g/cm\f$^3\f$ that the paper reproduces (\f$L = 33.22\f$ &Aring;; the paper
used 2000 monomers in a 100 &Aring; box, but the end-to-end distribution is an
intramolecular property and converges with the number of samples, not with
the box size). The 16 temperatures form a geometric ladder from 300 K to
450 K; with about 4800 degrees of freedom the constant-\f$C_v\f$ estimate gives
an exchange acceptance of roughly 20 % per pair. Exchanges are attempted every
200 steps (0.2 ps): an attempt costs nothing (the energies are known), and the
rate at which configurations travel the ladder is what limits the sampling at
300 K, where the coiling barrier itself is crossed only rarely.

Run from `examples/polymers/6_md_parallel_tempering_end_to_end_distance_hdda_liquid`:

```json
{
  "SimulationType" : "ParallelTemperingMolecularDynamics",
  "ParallelTemperingSwapEvery" : 200,
  "NumberOfInitializationCycles" : 5000,
  "NumberOfEquilibrationCycles" : 1000000,
  "NumberOfProductionCycles" : 10000000,
  "NumberOfThreads" : 16,
  "PrintEvery" : 100000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [33.22, 33.22, 33.22],
      "ExternalTemperatures" : [300.0, 308.2, 316.7, 325.3, 334.3, 343.4, 352.8, 362.5,
                                372.4, 382.6, 393.1, 403.9, 414.9, 426.3, 438.0, 450.0],
      "Ensemble" : "NVT",
      "TimeStep" : 0.001,
      "ChargeMethod" : "Ewald",
      "ComputeMoleculeProperties" : true,
      "SampleMoleculePropertiesEvery" : 10,
      "WriteMoleculePropertiesEvery" : 1000000,
      "NumberOfBinsMoleculeProperties" : 80,
      "EndToEndRangeMoleculeProperties" : 20.0
    }
  ],

  "Components" :
  [
    {
      "Name" : "HDDA",
      "TranslationProbability" : 1.0,
      "RotationProbability" : 1.0,
      "ReinsertionProbability" : 0.5,
      "PartialReinsertionProbability" : 1.0,
      "CrankshaftProbability" : 1.0,
      "BeadDisplacementProbability" : 1.0,
      "CreateNumberOfMolecules" : 100
    }
  ]
}
```

The 100 monomers are grown into the box with CBMC and relaxed by the Monte
Carlo moves of the initialization stage (the only stage in which these moves
are used), then integrated for 1 ns of equilibration and 10 ns of production
per replica (of the order of two days on 16 cores). The equilibration has to
do more than settle the packing and the thermostat: the CBMC-grown monomers
start with the conformational distribution of an isolated chain, and the
liquid one is reached through the exchanges (the equilibration stage also
attempts them) at about the round-trip time, so it is given roughly one to
three round trips. The length is set by the
coiled state: it is a minority population that the 300 K replica only receives
through the ladder, and estimating its weight to about 10 % needs of the order
of ten round trips, i.e. ten times the round-trip time of roughly 0.3–1 ns that
16 replicas at 20 % acceptance and 0.2 ps between attempts give (the
`round trips` count in the combined output shows the progress). A shorter run
already gives the shape of the semi-coiled main peak, and a run can be
extended from the binary restart file. The histogram of the 300 K replica,
`molecule_properties/end_to_end_HDDA_0_13.s0.txt` (column 1 the distance,
column 2 the normalized probability density, column 3 its 95 % confidence
error), is the quantity of Fig. S5; the header gives \f$\langle R\rangle\f$ and
\f$\langle R^2\rangle\f$. The files `.s1` … `.s15` are the same histogram at
the higher temperatures and show the coiled population grow with temperature.
The combined output `output/output.parallel_tempering_md.txt` reports the
exchange acceptance per pair, the round-trip statistics and the energy drift.
For this observable alone, Monte Carlo with regrowth moves (as in the
parallel-tempering MC example of the pHDDA-20 chain) is the cheaper route,
since a regrowth jumps between the coiled and uncoiled states without
crossing the barrier; the MD version is the one that is comparable with the
paper.

#### Spatial-decomposition molecular dynamics: 2000 HDDA monomers <a name="Example_parallel_11"></a>

The system of the paper behind the previous example at its original size:
2000 HDDA monomers (32 000 united atoms) at the experimental density of
1.025 g/cm\f$^3\f$ in a \f$90.2 \times 90.2 \times 90.2\f$ &Aring; box at
300 K, integrated with `"SimulationType" :
"MolecularDynamicsSpatialDecomposition"`. The replica-exchange example
parallelizes over temperatures and needs many small replicas; this driver
parallelizes a single system over space. The box is divided into cells
(a fraction of the cutoff plus the Verlet skin, \f$14 + 2 = 16\f$ &Aring;)
and the cells are distributed over
`"NumberOfThreads"` sub-domains, one persistent thread each; a thread builds
the Verlet lists of its own atoms, computes their short-range pair forces,
spreads their charges onto a private copy of the particle-mesh Ewald grid
and interpolates the reciprocal forces back, and handles the bonded terms and
Ewald corrections of its share of the molecules, with barriers between the
phases. The threads never write to the same atom, so the result is
independent of the thread count to rounding order. The cells are chosen a
quarter of the list cutoff wide here (\f$22 \times 22 \times 22\f$ cells,
\f$9^3\f$ stencil), because that lets the sub-domain boundaries fall so that
all 16 threads own nearly the same number of cells; with cells of the full
list cutoff (\f$5 \times 5 \times 5\f$) the largest sub-domain would hold four
times the atoms of the smallest. The reciprocal sum is a
particle-mesh (SPME/PPPM) solver with the same \f$\alpha\f$ as the Ewald
summation of the force field; with a 1 &Aring; mesh (\f$90^3\f$ points) and
fifth-order B-splines the reciprocal energy and forces agree with the exact
Ewald sum to about \f$10^{-4}\f$ relative, which the output reports at the
start of each MD stage (`Spatial-decomposition force engine check`).

Run from `examples/polymers/7_md_spatial_decomposition_end_to_end_distance_hdda_liquid_2000`:

```json
{
  "SimulationType" : "MolecularDynamicsSpatialDecomposition",
  "NumberOfThreads" : 12,
  "VerletSkin" : 2.0,
  "PPPMMeshSpacing" : 1.0,
  "PPPMInterpolationOrder" : 5,
  "NumberOfInitializationCycles" : 200,
  "NumberOfEquilibrationCycles" : 200000,
  "NumberOfProductionCycles" : 2000000,
  "PrintEvery" : 10000,
  "WriteBinaryRestartEvery" : 10000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [90.2, 90.2, 90.2],
      "ExternalTemperature" : 300.0,
      "Ensemble" : "NVT",
      "TimeStep" : 0.001,
      "ChargeMethod" : "Ewald",
      "ComputeMoleculeProperties" : true,
      "SampleMoleculePropertiesEvery" : 10,
      "WriteMoleculePropertiesEvery" : 100000,
      "NumberOfBinsMoleculeProperties" : 80,
      "EndToEndRangeMoleculeProperties" : 20.0
    }
  ],

  "Components" :
  [
    {
      "Name" : "HDDA",
      "TranslationProbability" : 1.0,
      "RotationProbability" : 1.0,
      "ReinsertionProbability" : 0.5,
      "PartialReinsertionProbability" : 1.0,
      "CrankshaftProbability" : 1.0,
      "BeadDisplacementProbability" : 1.0,
      "CreateNumberOfMolecules" : 2000
    }
  ]
}
```

The force field is that of the pHDDA examples with one change: the Coulomb
cutoff is set explicitly (`"CutOffCoulomb" : 14.0`) instead of `"auto"`,
because the automatic Coulomb cutoff is half the box (45 &Aring; here), which
would make the pair lists enormous and leave no room for the Verlet skin
(the neighbour lists need cutoff + skin at most half the box). The
initialization stage is serial Monte Carlo on 32 000 atoms and is kept short
(200 cycles, about a quarter of an hour, after about a minute of CBMC growth
of the 2000 monomers); it only has to remove the overlaps of the CBMC-grown
monomers, the 0.2 ns of MD equilibration does the rest. The
production stage of 2 ns samples the end-to-end histogram
`molecule_properties/end_to_end_HDDA_0_13.s0.txt` of the previous example
from 2000 monomers instead of 100 (a run can be extended from the binary
restart file; the engine settings are always taken from the input file, so
a restart may change the thread count).

Thread scaling of the force evaluation for this system (32 000 atoms,
cutoff 14 &Aring; plus 2 &Aring; skin, \f$96^3\f$ mesh, fifth-order splines;
Apple M4 Max, 12 performance and 4 efficiency cores; wall time per force
evaluation averaged over 50 MD steps including the neighbour-list rebuilds,
one rebuild per 17 steps):

| threads | sub-domains | ms per force evaluation | speed-up |
|--------:|:------------|------------------------:|---------:|
| exact code (all pairs + direct Ewald sum, serial) | – | 2350 | – |
| 1 | 1 × 1 × 1 | 96 | 1.0 |
| 2 | 1 × 1 × 2 | 51 | 1.9 |
| 4 | 1 × 2 × 2 | 27 | 3.6 |
| 8 | 2 × 2 × 2 | 15.1 | 6.3 |
| 12 | 2 × 2 × 3 | 11.2 | 8.6 |
| 16 | 2 × 2 × 4 | 14.0 | 6.9 |

The single-thread engine is 24 times faster than the exact code because the
Verlet lists replace the \f$O(N^2)\f$ pair loop, the mesh replaces the
direct k-space sum, and the pair kernel works on a compact copy of the atoms
with the periodic shifts resolved in the neighbour list (no minimum-image
operation per pair) using the specialised Lennard-Jones + tabulated-erfc
kernel (about 5.5 ns per neighbour-list pair). Every pair is evaluated once:
the force on a ghost image is stored in the thread's private force buffer and
collected by the owner after the pair phase, and the sub-domain cuts are
placed on the atom positions so every thread owns exactly the same number of
atoms. The parallel efficiency is then set by the barriers and the serial
parts of the list rebuild (the 4 efficiency cores do not add to the 12
performance cores because every phase waits for the slowest thread; set
`"NumberOfThreads"` to the number of equal cores). Per step on 12 threads the
pair forces take about 55 % of the time, the mesh (spreading, reduction,
FFTs, interpolation) 22 %, the list rebuilds 15 % and the bonded terms and
Ewald corrections 7 %. At 11 ms per step the 2 ns of production take 6 hours
on 12 cores; with `"SimulationType" : "MolecularDynamics"` the same
trajectory would take over a month.

For reference, LAMMPS (22 Jul 2025, FFTW3 build) on the same machine and a
system with the same pair and mesh workload (32 768 Lennard-Jones + charge
atoms in a 90.2 &Aring; box, `lj/cut/coul/long 14.0`, `neighbor 2.0 bin`,
`pppm 1e-4` with `mesh 96 96 96 order 5 diff ad`) needs 130 ms per step on
one core and 14.8 ms per step on 12 MPI ranks (172 ms with the default `ik`
PPPM differentiation); its pair phase takes 73 ms and 6.9 ms per step,
against 66 ms and 6.2 ms here.
