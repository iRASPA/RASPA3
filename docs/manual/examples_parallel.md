# Examples Parallel
\page examples_parallel Examples Parallel

These examples use the multithreaded replica-exchange drivers. Each run declares
a single system; RASPA replicates it internally onto a temperature, pressure, or
λ ladder and exchanges configurations between neighbouring replicas.

## Table of Contents
1. [Parallel tempering: methane in a box](#Example_parallel_1)
2. [Hyper-parallel tempering: methane in MFI](#Example_parallel_2)
3. [Parallel thermodynamic integration: NaCl in water](#Example_parallel_3)
4. [Parallel TMMC: methane vapour–liquid equilibrium](#Example_parallel_4)

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
