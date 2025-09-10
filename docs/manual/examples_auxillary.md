# Examples Auxillary
\page examples_auxillary Examples Auxillary

## Table of Contents
1. [Monte Carlo: Ideal gas Rosenbluth weight butane at 300K](#Example_auxillary_1)
2. [Monte Carlo: Ideal gas Rosenbluth weight C5-C9 at 573K](#Example_auxillary_2)
3. [Monte Carlo: Ideal gas Rosenbluth weight C6-isomers at 433K](#Example_auxillary_3)
4. [Charge-equilibration IRMOF-1](#Example_auxillary_4)
5. [Grid Interpolation: CO₂ in IRMOF-1](#Example_auxillary_5)


#### Monte Carlo: Ideal gas Rosenbluth weight butane at 300K<a name="Example_auxillary_1"></a>

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 100000,
  "NumberOfInitializationCycles" : 0,
  "PrintEvery" : 5000,

  "ForceField" : ".",

  "Systems" : [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 300.0,
      "ChargeMethod" : "None"
    }
  ],

  "Components" : [
    {
      "Name" : "butane",
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

```json
{
  "CriticalTemperature" : 425.125,
  "CriticalPressure" : 3796000.0,
  "AcentricFactor" : 0.201,
  "Type" : "flexible",
  "pseudoAtoms" :
    [
      ["CH3", [0.0, 0.0, 0.0]],
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH3", [0.0, 0.0, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3]
  ],
  "Bonds" : [
    [["CH3", "CH2"], "HARMONIC", [96500.0, 1.54]],
    [["CH2", "CH2"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH3", "CH2", "CH2"], "HARMONIC", [62500.0, 114]]
  ],
  "Torsions" : [
    [["CH3", "CH2", "CH2", "CH3"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ]
}
```

\f$1.304441\times10^{-1} \pm 1.109380\times10^{-5}\f$

#### Monte Carlo: Ideal gas Rosenbluth weight C5-C9 at 573K<a name="Example_auxillary_2"></a>

```
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 100000,
  "NumberOfInitializationCycles" : 0,
  "PrintEvery" : 5000,

  "ForceField" : ".",

  "Systems" : [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 573.0,
      "ChargeMethod" : "None"
    }
  ],

  "Components" : [
    {
      "Name" : "pentane",
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "hexane",
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "heptane",
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "octane",
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "nonane",
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

| Molecule        | Rosenbluth weight                             |
|-----------------|-----------------------------------------------|
|\f$n\f$-pentane  | \f$6.393805\times10^{-2} \pm  2.626340\times10^{-5}\f$ |
|\f$n\f$-hexane   | \f$1.642211\times10^{-2} \pm  7.127763\times10^{-6}\f$ |
|\f$n\f$-heptane  | \f$4.244256\times10^{-3} \pm  1.549234\times10^{-6}\f$ |
|\f$n\f$-octane   | \f$1.103900\times10^{-3} \pm  1.057183\times10^{-6}\f$ |
|\f$n\f$-nonane   | \f$2.885495\times10^{-4} \pm  2.165096\times10^{-7}\f$ |

#### Monte Carlo: Ideal gas Rosenbluth weight C6-isomers at 433K<a name="Example_auxillary_3"></a>

```
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 100000,
  "NumberOfInitializationCycles" : 0,
  "PrintEvery" : 5000,

  "Systems" : [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 433.0,
      "ChargeMethod" : "None"
    }
  ],

  "Components" : [
    {
      "Name" : "hexane",
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "2-methylpentane",
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "3-methylpentane",
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "22-methylbutane",
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

| Molecule        | Rosenbluth weight                             |
|-----------------|-----------------------------------------------|
|hexane           | \f$8.101508\times10^{-3} \pm  4.790084\times10^{-6}\f$ |
|2-methylpentane  | \f$3.746867\times10^{-2} \pm  5.200433\times10^{-5}\f$ |
|3-methylpentane  | \f$5.350863\times10^{-2} \pm  6.526268\times10^{-5}\f$ |
|22-methylbutane  | \f$1.347238\times10^{-1} \pm  1.493745\times10^{-4}\f$ |

#### Charge-equilibration IRMOF-1<a name="Example_auxillary_4"></a>

```
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 0,
  "NumberOfInitializationCycles" : 0,
  "NumberOfEquilibrationCycles" : 0,
  "PrintEvery" : 5000,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "IRMOF-1",
      "NumberOfUnitCells" : [1, 1, 1],
      "UseChargesFrom" : "ChargeEquilibration",
      "ChargeMethod" : "Ewald",
      "ExternalTemperature" : 300.0
    }
  ],

  "Components" : [
    {
      "Name" : "methane",
      "TranslationProbability" : 1.0,
      "CreateNumberOfMolecules" :0
    }
  ]
}
```

#### Grid Interpolation: CO₂ in IRMOF-1<a name="Example_auxillary_5"></a>

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 100000,
  "NumberOfInitializationCycles" : 20000,
  "NumberOfEquilibrationCycles" : 50000,
  "PrintEvery" : 5000,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "IRMOF-1",
      "NumberOfUnitCells" : [2, 2, 2],
      "ChargeMethod" : "Ewald",
      "CutOff" : 12.0,
      "ExternalTemperature" : 298.0,
      "ExternalPressure" : 1e5
    }
  ],

  "Components" : [
    {
      "Name" : "CO2",
      "ThermodynamicIntegration" : true,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "CFCMC_CBMC_SwapProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    }
  ]
}

```

In the `force_field.json` we can set additional options to switch on the use of grid-interpolation.
```json
{
  "UseInterpolationGrids" : ["C_co2", "O_co2"],
  "SpacingVDWGrid" : 0.15,
  "SpacingCoulombGrid" : 0.15,
  "NumberOfGridTestPoints" : 100000,
  "InterpolationScheme" : 3
}
```


