# Examples Auxillary
\page examples_auxillary Examples Auxillary

## Table of Contents
1. [Monte Carlo: Ideal gas Rosenbluth weight butane at 300K](#Example_auxillary_1)
2. [Monte Carlo: Ideal gas Rosenbluth weight C5-C9 at 573K](#Example_auxillary_2)
3. [Monte Carlo: Ideal gas Rosenbluth weight C6-isomers at 433K](#Example_auxillary_3)
4. [Charge-equilibration IRMOF-1](#Example_auxillary_4)
5. [Grid Interpolation: CO₂ in IRMOF-1](#Example_auxillary_5)


#### Monte Carlo: Ideal gas Rosenbluth weight butane at 300K<a name="Example_auxillary_1"></a>

To compare simulation values to experiments a reference state should be chosen. A convenient reference state is the ideal gas. The reference Rosenbluth value can be computed from a simulation of a single chain in an empty box at the desired temperature. 
```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 100000,
  "NumberOfInitializationCycles" : 0,
  "PrintEvery" : 5000,

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
Note that the ideal-gas Rosenbluth weight depends on the temperatures _and_ on the force field.
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
For this butane-model, we obtain at 300K an ideal-gas Rosenbluth weight of \f$1.304441\times10^{-1} \pm 1.109380\times10^{-5}\f$.

#### Monte Carlo: Ideal gas Rosenbluth weight C5-C9 at 573K<a name="Example_auxillary_2"></a>

Note that for Rosenbluth weights several chains can be computed simultaneously, since they are computed from Widom insertions where the molecule is never actually inserted in the system.

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

Similarly, we can obtain the ideal-gas Rosenbluth weights of hexane isomers.
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
    },
    {
      "Name" : "23-methylbutane",
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```
Note the ideal-gas Rosenbluth weights can vary a lot, even for isomers.

| Molecule        | Rosenbluth weight                             |
|-----------------|-----------------------------------------------|
|hexane           | \f$8.103214\times10^{-3} \pm  1.092925\times10^{-6}\f$ |
|2-methylpentane  | \f$4.703641\times10^{-2} \pm  1.562648\times10^{-5}\f$ |
|3-methylpentane  | \f$5.351406\times10^{-2} \pm  1.976970\times10^{-5}\f$ |
|22-methylbutane  | \f$2.265669\times10^{-1} \pm  2.553789\times10^{-5}\f$ |
|23-methylbutane  | \f$8.730600\times10^{-2} \pm  1.510154\times10^{-5}\f$ |

#### Charge-equilibration IRMOF-1<a name="Example_auxillary_4"></a>

TODO: not working yet ([charge equilibration]: no solution found').
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
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

#### Grid Interpolation: CO₂ in IRMOF-1<a name="Example_auxillary_5"></a>

Some simulations, especially with a large number of unitcells, by using grid-interpolation. Consider a simulation of CO₂, then we can make a pre-computed energy grid for the oxygen and the carbon of the molecule. An additional grid is needed for the electrostatics.
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
InterpolationScheme 1 is capable of interpolating only energies, InterpolationScheme 3 interpolates the energies and the forces, while InterpolationScheme 5 interpolates energies, forces, and the Hessian.

For each listed (pseudo-)atom, a grid will be made at the start of the simulation.
```
Generating an Ewald Real interpolation grid (172x172x172) for a unit charge
===============================================================================
(Using 1 1 1 periodic cells to create grid)
Percentage finished: 1
Percentage finished: 2
Percentage finished: 3
....
Percentage finished: 99
Percentage finished: 100
Grid done... (     20.537806 [s])

Generating an VDW interpolation grid (172x172x172) for C_co2
===============================================================================
(Using 1 1 1 periodic cells to create grid)
Percentage finished: 1
Percentage finished: 2
Percentage finished: 3
....
Percentage finished: 99
Percentage finished: 100
Grid done... (      6.685936 [s])
```

After grid creation, the grids are tested with random insertions
```
Testing VDW interpolation grid (172x172x172) for C_co2
-------------------------------------------------------------------------------
(Using 100000 points for testing)

Boltzmann average energy VDW (table):      -214.02840045285805
Boltzmann average energy VDW (full):       -214.01140960167118
Boltzmann relative error:                  0.00021947314459359044

Boltzmann average gradient(x) VDW (table): 3.712850362496449
Boltzmann average gradient(x) VDW (full):  3.713806990998235
Boltzmann relative error:                  0.002155151736540394

Boltzmann average gradient(y) VDW (table): 0.8850842781059334
Boltzmann average gradient(y) VDW (full):  0.8857628567145862
Boltzmann relative error:                  0.002131945045696959

Boltzmann average gradient(z) VDW (table): 0.6725652417255391
Boltzmann average gradient(z) VDW (full):  0.671549051697825
Boltzmann relative error:                  0.0020873926276597484

Testing Coulomb interpolation grid (172x172x172) for C_co2
-------------------------------------------------------------------------------
(Using 100000 points for testing)

Boltzmann average energy Real Ewald (table):      -192.95353813520873
Boltzmann average energy Real Ewald (full):       -192.95384459189873
Boltzmann relative error:                         5.643361809478482e-05

Boltzmann average gradient(x) Real Ewald (table): -1.115591762633454
Boltzmann average gradient(x) Real Ewald (full):  -1.1181301909005674
Boltzmann relative error:                         0.0009199018730806105

Boltzmann average gradient(y) Real Ewald (table): 2.994373232849875
Boltzmann average gradient(y) Real Ewald (full):  2.994950182057741
Boltzmann relative error:                         0.000927688382477166

Boltzmann average gradient(z) Real Ewald (table): 5.0800230953098415
Boltzmann average gradient(z) Real Ewald (full):  5.0797601170886315
Boltzmann relative error:                         0.0009238132876741674
```

Note that grids also work with the CFCMC techniques. Here, the fractional molecules are computed using full energy evaluations, while the energies of all integer molecule are interpolated.
