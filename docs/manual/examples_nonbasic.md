# Examples Non-Basic
\page examples_nonbasic Examples Non-Basic


## Table of Contents
1. [Monte Carlo: adsorption of binary mixture of CO₂ and methane in IRMOF-1](#Example_nonbasic_1)
2. [Monte Carlo: NPT CO₂](#Example_nonbasic_2)
3. [Monte Carlo: NPT propane](#Example_nonbasic_3)
4. [Monte Carlo: Gibbs CO₂](#Example_nonbasic_4)
5. [Molecular Dynamics: benzene in IRMOF-1](#Example_nonbasic_5)
6. [Monte Carlo: adsorption CO₂ in LTA-4A sodium](#Example_nonbasic_6)
7. [Molecular Dynamics: diffusion of CO₂ in LTA-4A sodium](#Example_nonbasic_7)
8. [Monte Carlo: adsorption of C6-isomers in MFI](#Example_nonbasic_8)



#### Monte Carlo: adsorption of binary mixture of CO₂ and methane in IRMOF-1<a name="Example_nonbasic_1"></a>

Adsorption of a binary mixture is specified at a total pressure and individual mol-fractions for the components.
```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 500000,
  "NumberOfInitializationCycles" : 50000,
  "PrintEvery" : 1000,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "IRMOF-1",
      "NumberOfUnitCells" : [1, 1, 1],
      "HeliumVoidFraction" : 0.81,
      "ChargeMethod" : "Ewald",
      "ExternalTemperature" : 300.0,
      "ExternalPressure" : 1e6
    }
  ],

  "Components" : [
    {
      "Name" : "CO2",
      "MolFraction" : 0.25,
      "FugacityCoefficient" : 1.0,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "SwapProbability" : 1.0,
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "methane",
      "MolFraction" : 0.75,
      "FugacityCoefficient" : 1.0,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "SwapProbability" : 1.0,
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

The fugacity \f$f_i\f$ of component \f$i\f$ in a gas-mixture is given by
$$f_i = \gamma_i p_i$$
where \f$\gamma_i\f$ is the fugacity coefficient and \f$p_i\f$ the partial pressure of component \f$i\f$.

The use of fugacity is prefered over using chemical potentials. For mixtures, the chemical potential diverges when the molar fraction of a component goes to zero. In contrast, the fugacity goes to zero. In addition, there is a very natural reference state for fugacity, namely an ideal gas. The fugacity and pressure converge at low pressure
$$\lim_{p\rightarrow0}\frac{f_i}{p_i} = 1$$

Using the Peng-Robinson equation of the state, the fugacity coefficients and the number of excess molecules are computed:
```
Component 0 [CO2]
    Mol-fraction:                 0.25 [-]
    Fugacity:                     237587.6177345715 [Pa]
    Fugacity coefficient:         0.9503504709382861 [-]
    Compressibility:              0.9714389724700717 [-]
    Excess molecules:             0.8675190741449798 [-]

Component 1 [methane]
    Mol-fraction:                 0.75 [-]
    Fugacity:                     734258.9620174329 [Pa]
    Fugacity coefficient:         0.9790119493565772 [-]
    Compressibility:              0.9714389724700717 [-]
    Excess molecules:             2.6025572224349394 [-]
```

Appreciable adsorption in MOF materials occurs at higher pressure than zeolites, usually in the range up to 10 bar. At these high pressures absolute and excess adsorption are different, and excess adsorption eventually even goes down. This is due to the fact that excess adsorption is relative to what would have been in the free pore volume at these conditions.

To compute the excess adsorption the void fraction of a structure needs to be specified using `HeliumVoidFraction'. RASPA automatically uses an equation of state (default: Peng-Robinson) to compute the fugacities from the pressure and mol-fraction as is done here for a mixture of CO₂ and methane. It also computes the amount of excess molecules from this equation of state. It also computes the amount of excess molecules from this equation of state.

```
Component 0 (CO2)
    Block[  0]  1.363929e+01
    Block[  1]  1.355632e+01
    Block[  2]  1.360550e+01
    Block[  3]  1.362871e+01
    Block[  4]  1.360001e+01
    ---------------------------------------------------------------------------
    Abs. loading average   1.360596e+01 +/-  3.988351e-02 [molecules/cell]
    Abs. loading average   1.360596e+01 +/-  3.988351e-02 [molecules/uc]
    Abs. loading average   2.209264e+00 +/-  6.476074e-03 [mol/kg-framework]
    Abs. loading average   9.720498e+01 +/-  2.849395e-01 [mg/g-framework]

    Block[  0]  1.277177e+01
    Block[  1]  1.268880e+01
    Block[  2]  1.273798e+01
    Block[  3]  1.276119e+01
    Block[  4]  1.273249e+01
    ---------------------------------------------------------------------------
    Excess loading average   1.273844e+01 +/-  3.988351e-02 [molecules/cell]
    Excess loading average   1.273844e+01 +/-  3.988351e-02 [molecules/uc]
    Excess loading average   2.068401e+00 +/-  6.476074e-03 [mol/kg-framework]
    Excess loading average   9.100717e+01 +/-  2.849395e-01 [mg/g-framework]


Component 1 (methane)
    Block[  0]  2.059514e+01
    Block[  1]  2.054633e+01
    Block[  2]  2.054283e+01
    Block[  3]  2.053066e+01
    Block[  4]  2.056720e+01
    ---------------------------------------------------------------------------
    Abs. loading average   2.055643e+01 +/-  3.143894e-02 [molecules/cell]
    Abs. loading average   2.055643e+01 +/-  3.143894e-02 [molecules/uc]
    Abs. loading average   3.337845e+00 +/-  5.104890e-03 [mol/kg-framework]
    Abs. loading average   5.354724e+01 +/-  8.189499e-02 [mg/g-framework]

    Block[  0]  1.799258e+01
    Block[  1]  1.794377e+01
    Block[  2]  1.794027e+01
    Block[  3]  1.792810e+01
    Block[  4]  1.796464e+01
    ---------------------------------------------------------------------------
    Excess loading average   1.795387e+01 +/-  3.143894e-02 [molecules/cell]
    Excess loading average   1.795387e+01 +/-  3.143894e-02 [molecules/uc]
    Excess loading average   2.915255e+00 +/-  5.104890e-03 [mol/kg-framework]
    Excess loading average   4.676786e+01 +/-  8.189499e-02 [mg/g-framework]
```

RASPA also computed the enthalpies of adsorption for the components in the mixture:
```
Component 0 [CO2]
-------------------------------------------------------------------------------
    Block[  0] -1.856496e+03
    Block[  1] -1.854692e+03
    Block[  2] -1.852168e+03
    Block[  3] -1.852296e+03
    Block[  4] -1.853725e+03
    ---------------------------------------------------------------------------
    Enthalpy of adsorption: -1.853875e+03 +/-  2.235383e+00 [K]
                            -1.541398e+01 +/-  1.858601e-02 [kJ/mol]
    Note: need to subtract the ideal-gas energy.

Component 1 [methane]
-------------------------------------------------------------------------------
    Block[  0] -1.404500e+03
    Block[  1] -1.406274e+03
    Block[  2] -1.407237e+03
    Block[  3] -1.409140e+03
    Block[  4] -1.402377e+03
    ---------------------------------------------------------------------------
    Enthalpy of adsorption: -1.405906e+03 +/-  3.213472e+00 [K]
                            -1.168935e+01 +/-  2.671830e-02 [kJ/mol]
    Note: need to subtract the ideal-gas energy.
```

Using Widom insertion we can verify that the measured fugacities are close to the imposed ones:
```
Component 0 [CO2]
    Imposed fugacity:                    2.375876e+05 [Pa]
    Measured fugacity:                   2.308130e+05 +/-  7.343665e+02 [Pa]

Component 1 [methane]
    Imposed fugacity:                    7.342590e+05 [Pa]
    Measured fugacity:                   7.143740e+05 +/-  1.022102e+03 [Pa]
```

#### Monte Carlo: NPT CO₂<a name="Example_nonbasic_2"></a>

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 50000,
  "NumberOfInitializationCycles" : 10000,
  "PrintEvery" : 1000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 250.0,
      "ExternalPressure" : 1e7,
      "ChargeMethod" : "Ewald",
      "VolumeMoveProbability" : 0.05
    }
  ],

  "Components" :
  [
    {
      "Name" : "CO2",
      "MoleculeDefinition" : "ExampleDefinitions",
      "TranslationProbability" : 1.0,
      "RotationProbability" : 1.0,
      "ReinsertionProbability" : 1.0,
      "CreateNumberOfMolecules" : 256
    }
  ]
}
```

```
    Block[  0] 256
    Block[  1] 256
    Block[  2] 256
    Block[  3] 256
    Block[  4] 256
    -----------------------------------------------------------------------
    Density average   2.560000e+02 +/-  0.000000e+00 [molecules]
    Density average   1.278004e-02 +/-  4.137859e-05 [molec/A^3]
    Density average   9.337315e+02 +/-  3.023191e+00 [kg.m⁻³]
```

```
Average pressure tensor:
-------------------------------------------------------------------------------
 9.9463e+01 -7.4457e-01 -2.1266e-01 +/- 3.2169e+00 2.6831e+00 3.3564e+00 [bar]
-7.4457e-01  9.9282e+01 -2.3006e-01 +/- 2.6831e+00 2.6834e+00 3.4452e+00 [bar]
-2.1266e-01 -2.3006e-01  1.0180e+02 +/- 3.3564e+00 3.4452e+00 2.9063e+00 [bar]

    Block[  0]  4.424693e+07
    Block[  1]  4.406916e+07
    Block[  2]  4.396276e+07
    Block[  3]  4.420751e+07
    Block[  4]  4.407319e+07
    ---------------------------------------------------------------------------
    Ideal gas pressure   4.411191e+07 +/-  8.601021e-03 [Pa]
                         4.411191e+02 +/-  1.428234e+00 [bar]


    Block[  0] -3.432308e+07
    Block[  1] -3.413620e+07
    Block[  2] -3.374362e+07
    Block[  3] -3.428751e+07
    Block[  4] -3.397906e+07
    ---------------------------------------------------------------------------
    Excess pressure  -3.409390e+07 +/-  2.961885e+05 [Pa]
                     -3.409390e+02 +/-  2.961885e+00 [bar]


    Block[  0]  9.923855e+06
    Block[  1]  9.932958e+06
    Block[  2]  1.021914e+07
    Block[  3]  9.919992e+06
    Block[  4]  1.009413e+07
    ---------------------------------------------------------------------------
    Pressure average   1.001802e+07 +/-  1.665082e+05 [Pa]
                       1.001802e+02 +/-  1.665082e+00 [bar]
```

```
Monte-Carlo moves statistics
===============================================================================

    Volume change        all:             4197072
    Volume change        total:           4197072
    Volume change        constructed:     4197072
    Volume change        accepted:        2099085
    Volume change        fraction:       0.500131
    Volume change        max-change:     0.027143

Component 0 [CO2]
    Reinsertion (CBMC)   all:            83934567
    Reinsertion (CBMC)   total:          83934567
    Reinsertion (CBMC)   constructed:    81465120
    Reinsertion (CBMC)   accepted:        6662598
    Reinsertion (CBMC)   fraction:       0.079378
    Reinsertion (CBMC)   max-change:     0.000000

    Translation          all:            83927109
    Translation          total:          27977763   27975663   27973683
    Translation          constructed:    27970649   27968514   27966458
    Translation          accepted:       10966567   10967397   10961668
    Translation          fraction:       0.391974   0.392033   0.391856
    Translation          max-change:     0.010000   0.010000   0.010000

    Rotation             all:            83941252
    Rotation             total:          27986886   27978440   27975926
    Rotation             constructed:    27979687   27971583   27968722
    Rotation             accepted:       10969893   10970338   10971701
    Rotation             fraction:       0.391965   0.392100   0.392184
    Rotation             max-change:     0.010000   0.010000   0.010000
```

#### Monte Carlo: NPT propane<a name="Example_nonbasic_3"></a>

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 50000,
  "NumberOfInitializationCycles" : 10000,
  "PrintEvery" : 1000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 250.0,
      "ExternalPressure" : 1e6,
      "ChargeMethod" : "None",
      "VolumeMoveProbability" : 0.05
    }
  ],

  "Components" :
  [
    {
      "Name" : "propane",
      "MoleculeDefinition" : "ExampleDefinitions",
      "TranslationProbability" : 1.0,
      "RotationProbability" : 1.0,
      "ReinsertionProbability" : 1.0,
      "CreateNumberOfMolecules" : 256
    }
  ]
}
```

```
Component 0 (propane)
    Block[  0] 256
    Block[  1] 256
    Block[  2] 256
    Block[  3] 256
    Block[  4] 256
    -----------------------------------------------------------------------
    Density average   2.560000e+02 +/-  0.000000e+00 [molecules]
    Density average   7.761091e-03 +/-  3.026467e-06 [molec/A^3]
    Density average   5.682864e+02 +/-  2.216054e-01 [kg.m⁻³]
```

```
Average pressure tensor:
-------------------------------------------------------------------------------
 9.0757e+00 -8.8946e-02  2.4468e-02 +/- 7.4469e-01 1.8894e+00 9.6681e-01 [bar]
-8.8946e-02  9.6946e+00 -8.2002e-02 +/- 1.8894e+00 2.2443e+00 3.2412e+00 [bar]
 2.4468e-02 -8.2002e-02  1.0125e+01 +/- 9.6681e-01 3.2412e+00 3.8172e+00 [bar]

    Block[  0]  2.678430e+07
    Block[  1]  2.678570e+07
    Block[  2]  2.678683e+07
    Block[  3]  2.680309e+07
    Block[  4]  2.678200e+07
    ---------------------------------------------------------------------------
    Ideal gas pressure   2.678838e+07 +/-  6.290863e-04 [Pa]
                         2.678838e+02 +/-  1.044623e-01 [bar]


    Block[  0] -2.581568e+07
    Block[  1] -2.558996e+07
    Block[  2] -2.607333e+07
    Block[  3] -2.588163e+07
    Block[  4] -2.576545e+07
    ---------------------------------------------------------------------------
    Excess pressure  -2.582521e+07 +/-  2.183639e+05 [Pa]
                     -2.582521e+02 +/-  2.183639e+00 [bar]


    Block[  0]  9.686217e+05
    Block[  1]  1.195745e+06
    Block[  2]  7.134989e+05
    Block[  3]  9.214610e+05
    Block[  4]  1.016545e+06
    ---------------------------------------------------------------------------
    Pressure average   9.631744e+05 +/-  2.159625e+05 [Pa]
                       9.631744e+00 +/-  2.159625e+00 [bar]
```

```
Monte-Carlo moves statistics
===============================================================================

    Volume change        all:             2536350
    Volume change        total:           2536350
    Volume change        constructed:     2536350
    Volume change        accepted:        1266056
    Volume change        fraction:       0.499165
    Volume change        max-change:     0.020876

Component 0 [propane]
    Reinsertion (CBMC)   all:            50693568
    Reinsertion (CBMC)   total:          50693568
    Reinsertion (CBMC)   constructed:    38269580
    Reinsertion (CBMC)   accepted:         301948
    Reinsertion (CBMC)   fraction:       0.005956
    Reinsertion (CBMC)   max-change:     0.000000

    Partial reinsertion (CBMC) all:            50699781
    Partial reinsertion (CBMC) total:          50699781
    Partial reinsertion (CBMC) constructed:    50689654
    Partial reinsertion (CBMC) accepted:       25051762
    Partial reinsertion (CBMC) fraction:       0.494120
    Partial reinsertion (CBMC) max-change:     0.000000

    Widom                all:            50691571
    Widom                total:          50691571
    Widom                constructed:    37774707
    Widom                accepted:              0
    Widom                fraction:       0.000000
    Widom                max-change:     0.000000

    Translation          all:            50686797
    Translation          total:          16893938   16891675   16901184
    Translation          constructed:    16893938   16891675   16901184
    Translation          accepted:        8444975    8442754    8454909
    Translation          fraction:       0.499882   0.499817   0.500255
    Translation          max-change:     0.759807   0.749178   0.739907

    Rotation             all:            50691933
    Rotation             total:          16896569   16894845   16900519
    Rotation             constructed:    16896542   16894825   16900495
    Rotation             accepted:        8444825    8453701    8457582
    Rotation             fraction:       0.499795   0.500372   0.500433
    Rotation             max-change:     0.604939   0.601181   0.601781
```

#### Monte Carlo: Gibbs CO₂<a name="Example_nonbasic_4"></a>

We can compute Vapor-Liquid Equilibrium (VLE) using the Gibbs-ensemble.
```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 200000,
  "NumberOfInitializationCycles" : 50000,
  "PrintEvery" : 1000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 240.0,
      "ChargeMethod" : "Ewald",
      "GibbsVolumeMoveProbability" : 0.01
    },
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 240.0,
      "ChargeMethod" : "Ewald",
      "GibbsVolumeMoveProbability" : 0.01
    }
  ],

  "Components" :
  [
    {
      "Name" : "CO2",
      "MoleculeDefinition" : "ExampleDefinitions",
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "GibbsSwapProbability" : 1.0,
      "CreateNumberOfMolecules" : [256, 256]
    }
  ]
}
```

For each box we obtain a computed density (both the number of molecules and the volume of the boxes change), computing to a vapor phase and a liquid phase. For the gas-phase we obtain
```
Component 0 (CO2)
    Block[  0] 8.581800000000001
    Block[  1] 11.114
    Block[  2] 9.5334
    Block[  3] 10.2088
    Block[  4] 8.6978
    -----------------------------------------------------------------------
    Density average   9.627160e+00 +/-  1.318860e+00 [molecules]
    Density average   4.605977e-04 +/-  5.480961e-05 [molec/A^3]
    Density average   3.365205e+01 +/-  4.004484e+00 [kg.m⁻³]
```
and for the liquid phase we obtain
```
    Block[  0] 503.4182
    Block[  1] 500.886
    Block[  2] 502.4666
    Block[  3] 501.7912
    Block[  4] 503.3022
    -----------------------------------------------------------------------
    Density average   5.023728e+02 +/-  1.318860e+00 [molecules]
    Density average   1.512513e-02 +/-  1.232200e-04 [molec/A^3]
    Density average   1.105068e+03 +/-  9.002663e+00 [kg.m⁻³]
```
Uisng the NIST chemical database, we can compare to the experimental values of 33.295 and 1088.9 kg/m³. The experimental pressure from NIST is 12.825 bar. From the simulations we obtain \f$13.08 \pm 1.3\f$ and \f$12.27 \pm 9.5\f$ bar for the vapor and the liquid phases, respectively. Note that in equilibrium, the temperatures, chemical potentials, and pressures of the phases is equal, but the vapor-phase data has significantly lower error bars.

#### Molecular Dynamics: benzene in IRMOF-1<a name="Example_nonbasic_5"></a>

```json
{
  "SimulationType" : "MolecularDynamics",
  "NumberOfCycles" : 500000,
  "NumberOfInitializationCycles" : 1000,
  "NumberOfEquilibrationCycles" : 100000,
  "PrintEvery" : 5000,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "IRMOF-1",
      "NumberOfUnitCells" : [1, 1, 1],
      "ChargeMethod" : "Ewald",
      "CutOff" : 12.0,
      "ExternalTemperature" : 300.0,
      "Ensemble" : "NVT",
      "ComputeMSD" : true,
      "SampleMSDEvery" : 10,
      "WriteMSDEvery" : 5000
    }
  ],

  "Components" : [
    {
      "Name" : "benzene",
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "CreateNumberOfMolecules" : 16
    }
  ]
}
```

#### Monte Carlo: adsorption CO₂ in LTA-4A sodium<a name="Example_nonbasic_6"></a>

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 25000,
  "NumberOfInitializationCycles" : 10000,
  "PrintEvery" : 1000,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "LTA4A",
      "NumberOfUnitCells" : [1, 1, 1],
      "ExternalTemperature" : 298.0,
      "ExternalPressure" : 1.0e4,
      "ChargeMethod" : "Ewald"
    }
  ],

  "Components" : [
    {
      "Name" : "sodium",
      "Type" : "Cation",
      "TranslationProbability" : 0.5,
      "RandomTranslationProbability" : 0.5,
      "CreateNuMberofmolecules" : 96
    },
    {
      "Name" : "CO2",
      "FugacityCoefficient" : 1.0,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "SwapProbability" : 1.0,
      "blockingPockets" : [
           [0.0,       0.0,        0.0,       4.0],
           [0.5,       0.0,        0.0,       4.0],
           [0.0,       0.5,        0.0,       4.0],
           [0.5,       0.5,        0.0,       4.0],
           [0.0,       0.0,        0.5,       4.0],
           [0.5,       0.0,        0.5,       4.0],
           [0.0,       0.5,        0.5,       4.0],
           [0.5,       0.5,        0.5,       4.0],
           [0.25,      0.0,        0.0,       1.0],
           [0.75,      0.0,        0.0,       1.0],
           [0.25,      0.5,        0.0,       1.0],
           [0.75,      0.5,        0.0,       1.0],
           [0.0,       0.25,       0.0,       1.0],
           [0.0,       0.75,       0.0,       1.0],
           [0.5,       0.25,       0.0,       1.0],
           [0.5,       0.75,       0.0,       1.0],
           [0.0,       0.0,        0.25,      1.0],
           [0.0,       0.0,        0.75,      1.0],
           [0.0,       0.5,        0.25,      1.0],
           [0.0,       0.5,        0.75,      1.0],
           [0.25,      0.0,        0.5,       1.0],
           [0.75,      0.0,        0.5,       1.0],
           [0.25,      0.5,        0.5,       1.0],
           [0.75,      0.5,        0.5,       1.0],
           [0.0,       0.25,       0.5,       1.0],
           [0.0,       0.75,       0.5,       1.0],
           [0.5,       0.25,       0.5,       1.0],
           [0.5,       0.75,       0.5,       1.0],
           [0.5,       0.0,        0.25,      1.0],
           [0.5,       0.0,        0.75,      1.0],
           [0.5,       0.5,        0.25,      1.0],
           [0.5,       0.5,        0.75,      1.0]
         ],
      "CreateNuMberofmolecules" : 0
    }
  ]
}
```

#### Molecular Dynamics: diffusion of CO₂ in LTA-4A sodium<a name="Example_nonbasic_7"></a>

```json
{
  "SimulationType" : "MolecularDynamics",
  "NumberOfCycles" : 250000,
  "NumberOfInitializationCycles" : 5000,
  "NumberOfEquilibrationCycles" : 10000,
  "PrintEvery" : 5000,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "LTA4A",
      "NumberOfUnitCells" : [1, 1, 1],
      "ExternalTemperature" : 600.0,
      "Ensemble" : "NVT",
      "ChargeMethod" : "Ewald",
      "ComputeMSD" : true,
      "SampleMSDEvery" : 10,
      "WriteMSDEvery" : 5000
    }
  ],

  "Components" : [
    {
      "Name" : "sodium",
      "Type" : "Cation",
      "TranslationProbability" : 0.5,
      "RandomTranslationProbability" : 0.5,
      "CreateNuMberofmolecules" : 96
    },
    {
      "Name" : "CO2",
      "FugacityCoefficient" : 1.0,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "blockingPockets" : [
           [0.0,       0.0,        0.0,       4.0],
           [0.5,       0.0,        0.0,       4.0],
           [0.0,       0.5,        0.0,       4.0],
           [0.5,       0.5,        0.0,       4.0],
           [0.0,       0.0,        0.5,       4.0],
           [0.5,       0.0,        0.5,       4.0],
           [0.0,       0.5,        0.5,       4.0],
           [0.5,       0.5,        0.5,       4.0],
           [0.25,      0.0,        0.0,       1.0],
           [0.75,      0.0,        0.0,       1.0],
           [0.25,      0.5,        0.0,       1.0],
           [0.75,      0.5,        0.0,       1.0],
           [0.0,       0.25,       0.0,       1.0],
           [0.0,       0.75,       0.0,       1.0],
           [0.5,       0.25,       0.0,       1.0],
           [0.5,       0.75,       0.0,       1.0],
           [0.0,       0.0,        0.25,      1.0],
           [0.0,       0.0,        0.75,      1.0],
           [0.0,       0.5,        0.25,      1.0],
           [0.0,       0.5,        0.75,      1.0],
           [0.25,      0.0,        0.5,       1.0],
           [0.75,      0.0,        0.5,       1.0],
           [0.25,      0.5,        0.5,       1.0],
           [0.75,      0.5,        0.5,       1.0],
           [0.0,       0.25,       0.5,       1.0],
           [0.0,       0.75,       0.5,       1.0],
           [0.5,       0.25,       0.5,       1.0],
           [0.5,       0.75,       0.5,       1.0],
           [0.5,       0.0,        0.25,      1.0],
           [0.5,       0.0,        0.75,      1.0],
           [0.5,       0.5,        0.25,      1.0],
           [0.5,       0.5,        0.75,      1.0]
         ],
      "CreateNuMberofmolecules" : 64
    }
  ]
}
```

#### Monte Carlo: adsorption of C6-isomers in MFI<a name="Example_nonbasic_8"></a>


```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 500000,
  "NumberOfInitializationCycles" : 100000,
  "PrintEvery" : 5000,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "MFI_SI",
      "NumberOfUnitCells" : [2, 2, 2],
      "ExternalTemperature" : 433.0,
      "ExternalPressure" : 1.0e5,
      "ChargeMethod" : "None"
    }
  ],

  "Components" : [
    {
      "Name" : "hexane",
      "MolFraction" : 0.2,
      "FugacityCoefficient" : 1.0,
      "IdealGasRosenbluthWeight" : 8.103901e-03,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "PartialReinsertionProbability" : 0.5,
      "SwapProbability" : 1.0,
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "2-methylpentane",
      "MolFraction" : 0.2,
      "FugacityCoefficient" : 1.0,
      "IdealGasRosenbluthWeight" : 4.704858e-02,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "PartialReinsertionProbability" : 0.5,
      "SwapProbability" : 1.0,
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "3-methylpentane",
      "MolFraction" : 0.2,
      "FugacityCoefficient" : 1.0,
      "IdealGasRosenbluthWeight" : 5.353616e-02,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "PartialReinsertionProbability" : 0.5,
      "SwapProbability" : 1.0,
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "22-methylbutane",
      "MolFraction" : 0.2,
      "FugacityCoefficient" : 1.0,
      "IdealGasRosenbluthWeight" : 2.265790e-01,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "PartialReinsertionProbability" : 0.5,
      "SwapProbability" : 1.0,
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    },
    {
      "Name" : "23-methylbutane",
      "MolFraction" : 0.2,
      "FugacityCoefficient" : 1.0,
      "IdealGasRosenbluthWeight" : 8.732375e-02,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "PartialReinsertionProbability" : 0.5,
      "SwapProbability" : 1.0,
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

```
Component 0 (hexane)
    Block[  0]  1.938018e+01
    Block[  1]  2.000694e+01
    Block[  2]  2.042880e+01
    Block[  3]  2.147259e+01
    Block[  4]  2.069375e+01
    ---------------------------------------------------------------------------
    Abs. loading average   2.039645e+01 +/-  9.680310e-01 [molecules/cell]
    Abs. loading average   2.549557e+00 +/-  1.210039e-01 [molecules/uc]
    Abs. loading average   4.420603e-01 +/-  2.098052e-02 [mol/kg-framework]
    Abs. loading average   3.812981e+01 +/-  1.809670e+00 [mg/g-framework]


Component 1 (2-methylpentane)
    Block[  0]  8.119860e+00
    Block[  1]  7.462680e+00
    Block[  2]  7.671530e+00
    Block[  3]  7.182510e+00
    Block[  4]  7.427390e+00
    ---------------------------------------------------------------------------
    Abs. loading average   7.572794e+00 +/-  4.365559e-01 [molecules/cell]
    Abs. loading average   9.465993e-01 +/-  5.456949e-02 [molecules/uc]
    Abs. loading average   1.641281e-01 +/-  9.461648e-03 [mol/kg-framework]
    Abs. loading average   1.415683e+01 +/-  8.161122e-01 [mg/g-framework]


Component 2 (3-methylpentane)
    Block[  0]  3.988190e+00
    Block[  1]  4.197100e+00
    Block[  2]  4.042350e+00
    Block[  3]  4.091960e+00
    Block[  4]  4.177710e+00
    ---------------------------------------------------------------------------
    Abs. loading average   4.099462e+00 +/-  1.099160e-01 [molecules/cell]
    Abs. loading average   5.124327e-01 +/-  1.373950e-02 [molecules/uc]
    Abs. loading average   8.884926e-02 +/-  2.382253e-03 [mol/kg-framework]
    Abs. loading average   7.663671e+00 +/-  2.054807e-01 [mg/g-framework]


Component 3 (22-methylbutane)
    Block[  0]  1.312820e+00
    Block[  1]  1.435550e+00
    Block[  2]  1.322000e+00
    Block[  3]  1.228440e+00
    Block[  4]  1.239790e+00
    ---------------------------------------------------------------------------
    Abs. loading average   1.307720e+00 +/-  1.028811e-01 [molecules/cell]
    Abs. loading average   1.634650e-01 +/-  1.286014e-02 [molecules/uc]
    Abs. loading average   2.834273e-02 +/-  2.229783e-03 [mol/kg-framework]
    Abs. loading average   2.444695e+00 +/-  1.923294e-01 [mg/g-framework]


Component 4 (23-methylbutane)
    Block[  0]  1.595000e+00
    Block[  1]  1.534460e+00
    Block[  2]  1.393040e+00
    Block[  3]  1.287950e+00
    Block[  4]  1.423150e+00
    ---------------------------------------------------------------------------
    Abs. loading average   1.446720e+00 +/-  1.499172e-01 [molecules/cell]
    Abs. loading average   1.808400e-01 +/-  1.873965e-02 [molecules/uc]
    Abs. loading average   3.135533e-02 +/-  3.249214e-03 [mol/kg-framework]
    Abs. loading average   2.704547e+00 +/-  2.802602e-01 [mg/g-framework]
```

```
Component 0 [hexane]
-------------------------------------------------------------------------------
    Block[  0] -6.442828e+03
    Block[  1] -6.490959e+03
    Block[  2] -6.486642e+03
    Block[  3] -6.398685e+03
    Block[  4] -6.427075e+03
    ---------------------------------------------------------------------------
    Enthalpy of adsorption: -6.449238e+03 +/-  4.898491e+01 [K]
                            -5.362196e+01 +/-  4.072833e-01 [kJ/mol]
    Note: need to subtract the ideal-gas energy.

Component 1 [2-methylpentane]
-------------------------------------------------------------------------------
    Block[  0] -6.345850e+03
    Block[  1] -6.512363e+03
    Block[  2] -6.322441e+03
    Block[  3] -6.341487e+03
    Block[  4] -6.339311e+03
    ---------------------------------------------------------------------------
    Enthalpy of adsorption: -6.372290e+03 +/-  9.783339e+01 [K]
                            -5.298219e+01 +/-  8.134323e-01 [kJ/mol]
    Note: need to subtract the ideal-gas energy.

Component 2 [3-methylpentane]
-------------------------------------------------------------------------------
    Block[  0] -5.846014e+03
    Block[  1] -6.077024e+03
    Block[  2] -5.809159e+03
    Block[  3] -5.762253e+03
    Block[  4] -5.786489e+03
    ---------------------------------------------------------------------------
    Enthalpy of adsorption: -5.856188e+03 +/-  1.579706e+02 [K]
                            -4.869107e+01 +/-  1.313441e+00 [kJ/mol]
    Note: need to subtract the ideal-gas energy.

Component 3 [22-methylbutane]
-------------------------------------------------------------------------------
    Block[  0] -5.300764e+03
    Block[  1] -5.438470e+03
    Block[  2] -5.733780e+03
    Block[  3] -5.454107e+03
    Block[  4] -5.464367e+03
    ---------------------------------------------------------------------------
    Enthalpy of adsorption: -5.478298e+03 +/-  1.954619e+02 [K]
                            -4.554911e+01 +/-  1.625161e+00 [kJ/mol]
    Note: need to subtract the ideal-gas energy.

Component 4 [23-methylbutane]
-------------------------------------------------------------------------------
    Block[  0] -5.318777e+03
    Block[  1] -5.611167e+03
    Block[  2] -5.072493e+03
    Block[  3] -5.249716e+03
    Block[  4] -5.261793e+03
    ---------------------------------------------------------------------------
    Enthalpy of adsorption: -5.302789e+03 +/-  2.427102e+02 [K]
                            -4.408985e+01 +/-  2.018005e+00 [kJ/mol]
    Note: need to subtract the ideal-gas energy.
```
