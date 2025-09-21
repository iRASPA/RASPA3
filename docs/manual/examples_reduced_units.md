# Examples Basic
\page examples_reduced_units Examples Reduced Units


## Table of Contents
1. [Monte Carlo: particle in box (NVT)](#Example_reduced_units_1)
2. [Monte Carlo: particle in box (CFCMC)](#Example_reduced_units_2)
3. [Monte Carlo: Gibbs Vapor-Liquid Equilibrium](#Example_reduced_units_3)
4. [Monte Carlo: Gibbs Vapor-Liquid Equilibrium (CFCMC)](#Example_reduced_units_4)


#### Monte Carlo: particles in box (NVT) <a name="Example_reduced_units_1"></a>

A lennard-Jones fluid is simple model of a fluid interacting via the reduced Lennard-Jones potential
$$U(r) = 4\left[\frac{1}{r^{12}} - \frac{1}{r^6}\right]$$

Here, we use a cutoff of 4.0 \f$\sigma\f$ and tail-corrections and compute via Widom insertion the chemical potential.

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 1000000,
  "NumberOfInitializationCycles" : 100000,
  "PrintEvery" : 5000,
  "Units" : "Reduced",

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [9.49914251593, 9.49914251593, 9.49914251593],
      "ExternalTemperature" : 1.2,
      "ChargeMethod" : "None",
      "OutputPDBMovie" : false,
      "SampleMovieEvery" : 10
    }
  ],

  "Components" :
  [
    {
      "Name" : "particle",
      "MoleculeDefinition" : "ExampleDefinitions",
      "TranslationProbability" : 1.0,
      "ReinsertionProbability" : 0.25,
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : 600
    }
  ]
}
```

The average of the total energy:
```
Total energy
-------------------------------------------------------------------------------
    Block[  0] -2.838221e+03
    Block[  1] -2.838781e+03
    Block[  2] -2.838190e+03
    Block[  3] -2.838883e+03
    Block[  4] -2.838060e+03
    ---------------------------------------------------------------------------
    Average  -2.838427e+03 +/-  4.670014e-01 [ε]
```

The average pressure is computed as:
```
Average pressure tensor:
-------------------------------------------------------------------------------
 7.0379e-01  3.4238e-04 -1.6270e-03 +/- 4.7707e-03 2.4653e-03 1.6942e-03 [ε.σ⁻³]
 3.4238e-04  7.0388e-01 -5.1091e-05 +/- 2.4653e-03 4.6656e-03 1.3557e-03 [ε.σ⁻³]
-1.6270e-03 -5.1091e-05  6.9996e-01 +/- 1.6942e-03 1.3557e-03 4.7725e-03 [ε.σ⁻³]

    Block[  0]  8.400000e-01
    Block[  1]  8.400000e-01
    Block[  2]  8.400000e-01
    Block[  3]  8.400000e-01
    Block[  4]  8.400000e-01
    ---------------------------------------------------------------------------
    Ideal gas pressure   8.400000e-01 +/-  0.000000e+00 [ε.σ⁻³]


    Block[  0] -1.351949e-01
    Block[  1] -1.408517e-01
    Block[  2] -1.367764e-01
    Block[  3] -1.402155e-01
    Block[  4] -1.342423e-01
    ---------------------------------------------------------------------------
    Excess pressure  -1.374562e-01 +/-  3.674814e-03 [ε.σ⁻³]


    Block[  0]  7.048051e-01
    Block[  1]  6.991483e-01
    Block[  2]  7.032236e-01
    Block[  3]  6.997845e-01
    Block[  4]  7.057577e-01
    ---------------------------------------------------------------------------
    Pressure average   7.025438e-01 +/-  3.674814e-03 [ε.σ⁻³]
```

Using Widom insertion, we find the chemical potential is given by:
```
    Widom insertion Rosenbluth weight statistics:
    ---------------------------------------------------------------------------
        Block[  0]  6.489275e+00
        Block[  1]  6.457452e+00
        Block[  2]  6.460533e+00
        Block[  3]  6.487067e+00
        Block[  4]  6.494927e+00
    ---------------------------------------------------------------------------
    Average Rosenbluth weight:    6.477852e+00 +/-  2.170840e-02 [-]


    Widom insertion chemical potential  statistics:
    ---------------------------------------------------------------------------
        Block[  0] -1.8701508415534316
        Block[  1] -1.8652347611211533
        Block[  2] -1.8657118770613563
        Block[  3] -1.8698104999091023
        Block[  4] -1.8710213520189722
    ---------------------------------------------------------------------------
    Beta * Excess chemical potential:          -1.868389e+00 +/-  3.352543e-03 [-]
    Beta * Ideal chemical potential:           -3.566749e-01 +/-  1.832108e-17 [-]
    Beta * Total chemical potential:           -2.225064e+00 +/-  3.352543e-03 [-]
    Beta * Imposed chemical potential:  -inf [-]
    ---------------------------------------------------------------------------
    Excess chemical potential:          -2.242067e+00 +/-  4.023052e-03 [ε]
    Ideal chemical potential:           -4.280099e-01 +/-  2.198529e-17 [ε]
    Total chemical potential:           -2.670077e+00 +/-  4.023052e-03 [ε]
    Imposed chemical potential:  -inf [ε]
    ---------------------------------------------------------------------------
    Imposed fugacity:             0.000000e+00 [ε.σ⁻³]
    Measured fugacity:            1.296726e-01 +/-  4.349127e-04 [ε.σ⁻³]
```


#### Monte Carlo: particles in a box (CFCMC)<a name="Example_reduced_units_2"></a>

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 1000000,
  "NumberOfInitializationCycles" : 100000,
  "NumberOfEquilibrationCycles" : 200000,
  "PrintEvery" : 5000,
  "Units" : "Reduced",

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [9.49914251593, 9.49914251593, 9.49914251593],
      "ExternalTemperature" : 1.2,
      "ChargeMethod" : "None"
    }
  ],

  "Components" :
  [
    {
      "Name" : "particle",
      "MoleculeDefinition" : "ExampleDefinitions",
      "ThermoDynamicIntegration" : true,
      "FugacityCoefficient" : 1.0,
      "TranslationProbability" : 1.0,
      "ReinsertionProbability" : 0.25,
      "WidomProbability" : 0.5,
      "CFCMC_CBMC_WidomProbability" : 0.5,
      "CreateNumberOfMolecules" : 600
    }
  ]
}
```

```
    Lambda histogram and bias:
    ---------------------------------------------------------------------------
     0-0.000000 (lambda) P:  1.30774e+00 +/- 7.14095e-02 bias:  2.21994e+00 [-]
     1-0.025000 (lambda) P:  1.30331e+00 +/- 6.80667e-02 bias:  2.21503e+00 [-]
     2-0.050000 (lambda) P:  1.28486e+00 +/- 4.81727e-02 bias:  2.21321e+00 [-]
     3-0.075000 (lambda) P:  1.29142e+00 +/- 6.08460e-02 bias:  2.24468e+00 [-]
     4-0.100000 (lambda) P:  1.25968e+00 +/- 5.46987e-02 bias:  2.33087e+00 [-]
     5-0.125000 (lambda) P:  1.32110e+00 +/- 7.34593e-02 bias:  2.67301e+00 [-]
     6-0.150000 (lambda) P:  1.24415e+00 +/- 4.60494e-02 bias:  3.20426e+00 [-]
     7-0.175000 (lambda) P:  1.14816e+00 +/- 6.80581e-02 bias:  3.71065e+00 [-]
     8-0.200000 (lambda) P:  1.04607e+00 +/- 5.85775e-02 bias:  3.78275e+00 [-]
     9-0.225000 (lambda) P:  9.88961e-01 +/- 5.01299e-02 bias:  3.56545e+00 [-]
    10-0.250000 (lambda) P:  9.60097e-01 +/- 5.06641e-03 bias:  3.18626e+00 [-]
    11-0.275000 (lambda) P:  9.93020e-01 +/- 1.36496e-02 bias:  2.73752e+00 [-]
    12-0.300000 (lambda) P:  9.58129e-01 +/- 1.76886e-02 bias:  2.19360e+00 [-]
    13-0.325000 (lambda) P:  9.51282e-01 +/- 1.46261e-02 bias:  1.66696e+00 [-]
    14-0.350000 (lambda) P:  9.38367e-01 +/- 2.47256e-02 bias:  1.18065e+00 [-]
    15-0.375000 (lambda) P:  9.27297e-01 +/- 3.68598e-02 bias:  7.80146e-01 [-]
    16-0.400000 (lambda) P:  9.21885e-01 +/- 2.43490e-02 bias:  4.63859e-01 [-]
    17-0.425000 (lambda) P:  9.25534e-01 +/- 1.99206e-02 bias:  2.43098e-01 [-]
    18-0.450000 (lambda) P:  9.13562e-01 +/- 1.29640e-02 bias:  9.40928e-02 [-]
    19-0.475000 (lambda) P:  9.21106e-01 +/- 3.09024e-02 bias:  3.13973e-02 [-]
    20-0.500000 (lambda) P:  9.27051e-01 +/- 2.34937e-02 bias:  1.18339e-02 [-]
    21-0.525000 (lambda) P:  9.25903e-01 +/- 2.87938e-02 bias:  2.85258e-02 [-]
    22-0.550000 (lambda) P:  9.16145e-01 +/- 8.62114e-03 bias:  1.48056e-02 [-]
    23-0.575000 (lambda) P:  9.35169e-01 +/- 2.43952e-02 bias:  9.52742e-03 [-]
    24-0.600000 (lambda) P:  9.28978e-01 +/- 2.50483e-02 bias:  1.15441e-02 [-]
    25-0.625000 (lambda) P:  9.19015e-01 +/- 1.88459e-02 bias:  0.00000e+00 [-]
    26-0.650000 (lambda) P:  9.14177e-01 +/- 1.72043e-02 bias:  2.82292e-03 [-]
    27-0.675000 (lambda) P:  9.02533e-01 +/- 1.27092e-02 bias:  6.07650e-04 [-]
    28-0.700000 (lambda) P:  9.21598e-01 +/- 1.68121e-02 bias:  8.97147e-03 [-]
    29-0.725000 (lambda) P:  9.35292e-01 +/- 2.23957e-02 bias:  1.08080e-02 [-]
    30-0.750000 (lambda) P:  9.13931e-01 +/- 2.73706e-02 bias:  3.08510e-03 [-]
    31-0.775000 (lambda) P:  9.32832e-01 +/- 3.29775e-02 bias:  1.53255e-02 [-]
    32-0.800000 (lambda) P:  9.26190e-01 +/- 1.82922e-02 bias:  1.63149e-02 [-]
    33-0.825000 (lambda) P:  9.23935e-01 +/- 2.22604e-02 bias:  1.12784e-02 [-]
    34-0.850000 (lambda) P:  9.17867e-01 +/- 1.99313e-02 bias:  9.48316e-03 [-]
    35-0.875000 (lambda) P:  9.21516e-01 +/- 2.50365e-02 bias:  1.29228e-02 [-]
    36-0.900000 (lambda) P:  9.27461e-01 +/- 2.29343e-02 bias:  7.83733e-03 [-]
    37-0.925000 (lambda) P:  9.24714e-01 +/- 1.77014e-02 bias:  1.15458e-02 [-]
    38-0.950000 (lambda) P:  9.25534e-01 +/- 3.27156e-02 bias:  6.93352e-03 [-]
    39-0.975000 (lambda) P:  9.24181e-01 +/- 3.17705e-02 bias:  1.63701e-02 [-]
    40-1.000000 (lambda) P:  9.30249e-01 +/- 1.97409e-02 bias:  3.75212e-03 [-]


    Lambda statistics:
    ---------------------------------------------------------------------------
     0-0.000000 (lambda) Free energy: 1.016884e-02 +/- 5.468281e-02 [-]
     1-0.025000 (lambda) Free energy: 1.356059e-02 +/- 5.115551e-02 [-]
     2-0.050000 (lambda) Free energy: 2.781803e-02 +/- 3.695043e-02 [-]
     3-0.075000 (lambda) Free energy: 2.272540e-02 +/- 4.665128e-02 [-]
     4-0.100000 (lambda) Free energy: 4.760534e-02 +/- 4.287716e-02 [-]
     5-0.125000 (lambda) Free energy: 0.000000e+00 +/- 5.537236e-02 [-]
     6-0.150000 (lambda) Free energy: 6.001769e-02 +/- 3.677614e-02 [-]
     7-0.175000 (lambda) Free energy: 1.403021e-01 +/- 5.927641e-02 [-]
     8-0.200000 (lambda) Free energy: 2.334221e-01 +/- 5.567939e-02 [-]
     9-0.225000 (lambda) Free energy: 2.895666e-01 +/- 5.020972e-02 [-]
    10-0.250000 (lambda) Free energy: 3.191872e-01 +/- 5.282070e-03 [-]
    11-0.275000 (lambda) Free energy: 2.854707e-01 +/- 1.375505e-02 [-]
    12-0.300000 (lambda) Free energy: 3.212391e-01 +/- 1.839641e-02 [-]
    13-0.325000 (lambda) Free energy: 3.284110e-01 +/- 1.536563e-02 [-]
    14-0.350000 (lambda) Free energy: 3.420804e-01 +/- 2.648885e-02 [-]
    15-0.375000 (lambda) Free energy: 3.539476e-01 +/- 3.963937e-02 [-]
    16-0.400000 (lambda) Free energy: 3.598010e-01 +/- 2.675731e-02 [-]
    17-0.425000 (lambda) Free energy: 3.558506e-01 +/- 2.152976e-02 [-]
    18-0.450000 (lambda) Free energy: 3.688703e-01 +/- 1.418260e-02 [-]
    19-0.475000 (lambda) Free energy: 3.606464e-01 +/- 3.405690e-02 [-]
    20-0.500000 (lambda) Free energy: 3.542129e-01 +/- 2.543852e-02 [-]
    21-0.525000 (lambda) Free energy: 3.554520e-01 +/- 3.124870e-02 [-]
    22-0.550000 (lambda) Free energy: 3.660469e-01 +/- 9.385207e-03 [-]
    23-0.575000 (lambda) Free energy: 3.454943e-01 +/- 2.611788e-02 [-]
    24-0.600000 (lambda) Free energy: 3.521365e-01 +/- 2.709814e-02 [-]
    25-0.625000 (lambda) Free energy: 3.629191e-01 +/- 2.065218e-02 [-]
    26-0.650000 (lambda) Free energy: 3.681973e-01 +/- 1.884513e-02 [-]
    27-0.675000 (lambda) Free energy: 3.810163e-01 +/- 1.412701e-02 [-]
    28-0.700000 (lambda) Free energy: 3.601124e-01 +/- 1.822148e-02 [-]
    29-0.725000 (lambda) Free energy: 3.453627e-01 +/- 2.403513e-02 [-]
    30-0.750000 (lambda) Free energy: 3.684664e-01 +/- 2.978111e-02 [-]
    31-0.775000 (lambda) Free energy: 3.479964e-01 +/- 3.583179e-02 [-]
    32-0.800000 (lambda) Free energy: 3.551421e-01 +/- 1.978602e-02 [-]
    33-0.825000 (lambda) Free energy: 3.575798e-01 +/- 2.398551e-02 [-]
    34-0.850000 (lambda) Free energy: 3.641690e-01 +/- 2.160428e-02 [-]
    35-0.875000 (lambda) Free energy: 3.602014e-01 +/- 2.741583e-02 [-]
    36-0.900000 (lambda) Free energy: 3.537708e-01 +/- 2.475013e-02 [-]
    37-0.925000 (lambda) Free energy: 3.567370e-01 +/- 1.905506e-02 [-]
    38-0.950000 (lambda) Free energy: 3.558506e-01 +/- 3.577566e-02 [-]
    39-0.975000 (lambda) Free energy: 3.573136e-01 +/- 3.427250e-02 [-]
    40-1.000000 (lambda) Free energy: 3.507692e-01 +/- 2.122766e-02 [-]
    ---------------------------------------------------------------------------
    Excess chemical potential: (ln(P(lambda=1))-ln(P(lambda=0)))
        Block[  0] -1.9425521171520197
        Block[  1] -1.9020588353830536
        Block[  2] -1.7945046752415006
        Block[  3] -1.864956751792043
        Block[  4] -1.877161937515691
    ---------------------------------------------------------------------------
    Beta * Excess chemical potential:    -1.875589e+00 +/-  6.764247e-02 [-]
    Beta * Ideal gas chemical potential: -3.566749e-01 +/-  4.313206e-14 [-]
    Beta * Total chemical potential:     -2.232264e+00 +/-  6.764247e-02 [-]
    Beta * Imposed chemical potential:   -inf [-]
    ---------------------------------------------------------------------------
    Excess chemical potential:    -2.250707e+00 +/-  8.117096e-02 [ε]
    Ideal gas chemical potential: -4.280099e-01 +/-  5.175848e-14 [ε]
    Total chemical potential:     -2.678716e+00 +/-  8.117096e-02 [ε]
    Imposed chemical potential:   -inf [ε]
    ---------------------------------------------------------------------------
    Imposed fugacity:              0.000000e+00 [ε.σ⁻³]
    Measured fugacity:             1.287423e-01 +/-  8.800834e-03 [ε.σ⁻³]
```

#### Monte Carlo: Gibbs Vapor-Liquid Equilibrium<a name="Example_reduced_units_3"></a>

Shifted potential at a  cutoff of \f$r_c=2.5\sigma\f$.

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 1000000,
  "NumberOfInitializationCycles" : 50000,
  "NumberOfEquilibrationCycles" : 200000,
  "PrintEvery" : 5000,
  "Units" : "Reduced",

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [12.5, 12.5, 12.5],
      "ExternalTemperature" : 0.6,
      "ChargeMethod" : "None",
      "GibbsVolumeMoveProbability" : 0.01
    },
    {
      "Type" : "Box",
      "BoxLengths" : [12.5, 12.5, 12.5],
      "ExternalTemperature" : 0.6,
      "ChargeMethod" : "None",
      "GibbsVolumeMoveProbability" : 0.01
    }
  ],

  "Components" :
  [
    {
      "Name" : "particle",
      "MoleculeDefinition" : "ExampleDefinitions",
      "ThermoDynamicIntegration" : true,
      "TranslationProbability" : 1.0,
      "ReinsertionProbability" : 0.25,
      "WidomProbability" : 0.25,
      "GibbsSwapProbability" : 1.0,
      "CreateNumberOfMolecules" : [128, 128]
    }
  ]
}
```

| Property           | Units | Vapor                                 | Liquid                                           |
|--------------------|-------|---------------------------------------|--------------------------------------------------|
| Density \f$\left\langle\rho\right\rangle\f$ | \f$N\cdot\sigma^{-3}\f$  | \f$0.001928 \pm 0.000059\f$ | \f$0.838342 \pm 0.000317\f$   |
| Molecules \f$\left\langle N\right\rangle\f$ | -  | \f$06.96\pm 0.21\f$ | \f$249.04 \pm 0.21\f$
| Energy \f$\left\langle U\right\rangle\f$ | \f$\epsilon\f$  | \f$-0.17 \pm  0.012\f$ | \f$-1321.05 \pm 0.71\f$ |
| Pressure \f$\left\langle p\right\rangle\f$ | \f$\epsilon\sigma^{-3}\f$   | \f$0.001132 \pm 0.000035\f$ | \f$0.001369 \pm 0.000728\f$ 
| Excess chemical potential \f$\left\langle \mu^{ex}\right\rangle\f$ | \f$\epsilon\f$  |  \f$-0.02528\pm0.00079\f$ | \f$-3.67556\pm 0.0261\f$ |
| Ideal chemical potential \f$\left\langle \mu^{id}\right\rangle\f$ | \f$\epsilon\f$  | \f$-3.7509\pm 0.0185\f$ | \f$-0.105798\pm  0.000227\f$ |
| Chemical potential \f$\left\langle \mu\right\rangle\f$ | \f$\epsilon\f$  | \f$-3.7762\pm 0.0177\f$ | \f$-3.7814\pm  0.0261\f$ |

Liquid
```
    Gibbs swap (CBMC)    all:           509958010
    Gibbs swap (CBMC)    total:        1019276656
    Gibbs swap (CBMC)    constructed:   659678014
    Gibbs swap (CBMC)    accepted:        1223173
    Gibbs swap (CBMC)    fraction:       0.001200
    Gibbs swap (CBMC)    max-change:     0.000000
```

Vapor
```
    Gibbs swap (CBMC)    all:           509911407
    Gibbs swap (CBMC)    total:        1019276656
    Gibbs swap (CBMC)    constructed:   659678014
    Gibbs swap (CBMC)    accepted:        1223173
    Gibbs swap (CBMC)    fraction:       0.001200
    Gibbs swap (CBMC)    max-change:     0.000000
```


#### Monte Carlo: Gibbs Vapor-Liquid Equilibrium (CFCMC)<a name="Example_reduced_units_4"></a>

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 1000000,
  "NumberOfInitializationCycles" : 50000,
  "NumberOfEquilibrationCycles" : 1000000,
  "RescaleWangLandauEvery" : 50000,
  "PrintEvery" : 5000,
  "Units" : "Reduced",

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [12.5, 12.5, 12.5],
      "ExternalTemperature" : 0.6,
      "ChargeMethod" : "None",
      "GibbsVolumeMoveProbability" : 0.01
    },
    {
      "Type" : "Box",
      "BoxLengths" : [12.5, 12.5, 12.5],
      "ExternalTemperature" : 0.6,
      "ChargeMethod" : "None",
      "GibbsVolumeMoveProbability" : 0.01
    }
  ],

  "Components" :
  [
    {
      "Name" : "particle",
      "MoleculeDefinition" : "ExampleDefinitions",
      "ThermoDynamicIntegration" : true,
      "TranslationProbability" : 1.0,
      "ReinsertionProbability" : 0.25,
      "WidomProbability" : 0.5,
      "Gibbs_CFCMC_SwapProbability" : 1.0,
      "CreateNumberOfMolecules" : [128, 128]
    }
  ]
}
```

| Property           | Units | Vapor                                  | Liquid                                           |
|--------------------|-------|----------------------------------------|--------------------------------------------------|
| Density \f$\left\langle\rho\right\rangle\f$ | \f$N\cdot\sigma^{-3}\f$ | \f$0.001946 \pm  0.000058\f$   | \f$ 0.836678 \pm 0.000218\f$ |
| Molecules \f$\left\langle N\right\rangle\f$ | - |  \f$7.02\pm 0.21\f$ |  \f$248.98 \pm 0.21\f$ |
| Energy \f$\left\langle U\right\rangle\f$ | \f$\epsilon\f$ | \f$-0.19 \pm 0.01\f$  | \f$-1322.98 \pm 1.48\f$|
| Pressure \f$\left\langle p\right\rangle\f$ | \f$\epsilon\sigma^{-3}\f$ | \f$0.001139 \pm  0.000034\f$ | \f$0.000494\pm 0.000366\f$ |
| Chemical potential Widom \f$\left\langle \mu\right\rangle\f$ | \f$\epsilon\f$ | \f$-3.7721 \pm 0.0171\f$ | \f$-3.76212 \pm 0.0358\f$|
| Chemical potential CFCMC \f$\left\langle \mu\right\rangle\f$ | \f$\epsilon\f$ | \f$-3.7676 \pm 0.0226\f$ | \f$-3.77704 \pm 0.0212\f$ |



Liquid
```
    Gibbs swap (CFCMC)   all:           433767283
    Gibbs swap (CFCMC)   total:         108435452  108433519  216898312
    Gibbs swap (CFCMC)   constructed:   104794521  108411050  160644912
    Gibbs swap (CFCMC)   accepted:       69696354   17629208  109049917
    Gibbs swap (CFCMC)   fraction:       0.642745   0.162581   0.502770
    Gibbs swap (CFCMC)   max-change:     1.000000   0.100000   0.466139
```

Vapor
```
    Gibbs swap (CFCMC)   all:           501025563
    Gibbs swap (CFCMC)   total:         125272078  125241165  250512320
    Gibbs swap (CFCMC)   constructed:   125267752  113611272  128364421
    Gibbs swap (CFCMC)   accepted:       69696354   17629207  126271863
    Gibbs swap (CFCMC)   fraction:       0.556360   0.140762   0.504055
    Gibbs swap (CFCMC)   max-change:     0.299122   0.100000   1.000000
```

The fractional molecule was 46.4058% in the liquid phase and 53.5942% in the vapor phase.
