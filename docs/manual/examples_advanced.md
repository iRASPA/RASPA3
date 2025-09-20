# Examples advanced
\page examples_advanced Examples Advanced


## Table of Contents
1. [Monte Carlo: CFCMC Gibbs CO₂](#Example_advanced_1)
2. [Monte Carlo: CFCMC adsorption of CO₂ in MFI](#Example_advanced_2)
3. [Monte Carlo: CFCMC binary mixture adsorption of CO₂ and N₂ in DMOF](#Example_advanced_3)
4. [Transition Matrix Monte Carlo: methane in Tobacco-667](#Example_advanced_4)
5. [HMC CFCMC co2 in MFI](#Example_advanced_5)


#### Monte Carlo: CFCMC Gibbs CO₂<a name="Example_advanced_1"></a>

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 100000,
  "NumberOfInitializationCycles" : 50000,
  "NumberOfEquilibrationCycles" : 100000,
  "PrintEvery" : 5000,

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
      "ThermodynamicIntegration" : true,
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "Gibbs_CFCMC_SwapProbability" : 1.0,
      "WidomProbability" : 1.0,
      "CreateNumberOfMolecules" : [256, 256]
    }
  ]
}
```

```
    Block[  0] 10.119293539391567
    Block[  1] 11.777381262659526
    Block[  2] 8.624489044751295
    Block[  3] 10.323746247042733
    Block[  4] 8.770005749072295
    -----------------------------------------------------------------------
    Density average   9.923058e+00 +/-  1.601127e+00 [molecules]
    Density average   4.743900e-04 +/-  7.312203e-05 [molec/A^3]
    Density average   3.465974e+01 +/-  5.342421e+00 [kg.m⁻³]
```

```
Component 0 (CO2)
    Block[  0] 501.8584437344306
    Block[  1] 500.24156263024753
    Block[  2] 503.4371500202725
    Block[  3] 501.690202583573
    Block[  4] 503.2028284882624
    -----------------------------------------------------------------------
    Density average   5.020940e+02 +/-  1.604760e+00 [molecules]
    Density average   1.513777e-02 +/-  6.821167e-05 [molec/A^3]
    Density average   1.105991e+03 +/-  4.983661e+00 [kg.m⁻³]
```

The pressure is \f$13.85 \pm 2.2\f$ and \f$13.80 \pm 4.1\f$  for the vapor and liquid phase, respectively.

[comment]: <>       liquid                                       vapor
[comment]: <>  P(q) 1.250459e+06 +/-  2.406476e+05 [Pa]          1.278280e+06 +/-  4.265266e+05 [Pa]
[comment]: <>  thermodynamic 1.261961e+06 +/-  2.928364e+05      1.190653e+06 +/-  1.326930e+05 [Pa]
[comment]: <>  WIdom 1.152688e+06 +/-  8.188189e+04 [Pa]         1.175143e+06 +/-  1.241933e+05 [Pa]


Liquid
```
    Gibbs swap (CFCMC)   all:             4312726
    Gibbs swap (CFCMC)   total:           1077767    1077914    2157045
    Gibbs swap (CFCMC)   constructed:      956937    1059902    1709417
    Gibbs swap (CFCMC)   accepted:         455521     159983    1079031
    Gibbs swap (CFCMC)   fraction:       0.422653   0.148419   0.500236
    Gibbs swap (CFCMC)   max-change:     0.100000   0.100000   0.417734
```

Vapor
```
    Gibbs swap (CFCMC)   total:           1055101    1056705    2111838
    Gibbs swap (CFCMC)   constructed:     1051789     568825    1158819
    Gibbs swap (CFCMC)   accepted:         455519     159986    1046301
    Gibbs swap (CFCMC)   fraction:       0.431730   0.151401   0.495446
    Gibbs swap (CFCMC)   max-change:     0.100000   0.100000   0.924593
```

#### Monte Carlo: Monte Carlo: CFCMC adsorption of  CO₂ in MFI<a name="Example_advanced_2"></a>

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 50000,
  "NumberOfInitializationCycles" : 25000,
  "NumberOfEquilibrationCycles" : 25000,
  "PrintEvery" : 5000,

  "Systems" :
  [
    {
      "Type" : "Framework",
      "Name" : "MFI_SI",
      "NumberOfUnitCells" : [2, 2, 2],
      "ExternalTemperature" : 353.0,
      "ExternalPressure" : 1.0e5,
      "ChargeMethod" : "Ewald"
    }
  ],

  "Components" :
  [
    {
      "Name" : "CO2",
      "MoleculeDefinition" : "ExampleDefinitions",
      "FugacityCoefficient" : 1.0,
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

```
Component 0 (CO2)
    Block[  0]  2.637955e+01
    Block[  1]  2.558354e+01
    Block[  2]  2.566131e+01
    Block[  3]  2.596787e+01
    Block[  4]  2.588552e+01
    ---------------------------------------------------------------------------
    Abs. loading average   2.589568e+01 +/-  3.885433e-01 [molecules/cell]
    Abs. loading average   3.236960e+00 +/-  4.856792e-02 [molecules/uc]
    Abs. loading average   4.001867e-01 +/-  6.004472e-03 [mol/kg-framework]
    Abs. loading average   1.760774e+01 +/-  2.641896e-01 [mg/g-framework]
```

```
    Lambda histogram and bias:
    ---------------------------------------------------------------------------
     0-0.000000 (lambda) P:  1.11971e+00 +/- 3.88645e-02 bias:  3.53844e+00 [-]
     1-0.025000 (lambda) P:  1.10844e+00 +/- 7.20769e-02 bias:  3.53094e+00 [-]
     2-0.050000 (lambda) P:  1.01721e+00 +/- 3.10996e-02 bias:  3.45437e+00 [-]
     3-0.075000 (lambda) P:  1.04776e+00 +/- 4.74884e-02 bias:  3.53375e+00 [-]
     4-0.100000 (lambda) P:  1.13734e+00 +/- 2.07070e-02 bias:  3.73688e+00 [-]
     5-0.125000 (lambda) P:  9.80310e-01 +/- 4.02058e-02 bias:  3.90750e+00 [-]
     6-0.150000 (lambda) P:  8.83140e-01 +/- 5.58786e-02 bias:  4.19781e+00 [-]
     7-0.175000 (lambda) P:  1.04366e+00 +/- 3.06908e-02 bias:  4.54750e+00 [-]
     8-0.200000 (lambda) P:  9.64730e-01 +/- 6.21992e-02 bias:  4.40281e+00 [-]
     9-0.225000 (lambda) P:  9.72725e-01 +/- 3.42513e-02 bias:  4.13125e+00 [-]
    10-0.250000 (lambda) P:  1.07953e+00 +/- 2.62395e-02 bias:  3.83875e+00 [-]
    11-0.275000 (lambda) P:  9.48535e-01 +/- 6.12627e-02 bias:  3.28375e+00 [-]
    12-0.300000 (lambda) P:  9.87485e-01 +/- 4.19983e-02 bias:  2.82594e+00 [-]
    13-0.325000 (lambda) P:  9.82360e-01 +/- 4.05467e-02 bias:  2.41813e+00 [-]
    14-0.350000 (lambda) P:  1.06580e+00 +/- 5.45650e-02 bias:  2.08813e+00 [-]
    15-0.375000 (lambda) P:  1.08281e+00 +/- 1.69916e-02 bias:  1.77156e+00 [-]
    16-0.400000 (lambda) P:  9.03230e-01 +/- 4.22367e-02 bias:  1.34000e+00 [-]
    17-0.425000 (lambda) P:  1.02152e+00 +/- 3.24750e-02 bias:  1.25656e+00 [-]
    18-0.450000 (lambda) P:  9.22500e-01 +/- 3.52187e-02 bias:  1.06906e+00 [-]
    19-0.475000 (lambda) P:  9.42590e-01 +/- 6.11106e-02 bias:  1.05844e+00 [-]
    20-0.500000 (lambda) P:  9.59605e-01 +/- 2.84455e-02 bias:  1.02781e+00 [-]
    21-0.525000 (lambda) P:  1.01024e+00 +/- 1.96642e-02 bias:  1.10344e+00 [-]
    22-0.550000 (lambda) P:  9.73750e-01 +/- 2.75137e-02 bias:  1.03906e+00 [-]
    23-0.575000 (lambda) P:  9.87485e-01 +/- 4.65690e-02 bias:  1.03969e+00 [-]
    24-0.600000 (lambda) P:  9.87075e-01 +/- 4.59072e-02 bias:  1.06125e+00 [-]
    25-0.625000 (lambda) P:  1.03279e+00 +/- 4.00545e-02 bias:  1.04188e+00 [-]
    26-0.650000 (lambda) P:  9.71700e-01 +/- 1.42554e-02 bias:  9.74687e-01 [-]
    27-0.675000 (lambda) P:  1.02603e+00 +/- 2.64331e-02 bias:  9.69375e-01 [-]
    28-0.700000 (lambda) P:  9.91380e-01 +/- 6.16172e-02 bias:  8.83750e-01 [-]
    29-0.725000 (lambda) P:  9.95890e-01 +/- 2.46353e-02 bias:  8.10937e-01 [-]
    30-0.750000 (lambda) P:  1.02008e+00 +/- 5.70174e-02 bias:  7.31562e-01 [-]
    31-0.775000 (lambda) P:  1.01803e+00 +/- 4.76585e-02 bias:  6.42188e-01 [-]
    32-0.800000 (lambda) P:  1.03402e+00 +/- 1.63555e-02 bias:  5.32500e-01 [-]
    33-0.825000 (lambda) P:  9.67190e-01 +/- 3.50874e-02 bias:  3.61563e-01 [-]
    34-0.850000 (lambda) P:  9.46895e-01 +/- 2.61931e-02 bias:  2.52500e-01 [-]
    35-0.875000 (lambda) P:  9.75595e-01 +/- 2.28767e-02 bias:  1.68125e-01 [-]
    36-0.900000 (lambda) P:  9.53455e-01 +/- 2.17147e-02 bias:  1.08125e-01 [-]
    37-0.925000 (lambda) P:  9.33775e-01 +/- 5.57364e-02 bias:  4.96875e-02 [-]
    38-0.950000 (lambda) P:  1.00942e+00 +/- 3.46322e-02 bias:  6.31250e-02 [-]
    39-0.975000 (lambda) P:  9.45460e-01 +/- 2.39014e-02 bias:  0.00000e+00 [-]
    40-1.000000 (lambda) P:  1.04878e+00 +/- 1.01561e-02 bias:  9.28125e-02 [-]


    Lambda statistics:
    ---------------------------------------------------------------------------
     0-0.000000 (lambda) Free energy: 5.514740e+00 +/- 1.224748e+01 [K]
     1-0.025000 (lambda) Free energy: 9.087320e+00 +/- 2.310895e+01 [K]
     2-0.050000 (lambda) Free energy: 3.940494e+01 +/- 1.076929e+01 [K]
     3-0.075000 (lambda) Free energy: 2.896100e+01 +/- 1.609857e+01 [K]
     4-0.100000 (lambda) Free energy: 1.592267e-13 +/- 6.409162e+00 [K]
     5-0.125000 (lambda) Free energy: 5.244830e+01 +/- 1.450673e+01 [K]
     6-0.150000 (lambda) Free energy: 8.929628e+01 +/- 2.239344e+01 [K]
     7-0.175000 (lambda) Free energy: 3.034505e+01 +/- 1.033438e+01 [K]
     8-0.200000 (lambda) Free energy: 5.810357e+01 +/- 2.204916e+01 [K]
     9-0.225000 (lambda) Free energy: 5.519021e+01 +/- 1.251106e+01 [K]
    10-0.250000 (lambda) Free energy: 1.841477e+01 +/- 8.617753e+00 [K]
    11-0.275000 (lambda) Free energy: 6.407972e+01 +/- 2.251353e+01 [K]
    12-0.300000 (lambda) Free energy: 4.987406e+01 +/- 1.493929e+01 [K]
    13-0.325000 (lambda) Free energy: 5.171089e+01 +/- 1.454625e+01 [K]
    14-0.350000 (lambda) Free energy: 2.293486e+01 +/- 1.782078e+01 [K]
    15-0.375000 (lambda) Free energy: 1.734385e+01 +/- 5.505586e+00 [K]
    16-0.400000 (lambda) Free energy: 8.135607e+01 +/- 1.660521e+01 [K]
    17-0.425000 (lambda) Free energy: 3.791413e+01 +/- 1.124422e+01 [K]
    18-0.450000 (lambda) Free energy: 7.390418e+01 +/- 1.322437e+01 [K]
    19-0.475000 (lambda) Free energy: 6.629913e+01 +/- 2.343526e+01 [K]
    20-0.500000 (lambda) Free energy: 5.998384e+01 +/- 1.041518e+01 [K]
    21-0.525000 (lambda) Free energy: 4.183205e+01 +/- 6.868973e+00 [K]
    22-0.550000 (lambda) Free energy: 5.481844e+01 +/- 9.961855e+00 [K]
    23-0.575000 (lambda) Free energy: 4.987406e+01 +/- 1.684993e+01 [K]
    24-0.600000 (lambda) Free energy: 5.002066e+01 +/- 1.644233e+01 [K]
    25-0.625000 (lambda) Free energy: 3.403923e+01 +/- 1.381230e+01 [K]
    26-0.650000 (lambda) Free energy: 5.556238e+01 +/- 5.186979e+00 [K]
    27-0.675000 (lambda) Free energy: 3.635906e+01 +/- 9.074194e+00 [K]
    28-0.700000 (lambda) Free energy: 4.848444e+01 +/- 2.157213e+01 [K]
    29-0.725000 (lambda) Free energy: 4.688221e+01 +/- 8.761448e+00 [K]
    30-0.750000 (lambda) Free energy: 3.841037e+01 +/- 1.954104e+01 [K]
    31-0.775000 (lambda) Free energy: 3.912049e+01 +/- 1.635473e+01 [K]
    32-0.800000 (lambda) Free energy: 3.361907e+01 +/- 5.604900e+00 [K]
    33-0.825000 (lambda) Free energy: 5.720459e+01 +/- 1.273141e+01 [K]
    34-0.850000 (lambda) Free energy: 6.469058e+01 +/- 9.734838e+00 [K]
    35-0.875000 (lambda) Free energy: 5.415023e+01 +/- 8.347045e+00 [K]
    36-0.900000 (lambda) Free energy: 6.225346e+01 +/- 8.031176e+00 [K]
    37-0.925000 (lambda) Free energy: 6.961589e+01 +/- 2.125460e+01 [K]
    38-0.950000 (lambda) Free energy: 4.211869e+01 +/- 1.212537e+01 [K]
    39-0.975000 (lambda) Free energy: 6.522595e+01 +/- 8.950700e+00 [K]
    40-1.000000 (lambda) Free energy: 2.861584e+01 +/- 3.422448e+00 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential: (ln(P(lambda=1))-ln(P(lambda=0)))/Beta
        Block[  0] -1205.778738016742
        Block[  1] -1200.258471899379
        Block[  2] -1182.001455454578
        Block[  3] -1193.677225369166
        Block[  4] -1184.820636933015
    ---------------------------------------------------------------------------
    Excess chemical potential:    -1.193206e+03 +/-  1.247486e+01 [K]
    Ideal gas chemical potential: -2.614621e+03 +/-  5.278695e+00 [K]
    Total chemical potential:     -3.807827e+03 +/-  1.069918e+01 [K]
    Imposed chemical potential:   -3.810353e+03 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential:    -9.920866e+00 +/-  1.037218e-01 [kJ/mol]
    Ideal gas chemical potential: -2.173917e+01 +/-  4.388953e-02 [kJ/mol]
    Total chemical potential:     -3.166004e+01 +/-  8.895799e-02 [kJ/mol]
    Imposed chemical potential:   -3.168105e+01 [kJ/mol]
    ---------------------------------------------------------------------------
    Imposed fugacity:              1.000000e+05 [Pa]
    Measured fugacity:             1.007184e+05 +/-  3.044534e+03 [Pa]


    Thermodynamic integration (dU/dlambda)
    ===========================================================================

     0-0.000000 (lambda) <dU/dlambda>:  0.000000e+00 +/- 0.000000e+00 [K/-]
     1-0.025000 (lambda) <dU/dlambda>:  4.877297e+01 +/- 1.819147e+01 [K/-]
     2-0.050000 (lambda) <dU/dlambda>:  2.971444e+02 +/- 2.128510e+01 [K/-]
     3-0.075000 (lambda) <dU/dlambda>:  9.924703e+02 +/- 7.221565e+01 [K/-]
     4-0.100000 (lambda) <dU/dlambda>:  2.819070e+03 +/- 1.607450e+02 [K/-]
     5-0.125000 (lambda) <dU/dlambda>:  5.359466e+03 +/- 5.263760e+02 [K/-]
     6-0.150000 (lambda) <dU/dlambda>:  4.997893e+03 +/- 5.687054e+02 [K/-]
     7-0.175000 (lambda) <dU/dlambda>:  5.139452e+02 +/- 5.078197e+02 [K/-]
     8-0.200000 (lambda) <dU/dlambda>: -2.509815e+03 +/- 1.619668e+02 [K/-]
     9-0.225000 (lambda) <dU/dlambda>: -4.599565e+03 +/- 2.370776e+02 [K/-]
    10-0.250000 (lambda) <dU/dlambda>: -5.894260e+03 +/- 1.999326e+02 [K/-]
    11-0.275000 (lambda) <dU/dlambda>: -6.463507e+03 +/- 9.885582e+01 [K/-]
    12-0.300000 (lambda) <dU/dlambda>: -6.522006e+03 +/- 1.565721e+02 [K/-]
    13-0.325000 (lambda) <dU/dlambda>: -6.034478e+03 +/- 7.350153e+01 [K/-]
    14-0.350000 (lambda) <dU/dlambda>: -5.162605e+03 +/- 5.610515e+01 [K/-]
    15-0.375000 (lambda) <dU/dlambda>: -4.112846e+03 +/- 3.090119e+01 [K/-]
    16-0.400000 (lambda) <dU/dlambda>: -3.032841e+03 +/- 3.325573e+01 [K/-]
    17-0.425000 (lambda) <dU/dlambda>: -2.061337e+03 +/- 1.305164e+01 [K/-]
    18-0.450000 (lambda) <dU/dlambda>: -1.208290e+03 +/- 9.943336e+00 [K/-]
    19-0.475000 (lambda) <dU/dlambda>: -5.272465e+02 +/- 3.772265e+00 [K/-]
    20-0.500000 (lambda) <dU/dlambda>:  0.000000e+00 +/- 0.000000e+00 [K/-]
    21-0.525000 (lambda) <dU/dlambda>: -2.851385e+01 +/- 3.300881e+00 [K/-]
    22-0.550000 (lambda) <dU/dlambda>: -6.441270e+01 +/- 8.153721e+00 [K/-]
    23-0.575000 (lambda) <dU/dlambda>: -1.210019e+02 +/- 2.078867e+01 [K/-]
    24-0.600000 (lambda) <dU/dlambda>: -1.944030e+02 +/- 2.297619e+01 [K/-]
    25-0.625000 (lambda) <dU/dlambda>: -3.147098e+02 +/- 1.717133e+01 [K/-]
    26-0.650000 (lambda) <dU/dlambda>: -4.828201e+02 +/- 4.213004e+01 [K/-]
    27-0.675000 (lambda) <dU/dlambda>: -7.343406e+02 +/- 2.439306e+01 [K/-]
    28-0.700000 (lambda) <dU/dlambda>: -1.036277e+03 +/- 5.651073e+01 [K/-]
    29-0.725000 (lambda) <dU/dlambda>: -1.283711e+03 +/- 4.723179e+01 [K/-]
    30-0.750000 (lambda) <dU/dlambda>: -1.542215e+03 +/- 1.131860e+02 [K/-]
    31-0.775000 (lambda) <dU/dlambda>: -1.642367e+03 +/- 5.587092e+01 [K/-]
    32-0.800000 (lambda) <dU/dlambda>: -1.623731e+03 +/- 7.325643e+01 [K/-]
    33-0.825000 (lambda) <dU/dlambda>: -1.522767e+03 +/- 5.862325e+01 [K/-]
    34-0.850000 (lambda) <dU/dlambda>: -1.324467e+03 +/- 3.446451e+01 [K/-]
    35-0.875000 (lambda) <dU/dlambda>: -1.048089e+03 +/- 2.831138e+01 [K/-]
    36-0.900000 (lambda) <dU/dlambda>: -7.894606e+02 +/- 1.822824e+01 [K/-]
    37-0.925000 (lambda) <dU/dlambda>: -5.301669e+02 +/- 1.683158e+01 [K/-]
    38-0.950000 (lambda) <dU/dlambda>: -3.125078e+02 +/- 7.076055e+00 [K/-]
    39-0.975000 (lambda) <dU/dlambda>: -1.379682e+02 +/- 4.741803e+00 [K/-]
    40-1.000000 (lambda) <dU/dlambda>:  0.000000e+00 +/- 0.000000e+00 [K/-]
    ---------------------------------------------------------------------------
    Excess chemical potential: integral du/dlambda over lambda (Simpson's rule)
        Block[  0] -1197.9879280304358
        Block[  1] -1237.7903793611556
        Block[  2] -1178.0811037066471
        Block[  3] -1197.6293268810336
        Block[  4] -1196.7803680639572
    ---------------------------------------------------------------------------
    Excess chemical potential:   -1.201365e+03 +/-  2.716600e+01 [K]
    Ideal chemical potential:    -2.614621e+03 +/-  5.278695e+00 [K]
    Total chemical potential:    -3.815986e+03 +/-  2.911666e+01 [K]
    Imposed chemical potential:  -3.810353e+03 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential:   -9.988711e+00 +/-  2.258708e-01 [kJ/mol]
    Ideal chemical potential:    -2.173917e+01 +/-  4.388953e-02 [kJ/mol]
    Total chemical potential:    -3.172789e+01 +/-  2.420894e-01 [kJ/mol]
    Imposed chemical potential:  -3.168105e+01 [kJ/mol]
    ---------------------------------------------------------------------------
    Imposed fugacity:             1.000000e+05 [Pa]
    Measured fugacity:            9.841693e+04 +/-  7.832774e+03 [Pa]


    Widom insertion Rosenbluth weight statistics:
    ---------------------------------------------------------------------------
        Block[  0]  3.027080e+01
        Block[  1]  3.060848e+01
        Block[  2]  3.123650e+01
        Block[  3]  3.070635e+01
        Block[  4]  3.051262e+01
    ---------------------------------------------------------------------------
    Average Rosenbluth weight:    3.066415e+01 +/-  4.433563e-01 [-]


    Henry coefficient based on Rosenbluth weight:
    ---------------------------------------------------------------------------
        Block[  0]  4.094336e-06
        Block[  1]  4.140009e-06
        Block[  2]  4.224953e-06
        Block[  3]  4.153246e-06
        Block[  4]  4.127043e-06
    ---------------------------------------------------------------------------
    Average Henry coefficient:    4.147538e-06 +/-  5.996701e-08 [mol/kg/Pa]
    Average Henry coefficient:    3.354788e-05 +/-  4.850506e-07 [molec./uc/Pa]


    Widom insertion chemical potential  statistics:
    ---------------------------------------------------------------------------
        Block[  0] -1203.7959155333265
        Block[  1] -1207.7118597927076
        Block[  2] -1214.881403748269
        Block[  3] -1208.8387408385317
        Block[  4] -1206.6045864239932
    ---------------------------------------------------------------------------
    Excess chemical potential:          -1.208353e+03 +/-  5.084363e+00 [K]
    Tail-correction chemical potential:  0.000000e+00 +/-  0.000000e+00 [K]
    Ideal chemical potential:           -2.614621e+03 +/-  5.278695e+00 [K]
    Total chemical potential:           -3.822974e+03 +/-  9.494304e+00 [K]
    Imposed chemical potential:         -3.810353e+03 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential:          -1.004681e+01 +/-  4.227376e-02 [kJ/mol]
    Tail-correction chemical potential:  0.000000e+00 +/-  0.000000e+00 [kj/mol]
    Ideal chemical potential:           -2.173917e+01 +/-  4.388953e-02 [kJ/mol]
    Total chemical potential:           -3.178599e+01 +/-  7.894006e-02 [kJ/mol]
    Imposed chemical potential:         -3.168105e+01 [kJ/mol]
    ---------------------------------------------------------------------------
    Imposed fugacity:                    1.000000e+05 [Pa]
    Measured fugacity:                   9.648786e+04 +/-  2.602012e+03 [Pa]
```

```
Component 0 [CO2]
    Reinsertion (CBMC)   all:              771188
    Reinsertion (CBMC)   total:            771188
    Reinsertion (CBMC)   constructed:      719332
    Reinsertion (CBMC)   accepted:         163465
    Reinsertion (CBMC)   fraction:       0.211965
    Reinsertion (CBMC)   max-change:     0.000000

    Widom                all:             1544154
    Widom                total:           1544154
    Widom                constructed:     1437086
    Widom                accepted:              0
    Widom                fraction:       0.000000
    Widom                max-change:     0.000000

    Translation          all:              771861
    Translation          total:            257253     257117     257491
    Translation          constructed:      257245     257109     257490
    Translation          accepted:         131318     131177     128845
    Translation          fraction:       0.510462   0.510184   0.500386
    Translation          max-change:     1.291610   1.303601   1.209865

    Rotation             all:              772857
    Rotation             total:            257292     258164     257401
    Rotation             constructed:      257292     258164     257400
    Rotation             accepted:         131138     130571     131353
    Rotation             fraction:       0.509685   0.505768   0.510305
    Rotation             max-change:     1.253359   1.277356   1.372233

    Swap (CB/CFCMC)      all:             1544699
    Swap (CB/CFCMC)      total:            362894     367725     814080
    Swap (CB/CFCMC)      constructed:      348494     367725     759652
    Swap (CB/CFCMC)      accepted:         174273     174264     450560
    Swap (CB/CFCMC)      fraction:       0.480231   0.473898   0.553459
    Swap (CB/CFCMC)      max-change:     0.100000   0.100000   1.000000
```


#### Monte Carlo: CFCMC binary mixture adsorption of CO₂ and N₂ in DMOF<a name="Example_advanced_3"></a>

```json
{
}
```

#### Transition Matrix Monte Carlo: methane in Tobacco-667<a name="Example_advanced_4"></a>

```json
{
  "SimulationType": "MonteCarlo",
  "NumberOfCycles": 10000,
  "NumberOfInitializationCycles": 5000,
  "NumberOfEquilibrationCycles": 5000,
  "PrintEvery": 1000,
  "Systems": [
    {
      "Type": "Framework",
      "Name": "MFI_SI",
      "NumberOfUnitCells": [
        2,
        2,
        2
      ],
      "ExternalTemperature": 353.0,
      "ExternalPressure": 1.0e5,
      "ChargeMethod": "Ewald",
      "HybridMCProbability": 0.01,
      "HybridMCMoveNumberOfSteps": 100,
      "TimeStep": 0.005,
      "Ensemble": "NVE"
    }
  ],
  "Components": [
    {
      "Name": "CO2",
      "MoleculeDefinition": "ExampleDefinitions",
      "FugacityCoefficient": 1.0,
      "ThermodynamicIntegration": true,
      "TranslationProbability": 0.5,
      "RotationProbability": 0.5,
      "CFCMC_CBMC_SwapProbability": 1.0,
      "CreateNumberOfMolecules": 0
    }
  ]
}
```

#### HMC CFCMC co2 in MFI<a name="Example_advanced_5"></a>


```json
{
  "SimulationType" : "MonteCarloTransitionMatrix",
  "NumberOfCycles" : 200000,
  "NumberOfInitializationCycles" : 100000,
  "PrintEvery" : 5000,

  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "tobacco-667",
      "HeliumVoidFraction" : 0.29,
      "NumberOfUnitCells" : [1, 1, 1],
      "ExternalTemperature" : 111.58,
      "ExternalPressure" : 36000,
      "ChargeMethod" : "None",
      "MacroStateUseBias" : true,
      "MacroStateMinimumNumberOfMolecules" : 84,
      "MacroStateMaximumNumberOfMolecules" : 168
    }
  ],

  "Components" : [
    {
      "Name" : "methane",
      "MoleculeDefinition" : "ExampleDefinitions",
      "TranslationProbability" : 1.0,
      "ReinsertionProbability" : 1.0,
      "SwapProbability" : 1.0,
      "CreateNuMberofmolecules" : 120
    }
  ]
}
```

