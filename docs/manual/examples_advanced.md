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
Component 0 (CO2)
    Block[  0] 8.48142863589921
    Block[  1] 10.984341035898728
    Block[  2] 8.17805265820302
    Block[  3] 11.649264803898573
    Block[  4] 11.717605184930896
    -----------------------------------------------------------------------
    Density average   1.020234e+01 +/-  2.155687e+00 [molecules]
    Density average   4.890088e-04 +/-  9.102441e-05 [molec/A^3]
    Density average   3.572782e+01 +/-  6.650399e+00 [kg.m⁻³]
```

```
Component 0 (CO2)
    Block[  0] 503.53255521364514
    Block[  1] 501.01922246001124
    Block[  2] 503.80987984438934
    Block[  3] 500.3507046083389
    Block[  4] 500.28150325039945
    -----------------------------------------------------------------------
    Density average   5.018102e+02 +/-  2.155462e+00 [molecules]
    Density average   1.508888e-02 +/-  1.991028e-04 [molec/A^3]
    Density average   1.102419e+03 +/-  1.454679e+01 [kg.m⁻³]
```

The pressure is \f$13.85 \pm 2.2\f$ and \f$13.80 \pm 4.1\f$  for the vapor and liquid phase, respectively.

[comment]: <>       liquid                                       vapor
[comment]: <>  P(q) 1.190208e+06 +/-  1.563101e+05               1.136575e+06 +/-  1.292557e+05
[comment]: <>  thermodynamic 1.104544e+06 +/-  2.513265e+05      1.228774e+06 +/-  1.884924e+05
[comment]: <>  WIdom 1.151731e+06 +/-  2.083492e+05              1.202304e+06 +/-  1.605820e+05



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
    Block[  0]  2.596590e+01
    Block[  1]  2.589989e+01
    Block[  2]  2.575433e+01
    Block[  3]  2.598268e+01
    Block[  4]  2.587816e+01
    ---------------------------------------------------------------------------
    Abs. loading average   2.589622e+01 +/-  1.124223e-01 [molecules/cell]
    Abs. loading average   3.237028e+00 +/-  1.405279e-02 [molecules/uc]
    Abs. loading average   4.001951e-01 +/-  1.737353e-03 [mol/kg-framework]
    Abs. loading average   1.760811e+01 +/-  7.644144e-02 [mg/g-framework]
```

```
    Lambda histogram and bias:
    ---------------------------------------------------------------------------
     0-0.000000 (lambda) P:  1.03509e+00 +/- 2.80282e-02 bias:  3.45717e+00 [-]
     1-0.025000 (lambda) P:  1.02570e+00 +/- 2.43641e-02 bias:  3.46100e+00 [-]
     2-0.050000 (lambda) P:  1.04554e+00 +/- 1.09987e-02 bias:  3.47367e+00 [-]
     3-0.075000 (lambda) P:  1.00905e+00 +/- 2.00211e-02 bias:  3.48293e+00 [-]
     4-0.100000 (lambda) P:  9.83836e-01 +/- 1.93164e-02 bias:  3.59941e+00 [-]
     5-0.125000 (lambda) P:  1.00614e+00 +/- 1.98527e-02 bias:  3.89707e+00 [-]
     6-0.150000 (lambda) P:  1.06198e+00 +/- 1.63620e-02 bias:  4.35492e+00 [-]
     7-0.175000 (lambda) P:  9.89617e-01 +/- 2.00399e-02 bias:  4.47791e+00 [-]
     8-0.200000 (lambda) P:  9.97079e-01 +/- 2.69905e-02 bias:  4.40930e+00 [-]
     9-0.225000 (lambda) P:  9.89125e-01 +/- 2.53263e-02 bias:  4.14570e+00 [-]
    10-0.250000 (lambda) P:  1.02566e+00 +/- 1.79124e-02 bias:  3.80521e+00 [-]
    11-0.275000 (lambda) P:  1.02967e+00 +/- 1.43307e-02 bias:  3.35303e+00 [-]
    12-0.300000 (lambda) P:  1.06432e+00 +/- 1.08445e-02 bias:  2.91754e+00 [-]
    13-0.325000 (lambda) P:  9.85558e-01 +/- 7.63416e-03 bias:  2.39111e+00 [-]
    14-0.350000 (lambda) P:  9.89822e-01 +/- 2.78330e-02 bias:  2.01068e+00 [-]
    15-0.375000 (lambda) P:  1.04234e+00 +/- 1.33687e-02 bias:  1.71383e+00 [-]
    16-0.400000 (lambda) P:  9.92938e-01 +/- 1.84908e-02 bias:  1.43910e+00 [-]
    17-0.425000 (lambda) P:  9.98637e-01 +/- 1.90948e-02 bias:  1.25029e+00 [-]
    18-0.450000 (lambda) P:  1.05821e+00 +/- 1.19575e-02 bias:  1.20377e+00 [-]
    19-0.475000 (lambda) P:  1.03361e+00 +/- 2.16408e-02 bias:  1.12168e+00 [-]
    20-0.500000 (lambda) P:  1.01459e+00 +/- 1.28761e-02 bias:  1.07840e+00 [-]
    21-0.525000 (lambda) P:  9.92364e-01 +/- 1.68753e-02 bias:  1.05789e+00 [-]
    22-0.550000 (lambda) P:  9.83426e-01 +/- 1.54087e-02 bias:  1.05881e+00 [-]
    23-0.575000 (lambda) P:  9.62188e-01 +/- 1.09326e-02 bias:  1.01568e+00 [-]
    24-0.600000 (lambda) P:  1.00622e+00 +/- 1.93342e-02 bias:  1.05246e+00 [-]
    25-0.625000 (lambda) P:  9.93266e-01 +/- 1.27839e-02 bias:  1.01189e+00 [-]
    26-0.650000 (lambda) P:  9.82524e-01 +/- 8.93184e-03 bias:  9.89023e-01 [-]
    27-0.675000 (lambda) P:  9.84738e-01 +/- 1.81184e-02 bias:  9.31699e-01 [-]
    28-0.700000 (lambda) P:  1.01762e+00 +/- 1.85933e-02 bias:  9.17266e-01 [-]
    29-0.725000 (lambda) P:  9.66493e-01 +/- 2.34966e-02 bias:  7.79727e-01 [-]
    30-0.750000 (lambda) P:  9.74570e-01 +/- 1.82541e-02 bias:  6.89746e-01 [-]
    31-0.775000 (lambda) P:  9.63131e-01 +/- 1.05647e-02 bias:  5.79590e-01 [-]
    32-0.800000 (lambda) P:  9.52430e-01 +/- 1.41047e-02 bias:  4.42109e-01 [-]
    33-0.825000 (lambda) P:  9.57965e-01 +/- 1.62243e-02 bias:  3.35156e-01 [-]
    34-0.850000 (lambda) P:  9.71782e-01 +/- 2.14812e-02 bias:  2.43906e-01 [-]
    35-0.875000 (lambda) P:  9.85599e-01 +/- 1.80528e-02 bias:  1.84102e-01 [-]
    36-0.900000 (lambda) P:  9.98637e-01 +/- 1.79721e-02 bias:  1.34004e-01 [-]
    37-0.925000 (lambda) P:  1.00651e+00 +/- 2.50035e-02 bias:  7.54492e-02 [-]
    38-0.950000 (lambda) P:  9.63664e-01 +/- 2.56777e-02 bias:  2.08984e-02 [-]
    39-0.975000 (lambda) P:  9.91216e-01 +/- 7.61675e-03 bias:  2.96680e-02 [-]
    40-1.000000 (lambda) P:  9.67149e-01 +/- 1.39261e-02 bias:  0.00000e+00 [-]


    Lambda statistics:
    ---------------------------------------------------------------------------
     0-0.000000 (lambda) Free energy: 9.831286e+00 +/- 9.651231e+00 [K]
     1-0.025000 (lambda) Free energy: 1.304787e+01 +/- 8.411487e+00 [K]
     2-0.050000 (lambda) Free energy: 6.283654e+00 +/- 3.717435e+00 [K]
     3-0.075000 (lambda) Free energy: 1.882370e+01 +/- 7.001442e+00 [K]
     4-0.100000 (lambda) Free energy: 2.775685e+01 +/- 6.913294e+00 [K]
     5-0.125000 (lambda) Free energy: 1.984354e+01 +/- 6.975695e+00 [K]
     6-0.150000 (lambda) Free energy: 7.759598e-01 +/- 5.427128e+00 [K]
     7-0.175000 (lambda) Free energy: 2.568869e+01 +/- 7.161410e+00 [K]
     8-0.200000 (lambda) Free energy: 2.303695e+01 +/- 9.543260e+00 [K]
     9-0.225000 (lambda) Free energy: 2.586424e+01 +/- 9.018092e+00 [K]
    10-0.250000 (lambda) Free energy: 1.306198e+01 +/- 6.140162e+00 [K]
    11-0.275000 (lambda) Free energy: 1.168181e+01 +/- 4.906686e+00 [K]
    12-0.300000 (lambda) Free energy: -3.753330e-14 +/- 3.583755e+00 [K]
    13-0.325000 (lambda) Free energy: 2.713953e+01 +/- 2.725112e+00 [K]
    14-0.350000 (lambda) Free energy: 2.561558e+01 +/- 9.799202e+00 [K]
    15-0.375000 (lambda) Free energy: 7.365032e+00 +/- 4.506251e+00 [K]
    16-0.400000 (lambda) Free energy: 2.450606e+01 +/- 6.541313e+00 [K]
    17-0.425000 (lambda) Free energy: 2.248580e+01 +/- 6.756117e+00 [K]
    18-0.450000 (lambda) Free energy: 2.031996e+00 +/- 3.990059e+00 [K]
    19-0.475000 (lambda) Free energy: 1.033501e+01 +/- 7.367644e+00 [K]
    20-0.500000 (lambda) Free energy: 1.689266e+01 +/- 4.455021e+00 [K]
    21-0.525000 (lambda) Free energy: 2.471018e+01 +/- 6.024291e+00 [K]
    22-0.550000 (lambda) Free energy: 2.790398e+01 +/- 5.513194e+00 [K]
    23-0.575000 (lambda) Free energy: 3.561088e+01 +/- 4.008431e+00 [K]
    24-0.600000 (lambda) Free energy: 1.981477e+01 +/- 6.819060e+00 [K]
    25-0.625000 (lambda) Free energy: 2.438947e+01 +/- 4.535359e+00 [K]
    26-0.650000 (lambda) Free energy: 2.822791e+01 +/- 3.211774e+00 [K]
    27-0.675000 (lambda) Free energy: 2.743336e+01 +/- 6.493289e+00 [K]
    28-0.700000 (lambda) Free energy: 1.583863e+01 +/- 6.434199e+00 [K]
    29-0.725000 (lambda) Free energy: 3.403501e+01 +/- 8.559175e+00 [K]
    30-0.750000 (lambda) Free energy: 3.109724e+01 +/- 6.612389e+00 [K]
    31-0.775000 (lambda) Free energy: 3.526508e+01 +/- 3.880640e+00 [K]
    32-0.800000 (lambda) Free energy: 3.920909e+01 +/- 5.209931e+00 [K]
    33-0.825000 (lambda) Free energy: 3.716359e+01 +/- 5.963098e+00 [K]
    34-0.850000 (lambda) Free energy: 3.210853e+01 +/- 7.796752e+00 [K]
    35-0.875000 (lambda) Free energy: 2.712485e+01 +/- 6.483978e+00 [K]
    36-0.900000 (lambda) Free energy: 2.248580e+01 +/- 6.367166e+00 [K]
    37-0.925000 (lambda) Free energy: 1.971410e+01 +/- 8.746533e+00 [K]
    38-0.950000 (lambda) Free energy: 3.506979e+01 +/- 9.342686e+00 [K]
    39-0.975000 (lambda) Free energy: 2.511878e+01 +/- 2.720472e+00 [K]
    40-1.000000 (lambda) Free energy: 3.379550e+01 +/- 5.085608e+00 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential: (ln(P(lambda=1))-ln(P(lambda=0)))/Beta
        Block[  0] -1194.6648919462536
        Block[  1] -1201.9045337402054
        Block[  2] -1205.4271448834616
        Block[  3] -1186.5025891427517
        Block[  4] -1193.8318164006112
    ---------------------------------------------------------------------------
    Excess chemical potential:    -1.196417e+03 +/-  9.193219e+00 [K]
    Ideal gas chemical potential: -2.614614e+03 +/-  1.534176e+00 [K]
    Total chemical potential:     -3.811031e+03 +/-  1.048346e+01 [K]
    Imposed chemical potential:   -3.810353e+03 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential:    -9.947569e+00 +/-  7.643670e-02 [kJ/mol]
    Ideal gas chemical potential: -2.173911e+01 +/-  1.275585e-02 [kJ/mol]
    Total chemical potential:     -3.168668e+01 +/-  8.716437e-02 [kJ/mol]
    Imposed chemical potential:   -3.168105e+01 [kJ/mol]
    ---------------------------------------------------------------------------
    Imposed fugacity:              1.000000e+05 [Pa]
    Measured fugacity:             9.980832e+04 +/-  2.963737e+03 [Pa]


    Thermodynamic integration (dU/dlambda)
    ===========================================================================

     0-0.000000 (lambda) <dU/dlambda>:  0.000000e+00 +/- 0.000000e+00 [K/-]
     1-0.025000 (lambda) <dU/dlambda>:  4.998471e+01 +/- 6.925114e+00 [K/-]
     2-0.050000 (lambda) <dU/dlambda>:  2.749367e+02 +/- 1.038544e+01 [K/-]
     3-0.075000 (lambda) <dU/dlambda>:  1.016009e+03 +/- 3.662323e+01 [K/-]
     4-0.100000 (lambda) <dU/dlambda>:  2.744508e+03 +/- 9.313486e+01 [K/-]
     5-0.125000 (lambda) <dU/dlambda>:  5.299892e+03 +/- 1.084978e+02 [K/-]
     6-0.150000 (lambda) <dU/dlambda>:  5.053869e+03 +/- 2.736328e+02 [K/-]
     7-0.175000 (lambda) <dU/dlambda>:  6.161639e+02 +/- 2.745291e+02 [K/-]
     8-0.200000 (lambda) <dU/dlambda>: -2.508732e+03 +/- 2.053790e+02 [K/-]
     9-0.225000 (lambda) <dU/dlambda>: -4.717010e+03 +/- 4.807733e+01 [K/-]
    10-0.250000 (lambda) <dU/dlambda>: -5.959086e+03 +/- 7.809899e+01 [K/-]
    11-0.275000 (lambda) <dU/dlambda>: -6.567734e+03 +/- 1.272796e+02 [K/-]
    12-0.300000 (lambda) <dU/dlambda>: -6.540853e+03 +/- 2.215682e+01 [K/-]
    13-0.325000 (lambda) <dU/dlambda>: -6.039769e+03 +/- 2.009045e+01 [K/-]
    14-0.350000 (lambda) <dU/dlambda>: -5.175836e+03 +/- 2.349653e+01 [K/-]
    15-0.375000 (lambda) <dU/dlambda>: -4.107353e+03 +/- 2.038599e+01 [K/-]
    16-0.400000 (lambda) <dU/dlambda>: -3.037654e+03 +/- 7.192137e+00 [K/-]
    17-0.425000 (lambda) <dU/dlambda>: -2.055914e+03 +/- 4.722887e+00 [K/-]
    18-0.450000 (lambda) <dU/dlambda>: -1.210938e+03 +/- 6.289587e+00 [K/-]
    19-0.475000 (lambda) <dU/dlambda>: -5.269930e+02 +/- 8.040971e-01 [K/-]
    20-0.500000 (lambda) <dU/dlambda>:  0.000000e+00 +/- 0.000000e+00 [K/-]
    21-0.525000 (lambda) <dU/dlambda>: -2.578384e+01 +/- 1.270255e+00 [K/-]
    22-0.550000 (lambda) <dU/dlambda>: -6.486676e+01 +/- 5.862348e+00 [K/-]
    23-0.575000 (lambda) <dU/dlambda>: -1.225703e+02 +/- 7.856852e+00 [K/-]
    24-0.600000 (lambda) <dU/dlambda>: -2.144323e+02 +/- 8.421214e+00 [K/-]
    25-0.625000 (lambda) <dU/dlambda>: -3.264465e+02 +/- 1.137434e+01 [K/-]
    26-0.650000 (lambda) <dU/dlambda>: -5.111359e+02 +/- 1.579978e+01 [K/-]
    27-0.675000 (lambda) <dU/dlambda>: -7.401240e+02 +/- 1.808912e+01 [K/-]
    28-0.700000 (lambda) <dU/dlambda>: -1.001663e+03 +/- 2.874093e+01 [K/-]
    29-0.725000 (lambda) <dU/dlambda>: -1.259176e+03 +/- 2.171075e+01 [K/-]
    30-0.750000 (lambda) <dU/dlambda>: -1.488642e+03 +/- 2.868165e+01 [K/-]
    31-0.775000 (lambda) <dU/dlambda>: -1.640020e+03 +/- 2.346251e+01 [K/-]
    32-0.800000 (lambda) <dU/dlambda>: -1.649601e+03 +/- 3.161446e+01 [K/-]
    33-0.825000 (lambda) <dU/dlambda>: -1.511252e+03 +/- 1.087026e+01 [K/-]
    34-0.850000 (lambda) <dU/dlambda>: -1.310901e+03 +/- 2.147970e+01 [K/-]
    35-0.875000 (lambda) <dU/dlambda>: -1.054397e+03 +/- 7.123572e+00 [K/-]
    36-0.900000 (lambda) <dU/dlambda>: -7.822825e+02 +/- 1.103420e+01 [K/-]
    37-0.925000 (lambda) <dU/dlambda>: -5.277253e+02 +/- 1.498797e+00 [K/-]
    38-0.950000 (lambda) <dU/dlambda>: -3.100845e+02 +/- 2.531547e+00 [K/-]
    39-0.975000 (lambda) <dU/dlambda>: -1.365752e+02 +/- 2.027126e+00 [K/-]
    40-1.000000 (lambda) <dU/dlambda>:  0.000000e+00 +/- 0.000000e+00 [K/-]
    ---------------------------------------------------------------------------
    Excess chemical potential: integral du/dlambda over lambda (Simpson's rule)
        Block[  0] -1209.381544910717
        Block[  1] -1205.975028545626
        Block[  2] -1198.5961422947423
        Block[  3] -1211.7348268092067
        Block[  4] -1211.3720695070717
    ---------------------------------------------------------------------------
    Excess chemical potential:   -1.207450e+03 +/-  6.744008e+00 [K]
    Ideal chemical potential:    -2.614614e+03 +/-  1.534176e+00 [K]
    Total chemical potential:    -3.822063e+03 +/-  5.486984e+00 [K]
    Imposed chemical potential:  -3.810353e+03 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential:   -1.003930e+01 +/-  5.607281e-02 [kJ/mol]
    Ideal chemical potential:    -2.173911e+01 +/-  1.275585e-02 [kJ/mol]
    Total chemical potential:    -3.177841e+01 +/-  4.562134e-02 [kJ/mol]
    Imposed chemical potential:  -3.168105e+01 [kJ/mol]
    ---------------------------------------------------------------------------
    Imposed fugacity:             1.000000e+05 [Pa]
    Measured fugacity:            9.673720e+04 +/-  1.510063e+03 [Pa]


    Widom insertion Rosenbluth weight statistics:
    ---------------------------------------------------------------------------
        Block[  0]  3.054430e+01
        Block[  1]  3.043135e+01
        Block[  2]  3.047939e+01
        Block[  3]  3.073296e+01
        Block[  4]  3.047557e+01
    ---------------------------------------------------------------------------
    Average Rosenbluth weight:    3.053288e+01 +/-  1.476984e-01 [-]


    Henry coefficient based on Rosenbluth weight:
    ---------------------------------------------------------------------------
        Block[  0]  4.131328e-06
        Block[  1]  4.116051e-06
        Block[  2]  4.122548e-06
        Block[  3]  4.156846e-06
        Block[  4]  4.122032e-06
    ---------------------------------------------------------------------------
    Average Henry coefficient:    4.129784e-06 +/-  1.997724e-08 [mol/kg/Pa]
    Average Henry coefficient:    3.340427e-05 +/-  1.615884e-07 [molec./uc/Pa]


    Widom insertion chemical potential  statistics:
    ---------------------------------------------------------------------------
        Block[  0] -1206.9709190397484
        Block[  1] -1205.6631608535224
        Block[  2] -1206.2199258243559
        Block[  3] -1209.144574641737
        Block[  4] -1206.1757538659074
    ---------------------------------------------------------------------------
    Excess chemical potential:          -1.206839e+03 +/-  1.704340e+00 [K]
    Tail-correction chemical potential:  0.000000e+00 +/-  0.000000e+00 [K]
    Ideal chemical potential:           -2.614614e+03 +/-  1.534176e+00 [K]
    Total chemical potential:           -3.821453e+03 +/-  1.437305e+00 [K]
    Imposed chemical potential:         -3.810353e+03 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential:          -1.003422e+01 +/-  1.417067e-02 [kJ/mol]
    Tail-correction chemical potential:  0.000000e+00 +/-  0.000000e+00 [kj/mol]
    Ideal chemical potential:           -2.173911e+01 +/-  1.275585e-02 [kJ/mol]
    Total chemical potential:           -3.177333e+01 +/-  1.195043e-02 [kJ/mol]
    Imposed chemical potential:         -3.168105e+01 [kJ/mol]
    ---------------------------------------------------------------------------
    Imposed fugacity:                    1.000000e+05 [Pa]
    Measured fugacity:                   9.690470e+04 +/-  3.944328e+02 [Pa]
```

```
Component 0 [CO2]
    Reinsertion (CBMC)   all:             3859555
    Reinsertion (CBMC)   total:           3859555
    Reinsertion (CBMC)   constructed:     3601038
    Reinsertion (CBMC)   accepted:         817601
    Reinsertion (CBMC)   fraction:       0.211838
    Reinsertion (CBMC)   max-change:     0.000000

    Widom                all:             7714051
    Widom                total:           7714051
    Widom                constructed:     7182081
    Widom                accepted:              0
    Widom                fraction:       0.000000
    Widom                max-change:     0.000000

    Translation          all:             3857673
    Translation          total:           1287453    1284355    1285865
    Translation          constructed:     1287410    1284299    1285858
    Translation          accepted:         644003     641837     645654
    Translation          fraction:       0.500215   0.499735   0.502116
    Translation          max-change:     1.308401   1.300085   1.211848

    Rotation             all:             3859666
    Rotation             total:           1285650    1287774    1286242
    Rotation             constructed:     1285649    1287774    1286239
    Rotation             accepted:         644634     643416     643138
    Rotation             fraction:       0.501407   0.499634   0.500013
    Rotation             max-change:     1.257012   1.248266   1.393558

    Swap (CB/CFCMC)      all:             7712618
    Swap (CB/CFCMC)      total:           1873620    1891529    3947469
    Swap (CB/CFCMC)      constructed:     1797038    1891529    3682717
    Swap (CB/CFCMC)      accepted:         882960     882967    2172416
    Swap (CB/CFCMC)      fraction:       0.471259   0.466801   0.550331
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

