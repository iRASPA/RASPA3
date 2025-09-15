# Examples advanced
\page examples_advanced Examples Advanced


## Table of Contents
1. [Monte Carlo: CFCMC Gibbs CO₂](#Example_advanced_1)
2. [Monte Carlo: CFCMC adsorption of CO₂ in MFI](#Example_advanced_2)
3. [Monte Carlo: CFCMC binary mixture adsorption of CO₂ and N₂ in DMOF](#Example_advanced_3)
4. [Transition Matrix Monte Carlo: methane in Tobacco-667](#Example_advanced_4)
5. [HMC CFCMC co2 in MFI](#Example_advanced_5)


#### Monte Carlo: CFCMC Gibbs CO₂<a name="Example_advanced_1"></a>

Status: Error in CFCMC Reinsertion move.

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
      "ReinsertionProbability" : 0.0,
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

Status: Error in CFCMC Reinsertion move.

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
Loadings
===============================================================================

Component 0 (CO2)
    Block[  0]  2.268176e+01
    Block[  1]  2.390155e+01
    Block[  2]  2.417783e+01
    Block[  3]  2.374630e+01
    Block[  4]  2.352465e+01
    ---------------------------------------------------------------------------
    Abs. loading average   2.360733e+01 +/-  7.063445e-01 [molecules/cell]
    Abs. loading average   2.950916e+00 +/-  8.829307e-02 [molecules/uc]
    Abs. loading average   3.648230e-01 +/-  1.091571e-02 [mol/kg-framework]
    Abs. loading average   1.605177e+01 +/-  4.802781e-01 [mg/g-framework]
```

```
    Lambda histogram and bias:
    ---------------------------------------------------------------------------
     0-0.000000 (lambda) P:  1.04427e+00 +/- 8.26754e-02 bias:  3.38016e+00 [-]
     1-0.025000 (lambda) P:  1.06477e+00 +/- 9.35717e-02 bias:  3.42941e+00 [-]
     2-0.050000 (lambda) P:  1.03771e+00 +/- 8.96983e-02 bias:  3.37832e+00 [-]
     3-0.075000 (lambda) P:  1.07133e+00 +/- 2.69698e-02 bias:  3.45141e+00 [-]
     4-0.100000 (lambda) P:  9.88510e-01 +/- 8.36838e-02 bias:  3.51119e+00 [-]
     5-0.125000 (lambda) P:  1.09060e+00 +/- 1.03676e-01 bias:  3.83504e+00 [-]
     6-0.150000 (lambda) P:  1.05124e+00 +/- 7.33517e-02 bias:  4.26115e+00 [-]
     7-0.175000 (lambda) P:  1.01803e+00 +/- 8.52975e-02 bias:  4.40352e+00 [-]
     8-0.200000 (lambda) P:  1.11028e+00 +/- 5.46020e-02 bias:  4.37732e+00 [-]
     9-0.225000 (lambda) P:  9.91380e-01 +/- 8.70265e-02 bias:  4.06430e+00 [-]
    10-0.250000 (lambda) P:  9.73340e-01 +/- 5.61460e-02 bias:  3.68947e+00 [-]
    11-0.275000 (lambda) P:  1.04550e+00 +/- 5.35056e-02 bias:  3.27412e+00 [-]
    12-0.300000 (lambda) P:  9.70060e-01 +/- 6.09049e-02 bias:  2.77873e+00 [-]
    13-0.325000 (lambda) P:  1.02377e+00 +/- 4.71204e-02 bias:  2.36098e+00 [-]
    14-0.350000 (lambda) P:  1.00614e+00 +/- 3.34936e-02 bias:  1.95043e+00 [-]
    15-0.375000 (lambda) P:  9.60220e-01 +/- 6.60378e-02 bias:  1.60172e+00 [-]
    16-0.400000 (lambda) P:  9.61860e-01 +/- 5.75418e-02 bias:  1.36879e+00 [-]
    17-0.425000 (lambda) P:  9.65140e-01 +/- 4.76942e-02 bias:  1.14223e+00 [-]
    18-0.450000 (lambda) P:  9.99580e-01 +/- 4.76330e-02 bias:  1.10270e+00 [-]
    19-0.475000 (lambda) P:  1.02787e+00 +/- 2.49098e-02 bias:  1.02848e+00 [-]
    20-0.500000 (lambda) P:  1.00450e+00 +/- 8.25656e-02 bias:  1.01609e+00 [-]
    21-0.525000 (lambda) P:  1.02459e+00 +/- 7.34400e-02 bias:  1.02617e+00 [-]
    22-0.550000 (lambda) P:  1.02992e+00 +/- 4.85756e-02 bias:  1.08035e+00 [-]
    23-0.575000 (lambda) P:  1.00819e+00 +/- 7.33076e-02 bias:  9.91621e-01 [-]
    24-0.600000 (lambda) P:  1.05780e+00 +/- 8.34629e-02 bias:  1.02906e+00 [-]
    25-0.625000 (lambda) P:  9.84820e-01 +/- 6.81377e-02 bias:  9.09668e-01 [-]
    26-0.650000 (lambda) P:  9.86050e-01 +/- 4.28514e-02 bias:  9.46465e-01 [-]
    27-0.675000 (lambda) P:  9.77030e-01 +/- 6.62581e-02 bias:  9.18867e-01 [-]
    28-0.700000 (lambda) P:  9.95890e-01 +/- 8.55023e-02 bias:  8.23457e-01 [-]
    29-0.725000 (lambda) P:  1.04099e+00 +/- 5.17018e-02 bias:  8.19375e-01 [-]
    30-0.750000 (lambda) P:  9.93020e-01 +/- 4.65673e-02 bias:  6.67949e-01 [-]
    31-0.775000 (lambda) P:  9.08560e-01 +/- 4.12495e-02 bias:  5.09863e-01 [-]
    32-0.800000 (lambda) P:  9.65960e-01 +/- 5.95063e-02 bias:  4.79512e-01 [-]
    33-0.825000 (lambda) P:  9.17170e-01 +/- 2.08006e-02 bias:  3.04355e-01 [-]
    34-0.850000 (lambda) P:  9.66780e-01 +/- 5.23925e-02 bias:  2.45059e-01 [-]
    35-0.875000 (lambda) P:  9.57760e-01 +/- 6.75122e-02 bias:  1.49746e-01 [-]
    36-0.900000 (lambda) P:  9.91380e-01 +/- 4.35486e-02 bias:  1.32578e-01 [-]
    37-0.925000 (lambda) P:  9.22500e-01 +/- 6.63900e-02 bias:  3.83203e-02 [-]
    38-0.950000 (lambda) P:  9.70060e-01 +/- 6.75601e-02 bias:  7.07227e-02 [-]
    39-0.975000 (lambda) P:  9.56530e-01 +/- 5.78841e-02 bias:  4.11133e-02 [-]
    40-1.000000 (lambda) P:  9.38900e-01 +/- 3.42868e-02 bias:  0.00000e+00 [-]


    Lambda statistics:
    ---------------------------------------------------------------------------
     0-0.000000 (lambda) Free energy: 2.163686e+01 +/- 2.812850e+01 [K]
     1-0.025000 (lambda) Free energy: 1.477427e+01 +/- 3.088778e+01 [K]
     2-0.050000 (lambda) Free energy: 2.386136e+01 +/- 2.981821e+01 [K]
     3-0.075000 (lambda) Free energy: 1.260612e+01 +/- 8.869530e+00 [K]
     4-0.100000 (lambda) Free energy: 4.100761e+01 +/- 3.020381e+01 [K]
     5-0.125000 (lambda) Free energy: 6.313140e+00 +/- 3.387329e+01 [K]
     6-0.150000 (lambda) Free energy: 1.928858e+01 +/- 2.435124e+01 [K]
     7-0.175000 (lambda) Free energy: 3.062025e+01 +/- 2.967910e+01 [K]
     8-0.200000 (lambda) Free energy: -6.385932e-14 +/- 1.752870e+01 [K]
     9-0.225000 (lambda) Free energy: 3.998421e+01 +/- 3.110420e+01 [K]
    10-0.250000 (lambda) Free energy: 4.646686e+01 +/- 2.035710e+01 [K]
    11-0.275000 (lambda) Free energy: 2.122132e+01 +/- 1.829352e+01 [K]
    12-0.300000 (lambda) Free energy: 4.765843e+01 +/- 2.163498e+01 [K]
    13-0.325000 (lambda) Free energy: 2.863551e+01 +/- 1.611392e+01 [K]
    14-0.350000 (lambda) Free energy: 3.476736e+01 +/- 1.160215e+01 [K]
    15-0.375000 (lambda) Free energy: 5.125744e+01 +/- 2.464157e+01 [K]
    16-0.400000 (lambda) Free energy: 5.065505e+01 +/- 2.112221e+01 [K]
    17-0.425000 (lambda) Free energy: 4.945335e+01 +/- 1.766779e+01 [K]
    18-0.450000 (lambda) Free energy: 3.707644e+01 +/- 1.675684e+01 [K]
    19-0.475000 (lambda) Free energy: 2.722463e+01 +/- 8.553036e+00 [K]
    20-0.500000 (lambda) Free energy: 3.534321e+01 +/- 2.876342e+01 [K]
    21-0.525000 (lambda) Free energy: 2.835288e+01 +/- 2.524720e+01 [K]
    22-0.550000 (lambda) Free energy: 2.652130e+01 +/- 1.666448e+01 [K]
    23-0.575000 (lambda) Free energy: 3.404885e+01 +/- 2.535401e+01 [K]
    24-0.600000 (lambda) Free energy: 1.709261e+01 +/- 2.782266e+01 [K]
    25-0.625000 (lambda) Free energy: 4.232778e+01 +/- 2.436669e+01 [K]
    26-0.650000 (lambda) Free energy: 4.188718e+01 +/- 1.538569e+01 [K]
    27-0.675000 (lambda) Free energy: 4.513115e+01 +/- 2.388537e+01 [K]
    28-0.700000 (lambda) Free energy: 3.838197e+01 +/- 2.979631e+01 [K]
    29-0.725000 (lambda) Free energy: 2.274736e+01 +/- 1.763891e+01 [K]
    30-0.750000 (lambda) Free energy: 3.940073e+01 +/- 1.649566e+01 [K]
    31-0.775000 (lambda) Free energy: 7.077889e+01 +/- 1.595799e+01 [K]
    32-0.800000 (lambda) Free energy: 4.915356e+01 +/- 2.149868e+01 [K]
    33-0.825000 (lambda) Free energy: 6.744942e+01 +/- 8.017267e+00 [K]
    34-0.850000 (lambda) Free energy: 4.885403e+01 +/- 1.915263e+01 [K]
    35-0.875000 (lambda) Free energy: 5.216296e+01 +/- 2.528374e+01 [K]
    36-0.900000 (lambda) Free energy: 3.998421e+01 +/- 1.549987e+01 [K]
    37-0.925000 (lambda) Free energy: 6.540395e+01 +/- 2.456294e+01 [K]
    38-0.950000 (lambda) Free energy: 4.765843e+01 +/- 2.475185e+01 [K]
    39-0.975000 (lambda) Free energy: 5.261659e+01 +/- 2.104246e+01 [K]
    40-1.000000 (lambda) Free energy: 5.918352e+01 +/- 1.298501e+01 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential: (ln(P(lambda=1))-ln(P(lambda=0)))/Beta
        Block[  0] -1144.4772539075184
        Block[  1] -1160.5884039951882
        Block[  2] -1188.6413089847767
        Block[  3] -1155.7622216340555
        Block[  4] -1131.0563727485423
    ---------------------------------------------------------------------------
    Excess chemical potential:    -1.155650e+03 +/-  2.664458e+01 [K]
    Ideal gas chemical potential: -2.647280e+03 +/-  1.066744e+01 [K]
    Total chemical potential:     -3.802930e+03 +/-  2.115599e+01 [K]
    Imposed chemical potential:   -3.810353e+03 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential:    -9.608608e+00 +/-  2.215354e-01 [kJ/mol]
    Ideal gas chemical potential: -2.201072e+01 +/-  8.869402e-02 [kJ/mol]
    Total chemical potential:     -3.161933e+01 +/-  1.759007e-01 [kJ/mol]
    Imposed chemical potential:   -3.168105e+01 [kJ/mol]
    ---------------------------------------------------------------------------
    Imposed fugacity:              1.000000e+05 [Pa]
    Measured fugacity:             1.021253e+05 +/-  6.113034e+03 [Pa]


    Thermodynamic integration (dU/dlambda)
    ===========================================================================

     0-0.000000 (lambda) <dU/dlambda>:  0.000000e+00 +/- 0.000000e+00 [K/-]
     1-0.025000 (lambda) <dU/dlambda>:  5.740556e+01 +/- 7.104575e+00 [K/-]
     2-0.050000 (lambda) <dU/dlambda>:  2.796506e+02 +/- 2.425887e+01 [K/-]
     3-0.075000 (lambda) <dU/dlambda>:  1.026811e+03 +/- 1.075461e+02 [K/-]
     4-0.100000 (lambda) <dU/dlambda>:  2.584142e+03 +/- 2.258609e+02 [K/-]
     5-0.125000 (lambda) <dU/dlambda>:  5.133233e+03 +/- 2.290718e+02 [K/-]
     6-0.150000 (lambda) <dU/dlambda>:  4.952365e+03 +/- 7.724317e+02 [K/-]
     7-0.175000 (lambda) <dU/dlambda>:  8.959663e+02 +/- 5.022221e+02 [K/-]
     8-0.200000 (lambda) <dU/dlambda>: -2.750024e+03 +/- 4.481720e+02 [K/-]
     9-0.225000 (lambda) <dU/dlambda>: -4.583193e+03 +/- 5.154787e+02 [K/-]
    10-0.250000 (lambda) <dU/dlambda>: -5.680666e+03 +/- 4.616630e+02 [K/-]
    11-0.275000 (lambda) <dU/dlambda>: -6.438317e+03 +/- 2.372826e+02 [K/-]
    12-0.300000 (lambda) <dU/dlambda>: -6.408192e+03 +/- 2.326307e+02 [K/-]
    13-0.325000 (lambda) <dU/dlambda>: -5.959944e+03 +/- 1.216547e+02 [K/-]
    14-0.350000 (lambda) <dU/dlambda>: -5.088885e+03 +/- 1.991323e+01 [K/-]
    15-0.375000 (lambda) <dU/dlambda>: -4.050472e+03 +/- 3.046127e+01 [K/-]
    16-0.400000 (lambda) <dU/dlambda>: -3.018268e+03 +/- 3.398578e+01 [K/-]
    17-0.425000 (lambda) <dU/dlambda>: -2.033371e+03 +/- 1.364816e+01 [K/-]
    18-0.450000 (lambda) <dU/dlambda>: -1.202849e+03 +/- 1.131192e+01 [K/-]
    19-0.475000 (lambda) <dU/dlambda>: -5.258241e+02 +/- 3.187243e+00 [K/-]
    20-0.500000 (lambda) <dU/dlambda>:  0.000000e+00 +/- 0.000000e+00 [K/-]
    21-0.525000 (lambda) <dU/dlambda>: -2.114840e+01 +/- 5.182677e+00 [K/-]
    22-0.550000 (lambda) <dU/dlambda>: -4.126118e+01 +/- 1.612878e+01 [K/-]
    23-0.575000 (lambda) <dU/dlambda>: -8.550023e+01 +/- 3.501469e+01 [K/-]
    24-0.600000 (lambda) <dU/dlambda>: -1.471067e+02 +/- 3.829334e+01 [K/-]
    25-0.625000 (lambda) <dU/dlambda>: -2.823943e+02 +/- 8.202454e+01 [K/-]
    26-0.650000 (lambda) <dU/dlambda>: -4.324543e+02 +/- 5.673625e+01 [K/-]
    27-0.675000 (lambda) <dU/dlambda>: -6.132271e+02 +/- 5.695995e+01 [K/-]
    28-0.700000 (lambda) <dU/dlambda>: -9.087684e+02 +/- 7.407049e+01 [K/-]
    29-0.725000 (lambda) <dU/dlambda>: -1.154457e+03 +/- 5.846130e+01 [K/-]
    30-0.750000 (lambda) <dU/dlambda>: -1.444578e+03 +/- 1.116919e+02 [K/-]
    31-0.775000 (lambda) <dU/dlambda>: -1.568963e+03 +/- 2.073155e+01 [K/-]
    32-0.800000 (lambda) <dU/dlambda>: -1.555382e+03 +/- 2.558521e+01 [K/-]
    33-0.825000 (lambda) <dU/dlambda>: -1.474193e+03 +/- 4.345040e+01 [K/-]
    34-0.850000 (lambda) <dU/dlambda>: -1.215539e+03 +/- 3.960011e+01 [K/-]
    35-0.875000 (lambda) <dU/dlambda>: -1.022996e+03 +/- 3.812844e+01 [K/-]
    36-0.900000 (lambda) <dU/dlambda>: -7.606471e+02 +/- 6.576672e+01 [K/-]
    37-0.925000 (lambda) <dU/dlambda>: -5.169291e+02 +/- 1.941765e+01 [K/-]
    38-0.950000 (lambda) <dU/dlambda>: -3.029754e+02 +/- 1.015839e+01 [K/-]
    39-0.975000 (lambda) <dU/dlambda>: -1.276039e+02 +/- 6.616096e+00 [K/-]
    40-1.000000 (lambda) <dU/dlambda>:  0.000000e+00 +/- 0.000000e+00 [K/-]
    ---------------------------------------------------------------------------
    Excess chemical potential: integral du/dlambda over lambda (Simpson's rule)
        Block[  0] -1147.6671622001415
        Block[  1] -1175.391099224101
        Block[  2] -1152.7591364848136
        Block[  3] -1162.2101184503545
        Block[  4] -1186.090657601686
    ---------------------------------------------------------------------------
    Excess chemical potential:   -1.163861e+03 +/-  1.977904e+01 [K]
    Ideal chemical potential:    -2.647280e+03 +/-  1.066744e+01 [K]
    Total chemical potential:    -3.811141e+03 +/-  1.964511e+01 [K]
    Imposed chemical potential:  -3.810353e+03 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential:   -9.676883e+00 +/-  1.644521e-01 [kJ/mol]
    Ideal chemical potential:    -2.201072e+01 +/-  8.869402e-02 [kJ/mol]
    Total chemical potential:    -3.168760e+01 +/-  1.633385e-01 [kJ/mol]
    Imposed chemical potential:  -3.168105e+01 [kJ/mol]
    ---------------------------------------------------------------------------
    Imposed fugacity:             1.000000e+05 [Pa]
    Measured fugacity:            9.977704e+04 +/-  5.513037e+03 [Pa]


    Widom insertion Rosenbluth weight statistics:
    ---------------------------------------------------------------------------
        Block[  0]  3.142435e+01
        Block[  1]  3.105714e+01
        Block[  2]  3.100723e+01
        Block[  3]  3.172060e+01
        Block[  4]  3.071428e+01
    ---------------------------------------------------------------------------
    Average Rosenbluth weight:    3.118273e+01 +/-  4.862780e-01 [-]


    Henry coefficient based on Rosenbluth weight:
    ---------------------------------------------------------------------------
        Block[  0]  4.250360e-06
        Block[  1]  4.200693e-06
        Block[  2]  4.193942e-06
        Block[  3]  4.290431e-06
        Block[  4]  4.154319e-06
    ---------------------------------------------------------------------------
    Average Henry coefficient:    4.217680e-06 +/-  6.577247e-08 [mol/kg/Pa]
    Average Henry coefficient:    3.411523e-05 +/-  5.320087e-07 [molec./uc/Pa]


    Widom insertion chemical potential  statistics:
    ---------------------------------------------------------------------------
        Block[  0] -1216.9978489634334
        Block[  1] -1212.8485787014306
        Block[  2] -1212.2808520460585
        Block[  3] -1220.3102373694749
        Block[  4] -1208.9299580646082
    ---------------------------------------------------------------------------
    Excess chemical potential:          -1.214273e+03 +/-  5.496857e+00 [K]
    Tail-correction chemical potential:  0.000000e+00 +/-  0.000000e+00 [K]
    Ideal chemical potential:           -2.647280e+03 +/-  1.066744e+01 [K]
    Total chemical potential:           -3.861553e+03 +/-  1.329125e+01 [K]
    Imposed chemical potential:         -3.810353e+03 [K]
    ---------------------------------------------------------------------------
    Excess chemical potential:          -1.009603e+01 +/-  4.570343e-02 [kJ/mol]
    Tail-correction chemical potential:  0.000000e+00 +/-  0.000000e+00 [kj/mol]
    Ideal chemical potential:           -2.201072e+01 +/-  8.869402e-02 [kJ/mol]
    Total chemical potential:           -3.210675e+01 +/-  1.105096e-01 [kJ/mol]
    Imposed chemical potential:         -3.168105e+01 [kJ/mol]
    ---------------------------------------------------------------------------
    Imposed fugacity:                    1.000000e+05 [Pa]
    Measured fugacity:                   8.649857e+04 +/-  3.223570e+03 [Pa]

```

```
Component 0 [CO2]
    Reinsertion (CBMC)   all:              355485
    Reinsertion (CBMC)   total:            355485
    Reinsertion (CBMC)   constructed:      332483
    Reinsertion (CBMC)   accepted:          79762
    Reinsertion (CBMC)   fraction:       0.224375
    Reinsertion (CBMC)   max-change:     0.000000

    Widom                all:              713540
    Widom                total:            713540
    Widom                constructed:      666233
    Widom                accepted:              0
    Widom                fraction:       0.000000
    Widom                max-change:     0.000000

    Translation          all:              355921
    Translation          total:            118702     118591     118628
    Translation          constructed:      118699     118587     118622
    Translation          accepted:          54751      53992      53886
    Translation          fraction:       0.461247   0.455279   0.454244
    Translation          max-change:     0.015069   0.018671   0.020785

    Rotation             all:              356118
    Rotation             total:            118633     118461     119024
    Rotation             constructed:      118629     118456     119021
    Rotation             accepted:          53837      53816      54012
    Rotation             fraction:       0.453811   0.454293   0.453791
    Rotation             max-change:     0.021357   0.020527   0.024613

    Swap (CB/CFCMC)      all:              712582
    Swap (CB/CFCMC)      total:            164445     167444     380693
    Swap (CB/CFCMC)      constructed:      158330     167444     356236
    Swap (CB/CFCMC)      accepted:          80079      80082     211619
    Swap (CB/CFCMC)      fraction:       0.486965   0.478261   0.555878
    Swap (CB/CFCMC)      max-change:     0.125781   0.113082   1.000000
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

