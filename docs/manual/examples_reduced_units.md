# Examples Basic
\page examples_reduced_units Examples Reduced Units


## Table of Contents
1. [Monte Carlo: particle in box (NVT)](#Example_reduced_units_1)
2. [Monte Carlo: particle in box (CFCMC)](#Example_reduced_units_2)
3. [Monte Carlo: Gibbs Vapor-Liquid Equilibrium](#Example_reduced_units_3)
4. [Monte Carlo: Gibbs Vapor-Liquid Equilibrium (CFCMC)](#Example_reduced_units_4)


#### Monte Carlo: particles in box (NVT) <a name="Example_reduced_units_1"></a>

#### Monte Carlo: particles in a box (CFCMC)<a name="Example_reduced_units_2"></a>

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
| Density \f$\left\langle\rho\right\rangle\f$ | \f$m\cdot\sigma^{-3}\f$  | \f$0.001928 \pm 0.000059\f$ | \f$0.838342 \pm 0.000317\f$   |
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
| Density \f$\left\langle\rho\right\rangle\f$ | \f$m\cdot\sigma^{-3}\f$ | \f$0.001946 \pm  0.000058\f$   | \f$ 0.836678 \pm 0.000218\f$ |
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
