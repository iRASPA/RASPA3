# Examples Basic
\page examples_basic Examples Basic


## Table of Contents
1. [Monte Carlo: methane in box](#Example_1)
2. [Monte Carlo: CO<sub>2</sub> and N<sub>2</sub> in two independent boxes](#Example_2)
3. [Monte Carlo: binary mixture CO<sub>2</sub> and N<sub>2</sub> in box](#Example_3)
4. [Monte Carlo: binary mixture propane and butane in box](#Example_4)
5. [Molecular Dynamics: methane in box (msd)](#Example_5)
6. [Monte Carlo: enthalpy of adsorption in MFI at zero loading](#Example_6)
7. [Monte Carlo: Henry coefficient of methane in MFI](#Example_7)
8. [Monte Carlo: adsorption of methane in MFI](#Example_8)
9. [Monte Carlo: adsorption of butane in MFI](#Example_9)
10. [Monte Carlo: adsorption of CO<sub>2</sub> in Cu-BTC](#Example_10)
11. [Monte Carlo: Henry coefficient of methane, CO<sub>2</sub> and N<sub>2</sub> in MFI](#Example_11)
12. [Monte Carlo: radial distribution function of water](#Example_12)
13. [Molecular Dynamics: radial distribution function of water](#Example_13)


#### Monte Carlo: methane in box <a name="Example_1"></a>

A Monte Carlo run of 100 methane molecules in a $30\times30\times30$ &Aring; box at 300K. After 1000 cycles of initialization the production run is started. A movie is written and every 10th configuration is appended to the movie. The movie is stored in `movies/', and can be viewed with iRASPA or VMD.

The inputs for the simulation are specified in a json-file called `simulation.json`:
```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 10000,
  "NumberOfInitializationCycles" : 1000,
  "PrintEvery" : 1000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 300.0,
      "ChargeMethod" : "None",
      "OutputPDBMovie" : true,
      "SampleMovieEvery" : 10
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
There are global settings, but also settings for `Systems` and `Components`. The latter are arrays of sections with options. In this example, we specify one system of type `Box` with box-lengths $30\times30\times30$ &Aring;. We also set the option to make movies to `true` and sample the movie-snapshots every 10 cycles.

In RASPA, the cycle is define as max(20,$N$) steps, where $N$ is the number of molecules in the system. In every cycle, each of the molecules has on average been used for a Monte Carlo move (accepted or rejected). There is a minimum of 20 steps to avoid that low-density systems or not sampled well. The definition of a cycle is less dependent on the system size. The number of Monte Carlo steps is roughly the number of cycles times the average number of molecules.

The forcefield is defined in `force_field.json`
```
{
  "MixingRule" : "Lorentz-Berthelot",
  "TruncationMethod" : "shifted",
  "TailCorrections" : false,
  "CutOffVDW" : 12.0,
  "PseudoAtoms" :
  [
    {
      "name" : "CH4",
      "framework": false,
      "print_to_output" : true,
      "element" : "C",
      "print_as" : "C",
      "mass" : 16.04246,
      "charge" :  0.0,
      "source" : "M. G. Martin et al., J. Chem. Phys. 2001, 114, 7174-7181"
    }
  ],
  "SelfInteractions" :
  [
    {
      "name" : "CH4",
      "type" : "lennard-jones",
      "parameters" : [158.5, 3.72],
      "source" : "M. G. Martin et al., J. Chem. Phys. 2001, 114, 7174-7181."
    }
  ]
}
```
where we defined the types of the atoms, and their parameters. Here we give them as self-interactions and a mixing rule.
We use a cutoff of 12 &Aring; shifted to zero at the cutoff, and omit tail-corrections.

The methane molecule is defined in `methane.json`
```
{
  "CriticalTemperature" : 190.564,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.01142,
  "Type" : "rigid",
  "pseudoAtoms" :
    [
      ["CH4",[0.0, 0.0, 0.0]]
    ]
}
```
which list the critical temperatures and acentric factors, used to convert pressure into fugacity for adsorption simulations, and the types of the atoms, along with their relative positions.


The output is written to the 'output' directory (one file per system), and the temperature and pressure are appended to all output filenames. In the output file, the simulation writes an important check to the file

```
Energy statistics           |    Energy [K] | Recomputed [K] |     Drift [K] |
=============================================================================
Total potential energy      | -1.892240e+04 |  -1.892240e+04 |  3.412879e-10 |
    molecule-molecule VDW   | -1.892240e+04 |  -1.892240e+04 |  3.412879e-10 |
-----------------------------------------------------------------------------
```

In Monte Carlo, only difference in energies are computed. These differences are continuously added to keep track of the current energies (from which average energies etc. are computed). Obviously, the current energy that is kept track off during the simulation should be equal to a full recalculation of the energies. The difference between the two signals an error. If the drift is higher than say $10^{-3}$ or $10^{-4}$ the results of the simulation are in error. This could be due to an error in one of the Monte Carlo moves or because the force field is ``wrong'' (a typical error is when one forgets to define required potentials).

The performance of Monte Carlo moves is monitored. Translation moves are usually scaled to achieve an acceptance rate of 50%.  Here, the move reached its upper limit of 1.5 &Aring; because of the low density of the system.

```
Component 0 [methane]
    Translation all:             1000000
    Translation total:            333611     334052     332337
    Translation constructed:      333340     333852     332065
    Translation accepted:         265415     266313     264230
    Translation fraction:       0.795582   0.797220   0.795066
    Translation max-change:     1.500000   1.500000   1.500000
```

Averages are computed along with an error bar. The error is computed by dividing the simulation in 5 blocks and calculating the standard deviation. The errors in RASPA are computed as the 95% confidence interval.

```
Total energy:
-------------------------------------------------------------------------------
    Block[  0] -1.821186e+04
    Block[  1] -1.813742e+04
    Block[  2] -1.829147e+04
    Block[  3] -1.836916e+04
    Block[  4] -1.833114e+04
    ---------------------------------------------------------------------------
    Average  -1.826821e+04 +/-  1.160851e+02 [K]
```

#### Monte Carlo: CO<sub>2</sub> and N<sub>2</sub> in two independent boxes <a name="Example_2"></a>

RASPA has a build-in structure of being able to simulate several systems at the same time. This has applications in Gibbs-ensembles and (hyper) parallel tempering for example. However, this capability can also be used for independent systems. The first box is $25\times25\times25$ &Aring; with 90 $^\circ$ angles, containing 100 N<sub>2</sub> and 0 CO<sub>2</sub> and molecules and moved around by translation, rotation and reinsertion. The second box is monoclinic and of size $30\times30\times30$ with $\beta=120^\circ,\alpha=\gamma=90^\circ$ containing 0 N<sub>2</sub> and 100 CO<sub>2</sub> molecules. The first system is at 300K, the second at 500K.

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 10000,
  "NumberOfInitializationCycles" : 1000,
  "PrintEvery" : 1000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [25.0, 25.0, 25.0],
      "ExternalTemperature" : 300.0,
      "ChargeMethod" : "Ewald"
    },
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "BoxAngles" : [90.0, 120.0, 90.0],
      "ExternalTemperature" : 500.0,
      "ChargeMethod" : "Ewald"
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
      "CreateNumberOfMolecules" : [100, 0]
    },
    {
      "Name" : "N2",
      "MoleculeDefinition" : "ExampleDefinitions",
      "TranslationProbability" : 1.0,
      "RotationProbability" : 1.0,
      "ReinsertionProbability" : 1.0,
      "CreateNumberOfMolecules" : [0, 100]
    }
  ]
}
```
with the N<sub>2</sub> defined as 
```json
{
  "CriticalTemperature": 126.192,
  "CriticalPressure": 3395800.0,
  "AcentricFactor": 0.0372,
  "Type": "rigid",
  "pseudoAtoms": [
    ["N_n2", [0.0, 0.0, 0.55]],
    ["N_com", [0.0, 0.0, 0.0]],
    ["N_n2", [0.0, 0.0, -0.55]]
  ]
}
```
and CO<sub>2</sub> defined as
```json
{
  "CriticalTemperature": 304.1282,
  "CriticalPressure": 7377300.0,
  "AcentricFactor": 0.22394,
  "Type": "rigid",
  "pseudoAtoms": [
    ["O_co2", [0.0, 0.0, 1.149]],
    ["C_co2", [0.0, 0.0, 0.0]],
    ["O_co2", [0.0, 0.0, -1.149]]
  ]
}
```

The force field is defined as
```json
{
  "MixingRule" : "Lorentz-Berthelot",
  "TruncationMethod" : "shifted",
  "TailCorrections" : false,
  "CutOffVDW" : 12.0,
  "CutOffCoulomb" : "auto",
  "PseudoAtoms" :
  [
    {
      "name" : "C_co2",
      "framework" : false,
      "print_to_output" : true,
      "element" : "C",
      "print_as" : "C",
      "mass" : 12.0,
      "charge" :  0.6512,
      "source" : "A. Garcia-Sanchez et al., J. Phys. Chem. C 2009, 113, 8814-8820"
    },
    {
      "name" : "O_co2", 
      "framework" : false,
      "print_to_output" : true,  
      "element" : "O",
      "print_as" : "O",
      "mass" : 15.9994,
      "charge" : -0.3256,
      "source" : "A. Garcia-Sanchez et al., J. Phys. Chem. C 2009, 113, 8814-8820"
    },
    {
      "name" : "N_n2",
      "framework" : false,
      "print_to_output" : true,
      "element" : "N",
      "print_as" : "N",
      "mass" : 14.00674,
      "charge" : -0.405,
      "source" : "A. Martin-Calvo et al. , Phys. Chem. Chem. Phys. 2011, 13, 11165-11174"
    },
    {
      "name" : "N_com", 
      "framework" : false,
      "print_to_output" : false, 
      "element" : "N",
      "print_as" : "-",
      "mass" : 0.0,
      "charge" :  0.810,
      "source" : "A. Martin-Calvo et al. , Phys. Chem. Chem. Phys. 2011, 13, 11165-11174"
    }
  ],
  "SelfInteractions" : 
  [
    {
      "name" : "O_co2",
      "type" : "lennard-jones",    
      "parameters" : [85.671, 3.017],
      "source" : "A. Garcia-Sanchez et al., J. Phys. Chem. C 2009, 113, 8814-8820"
    },
    {
      "name" : "C_co2",
      "type" : "lennard-jones",    
      "parameters" : [29.933, 2.745],
      "source" : "A. Garcia-Sanchez et al., J. Phys. Chem. C 2009, 113, 8814-8820"
    },
    {
      "name" : "N_n2",
      "type" : "lennard-jones",    
      "parameters" : [38.298, 3.306],
      "source" : "A. Martin-Calvo et al. , Phys. Chem. Chem. Phys. 2011, 13, 11165-11174"
    },
    {
      "name" : "N_com",
      "type" : "none",    
      "parameters" : [0.0, 1.0],
      "source" : "A. Martin-Calvo et al. , Phys. Chem. Chem. Phys. 2011, 13, 11165-11174"
    }
  ]
}
```
The CO<sub>2</sub> has the charges distrubted over the atoms to approximate the experimental quadrupole. For N<sub>2</sub> we can not do the same. However, we can add a "dummy" site in the center of the two nitrogen atoms to mimick the quadrupole. This `N_com` is placed at the center of mass and has no VDW interactions, and only acts as a charge-center.

Based on the self-interactions and the mixing rule, the cross-interactions are computed:
```
C_co2    - C_co2    Lennard-Jones p₀/kʙ:  29.93300 [K], p₁:  2.74500 [Å]
                                  shift:  -0.01715 [K], tailcorrections: false
C_co2    - O_co2    Lennard-Jones p₀/kʙ:  50.63981 [K], p₁:  2.88100 [Å]
                                  shift:  -0.03878 [K], tailcorrections: false
C_co2    - N_n2     Lennard-Jones p₀/kʙ:  33.85815 [K], p₁:  3.02550 [Å]
                                  shift:  -0.03478 [K], tailcorrections: false
C_co2    - N_com    Lennard-Jones p₀/kʙ:   0.00000 [K], p₁:  1.87250 [Å]
                                  shift:   0.00000 [K], tailcorrections: false
O_co2    - O_co2    Lennard-Jones p₀/kʙ:  85.67100 [K], p₁:  3.01700 [Å]
                                  shift:  -0.08653 [K], tailcorrections: false
O_co2    - N_n2     Lennard-Jones p₀/kʙ:  57.28026 [K], p₁:  3.16150 [Å]
                                  shift:  -0.07659 [K], tailcorrections: false
O_co2    - N_com    Lennard-Jones p₀/kʙ:   0.00000 [K], p₁:  2.00850 [Å]
                                  shift:   0.00000 [K], tailcorrections: false
N_n2     - N_n2     Lennard-Jones p₀/kʙ:  38.29800 [K], p₁:  3.30600 [Å]
                                  shift:  -0.06695 [K], tailcorrections: false
N_n2     - N_com    Lennard-Jones p₀/kʙ:   0.00000 [K], p₁:  2.15300 [Å]
                                  shift:   0.00000 [K], tailcorrections: false
N_com    - N_com    Lennard-Jones p₀/kʙ:   0.00000 [K], p₁:  1.00000 [Å]
                                  shift:   0.00000 [K], tailcorrections: false
```

There will be an output-file for each system
```
output/output_300_0.s0.txt
output/output_500_0.s1.txt
```

Note that we specify only relative probabilities of MC particle moves. They will be correctly rescaled as shown in the output-file:
```
Translation-move probability:             0.3333333333333333 [-]
Rotation-move probability:                0.3333333333333333 [-]
Reinsertion (CBMC) probability:           0.3333333333333333 [-]
```
At every MC-step, each move will be randomly selected with 1/3 probability.

#### Monte Carlo: binary mixture CO<sub>2</sub> and N<sub>2</sub> in box<a name="Example_3"></a>

#### Monte Carlo: binary mixture propane and butane in box<a name="Example_4"></a>

A Monte Carlo run of 50 propane and 50 butane molecules in a $30\times30\times30$ &Aring; box. The MC moves are translation, rotation, full reinsertion, and partial reinsertion. After 1000 steps of initialization the production run is started. We run for 20,000 cycles to get some decent statistics.

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 20000,
  "NumberOfInitializationCycles" : 5000,
  "PrintEvery" : 1000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 500.0,
      "ChargeMethod" : "None"
    }
  ],

  "Components" :
  [
    {
      "Name" : "propane",
      "TranslationProbability" : 1.0,
      "RotationProbability" : 1.0,
      "ReinsertionProbability" : 1.0,
      "PartialReinsertionProbability" : 1.0,
      "CreateNumberOfMolecules" : 50
    },
    {
      "Name" : "butane",
      "TranslationProbability" : 1.0,
      "RotationProbability" : 1.0,
      "ReinsertionProbability" : 1.0,
      "PartialReinsertionProbability" : 1.0,
      "CreateNumberOfMolecules" : 50
    }
  ]
}
```

The propane and butane molecules are modeled as flexible united-atom beads.
The intra-molecular force field contains bond and bend terms for `propane.json`
```json
{
  "CriticalTemperature" : 369.825,
  "CriticalPressure" : 4247660.0,
  "AcentricFactor" : 0.1524,
  "Type" : "flexible",
  "pseudoAtoms" :
    [
      ["CH3", [0.0, 0.0, 0.0]],
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH3", [0.0, 0.0, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2]
  ],
  "Bonds" : [
    [["CH3", "CH2"], "FIXED", [1.54]]
  ],
  "Bends" : [
    [["CH3", "CH2", "CH3"], "HARMONIC", [62500.0, 114]]
  ],
  "Partial-reinsertion" : [
    [0, 1],
    [1, 2],
    [0],
    [2]
  ]
}
```
and bond, bend, and torsion terms for `butane.json`
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
    [["CH3", "CH2"], "FIXED", [1.54]],
    [["CH2", "CH2"], "FIXED", [1.54]]
  ],
  "Bends" : [
    [["CH3", "CH2", "CH2"], "HARMONIC", [62500.0, 114]]
  ],
  "Torsions" : [
    [["CH3", "CH2", "CH2", "CH3"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ],
  "Partial-reinsertion" : [
    [0, 1],
    [2, 3],
    [0],
    [3]
  ]
}
```

The TraPPE forcefield is defined as
```json
{
  "MixingRule" : "Lorentz-Berthelot",
  "TruncationMethod" : "truncated",
  "TailCorrections" : true,
  "CutOff" : 12.0,
  "PseudoAtoms" :
  [
    {
      "name" : "CH3",
      "framework" : false,
      "print_to_output" : true,
      "element" : "C",
      "print_as" : "C",
      "mass" : 15.04,
      "charge" :  0.0,
      "source" : "M. G. Martin and J. I. Siepmann, J. Phys. Chem. B 1998, 102(14), 2569–2577."
    },
    {
      "name" : "CH2",
      "framework" : false,
      "print_to_output" : true,
      "element" : "C",
      "print_as" : "C",
      "mass" : 14.03,
      "charge" :  0.0,
      "source" : "M. G. Martin and J. I. Siepmann, J. Phys. Chem. B 1998, 102(14), 2569–2577."
    }
  ],
  "SelfInteractions" :
  [
    {
      "name" : "CH3",
      "type" : "lennard-jones",
      "parameters" : [98.0, 3.75],
      "source" : "M. G. Martin and J. I. Siepmann, J. Phys. Chem. B 1998, 102(14), 2569–2577."
    },
    {
      "name" : "CH2",
      "type" : "lennard-jones",
      "parameters" : [46.0, 3.95],
      "source" : "M. G. Martin and J. I. Siepmann, J. Phys. Chem. B 1998, 102(14), 2569–2577."
    }
  ]
}
```

```
Energy statistics             |  Energy [K]   | Recomputed [K]|  Drift [K]    |
===============================================================================
Total potential energy/kʙ     | -2.963848e+04 | -2.963848e+04 | -1.520919e-08 |
    molecule-molecule VDW/kʙ  | -8.633271e+04 | -8.633271e+04 | -8.505944e-09 |
    Van der Waals (Tail)/kʙ   | -5.017238e+03 | -5.017238e+03 |  0.000000e+00 |
    bond/kʙ                   |  0.000000e+00 |  0.000000e+00 |  0.000000e+00 |
    bend/kʙ                   |  3.391890e+04 |  3.391890e+04 | -7.468954e-09 |
    torsion/kʙ                |  2.779257e+04 |  2.779257e+04 |  7.657100e-10 |
-------------------------------------------------------------------------------

```

The translation and rotation moves leave the internal structure invariant.
```
Component 0 [propane]
    Translation          all:              500633
    Translation          total:            167364     167001     166268
    Translation          constructed:      167355     166992     166263
    Translation          accepted:         103946     104914     105966
    Translation          fraction:       0.621077   0.628224   0.637320
    Translation          max-change:     1.500000   1.500000   1.500000

    Rotation             all:              499525
    Rotation             total:            166142     166594     166789
    Rotation             constructed:      164548     165646     165967
    Rotation             accepted:         100422     109640     110938
    Rotation             fraction:       0.604435   0.658127   0.665140
    Rotation             max-change:     1.500000   1.500000   1.500000

Component 1 [butane]
    Translation          all:              500514
    Translation          total:            167059     166424     167031
    Translation          constructed:      167053     166417     167023
    Translation          accepted:         103464      97379      98048
    Translation          fraction:       0.619326   0.585126   0.587005
    Translation          max-change:     1.500000   1.500000   1.500000

    Rotation             all:              499580
    Rotation             total:            166131     166746     166703
    Rotation             constructed:      163313     164838     164828
    Rotation             accepted:          95163     102118     102226
    Rotation             fraction:       0.572819   0.612416   0.613222
    Rotation             max-change:     1.230832   1.160213   1.172102
```

The reinsertion-move regrows the molecule at a random position with a new internal structure.

```
Component 0 [propane]
    Reinsertion (CBMC)   all:              499698
    Reinsertion (CBMC)   total:            499698
    Reinsertion (CBMC)   constructed:      499642
    Reinsertion (CBMC)   accepted:         268113
    Reinsertion (CBMC)   fraction:       0.536550
    Reinsertion (CBMC)   max-change:     0.000000

Component 1 [butane]
    Reinsertion (CBMC)   all:              498890
    Reinsertion (CBMC)   total:            498890
    Reinsertion (CBMC)   constructed:      498706
    Reinsertion (CBMC)   accepted:         215229
    Reinsertion (CBMC)   fraction:       0.431416
    Reinsertion (CBMC)   max-change:     0.000000
```

The acceptance percentages are here high enough. But for dense systems, the insertion acceptance ratios become too small. In these cases, other moves (like partial-reinsertion or MC/MD hybrid moves) become essential to properly sample the internal structure of molecules.

```
Component 0 [propane]
    Partial reinsertion (CBMC) all:              500700
    Partial reinsertion (CBMC) total:            500700
    Partial reinsertion (CBMC) constructed:      500700
    Partial reinsertion (CBMC) accepted:         392298
    Partial reinsertion (CBMC) fraction:       0.783499
    Partial reinsertion (CBMC) max-change:     0.000000

Component 1 [butane]
    Partial reinsertion (CBMC) all:              500460
    Partial reinsertion (CBMC) total:            500460
    Partial reinsertion (CBMC) constructed:      500436
    Partial reinsertion (CBMC) accepted:         315821
    Partial reinsertion (CBMC) fraction:       0.631061
    Partial reinsertion (CBMC) max-change:     0.000000
```

The average energies of the internal potentials are computed as:
```
    Bend energy/kʙ 0 [propane]
    ---------------------------------------------------------------------------
        Block[  0]  1.240826e+04
        Block[  1]  1.237668e+04
        Block[  2]  1.251196e+04
        Block[  3]  1.243697e+04
        Block[  4]  1.245964e+04
        -----------------------------------------------------------------------
        Average   1.243870e+04 +/-  6.385477e+01 [K]

    Bend energy/kʙ 1 [butane]
    ---------------------------------------------------------------------------
        Block[  0]  2.468257e+04
        Block[  1]  2.482497e+04
        Block[  2]  2.495069e+04
        Block[  3]  2.489606e+04
        Block[  4]  2.459066e+04
        -----------------------------------------------------------------------
        Average   2.478899e+04 +/-  1.857683e+02 [K]

    Torsion energy/kʙ 1 [butane]
    ---------------------------------------------------------------------------
        Block[  0]  2.492220e+04
        Block[  1]  2.463722e+04
        Block[  2]  2.505684e+04
        Block[  3]  2.487121e+04
        Block[  4]  2.463080e+04
        -----------------------------------------------------------------------
        Average   2.482365e+04 +/-  2.308448e+02 [K]
```

In addition to the average energies, we also get information on the pressure:
```
Average pressure tensor:
-------------------------------------------------------------------------------
 1.7178e+02  4.2992e-02 -1.2965e-03 +/- 1.1175e+01 4.4986e+00 2.0788e+00 [bar]
 4.2992e-02  1.7332e+02 -2.1953e+00 +/- 4.4986e+00 5.9796e+00 3.6461e+00 [bar]
-1.2965e-03 -2.1953e+00  1.7461e+02 +/- 2.0788e+00 3.6461e+00 5.2866e+00 [bar]

    Block[  0]  2.556760e+07
    Block[  1]  2.556760e+07
    Block[  2]  2.556760e+07
    Block[  3]  2.556760e+07
    Block[  4]  2.556760e+07
    ---------------------------------------------------------------------------
    Ideal gas pressure   2.556760e+07 +/-  0.000000e+00 [Pa]
                         2.556760e+02 +/-  0.000000e+00 [bar]


    Block[  0] -8.290629e+06
    Block[  1] -8.220314e+06
    Block[  2] -8.465780e+06
    Block[  3] -8.732872e+06
    Block[  4] -7.510238e+06
    ---------------------------------------------------------------------------
    Excess pressure  -8.243966e+06 +/-  5.652840e+05 [Pa]
                     -8.243966e+01 +/-  5.652840e+00 [bar]


    Block[  0]  1.727697e+07
    Block[  1]  1.734729e+07
    Block[  2]  1.710182e+07
    Block[  3]  1.683473e+07
    Block[  4]  1.805736e+07
    ---------------------------------------------------------------------------
    Pressure average   1.732363e+07 +/-  5.652840e+05 [Pa]
                       1.732363e+02 +/-  5.652840e+00 [bar]
```


#### Molecular Dynamics: methane in box (msd)<a name="Example_5"></a>

#### Monte Carlo: enthalpy of adsorption of methane in MFI at zero loading<a name="Example_6"></a>

#### Monte Carlo: Henry coefficient of methane in MFI<a name="Example_7"></a>

#### Monte Carlo: adsorption of methane in MFI<a name="Example_8"></a>

#### Monte Carlo: adsorption of butane in MFI<a name="Example_9"></a>

#### Monte Carlo: adsorption of CO<sub>2</sub> in Cu-BTC<a name="Example_10"></a>

#### Monte Carlo: Henry coefficient of methane, CO<sub>2</sub> and N<sub>2</sub> in MFI<a name="Example_11"></a>

#### Monte Carlo: radial distribution function of water<a name="Example_12"></a>

#### Molecular Dynamics: radial distribution function of water<a name="Example_13"></a>
