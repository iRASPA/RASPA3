# Monte Carlo of Methane in a box
\page basic1 Monte Carlo of Methane in a box

Here we run a basic Lennard-Jones simulation, using a United-Atom representation of methane.

You can find the directory with the all input files at `examples/basic/1_mc_methane_in_box`.

As this is likely your first simulation, let me explain what is happening after the input file is read.

First, your `simulation.json` file is read by the InputReader. It parses the input in three sections. First of all, general simulation settings like `SimulationType` and `NumberOfCycles` is read and set to the right attributes. Also, it will read the `force_field.json`. 

Then, from the **list** of systems a System object is created for every entry. The System object is central to all simulations and holds the lists of Atom, Molecule, Thermostat, SimulationBox, MCMoveProbabilitiesParticles, and also a variety of other trackers and property objects. While in this run we don't include a framework, this is also the place where a Framework object can be created.

Next, it will create a Component object. This object holds information about the specific molecule, which holds thermodynamic and structural information. 

Next, it will create a MonteCarlo object, which is then initialized with the System objects and runs the simulation. This object loops over the cycles and performs the moves on the system. 

During your run a few things will be created. First of all, an `output/` directory will be created with a `.txt` file which holds logging statements. Here, you can find all information about the run. Before the run all the settings are logged. Here you can find, for example, all force field settings, information about your system, units, components, frameworks, etc. During the run, with a frequency of `PrintEvery` cycles, statistics about the system are logged. This includes the pressure, energy, number of molecules, etc. At the end all these properties are then wrapped up in reports that give the block averages and report the energy drift of the system. 

## Input file

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
      "TranslationProbability" : 0.5,
      "CreateNumberOfMolecules" : 100
    }
  ]
}
```

## Run with python

Very similarly to what happens behind the scene during a run initiated from a `simulation.json` file, we can also carry out these steps ourselves, using python code. 

We first create a `raspalib.ForceField` object, which without arguments, uses a default force field, which can be found in `data/forcefields/force_field.json`. Then, we have to create the MCMoveProbabilitiesParticles object ourselves, in this case we add just a `translationProbability` of 1.0. Then, we create the component, which in this case is the default methane component given by `raspa.Component`. We also create a SimulationBox. Then we add all these things to `raspa.System`. Afterwards, we pass the System to MonteCarlo and call run, after which the simulation is started. This allows us to interact with the simulation objects after the run.

```python
import raspalib
import numpy as np

ff = raspalib.ForceField("force_field.json")
ff.useCharge = False

mcmoves = raspalib.MCMoveProbabilities(translationProbability=1.0)

methane = raspalib.Component(
    componentId=0,
    forceField=ff,
    componentName="methane",
    fileName="methane.json",
    particleProbabilities=mcmoves,
)

box = raspalib.SimulationBox(30.0, 30.0, 30.0)

system = raspalib.System(
    systemId=0,
    externalTemperature=300.0,
    forceField=ff,
    components=[methane],
    initialNumberOfMolecules=[100],
    simulationBox=box,
    sampleMoviesEvery=10,
)

mc = raspalib.MonteCarlo(
    numberOfCycles=10000,
    numberOfInitializationCycles=1000,
    systems=[system],
    outputToFiles=True,
)
mc.run()

```

