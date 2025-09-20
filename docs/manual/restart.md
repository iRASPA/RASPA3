# Restart-files `RASPA`
\page restart Restart Files

## Writing restart-files

Restart-files are automatically written (and updated every 5000 cycles). Consider for example the following Gibbs-ensemble simulation
```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 25000,
  "NumberOfInitializationCycles" : 10000,
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
At the end of the simulation you will find in the `output`-directory the files
```
output/restart_240_0.s0.json
output/restart_240_0.s1.json
```


----------------------------------------------------------------------------------

## Using restart-files

To use these restart-files to start from the saved (and hopefully equilibrated) postions, 
first copy these files to the working directory (the location of the `simulation.json`).
```
cp output\restart_240_0.s0.json .
cp output\restart_240_0.s1.json .
```
This is because the restart-files in the output-directory are overwritten every run.
Next change the `simulation.json` to use the restart-files

The contents of the restart files look like
```
{
  "CO2": [
    [ 16.640068355882505, 9.764919807611305, 4.333198721377549 ],
    [ 17.15733718188185, 8.82872907970449, 4.75293778631024 ],
    [ 17.67460600788119, 7.892538351797674, 5.172676851242931 ],
    [ 3.815549021105568, 8.879311693833092, 17.616704775297364 ],
    [ 4.677375964586537, 8.549547910646638, 18.30132961182377 ],
    [ 5.539202908067507, 8.219784127460183, 18.985954448350178 ]
  ],
  "SimulationBox": {
    "angle-alpha": 90.0,
    "angle-beta": 90.0,
    "angle-gamma": 90.0,
    "length-a": 31.881622110825266,
    "length-b": 31.881622110825266,
    "length-c": 31.881622110825266
  }
}

```
So, the component-name as the key, and the value is an array of positions,
and information on the simulation-box (during NPT or Gibbs the box changes).
CO<sub>2</sub> molecules are usually modeled as `rigid` and implemented using quaternions.
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
The quaternions are computed from the positions using singular-value decompositions.
The advantage is that, even when the positions are input by hand and do not satisfy
the rigid constraints, an optimal mapping or best fit is obtained.
For output restart-files, reading these in, should result in the exact same
positions and energies.



```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfCycles" : 25000,
  "NumberOfInitializationCycles" : 1000,
  "PrintEvery" : 1000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 240.0,
      "ChargeMethod" : "Ewald",
      "GibbsVolumeMoveProbability" : 0.01,
      "RestartFileName" : "restart_240_0.s0"
    },
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 240.0,
      "ChargeMethod" : "Ewald",
      "GibbsVolumeMoveProbability" : 0.01,
      "RestartFileName" : "restart_240_0.s1"
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
      "CreateNumberOfMolecules" : [0, 0]
    }
  ]
}
```
Note because the positions are already equilibrated in this case, you can lower the amount of initialization cycles.
This is the main point of using restart-file. But if you use the files to run at different conditions (for example
a different temperature) you would still need a significant amount of initialization cycles.

The number of `CreateNumberOfMolecules` should be put to zero, since the molecules are read from the restart-files.
However, you still have the freedom to create additional molecules on top of what was read from the restart-files.
This might be usefull to construct very high density systems.

----------------------------------------------------------------------------------

## Using binary-restart-files

Binary restart-files are usefull in case of system crashes or limits on running-times in computer clusters.
By continuing from the binary restart-file you would get identical results are running the long job, i.e.
it contains the entire state of the simulation. Note that therefore whatever you change in the `simulation.json` 
is ignored since it uses the state from the binary restart-file.

During the simulation a file named `restart_data.bin` is written out.
The option to continue from this file is
```json
{
  "RestartFromBinaryFile" : "true"
}
```
Note not everything can be recovered, like written pdb-movies.
The setting
```json
{
  "WriteBinaryRestartEvery" : "5000"
}
```
controls how often the file is written out. The default is every 5000 cycles.

NOTE: the format of the binary-restart files will keep changing until version RASPA 3.1. 
That means that you probably cannot use restart-files from a different version.
After the release of RASPA 3.1 we will start the process of making them backwards-compatible.
