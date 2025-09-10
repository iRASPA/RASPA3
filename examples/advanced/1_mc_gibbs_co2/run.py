import raspalib
import numpy as np

ff = raspalib.ForceField("force_field.json")

mcmoves = raspalib.MCMoveProbabilities(
    translationProbability=0.5,
    rotationProbability=0.5,
    reinsertionCBMCProbability=0.5,
    gibbsSwapCFCMCProbability=1.0,
)

co2 = raspalib.Component(
    componentId=0,
    forceField=ff,
    componentName="CO2",
    fileName="CO2.json",
    particleProbabilities=mcmoves,
)


sysMoves = raspalib.MCMoveProbabilities(gibbsVolumeChangeProbability=0.01)

system0 = raspalib.System(
    systemId=0,
    externalTemperature=240.0,
    forceField=ff,
    components=[co2],
    initialNumberOfMolecules=[256],
    simulationBox=raspalib.SimulationBox(30.0, 30.0, 30.0),
    systemProbabilities=sysMoves,
)

system1 = raspalib.System(
    systemId=1,
    externalTemperature=240.0,
    forceField=ff,
    components=[co2],
    initialNumberOfMolecules=[256],
    simulationBox=raspalib.SimulationBox(30.0, 30.0, 30.0),
    systemProbabilities=sysMoves,
)

mc = raspalib.MonteCarlo(
    numberOfCycles=20000,
    numberOfInitializationCycles=10000,
    numberOfEquilibrationCycles=10000,
    printEvery=1000,
    systems=[system0, system1],
    outputToFiles=True,
)

mc.run()
