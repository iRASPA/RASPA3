import raspalib
import numpy as np

ff = raspalib.ForceField("force_field.json")

mcmoves = raspalib.MCMoveProbabilities(translationProbability=0.5, rotationProbability=0.5, reinsertionCBMCProbability=1.0)

water = raspalib.Component(
    componentId=0,
    forceField=ff,
    componentName="water",
    fileName="water.json",
    particleProbabilities=mcmoves,
)

box = raspalib.SimulationBox(24.83, 24.83, 24.83)

system = raspalib.System(
    systemId=0,
    externalTemperature=300.0,
    forceField=ff,
    components=[water],
    initialNumberOfMolecules=[512],
    simulationBox=box,
    sampleMoviesEvery=10,
)

mc = raspalib.MonteCarlo(
    numberOfCycles=10000,
    numberOfInitializationCycles=1000,
    printEvery=100,
    systems=[system],
    outputToFiles=True,
)
mc.run()
