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
    numberOfProductionCycles=10000,
    numberOfInitializationCycles=1000,
    systems=[system],
    outputToFiles=True,
)
mc.run()
