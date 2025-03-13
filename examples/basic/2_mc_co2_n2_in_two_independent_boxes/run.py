import raspalib
import numpy as np

ff = raspalib.ForceField("force_field.json")

mcmoves = raspalib.MCMoveProbabilities(
    translationProbability=1.0, reinsertionCBMCProbability=1.0, rotationProbability=1.0
)

box = raspalib.SimulationBox(30.0, 30.0, 30.0)

# define particles
co2 = raspalib.Component(
    componentId=0,
    forceField=ff,
    componentName="CO2",
    fileName="CO2.json",
    particleProbabilities=mcmoves,
)

n2 = raspalib.Component(
    componentId=1,
    forceField=ff,
    componentName="N2",
    fileName="N2.json",
    particleProbabilities=mcmoves,
)

# define systems
system0 = raspalib.System(
    systemId=0,
    externalTemperature=300.0,
    forceField=ff,
    components=[co2, n2],
    initialNumberOfMolecules=[100, 0],
    simulationBox=box,
)

system1 = raspalib.System(
    systemId=1,
    externalTemperature=300.0,
    forceField=ff,
    components=[co2, n2],
    initialNumberOfMolecules=[0, 100],
    simulationBox=box,
)

mc = raspalib.MonteCarlo(
    numberOfCycles=10000,
    numberOfInitializationCycles=1000,
    systems=[system0, system1],
    outputToFiles=True,
)

mc.run()
