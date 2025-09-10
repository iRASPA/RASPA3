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

# define system
system = raspalib.System(
    systemId=0,
    externalTemperature=500.0,
    forceField=ff,
    components=[co2, n2],
    initialNumberOfMolecules=[50, 50],
    simulationBox=box,
)


mc = raspalib.MonteCarlo(
    numberOfCycles=10000,
    numberOfInitializationCycles=1000,
    systems=[system],
    outputToFiles=True,
)
mc.run()

print(mc.systems[0].writeMCMoveStatistics())
