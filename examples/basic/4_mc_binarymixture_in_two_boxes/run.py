import raspalib
import numpy as np

ff = raspalib.ForceField("force_field.json")

mcmoves = raspalib.MCMoveProbabilities(
    translationProbability=1.0, reinsertionCBMCProbability=1.0, rotationProbability=1.0
)
box0 = raspalib.SimulationBox(25.0, 25.0, 25.0)
box1 = raspalib.SimulationBox(30.0, 30.0, 30.0)

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
system0 = raspalib.System(
    systemId=0,
    externalTemperature=500.0,
    forceField=ff,
    components=[co2, n2],
    initialNumberOfMolecules=[50, 25],
    simulationBox=box0,
)

system1 = raspalib.System(
    systemId=1,
    externalTemperature=500.0,
    forceField=ff,
    components=[co2, n2],
    initialNumberOfMolecules=[25, 50],
    simulationBox=box1,
)


mc = raspalib.MonteCarlo(
    numberOfCycles=10000,
    numberOfInitializationCycles=1000,
    systems=[system0, system1],
    outputToFiles=True,
)
mc.run()
