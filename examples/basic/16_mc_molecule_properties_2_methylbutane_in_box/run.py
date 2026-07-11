import raspalib

ff = raspalib.ForceField("force_field.json")
ff.useCharge = False

mcmoves = raspalib.MCMoveProbabilities(
    translationProbability=1.0,
    rotationProbability=1.0,
    reinsertionCBMCProbability=1.0,
    partialReinsertionCBMCProbability=1.0,
)

methylbutane = raspalib.Component(
    componentId=0,
    forceField=ff,
    componentName="2-methylbutane",
    fileName="2-methylbutane.json",
    particleProbabilities=mcmoves,
)

box = raspalib.SimulationBox(25.0, 25.0, 25.0)

system = raspalib.System(
    systemId=0,
    externalTemperature=298.0,
    forceField=ff,
    components=[methylbutane],
    initialNumberOfMolecules=[32],
    simulationBox=box,
)

mc = raspalib.MonteCarlo(
    numberOfProductionCycles=1000000,
    numberOfInitializationCycles=10000,
    printEvery=50000,
    systems=[system],
    outputToFiles=True,
)
mc.run()
