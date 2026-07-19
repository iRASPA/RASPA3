import raspalib

ff = raspalib.ForceField("force_field.json")
ff.useCharge = False

mcmoves = raspalib.MCMoveProbabilities(
    translationProbability=1.0,
    rotationProbability=1.0,
    reinsertionCBMCProbability=1.0,
)

cyclohexane = raspalib.Component(
    componentId=0,
    forceField=ff,
    componentName="cyclohexane",
    fileName="cyclohexane.json",
    particleProbabilities=mcmoves,
)

box = raspalib.SimulationBox(30.0, 30.0, 30.0)

system = raspalib.System(
    systemId=0,
    externalTemperature=300.0,
    forceField=ff,
    components=[cyclohexane],
    initialNumberOfMolecules=[30],
    simulationBox=box,
)

mc = raspalib.MonteCarlo(
    numberOfProductionCycles=10000,
    numberOfInitializationCycles=1000,
    printEvery=1000,
    systems=[system],
    outputToFiles=True,
)
mc.run()
