import raspalib

ff = raspalib.ForceField("force_field.json")
ff.useCharge = False

mcmoves = raspalib.MCMoveProbabilities(
    translationProbability=1.0,
    rotationProbability=1.0,
    reinsertionCBMCProbability=1.0,
    partialReinsertionCBMCProbability=1.0,
)

biphenyl = raspalib.Component(
    componentId=0,
    forceField=ff,
    componentName="diethyl-biphenyl",
    fileName="diethyl-biphenyl.json",
    particleProbabilities=mcmoves,
)

box = raspalib.SimulationBox(30.0, 30.0, 30.0)

system = raspalib.System(
    systemId=0,
    externalTemperature=298.0,
    forceField=ff,
    components=[biphenyl],
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
