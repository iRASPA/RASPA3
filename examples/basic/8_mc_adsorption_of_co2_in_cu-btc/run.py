import raspalib

ff = raspalib.ForceField("force_field.json")

mcmoves = raspalib.MCMoveProbabilities(
    translationProbability=0.5,
    rotationProbability=0.5,
    reinsertionCBMCProbability=0.5,
    swapProbability=1.0,
    widomProbability=1.0,
)

co2 = raspalib.Component(
    componentId=0,
    forceField=ff,
    componentName="CO2",
    fileName="CO2.json",
    particleProbabilities=mcmoves,
)

cubtc = raspalib.Framework(
    frameworkId=0,
    forceField=ff,
    componentName="Cu-BTC",
    fileName="Cu-BTC.cif",
    numberOfUnitCells=raspalib.int3(1, 1, 1),
)

system = raspalib.System(
    systemId=0,
    externalTemperature=300.0,
    externalPressure=1e4,
    forceField=ff,
    components=[co2],
    initialNumberOfMolecules=[0],
    frameworkComponents=[cubtc],
)

mc = raspalib.MonteCarlo(numberOfCycles=10000, numberOfInitializationCycles=2000, systems=[system], outputToFiles=True)

mc.run()
