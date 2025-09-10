import raspalib

ff = raspalib.ForceField("force_field.json")

mcmoves = raspalib.MCMoveProbabilities(
    translationProbability=0.5,
    rotationProbability=0.5,
    reinsertionCBMCProbability=0.5,
    swapCBCFCMCProbability=1.0,
)

co2 = raspalib.Component(
    componentId=0,
    forceField=ff,
    componentName="CO2",
    fileName="CO2.json",
    particleProbabilities=mcmoves,
)

mfi = raspalib.Framework(
    frameworkId=0,
    forceField=ff,
    componentName="MFI_SI",
    fileName="MFI_SI.cif",
    numberOfUnitCells=raspalib.int3(2, 2, 2),
)

system = raspalib.System(
    systemId=0,
    externalTemperature=353.0,
    externalPressure=1e5,
    forceField=ff,
    components=[co2],
    initialNumberOfMolecules=[0],
    frameworkComponents=[mfi],
)

mc = raspalib.MonteCarlo(
    numberOfCycles=100000,
    numberOfInitializationCycles=50000,
    numberOfEquilibrationCycles=50000,
    printEvery=1000,
    systems=[system],
    outputToFiles=True,
)

mc.run()
