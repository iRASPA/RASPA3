import raspalib

ff = raspalib.ForceField("force_field.json")
ff.useCharge = False

mcmoves = raspalib.MCMoveProbabilities(
    translationProbability=0.5,
    reinsertionCBMCProbability=0.5,
)

methane = raspalib.Component(
    componentId=0,
    forceField=ff,
    componentName="methane",
    fileName="methane.json",
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
    externalTemperature=300.0,
    forceField=ff,
    components=[methane],
    initialNumberOfMolecules=[1],
    frameworkComponents=[mfi],
)

mc = raspalib.MonteCarlo(numberOfCycles=5000, numberOfInitializationCycles=5000, systems=[system], outputToFiles=True)

mc.run()
