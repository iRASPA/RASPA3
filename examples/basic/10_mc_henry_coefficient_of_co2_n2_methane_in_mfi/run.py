import raspalib

ff = raspalib.ForceField("force_field.json")

mcmoves = raspalib.MCMoveProbabilities(widomProbability=1.0)

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

methane = raspalib.Component(
    componentId=2,
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

# define systems
system = raspalib.System(
    systemId=0,
    externalTemperature=300.0,
    forceField=ff,
    components=[co2, n2, methane],
    initialNumberOfMolecules=[0, 0, 0],
    frameworkComponents=[mfi],
)

mc = raspalib.MonteCarlo(
    numberOfCycles=20000,
    numberOfInitializationCycles=0,
    systems=[system],
    outputToFiles=True,
)

mc.run()
