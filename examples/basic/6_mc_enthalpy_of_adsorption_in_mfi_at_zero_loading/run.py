import raspa

ff = raspa.ForceField.exampleMoleculeForceField(useCharge=False)

mcmoves = raspa.MCMoveProbabilitiesParticles(
    translationProbability=0.5,
    reinsertionCBMCProbability=0.5,
)
methane = raspa.Component.exampleCH4(0, ff, particleProbabilities=mcmoves)
mfi = raspa.Framework(0, ff, "MFI_SI", "MFI_SI.cif", numberOfUnitCells=[2, 2, 2])

system = raspa.System(
    systemId=0,
    temperature=300.0,
    forceField=ff,
    components=[methane],
    initialNumberOfMolecules=[1],
    pressure=1e5,
    frameworkComponents=[mfi],
)

mc = raspa.MonteCarlo(
    numberOfCycles=5000,
    numberOfInitializationCycles=5000,
    numberOfEquilibrationCycles=0,
    systems=[system],
)

mc.run()
