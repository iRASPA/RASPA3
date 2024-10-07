import raspa

ff = raspa.ForceField.exampleMoleculeForceField()

mcmoves = raspa.MCMoveProbabilitiesParticles(
    translationProbability=0.5,
    reinsertionCBMCProbability=0.5,
    swapProbability=1.0,
)
methane = raspa.Component.exampleCH4(0, ff, particleProbabilities=mcmoves)
mfi = raspa.Framework(0, ff, "MFI_SI", "MFI_SI.cif", numberOfUnitCells=[2, 2, 2])

system = raspa.System(0, 300.0, ff, [methane], [0], pressure=1e5, frameworkComponents=[mfi])

mc = raspa.MonteCarlo(
    numberOfCycles=10000,
    numberOfInitializationCycles=2000,
    numberOfEquilibrationCycles=0,
    systems=[system],
)

mc.run()
