import raspa

ff = raspa.ForceField(fileName="force_field.json")

mcmoves = raspa.MCMoveProbabilitiesParticles(
    translationProbability=0.5,
    rotationProbability=0.5,
    reinsertionCBMCProbability=0.5,
    swapCBMCProbability=1.0,
    widomProbability=1.0,
)
co2 = raspa.Component.exampleCO2(0, ff, particleProbabilities=mcmoves)
cubtc = raspa.Framework(0, ff, "Cu-BTC", "Cu-BTC.cif", numberOfUnitCells=[1, 1, 1])

system = raspa.System(
    systemId=0,
    temperature=323.0,
    forceField=ff,
    components=[co2],
    initialNumberOfMolecules=[1],
    pressure=1e4,
    frameworkComponents=[cubtc],
)

mc = raspa.MonteCarlo(
    numberOfCycles=10000,
    numberOfInitializationCycles=2000,
    numberOfEquilibrationCycles=0,
    printEvery=500,
    systems=[system],
)

mc.run()
