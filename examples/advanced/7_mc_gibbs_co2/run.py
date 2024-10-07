import raspa
import numpy as np

ff = raspa.ForceField(fileName="force_field.json")

mcmoves = raspa.MCMoveProbabilitiesParticles(
    translationProbability=0.5,
    rotationProbability=0.5,
    reinsertionCBMCProbability=0.5,
    gibbsSwapCFCMCProbability=1.0,
)
co2 = raspa.Component.exampleCO2(0, ff, particleProbabilities=mcmoves)


sysMoves = raspa.MCMoveProbabilitiesSystem(volumeChangeProbability=0.01)
system0 = raspa.System(
    systemId=0,
    temperature=240.0,
    forceField=ff,
    components=[co2],
    initialNumberOfMolecules=[256],
    pressure=1e4,
    simulationBox=raspa.SimulationBox(30.0 * np.ones(3)),
    systemProbabilities=sysMoves,
)

system1 = raspa.System(
    systemId=1,
    temperature=240.0,
    forceField=ff,
    components=[co2],
    initialNumberOfMolecules=[256],
    simulationBox=raspa.SimulationBox(30.0 * np.ones(3)),
    systemProbabilities=sysMoves,
)

mc = raspa.MonteCarlo(
    numberOfCycles=20000,
    numberOfInitializationCycles=10000,
    numberOfEquilibrationCycles=10000,
    printEvery=1000,
    systems=[system0, system1],
)

mc.run()
