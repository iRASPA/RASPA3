import raspa
import numpy as np

ff = raspa.ForceField.exampleMoleculeForceField(useCharge=False)

mcmoves = raspa.MCMoveProbabilitiesParticles(translationProbability=0.5, reinsertionCBMCProbability=0.5)
methane = raspa.Component.exampleCH4(0, ff, particleProbabilities=mcmoves)

sysmoves = raspa.MCMoveProbabilitiesSystem(volumeChangeProbability=0.05)
system = raspa.System(
    systemId=0,
    temperature=300.0,
    pressure=1e6,
    systemProbabilities=sysmoves,
    forceField=ff,
    components=[methane],
    initialNumberOfMolecules=[256],
    simulationBox=raspa.SimulationBox(30.0 * np.ones(3)),
)

mc = raspa.MonteCarlo(numberOfCycles=50000, numberOfInitializationCycles=10000, systems=[system])

mc.run()
