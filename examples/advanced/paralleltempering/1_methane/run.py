import raspa
import numpy as np


ff = raspa.ForceField.exampleMoleculeForceField(useCharge=False)
mcmoves = raspa.MCMoveProbabilitiesParticles(1.0)
methane = raspa.Component.exampleCH4(0, ff, particleProbabilities=mcmoves)

sysmoves = raspa.MCMoveProbabilitiesSystem(parallelTemperingProbability=0.001)
systems = [
    raspa.System(
        systemId=i,
        temperature=200.0 + 50.0 * i,
        forceField=ff,
        components=[methane],
        initialNumberOfMolecules=[100],
        systemProbabilities=sysmoves,
        simulationBox=raspa.SimulationBox(30.0 * np.ones(3)),
    )
    for i in range(8)
]

mc = raspa.MonteCarlo(numberOfCycles=10000, numberOfInitializationCycles=10000, printEvery=100, systems=systems)
mc.run()
