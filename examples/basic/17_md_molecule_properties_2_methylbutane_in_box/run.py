import raspalib

ff = raspalib.ForceField("force_field.json")
ff.useCharge = False

methylbutane = raspalib.Component(
    componentId=0,
    forceField=ff,
    componentName="2-methylbutane",
    fileName="2-methylbutane.json",
)

box = raspalib.SimulationBox(25.0, 25.0, 25.0)

system = raspalib.System(
    systemId=0,
    externalTemperature=298.0,
    forceField=ff,
    components=[methylbutane],
    initialNumberOfMolecules=[32],
    simulationBox=box,
    ensemble=raspalib.Ensemble.NVT,
)

md = raspalib.MolecularDynamics(
    numberOfCycles=5000000,
    numberOfInitializationCycles=5000,
    numberOfEquilibrationCycles=10000,
    printEvery=10000,
    systems=[system],
    outputToFiles=True,
)
md.run()
