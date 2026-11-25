"""Monte Carlo: ideal gas rosenbluth weight butane.

The simulation of simulation.json, set up with raspalib: the force field and the molecule
are read from the files of this directory.

    python run.py
"""

import raspalib

force_field = raspalib.ForceField.from_file("force_field.json")
force_field.use_charge = False

butane = raspalib.Component.from_file(
    force_field,
    "butane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

system = raspalib.System(
    force_field=force_field,
    simulation_box=raspalib.SimulationBox(30.0, 30.0, 30.0),
    external_temperature=300.0,
    components=[butane],
    initial_number_of_molecules=[0],
)

mc = raspalib.MonteCarlo(
    number_of_production_cycles=1000000,
    print_every=5000,
    systems=[system],
    output_to_files=True,
)
mc.run()

energy, error = system.average_energies.result()
print(f"average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
