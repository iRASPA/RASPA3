"""Monte Carlo: binary mixture propane butane in box.

The simulation of simulation.json, set up with raspalib: the force field and the molecules
are read from the files of this directory.

    python run.py
"""

import raspalib

force_field = raspalib.ForceField.from_file("force_field.json")
force_field.use_charge = False

propane = raspalib.Component.from_file(
    force_field,
    "propane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(
        translation_probability=1.0,
        rotation_probability=1.0,
        reinsertion_cbmc_probability=1.0,
        partial_reinsertion_cbmc_probability=1.0,
    ),
)

butane = raspalib.Component.from_file(
    force_field,
    "butane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(
        translation_probability=1.0,
        rotation_probability=1.0,
        reinsertion_cbmc_probability=1.0,
        partial_reinsertion_cbmc_probability=1.0,
    ),
)

system = raspalib.System(
    force_field=force_field,
    simulation_box=raspalib.SimulationBox(30.0, 30.0, 30.0),
    external_temperature=500.0,
    components=[propane, butane],
    initial_number_of_molecules=[50, 50],
)

mc = raspalib.MonteCarlo(
    number_of_production_cycles=20000,
    number_of_initialization_cycles=5000,
    print_every=1000,
    systems=[system],
    output_to_files=True,
)
mc.run()

energy, error = system.average_energies.result()
print(f"average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
