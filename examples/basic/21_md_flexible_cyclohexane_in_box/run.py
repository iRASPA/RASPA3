"""Molecular dynamics: flexible cyclohexane in box.

The simulation of simulation.json, set up with raspalib: the force field and the molecule
are read from the files of this directory.

    python run.py
"""

import raspalib

force_field = raspalib.ForceField.from_file("force_field.json")
force_field.use_charge = False

cyclohexane = raspalib.Component.from_file(
    force_field,
    "cyclohexane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(
        translation_probability=1.0,
        rotation_probability=1.0,
        reinsertion_cbmc_probability=1.0,
    ),
)
cyclohexane.molecule_properties_settings = raspalib.MoleculePropertiesSettings()

system = raspalib.System(
    force_field=force_field,
    simulation_box=raspalib.SimulationBox(30.0, 30.0, 30.0),
    external_temperature=300.0,
    components=[cyclohexane],
    initial_number_of_molecules=[30],
)
system.set_thermostat(raspalib.Thermostat())

md = raspalib.MolecularDynamics(
    number_of_production_cycles=100000,
    number_of_initialization_cycles=5000,
    number_of_equilibration_cycles=10000,
    print_every=5000,
    systems=[system],
    output_to_files=True,
)
md.run()

energy, error = system.average_energies.result()
print(f"average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
