"""Monte Carlo: co2 n2 in two independent boxes.

The simulation of simulation.json, set up with raspalib: the force field and the molecules
are read from the files of this directory.

    python run.py
"""

import raspalib

force_field = raspalib.ForceField.from_file("force_field.json")

co2 = raspalib.Component.from_file(
    force_field,
    "CO2.json",
    particle_probabilities=raspalib.MCMoveProbabilities(
        translation_probability=1.0,
        rotation_probability=1.0,
        reinsertion_cbmc_probability=1.0,
    ),
)

n2 = raspalib.Component.from_file(
    force_field,
    "N2.json",
    particle_probabilities=raspalib.MCMoveProbabilities(
        translation_probability=1.0,
        rotation_probability=1.0,
        reinsertion_cbmc_probability=1.0,
    ),
)

system_0 = raspalib.System(
    force_field=force_field,
    simulation_box=raspalib.SimulationBox(25.0, 25.0, 25.0),
    external_temperature=300.0,
    components=[co2, n2],
    initial_number_of_molecules=[100, 0],
)
system_0.set_sample_pdb_movie(raspalib.SampleMovie(sample_every=10))

system_1 = raspalib.System(
    force_field=force_field,
    simulation_box=raspalib.SimulationBox(30.0, 30.0, 30.0, 90.0, 120.0, 90.0),
    external_temperature=500.0,
    components=[co2, n2],
    initial_number_of_molecules=[0, 100],
)
system_1.set_sample_pdb_movie(raspalib.SampleMovie(sample_every=10))

mc = raspalib.MonteCarlo(
    number_of_production_cycles=10000,
    number_of_initialization_cycles=1000,
    print_every=1000,
    systems=[system_0, system_1],
    output_to_files=True,
)
mc.run()

energy, error = system_0.average_energies.result()
print(f"system 0: average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
print(f"system 0: molecules per component: {system_0.number_of_molecules_per_component}")
energy, error = system_1.average_energies.result()
print(f"system 1: average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
print(f"system 1: molecules per component: {system_1.number_of_molecules_per_component}")
