"""Monte Carlo: gibbs co2.

The simulation of simulation.json, set up with raspalib: the force field and the molecule
are read from the files of this directory.

    python run.py
"""

import raspalib

force_field = raspalib.ForceField.from_file("force_field.json")

co2 = raspalib.Component.from_file(
    force_field,
    "CO2.json",
    particle_probabilities=raspalib.MCMoveProbabilities(
        translation_probability=0.5,
        rotation_probability=0.5,
        reinsertion_cbmc_probability=0.5,
        gibbs_swap_cfcmc_probability=1.0,
        widom_probability=0.5,
    ),
    thermodynamic_integration=True,
)

system_0 = raspalib.System(
    force_field=force_field,
    simulation_box=raspalib.SimulationBox(30.0, 30.0, 30.0),
    external_temperature=240.0,
    components=[co2],
    initial_number_of_molecules=[256],
    system_probabilities=raspalib.MCMoveProbabilities(gibbs_volume_probability=0.01),
)

system_1 = raspalib.System(
    force_field=force_field,
    simulation_box=raspalib.SimulationBox(30.0, 30.0, 30.0),
    external_temperature=240.0,
    components=[co2],
    initial_number_of_molecules=[256],
    system_probabilities=raspalib.MCMoveProbabilities(gibbs_volume_probability=0.01),
)

mc = raspalib.MonteCarlo(
    number_of_production_cycles=10000,
    number_of_initialization_cycles=10000,
    number_of_equilibration_cycles=5000,
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
