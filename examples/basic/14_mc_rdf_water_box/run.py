"""Monte Carlo: rdf water box.

The simulation of simulation.json, set up with raspalib: the force field and the molecule
are read from the files of this directory.

    python run.py
"""

import raspalib

force_field = raspalib.ForceField.from_file("force_field.json")

water = raspalib.Component.from_file(
    force_field,
    "water.json",
    particle_probabilities=raspalib.MCMoveProbabilities(
        translation_probability=0.5,
        rotation_probability=0.5,
        reinsertion_cbmc_probability=1.0,
    ),
)

system = raspalib.System(
    force_field=force_field,
    simulation_box=raspalib.SimulationBox(24.83, 24.83, 24.83),
    external_temperature=298.0,
    components=[water],
    initial_number_of_molecules=[512],
)
rdf = raspalib.PropertyConventionalRadialDistributionFunction(
    number_of_bins=128,
    range=12.0,
    sample_every=10,
    write_every=100,
)
system.set_property_conventional_rdf(rdf)

mc = raspalib.MonteCarlo(
    number_of_production_cycles=10000,
    number_of_initialization_cycles=10000,
    print_every=100,
    systems=[system],
    output_to_files=True,
)
mc.run()

energy, error = system.average_energies.result()
print(f"average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
