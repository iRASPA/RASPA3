"""Monte Carlo: ideal gas rosenbluth weights c6 isomers 433k.

The simulation of simulation.json, set up with raspalib: the force field and the molecules
are read from the files of this directory.

    python run.py
"""

import raspalib

force_field = raspalib.ForceField.from_file("force_field.json")
force_field.use_charge = False

hexane = raspalib.Component.from_file(
    force_field,
    "hexane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

c_2_methylpentane = raspalib.Component.from_file(
    force_field,
    "2-methylpentane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

c_3_methylpentane = raspalib.Component.from_file(
    force_field,
    "3-methylpentane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

c_22_dimethylbutane = raspalib.Component.from_file(
    force_field,
    "22-dimethylbutane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

c_23_dimethylbutane = raspalib.Component.from_file(
    force_field,
    "23-dimethylbutane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

system = raspalib.System(
    force_field=force_field,
    simulation_box=raspalib.SimulationBox(30.0, 30.0, 30.0),
    external_temperature=433.0,
    components=[hexane, c_2_methylpentane, c_3_methylpentane, c_22_dimethylbutane, c_23_dimethylbutane],
    initial_number_of_molecules=[0, 0, 0, 0, 0],
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
