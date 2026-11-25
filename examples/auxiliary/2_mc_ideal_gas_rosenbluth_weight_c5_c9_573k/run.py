"""Monte Carlo: ideal gas rosenbluth weight c5 c9 573k.

The simulation of simulation.json, set up with raspalib: the force field and the molecules
are read from the files of this directory.

    python run.py
"""

import raspalib

force_field = raspalib.ForceField.from_file("force_field.json")
force_field.use_charge = False

pentane = raspalib.Component.from_file(
    force_field,
    "pentane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

hexane = raspalib.Component.from_file(
    force_field,
    "hexane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

heptane = raspalib.Component.from_file(
    force_field,
    "heptane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

octane = raspalib.Component.from_file(
    force_field,
    "octane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

nonane = raspalib.Component.from_file(
    force_field,
    "nonane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

system = raspalib.System(
    force_field=force_field,
    simulation_box=raspalib.SimulationBox(30.0, 30.0, 30.0),
    external_temperature=573.0,
    components=[pentane, hexane, heptane, octane, nonane],
    initial_number_of_molecules=[0, 0, 0, 0, 0],
)

mc = raspalib.MonteCarlo(
    number_of_production_cycles=100000,
    print_every=5000,
    systems=[system],
    output_to_files=True,
)
mc.run()

energy, error = system.average_energies.result()
print(f"average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
