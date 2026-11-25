"""Monte Carlo: henry coefficient of co2 n2 methane in mfi.

The simulation of simulation.json, set up with raspalib: the force field, the molecules and the framework
are read from the files of this directory.

    python run.py
"""

import raspalib

force_field = raspalib.ForceField.from_file("force_field.json")

co2 = raspalib.Component.from_file(
    force_field,
    "CO2.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

n2 = raspalib.Component.from_file(
    force_field,
    "N2.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

methane = raspalib.Component.from_file(
    force_field,
    "methane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
)

mfi_si = raspalib.Framework.from_cif(force_field, "MFI_SI.cif", number_of_unit_cells=[2, 2, 2])

system = raspalib.System(
    force_field=force_field,
    framework_components=mfi_si,
    external_temperature=300.0,
    external_pressure=100000.0,
    components=[co2, n2, methane],
    initial_number_of_molecules=[0, 0, 0],
)

mc = raspalib.MonteCarlo(
    number_of_production_cycles=20000,
    print_every=1000,
    systems=[system],
    output_to_files=True,
)
mc.run()

energy, error = system.average_energies.result()
print(f"average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
print(f"molecules per component: {system.number_of_molecules_per_component}")
