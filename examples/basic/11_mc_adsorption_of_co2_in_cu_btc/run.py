"""Monte Carlo: adsorption of co2 in cu btc.

The simulation of simulation.json, set up with raspalib: the force field, the molecule and the framework
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
        swap_cbmc_probability=1.0,
        widom_probability=1.0,
    ),
    fugacity_coefficient=1.0,
)
co2.ideal_gas_rosenbluth_weight = 1.0

cu_btc = raspalib.Framework.from_cif(force_field, "Cu-BTC.cif", number_of_unit_cells=[1, 1, 1])

system = raspalib.System(
    force_field=force_field,
    framework_components=cu_btc,
    external_temperature=323.0,
    external_pressure=10000.0,
    components=[co2],
    initial_number_of_molecules=[0],
)

mc = raspalib.MonteCarlo(
    number_of_production_cycles=500000,
    number_of_initialization_cycles=50000,
    print_every=5000,
    systems=[system],
    output_to_files=True,
)
mc.run()

energy, error = system.average_energies.result()
print(f"average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
print(f"molecules per component: {system.number_of_molecules_per_component}")
