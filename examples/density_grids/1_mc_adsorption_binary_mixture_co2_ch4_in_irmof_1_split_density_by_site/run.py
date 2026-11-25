"""Monte Carlo: adsorption binary mixture co2 ch4 in irmof 1 split density by site.

The simulation of simulation.json, set up with raspalib: the force field, the molecules and the framework
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
)
co2.mol_fraction = 0.25
co2.ideal_gas_rosenbluth_weight = 1.0

methane = raspalib.Component.from_file(
    force_field,
    "methane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(
        translation_probability=0.5,
        rotation_probability=0.5,
        reinsertion_cbmc_probability=0.5,
        swap_cbmc_probability=1.0,
        widom_probability=1.0,
    ),
)
methane.mol_fraction = 0.75
methane.ideal_gas_rosenbluth_weight = 1.0

irmof_1 = raspalib.Framework.from_cif(force_field, "IRMOF-1.cif", number_of_unit_cells=[1, 1, 1])

system = raspalib.System(
    force_field=force_field,
    framework_components=irmof_1,
    external_temperature=300.0,
    external_pressure=1000000.0,
    helium_void_fraction=0.81,
    components=[co2, methane],
    initial_number_of_molecules=[0, 0],
)
density_grid = raspalib.PropertyDensityGrid(
    number_of_grid_points=(128, 128, 128),
    sample_every=10,
    write_every=10000,
    density_grid_pseudo_atoms_list=[force_field.find_pseudo_atom(name) for name in ("C_co2", "O_co2")],
)
system.set_property_density_grid(density_grid)

mc = raspalib.MonteCarlo(
    number_of_production_cycles=1000000,
    number_of_initialization_cycles=100000,
    print_every=10000,
    systems=[system],
    output_to_files=True,
)
mc.run()

energy, error = system.average_energies.result()
print(f"average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
print(f"molecules per component: {system.number_of_molecules_per_component}")
