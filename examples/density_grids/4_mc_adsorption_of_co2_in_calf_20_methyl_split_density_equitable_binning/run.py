"""Monte Carlo: adsorption of co2 in calf 20 methyl split density equitable binning.

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
        reinsertion_cbmc_probability=0.5,
        swap_cbmc_probability=1.0,
        widom_probability=1.0,
    ),
    fugacity_coefficient=1.0,
)
co2.ideal_gas_rosenbluth_weight = 1.0

calf_20_methyl = raspalib.Framework.from_cif(
    force_field,
    "CALF-20-methyl.cif",
    number_of_unit_cells=[2, 2, 2],
)

system = raspalib.System(
    force_field=force_field,
    framework_components=calf_20_methyl,
    external_temperature=296.0,
    external_pressure=1000000.0,
    components=[co2],
    initial_number_of_molecules=[0],
)
density_grid = raspalib.PropertyDensityGrid(
    number_of_grid_points=(118, 112, 138),
    sample_every=10,
    write_every=10000,
    density_grid_pseudo_atoms_list=[force_field.find_pseudo_atom(name) for name in ("Cx", "Ox")],
    binning_mode=raspalib.PropertyDensityGrid.Binning.EQUITABLE,
)
system.set_property_density_grid(density_grid)

mc = raspalib.MonteCarlo(
    number_of_production_cycles=3000000,
    number_of_initialization_cycles=100000,
    print_every=10000,
    systems=[system],
    output_to_files=True,
)
mc.run()

energy, error = system.average_energies.result()
print(f"average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
print(f"molecules per component: {system.number_of_molecules_per_component}")
