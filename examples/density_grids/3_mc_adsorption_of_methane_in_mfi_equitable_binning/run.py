"""Monte Carlo: adsorption of methane in mfi equitable binning.

The simulation of simulation.json, set up with raspalib: the force field, the molecule and the framework
are read from the files of this directory.

    python run.py
"""

import raspalib

force_field = raspalib.ForceField.from_file("force_field.json")
force_field.use_charge = False

methane = raspalib.Component.from_file(
    force_field,
    "methane.json",
    particle_probabilities=raspalib.MCMoveProbabilities(
        translation_probability=0.5,
        reinsertion_cbmc_probability=0.5,
        swap_cbmc_probability=1.0,
        widom_probability=1.0,
    ),
    fugacity_coefficient=1.0,
)
methane.ideal_gas_rosenbluth_weight = 1.0

mfi_si = raspalib.Framework.from_cif(force_field, "MFI_SI.cif", number_of_unit_cells=[2, 2, 2])

system = raspalib.System(
    force_field=force_field,
    framework_components=mfi_si,
    external_temperature=300.0,
    external_pressure=100000.0,
    components=[methane],
    initial_number_of_molecules=[0],
)
density_grid = raspalib.PropertyDensityGrid(
    number_of_grid_points=(268, 266, 180),
    sample_every=10,
    write_every=50000,
    binning_mode=raspalib.PropertyDensityGrid.Binning.EQUITABLE,
)
system.set_property_density_grid(density_grid)

mc = raspalib.MonteCarlo(
    number_of_production_cycles=2500000,
    number_of_initialization_cycles=10000,
    print_every=10000,
    systems=[system],
    output_to_files=True,
)
mc.run()

energy, error = system.average_energies.result()
print(f"average energy: {energy.total_energy.energy:.3f} +/- {error.total_energy.energy:.3f} K")
print(f"molecules per component: {system.number_of_molecules_per_component}")
