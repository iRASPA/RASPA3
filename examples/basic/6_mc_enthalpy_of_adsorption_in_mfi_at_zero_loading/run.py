"""Monte Carlo: enthalpy of adsorption in mfi at zero loading.

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
    ),
)

mfi_si = raspalib.Framework.from_cif(force_field, "MFI_SI.cif", number_of_unit_cells=[2, 2, 2])

system = raspalib.System(
    force_field=force_field,
    framework_components=mfi_si,
    external_temperature=300.0,
    helium_void_fraction=0.29,
    components=[methane],
    initial_number_of_molecules=[1],
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
