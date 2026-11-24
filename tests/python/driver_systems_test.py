"""The systems given to a driver hold its results when a phase returns.

A driver works on its own copy of the systems it is given. When 'run' or a phase returns, the systems of the
driver are assigned back to the objects that were given, so these can be read like 'driver.systems'.

The systems are set up the way the RASPA workshop notebooks do it, with the keywords of the 3.1.0 release.

Run with pytest, or with the standard library:

    python -m unittest discover -s tests/python -p "*_test.py"

The file name ends in _test because the repository ignores files that start with test_.
"""

import os
import unittest

os.environ.setdefault("RASPALIB_QUIET", "1")

import raspalib


def co2_force_field() -> raspalib.ForceField:
    # CO2: A. Garcia-Sanchez et al., J. Phys. Chem. C 2009, 113, 8814-8820
    return raspalib.ForceField(
        pseudo_atoms=[
            raspalib.PseudoAtom(name="C_co2", framework_type=False, mass=12.0, charge=0.6512, atomic_number=6),
            raspalib.PseudoAtom(name="O_co2", framework_type=False, mass=15.9994, charge=-0.3256, atomic_number=8),
        ],
        parameters=[raspalib.VDWParameters(29.933, 2.745), raspalib.VDWParameters(85.671, 3.017)],
        mixing_rule=raspalib.ForceField.MixingRule.LORENTZ_BERTHELOT,
        cutoff_framework_vdw=12.0,
        cutoff_molecule_vdw=12.0,
        cutoff_coulomb=12.0,
        shifted=True,
        tail_corrections=False,
        use_charge=True,
    )


def co2_system(
    number_of_molecules: int = 40,
    move_probabilities: raspalib.MCMoveProbabilities | None = None,
    system_probabilities: raspalib.MCMoveProbabilities | None = None,
) -> raspalib.System:
    force_field = co2_force_field()
    if move_probabilities is None:
        move_probabilities = raspalib.MCMoveProbabilities(
            translation_probability=0.5, rotation_probability=0.5, reinsertion_cbmc_probability=0.5
        )
    component = raspalib.Component(
        force_field=force_field,
        component_name="CO2",
        critical_temperature=304.1282,
        critical_pressure=7377300.0,
        acentric_factor=0.22394,
        defined_atoms=[
            raspalib.Atom(position=(0.0, 0.0, 1.149), charge=-0.3256, type=1),
            raspalib.Atom(position=(0.0, 0.0, 0.0), charge=0.6512, type=0),
            raspalib.Atom(position=(0.0, 0.0, -1.149), charge=-0.3256, type=1),
        ],
        number_of_blocks=5,
        number_of_lambda_bins=21,
        particle_probabilities=move_probabilities,
        fugacity_coefficient=1.0,
        thermodynamic_integration=False,
    )
    return raspalib.System(
        force_field=force_field,
        simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
        has_external_field=False,
        external_temperature=300.0,
        framework_components=None,
        components=[component],
        initial_number_of_molecules=[number_of_molecules],
        number_of_blocks=5,
        system_probabilities=system_probabilities or raspalib.MCMoveProbabilities(),
    )


def average_energy(system: raspalib.System) -> float:
    average, _ = system.average_energies.result()
    return average.total_energy.energy


def translation_counts(component: raspalib.Component) -> float:
    return sum(component.mc_moves_statistics[raspalib.Move.Types.TRANSLATION].total_counts)


class MonteCarloCopiesOut(unittest.TestCase):
    def test_run_leaves_the_results_in_the_given_system(self):
        system = co2_system()
        self.assertEqual(average_energy(system), 0.0)

        mc = raspalib.MonteCarlo(
            number_of_cycles=40, number_of_initialization_cycles=20, systems=[system], random_seed=12
        )
        mc.run()

        self.assertNotEqual(average_energy(system), 0.0)
        self.assertEqual(average_energy(system), average_energy(mc.systems[0]))
        self.assertEqual(
            system.number_of_molecules_per_component, mc.systems[0].number_of_molecules_per_component
        )

    def test_every_phase_copies_out(self):
        system = co2_system()
        mc = raspalib.MonteCarlo(
            number_of_cycles=40, number_of_initialization_cycles=20, systems=[system], random_seed=12
        )

        mc.setup()
        mc.initialize()
        self.assertGreater(translation_counts(system.components[0]), 0)
        self.assertEqual(average_energy(system), 0.0)

        mc.equilibrate()
        mc.production()
        self.assertNotEqual(average_energy(system), 0.0)
        self.assertEqual(average_energy(system), average_energy(mc.systems[0]))

        mc.tear_down()
        self.assertEqual(average_energy(system), average_energy(mc.systems[0]))

    def test_inside_a_call_back_the_systems_of_the_driver_are_the_live_ones(self):
        system = co2_system()
        mc = raspalib.MonteCarlo(
            number_of_cycles=40, number_of_initialization_cycles=20, systems=[system], random_seed=12
        )
        seen = []

        def call_back():
            seen.append((translation_counts(mc.systems[0].components[0]), translation_counts(system.components[0])))

        mc.setup()
        mc.initialize(call_back_function=call_back, call_back_every=10)
        mc.tear_down()

        self.assertGreater(len(seen), 0)
        live, given = seen[-1]
        self.assertGreater(live, 0)
        # the given system is brought up to date when the phase returns, not during it
        self.assertEqual(given, 0)
        self.assertGreater(translation_counts(system.components[0]), 0)

    def test_a_component_taken_before_the_run_shows_the_results(self):
        system = co2_system()
        component = system.components[0]
        mc = raspalib.MonteCarlo(
            number_of_cycles=40, number_of_initialization_cycles=20, systems=[system], random_seed=12
        )
        mc.run()

        self.assertGreater(translation_counts(component), 0)
        self.assertEqual(translation_counts(component), translation_counts(mc.systems[0].components[0]))
        self.assertEqual(component.name, "CO2")

    def test_two_systems_of_a_gibbs_ensemble(self):
        move_probabilities = raspalib.MCMoveProbabilities(
            translation_probability=0.5,
            rotation_probability=0.5,
            reinsertion_cbmc_probability=0.5,
            gibbs_swap_cbmc_probability=1.0,
        )
        system_probabilities = raspalib.MCMoveProbabilities(gibbs_volume_change_probability=0.01)
        system_a = co2_system(30, move_probabilities, system_probabilities)
        system_b = co2_system(30, move_probabilities, system_probabilities)

        mc = raspalib.MonteCarlo(
            number_of_cycles=40,
            number_of_initialization_cycles=20,
            systems=[system_a, system_b],
            random_seed=12,
        )
        mc.run()

        for given, own in zip((system_a, system_b), mc.systems):
            self.assertEqual(average_energy(given), average_energy(own))
            self.assertEqual(given.number_of_molecules_per_component, own.number_of_molecules_per_component)
        self.assertEqual(
            system_a.number_of_molecules_per_component[0] + system_b.number_of_molecules_per_component[0], 60
        )

    def test_a_system_given_twice_is_left_as_it_was(self):
        system = co2_system()
        mc = raspalib.MonteCarlo(
            number_of_cycles=20, number_of_initialization_cycles=10, systems=[system, system], random_seed=12
        )
        mc.run()

        self.assertEqual(average_energy(system), 0.0)
        self.assertNotEqual(average_energy(mc.systems[0]), 0.0)

    def test_the_systems_of_one_driver_given_to_the_next(self):
        system = co2_system()
        mc = raspalib.MonteCarlo(
            number_of_cycles=20, number_of_initialization_cycles=20, systems=[system], random_seed=12
        )
        mc.run()

        md = raspalib.MolecularDynamics(
            number_of_cycles=20, number_of_initialization_cycles=10, systems=mc.systems, random_seed=12
        )
        md.run()
        self.assertEqual(average_energy(mc.systems[0]), average_energy(md.systems[0]))


class MolecularDynamicsCopiesOut(unittest.TestCase):
    def test_sampled_properties_are_read_from_the_given_system(self):
        system = co2_system()
        system.time_step = 0.0005
        system.set_property_msd(raspalib.PropertyMeanSquaredDisplacement(sample_every=1))
        system.set_property_conserved_energy_evolution(
            raspalib.PropertyConservedEnergyEvolution(number_of_cycles=80, sample_every=1)
        )
        md = raspalib.MolecularDynamics(
            number_of_cycles=50,
            number_of_initialization_cycles=20,
            number_of_equilibration_cycles=10,
            systems=[system],
            random_seed=12,
        )
        md.setup()
        md.initialize()
        md.equilibrate()
        md.production()
        md.tear_down()

        given = [energy.conserved_energy() for energy in system.property_conserved_energy_evolution.result()]
        own = [energy.conserved_energy() for energy in md.systems[0].property_conserved_energy_evolution.result()]
        self.assertGreater(len(given), 0)
        self.assertEqual(given, own)
        self.assertTrue(any(value != 0.0 for value in given))

        given_msd = [(point.time, point.xyz) for point in system.property_msd.result()[0]]
        own_msd = [(point.time, point.xyz) for point in md.systems[0].property_msd.result()[0]]
        self.assertGreater(len(given_msd), 0)
        self.assertEqual(given_msd, own_msd)


class Keywords(unittest.TestCase):
    def test_the_production_cycles_go_by_both_keywords(self):
        for driver in (raspalib.MonteCarlo, raspalib.MolecularDynamics):
            driver(number_of_cycles=10, systems=[co2_system()])
            driver(number_of_production_cycles=10, systems=[co2_system()])
            driver(10)
        raspalib.PropertyConservedEnergyEvolution(number_of_cycles=10, sample_every=1)
        raspalib.PropertyConservedEnergyEvolution(number_of_production_cycles=10, sample_every=1)
        raspalib.PropertyNumberOfMoleculesEvolution(number_of_cycles=10, number_of_components=1, sample_every=1)
        raspalib.PropertyVolumeEvolution(number_of_cycles=10, sample_every=1)

    def test_systems_have_to_be_systems(self):
        with self.assertRaises(TypeError):
            raspalib.MonteCarlo(number_of_cycles=10, systems=[1])
        with self.assertRaises(TypeError):
            raspalib.MonteCarlo(number_of_cycles=10, systems=co2_system())


if __name__ == "__main__":
    unittest.main()
