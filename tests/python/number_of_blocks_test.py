"""The number of blocks is given once.

The averages of a run are kept in blocks, by the components (their lambda histograms and Rosenbluth weights), by
the systems (their averages and sampled properties) and by the driver, which says in which block a sample goes.
The driver is the place to give the number: it sizes its systems for it, as a system sizes its components. Where
the number is left out it is taken from below, and what is below then has to agree.

    python -m unittest discover -s tests/python -p "*_test.py"
"""

import os
import unittest

os.environ.setdefault("RASPALIB_QUIET", "1")

import raspalib

from driver_systems_test import co2_force_field


def co2_component(move_probabilities=None, **keywords) -> raspalib.Component:
    force_field = co2_force_field()
    return raspalib.Component(
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
        number_of_lambda_bins=21,
        particle_probabilities=move_probabilities
        or raspalib.MCMoveProbabilities(translation_probability=0.5, rotation_probability=0.5),
        fugacity_coefficient=1.0,
        **keywords,
    )


def co2_system(components=None, number_of_molecules=20, **keywords) -> raspalib.System:
    components = components if components is not None else [co2_component()]
    return raspalib.System(
        force_field=co2_force_field(),
        simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
        external_temperature=300.0,
        external_pressure=1.0e5,
        components=components,
        initial_number_of_molecules=[number_of_molecules] * len(components),
        **keywords,
    )


class TheNumberOfBlocksIsGivenOnce(unittest.TestCase):
    def test_given_to_the_driver_it_sizes_the_systems_and_their_components(self):
        system = co2_system()
        self.assertEqual(system.number_of_blocks, 5)

        mc = raspalib.MonteCarlo(
            number_of_cycles=40, number_of_initialization_cycles=10, systems=[system], number_of_blocks=8
        )
        self.assertEqual(mc.number_of_blocks, 8)
        self.assertEqual(mc.systems[0].number_of_blocks, 8)
        self.assertEqual(mc.systems[0].components[0].number_of_blocks, 8)
        mc.run()

        # copied out: the given system is that of the driver
        self.assertEqual(system.number_of_blocks, 8)
        self.assertEqual(system.components[0].number_of_blocks, 8)
        average, error = system.average_energies.result()
        self.assertNotEqual(average.total_energy.energy, 0.0)
        self.assertGreater(error.total_energy.energy, 0.0)

    def test_left_out_everywhere_it_is_five(self):
        mc = raspalib.MonteCarlo(number_of_cycles=10, systems=[co2_system()])
        self.assertEqual(mc.number_of_blocks, 5)
        self.assertEqual(raspalib.MonteCarlo(number_of_cycles=10).number_of_blocks, 5)

    def test_given_to_the_system_it_sizes_the_components_and_the_driver_takes_it(self):
        system = co2_system(number_of_blocks=7)
        self.assertEqual(system.number_of_blocks, 7)
        self.assertEqual(system.components[0].number_of_blocks, 7)

        md = raspalib.MolecularDynamics(number_of_cycles=10, systems=[system])
        self.assertEqual(md.number_of_blocks, 7)
        md.run()

    def test_given_to_the_components_the_system_takes_it(self):
        system = co2_system(components=[co2_component(number_of_blocks=6)])
        self.assertEqual(system.number_of_blocks, 6)
        self.assertEqual(raspalib.MonteCarlo(number_of_cycles=10, systems=[system]).number_of_blocks, 6)

    def test_components_that_differ_need_the_number_of_the_system(self):
        components = [co2_component(number_of_blocks=6), co2_component()]
        with self.assertRaises(ValueError):
            co2_system(components=components)
        system = co2_system(components=components, number_of_blocks=4)
        self.assertEqual([component.number_of_blocks for component in system.components], [4, 4])
        # the components that were given are not those of the system
        self.assertEqual([component.number_of_blocks for component in components], [6, 5])

    def test_systems_that_differ_need_the_number_of_the_driver(self):
        systems = [co2_system(number_of_blocks=4), co2_system(number_of_blocks=6)]
        with self.assertRaises(ValueError):
            raspalib.MonteCarlo(number_of_cycles=10, systems=systems)
        mc = raspalib.MonteCarlo(number_of_cycles=20, systems=systems, number_of_blocks=10)
        mc.run()
        self.assertEqual([system.number_of_blocks for system in systems], [10, 10])

    def test_zero_blocks_are_refused(self):
        with self.assertRaises(ValueError):
            co2_system(number_of_blocks=0)
        with self.assertRaises(ValueError):
            raspalib.MonteCarlo(number_of_cycles=10, systems=[co2_system()], number_of_blocks=0)

    def test_the_statistics_of_the_components_follow_the_driver(self):
        # a lambda histogram of five blocks took the samples of ten: the later blocks went past its end
        probabilities = raspalib.MCMoveProbabilities(
            translation_probability=0.5, rotation_probability=0.5, swap_cfcmc_probability=1.0, widom_probability=0.5
        )
        system = co2_system(components=[co2_component(probabilities)])
        mc = raspalib.MonteCarlo(
            number_of_cycles=200,
            number_of_initialization_cycles=50,
            number_of_equilibration_cycles=50,
            systems=[system],
            number_of_blocks=10,
            random_seed=12,
        )
        mc.run()

        component = system.components[0]
        self.assertEqual(component.number_of_blocks, 10)
        self.assertGreater(sum(component.lambda_histogram.result()[0]), 0.0)
        weight, _ = component.average_rosenbluth_weights.result()
        self.assertGreater(weight, 0.0)

    def test_the_sampled_properties_follow_the_driver(self):
        system = co2_system()
        molecules = raspalib.PropertyNumberOfMoleculesHistogram(value_range=(0, 50))
        energies = raspalib.PropertyEnergyHistogram(number_of_bins=20, value_range=(-20000.0, 0.0))
        rdf = raspalib.PropertyConventionalRadialDistributionFunction(number_of_bins=32, range=12.0, sample_every=1)
        system.set_number_of_molecules_histogram(molecules)
        system.set_average_energy_histogram(energies)
        system.set_property_conventional_rdf(rdf)

        mc = raspalib.MonteCarlo(number_of_cycles=40, systems=[system], number_of_blocks=10, random_seed=12)
        mc.run()

        _, histogram, _ = molecules.result()[0]
        self.assertEqual(max(range(50), key=lambda n: histogram[n]), 20)
        _, average, _ = energies.result()
        self.assertGreater(sum(bin.total_energy for bin in average), 0.0)

    def test_a_system_keeps_its_results_in_a_driver_of_the_same_number(self):
        system = co2_system()
        raspalib.MonteCarlo(number_of_cycles=20, systems=[system], number_of_blocks=8, random_seed=12).run()
        energy = system.average_energies.result()[0].total_energy.energy
        self.assertNotEqual(energy, 0.0)

        same = raspalib.MonteCarlo(number_of_cycles=20, systems=[system], random_seed=12)
        self.assertEqual(same.number_of_blocks, 8)
        self.assertEqual(same.systems[0].average_energies.result()[0].total_energy.energy, energy)

        # in another number of blocks the averages start anew
        other = raspalib.MonteCarlo(number_of_cycles=20, systems=[system], number_of_blocks=4, random_seed=12)
        self.assertEqual(other.systems[0].average_energies.result()[0].total_energy.energy, 0.0)
        self.assertEqual(system.average_energies.result()[0].total_energy.energy, energy)


if __name__ == "__main__":
    unittest.main()
