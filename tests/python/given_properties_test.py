"""The property given to a system holds the results when a phase returns.

A system builds its own property from the settings of the one it is given, sized for its molecules, components and
blocks. The object that was given is brought up to date when it is set and whenever the system is copied out of a
driver, so it can be read like 'system.property_...'.

A constructor asks for what the object uses: the sizes that the system or the driver know (the number of components,
blocks and pseudo-atoms, the cycles of the driver, the index of the system) need not be given.

    python -m unittest discover -s tests/python -p "*_test.py"
"""

import gc
import os
import tempfile
import unittest

os.environ.setdefault("RASPALIB_QUIET", "1")

import raspalib

from driver_systems_test import co2_system


def conserved(evolution) -> list[float]:
    return [energy.conserved_energy() for energy in evolution.result()]


def msd_points(msd) -> list[tuple[float, float]]:
    return [(point.time, point.xyz) for point in msd.result()[0]]


class GivenPropertiesAreCopiedOut(unittest.TestCase):
    def test_molecular_dynamics_leaves_the_results_in_the_given_properties(self):
        system = co2_system()
        system.time_step = 0.0005
        msd = raspalib.PropertyMeanSquaredDisplacement(sample_every=1)
        vacf = raspalib.PropertyVelocityAutoCorrelationFunction(
            number_of_buffers_vacf=4, buffer_length_vacf=16, sample_every=1
        )
        energy = raspalib.PropertyConservedEnergyEvolution(sample_every=1)
        system.set_property_msd(msd)
        system.set_property_vacf(vacf)
        system.set_property_conserved_energy_evolution(energy)

        md = raspalib.MolecularDynamics(
            number_of_cycles=50,
            number_of_initialization_cycles=20,
            number_of_equilibration_cycles=10,
            systems=[system],
            random_seed=12,
        )
        md.run()

        self.assertGreater(len(msd_points(msd)), 0)
        self.assertEqual(msd_points(msd), msd_points(md.systems[0].property_msd))
        self.assertEqual(
            [point.xyz for point in vacf.result()[0]],
            [point.xyz for point in md.systems[0].property_vacf.result()[0]],
        )
        self.assertEqual(conserved(energy), conserved(md.systems[0].property_conserved_energy_evolution))
        self.assertTrue(any(value != 0.0 for value in conserved(energy)))

    def test_every_phase_copies_out(self):
        system = co2_system()
        evolution = raspalib.PropertyNumberOfMoleculesEvolution(sample_every=1)
        system.set_property_number_of_molecules_evolution(evolution)
        mc = raspalib.MonteCarlo(
            number_of_cycles=20, number_of_initialization_cycles=10, systems=[system], random_seed=12
        )
        mc.setup()
        self.assertEqual(sum(evolution.result[0]), 0)
        mc.initialize()
        self.assertEqual(evolution.result[0][:10], [40] * 10)
        self.assertEqual(evolution.result[0][10:], [0] * 20)
        mc.equilibrate()
        mc.production()
        self.assertEqual(evolution.result[0], [40] * 30)

    def test_inside_a_call_back_the_given_property_is_that_of_the_last_phase(self):
        system = co2_system()
        evolution = raspalib.PropertyNumberOfMoleculesEvolution(sample_every=1)
        system.set_property_number_of_molecules_evolution(evolution)
        mc = raspalib.MonteCarlo(number_of_cycles=20, systems=[system], random_seed=12)
        seen = []

        def call_back():
            live = mc.systems[0].property_number_of_molecules_evolution.result[0]
            seen.append((sum(live), sum(evolution.result[0])))

        mc.setup()
        mc.production(call_back_function=call_back, call_back_every=10)

        live, given = seen[-1]
        self.assertGreater(live, 0)
        self.assertEqual(given, 0)
        self.assertEqual(sum(evolution.result[0]), 40 * 20)

    def test_a_property_given_to_two_systems_is_left_as_it_was(self):
        system_a, system_b = co2_system(30), co2_system(20)
        shared = raspalib.PropertyNumberOfMoleculesEvolution(sample_every=1)
        own = raspalib.PropertyVolumeEvolution(sample_every=1)
        system_a.set_property_number_of_molecules_evolution(shared)
        system_b.set_property_number_of_molecules_evolution(shared)
        system_a.set_property_volume_evolution(own)

        mc = raspalib.MonteCarlo(number_of_cycles=10, systems=[system_a, system_b], random_seed=12)
        mc.run()

        self.assertEqual(shared.result, [[]])
        self.assertEqual(system_a.property_number_of_molecules_evolution.result, [[30] * 10])
        self.assertEqual(system_b.property_number_of_molecules_evolution.result, [[20] * 10])
        self.assertEqual(own.result, [27000.0] * 10)

    def test_a_property_that_was_replaced_or_removed_is_left_as_it_was(self):
        system = co2_system()
        first = raspalib.PropertyVolumeEvolution(sample_every=1)
        second = raspalib.PropertyVolumeEvolution(sample_every=1)
        removed = raspalib.PropertyNumberOfMoleculesEvolution(sample_every=1)
        system.set_property_volume_evolution(first)
        system.set_property_volume_evolution(second)
        system.set_property_number_of_molecules_evolution(removed)
        system.set_property_number_of_molecules_evolution(None)

        raspalib.MonteCarlo(number_of_cycles=10, systems=[system], random_seed=12).run()

        self.assertEqual(first.result, [])
        self.assertEqual(second.result, [27000.0] * 10)
        self.assertEqual(removed.result, [[]])
        self.assertIsNone(system.property_number_of_molecules_evolution)

    def test_the_given_property_is_sized_by_the_system_when_it_is_set(self):
        system = co2_system()
        evolution = raspalib.PropertyNumberOfMoleculesEvolution(number_of_cycles=8, sample_every=2)
        self.assertEqual(evolution.result, [])
        system.set_property_number_of_molecules_evolution(evolution)
        self.assertEqual(evolution.result, [[0] * 4])

    def test_another_type_is_refused(self):
        system = co2_system()
        with self.assertRaises(TypeError):
            system.set_property_msd(raspalib.PropertyVolumeEvolution(sample_every=1))
        with self.assertRaises(TypeError):
            system.set_property_volume_evolution(1)

    def test_the_given_properties_do_not_keep_a_system_alive(self):
        import weakref

        system = co2_system()
        msd = raspalib.PropertyMeanSquaredDisplacement(sample_every=1)
        system.set_property_msd(msd)
        reference = weakref.ref(system)
        del system
        gc.collect()
        self.assertIsNone(reference())


class ConstructorsAskForWhatTheObjectUses(unittest.TestCase):
    def test_an_evolution_without_cycles_covers_all_phases_of_the_driver(self):
        system = co2_system()
        molecules = raspalib.PropertyNumberOfMoleculesEvolution(sample_every=2)
        volume = raspalib.PropertyVolumeEvolution()
        system.set_property_number_of_molecules_evolution(molecules)
        system.set_property_volume_evolution(volume)

        mc = raspalib.MonteCarlo(
            number_of_cycles=20,
            number_of_pre_initialization_cycles=4,
            number_of_initialization_cycles=10,
            number_of_equilibration_cycles=6,
            systems=[system],
            random_seed=12,
        )
        mc.run()

        self.assertEqual(molecules.result, [[40] * 20])
        self.assertEqual(volume.result, [27000.0] * 40)

    def test_cycles_that_are_given_are_kept(self):
        system = co2_system()
        volume = raspalib.PropertyVolumeEvolution(number_of_cycles=5, sample_every=1)
        system.set_property_volume_evolution(volume)
        raspalib.MonteCarlo(number_of_cycles=20, systems=[system], random_seed=12).run()
        self.assertEqual(volume.result, [27000.0] * 5)

    def test_the_keywords_of_the_release_are_still_taken(self):
        raspalib.PropertyNumberOfMoleculesEvolution(number_of_cycles=10, number_of_components=1, sample_every=1)
        raspalib.PropertyNumberOfMoleculesEvolution(number_of_production_cycles=10, sample_every=1)
        raspalib.PropertyNumberOfMoleculesEvolution(10, 1, 1)
        raspalib.PropertyVolumeEvolution(10, 1)
        raspalib.PropertyConservedEnergyEvolution(10, 1, None)
        raspalib.SampleMovie(system_id=0, sample_every=10)
        raspalib.WriteLammpsData(sample_every=10)
        raspalib.PropertyDensityGrid(
            number_of_frameworks=0,
            number_of_components=1,
            number_of_grid_points=(8, 8, 8),
            sample_every=1,
            write_every=100,
            density_grid_pseudo_atoms_list=[],
        )
        raspalib.PropertyEnergyHistogram(
            number_of_blocks=5, number_of_bins=10, value_range=(-100.0, 100.0), sample_every=1
        )
        raspalib.PropertyNumberOfMoleculesHistogram(number_of_blocks=5, value_range=(0, 50), sample_every=1)
        raspalib.PropertyRadialDistributionFunction(
            number_of_blocks=5, number_of_pseudo_atoms=2, number_of_bins=64, range=12.0, sample_every=1, write_every=100
        )

    def test_a_sampling_interval_of_zero_is_refused(self):
        for construct in (
            raspalib.PropertyNumberOfMoleculesEvolution,
            raspalib.PropertyVolumeEvolution,
            raspalib.PropertyConservedEnergyEvolution,
            raspalib.SampleMovie,
        ):
            with self.assertRaises(ValueError):
                construct(sample_every=0)

    def test_histograms_and_grids_are_sized_by_the_system(self):
        system = co2_system()
        energies = raspalib.PropertyEnergyHistogram(number_of_bins=20, value_range=(-20000.0, 0.0))
        molecules = raspalib.PropertyNumberOfMoleculesHistogram(value_range=(0, 50))
        rdf = raspalib.PropertyRadialDistributionFunction(number_of_bins=32, range=12.0, sample_every=1)
        grid = raspalib.PropertyDensityGrid(number_of_grid_points=(8, 8, 8), sample_every=1)
        system.set_average_energy_histogram(energies)
        system.set_number_of_molecules_histogram(molecules)
        system.set_property_rdf(rdf)
        system.set_property_density_grid(grid)

        mc = raspalib.MonteCarlo(
            number_of_cycles=20, number_of_initialization_cycles=10, systems=[system], random_seed=12
        )
        mc.run()

        _, average, _ = energies.result()
        self.assertGreater(sum(bin.total_energy for bin in average), 0.0)
        self.assertEqual(len(average), 20)
        _, histogram, _ = molecules.result()[0]
        self.assertEqual(max(range(50), key=lambda n: histogram[n]), 40)
        self.assertEqual(molecules.result(), mc.systems[0].property_number_of_molecules_histogram.result())

    def test_the_files_of_a_system_go_by_its_index_in_the_driver(self):
        previous = os.getcwd()
        with tempfile.TemporaryDirectory() as directory:
            os.chdir(directory)
            try:
                systems = [co2_system(10), co2_system(10)]
                writers = [raspalib.WriteLammpsData(sample_every=5) for _ in systems]
                for system, writer in zip(systems, writers):
                    self.assertIsNone(writer.system_id)
                    system.set_sample_pdb_movie(raspalib.SampleMovie(sample_every=5))
                    system.set_write_lammps_data(writer)
                self.assertFalse(os.path.exists("lammps"))

                raspalib.MonteCarlo(number_of_cycles=10, systems=systems, random_seed=12).run()

                self.assertEqual([writer.system_id for writer in writers], [0, 1])
                self.assertTrue(os.path.exists("lammps/s0.data"))
                self.assertTrue(os.path.exists("lammps/s1.data"))
                self.assertEqual(len(os.listdir("movies")), 2)
            finally:
                os.chdir(previous)


if __name__ == "__main__":
    unittest.main()
