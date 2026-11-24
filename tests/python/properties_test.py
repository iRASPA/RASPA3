"""The analyses of molecules, the elastic constants from fluctuations, and the remaining samplers.

A component asks for an analysis of its molecules with settings; the driver creates the analyses of its systems, and
the results are read from the system afterwards. The other properties are set on a system like the ones before.

Run with pytest, or with the standard library:

    python -m unittest discover -s tests/python -p "*_test.py"
"""

import math
import os
import tempfile
import unittest

os.environ.setdefault("RASPALIB_QUIET", "1")

import raspalib

from potentials_test import alkane, alkane_force_field
from settings_test import MolecularDynamicsEnsembles


def mean_of(histogram) -> float:
    values, density, _ = histogram
    return sum(value * weight for value, weight in zip(values, density)) / sum(density)


def area_of(histogram) -> float:
    values, density, _ = histogram
    return sum(density) * (values[1] - values[0])


def hexane(**settings) -> raspalib.Component:
    component = alkane(alkane_force_field(), 6)
    component.mc_moves_probabilities[raspalib.Move.Types.TRANSLATION] = 1.0
    component.mc_moves_probabilities[raspalib.Move.Types.REINSERTION_CBMC] = 1.0
    for name, value in settings.items():
        setattr(component, name, value)
    return component


def hexane_box(component: raspalib.Component) -> raspalib.System:
    return raspalib.System(
        force_field=alkane_force_field(),
        simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
        external_temperature=300.0,
        components=[component],
        initial_number_of_molecules=[20],
    )


def all_settings() -> dict:
    return dict(
        molecule_properties_settings=raspalib.MoleculePropertiesSettings(sample_every=1),
        molecule_shape_settings=raspalib.MoleculeShapeSettings(sample_every=1),
        molecule_backbone_settings=raspalib.MoleculeBackboneSettings(sample_every=1, number_of_wave_vectors=8),
        end_to_end_acf_settings=raspalib.EndToEndACFSettings(sample_every=1),
    )


class MoleculeAnalyses(unittest.TestCase):
    """Hexane in Monte Carlo, with every analysis of its molecules."""

    @classmethod
    def setUpClass(cls):
        cls.system = hexane_box(hexane(**all_settings()))
        cls.before = (
            cls.system.property_molecule_properties,
            cls.system.property_molecule_shape,
            cls.system.property_molecule_backbone,
            cls.system.property_end_to_end_acf,
        )
        mc = raspalib.MonteCarlo(
            number_of_cycles=400, number_of_initialization_cycles=100, systems=[cls.system], random_seed=5
        )
        mc.run()

    def test_the_driver_creates_the_analyses(self):
        self.assertEqual(self.before, (None, None, None, None))
        self.assertIsNotNone(self.system.property_molecule_properties)
        self.assertIsNotNone(self.system.property_molecule_shape)
        self.assertIsNotNone(self.system.property_molecule_backbone)
        self.assertIsNotNone(self.system.property_end_to_end_acf)

    def test_the_histograms_of_the_bonded_terms(self):
        properties = self.system.property_molecule_properties
        self.assertEqual(properties.number_of_bonds, [5])
        self.assertEqual(properties.number_of_bends, [4])
        self.assertEqual(properties.number_of_torsions, [3])

        for bond in range(5):
            histogram = properties.bond_histogram(component=0, bond=bond)
            self.assertAlmostEqual(area_of(histogram), 1.0, places=9)
            self.assertAlmostEqual(mean_of(histogram), 1.54, delta=0.03)
        for bend in range(4):
            histogram = properties.bend_histogram(component=0, bend=bend)
            self.assertAlmostEqual(area_of(histogram), 1.0, places=9)
            self.assertAlmostEqual(mean_of(histogram), 114.0, delta=2.0)
        for torsion in range(3):
            values, density, error = properties.torsion_histogram(component=0, torsion=torsion)
            self.assertEqual(len(values), 128)
            self.assertGreater(values[0], -180.0)
            self.assertLess(values[-1], 180.0)
            self.assertAlmostEqual(values[0], -values[-1], places=9)
            self.assertAlmostEqual(area_of((values, density, error)), 1.0, places=9)
            # the trans conformation, at plus or minus 180 degrees, is the most likely one
            self.assertGreater(density[0] + density[-1], density[len(density) // 2])

        end_to_end = properties.end_to_end_histogram(component=0)
        self.assertAlmostEqual(area_of(end_to_end), 1.0, places=9)
        # between the folded chain and the stretched one of five bonds
        self.assertGreater(mean_of(end_to_end), 4.0)
        self.assertLess(mean_of(end_to_end), 5 * 1.54)

    def test_the_shape_of_the_molecules(self):
        shape = self.system.property_molecule_shape
        moments = raspalib.PropertyMoleculeShape.Moment
        radius, radius_error = shape.result(component=0, moment=moments.RADIUS_OF_GYRATION)
        radius_squared, _ = shape.result(0, moments.RADIUS_OF_GYRATION_SQUARED)
        self.assertGreater(radius, 1.5)
        self.assertLess(radius, 3.0)
        self.assertGreater(radius_error, 0.0)
        # the mean of the square is at least the square of the mean
        self.assertGreaterEqual(radius_squared, radius**2)
        self.assertAlmostEqual(mean_of(shape.radius_of_gyration_histogram(0)), radius, delta=0.06)

        # the eigenvalues of the gyration tensor are ordered and add up to the squared radius of gyration
        eigenvalues = [shape.result(0, moment)[0] for moment in (moments.LAMBDA_1, moments.LAMBDA_2, moments.LAMBDA_3)]
        self.assertEqual(eigenvalues, sorted(eigenvalues, reverse=True))
        self.assertAlmostEqual(sum(eigenvalues), radius_squared, places=6)
        trace = sum(shape.result(0, moment)[0] for moment in (moments.TENSOR_XX, moments.TENSOR_YY, moments.TENSOR_ZZ))
        self.assertAlmostEqual(trace, radius_squared, places=6)

        anisotropy, _ = shape.result(0, moments.SHAPE_ANISOTROPY)
        self.assertGreater(anisotropy, 0.0)
        self.assertLess(anisotropy, 1.0)
        self.assertAlmostEqual(area_of(shape.shape_anisotropy_histogram(0)), 1.0, places=9)
        self.assertAlmostEqual(area_of(shape.prolateness_histogram(0)), 1.0, places=9)

    def test_the_backbone(self):
        backbone = self.system.property_molecule_backbone
        moments = raspalib.PropertyMoleculeBackbone.Moment
        self.assertEqual(backbone.backbone(0), [0, 1, 2, 3, 4, 5])
        self.assertAlmostEqual(backbone.contour_length(0), 5 * 1.54, delta=0.1)

        bond_length, _ = backbone.moment(component=0, moment=moments.BOND_LENGTH)
        self.assertAlmostEqual(bond_length, 1.54, delta=0.03)

        averages = backbone.result(0)
        # beads k bonds apart, for k = 1 to 5: further apart along the chain is further apart in space
        distances = averages.internal_distance_squared
        self.assertEqual(len(distances), 5)
        self.assertEqual(distances, sorted(distances))
        self.assertAlmostEqual(distances[0], backbone.moment(0, moments.BOND_LENGTH_SQUARED)[0], places=6)
        self.assertAlmostEqual(distances[-1], backbone.moment(0, moments.END_TO_END_SQUARED)[0], places=6)
        self.assertAlmostEqual(averages.bond_correlation[0], 1.0, places=9)

        wave_vectors = backbone.wave_vectors(0)
        self.assertEqual(len(wave_vectors), 8)
        self.assertAlmostEqual(wave_vectors[0], 0.01, places=9)
        self.assertAlmostEqual(wave_vectors[-1], 5.0, places=9)
        self.assertEqual(len(averages.form_factor), 8)
        # the form factor goes to one at long wavelengths and decays from there
        self.assertAlmostEqual(averages.form_factor[0], 1.0, places=3)
        self.assertLess(averages.form_factor[-1], 0.5)

        length, error = backbone.persistence_length_from_projection(0)
        self.assertGreater(length, 1.54)
        self.assertGreater(error, 0.0)

    def test_the_end_to_end_autocorrelation_function(self):
        function = self.system.property_end_to_end_acf
        self.assertEqual(function.time_unit, "cycles")
        self.assertTrue(function.has_data(0))
        points = function.result(component=0)
        self.assertGreater(len(points), 10)
        self.assertEqual([point.time for point in points[:3]], [0.0, 1.0, 2.0])
        self.assertAlmostEqual(points[0].normalized, 1.0, places=9)
        self.assertLess(points[10].normalized, points[0].normalized)
        self.assertGreater(points[0].number_of_samples, 0)

        times = function.relaxation_times(0)
        self.assertAlmostEqual(times.mean_squared_end_to_end, points[0].acf, places=6)
        self.assertGreater(times.longest_lag, 0.0)

    def test_the_three_analyses_agree_on_the_end_to_end_distance(self):
        # each of them keeps its own sums of the squared end-to-end distance
        shape, _ = self.system.property_molecule_shape.result(0, raspalib.PropertyMoleculeShape.Moment.END_TO_END_SQUARED)
        backbone, _ = self.system.property_molecule_backbone.moment(
            0, raspalib.PropertyMoleculeBackbone.Moment.END_TO_END_SQUARED
        )
        function = self.system.property_end_to_end_acf.relaxation_times(0).mean_squared_end_to_end
        self.assertAlmostEqual(shape, backbone, places=6)
        self.assertAlmostEqual(shape, function, places=6)
        # and the histogram has it to the width of its bins
        mean = mean_of(self.system.property_molecule_properties.end_to_end_histogram(0))
        self.assertLessEqual(mean**2, shape + 0.5)

    def test_what_is_not_there(self):
        properties = self.system.property_molecule_properties
        with self.assertRaises(IndexError):
            properties.bond_histogram(component=0, bond=5)
        with self.assertRaises(IndexError):
            properties.bend_histogram(component=1, bend=0)
        with self.assertRaises(IndexError):
            self.system.property_molecule_shape.radius_of_gyration_histogram(3)
        with self.assertRaises(IndexError):
            self.system.property_molecule_backbone.result(3)
        self.assertEqual(self.system.property_end_to_end_acf.result(3), [])
        with self.assertRaises(IndexError):
            self.system.property_end_to_end_acf.relaxation_times(3)


class AnalysesAndDrivers(unittest.TestCase):
    def test_only_what_is_asked_for_is_created(self):
        system = hexane_box(hexane(molecule_shape_settings=raspalib.MoleculeShapeSettings(sample_every=2)))
        mc = raspalib.MonteCarlo(number_of_cycles=20, number_of_initialization_cycles=10, systems=[system], random_seed=5)
        mc.run()
        self.assertIsNone(system.property_molecule_properties)
        self.assertIsNotNone(system.property_molecule_shape)
        self.assertIsNone(system.property_molecule_backbone)
        self.assertIsNone(system.property_end_to_end_acf)

    def test_without_settings_nothing_is_created(self):
        system = hexane_box(hexane())
        mc = raspalib.MonteCarlo(number_of_cycles=20, number_of_initialization_cycles=10, systems=[system], random_seed=5)
        mc.run()
        self.assertIsNone(system.property_molecule_properties)
        self.assertIsNone(system.property_end_to_end_acf)

    def test_molecular_dynamics_has_its_lags_in_picoseconds(self):
        system = hexane_box(hexane(end_to_end_acf_settings=raspalib.EndToEndACFSettings(sample_every=1)))
        system.time_step = 0.0005
        system.set_thermostat(raspalib.Thermostat())
        md = raspalib.MolecularDynamics(
            number_of_cycles=200,
            number_of_initialization_cycles=50,
            number_of_equilibration_cycles=50,
            systems=[system],
            random_seed=5,
        )
        md.run()
        function = system.property_end_to_end_acf
        self.assertEqual(function.time_unit, "ps")
        points = function.result(0)
        self.assertAlmostEqual(points[1].time, 0.0005, places=12)
        self.assertAlmostEqual(points[2].time, 0.0010, places=12)

    def test_an_analysis_a_molecule_cannot_have_is_refused(self):
        # a rigid molecule has no end-to-end atoms
        force_field = alkane_force_field()
        ethane = raspalib.Component(
            force_field=force_field,
            component_name="ethane",
            critical_temperature=305.33,
            critical_pressure=4871800.0,
            acentric_factor=0.0993,
            defined_atoms=[
                raspalib.Atom(position=(0.0, 0.0, 0.77), charge=0.0, type=0),
                raspalib.Atom(position=(0.0, 0.0, -0.77), charge=0.0, type=0),
            ],
        )
        ethane.end_to_end_acf_settings = raspalib.EndToEndACFSettings()
        system = raspalib.System(
            force_field=force_field,
            simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
            external_temperature=300.0,
            components=[ethane],
            initial_number_of_molecules=[5],
        )
        with self.assertRaises(ValueError):
            raspalib.MonteCarlo(number_of_cycles=10, systems=[system])

    def test_the_settings(self):
        properties = raspalib.MoleculePropertiesSettings(
            sample_every=5, write_every=100, number_of_bins=64, bond_range=3.0, end_to_end_range=20.0, end_to_end_bin_width=0.5
        )
        self.assertEqual(
            (properties.sample_every, properties.write_every, properties.number_of_bins),
            (5, 100, 64),
        )
        self.assertEqual((properties.bond_range, properties.end_to_end_range, properties.end_to_end_bin_width), (3.0, 20.0, 0.5))
        shape = raspalib.MoleculeShapeSettings(mass_weighted=True, radius_of_gyration_range=12.0)
        self.assertEqual((shape.mass_weighted, shape.radius_of_gyration_range, shape.write_every), (True, 12.0, None))
        backbone = raspalib.MoleculeBackboneSettings(wave_vector_lower_limit=0.1, wave_vector_upper_limit=2.0)
        self.assertEqual((backbone.wave_vector_lower_limit, backbone.wave_vector_upper_limit), (0.1, 2.0))
        self.assertEqual(raspalib.EndToEndACFSettings(number_of_block_elements=10).number_of_block_elements, 10)

        component = hexane(molecule_shape_settings=shape)
        self.assertTrue(component.molecule_shape_settings.mass_weighted)
        component.molecule_shape_settings = None
        self.assertIsNone(component.molecule_shape_settings)

        for refused in (
            lambda: raspalib.MoleculePropertiesSettings(sample_every=0),
            lambda: raspalib.MoleculePropertiesSettings(end_to_end_bin_width=0.0),
            lambda: raspalib.MoleculeShapeSettings(radius_of_gyration_bin_width=-0.1),
            lambda: raspalib.MoleculeBackboneSettings(wave_vector_lower_limit=0.0),
            lambda: raspalib.MoleculeBackboneSettings(wave_vector_lower_limit=2.0, wave_vector_upper_limit=1.0),
            lambda: raspalib.EndToEndACFSettings(number_of_block_elements=1),
        ):
            with self.assertRaises(ValueError):
                refused()


class ElasticConstants(unittest.TestCase):
    def fluid(self) -> raspalib.System:
        return MolecularDynamicsEnsembles().co2_fluid(pressure=None)

    def test_the_elastic_constants_of_a_fluid_at_constant_temperature(self):
        system = self.fluid()
        system.set_thermostat(raspalib.Thermostat())
        system.set_property_elastic_constants_fluctuation(
            raspalib.PropertyElasticConstantsFluctuation(number_of_blocks=5), sample_every=5
        )
        md = raspalib.MolecularDynamics(
            number_of_cycles=400,
            number_of_initialization_cycles=100,
            number_of_equilibration_cycles=100,
            systems=[system],
            random_seed=12,
        )
        md.run()

        result = system.property_elastic_constants_fluctuation.result()
        self.assertEqual(result["voigtOrder"], ["xx", "yy", "zz", "yz", "xz", "xy"])
        for tensor in ("born", "kinetic", "fluctuation", "stiffness", "tangentStiffness"):
            self.assertEqual(len(result[tensor]), 36, tensor)
            self.assertTrue(all(math.isfinite(value) for value in result[tensor]), tensor)
        self.assertEqual(len(result["configurationalStress"]), 6)
        # the Born term of the compression along an axis is positive
        self.assertGreater(result["born"][0], 0.0)
        self.assertIn("Samples: 80", system.property_elastic_constants_fluctuation.write_averages_statistics())

    def test_they_are_of_molecular_dynamics_at_constant_temperature(self):
        system = self.fluid()
        system.set_property_elastic_constants_fluctuation(raspalib.PropertyElasticConstantsFluctuation(), sample_every=5)
        md = raspalib.MolecularDynamics(number_of_cycles=10, systems=[system], random_seed=12)
        with self.assertRaises(RuntimeError):
            md.run()

    def test_setting_and_removing(self):
        system = self.fluid()
        self.assertIsNone(system.property_elastic_constants_fluctuation)
        system.set_property_elastic_constants_fluctuation(raspalib.PropertyElasticConstantsFluctuation())
        self.assertIsNotNone(system.property_elastic_constants_fluctuation)
        system.set_property_elastic_constants_fluctuation(None)
        self.assertIsNone(system.property_elastic_constants_fluctuation)
        with self.assertRaises(ValueError):
            system.set_property_elastic_constants_fluctuation(raspalib.PropertyElasticConstantsFluctuation(), sample_every=0)


class OtherSamplers(unittest.TestCase):
    def test_lammps_data_files_are_written(self):
        previous_directory = os.getcwd()
        with tempfile.TemporaryDirectory() as directory:
            os.chdir(directory)
            try:
                system = hexane_box(hexane())
                self.assertIsNone(system.write_lammps_data)
                system.set_write_lammps_data(raspalib.WriteLammpsData(system_id=0, sample_every=10))
                self.assertEqual(system.write_lammps_data.sample_every, 10)
                mc = raspalib.MonteCarlo(
                    number_of_cycles=20, number_of_initialization_cycles=10, systems=[system], random_seed=5
                )
                mc.run()
                self.assertTrue(os.path.exists(os.path.join("lammps", "s0.data")))
                self.assertGreater(os.path.getsize(os.path.join("lammps", "s0.data")), 0)

                system.set_write_lammps_data(None)
                self.assertIsNone(system.write_lammps_data)
            finally:
                os.chdir(previous_directory)
        with self.assertRaises(ValueError):
            raspalib.WriteLammpsData(sample_every=0)

    def test_the_waiting_times_of_radicals(self):
        system = hexane_box(hexane())
        self.assertIsNone(system.property_radical_waiting_times)
        system.set_property_radical_waiting_times(
            raspalib.PropertyRadicalWaitingTimes(bins_per_decade=5, write_every=100)
        )
        waiting_times = system.property_radical_waiting_times
        self.assertEqual((waiting_times.bins_per_decade, waiting_times.write_every), (5, 100))
        # nothing reacts in this system: the histogram is there and empty, from a femtosecond to a microsecond
        histogram = waiting_times.closed()
        self.assertEqual(histogram.count, 0)
        self.assertEqual(len(histogram.counts), 9 * 5)
        edges = histogram.bin_edges()
        self.assertEqual(len(edges), 9 * 5 + 1)
        self.assertAlmostEqual(edges[0], 1.0e-3, places=12)
        self.assertAlmostEqual(edges[-1], 1.0e6, places=3)

        system.set_property_radical_waiting_times(None)
        self.assertIsNone(system.property_radical_waiting_times)
        with self.assertRaises(ValueError):
            raspalib.PropertyRadicalWaitingTimes(bins_per_decade=0)


if __name__ == "__main__":
    unittest.main()
