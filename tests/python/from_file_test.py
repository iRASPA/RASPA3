"""A force field and a component are made from scratch or read from the files of the command-line program.

    python -m unittest discover -s tests/python -p "*_test.py"
"""

import json
import os
import pathlib
import tempfile
import unittest

os.environ.setdefault("RASPALIB_QUIET", "1")

import raspalib

EXAMPLE = pathlib.Path(__file__).resolve().parents[2] / "examples" / "basic" / "2_mc_co2_n2_in_two_independent_boxes"


def system_of(force_field: raspalib.ForceField, components: list[raspalib.Component], counts: list[int]):
    return raspalib.System(
        force_field=force_field,
        simulation_box=raspalib.SimulationBox(30.0, 30.0, 30.0),
        external_temperature=300.0,
        components=components,
        initial_number_of_molecules=counts,
    )


class FromFile(unittest.TestCase):
    def setUp(self):
        self.previous = os.getcwd()
        os.chdir(EXAMPLE)

    def tearDown(self):
        os.chdir(self.previous)

    def test_a_force_field_and_its_molecules_are_read(self):
        force_field = raspalib.ForceField.from_file("force_field.json")
        self.assertGreater(len(force_field.pseudo_atoms), 0)
        self.assertIsNotNone(force_field.find_pseudo_atom("C_co2"))

        moves = raspalib.MCMoveProbabilities(translation_probability=1.0, rotation_probability=1.0)
        co2 = raspalib.Component.from_file(force_field, "CO2.json", particle_probabilities=moves)
        n2 = raspalib.Component.from_file(force_field, "N2", component_name="nitrogen", particle_probabilities=moves)
        self.assertEqual(co2.name, "CO2")
        self.assertEqual(n2.name, "nitrogen")
        self.assertEqual(co2.type, raspalib.Component.Type.ADSORBATE)
        self.assertGreater(co2.critical_temperature, 0.0)

        system = system_of(force_field, [co2, n2], [20, 30])
        mc = raspalib.MonteCarlo(number_of_cycles=20, systems=[system], random_seed=12)
        mc.run()
        self.assertEqual(system.number_of_molecules_per_component, [20, 30])
        # the atoms carry the place of their component in the system, not an index given when the file was read
        self.assertEqual(sorted({atom.component_id for atom in system.atom_data}), [0, 1])
        self.assertNotEqual(system.average_energies.result()[0].total_energy.energy, 0.0)

    def test_a_file_that_is_missing_or_of_another_kind_is_refused(self):
        force_field = raspalib.ForceField.from_file("force_field.json")
        with self.assertRaises(ValueError):
            raspalib.ForceField.from_file("no_such_force_field.json")
        with self.assertRaises(ValueError):
            raspalib.ForceField.from_file("CO2.json")
        with self.assertRaises(ValueError):
            raspalib.Component.from_file(force_field, "no_such_molecule")
        with self.assertRaises(ValueError):
            raspalib.Component.from_file(force_field, "force_field.json")

    def test_a_molecule_of_another_force_field_is_refused(self):
        force_field = raspalib.ForceField.from_file("force_field.json")
        molecule = json.loads(pathlib.Path("CO2.json").read_text())
        text = json.dumps(molecule).replace("C_co2", "C_unknown")
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "other.json"
            path.write_text(text)
            with self.assertRaises(ValueError):
                raspalib.Component.from_file(force_field, str(path))

    def test_what_the_simulation_file_adds_is_set_on_the_objects(self):
        def energy(cutoff: float | None, use_charge: bool) -> float:
            force_field = raspalib.ForceField.from_file("force_field.json")
            if cutoff is not None:
                force_field.cutoff_molecule_vdw = cutoff
                force_field.cutoff_framework_vdw = cutoff
                self.assertEqual(force_field.cutoff_molecule_vdw, cutoff)
            force_field.use_charge = use_charge
            co2 = raspalib.Component.from_file(force_field, "CO2.json")
            mc = raspalib.MonteCarlo(number_of_cycles=1, systems=[system_of(force_field, [co2], [200])], random_seed=7)
            mc.setup()
            return mc.systems[0].compute_total_energies().potential_energy()

        reference = energy(None, True)
        self.assertNotEqual(reference, 0.0)
        self.assertNotAlmostEqual(energy(8.0, True), reference, places=3)
        self.assertNotAlmostEqual(energy(None, False), reference, places=3)
        self.assertAlmostEqual(energy(12.0, True), reference, places=6)

        force_field = raspalib.ForceField.from_file("force_field.json")
        with self.assertRaises(ValueError):
            force_field.cutoff_coulomb = 0.0


if __name__ == "__main__":
    unittest.main()
