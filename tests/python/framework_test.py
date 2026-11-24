"""Frameworks from CIF files, with the source of their charges, and flexible frameworks.

A framework is made from a CIF file or from the content of one; the zeolites that are built in need no file. A
framework is flexible when it is given a definition: potentials per combination of pseudo-atom types, which apply to
every such bond, bend or torsion of the framework, and groups of atoms that are fixed, rigid or flexible.

Run with pytest, or with the standard library:

    python -m unittest discover -s tests/python -p "*_test.py"
"""

import math
import os
import pathlib
import tempfile
import unittest

os.environ.setdefault("RASPALIB_QUIET", "1")

import raspalib

EXAMPLES = pathlib.Path(__file__).resolve().parents[2] / "examples" / "basic"
MFI_CIF = EXAMPLES / "8_mc_adsorption_of_methane_in_mfi" / "MFI_SI.cif"
CHAIN_DEFINITION = EXAMPLES / "22_md_mixed_framework_fixed_flexible" / "framework.json"

# three carbon atoms in a row, 1.5 Angstrom apart, in a box of 20 Angstrom
CHAIN_CIF = """data_chain
_cell_length_a    20.0000
_cell_length_b    20.0000
_cell_length_c    20.0000
_cell_angle_alpha 90.0000
_cell_angle_beta  90.0000
_cell_angle_gamma 90.0000
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
loop_
_symmetry_equiv_pos_as_xyz
 'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.40000 0.50000 0.50000
C2 C 0.47500 0.50000 0.50000
C3 C 0.55000 0.50000 0.50000
"""

# rock salt, with a charge on every atom in the file
SALT_CIF = """data_salt
_cell_length_a    5.6400
_cell_length_b    5.6400
_cell_length_c    5.6400
_cell_angle_alpha 90.0000
_cell_angle_beta  90.0000
_cell_angle_gamma 90.0000
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
loop_
_symmetry_equiv_pos_as_xyz
 'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_charge
Na Na 0.00000 0.00000 0.00000  0.80
Na Na 0.50000 0.50000 0.00000  0.80
Na Na 0.50000 0.00000 0.50000  0.80
Na Na 0.00000 0.50000 0.50000  0.80
Cl Cl 0.50000 0.50000 0.50000 -0.80
Cl Cl 0.00000 0.00000 0.50000 -0.80
Cl Cl 0.00000 0.50000 0.00000 -0.80
Cl Cl 0.50000 0.00000 0.00000 -0.80
"""


def zeolite_force_field() -> raspalib.ForceField:
    # zeolite: P. Bai, M. Tsapatsis, J. I. Siepmann, J. Phys. Chem. C 2013, 117, 24375-24387
    # CH4: M. G. Martin et al., J. Chem. Phys. 2001, 114, 7174-7181
    return raspalib.ForceField(
        pseudo_atoms=[
            raspalib.PseudoAtom(name="Si", framework_type=True, mass=28.0855, charge=2.05, atomic_number=14),
            raspalib.PseudoAtom(name="O", framework_type=True, mass=15.9994, charge=-1.025, atomic_number=8),
            raspalib.PseudoAtom(name="CH4", framework_type=False, mass=16.04246, charge=0.0, atomic_number=6),
        ],
        parameters=[
            raspalib.VDWParameters(22.0, 2.30),
            raspalib.VDWParameters(53.0, 3.30),
            raspalib.VDWParameters(158.5, 3.72),
        ],
        mixing_rule=raspalib.ForceField.MixingRule.LORENTZ_BERTHELOT,
        cutoff_framework_vdw=12.0,
        cutoff_molecule_vdw=12.0,
        cutoff_coulomb=12.0,
        shifted=True,
        tail_corrections=False,
        use_charge=True,
    )


def chain_force_field() -> raspalib.ForceField:
    return raspalib.ForceField(
        pseudo_atoms=[
            raspalib.PseudoAtom(name=name, framework_type=True, mass=12.0, charge=0.0, atomic_number=6)
            for name in ("C1", "C2", "C3")
        ],
        parameters=[raspalib.VDWParameters(50.0, 3.0)] * 3,
        mixing_rule=raspalib.ForceField.MixingRule.LORENTZ_BERTHELOT,
        cutoff_framework_vdw=10.0,
        cutoff_molecule_vdw=10.0,
        cutoff_coulomb=10.0,
        shifted=True,
        tail_corrections=False,
        use_charge=False,
    )


def salt_force_field() -> raspalib.ForceField:
    return raspalib.ForceField(
        pseudo_atoms=[
            raspalib.PseudoAtom(name="Na", framework_type=True, mass=22.98977, charge=1.0, atomic_number=11),
            raspalib.PseudoAtom(name="Cl", framework_type=True, mass=35.453, charge=-1.0, atomic_number=17),
        ],
        parameters=[raspalib.VDWParameters(15.0966, 2.65755), raspalib.VDWParameters(114.2126, 3.51932)],
        use_charge=True,
    )


def chain_definition() -> raspalib.FrameworkDefinition:
    """The definition of the example of a mixed fixed and flexible framework: the first atom is fixed."""
    return raspalib.FrameworkDefinition(
        bonds=[
            (("C1", "C2"), raspalib.BondType.HARMONIC, [50000.0, 1.5]),
            (("C2", "C3"), raspalib.BondType.HARMONIC, [50000.0, 1.5]),
        ],
        bends=[(("C1", "C2", "C3"), raspalib.BendType.HARMONIC, [20000.0, 180.0])],
        groups=[(raspalib.FrameworkGroupType.FIXED, [0]), (raspalib.FrameworkGroupType.FLEXIBLE, [1, 2])],
        exclude_intra_12_interactions=True,
        exclude_intra_13_interactions=True,
    )


def sorted_positions(atoms) -> list[tuple[float, float, float]]:
    return sorted(tuple(round(x, 3) for x in atom.position) for atom in atoms)


@unittest.skipUnless(MFI_CIF.exists(), "the MFI structure of the examples is not there")
class FrameworksFromFiles(unittest.TestCase):
    def test_a_framework_from_a_cif_file(self):
        force_field = zeolite_force_field()
        framework = raspalib.Framework.from_cif(
            force_field=force_field, file_name=str(MFI_CIF), number_of_unit_cells=(2, 2, 2)
        )
        self.assertEqual(framework.name, "MFI_SI")
        self.assertTrue(framework.rigid)
        self.assertEqual(framework.number_of_unit_cells, (2, 2, 2))
        # the unit cell of MFI has 96 silicon and 192 oxygen atoms
        self.assertEqual(len(framework.unit_cell_atoms), 288)
        self.assertEqual(framework.number_of_atoms, 8 * 288)
        self.assertEqual(len(framework.atoms), 8 * 288)
        self.assertAlmostEqual(framework.unit_cell_mass, 96 * 28.0855 + 192 * 15.9994, places=3)
        self.assertAlmostEqual(framework.mass, 8.0 * framework.unit_cell_mass, places=3)
        self.assertAlmostEqual(framework.net_charge, 0.0, places=6)
        box = framework.simulation_box
        self.assertEqual((box.length_a, box.length_b, box.length_c), (20.022, 19.899, 13.383))

    def test_the_extension_and_the_name(self):
        force_field = zeolite_force_field()
        without_extension = str(MFI_CIF)[: -len(".cif")]
        framework = raspalib.Framework.from_cif(
            force_field=force_field, file_name=without_extension, number_of_unit_cells=(1, 1, 1), component_name="MFI"
        )
        self.assertEqual(framework.name, "MFI")
        self.assertEqual(framework.number_of_atoms, 288)

    def test_the_content_of_a_file(self):
        force_field = zeolite_force_field()
        from_file = raspalib.Framework.from_cif(force_field, str(MFI_CIF), number_of_unit_cells=(1, 1, 2))
        from_content = raspalib.Framework.from_cif_string(
            force_field=force_field, cif=MFI_CIF.read_text(), component_name="MFI", number_of_unit_cells=(1, 1, 2)
        )
        self.assertEqual(from_content.name, "MFI")
        self.assertEqual(sorted_positions(from_content.atoms), sorted_positions(from_file.atoms))

    def test_without_a_number_of_unit_cells_the_cutoffs_fit(self):
        # a cutoff of 12 Angstrom needs 24 Angstrom in every direction: MFI is about 20 by 20 by 13
        framework = raspalib.Framework.from_cif(zeolite_force_field(), str(MFI_CIF))
        self.assertEqual(framework.number_of_unit_cells, (2, 2, 2))

    def test_the_file_and_the_built_in_zeolite_are_the_same_structure(self):
        force_field = zeolite_force_field()
        from_file = raspalib.Framework.from_cif(force_field, str(MFI_CIF), number_of_unit_cells=(2, 2, 2))
        built_in = raspalib.Framework.make_mfi(force_field=force_field, number_of_unit_cells=(2, 2, 2))
        self.assertEqual(built_in.number_of_atoms, from_file.number_of_atoms)
        for ours, theirs in zip(sorted_positions(built_in.atoms), sorted_positions(from_file.atoms)):
            for x, y in zip(ours, theirs):
                self.assertAlmostEqual(x, y, delta=0.01)

    def test_a_simulation_in_the_framework(self):
        force_field = zeolite_force_field()
        framework = raspalib.Framework.from_cif(force_field, str(MFI_CIF), number_of_unit_cells=(2, 2, 2))
        methane = raspalib.Component(
            force_field=force_field,
            component_name="CH4",
            critical_temperature=190.564,
            critical_pressure=4599200.0,
            acentric_factor=0.01142,
            defined_atoms=[raspalib.Atom(position=(0.0, 0.0, 0.0), charge=0.0, type=2)],
            particle_probabilities=raspalib.MCMoveProbabilities(
                translation_probability=0.5, reinsertion_cbmc_probability=1.0
            ),
        )
        system = raspalib.System(
            force_field=force_field,
            external_temperature=300.0,
            helium_void_fraction=0.3,
            framework_components=framework,
            components=[methane],
            initial_number_of_molecules=[1],
        )
        mc = raspalib.MonteCarlo(
            number_of_cycles=200, number_of_initialization_cycles=100, systems=[system], random_seed=12
        )
        mc.run()
        # the enthalpy of adsorption of methane in MFI at 300 K is about -20 kJ/mol: -2400 K
        energy, _ = system.average_energies.result()
        self.assertLess(energy.total_energy.energy, -1000.0)
        self.assertGreater(energy.total_energy.energy, -4000.0)


class BuiltInZeolites(unittest.TestCase):
    def test_the_zeolites_that_need_no_file(self):
        force_field = zeolite_force_field()
        # (chabazite in its rhombohedral cell)
        atoms_per_unit_cell = {"make_mfi": 288, "make_itq_29": 72, "make_fau": 576, "make_cha": 36}
        for maker, number_of_atoms in atoms_per_unit_cell.items():
            framework = getattr(raspalib.Framework, maker)(force_field=force_field)
            self.assertEqual(framework.number_of_atoms, number_of_atoms, maker)
            self.assertAlmostEqual(framework.net_charge, 0.0, places=6, msg=maker)
            doubled = getattr(raspalib.Framework, maker)(force_field=force_field, number_of_unit_cells=(1, 1, 2))
            self.assertEqual(doubled.number_of_atoms, 2 * number_of_atoms, maker)

    def test_they_need_the_pseudo_atoms_of_a_zeolite(self):
        with self.assertRaises(ValueError):
            raspalib.Framework.make_mfi(force_field=chain_force_field())


class Charges(unittest.TestCase):
    def framework(self, use_charges_from) -> raspalib.Framework:
        return raspalib.Framework.from_cif_string(
            force_field=salt_force_field(),
            cif=SALT_CIF,
            component_name="salt",
            number_of_unit_cells=(1, 1, 1),
            use_charges_from=use_charges_from,
        )

    def test_the_charges_of_the_pseudo_atoms(self):
        framework = self.framework(raspalib.Framework.UseChargesFrom.PSEUDO_ATOMS)
        self.assertEqual(sorted(atom.charge for atom in framework.atoms), [-1.0] * 4 + [1.0] * 4)

    def test_the_charges_of_the_file(self):
        framework = self.framework(raspalib.Framework.UseChargesFrom.CIF_FILE)
        self.assertEqual(sorted(atom.charge for atom in framework.atoms), [-0.8] * 4 + [0.8] * 4)
        self.assertAlmostEqual(framework.net_charge, 0.0, places=9)

    def test_charge_equilibration(self):
        framework = self.framework(raspalib.Framework.UseChargesFrom.CHARGE_EQUILIBRATION)
        charges = sorted(atom.charge for atom in framework.atoms)
        self.assertAlmostEqual(sum(charges), 0.0, places=6)
        # the chlorine atoms take charge from the sodium atoms, all four alike
        self.assertTrue(all(charge < -0.05 for charge in charges[:4]), charges)
        self.assertTrue(all(charge > 0.05 for charge in charges[4:]), charges)
        self.assertAlmostEqual(charges[0], charges[3], places=6)
        self.assertNotIn(1.0, charges)


class FlexibleFrameworks(unittest.TestCase):
    def framework(self, **keywords) -> raspalib.Framework:
        return raspalib.Framework.from_cif_string(
            force_field=chain_force_field(),
            cif=CHAIN_CIF,
            component_name="chain",
            number_of_unit_cells=(1, 1, 1),
            **keywords,
        )

    def test_without_a_definition_the_framework_is_rigid(self):
        framework = self.framework()
        self.assertTrue(framework.rigid)
        self.assertEqual(framework.number_of_mobile_atoms, 0)
        self.assertEqual(framework.intra_molecular_potentials.bonds, [])

    def test_a_definition_gives_every_bond_and_bend_its_potential(self):
        framework = self.framework(definition=chain_definition())
        self.assertFalse(framework.rigid)
        self.assertEqual(framework.number_of_atoms, 3)
        self.assertEqual(framework.number_of_fixed_atoms, 1)
        self.assertEqual(framework.number_of_mobile_atoms, 2)
        self.assertEqual(framework.connectivity_table.bonds, [[0, 1], [1, 2]])
        potentials = framework.intra_molecular_potentials
        self.assertEqual([bond.identifiers for bond in potentials.bonds], [[0, 1], [1, 2]])
        self.assertEqual({bond.type for bond in potentials.bonds}, {raspalib.BondType.HARMONIC})
        self.assertEqual([bend.identifiers for bend in potentials.bends], [[0, 1, 2]])

    def run_md(self, framework: raspalib.Framework) -> raspalib.System:
        system = raspalib.System(
            force_field=chain_force_field(),
            external_temperature=300.0,
            framework_components=framework,
        )
        system.time_step = 0.0005
        system.set_thermostat(raspalib.Thermostat())
        md = raspalib.MolecularDynamics(
            number_of_cycles=400, number_of_equilibration_cycles=200, systems=[system], random_seed=12
        )
        md.run()
        return system

    def test_molecular_dynamics_moves_the_flexible_atoms_only(self):
        system = self.run_md(self.framework(definition=chain_definition()))
        fixed, middle, end = (atom.position for atom in system.atom_data)
        # the first atom is where the file puts it; the others have moved, and are still bonded
        self.assertEqual(tuple(round(x, 9) for x in fixed), (8.0, 10.0, 10.0))
        self.assertGreater(math.dist(middle, (9.5, 10.0, 10.0)), 1.0e-4)
        self.assertGreater(math.dist(end, (11.0, 10.0, 10.0)), 1.0e-4)
        self.assertAlmostEqual(math.dist(fixed, middle), 1.5, delta=0.2)
        self.assertAlmostEqual(math.dist(middle, end), 1.5, delta=0.2)

    @unittest.skipUnless(CHAIN_DEFINITION.exists(), "the framework definition of the examples is not there")
    def test_a_definition_file_and_the_same_definition_in_python_are_the_same_simulation(self):
        from_file = self.framework(definition_file=str(CHAIN_DEFINITION))
        self.assertFalse(from_file.rigid)
        self.assertEqual(from_file.number_of_fixed_atoms, 1)

        in_python = self.run_md(self.framework(definition=chain_definition()))
        read = self.run_md(from_file)
        self.assertEqual(
            [atom.position for atom in in_python.atom_data], [atom.position for atom in read.atom_data]
        )

    def test_a_framework_that_is_defined_atom_by_atom(self):
        framework = raspalib.Framework(
            force_field=chain_force_field(),
            component_name="chain",
            simulation_box=raspalib.SimulationBox(a=20.0, b=20.0, c=20.0),
            space_group_hall_number=1,
            defined_atoms=[
                raspalib.Atom(position=(0.400, 0.5, 0.5), charge=0.0, type=0),
                raspalib.Atom(position=(0.475, 0.5, 0.5), charge=0.0, type=1),
                raspalib.Atom(position=(0.550, 0.5, 0.5), charge=0.0, type=2),
            ],
            number_of_unit_cells=(1, 1, 1),
            definition=chain_definition(),
        )
        self.assertFalse(framework.rigid)
        self.assertEqual(framework.number_of_mobile_atoms, 2)
        self.assertEqual(len(framework.intra_molecular_potentials.bonds), 2)

    def test_what_is_refused(self):
        with self.assertRaisesRegex(ValueError, "not both"):
            self.framework(definition=chain_definition(), definition_file="framework")
        with self.assertRaises(ValueError):
            self.framework(definition_file="no_such_definition")
        # a potential between pseudo-atoms the force field does not have
        unknown = raspalib.FrameworkDefinition(bonds=[(("C1", "C9"), raspalib.BondType.HARMONIC, [50000.0, 1.5])])
        with self.assertRaisesRegex(ValueError, "C9"):
            self.framework(definition=unknown)
        # groups are of a flexible framework, and name atoms the framework has
        rigid_with_groups = raspalib.FrameworkDefinition(
            flexible=False, groups=[(raspalib.FrameworkGroupType.FIXED, [0])]
        )
        with self.assertRaises(ValueError):
            self.framework(definition=rigid_with_groups)
        beyond = raspalib.FrameworkDefinition(groups=[(raspalib.FrameworkGroupType.FIXED, [7])])
        with self.assertRaises(ValueError):
            self.framework(definition=beyond)


class Refusals(unittest.TestCase):
    def test_what_is_refused(self):
        force_field = chain_force_field()
        with self.assertRaises(ValueError):
            raspalib.Framework.from_cif(force_field, "no_such_file.cif")
        # the atoms of the file have no pseudo-atom in this force field
        with self.assertRaises(ValueError):
            raspalib.Framework.from_cif_string(force_field, SALT_CIF, "salt", number_of_unit_cells=(1, 1, 1))
        with self.assertRaises(ValueError):
            raspalib.Framework.from_cif_string(force_field, CHAIN_CIF, "chain", number_of_unit_cells=(0, 1, 1))
        with self.assertRaises(ValueError):
            raspalib.Framework.from_cif_string(force_field, "this is not a cif file", "nothing")

    def test_a_file_in_another_directory_is_found_by_its_path(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "chain.cif"
            path.write_text(CHAIN_CIF)
            framework = raspalib.Framework.from_cif(chain_force_field(), str(path), number_of_unit_cells=(1, 1, 1))
            self.assertEqual(framework.name, "chain")
            self.assertEqual(framework.number_of_atoms, 3)


if __name__ == "__main__":
    unittest.main()
