"""Pair potentials other than Lennard-Jones, pairs given explicitly, and flexible molecules.

A flexible molecule is defined by its bonded terms, in the units of the molecule files: energies in Kelvin, lengths
in Angstrom, angles in degrees. The bonds are the connectivity of the molecule; every bond, bend and torsion they
form needs a potential.

Run with pytest, or with the standard library:

    python -m unittest discover -s tests/python -p "*_test.py"
"""

import math
import os
import unittest

os.environ.setdefault("RASPALIB_QUIET", "1")

import raspalib

# the energies of 'System.compute_total_energies' are in internal units of 10 J/mol (to six digits here)
KELVIN_TO_ENERGY = 0.831446261815324


def pseudo_atoms(*names: str) -> list[raspalib.PseudoAtom]:
    return [
        raspalib.PseudoAtom(name=name, framework_type=False, mass=16.0, charge=0.0, atomic_number=6) for name in names
    ]


def pair_energy(force_field: raspalib.ForceField, types: tuple[int, int], distance: float) -> float:
    """The potential energy [K] of two single-atom molecules at the given distance."""
    components = [
        raspalib.Component(
            force_field=force_field,
            component_name=f"atom_{index}",
            critical_temperature=190.564,
            critical_pressure=4599200.0,
            acentric_factor=0.01142,
            defined_atoms=[raspalib.Atom(position=(0.0, 0.0, 0.0), charge=0.0, type=type_)],
        )
        for index, type_ in enumerate(types)
    ]
    system = raspalib.System(
        force_field=force_field,
        simulation_box=raspalib.SimulationBox(a=40.0, b=40.0, c=40.0),
        external_temperature=300.0,
        components=components,
        # the molecules at these positions come on top of 'initial_number_of_molecules', which places at random
        initial_positions=[[(10.0, 10.0, 10.0)], [(10.0 + distance, 10.0, 10.0)]],
        initial_number_of_molecules=[0, 0],
    )
    assert system.number_of_molecules_per_component == [1, 1]
    return system.compute_total_energies().potential_energy() / KELVIN_TO_ENERGY


def lennard_jones(epsilon: float, sigma: float, distance: float) -> float:
    return 4.0 * epsilon * ((sigma / distance) ** 12 - (sigma / distance) ** 6)


class PairPotentials(unittest.TestCase):
    def assert_energy(self, energy: float, expected: float):
        # the conversion from the internal units is known here to six digits
        self.assertTrue(math.isclose(energy, expected, rel_tol=1.0e-5, abs_tol=1.0e-6), f"{energy} != {expected}")

    def force_field(self, parameters, binary_interactions=()) -> raspalib.ForceField:
        return raspalib.ForceField(
            pseudo_atoms=pseudo_atoms("A", "B"),
            parameters=parameters,
            mixing_rule=raspalib.ForceField.MixingRule.LORENTZ_BERTHELOT,
            shifted=False,
            tail_corrections=False,
            use_charge=False,
            binary_interactions=list(binary_interactions),
        )

    def test_lennard_jones_with_the_mixing_rule(self):
        force_field = self.force_field([raspalib.VDWParameters(100.0, 3.0), raspalib.VDWParameters(400.0, 4.0)])
        self.assert_energy(pair_energy(force_field, (0, 0), 3.5), lennard_jones(100.0, 3.0, 3.5))
        self.assert_energy(pair_energy(force_field, (1, 1), 4.5), lennard_jones(400.0, 4.0, 4.5))
        # Lorentz-Berthelot: the geometric mean of the strengths, the arithmetic mean of the sizes
        self.assert_energy(pair_energy(force_field, (0, 1), 4.0), lennard_jones(200.0, 3.5, 4.0))

    def test_a_potential_by_type_and_parameters(self):
        types = raspalib.VDWParameters.Type
        lennard_jones_by_type = raspalib.VDWParameters(type=types.LENNARD_JONES, parameters=[100.0, 3.0])
        self.assertEqual(lennard_jones_by_type.type, types.LENNARD_JONES)
        force_field = self.force_field([lennard_jones_by_type, raspalib.VDWParameters(400.0, 4.0)])
        self.assert_energy(pair_energy(force_field, (0, 0), 3.5), lennard_jones(100.0, 3.0, 3.5))

    def test_a_buckingham_potential(self):
        # u(r) = A exp(-B r) - C / r^6, with A and C in Kelvin
        a, b, c = 1.0e7, 3.5, 2.0e5
        buckingham = raspalib.VDWParameters(type=raspalib.VDWParameters.Type.BUCKINGHAM, parameters=[a, b, c])
        self.assertEqual(buckingham.type, raspalib.VDWParameters.Type.BUCKINGHAM)
        force_field = self.force_field([buckingham, raspalib.VDWParameters(400.0, 4.0)])

        for distance in (3.0, 4.0, 6.0):
            expected = a * math.exp(-b * distance) - c / distance**6
            self.assert_energy(pair_energy(force_field, (0, 0), distance), expected)
        # a form without a mixing rule does not interact with another unless the pair is given
        self.assertEqual(pair_energy(force_field, (0, 1), 4.0), 0.0)

    def test_a_pair_that_is_given_replaces_the_mixing_rule(self):
        selfs = [raspalib.VDWParameters(100.0, 3.0), raspalib.VDWParameters(400.0, 4.0)]
        force_field = self.force_field(selfs, [("A", "B", raspalib.VDWParameters(50.0, 3.2))])
        self.assert_energy(pair_energy(force_field, (0, 1), 4.0), lennard_jones(50.0, 3.2, 4.0))
        self.assert_energy(pair_energy(force_field, (1, 0), 4.0), lennard_jones(50.0, 3.2, 4.0))
        # the pairs that are not given keep theirs
        self.assert_energy(pair_energy(force_field, (0, 0), 3.5), lennard_jones(100.0, 3.0, 3.5))

        a, b, c = 1.0e7, 3.5, 2.0e5
        cross = raspalib.VDWParameters(type=raspalib.VDWParameters.Type.BUCKINGHAM, parameters=[a, b, c])
        force_field = self.force_field(selfs, [("B", "A", cross)])
        self.assert_energy(pair_energy(force_field, (0, 1), 4.0), a * math.exp(-b * 4.0) - c / 4.0**6)

    def test_a_given_pair_is_shifted_at_the_cutoff_like_the_others(self):
        force_field = raspalib.ForceField(
            pseudo_atoms=pseudo_atoms("A", "B"),
            parameters=[raspalib.VDWParameters(100.0, 3.0), raspalib.VDWParameters(400.0, 4.0)],
            shifted=True,
            use_charge=False,
            binary_interactions=[("A", "B", raspalib.VDWParameters(50.0, 3.2))],
        )
        expected = lennard_jones(50.0, 3.2, 4.0) - lennard_jones(50.0, 3.2, 12.0)
        self.assert_energy(pair_energy(force_field, (0, 1), 4.0), expected)

    def test_every_form_against_its_formula(self):
        """Each form with its parameters as in a force-field file, at distances inside the cutoff of 12 Angstrom."""
        types = raspalib.VDWParameters.Type
        cutoff = 12.0
        epsilon, sigma = 120.0, 3.4

        def u6(r):
            return (sigma / r) ** 6

        def lj(r):
            return 4.0 * epsilon * (u6(r) ** 2 - u6(r))

        def tang_toennies(order, x):
            return 1.0 - math.exp(-x) * sum(x**k / math.factorial(k) for k in range(order + 1))

        # hbar^2 / (24 mu k_B T) in Angstrom^2, for a reduced mass mu in atomic mass units
        hbar, atomic_mass, boltzmann = 1.054571817e-34, 1.66053906660e-27, 1.380649e-23
        reduced_mass, temperature = 2.0, 300.0
        feynman_hibbs = hbar**2 / (24.0 * reduced_mass * atomic_mass * boltzmann * temperature) * 1.0e20

        c6 = (sigma / cutoff) ** 6
        # the switched forms start switching 2 Angstrom below the cutoff unless the force field says otherwise
        switch = cutoff - 2.0

        def switched(r):
            if r <= switch:
                return lj(r)
            x = (r - switch) / (cutoff - switch)
            return lj(r) * (1.0 - 10.0 * x**3 + 15.0 * x**4 - 6.0 * x**5)

        def force_switched(r):
            # P. J. Steinbach and B. R. Brooks, J. Comput. Chem. 1994, 15, 667-683
            a, b = 4.0 * epsilon * sigma**12, 4.0 * epsilon * sigma**6
            if r <= switch:
                return a * (r**-12 - (cutoff * switch) ** -6) - b * (r**-6 - (cutoff * switch) ** -3)
            return (
                a * cutoff**6 / (cutoff**6 - switch**6) * (r**-6 - cutoff**-6) ** 2
                - b * cutoff**3 / (cutoff**3 - switch**3) * (r**-3 - cutoff**-3) ** 2
            )

        forms = {
            types.NONE: ([], lambda r: 0.0),
            types.LENNARD_JONES: ([epsilon, sigma], lj),
            types.BUCKINGHAM: ([1.0e7, 3.5, 2.0e5], lambda r: 1.0e7 * math.exp(-3.5 * r) - 2.0e5 / r**6),
            types.MORSE: (
                [150.0, 1.4, 3.8],
                lambda r: 150.0 * ((1.0 - math.exp(-1.4 * (r - 3.8))) ** 2 - 1.0),
            ),
            types.FEYNMAN_HIBBS: (
                [epsilon, sigma, reduced_mass],
                lambda r: lj(r) + 4.0 * epsilon * feynman_hibbs * (132.0 * u6(r) ** 2 - 30.0 * u6(r)) / r**2,
            ),
            types.MM3: (
                [epsilon, sigma],
                lambda r: epsilon * (1.84e5 * math.exp(-12.0 * r / sigma) - 2.25 * (sigma / r) ** 6),
            ),
            types.BORN_HUGGINS_MEYER: (
                [300.0, 3.2, 3.0, 1.5e5, 4.0e5],
                lambda r: 300.0 * math.exp(3.2 * (3.0 - r)) - 1.5e5 / r**6 - 4.0e5 / r**8,
            ),
            types.LENNARD_JONES_SHIFTED_FORCE: (
                [epsilon, sigma],
                lambda r: 4.0 * epsilon
                * (u6(r) ** 2 - u6(r) - c6 * (c6 - 1.0) + (12.0 * c6**2 - 6.0 * c6) * (r - cutoff) / cutoff),
            ),
            types.LENNARD_JONES_SECOND_ORDER_TAYLOR_SHIFTED: (
                [epsilon, sigma],
                lambda r: 4.0 * epsilon
                * (
                    u6(r) ** 2
                    - u6(r)
                    - c6 * (c6 - 1.0)
                    + (12.0 * c6**2 - 6.0 * c6) * (r - cutoff) / cutoff
                    - 0.5 * (156.0 * c6**2 - 42.0 * c6) * ((r - cutoff) / cutoff) ** 2
                ),
            ),
            types.LENNARD_JONES_SWITCHED: ([epsilon, sigma], switched),
            types.LENNARD_JONES_FORCE_SWITCHED: ([epsilon, sigma], force_switched),
            types.POTENTIAL_12_6: ([9.0e8, 6.0e5], lambda r: 9.0e8 / r**12 - 6.0e5 / r**6),
            types.POTENTIAL_12_6_2_0: (
                [9.0e8, -6.0e5, 40.0, 5.0],
                lambda r: 9.0e8 / r**12 - 6.0e5 / r**6 + 40.0 / r**2 + 5.0,
            ),
            types.CFF_9_6: ([3.0e7, 5.0e5], lambda r: 3.0e7 / r**9 - 5.0e5 / r**6),
            types.CFF_EPSILON_SIGMA: (
                [epsilon, sigma],
                lambda r: epsilon * (2.0 * (sigma / r) ** 9 - 3.0 * (sigma / r) ** 6),
            ),
            types.MATSUOKA_CLEMENTI_YOSHIMINE: (
                [8.0e6, 3.6, -9.0e4, 1.3],
                lambda r: 8.0e6 * math.exp(-3.6 * r) - 9.0e4 * math.exp(-1.3 * r),
            ),
            types.GENERIC: (
                [1.0e7, 3.5, 2.0e3, 2.0e5, 3.0e5, 4.0e5],
                lambda r: 1.0e7 * math.exp(-3.5 * r) - 2.0e3 / r**4 - 2.0e5 / r**6 - 3.0e5 / r**8 - 4.0e5 / r**10,
            ),
            types.PELLENQ_NICHOLSON: (
                [1.0e7, 3.5, 2.0e5, 3.0e5, 4.0e5],
                lambda r: 1.0e7 * math.exp(-3.5 * r)
                - tang_toennies(6, 3.5 * r) * 2.0e5 / r**6
                - tang_toennies(8, 3.5 * r) * 3.0e5 / r**8
                - tang_toennies(10, 3.5 * r) * 4.0e5 / r**10,
            ),
            types.HYDRATED_ION_WATER: (
                [1.0e7, 3.5, 2.0e3, 2.0e5, 4.0e7],
                lambda r: 1.0e7 * math.exp(-3.5 * r) - 2.0e3 / r**4 - 2.0e5 / r**6 - 4.0e7 / r**12,
            ),
            types.MIE: ([2.0e9, 14.0, 5.0e5, 6.0], lambda r: 2.0e9 / r**14.0 - 5.0e5 / r**6.0),
            types.WEEKS_CHANDLER_ANDERSEN: (
                [epsilon, sigma],
                lambda r: lj(r) + epsilon if r < 2.0 ** (1.0 / 6.0) * sigma else 0.0,
            ),
            types.REPULSIVE_HARMONIC: (
                [500.0, 4.0],
                lambda r: 0.5 * 500.0 * (1.0 - r / 4.0) ** 2 if r < 4.0 else 0.0,
            ),
        }
        self.assertEqual(set(forms), set(types))

        for type_, (parameters, formula) in forms.items():
            potential = raspalib.VDWParameters(type=type_, parameters=parameters)
            self.assertEqual(potential.type, type_)
            force_field = self.force_field([potential, raspalib.VDWParameters(400.0, 4.0)])
            for distance in (3.3, 3.7, 4.5, 7.0, 10.5, 11.5):
                with self.subTest(form=type_.name, distance=distance):
                    self.assert_energy(pair_energy(force_field, (0, 0), distance), formula(distance))

    def test_what_is_refused(self):
        types = raspalib.VDWParameters.Type
        with self.assertRaises(ValueError):
            raspalib.VDWParameters(type=types.BUCKINGHAM, parameters=[1.0e7])
        with self.assertRaises(ValueError):
            self.force_field([raspalib.VDWParameters(100.0, 3.0)])
        with self.assertRaises(ValueError):
            self.force_field(
                [raspalib.VDWParameters(100.0, 3.0), raspalib.VDWParameters(400.0, 4.0)],
                [("A", "C", raspalib.VDWParameters(50.0, 3.2))],
            )


def alkane_force_field() -> raspalib.ForceField:
    # D. Dubbeldam et al., J. Phys. Chem. B 2004, 108(33), 12301-12313
    return raspalib.ForceField(
        pseudo_atoms=[
            raspalib.PseudoAtom(name="CH3", framework_type=False, mass=15.04246, charge=0.0, atomic_number=6),
            raspalib.PseudoAtom(name="CH2", framework_type=False, mass=14.04246, charge=0.0, atomic_number=6),
        ],
        parameters=[raspalib.VDWParameters(108.0, 3.76), raspalib.VDWParameters(56.0, 3.96)],
        mixing_rule=raspalib.ForceField.MixingRule.LORENTZ_BERTHELOT,
        shifted=True,
        tail_corrections=False,
        use_charge=False,
    )


def butane_potentials() -> raspalib.IntraMolecularPotentials:
    # M. G. Martin and J. I. Siepmann, J. Phys. Chem. B 1998, 102, 2569-2577
    return raspalib.IntraMolecularPotentials(
        bonds=[
            raspalib.BondPotential(identifiers=(0, 1), type=raspalib.BondType.HARMONIC, parameters=[96500.0, 1.54]),
            raspalib.BondPotential(identifiers=(1, 2), type=raspalib.BondType.HARMONIC, parameters=[96500.0, 1.54]),
            raspalib.BondPotential(identifiers=(2, 3), type=raspalib.BondType.HARMONIC, parameters=[96500.0, 1.54]),
        ],
        bends=[
            raspalib.BendPotential(identifiers=(0, 1, 2), type=raspalib.BendType.HARMONIC, parameters=[62500.0, 114.0]),
            raspalib.BendPotential(identifiers=(1, 2, 3), type=raspalib.BendType.HARMONIC, parameters=[62500.0, 114.0]),
        ],
        torsions=[
            raspalib.TorsionPotential(
                identifiers=(0, 1, 2, 3),
                type=raspalib.TorsionType.TRAPPE,
                parameters=[0.0, 355.03, -68.19, 791.32],
            )
        ],
    )


def butane(
    force_field: raspalib.ForceField,
    potentials: raspalib.IntraMolecularPotentials,
    move_probabilities: raspalib.MCMoveProbabilities | None = None,
    **keywords,
) -> raspalib.Component:
    return raspalib.Component(
        force_field=force_field,
        component_name="butane",
        critical_temperature=425.125,
        critical_pressure=3796000.0,
        acentric_factor=0.201,
        defined_atoms=[
            raspalib.Atom(position=(0.0, 0.0, 0.0), charge=0.0, type=0),
            raspalib.Atom(position=(0.0, 0.0, 0.0), charge=0.0, type=1),
            raspalib.Atom(position=(0.0, 0.0, 0.0), charge=0.0, type=1),
            raspalib.Atom(position=(0.0, 0.0, 0.0), charge=0.0, type=0),
        ],
        intra_molecular_potentials=potentials,
        particle_probabilities=move_probabilities or raspalib.MCMoveProbabilities(),
        **keywords,
    )


def butane_box(move_probabilities: raspalib.MCMoveProbabilities, number_of_molecules: int = 20) -> raspalib.System:
    force_field = alkane_force_field()
    return raspalib.System(
        force_field=force_field,
        simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
        external_temperature=300.0,
        components=[butane(force_field, butane_potentials(), move_probabilities)],
        initial_number_of_molecules=[number_of_molecules],
    )


def geometry(system: raspalib.System, atoms_per_molecule: int = 4) -> tuple[list[float], list[float]]:
    """The bond lengths [Angstrom] and the bend angles [degrees] of the chain molecules of a system."""
    positions = [atom.position for atom in system.atom_data]
    lengths, angles = [], []
    for first in range(0, len(positions), atoms_per_molecule):
        molecule = positions[first : first + atoms_per_molecule]
        for a, b in zip(molecule, molecule[1:]):
            lengths.append(math.dist(a, b))
        for a, b, c in zip(molecule, molecule[1:], molecule[2:]):
            u = [x - y for x, y in zip(a, b)]
            v = [x - y for x, y in zip(c, b)]
            cosine = sum(x * y for x, y in zip(u, v)) / (math.hypot(*u) * math.hypot(*v))
            angles.append(math.degrees(math.acos(cosine)))
    return lengths, angles


class FlexibleMolecules(unittest.TestCase):
    def assert_butane_geometry(self, system: raspalib.System, number_of_molecules: int):
        lengths, angles = geometry(system)
        self.assertEqual(len(lengths), 3 * number_of_molecules)
        self.assertEqual(len(angles), 2 * number_of_molecules)
        # at 300 K the bond length spreads by about 0.06 Angstrom around 1.54, the angle by about 4 degrees around 114
        self.assertAlmostEqual(sum(lengths) / len(lengths), 1.54, delta=0.03)
        self.assertTrue(all(1.2 < length < 1.9 for length in lengths))
        self.assertAlmostEqual(sum(angles) / len(angles), 114.0, delta=3.0)
        self.assertTrue(all(95.0 < angle < 135.0 for angle in angles))

    def test_the_component_takes_its_connectivity_from_the_bonds(self):
        component = butane(alkane_force_field(), butane_potentials())
        self.assertFalse(component.rigid)
        self.assertEqual(component.connectivity_table.number_of_beads, 4)
        self.assertEqual(component.connectivity_table.bonds, [[0, 1], [1, 2], [2, 3]])
        self.assertEqual(len(component.connectivity_table.bends), 2)
        self.assertEqual(component.connectivity_table.torsions, [[0, 1, 2, 3]])
        potentials = component.intra_molecular_potentials
        self.assertEqual(len(potentials.bonds), 3)
        self.assertEqual(len(potentials.bends), 2)
        self.assertEqual(len(potentials.torsions), 1)
        self.assertEqual(potentials.bonds[1].identifiers, [1, 2])
        self.assertEqual(potentials.bonds[1].type, raspalib.BondType.HARMONIC)
        self.assertEqual(potentials.torsions[0].type, raspalib.TorsionType.TRAPPE)

    def test_a_molecule_without_bonds_is_rigid(self):
        force_field = alkane_force_field()
        component = raspalib.Component(
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
        self.assertTrue(component.rigid)
        self.assertEqual(component.connectivity_table.bonds, [])

    def test_monte_carlo_grows_the_molecule_with_its_bonds_and_bends(self):
        moves = raspalib.MCMoveProbabilities(
            translation_probability=1.0, rotation_probability=1.0, reinsertion_cbmc_probability=1.0
        )
        system = butane_box(moves)
        mc = raspalib.MonteCarlo(
            number_of_cycles=300, number_of_initialization_cycles=200, systems=[system], random_seed=5
        )
        mc.run()

        self.assert_butane_geometry(system, 20)
        statistics = system.components[0].mc_moves_statistics[raspalib.Move.Types.REINSERTION_CBMC]
        self.assertGreater(statistics.total_accepted, 0)

    def test_molecular_dynamics_of_the_flexible_molecule(self):
        system = butane_box(raspalib.MCMoveProbabilities(translation_probability=1.0, reinsertion_cbmc_probability=1.0))
        system.time_step = 0.0005
        system.set_thermostat(raspalib.Thermostat())
        md = raspalib.MolecularDynamics(
            number_of_cycles=400,
            number_of_initialization_cycles=100,
            number_of_equilibration_cycles=200,
            systems=[system],
            random_seed=5,
        )
        md.run()

        self.assert_butane_geometry(system, 20)
        temperature, _ = system.average_temperature.result()
        self.assertGreater(temperature, 200.0)
        self.assertLess(temperature, 400.0)

    def test_a_connectivity_table_that_is_given(self):
        table = raspalib.ConnectivityTable(number_of_beads=4, bonds=[(0, 1), (1, 2), (2, 3)])
        self.assertEqual(table.number_of_beads, 4)
        self.assertEqual(table.bonds, [[0, 1], [1, 2], [2, 3]])
        self.assertEqual(table.torsions, [[0, 1, 2, 3]])
        component = butane(alkane_force_field(), butane_potentials(), connectivity_table=table)
        self.assertEqual(component.connectivity_table.bonds, table.bonds)

        with self.assertRaises(ValueError):
            butane(
                alkane_force_field(),
                butane_potentials(),
                connectivity_table=raspalib.ConnectivityTable(number_of_beads=4, bonds=[(0, 1), (1, 2)]),
            )
        with self.assertRaises(ValueError):
            raspalib.ConnectivityTable(number_of_beads=2, bonds=[(0, 2)])

    def test_every_bend_and_torsion_of_the_bonds_needs_a_potential(self):
        complete = butane_potentials()
        without_a_bend = raspalib.IntraMolecularPotentials(
            bonds=complete.bonds, bends=complete.bends[:1], torsions=complete.torsions
        )
        with self.assertRaisesRegex(ValueError, "bend 1-2-3"):
            butane(alkane_force_field(), without_a_bend)

        without_the_torsion = raspalib.IntraMolecularPotentials(bonds=complete.bonds, bends=complete.bends)
        with self.assertRaisesRegex(ValueError, "torsion 0-1-2-3"):
            butane(alkane_force_field(), without_the_torsion)

        # the atoms of a term may be given in either direction
        reversed_terms = raspalib.IntraMolecularPotentials(
            bonds=complete.bonds,
            bends=[
                raspalib.BendPotential(identifiers=(2, 1, 0), type=raspalib.BendType.HARMONIC, parameters=[62500.0, 114.0]),
                raspalib.BendPotential(identifiers=(3, 2, 1), type=raspalib.BendType.HARMONIC, parameters=[62500.0, 114.0]),
            ],
            torsions=[
                raspalib.TorsionPotential(
                    identifiers=(3, 2, 1, 0),
                    type=raspalib.TorsionType.TRAPPE,
                    parameters=[0.0, 355.03, -68.19, 791.32],
                )
            ],
        )
        self.assertFalse(butane(alkane_force_field(), reversed_terms).rigid)

    def test_what_is_refused(self):
        with self.assertRaises(ValueError):
            raspalib.BondPotential(identifiers=(0, 0), type=raspalib.BondType.HARMONIC, parameters=[96500.0, 1.54])
        with self.assertRaises(ValueError):
            raspalib.BondPotential(identifiers=(0, 1), type=raspalib.BondType.HARMONIC, parameters=[96500.0])
        with self.assertRaises(ValueError):
            raspalib.TorsionPotential(identifiers=(0, 1, 2, 3), type=raspalib.TorsionType.TRAPPE, parameters=[0.0])
        with self.assertRaises(TypeError):
            raspalib.BendPotential(identifiers=(0, 1), type=raspalib.BendType.HARMONIC, parameters=[62500.0, 114.0])

        # a term that names an atom the component does not have
        too_long = raspalib.IntraMolecularPotentials(
            bonds=[raspalib.BondPotential(identifiers=(0, 4), type=raspalib.BondType.HARMONIC, parameters=[96500.0, 1.54])]
        )
        with self.assertRaises(ValueError):
            butane(alkane_force_field(), too_long)

        # an atom of a type the force field does not have
        with self.assertRaises(ValueError):
            raspalib.Component(
                force_field=alkane_force_field(),
                component_name="unknown",
                critical_temperature=190.564,
                critical_pressure=4599200.0,
                acentric_factor=0.01142,
                defined_atoms=[raspalib.Atom(position=(0.0, 0.0, 0.0), charge=0.0, type=2)],
            )

    def test_every_kind_of_bonded_term_is_constructed(self):
        terms = [
            (raspalib.BondPotential, raspalib.BondType, 2),
            (raspalib.UreyBradleyPotential, raspalib.UreyBradleyType, 2),
            (raspalib.BendPotential, raspalib.BendType, 3),
            (raspalib.InversionBendPotential, raspalib.InversionBendType, 4),
            (raspalib.TorsionPotential, raspalib.TorsionType, 4),
            (raspalib.BondBondPotential, raspalib.BondBondType, 3),
            (raspalib.BondBendPotential, raspalib.BondBendType, 3),
            (raspalib.BondTorsionPotential, raspalib.BondTorsionType, 4),
            (raspalib.BendBendPotential, raspalib.BendBendType, 4),
            (raspalib.BendTorsionPotential, raspalib.BendTorsionType, 4),
        ]
        for potential, types, number_of_atoms in terms:
            for type_ in types:
                # the number of parameters of a type is found from the refusal of one parameter too many
                constructed, refusals = None, []
                for number_of_parameters in range(0, 8):
                    try:
                        constructed = potential(
                            identifiers=tuple(range(number_of_atoms)), type=type_, parameters=[1.0] * number_of_parameters
                        )
                        break
                    except ValueError as error:
                        refusals.append(str(error))
                if constructed is None:
                    # a type that is listed but that the potential does not support says so
                    self.assertTrue(
                        any("not supported" in refusal for refusal in refusals), f"{potential.__name__} {type_.name}"
                    )
                    continue
                self.assertEqual(constructed.type, type_)
                self.assertEqual(constructed.identifiers, list(range(number_of_atoms)))
                self.assertTrue(repr(constructed))

        pair = raspalib.VanDerWaalsPotential(identifiers=(0, 3), parameters=raspalib.VDWParameters(108.0, 3.76), scaling=0.5)
        self.assertEqual((pair.identifiers, pair.scaling), ([0, 3], 0.5))
        charges = raspalib.CoulombPotential(identifiers=(0, 3), charge_a=0.5, charge_b=-0.5, scaling=0.8)
        self.assertEqual((charges.identifiers, charges.scaling), ([0, 3], 0.8))
        potentials = raspalib.IntraMolecularPotentials(van_der_waals=[pair], coulombs=[charges])
        self.assertEqual(len(potentials.van_der_waals), 1)
        self.assertEqual(len(potentials.coulombs), 1)


def alkane(
    force_field: raspalib.ForceField,
    number_of_atoms: int,
    potentials: raspalib.IntraMolecularPotentials | None = None,
    **keywords,
) -> raspalib.Component:
    """A linear alkane with the bonded terms of butane on every bond, bend and torsion of the chain."""
    if potentials is None:
        potentials = chain_potentials(number_of_atoms)
    types = [0] + [1] * (number_of_atoms - 2) + [0]
    return raspalib.Component(
        force_field=force_field,
        component_name="alkane",
        critical_temperature=425.125,
        critical_pressure=3796000.0,
        acentric_factor=0.201,
        defined_atoms=[raspalib.Atom(position=(0.0, 0.0, 0.0), charge=0.0, type=type_) for type_ in types],
        intra_molecular_potentials=potentials,
        **keywords,
    )


def chain_potentials(number_of_atoms: int, **keywords) -> raspalib.IntraMolecularPotentials:
    return raspalib.IntraMolecularPotentials(
        bonds=[
            raspalib.BondPotential(identifiers=(i, i + 1), type=raspalib.BondType.HARMONIC, parameters=[96500.0, 1.54])
            for i in range(number_of_atoms - 1)
        ],
        bends=[
            raspalib.BendPotential(
                identifiers=(i, i + 1, i + 2), type=raspalib.BendType.HARMONIC, parameters=[62500.0, 114.0]
            )
            for i in range(number_of_atoms - 2)
        ],
        torsions=[
            raspalib.TorsionPotential(
                identifiers=(i, i + 1, i + 2, i + 3),
                type=raspalib.TorsionType.TRAPPE,
                parameters=[0.0, 355.03, -68.19, 791.32],
            )
            for i in range(number_of_atoms - 3)
        ],
        **keywords,
    )


def unshifted_alkane_force_field(**keywords) -> raspalib.ForceField:
    return raspalib.ForceField(
        pseudo_atoms=[
            raspalib.PseudoAtom(name="CH3", framework_type=False, mass=15.04246, charge=0.0, atomic_number=6),
            raspalib.PseudoAtom(name="CH2", framework_type=False, mass=14.04246, charge=0.0, atomic_number=6),
        ],
        parameters=[raspalib.VDWParameters(108.0, 3.76), raspalib.VDWParameters(56.0, 3.96)],
        mixing_rule=raspalib.ForceField.MixingRule.LORENTZ_BERTHELOT,
        shifted=False,
        tail_corrections=False,
        use_charge=False,
        **keywords,
    )


def single_molecule(force_field: raspalib.ForceField, component: raspalib.Component) -> raspalib.System:
    return raspalib.System(
        force_field=force_field,
        simulation_box=raspalib.SimulationBox(a=40.0, b=40.0, c=40.0),
        external_temperature=300.0,
        components=[component],
        initial_number_of_molecules=[1],
    )


def subtract(a, b):
    return [x - y for x, y in zip(a, b)]


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def chain_energy(positions, pair_potential, scaling_14: float, pair_potential_14=None) -> float:
    """The energy [K] of one linear alkane: its bonds, bends and torsions, the 1-4 pairs with their scaling, and
    the pairs further apart in full. 'pair_potential(i, j, distance)' is the Van der Waals energy of two atoms."""
    pair_potential_14 = pair_potential_14 or pair_potential
    number_of_atoms = len(positions)
    energy = 0.0
    for i in range(number_of_atoms - 1):
        energy += 0.5 * 96500.0 * (math.dist(positions[i], positions[i + 1]) - 1.54) ** 2
    for i in range(number_of_atoms - 2):
        u = subtract(positions[i], positions[i + 1])
        v = subtract(positions[i + 2], positions[i + 1])
        angle = math.acos(dot(u, v) / (math.hypot(*u) * math.hypot(*v)))
        energy += 0.5 * 62500.0 * (angle - math.radians(114.0)) ** 2
    for i in range(number_of_atoms - 3):
        b1 = subtract(positions[i + 1], positions[i])
        b2 = subtract(positions[i + 2], positions[i + 1])
        b3 = subtract(positions[i + 3], positions[i + 2])
        n1, n2 = cross(b1, b2), cross(b2, b3)
        cosine = dot(n1, n2) / (math.hypot(*n1) * math.hypot(*n2))
        phi = math.acos(max(-1.0, min(1.0, cosine)))
        # the TraPPE torsion, zero for the trans conformation (phi of 180 degrees)
        energy += 355.03 * (1.0 + math.cos(phi)) - 68.19 * (1.0 - math.cos(2.0 * phi)) + 791.32 * (1.0 + math.cos(3.0 * phi))
    for i in range(number_of_atoms):
        for j in range(i + 3, number_of_atoms):
            distance = math.dist(positions[i], positions[j])
            if j == i + 3:
                energy += scaling_14 * pair_potential_14(i, j, distance)
            else:
                energy += pair_potential(i, j, distance)
    return energy


class IntraMolecularEnergies(unittest.TestCase):
    """The energy of one flexible molecule against the formulas of its terms, for the configuration it has."""

    def mixed(self, number_of_atoms: int, methyl: tuple[float, float] = (108.0, 3.76)):
        parameters = {0: methyl, 1: (56.0, 3.96)}

        def pair_potential(i: int, j: int, distance: float) -> float:
            (epsilon_i, sigma_i), (epsilon_j, sigma_j) = (
                parameters[0 if atom in (0, number_of_atoms - 1) else 1] for atom in (i, j)
            )
            return lennard_jones(math.sqrt(epsilon_i * epsilon_j), 0.5 * (sigma_i + sigma_j), distance)

        return pair_potential

    def assert_energy(self, system: raspalib.System, expected: float):
        energy = system.compute_total_energies().potential_energy() / KELVIN_TO_ENERGY
        self.assertTrue(math.isclose(energy, expected, rel_tol=1.0e-5, abs_tol=1.0e-3), f"{energy} != {expected}")

    def positions(self, system: raspalib.System):
        return [atom.position for atom in system.atom_data]

    def test_the_bonded_terms_of_butane(self):
        force_field = unshifted_alkane_force_field()
        system = single_molecule(force_field, alkane(force_field, 4))
        positions = self.positions(system)
        # the 1-4 pair does not interact unless it is given a scaling
        expected = chain_energy(positions, self.mixed(4), scaling_14=0.0)
        self.assertGreater(expected, 0.0)
        self.assert_energy(system, expected)

    def test_the_pairs_further_apart_than_three_bonds_interact_in_full(self):
        force_field = unshifted_alkane_force_field()
        system = single_molecule(force_field, alkane(force_field, 7))
        self.assert_energy(system, chain_energy(self.positions(system), self.mixed(7), scaling_14=0.0))

    def test_the_scaling_of_the_1_4_pairs(self):
        force_field = unshifted_alkane_force_field()
        component = alkane(force_field, 6, intra_14_van_der_waals_scaling=0.5)
        self.assertEqual(component.intra_14_van_der_waals_scaling, 0.5)
        system = single_molecule(force_field, component)
        positions = self.positions(system)
        with_scaling = chain_energy(positions, self.mixed(6), scaling_14=0.5)
        without = chain_energy(positions, self.mixed(6), scaling_14=0.0)
        self.assertNotAlmostEqual(with_scaling, without, places=2)
        self.assert_energy(system, with_scaling)

    def test_separate_parameters_for_the_1_4_pairs(self):
        # the 1-4 pairs of a methyl group have their own strength and size; the methylene groups keep theirs
        force_field = unshifted_alkane_force_field(parameters_14=[raspalib.VDWParameters(30.0, 3.4), None])
        system = single_molecule(force_field, alkane(force_field, 6, intra_14_van_der_waals_scaling=1.0))
        positions = self.positions(system)
        expected = chain_energy(
            positions, self.mixed(6), scaling_14=1.0, pair_potential_14=self.mixed(6, methyl=(30.0, 3.4))
        )
        regular = chain_energy(positions, self.mixed(6), scaling_14=1.0)
        self.assertNotAlmostEqual(expected, regular, places=2)
        self.assert_energy(system, expected)

    def test_a_1_4_pair_that_is_given(self):
        force_field = unshifted_alkane_force_field(
            binary_interactions_14=[("CH3", "CH2", raspalib.VDWParameters(20.0, 3.3))]
        )
        system = single_molecule(force_field, alkane(force_field, 5, intra_14_van_der_waals_scaling=1.0))
        positions = self.positions(system)

        # in pentane the two 1-4 pairs are the methyl-methylene pairs 0-3 and 1-4
        def given_14(i: int, j: int, distance: float) -> float:
            return lennard_jones(20.0, 3.3, distance)

        self.assert_energy(
            system, chain_energy(positions, self.mixed(5), scaling_14=1.0, pair_potential_14=given_14)
        )

    def test_a_cmap_term(self):
        # a map of one value adds that value, whatever the two dihedral angles are
        force_field = unshifted_alkane_force_field()
        flat = raspalib.CMAPMap(name="flat", resolution=24, energies=[150.0] * (24 * 24))
        self.assertEqual((flat.name, flat.resolution), ("flat", 24))
        potentials = chain_potentials(
            5, cmap_maps=[flat], cmaps=[raspalib.CMAPPotential(identifiers=(0, 1, 2, 3, 4), map_index=0)]
        )
        self.assertEqual(potentials.cmaps[0].identifiers, [0, 1, 2, 3, 4])
        system = single_molecule(force_field, alkane(force_field, 5, potentials))
        self.assert_energy(system, chain_energy(self.positions(system), self.mixed(5), scaling_14=0.0) + 150.0)

        with self.assertRaises(ValueError):
            raspalib.CMAPMap(name="short", resolution=24, energies=[0.0] * 10)
        with self.assertRaises(ValueError):
            alkane(
                force_field,
                5,
                chain_potentials(5, cmaps=[raspalib.CMAPPotential(identifiers=(0, 1, 2, 3, 4), map_index=0)]),
            )


def propyl_positions() -> list[tuple[float, float, float]]:
    """Three atoms with bonds of 1.54 Angstrom at an angle of 114 degrees, and a fourth continuing the chain."""
    half = math.radians(114.0) / 2.0
    x, z = 1.54 * math.sin(half), 1.54 * math.cos(half)
    return [(-x, 0.0, z), (0.0, 0.0, 0.0), (x, 0.0, z), (2.0 * x, 0.0, 0.0)]


class RigidBodies(unittest.TestCase):
    def component(self, force_field: raspalib.ForceField, **keywords) -> raspalib.Component:
        """A chain of four atoms of which the first three are one rigid body and the fourth is flexible."""
        potentials = raspalib.IntraMolecularPotentials(
            bonds=[raspalib.BondPotential(identifiers=(2, 3), type=raspalib.BondType.HARMONIC, parameters=[96500.0, 1.54])],
            bends=[raspalib.BendPotential(identifiers=(1, 2, 3), type=raspalib.BendType.HARMONIC, parameters=[62500.0, 114.0])],
            torsions=[
                raspalib.TorsionPotential(
                    identifiers=(0, 1, 2, 3), type=raspalib.TorsionType.TRAPPE, parameters=[0.0, 355.03, -68.19, 791.32]
                )
            ],
        )
        arguments = dict(
            force_field=force_field,
            component_name="semi_flexible",
            critical_temperature=425.125,
            critical_pressure=3796000.0,
            acentric_factor=0.201,
            defined_atoms=[
                raspalib.Atom(position=position, charge=0.0, type=type_)
                for position, type_ in zip(propyl_positions(), (0, 1, 1, 0))
            ],
            connectivity_table=raspalib.ConnectivityTable(number_of_beads=4, bonds=[(0, 1), (1, 2), (2, 3)]),
            intra_molecular_potentials=potentials,
            particle_probabilities=raspalib.MCMoveProbabilities(
                translation_probability=1.0, rotation_probability=1.0, reinsertion_cbmc_probability=1.0
            ),
            rigid_bodies=[[0, 1, 2]],
        )
        arguments.update(keywords)
        return raspalib.Component(**arguments)

    def test_a_rigid_body_keeps_its_geometry_and_the_rest_is_flexible(self):
        force_field = alkane_force_field()
        component = self.component(force_field)
        self.assertFalse(component.rigid)
        self.assertEqual(component.rigid_bodies, [[0, 1, 2]])
        self.assertEqual(component.connectivity_table.bonds, [[0, 1], [1, 2], [2, 3]])

        system = raspalib.System(
            force_field=force_field,
            simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
            external_temperature=300.0,
            components=[component],
            initial_number_of_molecules=[20],
        )
        mc = raspalib.MonteCarlo(
            number_of_cycles=300, number_of_initialization_cycles=200, systems=[system], random_seed=5
        )
        mc.run()

        reference = propyl_positions()
        positions = [atom.position for atom in system.atom_data]
        flexible_lengths = []
        for first in range(0, len(positions), 4):
            molecule = positions[first : first + 4]
            for i, j in ((0, 1), (1, 2), (0, 2)):
                self.assertAlmostEqual(
                    math.dist(molecule[i], molecule[j]), math.dist(reference[i], reference[j]), places=6
                )
            flexible_lengths.append(math.dist(molecule[2], molecule[3]))
        self.assertAlmostEqual(sum(flexible_lengths) / len(flexible_lengths), 1.54, delta=0.04)
        self.assertGreater(max(flexible_lengths) - min(flexible_lengths), 0.02)
        statistics = system.components[0].mc_moves_statistics[raspalib.Move.Types.REINSERTION_CBMC]
        self.assertGreater(statistics.total_accepted, 0)

    def test_all_atoms_in_one_rigid_body_is_a_rigid_molecule(self):
        force_field = alkane_force_field()
        component = self.component(
            force_field, intra_molecular_potentials=raspalib.IntraMolecularPotentials(), rigid_bodies=[[0, 1, 2, 3]]
        )
        self.assertTrue(component.rigid)

    def test_what_is_refused(self):
        force_field = alkane_force_field()
        # the connections inside a rigid body come from the table
        with self.assertRaisesRegex(ValueError, "connectivity_table"):
            self.component(force_field, connectivity_table=raspalib.ConnectivityTable())
        # the bond out of the rigid body needs its potential, and so do the bend and the torsion over it
        with self.assertRaisesRegex(ValueError, "bond 2-3"):
            self.component(force_field, intra_molecular_potentials=raspalib.IntraMolecularPotentials())
        with self.assertRaisesRegex(ValueError, "bend 1-2-3"):
            self.component(
                force_field,
                intra_molecular_potentials=raspalib.IntraMolecularPotentials(
                    bonds=[
                        raspalib.BondPotential(
                            identifiers=(2, 3), type=raspalib.BondType.HARMONIC, parameters=[96500.0, 1.54]
                        )
                    ]
                ),
            )
        # an atom in two rigid bodies, and an atom the component does not have
        with self.assertRaises(ValueError):
            self.component(force_field, rigid_bodies=[[0, 1], [1, 2]])
        with self.assertRaises(ValueError):
            self.component(force_field, rigid_bodies=[[0, 1, 7]])


class ChiralCenters(unittest.TestCase):
    def component(self, positions, chiral_centers) -> raspalib.Component:
        """A center atom with three bonded neighbors."""
        force_field = alkane_force_field()
        potentials = raspalib.IntraMolecularPotentials(
            bonds=[
                raspalib.BondPotential(identifiers=(0, i), type=raspalib.BondType.HARMONIC, parameters=[96500.0, 1.54])
                for i in (1, 2, 3)
            ],
            bends=[
                raspalib.BendPotential(identifiers=(i, 0, j), type=raspalib.BendType.HARMONIC, parameters=[62500.0, 112.0])
                for i, j in ((1, 2), (1, 3), (2, 3))
            ],
        )
        return raspalib.Component(
            force_field=force_field,
            component_name="branched",
            critical_temperature=407.8,
            critical_pressure=3640000.0,
            acentric_factor=0.184,
            defined_atoms=[
                raspalib.Atom(position=position, charge=0.0, type=type_) for position, type_ in zip(positions, (1, 0, 0, 0))
            ],
            intra_molecular_potentials=potentials,
            chiral_centers=chiral_centers,
        )

    def test_the_handedness_follows_from_the_positions(self):
        pyramid = [(0.0, 0.0, 0.0), (1.4, 0.0, -0.6), (-0.7, 1.2, -0.6), (-0.7, -1.2, -0.6)]
        mirrored = [(x, y, -z) for x, y, z in pyramid]

        component = self.component(pyramid, [(0, 1, 2, 3)])
        (center,) = component.chiral_centers
        self.assertEqual(center.identifiers, [0, 1, 2, 3])
        opposite = self.component(mirrored, [(0, 1, 2, 3)]).chiral_centers[0]
        chiralities = raspalib.ChiralCenter.Chirality
        self.assertEqual({center.chirality, opposite.chirality}, {chiralities.R, chiralities.S})
        # naming two of the neighbors in the other order names the other handedness
        self.assertEqual(self.component(pyramid, [(0, 2, 1, 3)]).chiral_centers[0].chirality, opposite.chirality)

        self.assertEqual(self.component(pyramid, []).chiral_centers, [])

    def test_monte_carlo_keeps_the_handedness(self):
        pyramid = [(0.0, 0.0, 0.0), (1.4, 0.0, -0.6), (-0.7, 1.2, -0.6), (-0.7, -1.2, -0.6)]

        def handedness(molecule) -> float:
            return dot(subtract(molecule[1], molecule[0]), cross(subtract(molecule[2], molecule[0]), subtract(molecule[3], molecule[0])))

        for positions in (pyramid, [(x, y, -z) for x, y, z in pyramid]):
            force_field = alkane_force_field()
            component = self.component(positions, [(0, 1, 2, 3)])
            component.mc_moves_probabilities[raspalib.Move.Types.TRANSLATION] = 1.0
            component.mc_moves_probabilities[raspalib.Move.Types.REINSERTION_CBMC] = 1.0
            system = raspalib.System(
                force_field=force_field,
                simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
                external_temperature=300.0,
                components=[component],
                initial_number_of_molecules=[20],
            )
            mc = raspalib.MonteCarlo(
                number_of_cycles=200, number_of_initialization_cycles=100, systems=[system], random_seed=5
            )
            mc.run()
            statistics = system.components[0].mc_moves_statistics[raspalib.Move.Types.REINSERTION_CBMC]
            self.assertGreater(statistics.total_accepted, 0)

            atoms = [atom.position for atom in system.atom_data]
            volumes = [handedness(atoms[first : first + 4]) for first in range(0, len(atoms), 4)]
            self.assertEqual(len(volumes), 20)
            self.assertTrue(all(volume * handedness(positions) > 0.0 for volume in volumes), volumes)

    def test_what_is_refused(self):
        pyramid = [(0.0, 0.0, 0.0), (1.4, 0.0, -0.6), (-0.7, 1.2, -0.6), (-0.7, -1.2, -0.6)]
        planar = [(0.0, 0.0, 0.0), (1.54, 0.0, 0.0), (-0.77, 1.33, 0.0), (-0.77, -1.33, 0.0)]
        with self.assertRaisesRegex(ValueError, "planar"):
            self.component(planar, [(0, 1, 2, 3)])
        # the first atom is the center: the others are bonded to it
        with self.assertRaisesRegex(ValueError, "not bonded"):
            self.component(pyramid, [(1, 0, 2, 3)])
        with self.assertRaises(ValueError):
            self.component(pyramid, [(0, 1, 1, 3)])
        with self.assertRaises(ValueError):
            self.component(pyramid, [(0, 1, 2, 4)])


class ChainMolecules(unittest.TestCase):
    """The descriptions of a chain that some moves need: the atoms a partial reinsertion keeps, and repeat units."""

    def run_moves(self, component: raspalib.Component, moves: dict) -> raspalib.System:
        force_field = alkane_force_field()
        for move, probability in moves.items():
            component.mc_moves_probabilities[move] = probability
        system = raspalib.System(
            force_field=force_field,
            simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
            external_temperature=300.0,
            components=[component],
            initial_number_of_molecules=[15],
        )
        mc = raspalib.MonteCarlo(
            number_of_cycles=200, number_of_initialization_cycles=100, systems=[system], random_seed=5
        )
        mc.run()
        return system

    def test_a_partial_reinsertion_regrows_what_is_not_fixed(self):
        types = raspalib.Move.Types
        moves = {types.TRANSLATION: 1.0, types.PARTIAL_REINSERTION_CBMC: 1.0}
        # either half of hexane, or its interior between the two ends
        fixed_atoms = [[0, 1, 2], [3, 4, 5], [0, 5]]
        component = alkane(alkane_force_field(), 6, partial_reinsertion_fixed_atoms=fixed_atoms)
        self.assertEqual(component.partial_reinsertion_fixed_atoms, fixed_atoms)

        system = self.run_moves(component, moves)
        statistics = system.components[0].mc_moves_statistics[types.PARTIAL_REINSERTION_CBMC]
        self.assertGreater(statistics.total_counts, 0)
        self.assertGreater(statistics.total_accepted, 0)
        lengths, angles = geometry(system, atoms_per_molecule=6)
        self.assertAlmostEqual(sum(lengths) / len(lengths), 1.54, delta=0.03)
        self.assertAlmostEqual(sum(angles) / len(angles), 114.0, delta=3.0)

        # without fixed atoms the move has nothing to do
        system = self.run_moves(alkane(alkane_force_field(), 6), moves)
        statistics = system.components[0].mc_moves_statistics[types.PARTIAL_REINSERTION_CBMC]
        self.assertEqual(statistics.total_counts, 0)

    def homopolymer(self, **keywords) -> raspalib.Component:
        """A chain of six equal beads."""
        return raspalib.Component(
            force_field=alkane_force_field(),
            component_name="chain",
            critical_temperature=425.125,
            critical_pressure=3796000.0,
            acentric_factor=0.201,
            defined_atoms=[raspalib.Atom(position=(0.0, 0.0, 0.0), charge=0.0, type=1) for _ in range(6)],
            intra_molecular_potentials=chain_potentials(6),
            **keywords,
        )

    def test_reptation_of_a_chain_with_repeat_units(self):
        types = raspalib.Move.Types
        component = self.homopolymer(repeat_units=[[0], [1], [2], [3], [4], [5]])
        self.assertEqual(component.repeat_units, [[0], [1], [2], [3], [4], [5]])
        # the ends of the chain follow from its repeat units
        self.assertEqual(component.end_to_end_atoms, [0, 5])

        system = self.run_moves(component, {types.TRANSLATION: 1.0, types.REPTATION: 1.0})
        statistics = system.components[0].mc_moves_statistics[types.REPTATION]
        self.assertGreater(statistics.total_accepted, 0)
        lengths, angles = geometry(system, atoms_per_molecule=6)
        self.assertAlmostEqual(sum(lengths) / len(lengths), 1.54, delta=0.03)
        self.assertAlmostEqual(sum(angles) / len(angles), 114.0, delta=3.0)

    def test_the_end_to_end_atoms(self):
        force_field = alkane_force_field()
        # the ends of a chain are the two atoms furthest apart along the bonds
        chain = alkane(force_field, 6)
        self.assertEqual(chain.end_to_end_atoms, [0, 5])
        self.assertEqual(chain.backbone_atoms(), [0, 1, 2, 3, 4, 5])

        given = alkane(force_field, 6, end_to_end_atoms=(1, 4))
        self.assertEqual(given.end_to_end_atoms, [1, 4])
        self.assertEqual(given.backbone_atoms(), [1, 2, 3, 4])

        rigid = raspalib.Component(
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
        self.assertIsNone(rigid.end_to_end_atoms)

    def test_what_is_refused(self):
        force_field = alkane_force_field()
        # the repeat units of hexane are not alike: its ends are methyl groups
        with self.assertRaisesRegex(ValueError, "shift-periodic"):
            alkane(force_field, 6, repeat_units=[[0], [1], [2], [3], [4], [5]])
        with self.assertRaises(ValueError):
            self.homopolymer(repeat_units=[[0, 1, 2], [3, 4]])
        with self.assertRaisesRegex(ValueError, "nothing to regrow"):
            alkane(force_field, 6, partial_reinsertion_fixed_atoms=[[0, 1, 2, 3, 4, 5]])
        with self.assertRaisesRegex(ValueError, "out of range"):
            alkane(force_field, 6, partial_reinsertion_fixed_atoms=[[9]])
        with self.assertRaises(ValueError):
            alkane(force_field, 6, end_to_end_atoms=(0, 0))
        # these describe a flexible molecule
        with self.assertRaisesRegex(ValueError, "rigid"):
            raspalib.Component(
                force_field=force_field,
                component_name="ethane",
                critical_temperature=305.33,
                critical_pressure=4871800.0,
                acentric_factor=0.0993,
                defined_atoms=[
                    raspalib.Atom(position=(0.0, 0.0, 0.77), charge=0.0, type=0),
                    raspalib.Atom(position=(0.0, 0.0, -0.77), charge=0.0, type=0),
                ],
                end_to_end_atoms=(0, 1),
            )


class ParameterDocumentation(unittest.TestCase):
    def test_every_potential_lists_the_parameters_of_its_types(self):
        description = raspalib.BondPotential.describe_parameters()
        self.assertIn("HARMONIC: p_0/k_B [K/Å^2], p_1 [Å]", description)
        self.assertIn("FIXED", description)
        self.assertEqual(raspalib.BondPotential.__doc__, description)
        for potential, types in (
            (raspalib.BondPotential, raspalib.BondType),
            (raspalib.UreyBradleyPotential, raspalib.UreyBradleyType),
            (raspalib.BendPotential, raspalib.BendType),
            (raspalib.InversionBendPotential, raspalib.InversionBendType),
            (raspalib.TorsionPotential, raspalib.TorsionType),
            (raspalib.BondBondPotential, raspalib.BondBondType),
            (raspalib.BondBendPotential, raspalib.BondBendType),
            (raspalib.BondTorsionPotential, raspalib.BondTorsionType),
            (raspalib.BendBendPotential, raspalib.BendBendType),
            (raspalib.BendTorsionPotential, raspalib.BendTorsionType),
        ):
            description = potential.describe_parameters()
            for type_ in types:
                self.assertIn(f"  {type_.name}: ", description, potential.__name__)
            # no values are left in the description
            self.assertNotIn("=", description, potential.__name__)

    def test_the_terms_that_have_no_energy_are_not_offered(self):
        self.assertFalse(hasattr(raspalib, "OutOfPlaneBendPotential"))
        self.assertFalse(hasattr(raspalib.BendType, "MM3_INPLANE"))


if __name__ == "__main__":
    unittest.main()
