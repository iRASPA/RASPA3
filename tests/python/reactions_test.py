"""Reactions in the reaction ensemble, and the partition functions of the molecules that react.

A reaction says how many molecules of each component it takes and gives, in the order of the components of the
system. It is driven by a reaction move of the system; the equilibrium follows from the partition functions of the
isolated molecules, which are given as numbers or taken from the thermochemical data of a species.

Run with pytest, or with the standard library:

    python -m unittest discover -s tests/python -p "*_test.py"
"""

import math
import os
import unittest

os.environ.setdefault("RASPALIB_QUIET", "1")

import raspalib

from settings_test import methane, methane_force_field

REACTION_MOVES = {
    "reaction_cbmc_probability": raspalib.Move.Types.REACTION_CBMC,
    "reaction_cfcmc_probability": raspalib.Move.Types.REACTION_CFCMC,
    "reaction_cbcfcmc_probability": raspalib.Move.Types.REACTION_CBCFCMC,
    "reaction_conventional_cfcmc_probability": raspalib.Move.Types.REACTION_CONVENTIONAL_CFCMC,
    "reaction_conventional_cbcfcmc_probability": raspalib.Move.Types.REACTION_CONVENTIONAL_CBCFCMC,
}


def species(names: list[str], ln_partition_functions: list[float]) -> list[raspalib.Component]:
    """Components that are alike but for their name and their partition function."""
    force_field = methane_force_field()
    components = []
    for name, ln_partition_function in zip(names, ln_partition_functions):
        component = methane(force_field, name, raspalib.MCMoveProbabilities(translation_probability=1.0))
        component.ln_partition_function = ln_partition_function
        components.append(component)
    return components


def reacting_system(
    components: list[raspalib.Component],
    number_of_molecules: list[int],
    reactions: list[raspalib.Reaction],
    **system_probabilities,
) -> raspalib.System:
    return raspalib.System(
        force_field=methane_force_field(),
        simulation_box=raspalib.SimulationBox(a=40.0, b=40.0, c=40.0),
        external_temperature=300.0,
        components=components,
        initial_number_of_molecules=number_of_molecules,
        system_probabilities=raspalib.MCMoveProbabilities(**system_probabilities),
        reactions=reactions,
    )


def run(system: raspalib.System, cycles: int = 2000, equilibration: int = 0) -> None:
    mc = raspalib.MonteCarlo(
        number_of_cycles=cycles,
        number_of_initialization_cycles=1000,
        number_of_equilibration_cycles=equilibration,
        systems=[system],
        random_seed=7,
    )
    mc.run()


class Isomerization(unittest.TestCase):
    """A <-> B for two species that interact alike: the composition follows from the partition functions alone,
    N_B / N_A = q_B / q_A."""

    def system(self, move_keyword: str, **reaction_keywords) -> raspalib.System:
        return reacting_system(
            species(["A", "B"], [0.0, math.log(3.0)]),
            [100, 0],
            [raspalib.Reaction(reactants=[1, 0], products=[0, 1], **reaction_keywords)],
            **{move_keyword: 1.0},
        )

    def test_the_equilibrium_composition(self):
        system = self.system("reaction_cbmc_probability")
        run(system, cycles=3000)

        loading, _ = system.average_loadings.result()
        number_a, number_b = loading.number_of_molecules
        self.assertAlmostEqual(number_a + number_b, 100.0, places=6)
        self.assertAlmostEqual(number_b / number_a, 3.0, delta=0.15)

        statistics = system.mc_moves_statistics[raspalib.Move.Types.REACTION_CBMC]
        self.assertGreater(statistics.total_accepted, 0)

    def test_every_reaction_move_reaches_it(self):
        for move_keyword, move in REACTION_MOVES.items():
            with self.subTest(move=move.name):
                system = self.system(move_keyword)
                run(system, cycles=2000, equilibration=1000)
                self.assertEqual(system.reactions[0].move, move)
                loading, _ = system.average_loadings.result()
                number_a, number_b = loading.number_of_molecules
                self.assertAlmostEqual(number_b / number_a, 3.0, delta=0.35)

    def test_the_other_way_around(self):
        # starting from the products, with the reaction written from B to A
        system = reacting_system(
            species(["A", "B"], [0.0, math.log(3.0)]),
            [0, 100],
            [raspalib.Reaction(reactants=[0, 1], products=[1, 0])],
            reaction_cbmc_probability=1.0,
        )
        run(system, cycles=3000)
        number_a, number_b = system.average_loadings.result()[0].number_of_molecules
        self.assertAlmostEqual(number_b / number_a, 3.0, delta=0.15)

    def test_the_lambda_histogram_of_a_fractional_reaction(self):
        system = self.system("reaction_cfcmc_probability", maximum_lambda_change=0.2, lambda_switch_point=0.25)
        reaction = system.reactions[0]
        self.assertEqual((reaction.maximum_lambda_change, reaction.lambda_switch_point), (0.2, 0.25))
        run(system, cycles=500, equilibration=500)

        reaction = system.reactions[0]
        self.assertGreaterEqual(reaction.current_lambda, 0.0)
        self.assertLessEqual(reaction.current_lambda, 1.0)
        histogram, _ = reaction.lambda_histogram.result()
        self.assertEqual(len(histogram), 21)
        self.assertGreater(sum(histogram), 0.0)


class Stoichiometry(unittest.TestCase):
    def test_a_dimerization_keeps_the_number_of_monomer_units(self):
        system = reacting_system(
            species(["monomer", "dimer"], [0.0, 2.0]),
            [100, 0],
            [raspalib.Reaction(reactants=[2, 0], products=[0, 1])],
            reaction_cbmc_probability=1.0,
        )
        run(system, cycles=500)
        monomers, dimers = system.number_of_molecules_per_component
        self.assertGreater(dimers, 0)
        self.assertEqual(monomers + 2 * dimers, 100)
        average_monomers, average_dimers = system.average_loadings.result()[0].number_of_molecules
        self.assertAlmostEqual(average_monomers + 2.0 * average_dimers, 100.0, places=6)

    def test_two_reactions_in_one_system(self):
        # 2 A <-> B + D and D <-> C, as the metathesis of propene with the isomerization of its product
        system = reacting_system(
            species(["A", "B", "C", "D"], [0.0, 0.0, 0.5, 0.0]),
            [100, 0, 0, 0],
            [
                raspalib.Reaction(reactants=[2, 0, 0, 0], products=[0, 1, 0, 1]),
                raspalib.Reaction(reactants=[0, 0, 0, 1], products=[0, 0, 1, 0]),
            ],
            reaction_cbmc_probability=1.0,
        )
        self.assertEqual([reaction.id for reaction in system.reactions], [0, 1])
        self.assertEqual(system.reactions[1].reactants, [0, 0, 0, 1])
        self.assertEqual(system.reactions[1].products, [0, 0, 1, 0])
        run(system, cycles=500)

        a, b, c, d = system.number_of_molecules_per_component
        self.assertEqual(a + b + c + d, 100)
        # every B comes with a D, which is a D still or has become a C
        self.assertEqual(b, c + d)
        self.assertGreater(b, 0)
        self.assertGreater(c, 0)


class ReactionsAndMoves(unittest.TestCase):
    def test_a_reaction_takes_the_reaction_move_of_its_system(self):
        reaction = raspalib.Reaction(reactants=[1, 0], products=[0, 1])
        self.assertIsNone(reaction.move)
        self.assertEqual((reaction.reactants, reaction.products), ([1, 0], [0, 1]))
        system = reacting_system(species(["A", "B"], [0.0, 0.0]), [10, 0], [reaction], reaction_cfcmc_probability=1.0)
        self.assertEqual(system.reactions[0].move, raspalib.Move.Types.REACTION_CFCMC)
        # the reaction that was given is a description; the system has its own
        self.assertIsNone(reaction.move)

    def test_with_several_reaction_moves_a_reaction_names_its_own(self):
        moves = dict(reaction_cbmc_probability=1.0, reaction_cfcmc_probability=1.0)
        named = raspalib.Reaction(reactants=[1, 0], products=[0, 1], move=raspalib.Move.Types.REACTION_CFCMC)
        system = reacting_system(species(["A", "B"], [0.0, 0.0]), [10, 0], [named], **moves)
        self.assertEqual(system.reactions[0].move, raspalib.Move.Types.REACTION_CFCMC)

        unnamed = raspalib.Reaction(reactants=[1, 0], products=[0, 1])
        with self.assertRaisesRegex(ValueError, "name the move"):
            reacting_system(species(["A", "B"], [0.0, 0.0]), [10, 0], [unnamed], **moves)

    def test_a_system_without_reactions(self):
        system = reacting_system(species(["A", "B"], [0.0, 0.0]), [10, 0], [])
        self.assertEqual(system.reactions, [])

    def test_the_reactions_are_in_the_system_a_driver_gives_back(self):
        system = reacting_system(
            species(["A", "B"], [0.0, 0.0]),
            [20, 0],
            [raspalib.Reaction(reactants=[1, 0], products=[0, 1])],
            reaction_cbmc_probability=1.0,
        )
        mc = raspalib.MonteCarlo(number_of_cycles=50, number_of_initialization_cycles=50, systems=[system], random_seed=7)
        mc.run()
        self.assertEqual(len(system.reactions), 1)
        self.assertEqual(len(mc.systems[0].reactions), 1)
        self.assertEqual(sum(system.number_of_molecules_per_component), 20)

    def test_what_is_refused(self):
        components = lambda: species(["A", "B"], [0.0, 0.0])
        isomerization = lambda **keywords: raspalib.Reaction(reactants=[1, 0], products=[0, 1], **keywords)

        # reactions without a reaction move that has a probability
        with self.assertRaisesRegex(ValueError, "no reaction move"):
            reacting_system(components(), [10, 0], [isomerization()])
        # a stoichiometry that is not of the components of the system
        with self.assertRaisesRegex(ValueError, "each of the 2 components"):
            reacting_system(
                components(),
                [10, 0],
                [raspalib.Reaction(reactants=[1, 0, 0], products=[0, 1, 0])],
                reaction_cbmc_probability=1.0,
            )
        # a move that has no probability in the system
        with self.assertRaises(ValueError):
            reacting_system(
                components(),
                [10, 0],
                [isomerization(move=raspalib.Move.Types.REACTION_CFCMC)],
                reaction_cbmc_probability=1.0,
            )
        # more lambda coordinates than can be followed: each of these reactions takes two
        with self.assertRaisesRegex(ValueError, "lambda coordinates"):
            reacting_system(
                components(), [10, 0], [isomerization()] * 3, reaction_conventional_cfcmc_probability=1.0
            )

        with self.assertRaises(ValueError):
            isomerization(move=raspalib.Move.Types.TRANSLATION)
        with self.assertRaises(ValueError):
            raspalib.Reaction(reactants=[1, 0], products=[0, 1, 0])
        with self.assertRaises(ValueError):
            raspalib.Reaction(reactants=[0, 0], products=[0, 1])
        with self.assertRaises(ValueError):
            isomerization(maximum_lambda_change=0.0)
        with self.assertRaises(ValueError):
            isomerization(lambda_switch_point=1.5)


class PartitionFunctions(unittest.TestCase):
    def component(self) -> raspalib.Component:
        return methane(methane_force_field(), "CO2", raspalib.MCMoveProbabilities(translation_probability=1.0))

    def system(self, component: raspalib.Component, temperature: float) -> raspalib.System:
        return raspalib.System(
            force_field=methane_force_field(),
            simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
            external_temperature=temperature,
            components=[component],
            initial_number_of_molecules=[5],
        )

    def test_a_partition_function_that_is_given(self):
        component = self.component()
        self.assertEqual(component.ln_partition_function, 0.0)
        component.ln_partition_function = 87.1384
        self.assertEqual(component.ln_partition_function, 87.1384)
        self.assertIsNone(component.ln_partition_function_species)
        self.assertEqual(self.system(component, 450.0).components[0].ln_partition_function, 87.1384)

    def test_the_partition_function_of_a_species(self):
        databases = raspalib.ThermochemicalDatabase
        at_300 = raspalib.ln_partition_function(species="CO2", temperature=300.0)
        at_600 = raspalib.ln_partition_function("CO2", 600.0)
        self.assertTrue(math.isfinite(at_300))
        self.assertNotEqual(at_300, at_600)
        # the two databases agree to a few parts in a hundred thousand
        from_janaf = raspalib.ln_partition_function("CO2", 300.0, database=databases.JANAF)
        self.assertNotEqual(from_janaf, at_300)
        self.assertAlmostEqual(from_janaf, at_300, delta=1.0e-4 * abs(at_300))

        self.assertTrue(raspalib.has_thermochemical_data("NH3"))
        self.assertFalse(raspalib.has_thermochemical_data("unobtainium"))
        with self.assertRaises(ValueError):
            raspalib.ln_partition_function("unobtainium", 300.0)
        with self.assertRaises(ValueError):
            raspalib.ln_partition_function("CO2", 0.0)

    def test_a_component_of_a_species_takes_the_temperature_of_its_system(self):
        component = self.component()
        component.set_ln_partition_function_species("CO2")
        self.assertEqual(component.ln_partition_function_species, "CO2")
        for temperature in (300.0, 600.0):
            in_system = self.system(component, temperature).components[0]
            self.assertEqual(in_system.ln_partition_function_species, "CO2")
            self.assertEqual(in_system.ln_partition_function, raspalib.ln_partition_function("CO2", temperature))

        janaf = self.component()
        janaf.set_ln_partition_function_species("CO2", database=raspalib.ThermochemicalDatabase.JANAF)
        self.assertEqual(
            self.system(janaf, 300.0).components[0].ln_partition_function,
            raspalib.ln_partition_function("CO2", 300.0, raspalib.ThermochemicalDatabase.JANAF),
        )

        # a number that is given afterwards is the partition function at every temperature
        component.ln_partition_function = 12.5
        self.assertIsNone(component.ln_partition_function_species)
        self.assertEqual(self.system(component, 600.0).components[0].ln_partition_function, 12.5)

        with self.assertRaises(ValueError):
            component.set_ln_partition_function_species("unobtainium")


if __name__ == "__main__":
    unittest.main()
