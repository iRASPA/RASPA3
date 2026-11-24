"""The settings of Monte Carlo and molecular dynamics that are made on the objects of a simulation.

The probabilities of all moves, the mole fraction and the type of a component, the charge method of a force field,
the ensemble and the barostat of molecular dynamics, the averages of the temperature and the box, restarting from a
binary restart file, and the number of threads.

A setting of a force field or a component is made before the object is handed on: a component and a system keep
their own copy of what they are constructed with.

Run with pytest, or with the standard library:

    python -m unittest discover -s tests/python -p "*_test.py"
"""

import os
import tempfile
import unittest

os.environ.setdefault("RASPALIB_QUIET", "1")

import raspalib

from driver_systems_test import average_energy, co2_force_field, co2_system, translation_counts


def methane_force_field() -> raspalib.ForceField:
    # CH4: M. G. Martin et al., J. Chem. Phys. 2001, 114, 7174-7181
    return raspalib.ForceField(
        pseudo_atoms=[
            raspalib.PseudoAtom(name="CH4", framework_type=False, mass=16.04246, charge=0.0, atomic_number=6)
        ],
        parameters=[raspalib.VDWParameters(158.5, 3.72)],
        mixing_rule=raspalib.ForceField.MixingRule.LORENTZ_BERTHELOT,
        cutoff_framework_vdw=12.0,
        cutoff_molecule_vdw=12.0,
        cutoff_coulomb=12.0,
        shifted=True,
        tail_corrections=False,
        use_charge=False,
    )


def methane(
    force_field: raspalib.ForceField, name: str, move_probabilities: raspalib.MCMoveProbabilities
) -> raspalib.Component:
    return raspalib.Component(
        force_field=force_field,
        component_name=name,
        critical_temperature=190.564,
        critical_pressure=4599200.0,
        acentric_factor=0.01142,
        defined_atoms=[raspalib.Atom(position=(0.0, 0.0, 0.0), charge=0.0, type=0)],
        number_of_blocks=5,
        number_of_lambda_bins=21,
        particle_probabilities=move_probabilities,
        fugacity_coefficient=1.0,
        thermodynamic_integration=False,
    )


def methane_box(
    components: list[raspalib.Component], number_of_molecules: list[int], pressure: float | None = None
) -> raspalib.System:
    return raspalib.System(
        force_field=methane_force_field(),
        simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
        has_external_field=False,
        external_temperature=300.0,
        external_pressure=pressure,
        framework_components=None,
        components=components,
        initial_number_of_molecules=number_of_molecules,
        number_of_blocks=5,
        system_probabilities=raspalib.MCMoveProbabilities(),
    )


class MoveProbabilities(unittest.TestCase):
    def test_every_move_has_a_keyword(self):
        for move in raspalib.Move.Types:
            keyword = move.name.lower() + "_probability"
            probabilities = raspalib.MCMoveProbabilities(**{keyword: 0.25})
            self.assertEqual(probabilities[move], 0.25, keyword)
            others = [other for other in raspalib.Move.Types if other != move]
            self.assertTrue(all(probabilities[other] == 0.0 for other in others), keyword)

    def test_the_keywords_of_the_notebooks_reach_their_moves(self):
        types = raspalib.Move.Types
        probabilities = raspalib.MCMoveProbabilities(
            translation_probability=0.1,
            reinsertion_cbmc_probability=0.2,
            identity_change_probability=0.3,
            gibbs_volume_change_probability=0.4,
            gibbs_swap_cfcmc_probability=0.5,
            widom_probability=0.6,
            widom_cfcmc_probability=0.7,
            widom_cbcfcmc_probability=0.8,
            parallel_tempering_probability=0.9,
            hybrid_mc_probability=1.0,
        )
        self.assertEqual(probabilities[types.TRANSLATION], 0.1)
        self.assertEqual(probabilities[types.REINSERTION_CBMC], 0.2)
        self.assertEqual(probabilities[types.IDENTITY_CHANGE_CBMC], 0.3)
        self.assertEqual(probabilities[types.GIBBS_VOLUME], 0.4)
        self.assertEqual(probabilities[types.GIBBS_SWAP_CFCMC], 0.5)
        self.assertEqual(probabilities[types.WIDOM], 0.6)
        self.assertEqual(probabilities[types.WIDOM_CFCMC], 0.7)
        self.assertEqual(probabilities[types.WIDOM_CBCFCMC], 0.8)
        self.assertEqual(probabilities[types.PARALLEL_TEMPERING], 0.9)
        self.assertEqual(probabilities[types.HYBRID_MC], 1.0)
        self.assertEqual(probabilities[types.GIBBS_IDENTITY_CHANGE_CBMC], 0.0)

    def test_a_probability_is_set_and_read_by_its_move(self):
        probabilities = raspalib.MCMoveProbabilities()
        probabilities[raspalib.Move.Types.PIVOT] = 0.5
        probabilities.set_probability(raspalib.Move.Types.CRANKSHAFT, 0.25)
        self.assertEqual(probabilities.get_probability(raspalib.Move.Types.PIVOT), 0.5)
        self.assertEqual(probabilities[raspalib.Move.Types.CRANKSHAFT], 0.25)

    def test_what_is_refused(self):
        with self.assertRaises(TypeError):
            raspalib.MCMoveProbabilities(no_such_move_probability=1.0)
        with self.assertRaises(TypeError):
            raspalib.MCMoveProbabilities(widom_probability="often")
        with self.assertRaises(TypeError):
            raspalib.MCMoveProbabilities(0.5)
        with self.assertRaises(ValueError):
            raspalib.MCMoveProbabilities(widom_probability=-1.0)
        with self.assertRaises(ValueError):
            raspalib.MCMoveProbabilities()[raspalib.Move.Types.WIDOM] = -1.0

    def test_the_widom_move_is_performed(self):
        moves = raspalib.MCMoveProbabilities(translation_probability=1.0, widom_probability=1.0)
        force_field = methane_force_field()
        system = methane_box([methane(force_field, "CH4", moves)], [40])
        mc = raspalib.MonteCarlo(
            number_of_cycles=100, number_of_initialization_cycles=50, systems=[system], random_seed=12
        )
        mc.run()

        statistics = system.components[0].mc_moves_statistics[raspalib.Move.Types.WIDOM]
        self.assertGreater(statistics.total_counts, 0)
        weight, _ = system.components[0].average_rosenbluth_weights.result()
        self.assertGreater(weight, 0.0)

    def test_a_multiple_try_swap_with_its_number_of_trial_molecules(self):
        moves = raspalib.MCMoveProbabilities(translation_probability=1.0, swap_multi_trial_cbmc_probability=1.0)
        force_field = methane_force_field()
        component = methane(force_field, "CH4", moves)
        component.number_of_trial_molecules = 4
        system = methane_box([component], [20], pressure=1.0e6)
        self.assertEqual(system.components[0].number_of_trial_molecules, 4)
        self.assertTrue(system.components[0].swappable)

        mc = raspalib.MonteCarlo(
            number_of_cycles=100, number_of_initialization_cycles=50, systems=[system], random_seed=12
        )
        mc.run()
        statistics = system.components[0].mc_moves_statistics[raspalib.Move.Types.SWAP_MULTI_TRIAL_CBMC]
        self.assertGreater(sum(statistics.total_counts), 0)


class ComponentSettings(unittest.TestCase):
    def test_the_mole_fractions_set_the_composition_of_a_mixture(self):
        # two components that differ in name only, in an empty box at a low pressure: the numbers of molecules
        # follow the mole fractions of the reservoir
        moves = raspalib.MCMoveProbabilities(translation_probability=0.5, swap_cbmc_probability=1.0)
        force_field = methane_force_field()
        first = methane(force_field, "CH4_a", moves)
        second = methane(force_field, "CH4_b", moves)
        first.mol_fraction = 0.25
        second.mol_fraction = 0.75
        system = methane_box([first, second], [0, 0], pressure=2.0e6)
        self.assertEqual([component.mol_fraction for component in system.components], [0.25, 0.75])

        mc = raspalib.MonteCarlo(
            number_of_cycles=4000, number_of_initialization_cycles=1000, systems=[system], random_seed=12
        )
        mc.run()

        loading, _ = system.average_loadings.result()
        number_first, number_second = loading.number_of_molecules
        self.assertGreater(number_first, 1.0)
        self.assertAlmostEqual(number_second / number_first, 3.0, delta=0.3)

    def test_settings_are_carried_into_the_system(self):
        force_field = co2_force_field()
        system = co2_system()
        component = system.components[0]

        component.type = raspalib.Component.Type.CATION
        component.pair_component = 0
        component.maximum_pair_distance = 6.0
        component.group_components = [0, 0]
        component.maximum_group_distance = 5.0
        component.identity_changes = [1]
        component.identity_switches = [1]
        component.gibbs_identity_changes = [1]
        component.tethered_proton_hop_sites = [[(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]]
        component.ideal_gas_rosenbluth_weight = 1.0
        component.ideal_gas_rosenbluth_weight_tolerance = 0.01
        component.ideal_gas_energy = 250.0
        component.starting_bead = 1
        component.pivot_randomization_fraction = 0.4
        component.pivot_cbmc_randomization_fraction = 0.6
        component.pivot_cbmc_number_of_trial_angles = 6
        component.crankshaft_randomization_fraction = 0.3
        component.crankshaft_max_segment_size = 3
        component.concerted_rotation_randomization_fraction = 0.1
        component.bead_flip_randomization_fraction = 0.7

        copy = raspalib.System(
            force_field=force_field,
            simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
            external_temperature=300.0,
            components=[component],
            initial_number_of_molecules=[5],
        ).components[0]

        self.assertEqual(copy.type, raspalib.Component.Type.CATION)
        self.assertEqual(copy.pair_component, 0)
        self.assertEqual(copy.maximum_pair_distance, 6.0)
        self.assertEqual(copy.group_components, [0, 0])
        self.assertEqual(copy.maximum_group_distance, 5.0)
        self.assertEqual(copy.identity_changes, [1])
        self.assertEqual(copy.identity_switches, [1])
        self.assertEqual(copy.gibbs_identity_changes, [1])
        self.assertEqual(copy.tethered_proton_hop_sites, [[(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]])
        self.assertEqual(copy.ideal_gas_rosenbluth_weight, 1.0)
        self.assertEqual(copy.ideal_gas_rosenbluth_weight_tolerance, 0.01)
        self.assertAlmostEqual(copy.ideal_gas_energy, 250.0)
        self.assertEqual(copy.starting_bead, 1)
        self.assertEqual(copy.pivot_randomization_fraction, 0.4)
        self.assertEqual(copy.pivot_cbmc_randomization_fraction, 0.6)
        self.assertEqual(copy.pivot_cbmc_number_of_trial_angles, 6)
        self.assertEqual(copy.crankshaft_randomization_fraction, 0.3)
        self.assertEqual(copy.crankshaft_max_segment_size, 3)
        self.assertEqual(copy.concerted_rotation_randomization_fraction, 0.1)
        self.assertEqual(copy.bead_flip_randomization_fraction, 0.7)
        self.assertEqual(copy.critical_temperature, 304.1282)
        self.assertEqual(copy.critical_pressure, 7377300.0)
        self.assertEqual(copy.acentric_factor, 0.22394)

    def test_none_clears_what_is_optional(self):
        component = co2_system().components[0]
        for name in (
            "pair_component",
            "maximum_pair_distance",
            "maximum_group_distance",
            "ideal_gas_rosenbluth_weight",
            "ideal_gas_rosenbluth_weight_tolerance",
            "ideal_gas_energy",
            "fugacity_coefficient",
        ):
            setattr(component, name, 1)
            self.assertIsNotNone(getattr(component, name), name)
            setattr(component, name, None)
            self.assertIsNone(getattr(component, name), name)

    def test_what_is_refused(self):
        component = co2_system().components[0]
        for name, value in (
            ("mol_fraction", -0.1),
            ("fugacity_coefficient", 0.0),
            ("starting_bead", 3),
            ("number_of_trial_molecules", 0),
            ("pivot_cbmc_number_of_trial_angles", 0),
            ("crankshaft_max_segment_size", 0),
            ("pivot_randomization_fraction", 1.5),
            ("bead_flip_randomization_fraction", -0.5),
            ("ideal_gas_rosenbluth_weight_tolerance", 0.0),
        ):
            with self.assertRaises(ValueError, msg=name):
                setattr(component, name, value)


class ForceFieldSettings(unittest.TestCase):
    def test_the_charge_method_changes_the_electrostatics(self):
        energies = {}
        for method in (raspalib.ForceField.ChargeMethod.EWALD, raspalib.ForceField.ChargeMethod.WOLF):
            force_field = co2_force_field()
            force_field.charge_method = method
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
                particle_probabilities=raspalib.MCMoveProbabilities(
                    translation_probability=0.5, rotation_probability=0.5
                ),
            )
            system = raspalib.System(
                force_field=force_field,
                simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
                external_temperature=300.0,
                components=[component],
                initial_number_of_molecules=[40],
            )
            mc = raspalib.MonteCarlo(
                number_of_cycles=40, number_of_initialization_cycles=20, systems=[system], random_seed=12
            )
            mc.run()
            energies[method] = average_energy(system)
        ewald, wolf = energies.values()
        self.assertNotEqual(ewald, 0.0)
        self.assertNotEqual(wolf, 0.0)
        self.assertNotEqual(ewald, wolf)

    def test_the_ewald_settings(self):
        force_field = co2_force_field()
        self.assertEqual(force_field.charge_method, raspalib.ForceField.ChargeMethod.EWALD)
        self.assertTrue(force_field.automatic_ewald)

        force_field.set_ewald_parameters(alpha=0.3, number_of_wave_vectors=(6, 7, 8))
        self.assertEqual(force_field.ewald_alpha, 0.3)
        self.assertEqual(force_field.number_of_wave_vectors, (6, 7, 8))
        self.assertFalse(force_field.automatic_ewald)

        force_field.ewald_precision = 1.0e-8
        self.assertEqual(force_field.ewald_precision, 1.0e-8)
        self.assertTrue(force_field.automatic_ewald)

        force_field.modified_shifted_force_beta = 0.25
        force_field.omit_ewald_fourier = True
        force_field.compute_polarization = True
        force_field.omit_inter_polarization = False
        force_field.omit_inter_interactions = True
        self.assertEqual(force_field.modified_shifted_force_beta, 0.25)
        self.assertTrue(force_field.omit_ewald_fourier)
        self.assertTrue(force_field.compute_polarization)
        self.assertFalse(force_field.omit_inter_polarization)
        self.assertTrue(force_field.omit_inter_interactions)
        self.assertEqual(
            (force_field.cutoff_framework_vdw, force_field.cutoff_molecule_vdw, force_field.cutoff_coulomb),
            (12.0, 12.0, 12.0),
        )

    def test_what_is_refused(self):
        force_field = co2_force_field()
        with self.assertRaises(ValueError):
            force_field.ewald_precision = 0.0
        with self.assertRaises(ValueError):
            force_field.modified_shifted_force_beta = -1.0
        with self.assertRaises(ValueError):
            force_field.set_ewald_parameters(alpha=0.0, number_of_wave_vectors=(6, 6, 6))


class MolecularDynamicsEnsembles(unittest.TestCase):
    def run_md(self, system: raspalib.System) -> raspalib.MolecularDynamics:
        md = raspalib.MolecularDynamics(
            number_of_cycles=300,
            number_of_initialization_cycles=100,
            number_of_equilibration_cycles=100,
            systems=[system],
            random_seed=12,
        )
        md.run()
        return md

    def thermostat(self) -> raspalib.Thermostat:
        return raspalib.Thermostat(
            thermostat_chain_length=5, number_of_yoshida_suzuki_steps=5, time_scale_parameter=0.15
        )

    def co2_fluid(self, pressure: float | None) -> raspalib.System:
        force_field = co2_force_field()
        component = co2_system().components[0]
        system = raspalib.System(
            force_field=force_field,
            simulation_box=raspalib.SimulationBox(a=30.0, b=30.0, c=30.0),
            external_temperature=300.0,
            external_pressure=pressure,
            components=[component],
            initial_number_of_molecules=[60],
        )
        system.time_step = 0.0005
        return system

    def test_the_volume_is_constant_with_a_thermostat_only(self):
        system = self.co2_fluid(pressure=None)
        system.set_thermostat(self.thermostat())
        self.assertEqual(system.molecular_dynamics_ensemble, raspalib.MolecularDynamicsEnsemble.NVE)
        self.assertIsNone(system.thermobarostat)
        self.run_md(system)

        self.assertAlmostEqual(system.simulation_box.volume, 27000.0, places=6)
        temperature, _ = system.average_temperature.result()
        self.assertGreater(temperature, 150.0)
        self.assertLess(temperature, 450.0)
        translational, _ = system.average_translational_temperature.result()
        rotational, _ = system.average_rotational_temperature.result()
        self.assertGreater(translational, 0.0)
        self.assertGreater(rotational, 0.0)

    def test_a_barostat_changes_the_volume(self):
        system = self.co2_fluid(pressure=5.0e7)
        system.set_thermostat(self.thermostat())
        system.set_thermobarostat(raspalib.Thermobarostat(ensemble=raspalib.MolecularDynamicsEnsemble.NPT))
        self.assertEqual(system.molecular_dynamics_ensemble, raspalib.MolecularDynamicsEnsemble.NPT)
        self.assertEqual(system.thermobarostat.cell_type, raspalib.CellType.ISOTROPIC)
        self.run_md(system)

        self.assertNotAlmostEqual(system.simulation_box.volume, 27000.0, places=3)
        average_box, _ = system.average_simulation_box.result()
        self.assertNotAlmostEqual(average_box.volume, 27000.0, places=3)

    def test_the_settings_of_a_barostat(self):
        ensembles = raspalib.MolecularDynamicsEnsemble
        barostat = raspalib.Thermobarostat(
            ensemble=ensembles.NPTPR,
            barostat_chain_length=3,
            number_of_yoshida_suzuki_steps=3,
            time_scale_parameter=2.0,
            number_of_respa_steps=2,
            coupling=raspalib.BarostatCoupling.ATOMIC,
        )
        self.assertEqual(barostat.ensemble, ensembles.NPTPR)
        self.assertEqual(barostat.cell_type, raspalib.CellType.REGULAR)
        self.assertEqual(barostat.coupling, raspalib.BarostatCoupling.ATOMIC)
        self.assertEqual(barostat.barostat_chain_length, 3)
        self.assertEqual(barostat.number_of_yoshida_suzuki_steps, 3)
        self.assertEqual(barostat.number_of_respa_steps, 2)
        self.assertEqual(barostat.time_scale_parameter, 2.0)
        self.assertEqual(
            raspalib.Thermobarostat(ensemble=ensembles.NPTPR, cell_type=raspalib.CellType.MONOCLINIC).cell_type,
            raspalib.CellType.MONOCLINIC,
        )

        for ensemble in (ensembles.NVE, ensembles.NVT, ensembles.MUVT):
            with self.assertRaises(ValueError):
                raspalib.Thermobarostat(ensemble=ensemble)
        with self.assertRaises(ValueError):
            raspalib.Thermobarostat(ensemble=ensembles.NPT, barostat_chain_length=0)
        with self.assertRaises(ValueError):
            raspalib.Thermobarostat(ensemble=ensembles.NPT, time_scale_parameter=0.0)

        system = self.co2_fluid(pressure=1.0e5)
        system.set_thermobarostat(barostat)
        self.assertEqual(system.thermobarostat.barostat_chain_length, 3)
        self.assertEqual(system.thermobarostat.coupling, raspalib.BarostatCoupling.ATOMIC)
        system.set_thermobarostat(None)
        self.assertIsNone(system.thermobarostat)


class SystemSettings(unittest.TestCase):
    def test_the_settings_of_a_system(self):
        system = co2_system()
        self.assertEqual(system.temperature, 300.0)

        system.external_pressure_diagonal = (1.0e5, 2.0e5, 3.0e5)
        self.assertEqual(system.external_pressure_diagonal, (1.0e5, 2.0e5, 3.0e5))

        system.molecular_dynamics_ensemble = raspalib.MolecularDynamicsEnsemble.MUVT
        self.assertEqual(system.molecular_dynamics_ensemble, raspalib.MolecularDynamicsEnsemble.MUVT)

        system.number_of_hybrid_mc_steps = 4
        self.assertEqual(system.number_of_hybrid_mc_steps, 4)
        with self.assertRaises(ValueError):
            system.number_of_hybrid_mc_steps = 0


class Restart(unittest.TestCase):
    def driver(self, number_of_cycles: int) -> raspalib.MonteCarlo:
        return raspalib.MonteCarlo(
            number_of_cycles=number_of_cycles,
            number_of_initialization_cycles=20,
            write_binary_restart_every=10,
            systems=[co2_system()],
            random_seed=12,
            output_to_files=True,
        )

    def test_a_run_continues_from_a_binary_restart_file(self):
        previous_directory = os.getcwd()
        with tempfile.TemporaryDirectory() as directory:
            os.chdir(directory)
            try:
                # restart files are written every 10 cycles: the last one is 5 cycles before the end of the run
                first = self.driver(number_of_cycles=45)
                first.run()
                self.assertTrue(os.path.exists("restart_data.bin"))

                second = self.driver(number_of_cycles=45)
                second.restart_from_binary_file("restart_data.bin")
                restored = second.systems[0].components[0]
                self.assertGreater(translation_counts(restored), 0)
                self.assertLess(translation_counts(restored), translation_counts(first.systems[0].components[0]))

                # the remaining cycles are the same ones: the random-number generator is part of the state
                second.run()
                self.assertNotEqual(average_energy(second.systems[0]), 0.0)
                self.assertEqual(average_energy(second.systems[0]), average_energy(first.systems[0]))
                self.assertEqual(
                    translation_counts(second.systems[0].components[0]),
                    translation_counts(first.systems[0].components[0]),
                )

                with self.assertRaises(Exception):
                    self.driver(number_of_cycles=40).restart_from_binary_file("no_such_file.bin")
            finally:
                os.chdir(previous_directory)


class Threads(unittest.TestCase):
    def test_the_number_of_threads(self):
        self.assertGreaterEqual(raspalib.number_of_threads(), 1)
        with self.assertRaises(ValueError):
            raspalib.set_number_of_threads(0)
        raspalib.set_number_of_threads(1)
        self.assertEqual(raspalib.number_of_threads(), 1)


if __name__ == "__main__":
    unittest.main()
