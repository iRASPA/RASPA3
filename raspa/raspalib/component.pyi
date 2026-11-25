from raspalib.reaction import *
from raspalib.molecule_analyses import *
import enum
import collections.abc
from raspalib.connectivity_table import *
from raspalib.intra_molecular_potentials import *
from raspalib.atom import *
from raspalib.force_field import *
from raspalib.mc_move_probabilities import *
from raspalib.mc_move_statistics import *
from raspalib.property_widom import *
from raspalib.property_lambda_probability_histogram import *

class Component():
    """
    Represent a simulation component in RASPA.

    A component groups force-field settings, molecular topology definitions,
    Monte Carlo move settings, and runtime statistics for a specific adsorbate
    or fluid species.
    """

    class Type(enum.IntEnum):
        """Whether the component is an adsorbate or a cation."""
        ADSORBATE: int = 0
        CATION: int = 1

    @staticmethod
    def from_file(
        force_field: ForceField,
        file_name: str,
        component_name: str | None = None,
        particle_probabilities: MCMoveProbabilities = ...,
        type: 'Component.Type' = ...,
        fugacity_coefficient: float | None = None,
        thermodynamic_integration: bool = False,
        number_of_blocks: int = 5,
        number_of_lambda_bins: int = 41,
    ) -> 'Component':
        """
        Read a component from a molecule file (JSON), as the command-line
        program reads the molecule of a component of ``simulation.json``.

        Args:
            force_field: The force field; it has to hold the pseudo-atoms of
                the molecule.
            file_name: The molecule file. The extension ``.json`` may be left
                out, and a file that is not found is looked for in the
                directory of ``RASPA_DIR``.
            component_name: The name of the component. Defaults to the name
                of the file without its extension.
            particle_probabilities: The probabilities of the Monte Carlo
                moves of the component.
            type: Adsorbate (the default) or cation.
            fugacity_coefficient: Optional fugacity coefficient; ``None``
                takes it from the equation of state.
            thermodynamic_integration: Whether dU/dlambda is sampled.
            number_of_blocks: Not needed: a system sizes its components for
                its number of blocks, and a driver its systems.
            number_of_lambda_bins: Number of bins in lambda histograms.

        What the component block of ``simulation.json`` sets besides (the mol
        fraction, the ideal-gas Rosenbluth weight, the analyses of the
        molecule) is set on the object before it is given to a system.

        Raises:
            ValueError: the file is not found or is not a molecule of this
                force field.
        """
        ...

    def __init__(
        self,
        force_field: ForceField,
        component_name: str,
        critical_temperature: float,
        critical_pressure: float,
        acentric_factor: float,
        defined_atoms: collections.abc.Sequence[Atom] = [],
        connectivity_table: ConnectivityTable = ConnectivityTable(),
        intra_molecular_potentials: IntraMolecularPotentials = IntraMolecularPotentials(),
        number_of_blocks: int = 5,
        number_of_lambda_bins: int = 41,
        particle_probabilities: MCMoveProbabilities = MCMoveProbabilities(),
        fugacity_coefficient: float | None = None,
        thermodynamic_integration: bool = False,
        blocking_pockets: collections.abc.Sequence[collections.abc.Sequence[float]] = [],
        rigid_bodies: collections.abc.Sequence[collections.abc.Sequence[int]] = [],
        chiral_centers: collections.abc.Sequence[collections.abc.Sequence[int]] = [],
        intra_14_van_der_waals_scaling: float = 0.0,
        intra_14_charge_charge_scaling: float = 0.0,
        partial_reinsertion_fixed_atoms: collections.abc.Sequence[collections.abc.Sequence[int]] = [],
        repeat_units: collections.abc.Sequence[collections.abc.Sequence[int]] = [],
        end_to_end_atoms: collections.abc.Sequence[int] | None = None,
    ) -> None:
        ...
        """
        Initialize a :class:`Component`.

        Args:
            force_field: Force-field definition used for this component.
            component_name: Human-readable component name.
            critical_temperature: Critical temperature of the component.
            critical_pressure: Critical pressure of the component.
            acentric_factor: Acentric factor of the component.
            defined_atoms: Atom definitions used to build the component
                topology. Defaults to an empty sequence.
            connectivity_table: Bond/angle/dihedral connectivity information.
                Defaults to an empty :class:`ConnectivityTable`.
            intra_molecular_potentials: Intra-molecular potential settings.
                Defaults to an empty :class:`IntraMolecularPotentials`.
            number_of_blocks: Number of sampling blocks used for statistics.
                Defaults to ``5``. Not needed: a system sizes its components
                for its number of blocks, and a driver its systems.
            number_of_lambda_bins: Number of bins in lambda histograms.
                Defaults to ``41``.
            particle_probabilities: Monte Carlo move probabilities for this
                component. Defaults to a fresh
                :class:`MCMoveProbabilities` instance.
        """

    @property
    def blocking_pockets(self) -> collections.abc.Sequence[tuple[float, float, float, float]]:
        """Return blocking-pocket definitions as ``(x, y, z, radius)`` tuples."""
        ...

    @blocking_pockets.setter
    def blocking_pockets(self, arg0: collections.abc.Sequence[tuple[float, float, float, float]]) -> None:
        """Set blocking-pocket definitions as ``(x, y, z, radius)`` tuples."""
        ...

    @property
    def average_rosenbluth_weights(self) -> PropertyWidom:
        """Return Widom insertion Rosenbluth-weight statistics."""
        ...

    @property
    def average_gibbs_rosenbluth_weights(self) -> PropertyWidom:
        """Return Gibbs-Widom insertion Rosenbluth-weight statistics."""
        ...

    @property
    def mc_moves_probabilities(self) -> MCMoveProbabilities:
        """Return configured Monte Carlo move probabilities."""
        ...

    @property
    def mc_moves_statistics(self) -> MCMoveStatistics:
        """Return Monte Carlo move acceptance/statistics counters."""
        ...

    @property
    def lambda_histogram(self) -> PropertyLambdaProbabilityHistogram:
        """Return lambda-probability histogram statistics."""
        ...

    # settings and results that are properties: a setting is made before the component is given to a system

    @property
    def type(self) -> 'Component.Type':
        """Adsorbate or cation."""
        ...

    @type.setter
    def type(self, value: 'Component.Type') -> None:
        ...

    @property
    def mol_fraction(self) -> float:
        """The mole fraction of the component in the reservoir [-]."""
        ...

    @mol_fraction.setter
    def mol_fraction(self, value: float) -> None:
        ...

    @property
    def fugacity_coefficient(self) -> float | None:
        """The fugacity coefficient [-]; None takes it from the equation of
        state."""
        ...

    @fugacity_coefficient.setter
    def fugacity_coefficient(self, value: float | None) -> None:
        ...

    @property
    def critical_temperature(self) -> float:
        """[K]"""
        ...

    @property
    def critical_pressure(self) -> float:
        """[Pa]"""
        ...

    @property
    def acentric_factor(self) -> float:
        """[-]"""
        ...

    @property
    def net_charge(self) -> float:
        """[e]"""
        ...

    @property
    def swappable(self) -> bool:
        """Whether a move of a system inserts or deletes the component."""
        ...

    @property
    def starting_bead(self) -> int:
        """The atom the growth of the molecule starts from."""
        ...

    @starting_bead.setter
    def starting_bead(self, value: int) -> None:
        ...

    @property
    def ideal_gas_rosenbluth_weight(self) -> float | None:
        """The Rosenbluth weight of the isolated molecule [-]; None computes
        it at the start of the run."""
        ...

    @ideal_gas_rosenbluth_weight.setter
    def ideal_gas_rosenbluth_weight(self, value: float | None) -> None:
        ...

    @property
    def ideal_gas_rosenbluth_weight_tolerance(self) -> float | None:
        """The relative standard error to which a weight that is not given
        is computed [-]."""
        ...

    @ideal_gas_rosenbluth_weight_tolerance.setter
    def ideal_gas_rosenbluth_weight_tolerance(self, value: float | None) -> None:
        ...

    @property
    def constraint_tolerance(self) -> float:
        """The relative error to which molecular dynamics restores a FIXED
        bond length with SHAKE/RATTLE [-]. Default 1e-8."""
        ...

    @constraint_tolerance.setter
    def constraint_tolerance(self, value: float) -> None:
        ...

    @property
    def constraint_maximum_iterations(self) -> int:
        """The iteration limit of SHAKE and of RATTLE per molecule and step.
        Default 500."""
        ...

    @constraint_maximum_iterations.setter
    def constraint_maximum_iterations(self, value: int) -> None:
        ...

    @property
    def distance_constraints(self) -> list[tuple[int, int, float]]:
        """The distance constraints (atom, atom, length [Å]) that molecular
        dynamics holds for the FIXED terms of the molecule."""
        ...

    @property
    def ideal_gas_energy(self) -> float | None:
        """The mean intramolecular energy of the isolated molecule [K], the
        reference of the enthalpy of adsorption."""
        ...

    @ideal_gas_energy.setter
    def ideal_gas_energy(self, value: float | None) -> None:
        ...

    @property
    def pair_component(self) -> int | None:
        """The index of the partner component of the pair-swap moves."""
        ...

    @pair_component.setter
    def pair_component(self, value: int | None) -> None:
        ...

    @property
    def maximum_pair_distance(self) -> float | None:
        """The largest separation of the pair of the pair-swap moves
        [Angstrom]."""
        ...

    @maximum_pair_distance.setter
    def maximum_pair_distance(self, value: float | None) -> None:
        ...

    @property
    def group_components(self) -> collections.abc.Sequence[int]:
        """The indices of the satellite components of the group-swap moves,
        one entry per molecule."""
        ...

    @group_components.setter
    def group_components(self, value: collections.abc.Sequence[int]) -> None:
        ...

    @property
    def maximum_group_distance(self) -> float | None:
        """The largest distance of a satellite from the central molecule of
        the group-swap moves [Angstrom]."""
        ...

    @maximum_group_distance.setter
    def maximum_group_distance(self, value: float | None) -> None:
        ...

    @property
    def tethered_proton_hop_sites(self) -> collections.abc.Sequence[collections.abc.Sequence[tuple[float, float, float]]]:
        """The sites of the tethered proton-hop move: per tethering site a
        list of fractional positions."""
        ...

    @tethered_proton_hop_sites.setter
    def tethered_proton_hop_sites(self, value: collections.abc.Sequence[collections.abc.Sequence[tuple[float, float, float]]]) -> None:
        ...

    @property
    def identity_changes(self) -> collections.abc.Sequence[int]:
        """The indices of the components the identity-change move changes
        the component into."""
        ...

    @identity_changes.setter
    def identity_changes(self, value: collections.abc.Sequence[int]) -> None:
        ...

    @property
    def identity_switches(self) -> collections.abc.Sequence[int]:
        """The indices of the partner components of the identity-switch
        move."""
        ...

    @identity_switches.setter
    def identity_switches(self, value: collections.abc.Sequence[int]) -> None:
        ...

    @property
    def gibbs_identity_changes(self) -> collections.abc.Sequence[int]:
        """The indices of the components the Gibbs identity-change move
        changes the component into."""
        ...

    @gibbs_identity_changes.setter
    def gibbs_identity_changes(self, value: collections.abc.Sequence[int]) -> None:
        ...

    @property
    def number_of_trial_molecules(self) -> int:
        """The number of complete trial molecules of the multiple-try moves."""
        ...

    @number_of_trial_molecules.setter
    def number_of_trial_molecules(self, value: int) -> None:
        ...

    @property
    def pivot_cbmc_number_of_trial_angles(self) -> int:
        """The number of trial angles of the configurational-bias pivot
        move."""
        ...

    @pivot_cbmc_number_of_trial_angles.setter
    def pivot_cbmc_number_of_trial_angles(self, value: int) -> None:
        ...

    @property
    def crankshaft_max_segment_size(self) -> int:
        """The largest number of atoms in the rotated segment of the
        crankshaft move."""
        ...

    @crankshaft_max_segment_size.setter
    def crankshaft_max_segment_size(self, value: int) -> None:
        ...

    @property
    def pivot_randomization_fraction(self) -> float:
        """The fraction of pivot attempts that draw the angle from the full
        circle."""
        ...

    @pivot_randomization_fraction.setter
    def pivot_randomization_fraction(self, value: float) -> None:
        ...

    @property
    def pivot_cbmc_randomization_fraction(self) -> float:
        """The fraction of configurational-bias pivot attempts that draw the
        trial angles from the full circle."""
        ...

    @pivot_cbmc_randomization_fraction.setter
    def pivot_cbmc_randomization_fraction(self, value: float) -> None:
        ...

    @property
    def crankshaft_randomization_fraction(self) -> float:
        """The fraction of crankshaft attempts that draw the angle from the
        full circle."""
        ...

    @crankshaft_randomization_fraction.setter
    def crankshaft_randomization_fraction(self, value: float) -> None:
        ...

    @property
    def concerted_rotation_randomization_fraction(self) -> float:
        """The fraction of concerted-rotation attempts that draw the driver
        angle from the full circle."""
        ...

    @concerted_rotation_randomization_fraction.setter
    def concerted_rotation_randomization_fraction(self, value: float) -> None:
        ...

    @property
    def bead_flip_randomization_fraction(self) -> float:
        """The fraction of bead-flip attempts that randomize the rotation
        completely."""
        ...

    @bead_flip_randomization_fraction.setter
    def bead_flip_randomization_fraction(self, value: float) -> None:
        ...

    @property
    def lambda_gibbs_histogram(self) -> PropertyLambdaProbabilityHistogram:
        """(self) -> raspalib._raspalib.PropertyLambdaProbabilityHistogram"""
        ...

    @property
    def lambda_pair_swap_histogram(self) -> PropertyLambdaProbabilityHistogram:
        """(self) -> raspalib._raspalib.PropertyLambdaProbabilityHistogram"""
        ...

    @property
    def lambda_group_swap_histogram(self) -> PropertyLambdaProbabilityHistogram:
        """(self) -> raspalib._raspalib.PropertyLambdaProbabilityHistogram"""
        ...

    @property
    def number_of_blocks(self) -> int:
        """The number of blocks the statistics of the component are kept in."""
        ...

    @property
    def rigid(self) -> bool:
        """Whether the molecule is a single rigid body."""
        ...

    @property
    def connectivity_table(self) -> ConnectivityTable:
        """The connectivity of the molecule, made from its bonds."""
        ...

    @property
    def intra_molecular_potentials(self) -> IntraMolecularPotentials:
        """The bonded terms of the molecule."""
        ...

    # A molecule without bonds is one rigid body. With bond potentials it is
    # flexible: the bonds are its connectivity, and every bond, bend and torsion
    # they form needs a potential. With ``rigid_bodies`` the listed groups of
    # atoms move as one unit and the other atoms are flexible; the
    # ``connectivity_table`` is then given, and the terms inside one rigid body
    # need no potential. ``chiral_centers`` are tuples of a center atom and
    # three of its bonded neighbors. The 1-4 pairs do not interact unless they
    # are given a scaling. ``partial_reinsertion_fixed_atoms`` are the sets of
    # atoms a partial reinsertion keeps in place (without them the move has
    # nothing to do); ``repeat_units`` are the monomers of a chain, for the
    # reptation move; ``end_to_end_atoms`` are the two atoms of the end-to-end
    # distance, by default the two furthest apart along the bonds.

    @property
    def rigid_bodies(self) -> collections.abc.Sequence[collections.abc.Sequence[int]]:
        """The groups of atoms that move as one rigid unit."""
        ...

    @property
    def chiral_centers(self) -> collections.abc.Sequence[ChiralCenter]:
        """The chiral centers, with the handedness they have in the positions
        of the atoms."""
        ...

    @property
    def intra_14_van_der_waals_scaling(self) -> float:
        """The scaling of the Van der Waals interaction of the 1-4 pairs."""
        ...

    @property
    def intra_14_charge_charge_scaling(self) -> float:
        """The scaling of the Coulomb interaction of the 1-4 pairs."""
        ...

    @property
    def partial_reinsertion_fixed_atoms(self) -> collections.abc.Sequence[collections.abc.Sequence[int]]:
        """The sets of atoms that stay where they are in a partial
        reinsertion; the others are grown again."""
        ...

    @property
    def repeat_units(self) -> collections.abc.Sequence[collections.abc.Sequence[int]]:
        """The repeat units of a chain molecule, for the reptation move."""
        ...

    @property
    def end_to_end_atoms(self) -> collections.abc.Sequence[int] | None:
        """The two atoms whose distance is the end-to-end distance of the
        molecule."""
        ...

    def backbone_atoms(self) -> collections.abc.Sequence[int]:
        """Return the atoms on the shortest path of bonds between the
        end-to-end atoms."""
        ...

    @property
    def molecule_properties_settings(self) -> MoleculePropertiesSettings | None:
        """Histograms of the bond lengths, bend angles, torsion angles and the end-to-end distance."""
        ...

    @molecule_properties_settings.setter
    def molecule_properties_settings(self, value: MoleculePropertiesSettings | None) -> None:
        ...

    @property
    def molecule_shape_settings(self) -> MoleculeShapeSettings | None:
        """The radius of gyration and the shape descriptors of the gyration tensor."""
        ...

    @molecule_shape_settings.setter
    def molecule_shape_settings(self, value: MoleculeShapeSettings | None) -> None:
        ...

    @property
    def molecule_backbone_settings(self) -> MoleculeBackboneSettings | None:
        """The statistics of the backbone of a chain."""
        ...

    @molecule_backbone_settings.setter
    def molecule_backbone_settings(self, value: MoleculeBackboneSettings | None) -> None:
        ...

    @property
    def end_to_end_acf_settings(self) -> EndToEndACFSettings | None:
        """The autocorrelation function of the end-to-end vector."""
        ...

    @end_to_end_acf_settings.setter
    def end_to_end_acf_settings(self, value: EndToEndACFSettings | None) -> None:
        ...

    @property
    def ln_partition_function(self) -> float:
        """The logarithm of the partition function of the isolated molecule,
        for the reaction moves."""
        ...

    @ln_partition_function.setter
    def ln_partition_function(self, value: float) -> None:
        ...

    def set_ln_partition_function_species(self, species: str, database: ThermochemicalDatabase = ...) -> None:
        """Take the partition function from the thermochemical data of a
        species, at the temperature of the system the component is given
        to."""
        ...

    @property
    def ln_partition_function_species(self) -> str | None: ...
