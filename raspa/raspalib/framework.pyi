import collections.abc
import enum

from raspalib.atom import *
from raspalib.simulation_box import *
from raspalib.force_field import *
from raspalib.connectivity_table import *
from raspalib.intra_molecular_potentials import *


class FrameworkGroupType(enum.IntEnum):
    """How the atoms of a group of a flexible framework move."""
    FIXED: int = 0
    RIGID: int = 1
    FLEXIBLE: int = 2


class FrameworkDefinition():
    """
    The definition of a flexible framework, as a framework-definition file
    holds it.

    A potential is given per combination of pseudo-atoms and applies to every
    bond, bend or torsion of those pseudo-atoms the framework has.
    """

    def __init__(
        self,
        flexible: bool = True,
        bonds: collections.abc.Sequence[tuple[collections.abc.Sequence[str], BondType, collections.abc.Sequence[float]]] = [],
        bends: collections.abc.Sequence[tuple[collections.abc.Sequence[str], BendType, collections.abc.Sequence[float]]] = [],
        torsions: collections.abc.Sequence[tuple[collections.abc.Sequence[str], TorsionType, collections.abc.Sequence[float]]] = [],
        improper_torsions: collections.abc.Sequence[tuple[collections.abc.Sequence[str], TorsionType, collections.abc.Sequence[float]]] = [],
        groups: collections.abc.Sequence[tuple[FrameworkGroupType, collections.abc.Sequence[int]]] = [],
        exclude_intra_12_interactions: bool = True,
        exclude_intra_13_interactions: bool = True,
        exclude_intra_bond_interactions: bool = False,
        exclude_intra_bend_interactions: bool = False,
        intra_14_van_der_waals_scaling: float = 0.0,
        intra_14_charge_charge_scaling: float = 0.0,
    ) -> None:
        """
        Initialize a :class:`FrameworkDefinition`.

        Args:
            flexible: Whether the atoms of the framework move.
            bonds: The bond potentials: the names of the two pseudo-atoms, the
                type and the parameters of that type.
            bends: The bend potentials, of three pseudo-atoms.
            torsions: The torsion potentials, of four pseudo-atoms.
            improper_torsions: The improper-torsion potentials.
            groups: Groups of atoms of the framework that are fixed, rigid or
                flexible: the type and the indices of the atoms.
            exclude_intra_12_interactions: Leave out the non-bonded interaction
                of bonded atoms.
            exclude_intra_13_interactions: Leave out the non-bonded interaction
                of atoms two bonds apart.
            exclude_intra_bond_interactions: Leave out the non-bonded
                interaction of the atoms of a bond potential.
            exclude_intra_bend_interactions: Leave out the non-bonded
                interaction of the outer atoms of a bend potential.
            intra_14_van_der_waals_scaling: Scaling of the Van der Waals
                interaction of atoms three bonds apart.
            intra_14_charge_charge_scaling: Scaling of the Coulomb interaction
                of atoms three bonds apart.
        """
        ...


class Framework():
    """
    Represent a (typically porous) host framework in RASPA.

    A framework combines structural information (unit-cell geometry, symmetry,
    atoms, and replication count) with the force-field settings used to model
    framework interactions. It is made atom by atom, from a CIF file
    (``from_cif``), from the content of one (``from_cif_string``), or it is one
    of the zeolites that are built in (``make_mfi`` and others). With a
    definition it is flexible.
    """

    class UseChargesFrom(enum.IntEnum):
        """Where the charges of the atoms of a CIF file come from."""
        PSEUDO_ATOMS: int = 0
        CIF_FILE: int = 1
        CHARGE_EQUILIBRATION: int = 2

    def __init__(
        self,
        force_field: ForceField,
        component_name: str,
        simulation_box: SimulationBox,
        space_group_hall_number: int,
        defined_atoms: collections.abc.Sequence[Atom],
        number_of_unit_cells: collections.abc.Sequence[int],
        definition: FrameworkDefinition | None = None,
        definition_file: str | None = None,
    ) -> None:
        """
        Initialize a :class:`Framework`.

        Args:
            force_field: Force-field definition used for the framework.
            component_name: Name/identifier of the framework component.
            simulation_box: Unit-cell box definition.
            space_group_hall_number: Space-group Hall number used for
                crystallographic symmetry.
            defined_atoms: The atoms of the asymmetric unit, in fractional
                positions.
            number_of_unit_cells: Number of replicated unit cells along each
                lattice direction.
            definition: The definition of a flexible framework.
            definition_file: The name of a framework-definition file, instead
                of a definition.
        """
        ...

    @staticmethod
    def from_cif(
        force_field: ForceField,
        file_name: str,
        number_of_unit_cells: collections.abc.Sequence[int] | None = None,
        use_charges_from: 'Framework.UseChargesFrom' = ...,
        component_name: str | None = None,
        definition: FrameworkDefinition | None = None,
        definition_file: str | None = None,
    ) -> 'Framework':
        """
        Make a framework from a CIF file.

        Args:
            force_field: Force-field definition used for the framework.
            file_name: The CIF file; ``.cif`` may be left out, and the file is
                also looked for in ``RASPA_DIR``.
            number_of_unit_cells: Number of unit cells along each lattice
                direction; ``None`` takes the smallest supercell in which the
                cutoffs of the force field fit.
            use_charges_from: The charges of the pseudo-atoms (the default),
                those of the file, or charge equilibration.
            component_name: The name of the framework; that of the file unless
                given.
            definition: The definition of a flexible framework.
            definition_file: The name of a framework-definition file, instead
                of a definition.
        """
        ...

    @staticmethod
    def from_cif_string(
        force_field: ForceField,
        cif: str,
        component_name: str,
        number_of_unit_cells: collections.abc.Sequence[int] | None = None,
        use_charges_from: 'Framework.UseChargesFrom' = ...,
        definition: FrameworkDefinition | None = None,
        definition_file: str | None = None,
    ) -> 'Framework':
        """Make a framework from the content of a CIF file, as ``from_cif``
        makes it from a file."""
        ...

    @staticmethod
    def make_mfi(force_field: ForceField, number_of_unit_cells: collections.abc.Sequence[int] = (1, 1, 1)) -> 'Framework':
        """The all-silica zeolite MFI; the force field needs the pseudo-atoms
        ``Si`` and ``O``."""
        ...

    @staticmethod
    def make_itq_29(force_field: ForceField, number_of_unit_cells: collections.abc.Sequence[int] = (1, 1, 1)) -> 'Framework':
        """The all-silica zeolite ITQ-29."""
        ...

    @staticmethod
    def make_fau(force_field: ForceField, number_of_unit_cells: collections.abc.Sequence[int] = (1, 1, 1)) -> 'Framework':
        """The all-silica zeolite FAU."""
        ...

    @staticmethod
    def make_cha(force_field: ForceField, number_of_unit_cells: collections.abc.Sequence[int] = (1, 1, 1)) -> 'Framework':
        """The all-silica zeolite CHA, in its rhombohedral cell."""
        ...

    @property
    def name(self) -> str: ...

    @property
    def rigid(self) -> bool:
        """Whether no atom of the framework moves."""
        ...

    @property
    def simulation_box(self) -> SimulationBox:
        """The unit cell."""
        ...

    @property
    def space_group_hall_number(self) -> int: ...

    @property
    def number_of_unit_cells(self) -> tuple[int, int, int]: ...

    @property
    def mass(self) -> float:
        """The mass of the framework [amu]."""
        ...

    @property
    def unit_cell_mass(self) -> float:
        """The mass of one unit cell [amu]."""
        ...

    @property
    def net_charge(self) -> float:
        """The net charge of the framework [e]."""
        ...

    @property
    def defined_atoms(self) -> collections.abc.Sequence[Atom]:
        """The atoms of the asymmetric unit, in fractional positions."""
        ...

    @property
    def unit_cell_atoms(self) -> collections.abc.Sequence[Atom]:
        """The atoms of one unit cell, in Cartesian positions."""
        ...

    @property
    def atoms(self) -> collections.abc.Sequence[Atom]:
        """The atoms of the supercell, in Cartesian positions."""
        ...

    @property
    def number_of_atoms(self) -> int: ...

    @property
    def number_of_fixed_atoms(self) -> int: ...

    @property
    def number_of_mobile_atoms(self) -> int: ...

    @property
    def number_of_rigid_groups(self) -> int: ...

    @property
    def connectivity_table(self) -> ConnectivityTable:
        """The bonds of a flexible framework."""
        ...

    @property
    def intra_molecular_potentials(self) -> IntraMolecularPotentials:
        """The bonded terms of a flexible framework, one for every bond, bend
        and torsion."""
        ...

    def print(self, force_field: ForceField) -> str: ...
