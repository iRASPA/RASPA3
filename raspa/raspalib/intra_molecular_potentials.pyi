import collections.abc
import enum

from raspalib.vdw_parameters import *

# The bonded terms of a molecule have one shape: the atoms of the term (indices
# into the atoms of the component), a type, and the parameters of that type in
# the units of the molecule files: energies in Kelvin, lengths in Angstrom,
# angles in degrees.

class BondType(enum.IntEnum):
    """The forms of a bond."""
    NONE: int = 0
    FIXED: int = 1
    HARMONIC: int = 2
    CORE_SHELL_SPRING: int = 3
    MORSE: int = 4
    LJ_12_6: int = 5
    LENNARD_JONES: int = 6
    BUCKINGHAM: int = 7
    RESTRAINED_HARMONIC: int = 8
    QUARTIC: int = 9
    CFF_QUARTIC: int = 10
    MM3: int = 11


class BondPotential():
    """A bond between 2 atoms of a molecule.

    The types and their parameters:
      NONE: no parameters
      FIXED: p_0 [Å]
      HARMONIC: p_0/k_B [K/Å^2], p_1 [Å]
      CORE_SHELL_SPRING: p_0/k_B [K/Å^2]
      MORSE: p_0/k_B [K/Å^2], p_1 [Å^-1], p_2 [Å]
      LJ_12_6: p_0/k_B [K Å^12], p_1/k_B [K Å^6]
      LENNARD_JONES: p_0/k_B [K], p_1 [Å]
      BUCKINGHAM: p_0/k_B [K], p_1 [Å^-1], p_2/k_B [K Å^6]
      RESTRAINED_HARMONIC: p_0/k_B [K/Å^2], p_1 [Å], p_2 [Å]
      QUARTIC: p_0/k_B [K/Å^2], p_1 [Å], p_2 [K/Å^3], p_3 [K/Å^4]
      CFF_QUARTIC: p_0/k_B [K/Å^2], p_1 [Å], p_2 [K/Å^3], p_3 [K/Å^4]
      MM3: p_0/k_B [mdyne/Å molecule], p_1 [Å]
    """

    def __init__(self, identifiers: collections.abc.Sequence[int], type: BondType,
                 parameters: collections.abc.Sequence[float] = []) -> None:
        """
        Initialize a :class:`BondPotential`.

        Args:
            identifiers: The 2 atoms of the term.
            type: The form of the potential.
            parameters: The parameters of that form.
        """
        ...

    @staticmethod
    def describe_parameters() -> str:
        """Return the parameters of every type, with their units."""
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]:
        """The atoms of the term, as indices into the atoms of the component."""
        ...

    @property
    def type(self) -> BondType:
        """The form of the potential."""
        ...


class UreyBradleyType(enum.IntEnum):
    """The forms of a Urey-Bradley term."""
    FIXED: int = 0
    HARMONIC: int = 1
    CORE_SHELL_SPRING: int = 2
    MORSE: int = 3
    LJ_12_6: int = 4
    LENNARD_JONES: int = 5
    BUCKINGHAM: int = 6
    RESTRAINED_HARMONIC: int = 7
    QUARTIC: int = 8
    CFF_QUARTIC: int = 9
    MM3: int = 10


class UreyBradleyPotential():
    """A Urey-Bradley term between 2 atoms of a molecule.

    The types and their parameters:
      FIXED: 1 parameters
      HARMONIC: p_0/k_B [K/Å^2], p_1 [Å]
      CORE_SHELL_SPRING: p_0/k_B [K/Å^2]
      MORSE: p_0/k_B [K/Å^2], p_1 [Å^-1], p_2 [Å]
      LJ_12_6: p_0/k_B [K Å^12], p_1/k_B [K Å^6]
      LENNARD_JONES: p_0/k_B [K], p_1 [Å]
      BUCKINGHAM: p_0/k_B [K], p_1 [Å^-1], p_2/k_B [K Å^6]
      RESTRAINED_HARMONIC: p_0/k_B [K/Å^2], p_1 [Å], p_2 [Å]
      QUARTIC: p_0/k_B [K/Å^2], p_1 [Å], p_2 [K/Å^3], p_3 [K/Å^4]
      CFF_QUARTIC: p_0/k_B [K/Å^2], p_1 [Å], p_2 [K/Å^3], p_3 [K/Å^4]
      MM3: p_0/k_B [mdyne/Å molecule], p_1 [Å]
    """

    def __init__(self, identifiers: collections.abc.Sequence[int], type: UreyBradleyType,
                 parameters: collections.abc.Sequence[float] = []) -> None:
        """
        Initialize a :class:`UreyBradleyPotential`.

        Args:
            identifiers: The 2 atoms of the term.
            type: The form of the potential.
            parameters: The parameters of that form.
        """
        ...

    @staticmethod
    def describe_parameters() -> str:
        """Return the parameters of every type, with their units."""
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]:
        """The atoms of the term, as indices into the atoms of the component."""
        ...

    @property
    def type(self) -> UreyBradleyType:
        """The form of the potential."""
        ...


class BendType(enum.IntEnum):
    """The forms of a bend."""
    RIGID: int = 0
    FIXED: int = 1
    HARMONIC: int = 2
    CORE_SHELL: int = 3
    QUARTIC: int = 4
    CFF_QUARTIC: int = 5
    HARMONIC_COSINE: int = 6
    COSINE: int = 7
    TAFIPOLSKY: int = 8
    MM3: int = 9


class BendPotential():
    """A bend between 3 atoms of a molecule.

    The types and their parameters:
      RIGID: 1 parameters
      FIXED: 1 parameters
      HARMONIC: p_0/k_B [K/rad^2], p_1 [degrees]
      CORE_SHELL: p_0/k_B [K/rad^2], p_1 [degrees]
      QUARTIC: p_0/k_B [K/rad^2], p_1 [degrees], p_2/k_B [K/rad^3], p_3/k_B [K/rad^4]
      CFF_QUARTIC: p_0/k_B [K/rad^2], p_1 [degrees], p_2/k_B [K/rad^3], p_3/k_B [K/rad^4]
      HARMONIC_COSINE: p_0/k_B [K], p_1 [degrees]
      COSINE: p_0/k_B [K], p_1 [-], p_2 [degrees]
      TAFIPOLSKY: p_0/k_B [K]
      MM3: p_0/k_B [mdyne A/rad^2], p_1 [degrees]
    """

    def __init__(self, identifiers: collections.abc.Sequence[int], type: BendType,
                 parameters: collections.abc.Sequence[float] = []) -> None:
        """
        Initialize a :class:`BendPotential`.

        Args:
            identifiers: The 3 atoms of the term.
            type: The form of the potential.
            parameters: The parameters of that form.
        """
        ...

    @staticmethod
    def describe_parameters() -> str:
        """Return the parameters of every type, with their units."""
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]:
        """The atoms of the term, as indices into the atoms of the component."""
        ...

    @property
    def type(self) -> BendType:
        """The form of the potential."""
        ...


class InversionBendType(enum.IntEnum):
    """The forms of a inversion bend."""
    HARMONIC: int = 0
    HARMONIC_COSINE: int = 1
    PLANAR: int = 2
    HARMONIC2: int = 3
    HARMONIC_COSINE2: int = 4
    PLANAR2: int = 5
    MM3: int = 6


class InversionBendPotential():
    """A inversion bend between 4 atoms of a molecule.

    The types and their parameters:
      HARMONIC: p_0/k_B [K/rad^2], p_1 [degrees]
      HARMONIC_COSINE: p_0/k_B [K/rad^2], p_1 [degrees]
      PLANAR: p_0/k_B [K]
      HARMONIC2: p_0/k_B [K/rad^2], p_1 [degrees]
      HARMONIC_COSINE2: p_0/k_B [K/rad^2], p_1 [degrees]
      PLANAR2: p_0/k_B [K]
      MM3: p_0/k_B [mdyne A/rad^2], p_1 [degrees]
    """

    def __init__(self, identifiers: collections.abc.Sequence[int], type: InversionBendType,
                 parameters: collections.abc.Sequence[float] = []) -> None:
        """
        Initialize a :class:`InversionBendPotential`.

        Args:
            identifiers: The 4 atoms of the term.
            type: The form of the potential.
            parameters: The parameters of that form.
        """
        ...

    @staticmethod
    def describe_parameters() -> str:
        """Return the parameters of every type, with their units."""
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]:
        """The atoms of the term, as indices into the atoms of the component."""
        ...

    @property
    def type(self) -> InversionBendType:
        """The form of the potential."""
        ...


class TorsionType(enum.IntEnum):
    """The forms of a torsion."""
    FIXED: int = 0
    HARMONIC: int = 1
    HARMONIC_COSINE: int = 2
    THREE_COSINE: int = 3
    RYCKAERT_BELLEMANS: int = 4
    TRAPPE: int = 5
    TRAPPE_EXTENDED: int = 6
    MOD_TRAPPE: int = 7
    CVFF: int = 8
    CFF: int = 9
    CFF2: int = 10
    OPLS: int = 11
    MM3: int = 12
    FOURIER_SERIES: int = 13
    FOURIER_SERIES2: int = 14
    CVFF_BLOCKED: int = 15
    POLYNOMIAL: int = 16


class TorsionPotential():
    """A torsion between 4 atoms of a molecule.

    The types and their parameters:
      FIXED: 1 parameters
      HARMONIC: p_0/k_B [K/rad^2], p_1 [degrees]
      HARMONIC_COSINE: p_0/k_B [K], p_1 [degrees]
      THREE_COSINE: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K]
      RYCKAERT_BELLEMANS: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K], p_3/k_B [K], p_4/k_B [K], p_5/k_B [K]
      TRAPPE: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K], p_3/k_B [K]
      TRAPPE_EXTENDED: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K], p_3/k_B [K], p_4/k_B [K]
      MOD_TRAPPE: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K], p_3/k_B [K], p_4 [degrees]
      CVFF: p_0/k_B [K], p_1 [-], p_2 [degrees]
      CFF: p_0/k_B [K], p_1 [-], p_2/k_B [K]
      CFF2: p_0/k_B [K], p_1 [-], p_2/k_B [K]
      OPLS: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K], p_3/k_B [K]
      MM3: p_0/k_B [Kcal/mol], p_1/k_B [Kcal/mol], p_2/k_B [Kcal/mol]
      FOURIER_SERIES: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K], p_3/k_B [K], p_4/k_B [K], p_5/k_B [K]
      FOURIER_SERIES2: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K], p_3/k_B [K], p_4/k_B [K], p_5/k_B [K]
      CVFF_BLOCKED: p_0 [rad], p_1/k_B [K], p_2 [-], p_3 [rad], p_4 [rad]
      POLYNOMIAL: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K], p_3/k_B [K], p_4/k_B [K], p_5/k_B [K]
    """

    def __init__(self, identifiers: collections.abc.Sequence[int], type: TorsionType,
                 parameters: collections.abc.Sequence[float] = []) -> None:
        """
        Initialize a :class:`TorsionPotential`.

        Args:
            identifiers: The 4 atoms of the term.
            type: The form of the potential.
            parameters: The parameters of that form.
        """
        ...

    @staticmethod
    def describe_parameters() -> str:
        """Return the parameters of every type, with their units."""
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]:
        """The atoms of the term, as indices into the atoms of the component."""
        ...

    @property
    def type(self) -> TorsionType:
        """The form of the potential."""
        ...


class BondBondType(enum.IntEnum):
    """The forms of a bond-bond cross term."""
    CVFF: int = 0
    CFF: int = 1


class BondBondPotential():
    """A bond-bond cross term between 3 atoms of a molecule.

    The types and their parameters:
      CVFF: p_0/k_B [K/Å^2], p_1 [Å], p_2 [Å]
      CFF: p_0/k_B [K/Å^2], p_1 [Å], p_2 [Å]
    """

    def __init__(self, identifiers: collections.abc.Sequence[int], type: BondBondType,
                 parameters: collections.abc.Sequence[float] = []) -> None:
        """
        Initialize a :class:`BondBondPotential`.

        Args:
            identifiers: The 3 atoms of the term.
            type: The form of the potential.
            parameters: The parameters of that form.
        """
        ...

    @staticmethod
    def describe_parameters() -> str:
        """Return the parameters of every type, with their units."""
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]:
        """The atoms of the term, as indices into the atoms of the component."""
        ...

    @property
    def type(self) -> BondBondType:
        """The form of the potential."""
        ...


class BondBendType(enum.IntEnum):
    """The forms of a bond-bend cross term."""
    CVFF: int = 0
    CFF: int = 1
    MM3: int = 2
    TRUNCATED_HARMONIC: int = 3
    SCREENED_HARMONIC: int = 4
    SCREENED_VESSAL: int = 5
    TRUNCATED_VESSAL: int = 6


class BondBendPotential():
    """A bond-bend cross term between 3 atoms of a molecule.

    The types and their parameters:
      CVFF: p_0 [degrees], p_1/k_B [K/A/rad], p_2 [A], p_3/k_B [K/A/rad], p_4 [A]
      CFF: p_0 [degrees], p_1/k_B [K/A/rad], p_2 [A], p_3/k_B [K/A/rad], p_4 [A]
      MM3: p_0 [mdyne/rad], p_1 [A], p_2 [A], p_3 [degrees]
      TRUNCATED_HARMONIC: p_0/k_B [K/A/rad], p_2 [degrees], p_3 [A]
      SCREENED_HARMONIC: p_0/k_B [K/rad^2], p_1 [degrees], p_2 [A], p_3 [A]
      SCREENED_VESSAL: p_0/k_B [K/rad^2], p_1 [degrees], p_2 [A], p_3 [A]
      TRUNCATED_VESSAL: p_0/k_B [K/rad^(4+p_2)], p_1 [degrees], p_2 [-], p_3 [A]
    """

    def __init__(self, identifiers: collections.abc.Sequence[int], type: BondBendType,
                 parameters: collections.abc.Sequence[float] = []) -> None:
        """
        Initialize a :class:`BondBendPotential`.

        Args:
            identifiers: The 3 atoms of the term.
            type: The form of the potential.
            parameters: The parameters of that form.
        """
        ...

    @staticmethod
    def describe_parameters() -> str:
        """Return the parameters of every type, with their units."""
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]:
        """The atoms of the term, as indices into the atoms of the component."""
        ...

    @property
    def type(self) -> BondBendType:
        """The form of the potential."""
        ...


class BondTorsionType(enum.IntEnum):
    """The forms of a bond-torsion cross term."""
    MM3: int = 0


class BondTorsionPotential():
    """A bond-torsion cross term between 4 atoms of a molecule.

    The types and their parameters:
      MM3: p_0 [kcal/A mole], p_1 [kcal/A mole], p_2 [kcal/A mole], p_3/k_B [A]
    """

    def __init__(self, identifiers: collections.abc.Sequence[int], type: BondTorsionType,
                 parameters: collections.abc.Sequence[float] = []) -> None:
        """
        Initialize a :class:`BondTorsionPotential`.

        Args:
            identifiers: The 4 atoms of the term.
            type: The form of the potential.
            parameters: The parameters of that form.
        """
        ...

    @staticmethod
    def describe_parameters() -> str:
        """Return the parameters of every type, with their units."""
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]:
        """The atoms of the term, as indices into the atoms of the component."""
        ...

    @property
    def type(self) -> BondTorsionType:
        """The form of the potential."""
        ...


class BendBendType(enum.IntEnum):
    """The forms of a bend-bend cross term."""
    CVFF: int = 0
    CFF: int = 1
    MM3: int = 2


class BendBendPotential():
    """A bend-bend cross term between 4 atoms of a molecule.

    The types and their parameters:
      CVFF: p_0/k_B [K/rad^2], p_1 [degrees], p_2 [degrees]
      CFF: p_0/k_B [K/rad^2], p_1 [degrees], p_2 [degrees]
      MM3: p_0 [mdyne A/rad^2], p_1 [degrees], p_2 [degrees]
    """

    def __init__(self, identifiers: collections.abc.Sequence[int], type: BendBendType,
                 parameters: collections.abc.Sequence[float] = []) -> None:
        """
        Initialize a :class:`BendBendPotential`.

        Args:
            identifiers: The 4 atoms of the term.
            type: The form of the potential.
            parameters: The parameters of that form.
        """
        ...

    @staticmethod
    def describe_parameters() -> str:
        """Return the parameters of every type, with their units."""
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]:
        """The atoms of the term, as indices into the atoms of the component."""
        ...

    @property
    def type(self) -> BendBendType:
        """The form of the potential."""
        ...


class BendTorsionType(enum.IntEnum):
    """The forms of a bend-torsion cross term."""
    SMOOTHED: int = 0
    SMOOTHED_THREE_COSINE: int = 1
    NICHOLAS: int = 2
    CFF: int = 3
    SMOOTHED_CFF: int = 4
    SMOOTHED_CFF2: int = 5
    SMOOTHED_CFF3: int = 6
    CVFF: int = 7


class BendTorsionPotential():
    """A bend-torsion cross term between 4 atoms of a molecule.

    The types and their parameters:
      SMOOTHED: p_0/k_B [K/rad^2], p_1 [-], p_2 [degrees]
      SMOOTHED_THREE_COSINE: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K]
      NICHOLAS: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K]
      CFF: p_0/k_B [K/rad^3], p_1 [degrees], p_2 [degrees]
      SMOOTHED_CFF: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K]
      SMOOTHED_CFF2: p_0/k_B [K], p_1/k_B [K], p_2/k_B [K]
      SMOOTHED_CFF3: p_0/k_B [K/rad^3], p_1 [degrees], p_2 [degrees]
      CVFF: p_0/k_B [K/rad^3], p_1 [degrees], p_2 [degrees]
    """

    def __init__(self, identifiers: collections.abc.Sequence[int], type: BendTorsionType,
                 parameters: collections.abc.Sequence[float] = []) -> None:
        """
        Initialize a :class:`BendTorsionPotential`.

        Args:
            identifiers: The 4 atoms of the term.
            type: The form of the potential.
            parameters: The parameters of that form.
        """
        ...

    @staticmethod
    def describe_parameters() -> str:
        """Return the parameters of every type, with their units."""
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]:
        """The atoms of the term, as indices into the atoms of the component."""
        ...

    @property
    def type(self) -> BendTorsionType:
        """The form of the potential."""
        ...


class VanDerWaalsPotential():
    """A pair of atoms inside a molecule whose Van der Waals interaction is given."""

    def __init__(self, identifiers: collections.abc.Sequence[int], parameters: VDWParameters,
                 scaling: float = 1.0) -> None:
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]: ...
    @property
    def scaling(self) -> float: ...


class CoulombPotential():
    """A pair of atoms inside a molecule whose Coulomb interaction is given."""

    def __init__(self, identifiers: collections.abc.Sequence[int], charge_a: float, charge_b: float,
                 scaling: float = 1.0) -> None:
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]: ...
    @property
    def scaling(self) -> float: ...


class CMAPMap():
    """
    An energy map over two consecutive backbone dihedral angles, the CMAP
    correction of protein force fields.
    """

    def __init__(self, name: str, resolution: int, energies: collections.abc.Sequence[float]) -> None:
        """
        Initialize a :class:`CMAPMap`.

        Args:
            name: The name of the map.
            resolution: The number of grid points per angle (at least 3).
            energies: The resolution x resolution energies [K]; the first angle
                is the slow index.
        """
        ...

    @property
    def name(self) -> str: ...
    @property
    def resolution(self) -> int: ...


class CMAPPotential():
    """A CMAP term of five atoms A, B, C, D, E with one of the maps."""

    def __init__(self, identifiers: collections.abc.Sequence[int], map_index: int) -> None:
        """
        Initialize a :class:`CMAPPotential`.

        Args:
            identifiers: The five atoms of the term.
            map_index: The index of the map in ``cmap_maps`` of the potentials.
        """
        ...

    @property
    def identifiers(self) -> collections.abc.Sequence[int]: ...
    @property
    def map_index(self) -> int: ...


class ChiralCenter():
    """
    A chiral center of a component: an atom with three of its bonded neighbors.

    Chiral centers are given to :class:`Component` as tuples of four atoms; the
    handedness follows from the positions of the atoms.
    """

    class Chirality(enum.IntEnum):
        """The handedness of a chiral center."""
        S: int = 0
        R: int = 1

    @property
    def identifiers(self) -> collections.abc.Sequence[int]:
        """The center atom and three of its bonded neighbors."""
        ...

    @property
    def chirality(self) -> 'ChiralCenter.Chirality':
        """The handedness the atoms have in the positions of the component."""
        ...


class IntraMolecularPotentials():
    """
    Store intra-molecular potential definitions for a component in RASPA.

    This container holds the bonded terms of a flexible molecule. The bonds are
    the connectivity of the molecule: every bond, bend and torsion they form
    needs a potential, which :class:`Component` checks.
    """

    def __init__(
        self,
        bonds: collections.abc.Sequence[BondPotential] = [],
        bends: collections.abc.Sequence[BendPotential] = [],
        torsions: collections.abc.Sequence[TorsionPotential] = [],
        improper_torsions: collections.abc.Sequence[TorsionPotential] = [],
        urey_bradleys: collections.abc.Sequence[UreyBradleyPotential] = [],
        inversion_bends: collections.abc.Sequence[InversionBendPotential] = [],
        bond_bonds: collections.abc.Sequence[BondBondPotential] = [],
        bond_bends: collections.abc.Sequence[BondBendPotential] = [],
        bond_torsions: collections.abc.Sequence[BondTorsionPotential] = [],
        bend_bends: collections.abc.Sequence[BendBendPotential] = [],
        bend_torsions: collections.abc.Sequence[BendTorsionPotential] = [],
        van_der_waals: collections.abc.Sequence[VanDerWaalsPotential] = [],
        coulombs: collections.abc.Sequence[CoulombPotential] = [],
        cmap_maps: collections.abc.Sequence[CMAPMap] = [],
        cmaps: collections.abc.Sequence[CMAPPotential] = [],
    ) -> None:
        """Initialize an :class:`IntraMolecularPotentials`."""
        ...

    @property
    def bonds(self) -> collections.abc.Sequence[BondPotential]: ...
    @property
    def bends(self) -> collections.abc.Sequence[BendPotential]: ...
    @property
    def torsions(self) -> collections.abc.Sequence[TorsionPotential]: ...
    @property
    def improper_torsions(self) -> collections.abc.Sequence[TorsionPotential]: ...
    @property
    def urey_bradleys(self) -> collections.abc.Sequence[UreyBradleyPotential]: ...
    @property
    def inversion_bends(self) -> collections.abc.Sequence[InversionBendPotential]: ...
    @property
    def bond_bonds(self) -> collections.abc.Sequence[BondBondPotential]: ...
    @property
    def bond_bends(self) -> collections.abc.Sequence[BondBendPotential]: ...
    @property
    def bond_torsions(self) -> collections.abc.Sequence[BondTorsionPotential]: ...
    @property
    def bend_bends(self) -> collections.abc.Sequence[BendBendPotential]: ...
    @property
    def bend_torsions(self) -> collections.abc.Sequence[BendTorsionPotential]: ...
    @property
    def van_der_waals(self) -> collections.abc.Sequence[VanDerWaalsPotential]: ...
    @property
    def coulombs(self) -> collections.abc.Sequence[CoulombPotential]: ...
    @property
    def cmap_maps(self) -> collections.abc.Sequence[CMAPMap]: ...
    @property
    def cmaps(self) -> collections.abc.Sequence[CMAPPotential]: ...
    @property
    def chiral_centers(self) -> collections.abc.Sequence[ChiralCenter]: ...
