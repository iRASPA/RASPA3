import collections.abc
import enum
import typing


class VDWParameters():
    """
    Store the parameters of a Van der Waals pair potential.

    A potential is a Lennard-Jones potential given by its strength and size, or
    a potential of another type given by the parameters of that type.
    """

    class Type(enum.IntEnum):
        """The form of the pair potential."""
        NONE: int = 0
        LENNARD_JONES: int = 1
        BUCKINGHAM: int = 2
        MORSE: int = 3
        FEYNMAN_HIBBS: int = 4
        MM3: int = 5
        BORN_HUGGINS_MEYER: int = 6
        LENNARD_JONES_SHIFTED_FORCE: int = 7
        POTENTIAL_12_6: int = 8
        POTENTIAL_12_6_2_0: int = 9
        CFF_9_6: int = 10
        CFF_EPSILON_SIGMA: int = 11
        MATSUOKA_CLEMENTI_YOSHIMINE: int = 12
        GENERIC: int = 13
        PELLENQ_NICHOLSON: int = 14
        HYDRATED_ION_WATER: int = 15
        MIE: int = 16
        WEEKS_CHANDLER_ANDERSEN: int = 17
        LENNARD_JONES_SECOND_ORDER_TAYLOR_SHIFTED: int = 18
        LENNARD_JONES_SWITCHED: int = 19
        LENNARD_JONES_FORCE_SWITCHED: int = 20
        REPULSIVE_HARMONIC: int = 100

    @typing.overload
    def __init__(self, epsilon: float, sigma: float) -> None:
        """
        Initialize a Lennard-Jones potential.

        Args:
            epsilon: Energy-well depth [K].
            sigma: Size/collision-diameter [Angstrom].
        """
        ...

    @typing.overload
    def __init__(self, type: 'VDWParameters.Type', parameters: collections.abc.Sequence[float]) -> None:
        """
        Initialize a potential of the given type.

        Args:
            type: The form of the potential.
            parameters: The parameters of that form, as in a force-field file
                (energies in Kelvin, lengths in Angstrom).
        """
        ...

    @property
    def type(self) -> 'VDWParameters.Type':
        """The form of the pair potential."""
        ...
