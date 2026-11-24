import collections.abc
import typing


class ConnectivityTable():
    """
    Store molecular connectivity information used by RASPA.

    A connectivity table describes which atoms are connected and provides the
    topological basis for bonded interactions (e.g., bonds, angles, and
    torsions) in a component. A component makes its table from the bonds of its
    intra-molecular potentials, so a table is seldom given.
    """

    @typing.overload
    def __init__(self) -> None:
        """Initialize an empty :class:`ConnectivityTable`."""
        ...

    @typing.overload
    def __init__(self, number_of_beads: int,
                 bonds: collections.abc.Sequence[collections.abc.Sequence[int]] = []) -> None:
        """
        Initialize a :class:`ConnectivityTable`.

        Args:
            number_of_beads: The number of atoms of the molecule.
            bonds: The pairs of connected atoms.
        """
        ...

    @property
    def number_of_beads(self) -> int:
        """The number of atoms of the molecule."""
        ...

    @property
    def bonds(self) -> collections.abc.Sequence[collections.abc.Sequence[int]]:
        """The pairs of connected beads."""
        ...

    @property
    def bends(self) -> collections.abc.Sequence[collections.abc.Sequence[int]]:
        """The triplets of beads that form a bend."""
        ...

    @property
    def torsions(self) -> collections.abc.Sequence[collections.abc.Sequence[int]]:
        """The quadruplets of beads that form a torsion."""
        ...
