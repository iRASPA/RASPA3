class Atom():
    """
    Represent a single atom in a simulation system.

    The class stores geometric, electrostatic, and bookkeeping metadata used by
    RASPA to identify atoms in molecules/components and to determine how an atom
    should be interpreted in coordinate space.
    """

    def __init__(
        self,
        position: tuple[float, float, float],
        charge: float = 0.0,
        scaling: float = 1.0,
        molecule_id: int = 0,
        type: int = 0,
        component_id: int = 0,
        group_id: bool = False,
        is_fractional: bool = False
    ) -> None:
        ...
        """
        Initialize an :class:`Atom`.

        Args:
            position: Cartesian coordinates ``(x, y, z)`` of the atom.
            charge: Partial charge assigned to the atom.
            scaling: Coupling/scaling factor (lambda) for this atom. Defaults
                to ``1.0``.
            molecule_id: Identifier of the parent molecule. Defaults to ``0``.
            type: Integer atom-type identifier used by force-field logic.
                Defaults to ``0``.
            component_id: Identifier of the parent component. Defaults to ``0``.
            group_id: Identifier/flag for group membership. Defaults to ``0``.
            is_fractional: Whether ``position`` is provided in fractional
                coordinates instead of Cartesian coordinates. Defaults to ``0``
                (``False``).
        """

    @property
    def charge(self) -> float:
        """The charge of the atom [e]."""
        ...

    @charge.setter
    def charge(self, value: float) -> None:
        ...

    @property
    def type(self) -> int:
        """The index of the pseudo-atom of the atom."""
        ...

    @type.setter
    def type(self, value: int) -> None:
        ...

    @property
    def molecule_id(self) -> int: ...

    @property
    def component_id(self) -> int: ...

    @property
    def scaling_vdw(self) -> float:
        """The scaling of the Van der Waals interactions of the atom."""
        ...

    @property
    def scaling_coulomb(self) -> float:
        """The scaling of the Coulomb interactions of the atom."""
        ...

    @property
    def is_fractional(self) -> bool:
        """Whether the atom is of a fractional molecule."""
        ...
