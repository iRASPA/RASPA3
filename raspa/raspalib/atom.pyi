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
        charge: float,
        scaling: float = 1.0,
        molecule_id: int = 0,
        type: int = 0,
        component_id: int = 0,
        group_id: bool = 0,
        is_fractional: bool = 0
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

