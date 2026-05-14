from typing import overload
import enum

class SimulationBox():
    """
    Represent simulation-cell geometry used by RASPA.

    The box can be rectangular (orthorhombic) or fully triclinic, depending on
    whether only edge lengths or both lengths and angles are provided.
    """

    class SimulationBoxType(enum.IntEnum):
        """Supported simulation-box geometry types."""

        Rectangular: int = 0
        Triclinic: int = 1

    @overload
    def __init__(self, a: float, b: float, c: float) -> None:
        ...
        """
        Initialize a rectangular simulation box from edge lengths.

        Args:
            a: ``a`` edge length of the unit cell.
            b: ``b`` edge length of the unit cell.
            c: ``c`` edge length of the unit cell.
        """

    @overload
    def __init__(self, a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> None:
        ...
        """
        Initialize a triclinic simulation box from lengths and angles.

        Args:
            a: ``a`` edge length of the unit cell.
            b: ``b`` edge length of the unit cell.
            c: ``c`` edge length of the unit cell.
            alpha: Angle (degrees) between ``b`` and ``c``.
            beta: Angle (degrees) between ``a`` and ``c``.
            gamma: Angle (degrees) between ``a`` and ``b``.
        """

    @property
    def volume(self) -> float:
        """Return the simulation-box volume."""
        ...

