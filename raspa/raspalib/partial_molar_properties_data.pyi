import collections.abc

class PartialMolarPropertiesData():
    """
    Store partial molar properties reported by RASPA.

    Holds, per swappable component, the partial molar internal energy and the
    partial molar volume obtained from energy(volume)-particle fluctuations.
    """

    def __init__(
        self,
        number_of_components: int
    ) -> None:
        ...
        """
        Initialize a :class:`PartialMolarPropertiesData`.

        Args:
            number_of_components: Number of (swappable) components.
        """

    @property
    def partial_molar_energy(self) -> collections.abc.Sequence[float]:
        """Return the partial molar internal energy per component (simulation units)."""
        ...

    @property
    def partial_molar_volume(self) -> collections.abc.Sequence[float]:
        """Return the partial molar volume per component (Angstrom^3)."""
        ...
