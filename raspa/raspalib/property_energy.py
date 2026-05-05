from raspalib.energy_status import *

class PropertyEnergy():
    """
    Configure and access average-energy sampling in RASPA.

    This property accumulates energy statistics and exposes paired
    :class:`EnergyStatus` results (e.g., mean values and associated uncertainty
    or block-averaged companion data).
    """

    def __init__(
        self,
        number_of_blocks: int,
        number_of_external_fields: int,
        number_of_frameworks: int,
        number_of_components: int
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyEnergy`.

        Args:
            number_of_blocks: Number of blocks used for block averaging.
            number_of_external_fields: Number of external fields in the system.
            number_of_frameworks: Number of frameworks in the system.
            number_of_components: Number of components in the system.
        """

    def result(self) -> tuple[EnergyStatus, EnergyStatus]:
        """Return the computed pair of energy-status aggregates."""
        ...

