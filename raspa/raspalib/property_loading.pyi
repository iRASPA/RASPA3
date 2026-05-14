from raspalib.loading_data import *

class PropertyLoading():
    """
    Configure and access loading-property sampling in RASPA.

    This property tracks component-wise loading metrics and exposes aggregate
    :class:`LoadingData` results.
    """

    def __init__(
        self,
        number_of_blocks: int,
        number_of_components: int
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyLoading`.

        Args:
            number_of_blocks: Number of blocks used for block averaging.
            number_of_components: Number of components in the system.
        """

    def average_loading_number_of_molecules(self, arg0: int) -> tuple[float, float]:
        """Return average molecule loading and companion statistic for a component."""
        ...
    def result(self) -> tuple[LoadingData, LoadingData]:
        """Return the computed pair of loading-data aggregates."""
        ...

