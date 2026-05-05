class WidomData:
    """
    Store Widom-insertion result components reported by RASPA.

    This container exposes total, excess, and ideal-gas contributions for
    quantities derived from Widom sampling.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize a :class:`WidomData`.
        """
    @property
    def total(self) -> float:
        """Return the total Widom-derived quantity."""
        ...
    @property
    def excess(self) -> float:
        """Return the excess (non-ideal) Widom contribution."""
        ...
    @property
    def ideal_gas(self) -> float:
        """Return the ideal-gas Widom contribution."""
        ...
