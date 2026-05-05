class PressureData():
    """
    Store pressure contributions reported by RASPA.

    This object provides total pressure and its decomposition into ideal-gas and
    excess components.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize a :class:`PressureData`.
        """

    @property
    def total_pressure(self) -> float:
        """Return the total pressure."""
        ...

    @property
    def excess_pressure(self) -> float:
        """Return the excess (non-ideal) pressure contribution."""
        ...

    @property
    def ideal_gas_pressure(self) -> float:
        """Return the ideal-gas pressure contribution."""
        ...


