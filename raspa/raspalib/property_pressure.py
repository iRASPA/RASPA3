from raspalib.pressure_data import *

class PropertyPressure():
    """
    Configure and access pressure-property sampling in RASPA.

    This property exposes paired :class:`PressureData` aggregates for pressure
    observables.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize a :class:`PropertyPressure`.
        """

    def result(self) -> tuple[PressureData, PressureData]:
        """Return the computed pair of pressure-data aggregates."""
        ...


