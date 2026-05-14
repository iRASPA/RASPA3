from raspalib.widom_data import *

class PropertyWidom:
    """
    Configure and access Widom-insertion sampling results in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize a :class:`PropertyWidom`.
        """

    def chemical_potential_result(self, temperature: float) -> tuple[WidomData, WidomData]:
        """Return Widom chemical-potential aggregates at a given temperature.

        Args:
            temperature: Temperature used to convert Widom statistics.

        Returns:
            Pair of :class:`WidomData` aggregates for chemical potential.
        """
        ...

    def fugacity_result(self, temperature: float) -> float:
        """Return fugacity derived from Widom statistics at a temperature.

        Args:
            temperature: Temperature used to convert Widom statistics.

        Returns:
            Fugacity in Pascal.
        """
        ...


