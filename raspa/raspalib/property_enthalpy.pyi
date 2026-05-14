from raspalib.enthalpy_of_adsorption_data import *

class PropertyEnthalpy():
    """
    Configure and access enthalpy-related sampling in RASPA.

    This property exposes paired enthalpy-of-adsorption results as
    :class:`EnthalpyOfAdsorptionData` objects.
    """

    def __init__(
        self,
        number_of_blocks: int,
        number_of_components: int
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyEnthalpy`.

        Args:
            number_of_blocks: Number of blocks used for block averaging.
            number_of_components: Number of components in the system.
        """

    def result(self) -> tuple[EnthalpyOfAdsorptionData, EnthalpyOfAdsorptionData]:
        """Return the computed pair of enthalpy-of-adsorption aggregates."""
        ...

