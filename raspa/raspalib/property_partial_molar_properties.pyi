from raspalib.partial_molar_properties_data import *

class PropertyPartialMolarProperties():
    """
    Configure and access partial-molar-property sampling in RASPA.

    Computes the partial molar internal energy and partial molar volume of each
    swappable component from grand-canonical/osmotic fluctuations, using the same
    particle-number covariance matrix as the enthalpy of adsorption.
    """

    def __init__(
        self,
        number_of_blocks: int,
        number_of_components: int
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyPartialMolarProperties`.

        Args:
            number_of_blocks: Number of blocks used for block averaging.
            number_of_components: Number of components in the system.
        """

    def average_properties(self) -> tuple[PartialMolarPropertiesData, PartialMolarPropertiesData]:
        """Return the (mean, 95% confidence interval) partial molar properties."""
        ...
