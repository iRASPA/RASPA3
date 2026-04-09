import raspalib.raspalib as raspalib
from .base import RaspaBase
from .energy_status import EnergyStatus
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class PropertyEnergy(RaspaBase):
    """
    A class representing the average energy property in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        numberOfBlocks: int,
        numberOfExternalFields: int,
        numberOfFrameworks: int,
        numberOfComponents: int
    ):
        """
        Initialize the PropertyEnergy object with provided parameters.

        Args:
            numberOfBlocks (int): The number of blocks.
            numberOfExternalFields (int): The number of external fields.
            numberOfFrameworks (int): The number of frameworks.
            numberOfComponents (int): The number of components.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.EnergyStatus(**self.cpp_args())


    def result(self) -> tuple[EnergyStatus, EnergyStatus]:
        """
        Returns the computed data.
        """
        self._cpp_obj.result()
