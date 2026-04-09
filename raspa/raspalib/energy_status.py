import raspalib.raspalib as raspalib
from .base import RaspaBase
from .energy_factor import EnergyFactor
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class EnergyStatus(RaspaBase):
    """
    A class representing energy-status-information in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self
    ):
        """
        Initialize the  object with provided parameters.

        Args:
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.EnergyStatus(**self.cpp_args())



    @property
    def totalEnergy(self) -> EnergyFactor:
        """
        Get the total energy of the energy status.

        Returns:
            float: The total energy.
        """
        return self._cpp_obj.totalEnergy
