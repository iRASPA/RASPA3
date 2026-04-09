import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class EnergyFactor(RaspaBase):
    """
    A class representing energy-factor-information in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        energy: float,
        dUdlambda: float
    ):
        """
        Initialize the  object with provided parameters.

        Args:
            energy (float): The energy.
            dUdlambda (float): The dUdlambda factor.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.EnergyStatus(**self.cpp_args())


    @property
    def energy(self) -> float:
        """
        Get the energy of the energy factor.

        Returns:
            float: The energy.
        """
        return self._cpp_obj.energy

    @property
    def dUdlambda(self) -> float:
        """
        Get the dUdlambda of the energy factor.

        Returns:
            float: The dUdlambda.
        """
        return self._cpp_obj.dUdlambda
