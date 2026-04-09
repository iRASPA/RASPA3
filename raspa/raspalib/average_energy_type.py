import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class AverageEnergyType(RaspaBase):
    """
    A class representing the various to the energies in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        totalEnergy: float,
        VanDerWaalsEnergy: float,
        CoulombEnergy: float,
        polarizationEnergy: float
    ):
        """
        Initialize object with provided parameters.

        Args:
            totalEnergy (float): The total energy.
            VanDerWaalsEnergy (float): The Van der Waals energy.
            CoulombEnergy (float): The Coulomb energy.
            polarizationEnergy (float): The polarization energy.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.AverageEnergyType(**self.cpp_args())
