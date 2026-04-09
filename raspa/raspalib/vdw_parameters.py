from typing import Literal

import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RASPA_DIR, SHARE_DIR, popSelf
import os
import json

class VDWParameters(RaspaBase):
    """
    A class representing Van der Waals parameters in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the Van der Waals parameters.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(self, epsilon: float, sigma: float):
        """
        Initialize the VDWParameter object with provided parameters.

        Args:
            epsilon (float): The epsilon parameter.
            sigma (float): The sigma parameter.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.VDWParameters(**self.cpp_args())
