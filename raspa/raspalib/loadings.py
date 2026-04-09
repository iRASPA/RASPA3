import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class Loadings(RaspaBase):
    """
    A class representing loading-information in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        numberOfMolecules: list[int],
        numberDensities: list[float],
        inverseNumberDensities: list[float]
    ):
        """
        Initialize the Loadings object with provided parameters.

        Args:
            numberOfMolecules (list[int]): The number of molecules for each component
            numberDensities (list[float]): The number density for each component
            inverseNumberDensities (list[float]): The inverse number-density for each component
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.Loadings(**self.cpp_args())
