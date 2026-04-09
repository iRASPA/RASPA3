from typing import Literal

import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class PropertyDensityGrid(RaspaBase):
    """
    A class representing a density-grid property in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        numberOfFrameworks: int,
        numberOfComponents: int,
        numberOfGridPoints: tuple[int, int, int] = (128, 128, 128),
        sampleEvery: int  = 1,
        writeEvery: int = 5000,
        densityGridPseudoAtomsList: list[str] = [],
        normalizationType: Literal["Max"] = "Max",
        binningMode: Literal["Standard"] = "Standard"
    ):
        """
        Initialize the PropertyDensityGrid object with provided parameters.

        Args:
            numberOfFrameworks (int): The number of frameworks.
            numberOfComponents (int): The number of components.
            numberOfGridPoints (tuple[int, int, int]): The number of grid points.
            sampleEvery (int): The sample frequency.
            densityGridPseudoAtomsList (list[str]): List of pseudo-atoms.
            normalizationType (): The normalization type
            binningMode (): The binning mode
        """
        super().__init__(**popSelf(locals()))
        self._settings["numberOfGridPoints"] = raspalib.int3(*self._settings["numberOfGridPoints"])
        self._settings["normalizationType"] = getattr(raspalib.PropertyDensityGrid.Normalization, self._settings["normalizationType"])
        self._settings["binningMode"] = getattr(raspalib.PropertyDensityGrid.Standard, self._settings["binningMode"])
        self._cpp_obj = raspalib.PropertyDensityGrid(**self.cpp_args())

        #self._cpp_obj = raspalib.PropertyDensityGrid(numberOfFrameworks, numberOfComponents, raspalib.int3(*numberOfGridPoints), sampleEvery, writeEvery, densityGridPseudoAtomsList, getattr(raspalib.PropertyDensityGrid.Normalization, normalizationType), getattr(raspalib.PropertyDensityGrid.Binning, binningMode))
