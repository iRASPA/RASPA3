import raspalib.raspalib as raspalib
from .base import RaspaBase
from .loadings import Loadings
from .average_energy_type import AverageEnergyType
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class PropertyEnergyHistogram(RaspaBase):
    """
    A class representing an energy-histogram property in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        numberOfBlocks: int,
        numberOfBins: int,
        valueRange: tuple[float, float],
        sampleEvery: int = 1,
        writeEvery: int = 5000
    ):
        """
        Initialize the PropertyEnergyHistogram object with provided parameters. 

        Args:
            numberOfBlocks (int): The number of blocks.
            numberOfBins (int): The number of bins.
            valueRange (int): The range of energy values to consider.
            sampleEvery (int): The sample frequency.
            writeEvery (int): The write frequency.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.PropertyEnergyHistogram(**self.cpp_args())

    def result(self) -> tuple[list[float], list[AverageEnergyType], list[AverageEnergyType]]:
        """
        Returns the computed data.
        """
        self._cpp_obj.result()
