import raspalib.raspalib as raspalib
from .base import RaspaBase
from .loadings import Loadings
from .average_energy_type import AverageEnergyType
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class PropertyLambdaProbabilityHistogram(RaspaBase):
    """
    A class representing a lambda-histogram property in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        numberOfBlocks: int,
        numberOfSamplePoints: int
    ):
        """
        Initialize the PropertyEnergyHistogram object with provided parameters. 

        Args:
            numberOfBlocks (int): The number of blocks.
            numberOfSamplePoints (int): The number of sample points.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.PropertyLambdaProbabilityHistogram(**self.cpp_args())
