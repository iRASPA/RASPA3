import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class MoveStatisticsDouble3(RaspaBase):
    """
    A class representing a component in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self
    ):
        """
        Initialize the PropertyLoading object with provided parameters.

        Args:
            numberOfBlocks (int): The number of blocks.
            numberOfComponents (int): The number of components.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.MoveStatisticsDouble3(**self.cpp_args())
