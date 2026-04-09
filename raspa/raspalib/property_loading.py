import raspalib.raspalib as raspalib
from .base import RaspaBase
from .loadings import Loadings
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class PropertyLoading(RaspaBase):
    """
    A class representing a component in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        numberOfBlocks: int,
        numberOfComponents: int
    ):
        """
        Initialize the PropertyLoading object with provided parameters.

        Args:
            numberOfBlocks (int): The number of blocks.
            numberOfComponents (int): The number of components.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.PropertyLoading(**self.cpp_args())

    def result(self) -> tuple[Loadings, Loadings]:
        """
        Returns the computed loadings and errors.
        """
        self._cpp_obj.result()
