import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class IntraMolecularPotentials(RaspaBase):
    """
    A class representing the intra-molecular potentials of a component in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self
    ):
        """
        Initialize the IntraMolecularPotentials object with provided parameters.

        Args:
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.IntraMolecularPotentials(**self.cpp_args())
