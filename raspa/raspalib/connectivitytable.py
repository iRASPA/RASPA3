import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RaspaError, SHARE_DIR, popSelf

import os

class ConnectivityTable(RaspaBase):
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
        Initialize the ConnectivityTable object.

        Args:
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.ConnectivityTable(**self.cpp_args())
