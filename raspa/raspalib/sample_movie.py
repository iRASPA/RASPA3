import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class SampleMovie(RaspaBase):
    """
    A class representing a component in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        systemId: int,
        sampleEvery: int = 1,
        restrictToBox: bool = True
    ):
        """
        Initialize the SampleMovie  object with provided parameters.

        Args:
            systemId (int): The ID of the system.
            sampleEvery (int, optional): The sample frequency. Default: 1.
            restrictToBox (bool, optional): Whether to resetrict particle to the box confinement. Default: true.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.SampleMovie(**self.cpp_args())
