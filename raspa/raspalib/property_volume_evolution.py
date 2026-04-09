import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class PropertyVolumeEvolution(RaspaBase):
    """
    A class representing a a volume-evoluion property in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        numberOfCycles: int,
        sampleEvery: int,
        writeEvery: int = None
    ):
        """
        Initialize the PropertyVolumeEvolution object with provided parameters.

        Args:
            numberOfCycles (int): The number of cycles.
            sampleEvery (int): The sample frequency.
            writeEvery (int): The write frequency.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.PropertyVolumeEvolution(**self.cpp_args())

    def result(self) -> list[list[float]]:
        """
        Returns the computed values.
        """
        self._cpp_obj.result()
