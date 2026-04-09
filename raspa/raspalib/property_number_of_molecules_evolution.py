import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RaspaError, SHARE_DIR, popSelf

import os

class PropertyNumberOfMoleculesEvolution(RaspaBase):
    """
    A class representing a component in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        numberOfCycles: int,
        numberOfComponents: int,
        sampleEvery: int,
        writeEvery: int = None
    ):
        """
        Initialize the PropertyNumberOfMoleculesEvolution object with provided parameters.

        Args:
            numberOfCycles (int): The number of cycles.
            numberOfComponents (int): The number of components.
            sampleEvery (int): The sample frequency.
            writeEvery (int): The write frequency.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.PropertyNumberOfMoleculesEvolution(**self.cpp_args())


    def result(self) -> list[list[int]]:
        """
        Runs just the equilibration.
        """
        self._cpp_obj.result()
