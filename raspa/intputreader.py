import raspa.raspalib as raspalib
from .base import RaspaBase
from .utils import popSelf


class InputReader(RaspaBase):
    """
    A class loading python simulation objects from conventional input file.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the RASPA object.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(self, fileName: str = "simulation.json"):
        """
        Initializes the InputReader object from file.

        Args:
            fileName (str, optional): _description_. Defaults to "simulation.json".
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.InputReader(**self.cpp_args())
