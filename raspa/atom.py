import numpy as np
import raspa.raspalib as raspalib
from .base import RaspaBase
from .utils import popSelf


class Atom(RaspaBase):
    """
    A class representing an atom in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the atom.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        position: np.ndarray,
        charge: float,
        lambda_: float = 1.0,
        moleculeId: int = 0,
        type: int = 0,
        componentId: int = 0,
        groupId: int = 0,
    ):
        """
        Initialize the Atom object with provided parameters.

        Args:
            position (np.ndarray): The position of the atom as a 3-element numpy array.
            charge (float): The charge of the atom.
            lambda_ (float, optional): The lambda value of the atom. Default is 1.0.
            moleculeId (int, optional): The molecule ID of the atom. Default is 0.
            type (int, optional): The type of the atom. Default is 0.
            componentId (int, optional): The component ID of the atom. Default is 0.
            groupId (int, optional): The group ID of the atom. Default is 0.
        """
        super().__init__(**popSelf(locals()))

        # special word "lambda" will crash in python
        self._settings["lambda"] = self._settings.pop("lambda_")
        self._cpp_obj = raspalib.Atom(**self.cpp_args())

    @property
    def position(self):
        """
        Get the position of the atom.

        Returns:
            np.ndarray: The position of the atom as a 3-element numpy array.
        """
        self._position = self._cpp_obj.position
        return np.array([self._position.x, self._position.y, self._position.z])

    @position.setter
    def position(self, position: np.ndarray):
        """
        Set the position of the atom.

        Args:
            position (np.ndarray): The new position of the atom as a 3-element numpy array.
        """
        assert position.shape == (3,)
        self._position = raspalib.double3(*position)
        self._cpp_obj.position = self._position
