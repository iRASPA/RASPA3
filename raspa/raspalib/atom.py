import numpy as np
import raspalib.raspalib as raspalib
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
        position: tuple[float, float, float],
        charge: float,
        scaling: float = 1.0,
        moleculeId: int = 0,
        type: int = 0,
        componentId: int = 0,
        groupId: bool = 0,
        isFractional: bool = 0
    ):
        """
        Initialize the Atom object with provided parameters.

        Args:
            position (tuple[float, float, float]): The position of the atom as a 3-element tuple.
            charge (float): The charge of the atom.
            scaling (float, optional): The lambda value of the atom. Default is 1.0.
            moleculeId (int, optional): The molecule ID of the atom. Default is 0.
            type (int, optional): The type of the atom. Default is 0.
            componentId (int, optional): The component ID of the atom. Default is 0.
            groupId (bool, optional): The group ID of the atom. Default is false.
            isFractional (bool, optional): The isFractional ID of the atom. Default is false.
        """
        super().__init__(**popSelf(locals()))

        self._settings["position"] = raspalib.double3(*self._settings["position"])
        self._cpp_obj = raspalib.Atom(**self.cpp_args())

    @property
    def position(self):
        """
        Get the position of the atom.

        Returns:
            tuple[float, float, float]: The position of the atom as a 3-element tuple.
        """
        self._position = self._cpp_obj.position
        return (self._position.x, self._position.y, self._position.z)

    @position.setter
    def position(self, position: tuple[float, float, float]):
        """
        Set the position of the atom.

        Args:
            position tuple(float, float, float)): The new position of the atom as a 3-element tuple.
        """
        self._position = raspalib.double3(*position)
        self._cpp_obj.position = self._position
