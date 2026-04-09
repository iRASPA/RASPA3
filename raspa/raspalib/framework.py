import raspalib.raspalib as raspalib
from .base import RaspaBase
from .forcefield import ForceField
from .simulationbox import SimulationBox
from .atom import Atom
from .utils import popSelf


class Framework(RaspaBase):
    """
    A class representing a framework in RASPA, managing the initialization and configuration
    of a simulation framework.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the framework.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        frameworkId: int,
        forceField: ForceField,
        componentName: str,
        simulationBox: SimulationBox,
        spaceGroupHallNumber: int = None,
        definedAtoms: list[Atom] = None,
        numberOfUnitCells: list[int] = [1, 1, 1],
    ):
        """
        Initialize the Framework object with provided parameters.

        Args:
            frameworkId (int): The ID of the framework.
            forceField (ForceField): The force field to be used.
            componentName (str): The name of the component.
            simulationBox (SimulationBox): The unit cell box. 
            spaceGroupHallNumber (int, optional): The space group Hall number. Default is None.
            definedAtoms (list[Atom], optional): A list of defined atoms. Default is None.
            numberOfUnitCells (list[int], optional): The number of unit cells in each dimension. Default is [1, 1, 1].
        """
        super().__init__(**popSelf(locals()))
        self._settings["numberOfUnitCells"] = raspalib.int3(*self._settings["numberOfUnitCells"])
        self._cpp_obj = raspalib.Framework(**self.cpp_args())
