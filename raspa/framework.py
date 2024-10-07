import raspa.raspalib as raspalib
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
        fileName: str = None,
        simulationBox: SimulationBox = None,
        spaceGroupHallNumber: int = None,
        definedAtoms: list[Atom] = None,
        numberOfUnitCells: list[int] = [1, 1, 1],
    ):
        """
        Initialize the Framework object with provided parameters. There are two methods to initialize the framework.
        Either from .cif file passing it through fileName or from explicit definition of simulation box and atoms.
        Either way a framework needs an id, force field and number of unit cells.

        Args:
            frameworkId (int): The ID of the framework.
            forceField (ForceField): The force field to be used.
            componentName (str): The name of the component.
            fileName (str, optional): The file name for framework initialization. Default is None.
            simulationBox (SimulationBox, optional): The simulation box. Default is None.
            spaceGroupHallNumber (int, optional): The space group Hall number. Default is None.
            definedAtoms (list[Atom], optional): A list of defined atoms. Default is None.
            numberOfUnitCells (list[int], optional): The number of unit cells in each dimension. Default is [1, 1, 1].
        """
        super().__init__(**popSelf(locals()))
        self._settings["numberOfUnitCells"] = raspalib.int3(*self._settings["numberOfUnitCells"])

        # pick either file or manual init
        if self._settings["fileName"] is None:
            self.drop_args("fileName")
        else:
            self.drop_args("simulationBox", "spaceGroupHallNumber", "definedAtoms")

        self._cpp_obj = raspalib.Framework(**self.cpp_args())
