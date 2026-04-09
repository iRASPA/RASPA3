import raspalib.raspalib as raspalib

from .base import RaspaBase
from .system import System
from .utils import popSelf
from .system import System


class MonteCarlo(RaspaBase):
    """

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the Van der Waals parameters.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        numberOfCycles: int = 0,
        numberOfInitializationCycles: int = 0,
        numberOfEquilibrationCycles: int = 0,
        printEvery: int = 1000,
        writeBinaryRestartEvery: int = 5000,
        rescaleWangLandauEvery: int = 5000,
        optimizeMCMovesEvery: int = 5000,
        systems: list[System] = [],
        randomSeed: int = None,
        numberOfBlocks: int = 5,
        outputToFiles: bool = False
    ):
        """
        Initializes a Monte Carlo object.

        Args:
            numberOfCycles (int, optional): _description_. Defaults to 0.
            numberOfInitializationCycles (int, optional): _description_. Defaults to 0.
            numberOfEquilibrationCycles (int, optional): _description_. Defaults to 0.
            numberOfEquilibrationCycles (int, optional): _description_. Defaults to 0.
            printEvery (int, optional): _description_. Defaults to 1000.
            writeBinaryRestartEvery (int, optional): _description_. Defaults to 100.
            rescaleWangLandauEvery (int, optional): _description_. Defaults to 100.
            optimizeMCMovesEvery (int, optional): _description_. Defaults to 100.
            systems (list[System], optional): _description_. Defaults to None.
            randomSeed (int, optional): _description_. Defaults to a random integer.
            numberOfBlocks (int, optional): _description_. Defaults to 5.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.MonteCarlo(**self.cpp_args())

    @property
    def systems(self) -> list[System]:
        """
        Get the positions of the atoms in the system.

        Returns:
            np.ndarray: An array of atom positions.
        """
        return self._cpp_obj.systems

    def run(self):
        """
        Runs the Monte Carlo simulation.
        """
        self._cpp_obj.run()

    def initialize(self,
                   call_back_function: callable = None):
        """
        Runs just the initialization.
        """
        self._cpp_obj.initialize(call_back_function)

    def equilibrate(self, 
                    call_back_function: callable = None):
        """
        Runs just the equilibration.
        """
        self._cpp_obj.equilibrate(call_back_function)

    def production(self,
                   call_back_function: callable = None):
        """
        Runs just the production.
        """
        self._cpp_obj.production(call_back_function)
