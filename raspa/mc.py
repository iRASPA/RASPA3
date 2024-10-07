import raspa.raspalib as raspalib

from .base import RaspaBase
from .system import System
from .intputreader import InputReader
from .utils import popSelf


class RandomSeed(RaspaBase):
    """
    A class holding the random seed cpp object.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the pseudo atom.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(self, seed=12):
        """
        Initialize random seed.

        Args:
            seed (int, optional): _description_. Defaults to 12.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.random(**self.cpp_args())


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
        systems: list[System] = None,
        numberOfEquilibrationCycles: int = 0,
        printEvery: int = 1000,
        writeBinaryRestartEvery: int = 100,
        rescaleWangLandauEvery: int = 100,
        optimizeMCMovesEvery: int = 100,
        randomSeed: RandomSeed = RandomSeed(),
        numberOfBlocks: int = 5,
        inputReader: InputReader = None,
    ):
        """
        Initializes a Monte Carlo object. There are two different options when it comes to initalization. First, there
        is the initialization via InputReader class. If an InputReader is specified all other parameters are ignored.
        Second, there is the option of initializing using explicit variables (cycles, systems, etc.).

        Args:
            numberOfCycles (int, optional): _description_. Defaults to 0.
            numberOfInitializationCycles (int, optional): _description_. Defaults to 0.
            systems (list[System], optional): _description_. Defaults to None.
            numberOfEquilibrationCycles (int, optional): _description_. Defaults to 0.
            printEvery (int, optional): _description_. Defaults to 1000.
            writeBinaryRestartEvery (int, optional): _description_. Defaults to 100.
            rescaleWangLandauEvery (int, optional): _description_. Defaults to 100.
            optimizeMCMovesEvery (int, optional): _description_. Defaults to 100.
            randomSeed (RandomSeed, optional): _description_. Defaults to RandomSeed().
            numberOfBlocks (int, optional): _description_. Defaults to 5.
            inputReader (InputReader, optional): _description_. Defaults to None.
        """
        super().__init__(**popSelf(locals()))

        if self._settings["inputReader"] is None:
            self.drop_args("inputReader")
        else:
            self.drop_args(
                "numberOfCycles",
                "numberOfInitializationCycles",
                "systems",
                "numberOfEquilibrationCycles",
                "printEvery",
                "writeBinaryRestartEvery",
                "rescaleWangLandauEvery",
                "optimizeMCMovesEvery",
                "randomSeed",
                "numberOfBlocks",
            )

        self._cpp_obj = raspalib.MonteCarlo(**self.cpp_args())

    def run(self):
        """
        Runs the Monte Carlo simulation.
        """
        self._cpp_obj.run()

    def initialize(self):
        """
        Runs just the initialization.
        """
        self._cpp_obj.initialize()

    def equilibrate(self):
        """
        Runs just the equilibration.
        """
        self._cpp_obj.equilibrate()

    def production(self):
        """
        Runs just the production.
        """
        self._cpp_obj.production()
