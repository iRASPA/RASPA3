import raspa.raspalib as raspalib
from .base import RaspaBase
from .utils import popSelf
from .simulationbox import SimulationBox
from .forcefield import ForceField
from .framework import Framework
from .component import Component
from .mcmoveprobabilities import MCMoveProbabilitiesSystem

import numpy as np


class System(RaspaBase):
    """
    A class representing a system in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the system.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        systemId: int,
        temperature: float,
        forceField: ForceField,
        components: list[Component],
        initialNumberOfMolecules: list[int],
        numberOfBlocks: int = 5,
        pressure: float = None,
        heliumVoidFraction: float = 0.29,
        frameworkComponents: list[Framework] = [],
        simulationBox: SimulationBox = None,
        systemProbabilities: MCMoveProbabilitiesSystem = MCMoveProbabilitiesSystem(),
        sampleMoviesEvery: int = None,
    ):
        """
        Initialize the System object with provided parameters.

        Args:
            systemId (int): The ID of the system.
            temperature (float): The temperature of the system.
            forceField (ForceField): The force field to be used.
            components (list[Component]): A list of components in the system.
            initialNumberOfMolecules (list[int]): A list of initial number of molecules for each component.
            numberOfBlocks (int, optional): The number of blocks for the simulation. Default is 5.
            pressure (float, optional): The pressure of the system. Default is None.
            frameworkComponents (list[Framework], optional): A list of framework components. Default is an empty list.
            simulationBox (SimulationBox, optional): The simulation box. Default is None.
            systemProbabilities (MCMoveProbabilitiesSystem, optional): The move probabilities system. Default is a new instance of MCMoveProbabilitiesSystem.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.System(**self.cpp_args())

    @property
    def atomPositions(self):
        """
        Get the positions of the atoms in the system.

        Returns:
            np.ndarray: An array of atom positions.
        """
        return np.array([[atom.position.x, atom.position.y, atom.position.z] for atom in self._cpp_obj.atomPositions])

    @atomPositions.setter
    def atomPositions(self, index_position_tuple: tuple[np.ndarray, np.ndarray]):
        """
        Set the positions of the atoms in the system.

        Args:
            index_position_tuple (tuple[np.ndarray, np.ndarray]): A tuple containing an array of indices and an array of positions.
        """
        indices, position = index_position_tuple
        for i, idx in enumerate(indices):
            setattr(self.atomPositions[idx], "position", raspalib.double3(*position[i]))

    def computeTotalEnergies(self):
        """
        Compute the total energies of the system.

        Returns:
            The total energies of the system.
        """
        return self._cpp_obj.computeTotalEnergies()
