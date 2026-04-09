import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import popSelf
from .simulationbox import SimulationBox
from .forcefield import ForceField
from .framework import Framework
from .component import Component
from .mc_move_probabilities import MCMoveProbabilities
from .property_energy import PropertyEnergy
from .property_loading import PropertyLoading
from .property_number_of_molecules_evolution import PropertyNumberOfMoleculesEvolution
from .property_volume_evolution import PropertyVolumeEvolution

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
        forceField: ForceField,
        simulationBox: SimulationBox = None,
        hasExternalField: bool = False,
        externalTemperature: float = 300.0,
        externalPressure: float = None,
        heliumVoidFraction: float = 0.29,
        frameworkComponents: Framework = None,
        components: list[Component] = [],
        initialPositions: list[tuple[float, float, float]] = [],
        initialNumberOfMolecules: list[int] = [],
        numberOfBlocks: int = 5,
        systemProbabilities: MCMoveProbabilities = MCMoveProbabilities()
    ):
        """
        Initialize the System object with provided parameters.

        Args:
            systemId (int): The ID of the system.
            forceField (ForceField): The force field to be used.
            simulationBox (SimulationBox, optional): The simulation box. Default is None.
            externalTemperature (float): The temperature of the system.
            externalPressure (float | None): The pressure of the system. Default is None.
            heliumVoidFraction (float): The helium void-fraction of the system.
            frameworkComponents (Framework | None, optional): The framework component if present. Default is None.
            components (list[Component]): A list of components in the system.
            initialPositions (list[int]): A list of initial positions. Default is empty list.
            initialNumberOfMolecules (list[int]): A list of initial number of molecules for each component.
            numberOfBlocks (int, optional): The number of blocks for the simulation. Default is 5.
            systemProbabilities (MCMoveProbabilitiesSystem, optional): The move probabilities system. Default is a new instance of MCMoveProbabilitiesSystem.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.System(**self.cpp_args())

    @property
    def averageEnergies(self) -> PropertyEnergy:
        """
        Get the average-energies property.

        Returns:
            PropertyEnergy: The average-energies property.
        """
        return self._cpp_obj.averageEnergies

    @property
    def averageLoadings(self) -> PropertyLoading:
        """
        Get the avearge-loading property.

        Returns:
            PropertyLoading: The avearge-loading property.
        """
        return self._cpp_obj.averageLoadings

    @property
    def propertyNumberOfMoleculesEvolution(self) -> PropertyNumberOfMoleculesEvolution:
        """
        Get the number of molecule evolution property.

        Returns:
            PropertyNumberOfMoleculesEvolution: The number of molecule evolution property.
        """
        return self._cpp_obj.propertyNumberOfMoleculesEvolution

    @property
    def propertyVolumeEvolution(self) -> PropertyVolumeEvolution:
        """
        Get the volume evolution property.

        Returns:
            PropertyVolumeEvolution: The volume evolution property.
        """
        return self._cpp_obj.propertyVolumeEvolution
