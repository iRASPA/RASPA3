import raspalib.raspalib as raspalib
import numpy as np
from .base import RaspaBase
from .forcefield import ForceField
from .atom import Atom
from .connectivitytable import ConnectivityTable
from .intra_molecular_potentials import IntraMolecularPotentials
from .mc_move_probabilities import MCMoveProbabilities
from .utils import RaspaError, SHARE_DIR, popSelf

import os


class Component(RaspaBase):
    """
    A class representing a component in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the component.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        componentId: int,
        forceField: ForceField,
        componentName: str,
        criticalTemperature: float,
        criticalPressure: float,
        acentricFactor: float,
        definedAtoms: list[Atom] = [],
        connectivityTable: ConnectivityTable = ConnectivityTable(),
        intraMolecularPotentials: IntraMolecularPotentials = IntraMolecularPotentials(),
        numberOfBlocks: int = 5,
        numberOfLambdaBins: int = 21,
        particleProbabilities: MCMoveProbabilities = MCMoveProbabilities(),
        fugacityCoefficient: float = None,
        thermodynamicIntegration: bool = False,
        blockingPockets: list[tuple[float, float, float, float]] = []
    ):
        """
        Initialize the Component object with provided parameters. 

        Args:
            componentId (int): The ID of the component.
            forceField (ForceField): The force field to be used.
            componentName (str): The name of the component.
            criticalTemperature (float, optional): The critical temperature. Default is None.
            criticalPressure (float, optional): The critical pressure. Default is None.
            acentricFactor (float, optional): The acentric factor. Default is None.
            definedAtoms (list[Atom], optional): A list of defined atoms. Default is None.
            connectivityTable (ConnectivityTable, optional): The connectivity table of the atoms. Default is empty.
            intraMolecularPotentials (IntraMolecularPotentials, optional): The intra molecular potentials. Default is none.
            numberOfBlocks (int, optional): The number of blocks for the simulation. Default is 5.
            numberOfLambdaBins (int, optional): The number of lambda bins. Default is 21.
            particleProbabilities (MCMoveProbabilitiesParticles, optional): The particle move probabilities. Default is a new instance of MCMoveProbabilitiesParticles.
            fugacityCoefficient (float, optional): The fugacity coefficient. Default is None.
            thermodynamicIntegration (bool, optional): Whether to use thermodynamic integration. Default is False.
            blockingPockets (list[tuple[float, float, float, float]]): List of blocking-pockets
        """
        super().__init__(**popSelf(locals()))
        self._settings["blockingPockets"] = [raspalib.double4(*p) for p in self._settings["blockingPockets"]]
        self._cpp_obj = raspalib.Component(**self.cpp_args())

    @property
    def blockingPockets(self):
        """
        Get the blocking-pockets

        Returns:
            (list[tuple[float, float, float, float]]): The positions and radius of the pockets as a list of 4-element tuples.
        """
        self._blockingPockets = self._cpp_obj.blockingPockets
        return [(p.x, p.y, p.z, p.w) for p in self._blockingPockets]

    @blockingPockets.setter
    def blockingPockets(self, positions: list[tuple[float, float, float, float]]):
        """
        Set the blocking-pockets

        Args:
            positions (list[tuple[float, float, float, float]]): The blocking-pockets a list of 4-element tuples.
        """
        self._cpp_obj.blockingPockets = [raspalib.double4(*p) for p in positions]

