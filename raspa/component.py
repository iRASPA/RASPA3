import raspa.raspalib as raspalib
from .base import RaspaBase
from .forcefield import ForceField
from .atom import Atom
from .mcmoveprobabilities import MCMoveProbabilitiesParticles
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
        fileName: str = None,
        criticalTemperature: float = None,
        criticalPressure: float = None,
        acentricFactor: float = None,
        definedAtoms: list[Atom] = None,
        numberOfBlocks: int = 5,
        numberOfLambdaBins: int = 21,
        particleProbabilities: MCMoveProbabilitiesParticles = MCMoveProbabilitiesParticles(),
        fugacityCoefficient: float = None,
        thermodynamicIntegration: bool = False,
    ):
        """
        Initialize the Component object with provided parameters. There are two types of initialization. Both
        initializations require at least an id, a force field and a move probabilities objects. The first init requires
        explicit specification of the object using thermodynamic constants (criticalTemperature, criticalPressure,
        acentricFactor) and definedAtoms. The second initializes from file (example CO2 in data/molecules/co2.json).

        Args:
            componentId (int): The ID of the component.
            forceField (ForceField): The force field to be used.
            componentName (str): The name of the component.
            fileName (str, optional): The file name for component initialization. Default is None.
            criticalTemperature (float, optional): The critical temperature. Default is None.
            criticalPressure (float, optional): The critical pressure. Default is None.
            acentricFactor (float, optional): The acentric factor. Default is None.
            definedAtoms (list[Atom], optional): A list of defined atoms. Default is None.
            numberOfBlocks (int, optional): The number of blocks for the simulation. Default is 5.
            numberOfLambdaBins (int, optional): The number of lambda bins. Default is 21.
            particleProbabilities (MCMoveProbabilitiesParticles, optional): The particle move probabilities. Default is a new instance of MCMoveProbabilitiesParticles.
            fugacityCoefficient (float, optional): The fugacity coefficient. Default is None.
            thermodynamicIntegration (bool, optional): Whether to use thermodynamic integration. Default is False.
        """
        super().__init__(**popSelf(locals()))

        if self._settings["fileName"] is None:
            self.drop_args("fileName")
        else:
            self._settings["type"] = raspalib.Component.Adsorbate
            self.drop_args("criticalTemperature", "criticalPressure", "acentricFactor", "definedAtoms")

        self._cpp_obj = raspalib.Component(**self.cpp_args())

    @classmethod
    def exampleCO2(
        cls,
        componentId: int,
        forceField: ForceField,
        particleProbabilities: MCMoveProbabilitiesParticles = MCMoveProbabilitiesParticles(),
        numberOfBlocks: int = 5,
        numberOfLambdaBins: int = 21,
        fugacityCoefficient: float = None,
        thermodynamicIntegration: bool = False,
    ):
        """
        Create an example CO2 component.

        Args:
            componentId (int): The ID of the component.
            forceField (ForceField): The force field to be used.
            particleProbabilities (MCMoveProbabilitiesParticles, optional): The particle move probabilities. Default is a new instance of MCMoveProbabilitiesParticles.
            numberOfBlocks (int, optional): The number of blocks for the simulation. Default is 5.
            numberOfLambdaBins (int, optional): The number of lambda bins. Default is 21.
            fugacityCoefficient (float, optional): The fugacity coefficient. Default is None.
            thermodynamicIntegration (bool, optional): Whether to use thermodynamic integration. Default is False.

        Returns:
            Component: An instance of Component representing CO2.
        """
        settings = locals()
        settings.pop("cls")
        return cls(
            componentName="CO2",
            fileName=os.path.join(SHARE_DIR, "molecules", "example_definitions", "co2.json"),
            **settings,
        )

    @classmethod
    def exampleCH4(
        cls,
        componentId: int,
        forceField: ForceField,
        particleProbabilities: MCMoveProbabilitiesParticles = MCMoveProbabilitiesParticles(),
        numberOfBlocks: int = 5,
        numberOfLambdaBins: int = 21,
        fugacityCoefficient: float = None,
        thermodynamicIntegration: bool = False,
    ):
        """
        Create an example CH4 component.

        Args:
            componentId (int): The ID of the component.
            forceField (ForceField): The force field to be used.
            particleProbabilities (MCMoveProbabilitiesParticles, optional): The particle move probabilities. Default is a new instance of MCMoveProbabilitiesParticles.
            numberOfBlocks (int, optional): The number of blocks for the simulation. Default is 5.
            numberOfLambdaBins (int, optional): The number of lambda bins. Default is 21.
            fugacityCoefficient (float, optional): The fugacity coefficient. Default is None.
            thermodynamicIntegration (bool, optional): Whether to use thermodynamic integration. Default is False.

        Returns:
            Component: An instance of Component representing CH4.
        """
        settings = locals()
        settings.pop("cls")
        return cls(
            componentName="CH4",
            fileName=os.path.join(SHARE_DIR, "molecules", "example_definitions", "methane.json"),
            **settings,
        )

    @classmethod
    def exampleN2(
        cls,
        componentId: int,
        forceField: ForceField,
        particleProbabilities: MCMoveProbabilitiesParticles = MCMoveProbabilitiesParticles(),
        numberOfBlocks: int = 5,
        numberOfLambdaBins: int = 21,
        fugacityCoefficient: float = None,
        thermodynamicIntegration: bool = False,
    ):
        """
        Create an example N2 component.

        Args:
            componentId (int): The ID of the component.
            forceField (ForceField): The force field to be used.
            particleProbabilities (MCMoveProbabilitiesParticles, optional): The particle move probabilities. Default is a new instance of MCMoveProbabilitiesParticles.
            numberOfBlocks (int, optional): The number of blocks for the simulation. Default is 5.
            numberOfLambdaBins (int, optional): The number of lambda bins. Default is 21.
            fugacityCoefficient (float, optional): The fugacity coefficient. Default is None.
            thermodynamicIntegration (bool, optional): Whether to use thermodynamic integration. Default is False.

        Returns:
            Component: An instance of Component representing N2.
        """
        settings = locals()
        settings.pop("cls")
        return cls(
            componentName="N2",
            fileName=os.path.join(SHARE_DIR, "molecules", "example_definitions", "n2.json"),
            **settings,
        )
