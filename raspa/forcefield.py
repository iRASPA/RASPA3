from typing import Literal

import raspa.raspalib as raspalib
from .base import RaspaBase
from .utils import RASPA_DIR, SHARE_DIR, popSelf
import os
import json


class PseudoAtom(RaspaBase):
    """
    A class representing a pseudo atom in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the pseudo atom.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        name: str = "C",
        mass: float = 1.0,
        charge: float = 0.0,
        polarizability: float = 0.0,
        atomicNumber: int = 8,
        printToPDB: bool = False,
        source: str = "-",
    ):
        """
        Initialize the PseudoAtom object with provided parameters.

        Args:
            name (str): The name of the pseudo atom. Default is "C".
            mass (float): The mass of the pseudo atom. Default is 1.0.
            charge (float): The charge of the pseudo atom. Default is 0.0.
            polarizability (float): The polarizability of the pseudo atom. Default is 0.0.
            atomicNumber (int): The atomic number of the pseudo atom. Default is 8.
            printToPDB (bool): Whether to print to PDB. Default is False.
            source (str): The source of the pseudo atom. Default is "-".
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.PseudoAtom(**self.cpp_args())


class VDWParameter(RaspaBase):
    """
    A class representing Van der Waals parameters in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the Van der Waals parameters.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(self, epsilon: float, sigma: float):
        """
        Initialize the VDWParameter object with provided parameters.

        Args:
            epsilon (float): The epsilon parameter.
            sigma (float): The sigma parameter.
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.VDWParameters(**self.cpp_args())


class ForceField(RaspaBase):
    """
    A class representing a force field in RASPA.

    Inherits from RaspaBase.

    Attributes:
        _settings (dict): A dictionary storing the settings for the force field.
        _cpp_obj: A reference to the associated C++ object.
    """

    def __init__(
        self,
        fileName: str = None,
        pseudoAtoms: list[PseudoAtom] = None,
        parameters: list[VDWParameter] = None,
        mixingRule: Literal["Lorentz_Berthelot"] = "Lorentz_Berthelot",
        cutOff: float = 12.0,
        shifted: bool = False,
        tailCorrections: bool = False,
        useCharge=True,
    ):
        """
        Initialize the ForceField object with provided parameters. There are two methods of initializing an object.
        First, there is the initialization via a force field json file. When a file is used in the constructor all
        other parameters are ignored. The other method is via explicit definition of pseudoAtoms, parameters and mixing
        rule.

        Args:
            fileName (str, optional): The file name for force field initialization. Default is None.
            pseudoAtoms (list[PseudoAtom], optional): A list of pseudo atoms. Default is None.
            parameters (list[VDWParameter], optional): A list of Van der Waals parameters. Default is None.
            mixingRule (Literal["Lorentz_Berthelot"], optional): The mixing rule. Default is "Lorentz_Berthelot".
            cutOff (float, optional): The cut-off distance. Default is 12.0.
            shifted (bool, optional): Whether to use shifted potential. Default is False.
            tailCorrections (bool, optional): Whether to apply tail corrections. Default is False.
        """
        super().__init__(**popSelf(locals()))

        # handle double init
        if self._settings["fileName"] is None:
            self.drop_args("fileName")
            self._settings["mixingRule"] = getattr(raspalib.ForceField.MixingRule, self._settings["mixingRule"])
        else:
            self.drop_args(
                "pseudoAtoms", "parameters", "mixingRule", "cutOff", "shifted", "tailCorrections", "useCharge"
            )
        self._cpp_obj = raspalib.ForceField(**self.cpp_args())
        self.useCharge = useCharge

    @classmethod
    def exampleMoleculeForceField(cls, useCharge=True):
        """
        Create an example molecule force field.

        Returns:
            ForceField: An instance of ForceField with example molecule force field settings.
        """
        return cls(
            fileName=os.path.join(SHARE_DIR, "forcefields", "example_molecule_forcefield", "force_field.json"),
            useCharge=useCharge,
        )
