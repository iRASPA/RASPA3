from typing import Literal

import raspalib.raspalib as raspalib
from .base import RaspaBase
from .utils import RASPA_DIR, SHARE_DIR, popSelf
from .pseudo_atom import PseudoAtom
from .vdw_parameters import VDWParameters
import os
import json

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
        pseudoAtoms: list[PseudoAtom] = None,
        parameters: list[VDWParameters] = None,
        mixingRule: Literal["Lorentz_Berthelot"] = "Lorentz_Berthelot",
        cutOffFrameworkVDW: float = 12.0,
        cutOffMoleculeVDW: float = 12.0,
        cutOffCoulomb: float = 12.0,
        shifted: bool = False,
        tailCorrections: bool = False,
        useCharge=True,
    ):
        """
        Initialize the ForceField object with provided parameters.

        Args:
            pseudoAtoms (list[PseudoAtom], optional): A list of pseudo atoms. Default is None.
            parameters (list[VDWParameter], optional): A list of Van der Waals parameters. Default is None.
            mixingRule (Literal["Lorentz_Berthelot"], optional): The mixing rule. Default is "Lorentz_Berthelot".
            cutOffFrameworkVDW (float, optional): The framework-molecule Van der Waals cut-off distance. Default is 12.0.
            cutOffMoleculeVDW (float, optional): The molecule-molecule Van der Waaks cut-off distance. Default is 12.0.
            cutOffCoulomb (float, optional): The cut-off distance for the real-part of the Ewald-summation . Default is 12.0.
            shifted (bool, optional): Whether to use shifted potential. Default is False.
            tailCorrections (bool, optional): Whether to apply tail corrections. Default is False.
            useCharge (bool, optional): Whether to compute electrostatics or not.
        """
        super().__init__(**popSelf(locals()))
        self._settings["mixingRule"] = getattr(raspalib.ForceField.MixingRule, self._settings["mixingRule"])
        self._cpp_obj = raspalib.ForceField(**self.cpp_args())
