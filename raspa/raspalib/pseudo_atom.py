from typing import Literal

import raspalib.raspalib as raspalib
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
        frameworkType: bool = False,
        mass: float = 1.0,
        charge: float = 0.0,
        polarizability: float = 0.0,
        atomicNumber: int = 8,
        printToPDB: bool = True,
        source: str = "-"
    ):
        """
        Initialize the PseudoAtom object with provided parameters.

        Args:
            name (str): The name of the pseudo atom. Default is "C".
            frameworkType (bool): Whether it is a framework-atom or not
            mass (float): The mass of the pseudo atom. Default is 1.0.
            charge (float): The charge of the pseudo atom. Default is 0.0.
            polarizability (float): The polarizability of the pseudo atom. Default is 0.0.
            atomicNumber (int): The atomic number of the pseudo atom. Default is 8.
            printToPDB (bool): Whether to print to PDB. Default is False.
            source (str): The source of the pseudo atom. Default is "-".
        """
        super().__init__(**popSelf(locals()))
        self._cpp_obj = raspalib.PseudoAtom(**self.cpp_args())
