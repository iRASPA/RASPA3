import numpy as np

import raspalib.raspalib as raspalib
from .utils import RaspaError, popSelf
from .base import RaspaBase


class SimulationBox(RaspaBase):
    def __init__(self, a: float, b: float, c: float):
        super().__init__(**popSelf(locals()))

        # can not init with kwargs, as "box" is given as either 3 or 6 floats
        self._cpp_obj = raspalib.SimulationBox(a, b, c)

    @property
    def boxType(self):
        return self._cpp_obj.type.name
