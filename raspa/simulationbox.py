import numpy as np

import raspa.raspalib as raspalib
from .utils import RaspaError, popSelf
from .base import RaspaBase


class SimulationBox(RaspaBase):
    def __init__(self, box: np.ndarray):
        if not (box.shape == (3,) or box.shape == (6,)):
            raise RaspaError(f"Box shape should be either (3,) or (6,) not {self.box.shape}")
        super().__init__(**popSelf(locals()))

        # can not init with kwargs, as "box" is given as either 3 or 6 floats
        self._cpp_obj = raspalib.SimulationBox(*self._settings["box"])

    @property
    def boxType(self):
        return self._cpp_obj.type.name
