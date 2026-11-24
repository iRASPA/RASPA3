from raspalib.simulation_box import *


class PropertySimulationBox:
    """The block-averaged simulation box."""

    def result(self) -> tuple[SimulationBox, SimulationBox]:
        """Return the average simulation box and its error."""
        ...
