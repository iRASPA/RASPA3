from raspalib.energy_status import *

class PropertyEnergy():
    """
    Configure and access average-energy sampling in RASPA.

    This property accumulates energy statistics and exposes paired
    :class:`EnergyStatus` results (e.g., mean values and associated uncertainty
    or block-averaged companion data).
    """

    # obtained from a System; it has no constructor

    def result(self) -> tuple[EnergyStatus, EnergyStatus]:
        """Return the computed pair of energy-status aggregates."""
        ...

