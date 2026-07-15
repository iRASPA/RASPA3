from raspalib.energy_dudlambda import *

class EnergyStatus():
    """
    Represent the current energy-state summary reported by RASPA.

    This object exposes aggregate energy quantities for a configuration,
    including the total energy packaged as an :class:`EnergyDuDlambda`.
    """

    def __init__(
        self
    ) -> None:
        ...
        """
        Initialize an :class:`EnergyStatus`.
        """

    @property
    def total_energy(self) -> EnergyDuDlambda:
        """Return the total energy and corresponding ``dU/dlambda`` factor."""
        ...


