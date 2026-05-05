import collections.abc
from raspalib.running_energy import *

class PropertyConservedEnergyEvolution():
    """
    Configure and access conserved-energy evolution sampling in RASPA.

    This property tracks conserved-energy values over simulation cycles and
    exposes the sampled time-series as :class:`RunningEnergy` entries.
    """

    def __init__(
        self,
        number_of_cycles: int,
        sample_every: int,
        write_every: int | None = None
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyConservedEnergyEvolution`.

        Args:
            number_of_cycles: Number of simulation cycles to monitor.
            sample_every: Sampling interval in cycles.
            write_every: Optional output-write interval in cycles.
        """

    @property
    def result(self) -> collections.abc.Sequence[RunningEnergy]:
        """Return sampled conserved-energy evolution data."""
        ...

