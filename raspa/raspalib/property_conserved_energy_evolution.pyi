import typing
import collections.abc
from raspalib.running_energy import *

class PropertyConservedEnergyEvolution():
    """
    Configure and access conserved-energy evolution sampling in RASPA.

    This property tracks conserved-energy values over simulation cycles and
    exposes the sampled time-series as :class:`RunningEnergy` entries.
    """

    @typing.overload
    def __init__(self, sample_every: int = 1, write_every: int | None = None) -> None:
        ...

    @typing.overload
    def __init__(
        self,
        number_of_cycles: int,
        sample_every: int = 1,
        write_every: int | None = None
    ) -> None:
        ...

    @typing.overload
    def __init__(
        self,
        number_of_production_cycles: int,
        sample_every: int = 1,
        write_every: int | None = None
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyConservedEnergyEvolution`.

        Args:
            number_of_cycles: The keyword of the 3.1.0 release for
                ``number_of_production_cycles``.
            number_of_production_cycles: Number of cycles to record, counted
                over all phases of the driver. Not needed: without it the
                driver records the cycles of all its phases.
            sample_every: Sampling interval in cycles.
            write_every: Optional output-write interval in cycles.

        The object given to the setter of a system holds the results when a
        phase of the driver returns.
        """

    @property
    def sample_every(self) -> int: ...
    @property
    def write_every(self) -> int | None: ...

    def result(self) -> collections.abc.Sequence[RunningEnergy]:
        """Return the sampled running energies, in internal units like every :class:`RunningEnergy`."""
        ...

