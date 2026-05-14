import collections.abc

class PropertyVolumeEvolution():
    """
    Configure and access volume-evolution sampling in RASPA.

    This property tracks simulation-box volume as a function of cycle.
    """

    def __init__(
        self,
        number_of_cycles: int,
        sample_every: int,
        write_every: int | None = None
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyVolumeEvolution`.

        Args:
            number_of_cycles: Number of simulation cycles to monitor.
            sample_every: Sampling interval in cycles.
            write_every: Optional output-write interval in cycles.
        """

    @property
    def result(self) -> collections.abc.Sequence[float]:
        """Return sampled volume values over the simulation trajectory."""
        ...

