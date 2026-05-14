import collections.abc

class PropertyNumberOfMoleculesEvolution():
    """
    Configure and access molecule-count evolution sampling in RASPA.

    This property tracks per-component molecule counts over simulation cycles.
    """

    def __init__(
        self,
        number_of_cycles: int,
        number_of_components: int,
        sample_every: int,
        write_every: int | None = None
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyNumberOfMoleculesEvolution`.

        Args:
            number_of_cycles: Number of simulation cycles to monitor.
            number_of_components: Number of components in the system.
            sample_every: Sampling interval in cycles.
            write_every: Optional output-write interval in cycles.
        """

    @property
    def result(self) -> collections.abc.Sequence[collections.abc.Sequence[int]]:
        """Return sampled molecule counts indexed by cycle and component."""
        ...


