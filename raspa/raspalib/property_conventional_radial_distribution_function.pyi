import collections.abc

class PropertyConventionalRadialDistributionFunction:
    """
    Configure and access conventional radial distribution functions in RASPA.

    This property controls RDF sampling and provides per-pair distribution data
    accumulated during simulation.
    """

    def __init__(self,
                 number_of_bins: int,
                 range: float,
                 sample_every: int,
                 write_every: int | None = None) -> None:
        ...
        """
        Initialize a :class:`PropertyConventionalRadialDistributionFunction`.

        Args:
            number_of_bins: Number of histogram bins in the RDF.
            range: Maximum distance range of the RDF histogram.
            sample_every: Sampling interval in cycles.
            write_every: Optional output-write interval in cycles.
        """

    def result(self) -> list[list[tuple[collections.abc.Sequence[float],
                                        collections.abc.Sequence[float],
                                        collections.abc.Sequence[float]]]]:
        """Return sampled RDF data for all pseudo-atom pair combinations."""
        ...

