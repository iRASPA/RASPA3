import collections.abc

class PropertyConventionalRadialDistributionFunction:
    """
    Configure and access conventional radial distribution functions in RASPA.

    This property controls RDF sampling and provides per-pair distribution data
    accumulated during simulation.
    """

    def __init__(self,
                 numberOf_blocks: int,
                 number_of_pseudo_atoms: int,
                 number_of_bins: int,
                 range: float,
                 sample_every: int,
                 write_every: int) -> None:
        ...
        """
        Initialize a :class:`PropertyConventionalRadialDistributionFunction`.

        Args:
            numberOf_blocks: Number of blocks used for block averaging.
            number_of_pseudo_atoms: Number of pseudo-atom types considered.
            number_of_bins: Number of histogram bins in the RDF.
            range: Maximum distance range of the RDF histogram.
            sample_every: Sampling interval in cycles.
            write_every: Output-write interval in cycles.
        """

    @property
    def result(self) -> list[list[tuple[collections.abc.Sequence[float],
                                        collections.abc.Sequence[float],
                                        collections.abc.Sequence[float]]]]:
        """Return sampled RDF data for all pseudo-atom pair combinations."""
        ...

