import typing

class PropertyRadialDistributionFunction:
    """
    Configure and access force-based radial distribution functions in RASPA.

    The estimator uses the forces on the sites next to their distances, which
    lowers the noise at small bin widths.
    """

    @typing.overload
    def __init__(
        self,
        number_of_bins: int,
        range: float,
        sample_every: int = 10,
        write_every: int = 5000
    ) -> None:
        ...

    @typing.overload
    def __init__(
        self,
        number_of_blocks: int,
        number_of_pseudo_atoms: int,
        number_of_bins: int,
        range: float,
        sample_every: int,
        write_every: int
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyRadialDistributionFunction`.

        Args:
            number_of_blocks: Not needed: the system sizes the histograms for
                its blocks when the property is set.
            number_of_pseudo_atoms: Not needed: the system sizes the histograms
                for the pseudo-atoms of its force field when the property is set.
            number_of_bins: Number of histogram bins in the RDF.
            range: Maximum distance of the RDF histogram [Å].
            sample_every: Sampling interval in cycles.
            write_every: Output-write interval in cycles.
        """

    def result(self, atom_type_a: int, atom_type_b: int) -> tuple[list[float], list[float]]:
        """Return the distances and the RDF of a pair of pseudo-atom types."""
        ...
