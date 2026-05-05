import collections.abc

class PropertyMeanSquaredDisplacement:
    """
    Configure and access mean-squared-displacement sampling in RASPA.

    This property tracks MSD evolution and exposes sampled time/displacement
    series used for diffusion analysis.
    """

    def __init__(self,
                 sample_every: int,
                 write_every: int | None) -> None:
        ...
        """
        Initialize a :class:`PropertyMeanSquaredDisplacement`.

        Args:
            sample_every: Sampling interval in cycles.
            write_every: Optional output-write interval in cycles.
        """

    @property
    def result(self) -> tuple[collections.abc.Sequence[float],
                              collections.abc.Sequence[float]]:
        """Return sampled time points and mean-squared-displacement values."""
        ...

