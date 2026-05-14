import collections.abc

class PropertyNumberOfMoleculesHistogram():
    """
    Store number-of-molecules histogram statistics from RASPA.

    This container represents histogram-style occupancy/count statistics over
    molecule-number states.
    """

    def __init__(
        self,
        number_of_blocks: int,
        value_range: tuple[float, float],
        sample_every: int = 1,
        write_every: int | None = None
    ) -> None:
        ...


    def result(self) -> tuple[list[collections.abc.Sequence[float]], \
                              list[collections.abc.Sequence[float]]]:
        """Return histogram bin positions and averaged number of molecules contributions."""
        ...


