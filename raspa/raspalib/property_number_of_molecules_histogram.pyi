class PropertyNumberOfMoleculesHistogram():
    """
    Store number-of-molecules histogram statistics from RASPA.

    This container represents histogram-style occupancy/count statistics over
    molecule-number states.
    """

    def __init__(
        self,
        number_of_blocks: int,
        value_range: tuple[int, int],
        sample_every: int,
        write_every: int | None = None
    ) -> None:
        ...


    def result(self) -> list[tuple[list[float], list[float], list[float]]]:
        """Return, per component, the molecule counts of the bins, the averaged histogram and its error."""
        ...


