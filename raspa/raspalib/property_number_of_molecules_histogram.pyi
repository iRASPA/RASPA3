import typing

class PropertyNumberOfMoleculesHistogram():
    """
    Store number-of-molecules histogram statistics from RASPA.

    This container represents histogram-style occupancy/count statistics over
    molecule-number states.
    """

    @typing.overload
    def __init__(
        self,
        value_range: tuple[int, int],
        sample_every: int = 1,
        write_every: int | None = None
    ) -> None:
        ...

    @typing.overload
    def __init__(
        self,
        number_of_blocks: int,
        value_range: tuple[int, int],
        sample_every: int,
        write_every: int | None = None
    ) -> None:
        ...
        """
        ``number_of_blocks`` is not needed: the system sizes the histogram for
        its blocks and components when the property is set.
        """


    def result(self) -> list[tuple[list[float], list[float], list[float]]]:
        """Return, per component, the molecule counts of the bins, the averaged histogram and its error."""
        ...


