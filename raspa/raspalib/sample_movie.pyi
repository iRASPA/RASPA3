import typing

class SampleMovie():
    """
    Configure trajectory/movie sampling for a RASPA system.

    This object defines how often snapshots are sampled and whether sampled
    coordinates are wrapped to the simulation box.
    """

    @typing.overload
    def __init__(
        self,
        sample_every: int = 1,
        restrict_to_box: bool = True,
        tag: str | None = None
    ) -> None:
        ...

    @typing.overload
    def __init__(
        self,
        system_id: int,
        sample_every: int,
        restrict_to_box: bool = True,
        tag: str | None = None
    ) -> None:
        ...
        """
        Initialize a :class:`SampleMovie`.

        Args:
            system_id: Not used: the movie is written under the index the
                system has in its driver.
            sample_every: Sampling interval in cycles. Defaults to ``1``.
            restrict_to_box: Whether sampled coordinates are restricted/wrapped
                to the simulation box. Defaults to ``True``.
            tag: Optionally concat filename
        """

    @property
    def sample_every(self) -> int: ...
    @property
    def restrict_to_box(self) -> bool: ...
    @property
    def tag(self) -> str | None: ...
