class SampleMovie():
    """
    Configure trajectory/movie sampling for a RASPA system.

    This object defines how often snapshots are sampled and whether sampled
    coordinates are wrapped to the simulation box.
    """

    def __init__(
        self,
        system_id: int,
        sample_every: int = 1,
        restrict_to_box: bool = True,
        tag: str | None
    ) -> None:
        ...
        """
        Initialize a :class:`SampleMovie`.

        Args:
            system_id: Identifier of the system to sample.
            sample_every: Sampling interval in cycles. Defaults to ``1``.
            restrict_to_box: Whether sampled coordinates are restricted/wrapped
                to the simulation box. Defaults to ``True``.
            tag: Optionally concat filename
        """
