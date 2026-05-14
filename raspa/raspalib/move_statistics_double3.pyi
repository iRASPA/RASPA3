class MoveStatisticsDouble3():
    """
    Track Monte Carlo move statistics for vector-valued move metrics.

    This container stores accepted/constructed counters, running totals, and
    adaptive bounds/targets for three-component move parameters.
    """
    def __init__(
        self
    ) -> None:
        ...
        """
        Initialize a :class:`MoveStatisticsDouble3`.
        """

    @property
    def accepted(self) -> tuple[float, float, float]:
        """Return accepted values/counts in the current window."""
        ...
    @property
    def all_counts(self) -> int:
        """Return the total number of attempted updates counted."""
        ...
    @property
    def constructed(self) -> tuple[float, float, float]:
        """Return constructed/proposed values/counts in the current window."""
        ...
    @property
    def counts(self) -> tuple[float, float, float]:
        """Return sample counts in the current averaging window."""
        ...
    @property
    def lower_limit(self) -> tuple[float, float, float]:
        """Return the lower adaptive bounds for the move parameters."""
        ...
    @lower_limit.setter
    def lower_limit(self, arg0: tuple[float, float, float]) -> None:
        """Set the lower adaptive bounds for the move parameters."""
        ...
    @property
    def max_change(self) -> tuple[float, float, float]:
        """Return maximum parameter changes allowed per adaptation step."""
        ...
    @max_change.setter
    def max_change(self, arg0: tuple[float, float, float]) -> None:
        """Set maximum parameter changes allowed per adaptation step."""
        ...
    @property
    def target_acceptance(self) -> tuple[float, float, float]:
        """Return target acceptance ratios for adaptive tuning."""
        ...
    @target_acceptance.setter
    def target_acceptance(self, arg0: tuple[float, float, float]) -> None:
        """Set target acceptance ratios for adaptive tuning."""
        ...
    @property
    def total_accepted(self) -> tuple[float, float, float]:
        """Return cumulative accepted values/counts over the full run."""
        ...
    @property
    def total_constructed(self) -> tuple[float, float, float]:
        """Return cumulative constructed/proposed values/counts over the run."""
        ...
    @property
    def total_counts(self) -> tuple[float, float, float]:
        """Return cumulative sample counts over the full run."""
        ...
    @property
    def upper_limit(self) -> tuple[float, float, float]:
        """Return the upper adaptive bounds for the move parameters."""
        ...
    @upper_limit.setter
    def upper_limit(self, arg0: tuple[float, float, float]) -> None:
        """Set the upper adaptive bounds for the move parameters."""
        ...


