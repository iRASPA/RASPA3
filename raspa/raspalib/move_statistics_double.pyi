class MoveStatisticsDouble():
    """
    Track Monte Carlo move statistics for scalar-valued move metrics.

    This container stores accepted/constructed counters, running totals, and
    adaptive bounds/targets used during move-step tuning.
    """
    def __init__(
        self
    ) -> None:
        ...
        """
        Initialize a :class:`MoveStatisticsDouble`.
        """

    @property
    def accepted(self) -> float:
        """Return the accepted value/count in the current window."""
        ...
    @property
    def all_counts(self) -> int:
        """Return the total number of attempted updates counted."""
        ...
    @property
    def constructed(self) -> float:
        """Return the constructed/proposed value/count in the current window."""
        ...
    @property
    def counts(self) -> float:
        """Return the number of samples in the current averaging window."""
        ...
    @property
    def lower_limit(self) -> float:
        """Return the lower adaptive bound for the move parameter."""
        ...
    @lower_limit.setter
    def lower_limit(self, arg0: float) -> None:
        """Set the lower adaptive bound for the move parameter."""
        ...
    @property
    def max_change(self) -> float:
        """Return the maximum allowed parameter change per adaptation step."""
        ...
    @max_change.setter
    def max_change(self, arg0: float) -> None:
        """Set the maximum allowed parameter change per adaptation step."""
        ...
    @property
    def target_acceptance(self) -> float:
        """Return the target acceptance ratio for adaptive tuning."""
        ...
    @target_acceptance.setter
    def target_acceptance(self, arg0: float) -> None:
        """Set the target acceptance ratio for adaptive tuning."""
        ...
    @property
    def total_accepted(self) -> float:
        """Return cumulative accepted value/count over the full run."""
        ...
    @property
    def total_constructed(self) -> float:
        """Return cumulative constructed/proposed value/count over the full run."""
        ...
    @property
    def total_counts(self) -> float:
        """Return cumulative number of samples over the full run."""
        ...
    @property
    def upper_limit(self) -> float:
        """Return the upper adaptive bound for the move parameter."""
        ...
    @upper_limit.setter
    def upper_limit(self, arg0: float) -> None:
        """Set the upper adaptive bound for the move parameter."""
        ...


