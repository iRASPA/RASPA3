import collections.abc
import typing

class PropertyLambdaProbabilityHistogram():
    """
    Configure and access lambda-probability histogram statistics in RASPA.

    This property stores raw and normalized lambda occupancy information used in
    expanded-ensemble or free-energy style workflows.
    """

    @typing.overload
    def __init__(self) -> None:
        """Initialize an empty :class:`PropertyLambdaProbabilityHistogram`."""
        ...

    @typing.overload
    def __init__(
        self,
        number_of_blocks: int,
        number_of_sample_points: int
    ) -> None:
        """
        Initialize a :class:`PropertyLambdaProbabilityHistogram`.

        Args:
            number_of_blocks: Number of blocks used for block averaging.
            number_of_sample_points: Number of lambda sample points/bins.
        """
        ...

    def result(self) -> tuple[list[float], list[float]]:
        """Return the normalized lambda probability per bin and its confidence interval (dimensionless)."""
        ...

    def average_dudlambda(self) -> tuple[list[float], list[float]]:
        """Return the average dU/dlambda per lambda bin and its confidence interval, in Kelvin."""
        ...

    @property
    def bias_factor(self) -> collections.abc.Sequence[float]:
        """Return the lambda-dependent bias factors."""
        ...

    @bias_factor.setter
    def bias_factor(self, arg0: collections.abc.Sequence[float]) -> None:
        """Set the lambda-dependent bias factors."""
        ...

    @property
    def histogram(self) -> collections.abc.Sequence[float]:
        """Return raw histogram counts/probabilities per lambda bin."""
        ...

    @property
    def occupancy_count(self) -> float:
        """Return the counted number of occupied histogram samples."""
        ...

    @property
    def occupancy_total(self) -> float:
        """Return the total occupancy accumulator for normalization."""
        ...

