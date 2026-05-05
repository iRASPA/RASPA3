import collections.abc

from raspalib.average_energy_type import *

class PropertyEnergyHistogram():
    """
    Configure and access energy-histogram sampling in RASPA.

    This property collects histogrammed energy statistics over simulation
    cycles and exposes averaged energy contributions per histogram bin.
    """

    def __init__(
        self,
        number_of_blocks: int,
        number_of_bins: int,
        value_range: tuple[float, float],
        sample_every: int = 1,
        write_every: int = 5000
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyEnergyHistogram`.

        Args:
            number_of_blocks: Number of blocks used for block averaging.
            number_of_bins: Number of histogram bins.
            value_range: Inclusive energy range ``(min_value, max_value)``
                covered by the histogram.
            sample_every: Sampling interval in cycles. Defaults to ``1``.
            write_every: Output-write interval in cycles. Defaults to ``5000``.
        """

    def result(self) -> tuple[collections.abc.Sequence[float], \
                              collections.abc.Sequence[AverageEnergyType], \
                              collections.abc.Sequence[AverageEnergyType]]:
        """Return histogram bin positions and averaged energy contributions."""
        ...
