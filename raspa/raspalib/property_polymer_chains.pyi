import collections.abc
from typing import Any


class PropertyPolymerChains():
    """The size and shape of the chains of the polymerization driver, per
    chain length (the number of molecules of a chain), its topology, whether
    it is living and the conversion at the sample. A chain is measured once
    its topology has been constant for the force ramp plus ``lag_steps``
    steps; the raw sums per key are combined over independent runs with
    ``raspalib.polymer_chains.combine``."""

    def __init__(self, sample_every: int = 1000, write_every: int = 5000, sampling: str = 'periodic',
                 lag_steps: int = 0, conversion_bins: int = 10, minimum_length: int = 2,
                 mass_weighted: bool = False) -> None:
        """``sampling``: 'periodic' (every chain every ``sample_every`` MD
        steps) or 'before_reaction' (a chain just before a reaction changes
        it, and every chain once more at the end of the run)."""
        ...

    @property
    def sample_every(self) -> int: ...
    @property
    def write_every(self) -> int: ...
    @property
    def sampling(self) -> str: ...
    @property
    def lag_steps(self) -> int: ...
    @property
    def conversion_bins(self) -> int: ...
    @property
    def minimum_length(self) -> int: ...
    @property
    def mass_weighted(self) -> bool: ...
    @property
    def sampling_events(self) -> int: ...
    @property
    def skipped_ramping(self) -> int: ...
    @property
    def skipped_spanning(self) -> int: ...

    @property
    def moments(self) -> collections.abc.Sequence[str]:
        """The names of the accumulated per-chain values, in the order of
        'sum' and 'sumOfSquares'."""
        ...

    def entries(self) -> list[dict[str, Any]]:
        """The raw sums per key (length, topology, living, conversion bin), as
        the 'entries' of the JSON output; combine runs with
        raspalib.polymer_chains.combine."""
        ...
