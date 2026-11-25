import typing
import collections.abc

class PropertyNumberOfMoleculesEvolution():
    """
    Configure and access molecule-count evolution sampling in RASPA.

    This property tracks per-component molecule counts over simulation cycles.
    """

    @typing.overload
    def __init__(self, sample_every: int = 1, write_every: int | None = None) -> None:
        ...

    @typing.overload
    def __init__(
        self,
        number_of_cycles: int,
        number_of_components: int | None = None,
        sample_every: int = 1,
        write_every: int | None = None
    ) -> None:
        ...

    @typing.overload
    def __init__(
        self,
        number_of_production_cycles: int,
        number_of_components: int | None = None,
        sample_every: int = 1,
        write_every: int | None = None
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyNumberOfMoleculesEvolution`.

        Args:
            number_of_cycles: The keyword of the 3.1.0 release for
                ``number_of_production_cycles``.
            number_of_production_cycles: Number of cycles to record, counted
                over all phases of the driver. Not needed: without it the
                driver records the cycles of all its phases.
            number_of_components: Not needed: the system sizes the property
                for its components when it is set.
            sample_every: Sampling interval in cycles.
            write_every: Optional output-write interval in cycles.

        The object given to the setter of a system holds the results when a
        phase of the driver returns.
        """

    @property
    def sample_every(self) -> int: ...
    @property
    def write_every(self) -> int | None: ...

    @property
    def result(self) -> collections.abc.Sequence[collections.abc.Sequence[int]]:
        """Return sampled molecule counts indexed by cycle and component."""
        ...


