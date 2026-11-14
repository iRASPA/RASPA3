import collections.abc

class LoadingData():
    """
    Store component-wise loading metrics produced by RASPA.

    The object bundles molecule counts and associated density representations
    for each simulated component.
    """

    def __init__(self, number_of_components: int) -> None:
        ...
        """
        Initialize a :class:`LoadingData` with zero loadings for ``number_of_components`` components.
        """

    @property
    def inverse_number_densities(self) -> collections.abc.Sequence[float]:
        """Return component-wise inverse number densities."""
        ...

    @property
    def number_densities(self) -> collections.abc.Sequence[float]:
        """Return component-wise number densities."""
        ...

    @property
    def number_of_molecules(self) -> collections.abc.Sequence[int]:
        """Return component-wise molecule counts."""
        ...
