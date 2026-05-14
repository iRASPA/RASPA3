import collections.abc

class LoadingData():
    """
    Store component-wise loading metrics produced by RASPA.

    The object bundles molecule counts and associated density representations
    for each simulated component.
    """

    def __init__(
        self,
        number_of_molecules: collections.abc.Sequence[int],
        number_densities: collections.abc.Sequence[float],
        inverse_number_densities: collections.abc.Sequence[float]
    ) -> None:
        ...
        """
        Initialize a :class:`LoadingData`.

        Args:
            number_of_molecules: Number of molecules for each component.
            number_densities: Number density for each component.
            inverse_number_densities: Inverse number density for each component.
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
