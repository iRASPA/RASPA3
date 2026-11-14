import enum
import collections.abc

class PropertyDensityGrid():
    """
    Configure and access spatial density-grid sampling in RASPA.

    This property controls three-dimensional density histogram collection for
    selected pseudo-atoms across frameworks/components.
    """

    class Binning(enum.IntEnum):
        """Available strategies for assigning samples to grid bins."""

        STANDARD: int = 0
        EQUITABLE: int = 1

    class Normalization(enum.IntEnum):
        """Available normalization modes for exported density grids."""

        MAX: int = 0
        NUMBER_DENSITY: int = 1

    def __init__(
        self,
        number_of_frameworks: int,
        number_of_components: int,
        number_of_grid_points: tuple[int, int, int],
        sample_every: int,
        write_every: int | None,
        density_grid_pseudo_atoms_list: collections.abc.Sequence[str],
        normalization_type: Normalization = Normalization.MAX,
        binning_mode: Binning = Binning.STANDARD
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyDensityGrid`.

        Args:
            number_of_frameworks: Number of frameworks in the system.
            number_of_components: Number of components in the system.
            number_of_grid_points: Grid resolution as ``(nx, ny, nz)``.
            sample_every: Sampling interval in cycles.
            write_every: Optional output-write interval in cycles.
            density_grid_pseudo_atoms_list: Pseudo-atom names for which density
                grids are accumulated.
            normalization_type: Density-grid normalization mode. Defaults to
                :attr:`Normalization.MAX`.
            binning_mode: Density-grid binning mode. Defaults to
                :attr:`Binning.STANDARD`.
        """

