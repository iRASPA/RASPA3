import enum
import typing
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

    @typing.overload
    def __init__(
        self,
        number_of_grid_points: tuple[int, int, int],
        sample_every: int = 10,
        write_every: int = 5000,
        density_grid_pseudo_atoms_list: collections.abc.Sequence[int] = [],
        normalization_type: Normalization = Normalization.MAX,
        binning_mode: Binning = Binning.STANDARD
    ) -> None:
        ...

    @typing.overload
    def __init__(
        self,
        number_of_frameworks: int,
        number_of_components: int,
        number_of_grid_points: tuple[int, int, int],
        sample_every: int,
        write_every: int,
        density_grid_pseudo_atoms_list: collections.abc.Sequence[int],
        normalization_type: Normalization = Normalization.MAX,
        binning_mode: Binning = Binning.STANDARD
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyDensityGrid`.

        Args:
            number_of_frameworks: Not needed: the system sizes the grid for
                its framework when the property is set.
            number_of_components: Not needed: the system sizes the grid for
                its components when the property is set.
            number_of_grid_points: Grid resolution as ``(nx, ny, nz)``.
            sample_every: Sampling interval in cycles.
            write_every: Output-write interval in cycles.
            density_grid_pseudo_atoms_list: Pseudo-atom indices for which density
                grids are accumulated.
            normalization_type: Density-grid normalization mode. Defaults to
                :attr:`Normalization.MAX`.
            binning_mode: Density-grid binning mode. Defaults to
                :attr:`Binning.STANDARD`.
        """

