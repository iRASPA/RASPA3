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
        number_of_grid_points: tuple[int, int, int] = (128, 128, 128),
        sample_every: int  = 1,
        write_every: int = 5000,
        density_grid_pseudo_atoms_list: collections.abc.Sequence[str] = [],
        normalization_type: PropertyDensityGrid.Normalization = Normalization.MAX,
        binning_mode: PropertyDensityGrid.Binning = Binning.STANDARD
    ) -> None:
        ...
        """
        Initialize a :class:`PropertyDensityGrid`.

        Args:
            number_of_frameworks: Number of frameworks in the system.
            number_of_components: Number of components in the system.
            number_of_grid_points: Grid resolution as ``(nx, ny, nz)``.
                Defaults to ``(128, 128, 128)``.
            sample_every: Sampling interval in cycles. Defaults to ``1``.
            write_every: Output-write interval in cycles. Defaults to ``5000``.
            density_grid_pseudo_atoms_list: Pseudo-atom names for which density
                grids are accumulated. Defaults to an empty sequence.
            normalization_type: Density-grid normalization mode. Defaults to
                :attr:`Normalization.MAX`.
            binning_mode: Density-grid binning mode. Defaults to
                :attr:`Binning.STANDARD`.
        """

