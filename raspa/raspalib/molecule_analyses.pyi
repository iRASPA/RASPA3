import collections.abc
import enum

# The analyses of the molecules of a component. A component asks for one with
# settings (``component.molecule_shape_settings = MoleculeShapeSettings(...)``);
# the driver creates the analyses of its systems, and the results are read from
# ``property_molecule_properties``, ``property_molecule_shape``,
# ``property_molecule_backbone`` and ``property_end_to_end_acf`` of a system.
# A histogram comes as the centers of its bins, the probability density, and
# the error of it.


class MoleculePropertiesSettings():
    """Histograms of the bond lengths, the bend and torsion angles, and the
    end-to-end distance of the molecules of a component."""

    def __init__(self, sample_every: int = 10, write_every: int | None = None, number_of_bins: int = 128,
                 bond_range: float = 4.0, end_to_end_range: float | None = None,
                 end_to_end_bin_width: float = 0.25) -> None:
        """
        Args:
            sample_every: Sampling interval in cycles.
            write_every: Interval of writing the histograms to files; ``None``
                writes nothing.
            number_of_bins: Bins of the bond, bend and torsion histograms.
            bond_range: Upper limit of the bond-length histogram [Angstrom].
            end_to_end_range: Upper limit of the end-to-end histogram
                [Angstrom]; by default the contour length.
            end_to_end_bin_width: Bin width of that histogram [Angstrom].
        """
        ...

    @property
    def sample_every(self) -> int: ...
    @property
    def write_every(self) -> int | None: ...
    @property
    def number_of_bins(self) -> int: ...
    @property
    def bond_range(self) -> float: ...
    @property
    def end_to_end_range(self) -> float | None: ...
    @property
    def end_to_end_bin_width(self) -> float: ...


class MoleculeShapeSettings():
    """The radius of gyration and the shape descriptors of the gyration
    tensor of the molecules of a component."""

    def __init__(self, sample_every: int = 10, write_every: int | None = None, number_of_bins: int = 128,
                 mass_weighted: bool = False, radius_of_gyration_range: float | None = None,
                 radius_of_gyration_bin_width: float = 0.1) -> None:
        """
        Args:
            sample_every: Sampling interval in cycles.
            write_every: Interval of writing to files; ``None`` writes nothing.
            number_of_bins: Bins of the shape-anisotropy and prolateness
                histograms.
            mass_weighted: Weigh the atoms by their mass instead of equally.
            radius_of_gyration_range: Upper limit of the radius-of-gyration
                histogram [Angstrom]; by default from the contour length.
            radius_of_gyration_bin_width: Bin width of it [Angstrom].
        """
        ...

    @property
    def sample_every(self) -> int: ...
    @property
    def write_every(self) -> int | None: ...
    @property
    def number_of_bins(self) -> int: ...
    @property
    def mass_weighted(self) -> bool: ...
    @property
    def radius_of_gyration_range(self) -> float | None: ...
    @property
    def radius_of_gyration_bin_width(self) -> float: ...


class MoleculeBackboneSettings():
    """The statistics of the backbone of a chain: internal distances, bond
    correlations and the single-chain form factor."""

    def __init__(self, sample_every: int = 10, write_every: int | None = None, number_of_wave_vectors: int = 64,
                 wave_vector_lower_limit: float = 0.01, wave_vector_upper_limit: float = 5.0) -> None:
        """
        Args:
            sample_every: Sampling interval in cycles.
            write_every: Interval of writing to files; ``None`` writes nothing.
            number_of_wave_vectors: Points of the logarithmic grid of wave
                vectors of the form factor.
            wave_vector_lower_limit: Smallest wave vector [1/Angstrom].
            wave_vector_upper_limit: Largest wave vector [1/Angstrom].
        """
        ...

    @property
    def sample_every(self) -> int: ...
    @property
    def write_every(self) -> int | None: ...
    @property
    def number_of_wave_vectors(self) -> int: ...
    @property
    def wave_vector_lower_limit(self) -> float: ...
    @property
    def wave_vector_upper_limit(self) -> float: ...


class EndToEndACFSettings():
    """The autocorrelation function of the end-to-end vector of the molecules
    of a component."""

    def __init__(self, sample_every: int = 10, write_every: int | None = None,
                 number_of_block_elements: int = 25) -> None:
        """
        Args:
            sample_every: Sampling interval in cycles.
            write_every: Interval of writing to files; ``None`` writes nothing.
            number_of_block_elements: Elements per block of the order-N scheme.
        """
        ...

    @property
    def sample_every(self) -> int: ...
    @property
    def write_every(self) -> int | None: ...
    @property
    def number_of_block_elements(self) -> int: ...


class PropertyMoleculeProperties():
    """The histograms of the bonded terms of the molecules."""

    @property
    def number_of_bonds(self) -> collections.abc.Sequence[int]:
        """The number of bonds of a molecule, per component."""
        ...

    @property
    def number_of_bends(self) -> collections.abc.Sequence[int]: ...

    @property
    def number_of_torsions(self) -> collections.abc.Sequence[int]: ...

    def bond_histogram(self, component: int, bond: int) -> tuple[collections.abc.Sequence[float], collections.abc.Sequence[float], collections.abc.Sequence[float]]:
        """The histogram of the length of a bond [Angstrom], by the index of
        its bond potential."""
        ...

    def bend_histogram(self, component: int, bend: int) -> tuple[collections.abc.Sequence[float], collections.abc.Sequence[float], collections.abc.Sequence[float]]:
        """The histogram of a bend angle [degrees]."""
        ...

    def torsion_histogram(self, component: int, torsion: int) -> tuple[collections.abc.Sequence[float], collections.abc.Sequence[float], collections.abc.Sequence[float]]:
        """The histogram of a torsion angle, from -180 to 180 [degrees]."""
        ...

    def end_to_end_histogram(self, component: int) -> tuple[collections.abc.Sequence[float], collections.abc.Sequence[float], collections.abc.Sequence[float]]:
        """The histogram of the end-to-end distance [Angstrom]."""
        ...


class PropertyMoleculeShape():
    """The radius of gyration and the shape of the molecules."""

    class Moment(enum.IntEnum):
        """The averages that are kept: lengths in Angstrom, the eigenvalues of the gyration tensor, the asphericity and the acylindricity in Angstrom^2."""
        RADIUS_OF_GYRATION: int = 0
        RADIUS_OF_GYRATION_SQUARED: int = 1
        RADIUS_OF_GYRATION_FOURTH: int = 2
        LAMBDA_1: int = 3
        LAMBDA_2: int = 4
        LAMBDA_3: int = 5
        ASPHERICITY: int = 6
        ACYLINDRICITY: int = 7
        SHAPE_ANISOTROPY: int = 8
        SECOND_INVARIANT: int = 9
        PROLATENESS: int = 10
        INVERSE_HYDRODYNAMIC_RADIUS: int = 11
        END_TO_END_SQUARED: int = 12
        TENSOR_XX: int = 13
        TENSOR_YY: int = 14
        TENSOR_ZZ: int = 15
        TENSOR_XY: int = 16
        TENSOR_XZ: int = 17
        TENSOR_YZ: int = 18

    def result(self, component: int, moment: 'PropertyMoleculeShape.Moment') -> tuple[float, float]:
        """The average of a moment over the molecules of the component, and
        its error."""
        ...

    def radius_of_gyration_histogram(self, component: int) -> tuple[collections.abc.Sequence[float], collections.abc.Sequence[float], collections.abc.Sequence[float]]:
        """The histogram of the radius of gyration [Angstrom]."""
        ...

    def shape_anisotropy_histogram(self, component: int) -> tuple[collections.abc.Sequence[float], collections.abc.Sequence[float], collections.abc.Sequence[float]]:
        """The histogram of the relative shape anisotropy, from 0 to 1."""
        ...

    def prolateness_histogram(self, component: int) -> tuple[collections.abc.Sequence[float], collections.abc.Sequence[float], collections.abc.Sequence[float]]:
        """The histogram of the prolateness, from -1/4 to 2."""
        ...


class PropertyMoleculeBackbone():
    """The statistics of the backbone of chain molecules."""

    class Moment(enum.IntEnum):
        """The averages that are kept of the backbone."""
        BOND_LENGTH: int = 0
        BOND_LENGTH_SQUARED: int = 1
        END_TO_END_SQUARED: int = 2
        PROJECTION: int = 3
        RADIUS_OF_GYRATION_SQUARED: int = 4

    class Averages():
        """The averages over the molecules of a component and over the run."""

        @property
        def internal_distance_squared(self) -> collections.abc.Sequence[float]:
            """The mean squared distance of backbone beads k bonds apart
            [Angstrom^2], for k = 1, 2, ..."""
            ...

        @property
        def bond_correlation(self) -> collections.abc.Sequence[float]:
            """The correlation of the directions of backbone bonds k bonds
            apart, for k = 0, 1, ..."""
            ...

        @property
        def form_factor(self) -> collections.abc.Sequence[float]:
            """The single-chain form factor at the wave vectors of
            ``wave_vectors``."""
            ...

    def result(self, component: int) -> 'PropertyMoleculeBackbone.Averages': ...

    def moment(self, component: int, moment: 'PropertyMoleculeBackbone.Moment') -> tuple[float, float]:
        """The average of a moment and its error."""
        ...

    def persistence_length_from_projection(self, component: int) -> tuple[float, float]:
        """The persistence length from the projection of the chain on its end
        bonds, and its error [Angstrom]."""
        ...

    def persistence_length_from_fit(self, component: int) -> tuple[float, float]:
        """The persistence length from the decay of the bond correlation, and
        its error [Angstrom]."""
        ...

    def flory_exponent(self, component: int) -> tuple[float, float]:
        """The Flory exponent from the internal distances, and its error."""
        ...

    def wave_vectors(self, component: int) -> collections.abc.Sequence[float]:
        """The wave vectors of the form factor [1/Angstrom]."""
        ...

    def backbone(self, component: int) -> collections.abc.Sequence[int]:
        """The atoms of the backbone of the component."""
        ...

    def contour_length(self, component: int) -> float:
        """The contour length of the backbone [Angstrom]."""
        ...


class EndToEndAutoCorrelationFunctionData():
    """One lag of the end-to-end autocorrelation function."""

    @property
    def time(self) -> float:
        """The lag: picoseconds, or cycles for Monte Carlo."""
        ...

    @property
    def acf(self) -> float:
        """<R(0).R(t)> [Angstrom^2]."""
        ...

    @property
    def normalized(self) -> float:
        """<R(0).R(t)> / <R^2>."""
        ...

    @property
    def number_of_samples(self) -> float: ...
    @property
    def end_to_end_squared_acf(self) -> float | None: ...
    @property
    def radius_of_gyration_squared_acf(self) -> float | None: ...


class RelaxationTimeEstimates():
    """Estimates of a relaxation time, in the unit of the lags."""

    @property
    def integrated(self) -> float | None: ...
    @property
    def integrated_is_lower_bound(self) -> bool: ...
    @property
    def one_over_e(self) -> float | None: ...
    @property
    def exponential_fit(self) -> float | None: ...


class ScalarRelaxation():
    """The relaxation of a scalar of the molecule."""

    @property
    def available(self) -> bool: ...
    @property
    def mean(self) -> float: ...
    @property
    def variance(self) -> float: ...
    @property
    def times(self) -> RelaxationTimeEstimates: ...


class EndToEndRelaxationTimes():
    """The relaxation times of the end-to-end vector."""

    @property
    def mean_squared_end_to_end(self) -> float:
        """<R^2> [Angstrom^2]."""
        ...

    @property
    def longest_lag(self) -> float: ...

    @property
    def integrated(self) -> float | None:
        """The integral of the normalized function."""
        ...

    @property
    def integrated_is_lower_bound(self) -> bool: ...

    @property
    def one_over_e(self) -> float | None:
        """The lag at which the function has decayed to 1/e."""
        ...

    @property
    def exponential_fit(self) -> float | None: ...
    @property
    def end_to_end_squared(self) -> ScalarRelaxation: ...
    @property
    def radius_of_gyration_squared(self) -> ScalarRelaxation: ...


class PropertyEndToEndAutoCorrelationFunction():
    """The autocorrelation function of the end-to-end vector."""

    def has_data(self, component: int) -> bool: ...

    def result(self, component: int) -> collections.abc.Sequence[EndToEndAutoCorrelationFunctionData]:
        """The autocorrelation function of the component, per lag."""
        ...

    def relaxation_times(self, component: int) -> EndToEndRelaxationTimes:
        """The relaxation times of the end-to-end vector, in the unit of the
        lags."""
        ...

    @property
    def time_unit(self) -> str:
        """The unit of the lags: ``ps`` for molecular dynamics, ``cycles`` for
        Monte Carlo."""
        ...
