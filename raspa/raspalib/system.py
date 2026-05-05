import collections.abc
from raspalib.force_field import *
from raspalib.mc_move_probabilities import *
from raspalib.mc_move_statistics import *
from raspalib.property_energy import *
from raspalib.property_conserved_energy_evolution import *
from raspalib.framework import *
from raspalib.component import *
from raspalib.simulation_box import *
from raspalib.loading_data import *
from raspalib.property_loading import *
from raspalib.property_enthalpy import *
from raspalib.sample_movie import *
from raspalib.property_energy import *
from raspalib.property_pressure import *
from raspalib.property_density_grid import *
from raspalib.property_energy_histogram import *
from raspalib.property_number_of_molecules_histogram import *
from raspalib.property_conventional_radial_distribution_function import *
from raspalib.property_mean_squared_displacement import *
from raspalib.property_velocity_autocorrelation_function import *
from raspalib.thermostat import *
from raspalib.property_number_of_molecules_evolution import *
from raspalib.property_volume_evolution import *
from raspalib.property_conserved_energy_evolution import *


class System():
    """
    Represent a full simulation system configuration in RASPA.

    A system combines force-field settings, box geometry, framework/component
    definitions, initial state, move probabilities, and optional sampled
    properties.
    """

    def __init__(
        self,
        force_field: ForceField,
        simulation_box: SimulationBox | None = None,
        has_external_field: bool = False,
        external_temperature: float = 300.0,
        external_pressure: float | None = None,
        helium_void_fraction: float = 0.29,
        framework_components: Framework | None = None,
        components: collections.abc.Sequence[Component] = [],
        initial_positions: collections.abc.Sequence[tuple[float, float, float]] = [],
        initial_number_of_molecules: collections.abc.Sequence[int] = [],
        number_of_blocks: int = 5,
        system_probabilities: MCMoveProbabilities = MCMoveProbabilities()
    ) -> None:
        ...
        """
        Initialize a :class:`System`.

        Args:
            force_field: Force-field definition used by the system.
            simulation_box: Optional simulation-box geometry.
            has_external_field: Whether an external field is enabled.
            external_temperature: External-field temperature (K).
            external_pressure: Optional external pressure.
            helium_void_fraction: Helium void fraction used for related
                corrections/analysis.
            framework_components: Optional framework component definition.
            components: Sequence of mobile components in the system.
            initial_positions: Optional initial positions.
            initial_number_of_molecules: Initial molecule counts per component.
            number_of_blocks: Number of blocks used for block averaging.
            system_probabilities: Monte Carlo move probabilities at the system
                level.
        """

    @property
    def components(self) -> collections.abc.Sequence[Component]:
        """Return the configured list of components."""
        ...

    @property
    def number_of_molecules_per_component(self) -> collections.abc.Sequence[int]:
        """Return molecule counts for each configured component."""
        ...

    def framework_mass(self) -> float | None:
        """Return the total framework mass, if a framework is present."""
        ...

    @property
    def mc_moves_statistics(self) -> MCMoveStatistics:
        """Return Monte Carlo move acceptance/statistics counters."""
        ...

    @property
    def simulation_box(self) -> SimulationBox:
        """Return the system simulation box."""
        ...

    @property
    def loadings(self) -> LoadingData:
        """Return current loading data."""
        ...

    @property
    def average_loadings(self) -> PropertyLoading:
        """Return the average-loading property accessor."""
        ...

    @property
    def average_enthalpies_of_adsorption(self) -> PropertyEnthalpy:
        """Return the average enthalpy-of-adsorption property accessor."""
        ...

    @property
    def sample_pdb_movie(self) -> SampleMovie | None:
        """Return configured PDB-movie sampling settings, if enabled."""
        ...

    @property
    def average_energies(self) -> PropertyEnergy:
        """Return the average-energy property accessor."""
        ...
    @property
    def average_pressure(self) -> PropertyPressure:
        """Return the average-pressure property accessor."""
        ...

    @property
    def property_density_grid(self) -> PropertyDensityGrid | None:
        """Return the density-grid property, if configured."""
        ...

    @property
    def average_energy_histogram(self) -> PropertyEnergyHistogram | None:
        """Return the energy-histogram property, if configured."""
        ...

    @property
    def average_number_of_molecules_histogram(self) -> PropertyNumberOfMoleculesHistogram | None:
        """Return the molecule-count histogram property, if configured."""
        ...

    @property
    def property_conventional_rdf(self) -> PropertyConventionalRadialDistributionFunction | None:
        """Return the conventional radial distribution function property."""
        ...

    @property
    def property_msd(self) -> PropertyMeanSquaredDisplacement | None:
        """Return the mean-squared-displacement property, if configured."""
        ...

    @property
    def property_vacf(self) -> PropertyVelocityAutoCorrelationFunction | None:
        """Return the velocity-autocorrelation property, if configured."""
        ...

    @property
    def property_number_of_molecules_evolution(self) -> PropertyNumberOfMoleculesEvolution | None:
        """Return molecule-count evolution property, if configured."""
        ...

    @property
    def property_volume_evolution(self) -> PropertyVolumeEvolution | None:
        """Return volume-evolution property, if configured."""
        ...

    @property
    def property_conserved_energy_evolution(self) -> PropertyConservedEnergyEvolution | None:
        """Return conserved-energy evolution property, if configured."""
        ...

    def set_thermostat(self, rdf: Thermostat) -> None:
        """Attach thermostat settings to the system."""
        ...

    def set_sample_pdb_movie(self, sample_movie: SampleMovie) -> None:
        """Enable/configure PDB-movie sampling."""
        ...

    def set_average_energy_histogram(self, property: PropertyEnergyHistogram) -> None:
        """Attach an energy-histogram property configuration."""
        ...

    def set_property_density_grid(self, property: PropertyDensityGrid) -> None:
        """Attach a density-grid property configuration."""
        ...

    def set_property_number_of_molecules_evolution(self, property: PropertyNumberOfMoleculesEvolution) -> None:
        """Attach molecule-count evolution property configuration."""
        ...

    def set_property_volume_evolution(self, property: PropertyVolumeEvolution) -> None:
        """Attach volume-evolution property configuration."""
        ...

    def set_property_conserved_energy_evolution(self, property: PropertyConservedEnergyEvolution) -> None:
        """Attach conserved-energy evolution property configuration."""
        ...

    def set_property_conventional_rdf(self, property: PropertyConventionalRadialDistributionFunction) -> None:
        """Attach conventional RDF property configuration."""
        ...

    def set_property_msd(self, property: PropertyMeanSquaredDisplacement) -> None:
        """Attach mean-squared-displacement property configuration."""
        ...

    def set_property_msd(self, property: PropertyVelocityAutoCorrelationFunction) -> None:
        """Attach velocity-autocorrelation property configuration."""
        ...

