from raspalib.reaction import *
from raspalib.molecule_analyses import *
from raspalib.property_elastic_constants_fluctuation import *
from raspalib.property_radical_waiting_times import *
from raspalib.write_lammps_data import *
from raspalib.running_energy import *
from raspalib.thermobarostat import *
from raspalib.property_temperature import *
from raspalib.property_simulation_box import *
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
from raspalib.property_partial_molar_properties import *
from raspalib.sample_movie import *
from raspalib.property_energy import *
from raspalib.property_pressure import *
from raspalib.property_density_grid import *
from raspalib.property_energy_histogram import *
from raspalib.property_number_of_molecules_histogram import *
from raspalib.property_conventional_radial_distribution_function import *
from raspalib.property_radial_distribution_function import *
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
        external_temperature: float = 298.0,
        external_pressure: float | None = None,
        helium_void_fraction: float = 0.0,
        framework_components: Framework | None = None,
        components: collections.abc.Sequence[Component] = [],
        initial_positions: collections.abc.Sequence[collections.abc.Sequence[tuple[float, float, float]]] = [],
        initial_number_of_molecules: collections.abc.Sequence[int] = [],
        number_of_blocks: int | None = None,
        system_probabilities: MCMoveProbabilities = MCMoveProbabilities(),
        reactions: collections.abc.Sequence[Reaction] = []
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
            initial_positions: Optional initial positions, per component the
                atom positions of its initial molecules.
            initial_number_of_molecules: Initial molecule counts per component.
            number_of_blocks: Number of blocks used for block averaging. Not
                needed: without it the system takes the number of its
                components (which then have to agree), and a driver sizes its
                systems for its own. When given, the components are sized
                for it.
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
    def average_partial_molar_properties(self) -> PropertyPartialMolarProperties:
        """Return the average partial-molar-properties accessor."""
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
    def property_number_of_molecules_histogram(self) -> PropertyNumberOfMoleculesHistogram | None:
        """Return the molecule-count histogram property, if configured."""
        ...

    @property
    def average_energy_histogram(self) -> PropertyEnergyHistogram | None:
        """Return the energy-histogram property, if configured."""
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

    def set_thermostat(self, thermostat: Thermostat | None = None) -> None:
        """Attach thermostat settings to the system."""
        ...

    # The setters of the properties. The system builds its own property from the
    # settings of the one given, sized for its molecules, components and blocks.
    # The given object is brought up to date when it is set and whenever a phase
    # of a driver returns, so the results can be read from it as well as from
    # ``system.property_...``. An object given to more than one system of a
    # driver is left as it was. None removes the property.

    def set_sample_pdb_movie(self, property: SampleMovie | None = None) -> None:
        """Enable/configure PDB-movie sampling."""
        ...

    def set_number_of_molecules_histogram(self, property: PropertyNumberOfMoleculesHistogram | None = None) -> None:
        """Attach a number of molecules-histogram property configuration."""
        ...

    def set_average_energy_histogram(self, property: PropertyEnergyHistogram | None = None) -> None:
        """Attach an energy-histogram property configuration."""
        ...

    def set_property_density_grid(self, property: PropertyDensityGrid | None = None) -> None:
        """Attach a density-grid property configuration."""
        ...

    def set_property_number_of_molecules_evolution(self, property: PropertyNumberOfMoleculesEvolution | None = None) -> None:
        """Attach molecule-count evolution property configuration."""
        ...

    def set_property_volume_evolution(self, property: PropertyVolumeEvolution | None = None) -> None:
        """Attach volume-evolution property configuration."""
        ...

    def set_property_conserved_energy_evolution(self, property: PropertyConservedEnergyEvolution | None = None) -> None:
        """Attach conserved-energy evolution property configuration."""
        ...

    def set_property_conventional_rdf(self, property: PropertyConventionalRadialDistributionFunction | None = None) -> None:
        """Attach conventional RDF property configuration."""
        ...

    def set_property_msd(self, property: PropertyMeanSquaredDisplacement | None = None) -> None:
        """Attach mean-squared-displacement property configuration."""
        ...

    def set_property_vacf(self, property: PropertyVelocityAutoCorrelationFunction | None = None) -> None:
        """Attach velocity-autocorrelation property configuration."""
        ...

    @property
    def property_rdf(self) -> PropertyRadialDistributionFunction | None:
        """Return the force-based RDF property, if configured."""
        ...

    def set_property_rdf(self, property: PropertyRadialDistributionFunction | None = None) -> None:
        """Attach force-based RDF property configuration."""
        ...

    @property
    def temperature(self) -> float:
        """The external temperature [K]."""
        ...

    @property
    def external_pressure_diagonal(self) -> tuple[float, float, float]:
        """The external pressure along x, y and z [Pa], for the anisotropic
        volume move."""
        ...

    @external_pressure_diagonal.setter
    def external_pressure_diagonal(self, value: tuple[float, float, float]) -> None:
        ...

    @property
    def molecular_dynamics_ensemble(self) -> MolecularDynamicsEnsemble:
        """The ensemble of molecular dynamics. 'set_thermobarostat' sets it
        to that of the barostat."""
        ...

    @molecular_dynamics_ensemble.setter
    def molecular_dynamics_ensemble(self, value: MolecularDynamicsEnsemble) -> None:
        ...

    @property
    def thermobarostat(self) -> Thermobarostat | None:
        """(self) -> raspalib._raspalib.Thermobarostat | None"""
        ...

    @property
    def number_of_hybrid_mc_steps(self) -> int:
        """The number of molecular-dynamics steps of a hybrid Monte Carlo
        move."""
        ...

    @number_of_hybrid_mc_steps.setter
    def number_of_hybrid_mc_steps(self, value: int) -> None:
        ...

    @property
    def average_temperature(self) -> PropertyTemperature:
        """(self) -> raspalib._raspalib.PropertyTemperature"""
        ...

    @property
    def average_translational_temperature(self) -> PropertyTemperature:
        """(self) -> raspalib._raspalib.PropertyTemperature"""
        ...

    @property
    def average_rotational_temperature(self) -> PropertyTemperature:
        """(self) -> raspalib._raspalib.PropertyTemperature"""
        ...

    @property
    def average_simulation_box(self) -> PropertySimulationBox:
        """(self) -> raspalib._raspalib.PropertySimulationBox"""
        ...

    def set_thermobarostat(self, thermobarostat: Thermobarostat | None = None) -> None:
        """Build the barostat of the system from the given settings, with the
        temperature, the pressure and the time step the system has at this
        moment; ``None`` removes it. The ensembles with a barostat are
        thermostatted as well: see ``set_thermostat``."""
        ...

    def compute_total_energies(self) -> RunningEnergy:
        """Return the energies of the configuration as it is, in internal
        units of 10 J/mol."""
        ...

    @property
    def property_molecule_properties(self) -> PropertyMoleculeProperties | None:
        """The histograms of the molecules of the components that have
        ``molecule_properties_settings``; created by the driver."""
        ...

    @property
    def property_molecule_shape(self) -> PropertyMoleculeShape | None: ...

    @property
    def property_molecule_backbone(self) -> PropertyMoleculeBackbone | None: ...

    @property
    def property_end_to_end_acf(self) -> PropertyEndToEndAutoCorrelationFunction | None: ...

    @property
    def property_elastic_constants_fluctuation(self) -> PropertyElasticConstantsFluctuation | None: ...

    def set_property_elastic_constants_fluctuation(
        self, property: PropertyElasticConstantsFluctuation | None = None, sample_every: int = 100
    ) -> None:
        """Sample the elastic constants from the fluctuations of the stress,
        every ``sample_every`` steps of molecular dynamics at constant volume
        and temperature."""
        ...

    @property
    def property_radical_waiting_times(self) -> PropertyRadicalWaitingTimes | None: ...

    def set_property_radical_waiting_times(self, property: PropertyRadicalWaitingTimes | None = None) -> None: ...

    @property
    def write_lammps_data(self) -> WriteLammpsData | None: ...

    def set_write_lammps_data(self, property: WriteLammpsData | None = None) -> None: ...

    @property
    def reactions(self) -> collections.abc.Sequence[Reaction]:
        """The reactions of the system."""
        ...

    @property
    def number_of_blocks(self) -> int:
        """The number of blocks the averages of the system are kept in."""
        ...

    @property
    def mc_moves_probabilities(self) -> MCMoveProbabilities:
        """The probabilities of the moves of the system (volume, Gibbs
        volume, reactions)."""
        ...
