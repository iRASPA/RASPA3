import enum
import collections.abc

class PseudoAtom():
    """
    A class representing a pseudo atom in RASPA.
    """

    def __init__(
        self,
        name: str = "C",
        framework_type: bool = False,
        mass: float = 1.0,
        charge: float = 0.0,
        polarizability: float = 0.0,
        atomic_number: int = 8,
        print_to_pdb: bool = True,
        source: str = "-"
    ) -> None:
        ...
        """
        Initialize the PseudoAtom object with provided parameters.

        Args:
            name (str): The name of the pseudo atom. Default is "C".
            frameworkType (bool): Whether it is a framework-atom or not
            mass (float): The mass of the pseudo atom. Default is 1.0.
            charge (float): The charge of the pseudo atom. Default is 0.0.
            polarizability (float): The polarizability of the pseudo atom. Default is 0.0.
            atomicNumber (int): The atomic number of the pseudo atom. Default is 8.
            printToPDB (bool): Whether to print to PDB. Default is False.
            source (str): The source of the pseudo atom. Default is "-".
        """


class VDWParameters():
    """
    A class representing Van der Waals parameters in RASPA.
    """

    def __init__(self, epsilon: float, sigma: float) -> None:
        ...
        """
        Initialize the VDWParameter object with provided parameters.

        Args:
            epsilon (float): The epsilon parameter.
            sigma (float): The sigma parameter.
        """

class Atom():
    """
    A class representing an atom in RASPA.
    """

    def __init__(
        self,
        position: tuple[float, float, float],
        charge: float,
        scaling: float = 1.0,
        molecule_id: int = 0,
        type: int = 0,
        component_id: int = 0,
        group_id: bool = 0,
        is_fractional: bool = 0
    ) -> None:
        ...
        """
        Initialize the Atom object with provided parameters.

        Args:
            position (tuple[float, float, float]): The position of the atom as a 3-element tuple.
            charge (float): The charge of the atom.
            scaling (float, optional): The lambda value of the atom. Default is 1.0.
            molecule_id (int, optional): The molecule ID of the atom. Default is 0.
            type (int, optional): The type of the atom. Default is 0.
            component_id (int, optional): The component ID of the atom. Default is 0.
            group_id (bool, optional): The group ID of the atom. Default is false.
            is_fractional (bool, optional): The isFractional ID of the atom. Default is false.
        """


class ForceField():
    """
    A class representing a force field in RASPA.
    """

    class MixingRule:
        """
        Members:

          Lorentz_Berthelot
        """
        LORENTZ_BERTHELOT: int = 0
        JORGENSEN: int = 1

    def __init__(
        self,
        pseudo_atoms: collections.abc.Sequence[PseudoAtom] = None,
        parameters: collections.abc.Sequence[VDWParameters] = None,
        mixing_rule: MixingRule = MixingRule.Lorentz_Berthelot,
        cutoff_framework_vdw: float = 12.0,
        cutoff_molecule_vdw: float = 12.0,
        cutoff_coulomb: float = 12.0,
        shifted: bool = False,
        tail_corrections: bool = False,
        use_charge=True,
    ) -> None:
        """
        Initialize the ForceField object with provided parameters.

        Args:
            pseudo_atoms (Sequence[PseudoAtom], optional): A list of pseudo atoms. Default is None.
            parameters (Sequence[VDWParameter], optional): A list of Van der Waals parameters. Default is None.
            mixing_rule (enum(MixingRule), optional): The mixing rule. Default is "Lorentz_Berthelot".
            cutoff_framework_vdw (float, optional): The framework-molecule Van der Waals cut-off distance. Default is 12.0.
            cutoff_molecule_vdw (float, optional): The molecule-molecule Van der Waaks cut-off distance. Default is 12.0.
            cutoff_coulomb (float, optional): The cut-off distance for the real-part of the Ewald-summation . Default is 12.0.
            shifted (bool, optional): Whether to use shifted potential. Default is False.
            tail_corrections (bool, optional): Whether to apply tail corrections. Default is False.
            use_charge (bool, optional): Whether to compute electrostatics or not.
        """

    @property
    def pseudo_atoms(self) -> collections.abc.Sequence[PseudoAtom]:
        ...

    @property
    def vdw_parameters(self) -> collections.abc.Sequence[VDWParameters]:
        ...


class MCMoveProbabilities():
    """
    A class representing all move probabilities in RASPA.

    """

    def __init__(
        self,
        translation_probability: float = 0.0, 
        random_translation_probability: float = 0.0,
        rotation_probability: float = 0.0, 
        random_rotation_probability: float = 0.0,
        volume_change_probability: float = 0.0, 
        reinsertion_cbmc_probability: float = 0.0,
        partial_reinsertion_cbmc_probability: float = 0.0, 
        identity_change_probability: float = 0.0,
        swap_probability: float = 0.0, 
        swap_cbmc_probability: float = 0.0, 
        swap_cfcmc_probability: float = 0.0,
        swap_cbcfcmc_probability: float = 0.0, 
        gibbs_volume_change_probability: float = 0.0,
        gibbs_swap_cbmc_probability: float = 0.0,
        gibbs_swap_cfcmc_probability: float = 0.0,
        widom_probability: float = 0.0,
        widom_cfcmc_probability: float = 0.0,
        widom_cbcfcmc_probability: float = 0.0,
        parallel_tempering_probability: float = 0.0,
        hybrid_mc_probability: float = 0.0
    ) -> None:
        ...
        """
        Initialize a particle mc moves object that holds all probabilities for moves. It will be normalized after init.

        Args:
            probability_translation_move (float, optional): _description_. Defaults to 0.0.
            probability_random_translation_move (float, optional): _description_. Defaults to 0.0.
            probability_rotation_move (float, optional): _description_. Defaults to 0.0.
            probability_random_rotation_move (float, optional): _description_. Defaults to 0.0.
            probability_volume_move (float, optional): _description_. Defaults to 0.0.
            probability_reinsertion_move_cbmc (float, optional): _description_. Defaults to 0.0.
            probability_identity_change_move_cbmc (float, optional): _description_. Defaults to 0.0.
            probability_swap_move (float, optional): _description_. Defaults to 0.0.
            probability_swap_move_cbmc (float, optional): _description_. Defaults to 0.0.
            probability_swap_move_cfcmc (float, optional): _description_. Defaults to 0.0.
            probability_swap_move_cfcmc_cbmc (float, optional): _description_. Defaults to 0.0.
            probability_gibbs_volume_move (float, optional): _description_. Defaults to 0.0.
            probability_gibbs_swap_move_cbmc (float, optional): _description_. Defaults to 0.0.
            probability_gibbs_swap_move_cfcmc (float, optional): _description_. Defaults to 0.0.
            probability_gibbs_swap_move_cfcmc_cbmc (float, optional): _description_. Defaults to 0.0.
            probability_widom_move (float, optional): _description_. Defaults to 0.0.
            probability_widom_move_cfcmc (float, optional): _description_. Defaults to 0.0.
            probability_widom_move_cfcmc_cbmc (float, optional): _description_. Defaults to 0.0.
            probability_parallel_tempering_swap (float, optional): _description_. Defaults to 0.0.
        """

class ConnectivityTable():
    """
    A class representing a component in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the ConnectivityTable object.
        """

class IntraMolecularPotentials():
    """
    A class representing the intra-molecular potentials of a component in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the IntraMolecularPotentials object with provided parameters.
        """

class Move:
    """
    This is my move class.
    """
    def __init__(self) -> None:
        ...

    class Types(enum.IntEnum):
        """
        An enumeration for the supported Monte Carlo moves.

        - TRANSLATION: A randomly picked particle will be translated [0...maxChange) randomly in x, y, or z.
        - RANDOM_TRANSLATION:  A randomly picked particle will be translated [0...box-length) randomly in x, y, or z.
        - ROTATION:
        - RANDOM_ROTATION:
        - VOLUME_CHANGE:
        - REINSERTION_CBMC:
        - PARTIAL_REINSERTION_CBMC:
        - IDENTITY_CHANGE_CBMC:
        - SWAP:
        - SWAP_CBMC:
        - SWAP_CFCMC:
        - SWAP_CBCFCMC:
        - GIBBS_VOLUME:
        - GIBBS_SWAP_CBMC:
        - GIBBS_SWAP_CFCMC:
        - WIDOM:
        - WIDOM_CFCMC:
        - WIDOM_CBCFCMC:
        - PARALELL_TEMPERING:
        - HYBRID_MC:
        """
        def __init__(self) -> None:
            ...

        TRANSLATION = 0
        RANDOM_TRANSLATION = 1
        ROTATION = 2
        RANDOM_ROTATION = 3
        VOLUME_CHANGE = 4
        REINSERTION_CBMC = 5
        PARTIAL_REINSERTION_CBMC = 6
        IDENTITY_CHANGE_CBMC = 7
        SWAP = 8
        SWAP_CBMC = 9
        SWAP_CFCMC = 10
        SWAP_CBCFCMC = 11
        GIBBS_VOLUME = 12
        GIBBS_SWAP_CBMC = 13
        GIBBS_SWAP_CFCMC = 14
        WIDOM = 15
        WIDOM_CFCMC = 16
        WIDOM_CBCFCMC = 17
        PARALLEL_TEMPERING = 18
        HYBRID_MC = 19

class MoveStatisticsDouble():
    """
    A class representing a move-statistics<double> in RASPA.
    """
    def __init__(
        self
    ) -> None:
        ...
        """
        Initialize the MoveStatisticsDouble object.
        """

    @property
    def accepted(self) -> float:
        ...
    @property
    def all_counts(self) -> int:
        ...
    @property
    def constructed(self) -> float:
        ...
    @property
    def counts(self) -> float:
        ...
    @property
    def lower_limit(self) -> float:
        ...
    @lowerLimit.setter
    def lower_limit(self, arg0: float) -> None:
        ...
    @property
    def max_change(self) -> float:
        ...
    @maxChange.setter
    def max_change(self, arg0: float) -> None:
        ...
    @property
    def target_acceptance(self) -> float:
        ...
    @target_acceptance.setter
    def target_acceptance(self, arg0: float) -> None:
        ...
    @property
    def total_accepted(self) -> float:
        ...
    @property
    def total_constructed(self) -> float:
        ...
    @property
    def total_counts(self) -> float:
        ...
    @property
    def upper_limit(self) -> float:
        ...
    @upperLimit.setter
    def upper_limit(self, arg0: float) -> None:
        ...

class MoveStatisticsDouble3():
    """
    A class representing a move-statistics<double3> in RASPA.
    """
    def __init__(
        self
    ) -> None:
        ...
        """
        Initialize the MoveStatisticsDouble3 object.
        """

    @property
    def accepted(self) -> tuple[float, float, float]:
        ...
    @property
    def all_counts(self) -> int:
        ...
    @property
    def constructed(self) -> tuple[float, float, float]:
        ...
    @property
    def counts(self) -> tuple[float, float, float]:
        ...
    @property
    def lower_limit(self) -> tuple[float, float, float]:
        ...
    @lowerLimit.setter
    def lower_limit(self, arg0: tuple[float, float, float]) -> None:
        ...
    @property
    def max_change(self) -> tuple[float, float, float]:
        ...
    @maxChange.setter
    def max_change(self, arg0: tuple[float, float, float]) -> None:
        ...
    @property
    def target_acceptance(self) -> tuple[float, float, float]:
        ...
    @targetAcceptance.setter
    def target_acceptance(self, arg0: tuple[float, float, float]) -> None:
        ...
    @property
    def total_accepted(self) -> tuple[float, float, float]:
        ...
    @property
    def total_constructed(self) -> tuple[float, float, float]:
        ...
    @property
    def total_counts(self) -> tuple[float, float, float]:
        ...
    @property
    def upper_limit(self) -> tuple[float, float, float]:
        ...
    @upperLimit.setter
    def upper_limit(self, arg0: tuple[float, float, float]) -> None:
        ...

class MCMoveStatistics:
    """
    A class representing a the move statistics in RASPA.
    """
    def __getitem__(self, arg0: Move.Types) -> MoveStatisticsDouble | MoveStatisticsDouble3:
        ...

class WidomData:
    """
    A class representing a Widom-data property in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the Widom-data object.
        """
    @property
    def total(self) -> float:
        ...
    @property
    def excess(self) -> float:
        ...
    @property
    def ideal_gas(self) -> float:
        ...

class LoadingData():
    """
    A class representing loading-information in RASPA.
    """

    def __init__(
        self,
        number_of_molecules: collections.abc.Sequence[int],
        number_densities: collections.abc.Sequence[float],
        inverse_number_densities: collections.abc.Sequence[float]
    ) -> None:
        ...
        """
        Initialize the LoadingData object with provided parameters.

        Args:
            number_of_molecules (Sequence[int]): The number of molecules for each component
            number_densities (Sequence[float]): The number density for each component
            inverse_number_densities (Sequence[float]): The inverse number-density for each component
        """

    @property
    def inverse_number_densities(self) -> collections.abc.Sequence[float]:
        ...

    @property
    def number_densities(self) -> collections.abc.Sequence[float]:
        ...

    @property
    def number_of_molecules(self) -> collections.abc.Sequence[float]:
        ...

class PressureData():
    """
    A class representing pressure-information in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the PressureData object with provided parameters.
        """

    @property
    def total_pressure(self) -> float:
        ...

    @property
    def excess_pressure(self) -> float:
        ...

    @property
    def ideal_gas_pressure(self) -> float:
        ...



class SampleMovie():
    """
    A class representing a movie in RASPA.
    """

    def __init__(
        self,
        system_id: int,
        sample_every: int = 1,
        restrict_to_box: bool = True
    ) -> None:
        ...
        """
        Initialize the SampleMovie  object with provided parameters.

        Args:
            system_id (int): The ID of the system.
            sample_every (int, optional): The sample frequency. Default: 1.
            restrict_to_box (bool, optional): Whether to resetrict particle to the box confinement. Default: true.
        """

class EnergyFactor():
    """
    A class representing energy-factor-information in RASPA.
    """

    def __init__(
        self,
        energy: float,
        dudlambda: float
    ) -> None:
        ...
        """
        Initialize the  object with provided parameters.

        Args:
            energy (float): The energy.
            dudlambda (float): The dUdlambda factor.
        """

    @property
    def energy(self) -> float:
        ...
        """
        Get the energy of the energy factor.

        Returns:
            float: The energy.
        """

    @property
    def dudlambda(self) -> float:
        ...
        """
        Get the dudlambda of the energy factor.

        Returns:
            float: The dUdlambda.
        """

class EnergyStatus():
    """
    A class representing energy-status-information in RASPA.
    """

    def __init__(
        self
    ) -> None:
        ...
        """
        Initialize the  object with provided parameters.

        Args:
        """

    @property
    def total_energy(self) -> EnergyFactor:
        ...
        """
        Get the total energy of the energy status.

        Returns:
            float: The total energy.
        """


class AverageEnergyType():
    """
    A class representing the various to the energies in RASPA.
    """

    def __init__(
        self,
        total_energy: float,
        van_der_waals_energy: float,
        coulomb_energy: float,
        polarization_energy: float
    ) -> None:
        ...
        """
        Initialize object with provided parameters.

        Args:
            total_energy (float): The total energy.
            van_der_waals_energy (float): The Van der Waals energy.
            coulomb_energy (float): The Coulomb energy.
            polarization_energy (float): The polarization energy.
        """
    @property
    def coulomb_energy(self) -> float:
        ...
    @property
    def van_der_waals_energy(self) -> float:
        ...
    @property
    def polarization_energy(self) -> float:
        ...
    @property
    def total_energy(self) -> float:
        ...

class PropertyEnergy():
    """
    A class representing the average energy property in RASPA.
    """

    def __init__(
        self,
        number_of_blocks: int,
        number_of_external_fields: int,
        number_of_frameworks: int,
        number_of_components: int
    ) -> None:
        ...
        """
        Initialize the PropertyEnergy object with provided parameters.

        Args:
            number_of_blocks (int): The number of blocks.
            number_of_external_fields (int): The number of external fields.
            number_of_frameworks (int): The number of frameworks.
            number_of_components (int): The number of components.
        """

    def result(self) -> tuple[EnergyStatus, EnergyStatus]:
        ...
        """
        Returns the computed data.
        """


class PropertyEnergyHistogram():
    """
    A class representing an energy-histogram property in RASPA.
    """

    def __init__(
        self,
        number_of_blocks: int,
        number_of_bins: int,
        value_range: tuple[float, float],
        sample_every: int = 1,
        write_every: int = 5000
    ) -> None:
        ...
        """
        Initialize the PropertyEnergyHistogram object with provided parameters.

        Args:
            number_of_blocks (int): The number of blocks.
            number_of_bins (int): The number of bins.
            value_range (int): The range of energy values to consider.
            sample_every (int): The sample frequency.
            write_every (int): The write frequency.
        """


class PropertyLoading():
    """
    A class representing a loading property in RASPA.
    """

    def __init__(
        self,
        number_of_blocks: int,
        number_of_components: int
    ) -> None:
        ...
        """
        Initialize the PropertyLoading object with provided parameters.

        Args:
            number_of_blocks (int): The number of blocks.
            number_of_components (int): The number of components.
        """

    def average_loading_number_of_molecules(self, arg0: int) -> tuple[float, float]:
        ...
    def result(self) -> tuple[LoadingData, LoadingData]:
        ...

class EnthalpyOfAdsorptionData():
    """
    A class representing pressure-information in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the PressureData object with provided parameters.
        """

    @property
    def values(self) -> collections.abc.Sequence[float]:
        ...

class PropertyEnthalpy():
    """
    A class representing a loading property in RASPA.
    """

    def __init__(
        self,
        number_of_blocks: int,
        number_of_components: int
    ) -> None:
        ...
        """
        Initialize the PropertyEnthalpy object with provided parameters.

        Args:
            number_of_blocks (int): The number of blocks.
            number_of_components (int): The number of components.
        """

    def result(self) -> tuple[EnthalpyOfAdsorptionData, EnthalpyOfAdsorptionData]:
        ...

class PropertyPressure():
    """
    A class representing the pressure-property property in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the PropertyPressure object.
        """

    def result(self) -> tuple[PressureData, PressuresData]:
        ...

class PropertyNumberOfMoleculesEvolution():
    """
    A class representing a number of molecules evolution property in RASPA.
    """

    def __init__(
        self,
        number_of_cycles: int,
        number_of_components: int,
        sample_every: int,
        write_every: int | None = None
    ) -> None:
        ...
        """
        Initialize the PropertyNumberOfMoleculesEvolution object with provided parameters.

        Args:
            number_of_cycles (int): The number of cycles.
            number_of_components (int): The number of components.
            sample_every (int): The sample frequency.
            write_every (int): The write frequency.
        """

    @property
    def result(self) -> collections.abc.Sequence[collections.abc.Sequence[int]]:
        ...


class PropertyVolumeEvolution():
    """
    A class representing a a volume-evolution property in RASPA.
    """

    def __init__(
        self,
        number_of_cycles: int,
        sample_every: int,
        write_every: int | None = None
    ) -> None:
        ...
        """
        Initialize the PropertyVolumeEvolution object with provided parameters.

        Args:
            number_of_cycles (int): The number of cycles.
            sample_every (int): The sample frequency.
            write_every (int): The write frequency.
        """

    @property
    def result(self) -> collections.abc.Sequence[float]:
        ...

class RunningEnergy:
    """
    A class representing a Widom-data property in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the Widom-data object.
        """
    @property
    def conserved_energy(self) -> float:
        ...
    @property
    def potential_energy(self) -> float:
        ...
    @property
    def kinetic_energy(self) -> float:
        ...
    @property
    def nose_hoover_energy(self) -> float:
        ...

class PropertyConservedEnergyEvolution():
    """
    A class representing a number of molecules evolution property in RASPA.
    """

    def __init__(
        self,
        number_of_cycles: int,
        sample_every: int,
        write_every: int | None = None
    ) -> None:
        ...
        """
        Initialize the PropertyNumberOfMoleculesEvolution object with provided parameters.

        Args:
            number_of_cycles (int): The number of cycles.
            sample_every (int): The sample frequency.
            write_every (int): The write frequency.
        """

    @property
    def result(self) -> collections.abc.Sequence[collections.abc.Sequence[int]]:
        ...


class PropertyDensityGrid():
    """
    A class representing a density-grid property in RASPA.
    """

    class Binning(enum.IntEnum):
        STANDARD: int = 0
        EQUITABLE: int = 1

    class Normalization(enum.IntEnum):
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
        normalization_type: PropertyDensityGrid.Normalization = Normalization.MAX
        binning_mode: PropertyDensityGrid.Binning = Binning.STANDARD
    ) -> None:
        ...
        """
        Initialize the PropertyDensityGrid object with provided parameters.

        Args:
            number_of_frameworks (int): The number of frameworks.
            number_of_components (int): The number of components.
            number_of_grid_points (tuple[int, int, int]): The number of grid points.
            sample_every (int): The sample frequency.
            density_grid_pseudo_atoms_list (Sequence[str]): List of pseudo-atoms.
            normalization_type (): The normalization type
            binning_mode (): The binning mode
        """

class PropertyLambdaProbabilityHistogram():
    """
    A class representing a lambda-histogram property in RASPA.
    """

    def __init__(
        self,
        number_of_blocks: int,
        number_of_sample_points: int
    ) -> None:
        ...
        """
        Initialize the PropertyEnergyHistogram object with provided parameters.

        Args:
            number_of_blocks (int): The number of blocks.
            number_of_sample_points (int): The number of sample points.
        """

    def normalized_average_probability_histogram(self) -> tuple[collections.abc.Sequence[float], collections.abc.Sequence[float]]:
        ...

    @property
    def bias_factor(self) -> collections.abc.Sequence[float]:
        ...

    @biasFactor.setter
    def bias_factor(self, arg0: collections.abc.Sequence[float]) -> None:
        ...

    @property
    def histogram(self) -> collections.abc.Sequence[float]:
        ...

    @property
    def occupancy_count(self) -> float:
        ...

    @property
    def occupancy_total(self) -> float:
        ...


class PropertyWidom:
    """
    A class representing a Widom-insertion property in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize the PropertyWidom object.
        """

    @property
    def chemical_potential_result(self, temperature: float) -> WidomData:
        ...
    """
        Get the chemical potentials.

        Returns:
            WidomData: The chemical potentials in units of Kelvin.
        """

    @property
    def fugacity_result(self, temperature: float) -> float:
        ...
    """
        Get the fugacity.

        Returns:
            float: The fugacity in units of Pascals.
        """

class PropertyConventionalRadialDistributionFunction:
    """
    A class representing a PropertyConventionalRadialDistributionFunction property in RASPA.
    """

    def __init__(self,
                 numberOf_blocks: int,
                 number_of_pseudo_atoms: int,
                 number_of_bins: int,
                 range: double,
                 sample_every: int,
                 write_every: int) -> None:
        ...
        """
        Initialize the PropertyWidom object.
        """

    @property
    def result(self) -> list[list[tuple[collections.abc.Sequence[double],
                                        collections.abc.Sequence[double],
                                        collections.abc.Sequence[double]]]]:
        ...


class PropertyMeanSquaredDisplacement:
    """
    A class representing a PropertyMeanSquaredDisplacement property in RASPA.
    """

    def __init__(self,
                 sample_every: int,
                 write_every: int | None) -> None:
        ...
        """
        Initialize the PropertyMeanSquaredDisplacement object.
        """

    @property
    def result(self) -> list[collections.abc.Sequence[double],
                             collections.abc.Sequence[double]):
        ...

class PropertyVelocityAutoCorrelationFunction:
    """
    A class representing a PropertyVelocityAutoCorrelationFunction property in RASPA.
    """

    def __init__(self,
                 number_of_components: int,
                 number_of_molecules_per_component: collections.abc.Sequence[int],
                 number_of_particles: int,
                 number_of_buffers_vacf: int,
                 buffer_length_vacf: int,
                 sample_every: int,
                 write_every: int) -> None:
        ...
        """
        Initialize the PropertyWidom object.
        """

    @property
    def result(self) -> list[collections.abc.Sequence[double],
                             collections.abc.Sequence[double]):
        ...

class Thermostat:
    """
    A class representing a Thermostat in RASPA.
    """
    def __init__(self,
                 thermostat_chain_length: int,
                 number_of_yoshida_suzuki_steps: int) -> None:
        ...
        """
        Initialize the Thermostat object.
        """

class Component():
    """
    A class representing a component in RASPA.
    """

    def __init__(
        self,
        force_field: ForceField,
        component_name: str,
        critical_temperature: float,
        critical_pressure: float,
        acentric_factor: float,
        defined_atoms: collections.abc.Sequence[Atom] = [],
        connectivity_table: ConnectivityTable = ConnectivityTable(),
        intra_molecular_potentials: IntraMolecularPotentials = IntraMolecularPotentials(),
        number_of_blocks: int = 5,
        number_of_lambda_bins: int = 21,
        particle_probabilities: MCMoveProbabilities = MCMoveProbabilities(),
        fugacity_coefficient: float | None = None,
        thermodynamic_integration: bool = False,
        blocking_pockets: collections.abc.Sequence[tuple[float, float, float, float]] = []
    ) -> None:
        ...
        """
        Initialize the Component object with provided parameters. 

        Args:
            force_field (ForceField): The force field to be used.
            component_name (str): The name of the component.
            critical_temperature (float, optional): The critical temperature. Default is None.
            critical_pressure (float, optional): The critical pressure. Default is None.
            acentric_factor (float, optional): The acentric factor. Default is None.
            defined_atoms (Sequence[Atom], optional): A list of defined atoms. Default is None.
            connectivity_table (ConnectivityTable, optional): The connectivity table of the atoms. Default is empty.
            intra_molecular_potentials (IntraMolecularPotentials, optional): The intra molecular potentials. Default is none.
            number_of_blocks (int, optional): The number of blocks for the simulation. Default is 5.
            number_of_lambda_bins (int, optional): The number of lambda bins. Default is 21.
            particle_probabilities (MCMoveProbabilities, optional): The particle move probabilities. Default is a new instance of MCMoveProbabilities
            fugacity_coefficient (float, optional): The fugacity coefficient. Default is None.
            thermodynamic_integration (bool, optional): Whether to use thermodynamic integration. Default is False.
            blocking_pockets (Sequence[tuple[float, float, float, float]]): List of blocking-pockets
        """

    @property
    def blocking_pockets(self) -> collections.abc.Sequence[tuple[float, float, float, float]]:
        ...

    @blockingPockets.setter
    def blocking_pockets(self, arg0: collections.abc.Sequence[tuple[float, float, float, float]]) -> None:
        ...

    @property
    def average_rosenbluth_weights(self) -> PropertyWidom:
        ...
    """
        Get the Widom-move data.

        Returns:
            PropertyWidom: The results of the Widom insertions.
        """

    @property
    def average_gibbs_rosenbluth_weights(self) -> PropertyWidom:
        ...
    """
        Get the Gibbs-Widom-move data.

        Returns:
            PropertyWidom: The results of the Gibbs Widom insertions.
        """

    @property
    def mc_moves_probabilities(self) -> MCMoveProbabilities:
        ...
    """
        Get the move-probabilities

        Returns:
            MCMoveProbabilities: The probabilties for each of the moves.
        """

    @property
    def mc_moves_statistics(self) -> MCMoveStatistics:
        ...
    """
        Get the move-statistics

        Returns:
            MCMoveStatistics: The statistics for each of the moves.
        """

    @property
    def lambda_histogram(self) -> PropertyLambdaProbabilityHistogram:
        ...
    """
        Get the lambda-histogram-statistics

        Returns:
            PropertyLambdaProbabilityHistogram: The lambda histogram.
        """

    @property
    def average_dudlambda(self) -> tuple[collections.abc.Sequence[tuple[float, float, float]], collections.abc.Sequence[tuple[float, float, float]]]:
        ...
    """ 
        Get the thermodynamic integration data.
        
        Returns:
            tuple[]: The measured thermodynamic integration data.
        """


class SimulationBox():
    """
    A class representing simulation box in RASPA.
    """

    class SimulationBoxType(enum.IntEnum):
        Rectangular: int = 0
        Triclinic: int = 1

    @typing.overload
    def __init__(self, a: float, b: float, c: float) -> None:
        ...
        """
        Initialize a particle mc moves object that holds all probabilities for moves. It will be normalized after init.

        Args:
            a (float): The a-length of the box cell.
            b (float): The b-length of the box cell.
            c (float): The c-length of the box cell.
        """

    @typing.overload
    def __init__(self, a: float, b: float, c: float. alpha: float, beta: float, gamma: float) -> None:
        ...
        """
        Initialize a particle mc moves object that holds all probabilities for moves. It will be normalized after init.

        Args:
            a (float): The a-length of the box cell.
            b (float): The b-length of the box cell.
            c (float): The c-length of the box cell.
            alpha (float): The alpha-angle in degrees of the box cell.
            beta (float): The beta-angle in degrees of the box cell.
            gamma (float): The gamma-angle in degrees of the box cell.
        """

    @property
    def volume(self) -> float
        ...
        """
        Get the volume.

        Returns:
            float: The volume.
        """


class Framework():
    """
    A class representing a framework in RASPA, managing the initialization and configuration
    of a simulation framework.
    """

    def __init__(
        self,
        force_field: ForceField,
        component_name: str,
        simulation_box: SimulationBox,
        space_group_hall_number: int = None,
        defined_atoms: collections.abc.Sequence[Atom] = None,
        number_of_unit_cells: collections.abc.Sequence[int] = [1, 1, 1],
    ) -> None:
        ...
        """
        Initialize the Framework object with provided parameters.

        Args:
            force_field (ForceField): The force field to be used.
            component_name (str): The name of the component.
            simulation_box (SimulationBox): The unit cell box. 
            space_group_hall_number (int, optional): The space group Hall number. Default is None.
            defined_atoms (Sequence[Atom], optional): A list of defined atoms. Default is None.
            number_of_unit_cells (Sequence[int], optional): The number of unit cells in each dimension. Default is [1, 1, 1].
        """


class System():
    """
    A class representing a system in RASPA.
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
        Initialize the System object with provided parameters.

        Args:
            force_field (ForceField): The force field to be used.
            simulation_box (SimulationBox, optional): The simulation box. Default is None.
            external_temperature (float): The temperature of the system.
            external_pressure (float | None): The pressure of the system. Default is None.
            helium_void_fraction (float): The helium void-fraction of the system.
            framework_components (Framework | None, optional): The framework component if present. Default is None.
            components (Sequence[Component]): A list of components in the system.
            initial_positions (Sequence[int]): A list of initial positions. Default is empty list.
            initial_number_of_molecules (Sequence[int]): A list of initial number of molecules for each component.
            number_of_blocks (int, optional): The number of blocks for the simulation. Default is 5.
            system_probabilities (MCMoveProbabilitiesSystem, optional): The move probabilities system. Default is a new instance of MCMoveProbabilitiesSystem.
        """

    @property
    def components(self) -> collections.abc.Sequence[Component]:
        ...
        """
        Get the list of components.

        Returns:
            Sequence[Component]: The list of components.
        """

    @property
    def number_of_molecules_per_component(self) -> collections.abc.Sequence[int]:
        ...
        """
        Get the list of components.

        Returns:
            Sequence[Component]: The list of components.
        """

    def framework_mass(self) -> float | None:
        ...

    @property
    def mc_moves_statistics(self) -> MCMoveStatistics:
        ...
    """ 
        Get the move-statistics
        
        Returns:
            MCMoveStatistics: The statistics for each of the moves.
        """
    
    @property
    def simulation_box(self) -> SimulationBox:
        ...

    @property
    def loadings(self) -> LoadingData:
        ...

    @property
    def average_loadings(self) -> PropertyLoading:
        ...

    @property
    def average_enthalpies_of_adsorption(self) -> PropertyEnthalpy:
        ...

    @property
    def sample_pdb_movie(self) -> SampleMovie | None:
        ...

    @property
    def average_energies(self) -> PropertyEnergy:
        ...
        """
        Get the average-energies property.

        Returns:
            PropertyEnergy: The average-energies property.
        """
    @property
    def average_pressure(self) -> PropertyPressure:
        ...
        """
        Get the average-pressure property.

        Returns:
            PropertyPressure: The average-pressure property.
        """

    @property
    def property_density_grid(self) -> PropertyDensityGrid | None:
        ...
        """
        Get the rdf property.

        Returns:
            PropertyDensityGrid | None: The rdf property.
        """

    @property
    def average_energy_histogram(self) -> PropertyEnergyHistogram | None:
        ...
        """
        Get the rdf property.

        Returns:
            PropertyEnergyHistogram | None: The rdf property.
        """

    @property
    def average_number_of_molecules_histogram(self) -> PropertyNumberOfMoleculesHistogram | None:
        ...
        """
        Get the rdf property.

        Returns:
            PropertyNumberOfMoleculesHistogram | None: The rdf property.
        """

    @property
    def property_conventional_rdf(self) -> PropertyConventionalRadialDistributionFunction | None:
        ...
        """
        Get the rdf property.

        Returns:
            PropertyConventionalRadialDistributionFunction | None: The rdf property.
        """

    @property
    def property_msd(self) -> PropertyMeanSquaredDisplacement | None:
        ...
        """
        Get the msd property.

        Returns:
            PropertyMeanSquaredDisplacement | None: The msd property.
        """

    @property
    def property_vacf(self) -> PropertyVelocityAutoCorrelationFunction | None:
        ...
        """
        Get the vacf property.

        Returns:
            PropertyVelocityAutoCorrelationFunction | None: The vacf property.
        """
    
    @property
    def property_number_of_molecules_evolution(self) -> PropertyNumberOfMoleculesEvolution | None:
        ...
        """
        Get the average-pressure property.

        Returns:
            PropertyNumberOfMoleculesEvolution: The average-pressure property.
        """
    
    @property
    def property_volume_evolution(self) -> PropertyVolumeEvolution | None:
        ...
        """
        Get the average-pressure property.

        Returns:
            PropertyVolumeEvolution: The average-pressure property.
        """

    @property
    def property_conserved_energy_evolution(self) -> PropertyConservedEnergyEvolution | None:
        ...
        """
        Get the average-pressure property.

        Returns:
            PropertyConservedEnergyEvolution: The average-pressure property.
        """

    @set_thermostat.setter
    def set_thermostat(self, rdf: Thermostat) -> None:
        ...

    @set_sample_pdb_movie.setter
    def set_sample_pdb_movie(self, sample_movie: SampleMovie) -> None:
        ...

    @set_average_energy_histogram.setter
    def set_average_energy_histogram(self, property: PropertyEnergyHistogram) -> None:
        ...

    @set_property_density_grid.setter
    def set_property_density_grid(self, property: PropertyDensityGrid) -> None:
        ...

    @set_property_number_of_molecules_evolution.setter
    def set_property_number_of_molecules_evolution(self, property: PropertyNumberOfMoleculesEvolution) -> None:
        ...

    @set_property_volume_evolution.setter
    def set_property_volume_evolution(self, property: PropertyVolumeEvolution) -> None:
        ...

    @set_property_conserved_energy_evolution.setter
    def set_property_conserved_energy_evolution(self, property: PropertyConservedEnergyEvolution) -> None:
        ...

    @set_property_conventional_rdf.setter
    def set_property_conventional_rdf(self, property: PropertyConventionalRadialDistributionFunction) -> None:
        ...

    @set_property_msd.setter
    def set_property_msd(self, property: PropertyMeanSquaredDisplacement) -> None:
        ...

    @set_property_vacf.setter
    def set_property_msd(self, property: PropertyVelocityAutoCorrelationFunction) -> None:
        ...

class MonteCarlo():
    """
    A class representing a Monte Carlo simulation in RASPA.
    """

    def __init__(
        self,
        number_of_cycles: int = 0,
        number_of_initialization_cycles: int = 0,
        number_of_equilibration_cycles: int = 0,
        printEvery: int = 5000,
        write_binary_restart_every: int = 5000,
        rescale_wang_landau_every: int = 5000,
        optimize_mc_moves_every: int = 5000,
        systems: collections.abc.Sequence[System] = [],
        random_seed: int | None = None,
        number_of_blocks: int = 5,
        output_to_files: bool = False
    ) -> None:
        ...
        """
        Initializes a Monte Carlo object.

        Args:
            number_of_cycles (int, optional): _description_. Defaults to 0.
            number_of_initialization_cycles (int, optional): _description_. Defaults to 0.
            number_of_equilibration_cycles (int, optional): _description_. Defaults to 0.
            number_of_equilibration_cycles (int, optional): _description_. Defaults to 0.
            print_every (int, optional): _description_. Defaults to 5000.
            write_binary_restart_every (int, optional): _description_. Defaults to 5000.
            rescale_wang_landau_every (int, optional): _description_. Defaults to 5000.
            optimize_mc_moves_every (int, optional): _description_. Defaults to 5000.
            systems (Sequence[System], optional): _description_. Defaults to None.
            random_seed (int, optional): _description_. Defaults to a random integer.
            number_of_blocks (int, optional): _description_. Defaults to 5.
        """
    def setup(self) -> None:
        ...
    def tear_down(self) -> None:
        ...
    def equilibrate(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...
    def initialize(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...
    def production(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...
    def run(self) -> None:
        ...

    @property
    def systems(self) -> collections.abc.Sequence[System]:
        ...
        """
        Get the list of systems.

        Returns:
            Sequence[System]: The list of systems.
        """

class EquationOfState():
    """
    A class representing a lambda-histogram property in RASPA.
    """

    class EquationOfStateType(enum.IntEnum):
        PENG_ROBINSON: int = 0
        PENG_ROBINSON_GASEM: int = 1
        SOAVE_REDLICH_KWONG: int = 2

    class MixingRules(enum.IntEnum):
        VAN_DER_WAALS: int = 0

    class FluidState(enum.IntEnum):
        UNKNOWN: int = 0
        SUPER_CRITICAL_FLUID: int = 1
        VAPOR: int = 2
        LIQUID: int = 3
        VAPOR_LIQUID: int = 4

    class FluidInput():
        """
        A class representing a lambda-histogram property in RASPA.
        """
        def __init__(
            self,
            critical_temperature: float,
            critical_pressure: float,
            acentric_factor: float,
            mol_fraction: float = 1.0,
            swappable: bool = True
        ) -> None:
            ...

    @property
    def critical_temperature(self) -> float:
        ...

    @property
    def critical_pressure(self) -> float:
        ...

    @property
    def acentric_factor(self) -> float:
        ...

    @property
    def mol_fraction(self) -> float:
        ...

    @property
    def swappable(self) -> bool:
        ...

    class FluidResult():
        """
        A class representing a lambda-histogram property in RASPA.
        """

        @property
        def compressibility(self) -> float:
            ...

        @property
        def fugacity_coefficient(self) -> float | None:
            ...

        @property
        def fluid_state(self) -> FluidState:
            ...

    def __init__(self) -> None:
        ...

    @staticmethod
    def computeFluidProperties(temperature: float,
                               pressure: float,
                               properties: collections.abc.Sequence[FluidInput],
                               type: EquationOfStateType,
                               mixing_rules: EquationOfStateMultiComponentMixingRules,
                               ) -> collections.abc.Sequence[EquationOfState.FluidResult]:
        ...

class MolecularDynamics():
    """
    A class representing a Molecular Dynamics simulation in RASPA.
    """

    def __init__(
        self,
        number_of_cycles: int = 0,
        number_of_initialization_cycles: int = 0,
        number_of_equilibration_cycles: int = 0,
        print_every: int = 5000,
        write_binary_restart_every: int = 5000,
        rescale_wang_landau_every: int = 5000,
        optimize_mc_moves_every: int = 5000,
        systems: collections.abc.Sequence[System] = [],
        random_seed: int | None = None,
        number_of_blocks: int = 5,
        output_to_files: bool = False
    ) -> None:
        ...
        """
        Initializes a Moelcular Dynamics object.

        Args:
            number_of_cycles (int, optional): _description_. Defaults to 0.
            number_of_initialization_cycles (int, optional): _description_. Defaults to 0.
            number_of_equilibration_cycles (int, optional): _description_. Defaults to 0.
            number_of_equilibration_cycles (int, optional): _description_. Defaults to 0.
            print_every (int, optional): _description_. Defaults to 5000.
            write_binary_restart_every (int, optional): _description_. Defaults to 5000.
            rescale_wang_landau_every (int, optional): _description_. Defaults to 5000.
            optimize_mc_moves_every (int, optional): _description_. Defaults to 5000.
            systems (Sequence[System], optional): _description_. Defaults to None.
            random_seed (int, optional): _description_. Defaults to a random integer.
            number_of_blocks (int, optional): _description_. Defaults to 5.
        """
    def setup(self) -> None:
        ...
    def tear_down(self) -> None:
        ...
    def equilibrate(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...
    def initialize(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...
    def production(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        ...

    @property
    def systems(self) -> collections.abc.Sequence[System]:
        ...
        """
        Get the list of systems.

        Returns:
            Sequence[System]: The list of systems.
        """


