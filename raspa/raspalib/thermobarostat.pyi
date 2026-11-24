import enum


class MolecularDynamicsEnsemble(enum.IntEnum):
    """The ensemble of a molecular-dynamics simulation."""
    NVE: int = 0
    NVT: int = 1
    NPT: int = 2
    NPTPR: int = 3
    MUVT: int = 4
    MUPT: int = 5
    MUPTPR: int = 6


class BarostatCoupling(enum.IntEnum):
    """Whether the barostat scales the centers of mass of the molecules or the atoms."""
    MOLECULAR: int = 0
    ATOMIC: int = 1


class CellType(enum.IntEnum):
    """Which degrees of freedom of the cell change."""
    FIXED: int = 0
    ISOTROPIC: int = 1
    ANISOTROPIC: int = 2
    MONOCLINIC: int = 3
    REGULAR: int = 4
    REGULAR_UPPER_TRIANGLE: int = 5
    MONOCLINIC_UPPER_TRIANGLE: int = 6


class MonoclinicAngleType(enum.IntEnum):
    """The angle that changes in a monoclinic cell."""
    ALPHA: int = 0
    BETA: int = 1
    GAMMA: int = 2


class Thermobarostat:
    """
    Define barostat settings for molecular dynamics.

    ``System.set_thermobarostat`` builds the barostat of a system from these
    settings, with the temperature, the pressure and the time step of the
    system.
    """
    def __init__(self,
                 ensemble: MolecularDynamicsEnsemble,
                 barostat_chain_length: int = 5,
                 number_of_yoshida_suzuki_steps: int = 5,
                 time_scale_parameter: float = 1.0,
                 number_of_respa_steps: int = 5,
                 coupling: BarostatCoupling = BarostatCoupling.MOLECULAR,
                 cell_type: CellType | None = None,
                 monoclinic_angle: MonoclinicAngleType = MonoclinicAngleType.BETA) -> None:
        """
        Initialize a :class:`Thermobarostat`.

        Args:
            ensemble: An ensemble with a barostat: NPT, NPTPR, MUPT or MUPTPR.
            barostat_chain_length: Length of the barostat chain.
            number_of_yoshida_suzuki_steps: Number of Yoshida-Suzuki integration
                steps.
            time_scale_parameter: Barostat time-scale parameter [ps].
            number_of_respa_steps: Number of RESPA steps.
            coupling: Molecular or atomic scaling.
            cell_type: The cell changes of a flexible barostat (NPTPR, MUPTPR);
                ``None`` changes the whole cell. An isotropic barostat scales
                the volume.
            monoclinic_angle: The angle that changes in a monoclinic cell.
        """
        ...

    @property
    def ensemble(self) -> MolecularDynamicsEnsemble: ...
    @property
    def cell_type(self) -> CellType: ...
    @property
    def coupling(self) -> BarostatCoupling: ...
    @property
    def barostat_chain_length(self) -> int: ...
    @property
    def number_of_yoshida_suzuki_steps(self) -> int: ...
    @property
    def number_of_respa_steps(self) -> int: ...
    @property
    def time_scale_parameter(self) -> float: ...
