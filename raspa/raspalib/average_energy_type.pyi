class AverageEnergyType():
    """
    Container for average interaction-energy contributions in RASPA.

    This type groups the total average energy and its primary decomposed terms
    (Van der Waals, Coulomb, and polarization) for convenient access from the
    Python API.
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
        Initialize an :class:`AverageEnergyType`.

        Args:
            total_energy: Total average energy contribution.
            van_der_waals_energy: Average Van der Waals energy contribution.
            coulomb_energy: Average Coulomb (electrostatic) energy contribution.
            polarization_energy: Average polarization energy contribution.
        """
    @property
    def coulomb_energy(self) -> float:
        """Return the average Coulomb (electrostatic) energy contribution."""
        ...
    @property
    def van_der_waals_energy(self) -> float:
        """Return the average Van der Waals energy contribution."""
        ...
    @property
    def polarization_energy(self) -> float:
        """Return the average polarization energy contribution."""
        ...
    @property
    def total_energy(self) -> float:
        """Return the total average energy contribution."""
        ...
