class RunningEnergy:
    """
    Provide per-sample running energy terms reported by RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize a :class:`RunningEnergy`.
        """
    def conserved_energy(self) -> float:
        """Return the conserved total energy."""
        ...
    def potential_energy(self) -> float:
        """Return the potential-energy contribution."""
        ...
    def kinetic_energy(self) -> float:
        """Return the total kinetic-energy contribution."""
        ...
    def translational_kinetic_energy(self) -> float:
        """Return the translational kinetic-energy contribution."""
        ...
    def rotational_kinetic_energy(self) -> float:
        """Return the rotational kinetic-energy contribution."""
        ...
    def thermostat_energy(self) -> float:
        """Return the thermostat energy contribution."""
        ...
    def coulomb_energy(self) -> float:
        """Return the Coulomb (electrostatic) energy contribution."""
        ...
    def van_der_waals_energy(self) -> float:
        """Return the Van der Waals energy contribution."""
        ...

