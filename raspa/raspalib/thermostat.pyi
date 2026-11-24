class Thermostat:
    """
    Define thermostat integration settings for RASPA simulations.

    The thermostat parameters configure chain length and integrator detail for
    temperature control during molecular dynamics.
    """
    def __init__(self,
                 thermostat_chain_length: int = 5,
                 number_of_yoshida_suzuki_steps: int = 5,
                 time_scale_parameter: float = 0.15) -> None:
        ...
        """
        Initialize a :class:`Thermostat`.

        Args:
            thermostat_chain_length: Length of the thermostat chain.
            number_of_yoshida_suzuki_steps: Number of Yoshida-Suzuki integration
                steps.
            time_scale_parameter: Thermostat time-scale parameter.
        """

    @property
    def langevin(self) -> bool:
        """Whether the temperature is held by a Langevin thermostat (friction
        and random kicks on every thermostatted velocity) instead of the
        Nose-Hoover chain. Default False."""
        ...

    @langevin.setter
    def langevin(self, value: bool) -> None:
        ...

    @property
    def langevin_friction(self) -> float:
        """The friction coefficient of the Langevin thermostat [1/ps].
        Default 1.0."""
        ...

    @langevin_friction.setter
    def langevin_friction(self, value: float) -> None:
        ...
