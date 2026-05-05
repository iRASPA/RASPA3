class Thermostat:
    """
    Define thermostat integration settings for RASPA simulations.

    The thermostat parameters configure chain length and integrator detail for
    temperature control during molecular dynamics.
    """
    def __init__(self,
                 thermostat_chain_length: int,
                 number_of_yoshida_suzuki_steps: int,
                 time_scale_parameter: float) -> None:
        ...
        """
        Initialize a :class:`Thermostat`.

        Args:
            thermostat_chain_length: Length of the thermostat chain.
            number_of_yoshida_suzuki_steps: Number of Yoshida-Suzuki integration
                steps.
            time_scale_parameter: Thermostat time-scale parameter.
        """
