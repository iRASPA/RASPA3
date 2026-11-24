import typing


class PropertyElasticConstantsFluctuation():
    """
    The elastic constants from the fluctuations of the stress, in molecular
    dynamics at constant volume and temperature.
    """

    def __init__(self, number_of_blocks: int = 5) -> None:
        ...

    def result(self, eigenvalue_tolerance: float = 1.0e-8) -> dict[str, typing.Any]:
        """The averages as a dictionary: the stress, the Born, kinetic and
        fluctuation terms and the stiffness in Voigt order [GPa], with the
        compliance, the moduli and the Poisson ratios derived from them when
        the stiffness can be inverted."""
        ...

    def write_averages_statistics(self, eigenvalue_tolerance: float = 1.0e-8) -> str: ...
