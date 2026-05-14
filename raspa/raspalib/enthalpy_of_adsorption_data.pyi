import collections.abc

class EnthalpyOfAdsorptionData():
    """
    Store enthalpy-of-adsorption data reported by RASPA.

    The object provides access to sampled or averaged enthalpy values produced
    during simulation analysis.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize an :class:`EnthalpyOfAdsorptionData`.
        """

    @property
    def values(self) -> collections.abc.Sequence[float]:
        """Return the enthalpy-of-adsorption values."""
        ...

