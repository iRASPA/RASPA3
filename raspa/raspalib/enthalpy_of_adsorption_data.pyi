import collections.abc

class EnthalpyOfAdsorptionData():
    """
    Store enthalpy-of-adsorption data reported by RASPA.

    The object provides access to sampled or averaged enthalpy values produced
    during simulation analysis.
    """

    def __init__(self, number_of_components: int) -> None:
        ...
        """
        Initialize an :class:`EnthalpyOfAdsorptionData` with one (zero) value per component.
        """

    def __getitem__(self, index: int) -> float:
        """Return the value of component ``index``; raises ``IndexError`` past the last component."""
        ...

    def __len__(self) -> int:
        """Return the number of components."""
        ...

    @property
    def values(self) -> collections.abc.Sequence[float]:
        """Return the enthalpy-of-adsorption values."""
        ...

