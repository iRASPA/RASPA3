class EnergyFactor():
    """
    Store an energy term and its thermodynamic-integration derivative.

    This lightweight container is used by RASPA to pass both an energy value
    and the corresponding ``dU/dlambda`` contribution together.
    """

    def __init__(
        self,
        energy: float,
        dudlambda: float
    ) -> None:
        ...
        """
        Initialize an :class:`EnergyFactor`.

        Args:
            energy: Energy contribution.
            dudlambda: Derivative of energy with respect to lambda
                (``dU/dlambda``).
        """

    @property
    def energy(self) -> float:
        """Return the stored energy contribution."""
        ...

    @property
    def dudlambda(self) -> float:
        """Return the stored ``dU/dlambda`` contribution."""
        ...

