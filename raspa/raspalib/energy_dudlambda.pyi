class EnergyDuDlambda():
    """
    Store an energy term and its thermodynamic-integration derivative.

    This lightweight container is used by RASPA to pass both an energy value
    and the corresponding ``dU/dlambda`` contribution together.
    """

    # returned by the library only; it has no constructor

    @property
    def energy(self) -> float:
        """Return the stored energy contribution."""
        ...

    @property
    def dudlambda(self) -> float:
        """Return the stored ``dU/dlambda`` contribution."""
        ...
