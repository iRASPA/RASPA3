class WriteLammpsData():
    """Write the configuration as LAMMPS data files (``lammps/s<system>.data``
    and its tables) every ``sample_every`` cycles.

    Without ``system_id`` the files go by the index the system has in its
    driver; the files are then made when the driver is constructed."""

    def __init__(self, system_id: int | None = None, sample_every: int = 10) -> None:
        ...

    @property
    def system_id(self) -> int | None:
        """The index in the names of the files; None until a driver gives the
        index of the system."""
        ...
    @property
    def sample_every(self) -> int: ...
