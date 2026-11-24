class WriteLammpsData():
    """Write the configuration as LAMMPS data files (``lammps/s<system>.data``
    and its tables) every ``sample_every`` cycles."""

    def __init__(self, system_id: int = 0, sample_every: int = 10) -> None:
        ...

    @property
    def system_id(self) -> int: ...
    @property
    def sample_every(self) -> int: ...
