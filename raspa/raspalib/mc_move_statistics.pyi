from raspalib.move import *
from raspalib.move_statistics_double import *
from raspalib.move_statistics_double3 import *

class MCMoveStatistics:
    """
    Mapping-style access to Monte Carlo move statistics in RASPA.

    Instances expose per-move acceptance/statistics records and are indexed by
    :class:`Move.Types`.
    """
    def __getitem__(self, arg0: Move.Types) -> MoveStatisticsDouble | MoveStatisticsDouble3:
        """Return statistics for the requested move type."""
        ...
