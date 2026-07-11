import collections.abc

from raspalib.system import *

class MonteCarlo():
    """
    Configure and run Monte Carlo workflows in RASPA.

    This class stores run-control settings and orchestrates initialization,
    equilibration, and production stages for one or more simulation systems.
    """

    def __init__(
        self,
        number_of_production_cycles: int = 0,
        number_of_pre_initialization_cycles: int = 0,
        number_of_initialization_cycles: int = 0,
        number_of_equilibration_cycles: int = 0,
        printEvery: int = 5000,
        write_binary_restart_every: int = 5000,
        rescale_wang_landau_every: int = 5000,
        optimize_mc_moves_every: int = 5000,
        systems: collections.abc.Sequence[System] = [],
        random_seed: int | None = None,
        number_of_blocks: int = 5,
        output_to_files: bool = False
    ) -> None:
        ...
        """
        Initialize a :class:`MonteCarlo` simulation controller.

        Args:
            number_of_production_cycles: Number of production cycles. Defaults to ``0``.
            number_of_initialization_cycles: Number of initialization cycles.
                Defaults to ``0``.
            number_of_equilibration_cycles: Number of equilibration cycles.
                Defaults to ``0``.
            printEvery: Reporting interval (in cycles). Defaults to ``5000``.
            write_binary_restart_every: Interval for writing binary restart
                files. Defaults to ``5000``.
            rescale_wang_landau_every: Interval for Wang-Landau rescaling
                operations. Defaults to ``5000``.
            optimize_mc_moves_every: Interval for Monte Carlo move-probability
                optimization. Defaults to ``5000``.
            systems: Sequence of simulation systems to run. Defaults to an
                empty sequence.
            random_seed: Optional random seed. If ``None``, the backend selects
                a seed.
            number_of_blocks: Number of blocks used for block averaging.
                Defaults to ``5``.
            output_to_files: Whether output should be written to files.
                Defaults to ``False``.
            number_of_pre_initialization_cycles: Number of pre-initialization
                cycles, run with only translation, rotation, reinsertion and
                partial-reinsertion moves. Defaults to ``0``.
        """
    def setup(self) -> None:
        """Allocate and prepare simulation resources before running."""
        ...
    def tear_down(self) -> None:
        """Release simulation resources after completion."""
        ...
    def pre_initialize(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        """Run the pre-initialization stage.

        Args:
            call_back_function: Optional callback invoked during execution.
        """
        ...
    def equilibrate(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        """Run the equilibration stage.

        Args:
            call_back_function: Optional callback invoked during execution.
        """
        ...
    def initialize(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        """Run the initialization stage.

        Args:
            call_back_function: Optional callback invoked during execution.
        """
        ...
    def production(self, call_back_function: collections.abc.Callable[[], None] | None = None) -> None:
        """Run the production stage.

        Args:
            call_back_function: Optional callback invoked during execution.
        """
        ...
    def run(self) -> None:
        """Execute the full configured Monte Carlo workflow."""
        ...

    @property
    def systems(self) -> collections.abc.Sequence[System]:
        """Return the configured list of simulation systems."""
        ...


