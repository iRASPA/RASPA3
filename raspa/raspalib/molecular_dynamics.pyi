import typing
import collections.abc

from raspalib.system import *

class MolecularDynamics():
    """
    Configure and run molecular-dynamics workflows in RASPA.

    This class stores run-control settings and orchestrates initialization,
    equilibration, and production stages for one or more simulation systems.
    """

    @typing.overload
    def __init__(
        self,
        number_of_cycles: int,
        number_of_pre_initialization_cycles: int = 0,
        number_of_initialization_cycles: int = 0,
        number_of_equilibration_cycles: int = 0,
        print_every: int = 5000,
        write_binary_restart_every: int = 5000,
        rescale_wang_landau_every: int = 1000,
        optimize_mc_moves_every: int = 100,
        systems: collections.abc.Sequence[System] = [],
        random_seed: int | None = None,
        number_of_blocks: int | None = None,
        output_to_files: bool = False
    ) -> None:
        ...

    @typing.overload
    def __init__(
        self,
        number_of_production_cycles: int,
        number_of_pre_initialization_cycles: int = 0,
        number_of_initialization_cycles: int = 0,
        number_of_equilibration_cycles: int = 0,
        print_every: int = 5000,
        write_binary_restart_every: int = 5000,
        rescale_wang_landau_every: int = 1000,
        optimize_mc_moves_every: int = 100,
        systems: collections.abc.Sequence[System] = [],
        random_seed: int | None = None,
        number_of_blocks: int | None = None,
        output_to_files: bool = False
    ) -> None:
        ...
        """
        Initialize a :class:`MolecularDynamics` simulation controller.

        Args:
            number_of_cycles: The keyword of the 3.1.0 release for
                ``number_of_production_cycles``.
            number_of_production_cycles: Number of production cycles.
            number_of_pre_initialization_cycles: Number of pre-initialization cycles.
                Defaults to ``0``.
            number_of_initialization_cycles: Number of initialization cycles.
                Defaults to ``0``.
            number_of_equilibration_cycles: Number of equilibration cycles.
                Defaults to ``0``.
            print_every: Reporting interval (in cycles). Defaults to ``5000``.
            write_binary_restart_every: Interval for writing binary restart
                files. Defaults to ``5000``.
            rescale_wang_landau_every: Interval for Wang-Landau rescaling
                operations. Defaults to ``1000``.
            optimize_mc_moves_every: Interval for Monte Carlo move-probability
                optimization. Defaults to ``100``.
            systems: Sequence of simulation systems to run. Defaults to an
                empty sequence. The driver works on its own copy of them;
                when ``run`` or a phase returns, its systems are assigned
                back to the objects given here (not to an object that
                was given more than once). The properties that were
                given to the setters of those systems are brought up to
                date as well. An evolution property made without a number
                of cycles records the cycles of all phases of this driver.
            random_seed: Optional random seed. If ``None``, the backend selects
                a seed.
            number_of_blocks: Number of blocks used for block averaging.
                This is the one place to give it: the driver sizes its
                copies of the systems, their components and their sampled
                properties for it (statistics that were kept in another
                number of blocks start anew). Without it the driver takes
                the number of its systems, which then have to agree.
            output_to_files: Whether output should be written to files.
                Defaults to ``False``.
        """
    def setup(self) -> None:
        """Allocate and prepare simulation resources before running."""
        ...
    def tear_down(self) -> None:
        """Release simulation resources after completion."""
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

    @property
    def number_of_blocks(self) -> int:
        """The number of blocks the averages and their errors are computed over."""
        ...

    @property
    def systems(self) -> collections.abc.Sequence[System]:
        """Return the systems of the driver.

        These are the systems the driver works on. Inside a call-back
        they are the ones to read: the objects given to the constructor
        are brought up to date when the phase returns.
        """
        ...

    def restart_from_binary_file(self, file_name: str = "restart_data.bin") -> None:
        """Read the state of the driver from a binary restart file.

        ``run`` then continues where that run stopped. A restart file is
        written every ``write_binary_restart_every`` cycles when
        ``output_to_files`` is set.
        """
        ...
