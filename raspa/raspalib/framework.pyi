import collections.abc

from raspalib.atom import *
from raspalib.simulation_box import *
from raspalib.force_field import *

class Framework():
    """
    Represent a (typically porous) host framework in RASPA.

    A framework combines structural information (unit-cell geometry, symmetry,
    atoms, and replication count) with the force-field settings used to model
    framework interactions.
    """

    def __init__(
        self,
        force_field: ForceField,
        component_name: str,
        simulation_box: SimulationBox,
        space_group_hall_number: int = None,
        defined_atoms: collections.abc.Sequence[Atom] = None,
        number_of_unit_cells: collections.abc.Sequence[int] = [1, 1, 1],
    ) -> None:
        ...
        """
        Initialize a :class:`Framework`.

        Args:
            force_field: Force-field definition used for the framework.
            component_name: Name/identifier of the framework component.
            simulation_box: Unit-cell box definition.
            space_group_hall_number: Optional space-group Hall number used for
                crystallographic symmetry. Defaults to ``None``.
            defined_atoms: Optional atom definitions for the framework.
                Defaults to ``None``.
            number_of_unit_cells: Number of replicated unit cells along each
                lattice direction. Defaults to ``[1, 1, 1]``.
        """
