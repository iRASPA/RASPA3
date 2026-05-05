import enum
import collections.abc

from raspalib.pseudo_atom import *
from raspalib.vdw_parameters import *

class ForceField():
    """
    Define non-bonded force-field settings used by RASPA.

    This object combines pseudo-atom definitions, Van der Waals parameter
    tables, mixing-rule selection, and electrostatics/cutoff options.
    """

    class MixingRule(enum.IntEnum):
        """
        Mixing rule used to combine unlike Van der Waals parameters.
        """
        LORENTZ_BERTHELOT: int = 0
        JORGENSEN: int = 1

    def __init__(
        self,
        pseudo_atoms: collections.abc.Sequence[PseudoAtom] = None,
        parameters: collections.abc.Sequence[VDWParameters] = None,
        mixing_rule: MixingRule = MixingRule.LORENTZ_BERTHELOT,
        cutoff_framework_vdw: float = 12.0,
        cutoff_molecule_vdw: float = 12.0,
        cutoff_coulomb: float = 12.0,
        shifted: bool = False,
        tail_corrections: bool = False,
        use_charge=True,
    ) -> None:
        """
        Initialize a :class:`ForceField`.

        Args:
            pseudo_atoms: Sequence of pseudo-atom definitions. Defaults to
                ``None``.
            parameters: Sequence of Van der Waals parameters. Defaults to
                ``None``.
            mixing_rule: Mixing rule for unlike interactions. Defaults to
                :attr:`MixingRule.LORENTZ_BERTHELOT`.
            cutoff_framework_vdw: Framework-molecule Van der Waals cutoff
                distance. Defaults to ``12.0``.
            cutoff_molecule_vdw: Molecule-molecule Van der Waals cutoff
                distance. Defaults to ``12.0``.
            cutoff_coulomb: Real-space cutoff distance for Ewald electrostatic
                interactions. Defaults to ``12.0``.
            shifted: Whether to apply shifted non-bonded potentials. Defaults
                to ``False``.
            tail_corrections: Whether to apply long-range tail corrections.
                Defaults to ``False``.
            use_charge: Whether electrostatic interactions are enabled.
                Defaults to ``True``.
        """

    @property
    def pseudo_atoms(self) -> collections.abc.Sequence[PseudoAtom]:
        """Return the configured pseudo-atom definitions."""
        ...

    @property
    def vdw_parameters(self) -> collections.abc.Sequence[VDWParameters]:
        """Return the configured Van der Waals interaction parameters."""
        ...

