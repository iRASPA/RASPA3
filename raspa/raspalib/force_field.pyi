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
        SIXTH_POWER: int = 2

    class ChargeMethod(enum.IntEnum):
        """
        How the electrostatic interactions are computed.
        """
        EWALD: int = 0
        COULOMB: int = 1
        WOLF: int = 2
        DAMPED_SHIFTED_FORCE: int = 3
        MODIFIED_SHIFTED_FORCE: int = 4
        ZERO_DIPOLE: int = 5

    def __init__(
        self,
        pseudo_atoms: collections.abc.Sequence[PseudoAtom],
        parameters: collections.abc.Sequence[VDWParameters],
        mixing_rule: MixingRule = MixingRule.LORENTZ_BERTHELOT,
        cutoff_framework_vdw: float = 12.0,
        cutoff_molecule_vdw: float = 12.0,
        cutoff_coulomb: float = 12.0,
        shifted: bool = True,
        tail_corrections: bool = False,
        use_charge=True,
        binary_interactions: collections.abc.Sequence[tuple[str, str, VDWParameters]] = [],
        parameters_14: collections.abc.Sequence[VDWParameters | None] = [],
        binary_interactions_14: collections.abc.Sequence[tuple[str, str, VDWParameters]] = [],
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
            binary_interactions: Pairs of pseudo-atoms whose potential is given
                instead of mixed: the two names and the potential of the pair.
            parameters_14: Separate parameters for the intramolecular 1-4
                pairs, one entry for every pseudo-atom; ``None`` where the 1-4
                pairs have the regular parameters.
            binary_interactions_14: Pairs of pseudo-atoms whose 1-4 potential
                is given.
        """

    @property
    def pseudo_atoms(self) -> collections.abc.Sequence[PseudoAtom]:
        """Return the configured pseudo-atom definitions."""
        ...

    @property
    def vdw_parameters(self) -> collections.abc.Sequence[VDWParameters]:
        """Return the configured Van der Waals interaction parameters."""
        ...

    # settings that are properties: a setting is made before the force field is given to a component or a system

    @property
    def use_charge(self) -> bool:
        """(self) -> bool"""
        ...

    @use_charge.setter
    def use_charge(self, value: bool) -> None:
        ...

    @property
    def charge_method(self) -> 'ForceField.ChargeMethod':
        """How the electrostatic interactions are computed."""
        ...

    @charge_method.setter
    def charge_method(self, value: 'ForceField.ChargeMethod') -> None:
        ...

    @property
    def ewald_precision(self) -> float:
        """The relative precision the Ewald parameters are chosen for [-];
        setting it makes the choice automatic."""
        ...

    @ewald_precision.setter
    def ewald_precision(self, value: float) -> None:
        ...

    @property
    def ewald_alpha(self) -> float:
        """The Ewald damping parameter alpha [1/Angstrom]."""
        ...

    @property
    def number_of_wave_vectors(self) -> tuple[int, int, int]:
        """The number of wave vectors of the Ewald sum in each direction."""
        ...

    @property
    def automatic_ewald(self) -> bool:
        """Whether alpha and the wave vectors are chosen from
        'ewald_precision'."""
        ...

    @property
    def modified_shifted_force_beta(self) -> float:
        """The switching parameter beta of the modified shifted-force charge
        method [1/Angstrom]."""
        ...

    @modified_shifted_force_beta.setter
    def modified_shifted_force_beta(self, value: float) -> None:
        ...

    @property
    def omit_ewald_fourier(self) -> bool:
        """Leaves out the Fourier part of the Ewald sum."""
        ...

    @omit_ewald_fourier.setter
    def omit_ewald_fourier(self, value: bool) -> None:
        ...

    @property
    def compute_polarization(self) -> bool:
        """Whether the polarization energy is computed."""
        ...

    @compute_polarization.setter
    def compute_polarization(self, value: bool) -> None:
        ...

    @property
    def omit_inter_polarization(self) -> bool:
        """Leaves out the polarization between molecules."""
        ...

    @omit_inter_polarization.setter
    def omit_inter_polarization(self, value: bool) -> None:
        ...

    @property
    def omit_inter_interactions(self) -> bool:
        """Leaves out the interactions between molecules."""
        ...

    @omit_inter_interactions.setter
    def omit_inter_interactions(self, value: bool) -> None:
        ...

    @property
    def cutoff_framework_vdw(self) -> float:
        """(self) -> float"""
        ...

    @property
    def cutoff_molecule_vdw(self) -> float:
        """(self) -> float"""
        ...

    @property
    def cutoff_coulomb(self) -> float:
        """(self) -> float"""
        ...

    def set_ewald_parameters(self, alpha: float, number_of_wave_vectors: collections.abc.Sequence[int]) -> None:
        """Fix the Ewald damping parameter alpha [1/Angstrom] and the number
        of wave vectors, instead of choosing them from ``ewald_precision``."""
        ...
