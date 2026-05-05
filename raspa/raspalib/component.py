import collections.abc
from raspalib.connectivity_table import *
from raspalib.intra_molecular_potentials import *
from raspalib.atom import *
from raspalib.force_field import *
from raspalib.mc_move_probabilities import *
from raspalib.mc_move_statistics import *
from raspalib.property_widom import *
from raspalib.property_lambda_probability_histogram import *

class Component():
    """
    Represent a simulation component in RASPA.

    A component groups force-field settings, molecular topology definitions,
    Monte Carlo move settings, and runtime statistics for a specific adsorbate
    or fluid species.
    """

    def __init__(
        self,
        force_field: ForceField,
        component_name: str,
        critical_temperature: float,
        critical_pressure: float,
        acentric_factor: float,
        defined_atoms: collections.abc.Sequence[Atom] = [],
        connectivity_table: ConnectivityTable = ConnectivityTable(),
        intra_molecular_potentials: IntraMolecularPotentials = IntraMolecularPotentials(),
        number_of_blocks: int = 5,
        number_of_lambda_bins: int = 21,
        particle_probabilities: MCMoveProbabilities = MCMoveProbabilities(),
        fugacity_coefficient: float | None = None,
        thermodynamic_integration: bool = False,
        blocking_pockets: collections.abc.Sequence[tuple[float, float, float, float]] = []
    ) -> None:
        ...
        """
        Initialize a :class:`Component`.

        Args:
            force_field: Force-field definition used for this component.
            component_name: Human-readable component name.
            critical_temperature: Critical temperature of the component.
            critical_pressure: Critical pressure of the component.
            acentric_factor: Acentric factor of the component.
            defined_atoms: Atom definitions used to build the component
                topology. Defaults to an empty sequence.
            connectivity_table: Bond/angle/dihedral connectivity information.
                Defaults to an empty :class:`ConnectivityTable`.
            intra_molecular_potentials: Intra-molecular potential settings.
                Defaults to an empty :class:`IntraMolecularPotentials`.
            number_of_blocks: Number of sampling blocks used for statistics.
                Defaults to ``5``.
            number_of_lambda_bins: Number of bins in lambda histograms.
                Defaults to ``21``.
            particle_probabilities: Monte Carlo move probabilities for this
                component. Defaults to a fresh
                :class:`MCMoveProbabilities` instance.
            fugacity_coefficient: Optional fugacity coefficient override.
                Defaults to ``None``.
            thermodynamic_integration: Enable thermodynamic integration
                bookkeeping. Defaults to ``False``.
            blocking_pockets: Sequence of blocking-pocket definitions as
                ``(x, y, z, radius)`` tuples. Defaults to an empty sequence.
        """

    @property
    def blocking_pockets(self) -> collections.abc.Sequence[tuple[float, float, float, float]]:
        """Return blocking-pocket definitions as ``(x, y, z, radius)`` tuples."""
        ...

    @blocking_pockets.setter
    def blocking_pockets(self, arg0: collections.abc.Sequence[tuple[float, float, float, float]]) -> None:
        """Set blocking-pocket definitions as ``(x, y, z, radius)`` tuples."""
        ...

    @property
    def average_rosenbluth_weights(self) -> PropertyWidom:
        """Return Widom insertion Rosenbluth-weight statistics."""
        ...

    @property
    def average_gibbs_rosenbluth_weights(self) -> PropertyWidom:
        """Return Gibbs-Widom insertion Rosenbluth-weight statistics."""
        ...

    @property
    def mc_moves_probabilities(self) -> MCMoveProbabilities:
        """Return configured Monte Carlo move probabilities."""
        ...

    @property
    def mc_moves_statistics(self) -> MCMoveStatistics:
        """Return Monte Carlo move acceptance/statistics counters."""
        ...

    @property
    def lambda_histogram(self) -> PropertyLambdaProbabilityHistogram:
        """Return lambda-probability histogram statistics."""
        ...

    @property
    def average_dudlambda(self) -> tuple[collections.abc.Sequence[tuple[float, float, float]], collections.abc.Sequence[tuple[float, float, float]]]:
        """Return average dU/dlambda data used in thermodynamic integration."""
        ...


