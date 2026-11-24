from raspalib.move import *

class MCMoveProbabilities():
    """
    Store Monte Carlo move-selection probabilities for a component or a system.

    This container holds the relative probabilities of the moves; they are
    normalized when a move is drawn. Every move of :class:`Move.Types` has a
    keyword: its name in lower case with ``_probability`` appended.
    """

    def __init__(
        self,
        *,
        translation_probability: float = 0.0,
        random_translation_probability: float = 0.0,
        rotation_probability: float = 0.0,
        random_rotation_probability: float = 0.0,
        translation_smart_mc_probability: float = 0.0,
        translation_smart_mc_all_probability: float = 0.0,
        rotation_smart_mc_probability: float = 0.0,
        rotation_smart_mc_all_probability: float = 0.0,
        translation_rotation_smart_mc_probability: float = 0.0,
        volume_change_probability: float = 0.0,
        anisotropic_volume_change_probability: float = 0.0,
        reinsertion_cbmc_probability: float = 0.0,
        partial_reinsertion_cbmc_probability: float = 0.0,
        identity_change_cbmc_probability: float = 0.0,
        identity_switch_cbmc_probability: float = 0.0,
        swap_probability: float = 0.0,
        swap_cbmc_probability: float = 0.0,
        swap_cfcmc_probability: float = 0.0,
        swap_cbcfcmc_probability: float = 0.0,
        gibbs_volume_probability: float = 0.0,
        gibbs_swap_cbmc_probability: float = 0.0,
        gibbs_conventional_cfcmc_probability: float = 0.0,
        gibbs_conventional_cbcfcmc_probability: float = 0.0,
        gibbs_swap_cfcmc_probability: float = 0.0,
        gibbs_swap_cbcfcmc_probability: float = 0.0,
        gibbs_identity_change_cbmc_probability: float = 0.0,
        widom_probability: float = 0.0,
        widom_cfcmc_probability: float = 0.0,
        widom_cbcfcmc_probability: float = 0.0,
        parallel_tempering_probability: float = 0.0,
        hybrid_mc_probability: float = 0.0,
        reaction_cbmc_probability: float = 0.0,
        reaction_conventional_cfcmc_probability: float = 0.0,
        reaction_conventional_cbcfcmc_probability: float = 0.0,
        reaction_cfcmc_probability: float = 0.0,
        reaction_cbcfcmc_probability: float = 0.0,
        pair_swap_probability: float = 0.0,
        pair_swap_cbmc_probability: float = 0.0,
        pair_swap_cfcmc_probability: float = 0.0,
        pair_swap_cbcfcmc_probability: float = 0.0,
        group_swap_probability: float = 0.0,
        group_swap_cbmc_probability: float = 0.0,
        group_swap_cfcmc_probability: float = 0.0,
        group_swap_cbcfcmc_probability: float = 0.0,
        tethered_proton_hop_probability: float = 0.0,
        pivot_probability: float = 0.0,
        crankshaft_probability: float = 0.0,
        reptation_probability: float = 0.0,
        concerted_rotation_probability: float = 0.0,
        bead_displacement_probability: float = 0.0,
        bead_flip_probability: float = 0.0,
        double_bridging_probability: float = 0.0,
        intramolecular_double_rebridging_probability: float = 0.0,
        cross_link_swap_probability: float = 0.0,
        cross_link_formation_scission_probability: float = 0.0,
        cross_link_exchange_probability: float = 0.0,
        pivot_cbmc_probability: float = 0.0,
        swap_multi_trial_cbmc_probability: float = 0.0,
        reinsertion_multi_trial_cbmc_probability: float = 0.0,
        gibbs_swap_multi_trial_cbmc_probability: float = 0.0,
        identity_change_probability: float = 0.0,
        gibbs_volume_change_probability: float = 0.0,
    ) -> None:
        """
        Initialize an :class:`MCMoveProbabilities`.

        The probabilities are given by keyword. ``identity_change_probability``
        and ``gibbs_volume_change_probability`` are the keywords of earlier
        releases for ``identity_change_cbmc_probability`` and
        ``gibbs_volume_probability``.
        """
        ...

    def __getitem__(self, move: Move.Types) -> float:
        """Return the probability of a move."""
        ...

    def __setitem__(self, move: Move.Types, probability: float) -> None:
        """Set the probability of a move."""
        ...

    def get_probability(self, move: Move.Types) -> float:
        """Return the probability of a move."""
        ...

    def set_probability(self, move: Move.Types, probability: float) -> None:
        """Set the probability of a move."""
        ...

    def join(self, other: 'MCMoveProbabilities') -> None:
        """Add the probabilities of another set to this one."""
        ...
