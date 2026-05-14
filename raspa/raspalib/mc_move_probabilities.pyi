class MCMoveProbabilities():
    """
    Store Monte Carlo move-selection probabilities for a component.

    This container holds the relative probabilities for the supported MC moves.
    Values are typically normalized internally by the backend after
    initialization.
    """

    def __init__(
        self,
        translation_probability: float = 0.0,
        random_translation_probability: float = 0.0,
        rotation_probability: float = 0.0,
        random_rotation_probability: float = 0.0,
        volume_change_probability: float = 0.0,
        reinsertion_cbmc_probability: float = 0.0,
        partial_reinsertion_cbmc_probability: float = 0.0,
        identity_change_probability: float = 0.0,
        swap_probability: float = 0.0,
        swap_cbmc_probability: float = 0.0,
        swap_cfcmc_probability: float = 0.0,
        swap_cbcfcmc_probability: float = 0.0,
        gibbs_volume_change_probability: float = 0.0,
        gibbs_swap_cbmc_probability: float = 0.0,
        gibbs_swap_cfcmc_probability: float = 0.0,
        widom_probability: float = 0.0,
        widom_cfcmc_probability: float = 0.0,
        widom_cbcfcmc_probability: float = 0.0,
        parallel_tempering_probability: float = 0.0,
        hybrid_mc_probability: float = 0.0
    ) -> None:
        ...
        """
        Initialize an :class:`MCMoveProbabilities`.

        Args:
            translation_probability: Probability weight for translation moves.
            random_translation_probability: Probability weight for random
                translation moves.
            rotation_probability: Probability weight for rotation moves.
            random_rotation_probability: Probability weight for random rotation
                moves.
            volume_change_probability: Probability weight for volume-change
                moves.
            reinsertion_cbmc_probability: Probability weight for CBMC
                reinsertion moves.
            partial_reinsertion_cbmc_probability: Probability weight for
                partial CBMC reinsertion moves.
            identity_change_probability: Probability weight for identity-change
                moves.
            swap_probability: Probability weight for particle swap moves.
            swap_cbmc_probability: Probability weight for CBMC swap moves.
            swap_cfcmc_probability: Probability weight for CFCMC swap moves.
            swap_cbcfcmc_probability: Probability weight for combined CBMC-CFCMC
                swap moves.
            gibbs_volume_change_probability: Probability weight for Gibbs
                volume-change moves.
            gibbs_swap_cbmc_probability: Probability weight for Gibbs CBMC swap
                moves.
            gibbs_swap_cfcmc_probability: Probability weight for Gibbs CFCMC
                swap moves.
            widom_probability: Probability weight for Widom insertion moves.
            widom_cfcmc_probability: Probability weight for Widom-CFCMC moves.
            widom_cbcfcmc_probability: Probability weight for combined Widom
                CBMC-CFCMC moves.
            parallel_tempering_probability: Probability weight for parallel
                tempering exchanges.
            hybrid_mc_probability: Probability weight for hybrid MC moves.
        """


