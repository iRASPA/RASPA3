import enum

class Move:
    """
    Represent a Monte Carlo move descriptor used by RASPA.

    The nested :class:`Types` enumeration lists all supported move families
    used in move-probability and move-statistics APIs.
    """
    def __init__(self) -> None:
        """Initialize a :class:`Move` helper object."""
        ...

    class Types(enum.IntEnum):
        """
        Enumerate supported Monte Carlo move types.

        Members:
            TRANSLATION: Local particle translation move.
            RANDOM_TRANSLATION: Full-range random particle translation move.
            ROTATION: Local particle rotation move.
            RANDOM_ROTATION: Full-range random particle rotation move.
            VOLUME_CHANGE: Simulation-box volume-change move.
            ANISOTROPIC_VOLUME_CHANGE: Anisotropic simulation-box volume-change move.
            REINSERTION_CBMC: Configurational-bias reinsertion move.
            PARTIAL_REINSERTION_CBMC: Partial configurational-bias reinsertion.
            IDENTITY_CHANGE_CBMC: Configurational-bias identity-change move.
            SWAP: Particle swap move.
            SWAP_CBMC: Configurational-bias particle swap move.
            SWAP_CFCMC: Continuous fractional-component swap move.
            SWAP_CBCFCMC: Combined CBMC/CFCMC swap move.
            GIBBS_VOLUME: Gibbs-ensemble volume-change move.
            GIBBS_SWAP_CBMC: Gibbs-ensemble CBMC swap move.
            GIBBS_CONVENTIONAL_CFCMC: Parallel CFCMC Gibbs swap move.
            GIBBS_CONVENTIONAL_CBCFCMC: Parallel CB/CFCMC Gibbs swap move.
            GIBBS_SWAP_CFCMC: Serial CFCMC Gibbs swap move.
            GIBBS_SWAP_CBCFCMC: Serial CFCMC/CBMC Gibbs swap move.
            GIBBS_IDENTITY_CHANGE_CBMC: Gibbs-ensemble CBMC identity-change move.
            WIDOM: Widom insertion move.
            WIDOM_CFCMC: Widom insertion using CFCMC.
            WIDOM_CBCFCMC: Widom insertion using combined CBMC/CFCMC.
            PARALLEL_TEMPERING: Parallel-tempering exchange move.
            HYBRID_MC: Hybrid Monte Carlo move.
            REACTION_CBMC: Reaction move using CBMC.
            REACTION_CONVENTIONAL_CFCMC: Parallel CFCMC reaction move.
            REACTION_CONVENTIONAL_CBCFCMC: Parallel CB/CFCMC reaction move.
            REACTION_CFCMC: Serial CFCMC reaction move.
            REACTION_CBCFCMC: Serial CB/CFCMC reaction move.
            PAIR_SWAP: Pair swap move.
            PAIR_SWAP_CBMC: Configurational-bias pair swap move.
            PAIR_SWAP_CFCMC: CFCMC pair swap move.
            PAIR_SWAP_CBCFCMC: Combined CBMC/CFCMC pair swap move.
        """

        TRANSLATION = 0
        RANDOM_TRANSLATION = 1
        ROTATION = 2
        RANDOM_ROTATION = 3
        VOLUME_CHANGE = 4
        ANISOTROPIC_VOLUME_CHANGE = 5
        REINSERTION_CBMC = 6
        PARTIAL_REINSERTION_CBMC = 7
        IDENTITY_CHANGE_CBMC = 8
        SWAP = 9
        SWAP_CBMC = 10
        SWAP_CFCMC = 11
        SWAP_CBCFCMC = 12
        GIBBS_VOLUME = 13
        GIBBS_SWAP_CBMC = 14
        GIBBS_CONVENTIONAL_CFCMC = 15
        GIBBS_CONVENTIONAL_CBCFCMC = 16
        GIBBS_SWAP_CFCMC = 17
        GIBBS_SWAP_CBCFCMC = 18
        GIBBS_IDENTITY_CHANGE_CBMC = 19
        WIDOM = 20
        WIDOM_CFCMC = 21
        WIDOM_CBCFCMC = 22
        PARALLEL_TEMPERING = 23
        HYBRID_MC = 24
        REACTION_CBMC = 25
        REACTION_CONVENTIONAL_CFCMC = 26
        REACTION_CONVENTIONAL_CBCFCMC = 27
        REACTION_CFCMC = 28
        REACTION_CBCFCMC = 29
        PAIR_SWAP = 30
        PAIR_SWAP_CBMC = 31
        PAIR_SWAP_CFCMC = 32
        PAIR_SWAP_CBCFCMC = 33
