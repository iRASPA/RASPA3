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
            REINSERTION_CBMC: Configurational-bias reinsertion move.
            PARTIAL_REINSERTION_CBMC: Partial configurational-bias reinsertion.
            IDENTITY_CHANGE_CBMC: Configurational-bias identity-change move.
            SWAP: Particle swap move.
            SWAP_CBMC: Configurational-bias particle swap move.
            SWAP_CFCMC: Continuous fractional-component swap move.
            SWAP_CBCFCMC: Combined CBMC/CFCMC swap move.
            GIBBS_VOLUME: Gibbs-ensemble volume-change move.
            GIBBS_SWAP_CBMC: Gibbs-ensemble CBMC swap move.
            GIBBS_SWAP_CFCMC: Gibbs-ensemble CFCMC swap move.
            WIDOM: Widom insertion move.
            WIDOM_CFCMC: Widom insertion using CFCMC.
            WIDOM_CBCFCMC: Widom insertion using combined CBMC/CFCMC.
            PARALLEL_TEMPERING: Parallel-tempering exchange move.
            HYBRID_MC: Hybrid Monte Carlo move.
        """

        TRANSLATION = 0
        RANDOM_TRANSLATION = 1
        ROTATION = 2
        RANDOM_ROTATION = 3
        VOLUME_CHANGE = 4
        REINSERTION_CBMC = 5
        PARTIAL_REINSERTION_CBMC = 6
        IDENTITY_CHANGE_CBMC = 7
        SWAP = 8
        SWAP_CBMC = 9
        SWAP_CFCMC = 10
        SWAP_CBCFCMC = 11
        GIBBS_VOLUME = 12
        GIBBS_SWAP_CBMC = 13
        GIBBS_SWAP_CFCMC = 14
        WIDOM = 15
        WIDOM_CFCMC = 16
        WIDOM_CBCFCMC = 17
        PARALLEL_TEMPERING = 18
        HYBRID_MC = 19
