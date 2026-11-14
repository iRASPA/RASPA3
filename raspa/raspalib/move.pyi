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
            TRANSLATION_SMART_MC: Force-biased (smart MC) translation of a single molecule.
            TRANSLATION_SMART_MC_ALL: Force-biased (smart MC) translation of all molecules simultaneously.
            ROTATION_SMART_MC: Torque-biased (smart MC) rotation of a single molecule.
            ROTATION_SMART_MC_ALL: Torque-biased (smart MC) rotation of all molecules simultaneously.
            TRANSLATION_ROTATION_SMART_MC: Combined force- and torque-biased (smart MC) move of a single molecule.
            VOLUME_CHANGE: Simulation-box volume-change move.
            ANISOTROPIC_VOLUME_CHANGE: Anisotropic simulation-box volume-change move.
            REINSERTION_CBMC: Configurational-bias reinsertion move.
            PARTIAL_REINSERTION_CBMC: Partial configurational-bias reinsertion.
            IDENTITY_CHANGE_CBMC: Configurational-bias identity-change move.
            IDENTITY_SWITCH_CBMC: Canonical pair identity exchange (MCCCS-MN 'swatch').
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
            GROUP_SWAP: Insertion/deletion of a neutral group of molecules (conventional).
            GROUP_SWAP_CBMC: Insertion/deletion of a neutral group of molecules (CBMC).
            GROUP_SWAP_CFCMC: CFCMC insertion/deletion of a neutral group of molecules.
            GROUP_SWAP_CBCFCMC: CB/CFCMC insertion/deletion of a neutral group of molecules.
            TETHERED_PROTON_HOP: Relocation of a tethered proton among its discrete candidate sites.
            PIVOT: Rotation of a chain part about a randomly selected non-ring bond.
            CRANKSHAFT: Rotation of a small chain segment about the axis through its two anchors.
            REPTATION: Slithering-snake move: one repeat unit removed at one chain end, regrown at the other.
            CONCERTED_ROTATION: ConRot: driver torsion change with re-bridging of the following trimer (Dodd, Boone, Theodorou).
            BEAD_DISPLACEMENT: Random displacement of a single flexible bead.
            BEAD_FLIP: Bond-length-preserving rotation of a single bead (kink jump / end rotation).
            DOUBLE_BRIDGING: Connectivity-altering exchange of chain tails between two chains (Karayiannis, Mavrantzas, Theodorou).
            INTRAMOLECULAR_DOUBLE_REBRIDGING: Connectivity-altering reversal of a chain segment via two closures (IDR).
            CROSS_LINK_SWAP: Bond-count-conserving relocation of one end of a cross-link (system move).
            CROSS_LINK_FORMATION_SCISSION: Formation or scission of a cross-link between two reactive sites (system move).
            CROSS_LINK_EXCHANGE: Two cross-links trade partners, (a,b)+(c,d) -> (a,c)+(b,d) (system move).
            PIVOT_CBMC: Pivot with Rosenbluth selection of the rotation angle among k trial angles.
            SWAP_MULTI_TRIAL_CBMC: Multiple-try CBMC insertion/deletion: M trial molecules, selection by weight.
            REINSERTION_MULTI_TRIAL_CBMC: Multiple-try CBMC reinsertion.
            GIBBS_SWAP_MULTI_TRIAL_CBMC: Multiple-try CBMC Gibbs swap (cross-system move).
        """

        TRANSLATION = 0
        RANDOM_TRANSLATION = 1
        ROTATION = 2
        RANDOM_ROTATION = 3
        TRANSLATION_SMART_MC = 4
        TRANSLATION_SMART_MC_ALL = 5
        ROTATION_SMART_MC = 6
        ROTATION_SMART_MC_ALL = 7
        TRANSLATION_ROTATION_SMART_MC = 8
        VOLUME_CHANGE = 9
        ANISOTROPIC_VOLUME_CHANGE = 10
        REINSERTION_CBMC = 11
        PARTIAL_REINSERTION_CBMC = 12
        IDENTITY_CHANGE_CBMC = 13
        IDENTITY_SWITCH_CBMC = 14
        SWAP = 15
        SWAP_CBMC = 16
        SWAP_CFCMC = 17
        SWAP_CBCFCMC = 18
        GIBBS_VOLUME = 19
        GIBBS_SWAP_CBMC = 20
        GIBBS_CONVENTIONAL_CFCMC = 21
        GIBBS_CONVENTIONAL_CBCFCMC = 22
        GIBBS_SWAP_CFCMC = 23
        GIBBS_SWAP_CBCFCMC = 24
        GIBBS_IDENTITY_CHANGE_CBMC = 25
        WIDOM = 26
        WIDOM_CFCMC = 27
        WIDOM_CBCFCMC = 28
        PARALLEL_TEMPERING = 29
        HYBRID_MC = 30
        REACTION_CBMC = 31
        REACTION_CONVENTIONAL_CFCMC = 32
        REACTION_CONVENTIONAL_CBCFCMC = 33
        REACTION_CFCMC = 34
        REACTION_CBCFCMC = 35
        PAIR_SWAP = 36
        PAIR_SWAP_CBMC = 37
        PAIR_SWAP_CFCMC = 38
        PAIR_SWAP_CBCFCMC = 39
        GROUP_SWAP = 40
        GROUP_SWAP_CBMC = 41
        GROUP_SWAP_CFCMC = 42
        GROUP_SWAP_CBCFCMC = 43
        TETHERED_PROTON_HOP = 44
        PIVOT = 45
        CRANKSHAFT = 46
        REPTATION = 47
        CONCERTED_ROTATION = 48
        BEAD_DISPLACEMENT = 49
        BEAD_FLIP = 50
        DOUBLE_BRIDGING = 51
        INTRAMOLECULAR_DOUBLE_REBRIDGING = 52
        CROSS_LINK_SWAP = 53
        CROSS_LINK_FORMATION_SCISSION = 54
        CROSS_LINK_EXCHANGE = 55
        PIVOT_CBMC = 56
        SWAP_MULTI_TRIAL_CBMC = 57
        REINSERTION_MULTI_TRIAL_CBMC = 58
        GIBBS_SWAP_MULTI_TRIAL_CBMC = 59
