module;

export module mc_moves_move_types;

import std;

export struct Move
{
  enum class Types : std::size_t
  {
    Translation = 0,
    RandomTranslation = 1,
    Rotation = 2,
    RandomRotation = 3,
    VolumeChange = 4,
    AnisotropicVolumeChange = 5,
    ReinsertionCBMC = 6,
    PartialReinsertionCBMC = 7,
    IdentityChangeCBMC = 8,
    Swap = 9,
    SwapCBMC = 10,
    SwapCFCMC = 11,
    SwapCBCFCMC = 12,
    GibbsVolume = 13,
    GibbsSwapCBMC = 14,
    GibbsConventionalCFCMC = 15,      // parallel CFCMC Gibbs
    GibbsConventionalCBCFCMC = 16,    // parallel CFCMC/CBMC Gibbs
    GibbsSwapCFCMC = 17,              // serial CFCMC Gibbs
    GibbsSwapCBCFCMC = 18,            // serial CFCMC/CBMC Gibbs
    GibbsIdentityChangeCBMC = 19,
    Widom = 20,
    WidomCFCMC = 21,
    WidomCBCFCMC = 22,
    ParallelTempering = 23,
    HybridMC = 24,
    ReactionCBMC = 25,
    ReactionConventionalCFCMC = 26,    // parallel CFCMC reaction MC
    ReactionConventionalCBCFCMC = 27,  // parallel CFCMC/CBMC reaction MC
    ReactionCFCMC = 28,                // serial CFCMC reaction MC
    ReactionCBCFCMC = 29,              // serial CFCMC/CBMC reaction MC
    PairSwap = 30,
    PairSwapCBMC = 31,
    PairSwapCFCMC = 32,
    PairSwapCBCFCMC = 33,
    TranslationSmartMC = 34,     // force-biased (smart MC) translation of a single molecule
    TranslationSmartMCAll = 35,  // force-biased (smart MC) translation of all molecules simultaneously
    RotationSmartMC = 36,        // torque-biased (smart MC) rotation of a single molecule
    RotationSmartMCAll = 37,     // torque-biased (smart MC) rotation of all molecules simultaneously
    Count = 38
  };

  /**
   * \brief Sub-timing labels recorded per Monte Carlo move.
   *
   * This is the union of all timing components any move can report. Each move stores a
   * fixed-size array indexed by these labels; the subset that a given move actually uses
   * is described by Move::timingKeys.
   */
  enum class Timing : std::size_t
  {
    Total = 0,
    ExternalField,
    Framework,
    Molecule,
    NonEwald,
    Ewald,
    Tail,
    ExternalFieldMolecule,
    FrameworkMolecule,
    MoleculeMolecule,
    InsertionTotal,
    InsertionExternalField,
    InsertionFramework,
    InsertionMolecule,
    InsertionNonEwald,
    InsertionEwald,
    InsertionTail,
    DeletionTotal,
    DeletionExternalField,
    DeletionFramework,
    DeletionMolecule,
    DeletionNonEwald,
    DeletionEwald,
    DeletionTail,
    LambdaExternalField,
    LambdaFramework,
    LambdaMolecule,
    LambdaNonEwald,
    LambdaEwald,
    LambdaTail,
    LambdaInterchangeNonEwald,
    LambdaInterchangeEwald,
    LambdaInterchangeTail,
    LambdaChangeNonEwald,
    LambdaChangeEwald,
    LambdaChangeTail,
    LambdaShuffleNonEwald,
    LambdaShuffleEwald,
    LambdaShuffleTail,
    Energy,
    Fugacity,
    Integration,
    Count
  };

  static std::array<std::string, std::to_underlying(Move::Types::Count)> moveNames;

  ///< Display names for each Timing label (index by std::to_underlying(Timing)).
  static std::array<std::string, std::to_underlying(Move::Timing::Count)> timingNames;

  ///< Ordered list of sub-timings (excluding Total) that each move reports; drives output.
  static std::array<std::vector<Move::Timing>, std::to_underlying(Move::Types::Count)> timingKeys;
};


export inline std::unordered_set<Move::Types> componentMoves = {Move::Types::Translation,        Move::Types::RandomTranslation,
                                                                Move::Types::Rotation,           Move::Types::RandomRotation,
                                                                Move::Types::TranslationSmartMC,
                                                                Move::Types::RotationSmartMC,
                                                                Move::Types::ReinsertionCBMC,    Move::Types::PartialReinsertionCBMC,
                                                                Move::Types::IdentityChangeCBMC, Move::Types::Swap,
                                                                Move::Types::SwapCBMC,           Move::Types::PairSwapCBMC,
                                                                Move::Types::PairSwap,
                                                                Move::Types::PairSwapCFCMC,      Move::Types::PairSwapCBCFCMC,
                                                                Move::Types::SwapCFCMC,          Move::Types::SwapCBCFCMC,
                                                                Move::Types::Widom,
                                                                Move::Types::WidomCFCMC,         Move::Types::WidomCBCFCMC};

export inline std::unordered_set<Move::Types> systemMoves = {Move::Types::VolumeChange, Move::Types::AnisotropicVolumeChange,
                                                           Move::Types::HybridMC, Move::Types::TranslationSmartMCAll,
                                                           Move::Types::RotationSmartMCAll};

export inline std::unordered_set<Move::Types> crossSystemMoves = {Move::Types::GibbsVolume, Move::Types::GibbsSwapCBMC,
                                                                  Move::Types::GibbsSwapCFCMC,
                                                                  Move::Types::GibbsSwapCBCFCMC,
                                                                  Move::Types::GibbsConventionalCFCMC,
                                                                  Move::Types::GibbsConventionalCBCFCMC,
                                                                  Move::Types::GibbsIdentityChangeCBMC,
                                                                  Move::Types::ParallelTempering};

export inline std::unordered_set<Move::Types> groupMoves = {Move::Types::ReactionCBMC, Move::Types::ReactionConventionalCFCMC,
                                                            Move::Types::ReactionConventionalCBCFCMC,
                                                            Move::Types::ReactionCFCMC,
                                                            Move::Types::ReactionCBCFCMC};

