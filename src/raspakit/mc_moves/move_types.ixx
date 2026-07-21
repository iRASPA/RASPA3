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
    TranslationSmartMC = 4,           // force-biased (smart MC) translation of a single molecule
    TranslationSmartMCAll = 5,        // force-biased (smart MC) translation of all molecules simultaneously
    RotationSmartMC = 6,              // torque-biased (smart MC) rotation of a single molecule
    RotationSmartMCAll = 7,           // torque-biased (smart MC) rotation of all molecules simultaneously
    TranslationRotationSmartMC = 8,   // combined force- and torque-biased (smart MC) move of a single molecule
    VolumeChange = 9,
    AnisotropicVolumeChange = 10,
    ReinsertionCBMC = 11,
    PartialReinsertionCBMC = 12,
    IdentityChangeCBMC = 13,
    Swap = 14,
    SwapCBMC = 15,
    SwapCFCMC = 16,
    SwapCBCFCMC = 17,
    GibbsVolume = 18,
    GibbsSwapCBMC = 19,
    GibbsConventionalCFCMC = 20,      // parallel CFCMC Gibbs
    GibbsConventionalCBCFCMC = 21,    // parallel CFCMC/CBMC Gibbs
    GibbsSwapCFCMC = 22,              // serial CFCMC Gibbs
    GibbsSwapCBCFCMC = 23,            // serial CFCMC/CBMC Gibbs
    GibbsIdentityChangeCBMC = 24,
    Widom = 25,
    WidomCFCMC = 26,
    WidomCBCFCMC = 27,
    ParallelTempering = 28,
    HybridMC = 29,
    ReactionCBMC = 30,
    ReactionConventionalCFCMC = 31,    // parallel CFCMC reaction MC
    ReactionConventionalCBCFCMC = 32,  // parallel CFCMC/CBMC reaction MC
    ReactionCFCMC = 33,                // serial CFCMC reaction MC
    ReactionCBCFCMC = 34,              // serial CFCMC/CBMC reaction MC
    PairSwap = 35,
    PairSwapCBMC = 36,
    PairSwapCFCMC = 37,
    PairSwapCBCFCMC = 38,
    GroupSwap = 39,                   // insertion/deletion of a neutral group of molecules (conventional)
    GroupSwapCBMC = 40,               // insertion/deletion of a neutral group of molecules (CBMC)
    GroupSwapCFCMC = 41,              // CFCMC insertion/deletion of a neutral group of molecules
    GroupSwapCBCFCMC = 42,            // CB/CFCMC insertion/deletion of a neutral group of molecules
    Count = 43
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
                                                                Move::Types::TranslationRotationSmartMC,
                                                                Move::Types::ReinsertionCBMC,    Move::Types::PartialReinsertionCBMC,
                                                                Move::Types::IdentityChangeCBMC, Move::Types::Swap,
                                                                Move::Types::SwapCBMC,           Move::Types::PairSwapCBMC,
                                                                Move::Types::PairSwap,
                                                                Move::Types::GroupSwap,          Move::Types::GroupSwapCBMC,
                                                                Move::Types::GroupSwapCFCMC,     Move::Types::GroupSwapCBCFCMC,
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

