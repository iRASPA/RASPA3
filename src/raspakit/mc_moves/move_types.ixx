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
    Count = 34
  };

  static std::array<std::string, std::to_underlying(Move::Types::Count)> moveNames;
};


export inline std::unordered_set<Move::Types> componentMoves = {Move::Types::Translation,        Move::Types::RandomTranslation,
                                                                Move::Types::Rotation,           Move::Types::RandomRotation,
                                                                Move::Types::ReinsertionCBMC,    Move::Types::PartialReinsertionCBMC,
                                                                Move::Types::IdentityChangeCBMC, Move::Types::Swap,
                                                                Move::Types::SwapCBMC,           Move::Types::PairSwapCBMC,
                                                                Move::Types::PairSwap,
                                                                Move::Types::PairSwapCFCMC,      Move::Types::PairSwapCBCFCMC,
                                                                Move::Types::SwapCFCMC,          Move::Types::SwapCBCFCMC,
                                                                Move::Types::Widom,
                                                                Move::Types::WidomCFCMC,         Move::Types::WidomCBCFCMC};

export inline std::unordered_set<Move::Types> systemMoves = {Move::Types::VolumeChange, Move::Types::AnisotropicVolumeChange,
                                                           Move::Types::HybridMC};

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

