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
    ReinsertionCBMC = 5,
    PartialReinsertionCBMC = 6,
    IdentityChangeCBMC = 7,
    Swap = 8,
    SwapCBMC = 9,
    SwapCFCMC = 10,
    SwapCBCFCMC = 11,
    GibbsVolume = 12,
    GibbsSwapCBMC = 13,
    GibbsSwapCFCMC = 14,
    GibbsIdentityChangeCBMC = 15,
    Widom = 16,
    WidomCFCMC = 17,
    WidomCBCFCMC = 18,
    ParallelTempering = 19,
    HybridMC = 20,
    ReactionCBMC = 21,
    ReactionConventionalCFCMC = 22,
    ReactionConventionalCFCMCCBMC = 23,
    ReactionCFCMC = 24,
    ReactionCFCMCCBMC = 25,
    GibbsSwapCBCFCMC = 26,
    GibbsConventionalCFCMC = 27,
    GibbsConventionalCFCMCCBMC = 28,
    PairSwap = 29,
    PairSwapCBMC = 30,
    AnisotropicVolumeChange = 31,
    Count = 32
  };

  static std::array<std::string, std::to_underlying(Move::Types::Count)> moveNames;
};


export inline std::unordered_set<Move::Types> componentMoves = {Move::Types::Translation,        Move::Types::RandomTranslation,
                                                                Move::Types::Rotation,           Move::Types::RandomRotation,
                                                                Move::Types::ReinsertionCBMC,    Move::Types::PartialReinsertionCBMC,
                                                                Move::Types::IdentityChangeCBMC, Move::Types::Swap,
                                                                Move::Types::SwapCBMC,           Move::Types::PairSwapCBMC,
                                                                Move::Types::PairSwap,
                                                                Move::Types::SwapCFCMC,          Move::Types::SwapCBCFCMC,
                                                                Move::Types::Widom,
                                                                Move::Types::WidomCFCMC,         Move::Types::WidomCBCFCMC};

export inline std::unordered_set<Move::Types> systemMoves = {Move::Types::VolumeChange, Move::Types::AnisotropicVolumeChange,
                                                           Move::Types::HybridMC};

export inline std::unordered_set<Move::Types> crossSystemMoves = {Move::Types::GibbsVolume, Move::Types::GibbsSwapCBMC,
                                                                  Move::Types::GibbsSwapCFCMC,
                                                                  Move::Types::GibbsSwapCBCFCMC,
                                                                  Move::Types::GibbsConventionalCFCMC,
                                                                  Move::Types::GibbsConventionalCFCMCCBMC,
                                                                  Move::Types::GibbsIdentityChangeCBMC,
                                                                  Move::Types::ParallelTempering};

export inline std::unordered_set<Move::Types> groupMoves = {Move::Types::ReactionCBMC, Move::Types::ReactionConventionalCFCMC,
                                                            Move::Types::ReactionConventionalCFCMCCBMC,
                                                            Move::Types::ReactionCFCMC,
                                                            Move::Types::ReactionCFCMCCBMC};

