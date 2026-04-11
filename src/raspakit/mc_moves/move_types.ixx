module;

export module mc_moves_move_types;

import std;

export enum class MoveTypes : std::size_t {
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
  Widom = 15,
  WidomCFCMC = 16,
  WidomCBCFCMC = 17,
  ParallelTempering = 18,
  HybridMC = 19,
  Count = 20
};

export inline std::array<std::string, std::to_underlying(MoveTypes::Count)> moveNames = {
    "Translation",
    "Random translation",
    "Rotation",
    "Random rotation",
    "Volume change",
    "Reinsertion (CBMC)",
    "Partial reinsertion (CBMC)",
    "Identity change (CBMC)",
    "Swap",
    "Swap (CBMC)",
    "Swap (CFCMC)",
    "Swap (CB/CFCMC)",
    "Gibbs volume",
    "Gibbs swap (CBMC)",
    "Gibbs swap (CFCMC)",
    "Widom",
    "Widom (CFCMC)",
    "Widom (CB/CFCMC)",
    "Parallel tempering",
    "Hybrid MC"
};

export inline std::unordered_set<MoveTypes> componentMoves = {MoveTypes::Translation,        MoveTypes::RandomTranslation,
                                                              MoveTypes::Rotation,           MoveTypes::RandomRotation,
                                                              MoveTypes::ReinsertionCBMC,    MoveTypes::PartialReinsertionCBMC,
                                                              MoveTypes::IdentityChangeCBMC, MoveTypes::Swap,
                                                              MoveTypes::SwapCBMC,           MoveTypes::SwapCFCMC,
                                                              MoveTypes::SwapCBCFCMC,        MoveTypes::Widom,
                                                              MoveTypes::WidomCFCMC,         MoveTypes::WidomCBCFCMC};

export inline std::unordered_set<MoveTypes> systemMoves = {MoveTypes::VolumeChange, MoveTypes::HybridMC};

export inline std::unordered_set<MoveTypes> crossSystemMoves = {MoveTypes::GibbsVolume, MoveTypes::GibbsSwapCBMC,
                                                                MoveTypes::GibbsSwapCFCMC, MoveTypes::ParallelTempering};

export inline std::unordered_set<MoveTypes> groupMoves = {};

