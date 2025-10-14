module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <map>
#include <string>
#include <unordered_set>
#endif

export module mc_moves_move_types;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

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
  Count
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

export inline std::map<MoveTypes, std::string> moveNames = {
    {MoveTypes::Translation, "Translation"},
    {MoveTypes::RandomTranslation, "Random translation"},
    {MoveTypes::Rotation, "Rotation"},
    {MoveTypes::RandomRotation, "Random rotation"},
    {MoveTypes::VolumeChange, "Volume change"},
    {MoveTypes::ReinsertionCBMC, "Reinsertion (CBMC)"},
    {MoveTypes::PartialReinsertionCBMC, "Partial reinsertion (CBMC)"},
    {MoveTypes::IdentityChangeCBMC, "Identity change (CBMC)"},
    {MoveTypes::Swap, "Swap"},
    {MoveTypes::SwapCBMC, "Swap (CBMC)"},
    {MoveTypes::SwapCFCMC, "Swap (CFCMC)"},
    {MoveTypes::SwapCBCFCMC, "Swap (CB/CFCMC)"},
    {MoveTypes::GibbsVolume, "Gibbs volume"},
    {MoveTypes::GibbsSwapCBMC, "Gibbs swap (CBMC)"},
    {MoveTypes::GibbsSwapCFCMC, "Gibbs swap (CFCMC)"},
    {MoveTypes::Widom, "Widom"},
    {MoveTypes::WidomCFCMC, "Widom (CFCMC)"},
    {MoveTypes::WidomCBCFCMC, "Widom (CB/CFCMC)"},
    {MoveTypes::ParallelTempering, "Parallel tempering"},
    {MoveTypes::HybridMC, "Hybrid MC"},
};
