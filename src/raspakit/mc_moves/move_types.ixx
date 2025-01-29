module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#include <map>
#endif

export module mc_moves_move_types;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <map>;
#endif

export enum class MoveTypes : size_t {
  Translation = 0,
  RandomTranslation = 1,
  Rotation = 2,
  RandomRotation = 3,
  VolumeChange = 4,
  ReinsertionCBMC = 5,
  IdentityChangeCBMC = 6,
  Swap = 7,
  SwapCBMC = 8,
  SwapCFCMC = 9,
  SwapCBCFCMC = 10,
  GibbsVolume = 11,
  GibbsSwapCBMC = 12,
  GibbsSwapCFCMC = 13,
  Widom = 14,
  WidomCFCMC = 15,
  WidomCBCFCMC = 16,
  ParallelTempering = 17,
  HybridMC = 18,
  Count
};

export std::map<MoveTypes, std::string> moveNames = 
{
  {MoveTypes::Translation, "Translation"},
  {MoveTypes::RandomTranslation, "Random translation"},
  {MoveTypes::Rotation, "Rotation"},
  {MoveTypes::RandomRotation, "Random rotation"},
  {MoveTypes::VolumeChange, "Volume change"},
  {MoveTypes::ReinsertionCBMC, "Reinsertion (CBMC)"},
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