module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#endif

export module mc_moves_move_types;

#ifndef USE_LEGACY_HEADERS
import <string>;
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
