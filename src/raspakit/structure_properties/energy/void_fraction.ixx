module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <optional>
#include <vector>
#endif

export module energy_void_fraction;

#ifdef USE_STD_IMPORT
import std;
#endif

import int3;
import double2;
import double3;
import double3x3;
import forcefield;
import framework;

export struct EnergyVoidFraction
{
  EnergyVoidFraction();
  ~EnergyVoidFraction();

  void run(const ForceField &forceField, const Framework &framework);
};
