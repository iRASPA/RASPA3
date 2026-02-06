module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <optional>
#include <vector>
#include <string>
#endif

export module energy_surface_area;

#ifdef USE_STD_IMPORT
import std;
#endif

import int3;
import double2;
import double3;
import double3x3;
import forcefield;
import framework;

export struct EnergySurfaceArea
{
  EnergySurfaceArea();

  void run(const ForceField &forceField, const Framework &framework, double isoValue,
                         std::string probePseudoAtom);
};
