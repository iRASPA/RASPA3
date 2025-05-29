module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <array>
#include <vector>
#include <optional>
#endif

export module energy_surface_area;

import int3;
import double2;
import double3;
import double3x3;
import forcefield;
import framework;

export struct EnergySurfaceArea
{
  EnergySurfaceArea();

  void run(const ForceField &forceField, const Framework &framework);
};
