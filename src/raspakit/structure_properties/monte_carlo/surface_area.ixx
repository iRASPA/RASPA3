module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>
#include <optional>
#include <tuple>
#endif

export module mc_surface_area;

#ifdef USE_STD_IMPORT
import std;
#endif

import framework;
import forcefield;

export struct MC_SurfaceArea
{
  std::vector<double> data;

  MC_SurfaceArea() {};

  void run(const ForceField &forceField, const Framework &framework, double wellDepthFactor,
           std::string probePseudoAtom, std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps) const;
};
