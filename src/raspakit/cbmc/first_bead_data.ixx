module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <vector>
#endif

export module cbmc_first_bead_data;

#ifdef USE_STD_IMPORT
import std;
#endif

import atom;
import double3x3;
import double3;
import randomnumbers;
import running_energy;

export struct FirstBeadData
{
  Atom atom;
  RunningEnergy energies;
  double RosenbluthWeight;
  double storedR;

  FirstBeadData() noexcept = delete;
  FirstBeadData(Atom atom, RunningEnergy energies, double RosenbluthWeight, double storedR) noexcept
      : atom(atom), energies(energies), RosenbluthWeight(RosenbluthWeight), storedR(storedR)
  {
  }
};
