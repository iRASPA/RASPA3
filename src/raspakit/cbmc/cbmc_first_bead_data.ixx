module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#endif

export module cbmc_first_bead_data;

#ifndef USE_LEGACY_HEADERS
import <vector>;
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
  FirstBeadData(Atom atom, RunningEnergy energies, double RosenbluthWeight, double storedR) noexcept :
      atom(atom), energies(energies), RosenbluthWeight(RosenbluthWeight), storedR(storedR)
  {
  }
};

