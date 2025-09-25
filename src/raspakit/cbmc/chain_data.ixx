module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <vector>
#endif

export module cbmc_chain_data;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import atom;
import molecule;
import double3x3;
import double3;
import randomnumbers;
import running_energy;

export struct ChainGrowData
{
  Molecule molecule;
  std::vector<Atom> atom;
  RunningEnergy energies;
  double RosenbluthWeight;
  double storedR;

  ChainGrowData() : molecule(), atom(), energies(), RosenbluthWeight(), storedR() {}

  ChainGrowData(const Molecule &molecule, std::vector<Atom> atom, RunningEnergy energies, double RosenbluthWeight,
                double storedR) noexcept
      : molecule(molecule), atom(atom), energies(energies), RosenbluthWeight(RosenbluthWeight), storedR(storedR)
  {
  }
};

export struct ChainRetraceData
{
  RunningEnergy energies;
  double RosenbluthWeight;
  double storedR;

  ChainRetraceData() : energies(), RosenbluthWeight(), storedR() {}

  ChainRetraceData(RunningEnergy energies, double RosenbluthWeight, double storedR) noexcept
      : energies(energies), RosenbluthWeight(RosenbluthWeight), storedR(storedR)
  {
  }
};
