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

export struct ChainData
{
  Molecule molecule;
  std::vector<Atom> atom;
  std::vector<double3> electricField;
  RunningEnergy energies;
  double RosenbluthWeight;
  double storedR;

  ChainData(const Molecule &molecule, std::vector<Atom> atom, RunningEnergy energies, double RosenbluthWeight,
            double storedR) noexcept
      : molecule(molecule),
        atom(atom),
        electricField(atom.size()),
        energies(energies),
        RosenbluthWeight(RosenbluthWeight),
        storedR(storedR)
  {
  }
};
