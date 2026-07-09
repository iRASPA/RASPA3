module;

export module cbmc_chain_data;

import std;

import atom;
import molecule;
import double3x3;
import double3;
import randomnumbers;
import running_energy;

export struct ChainGrowData
{
  Molecule molecule;
  std::vector<Atom> atoms;  ///< All atoms of the grown molecule.
  RunningEnergy energies;
  double RosenbluthWeight;
  /// Retained partial Rosenbluth weight for the multiple-first-bead reinsertion scheme (Esselink et
  /// al., 'r' in Eq. 16-18): the Rosenbluth weight minus the Boltzmann factor of the selected trial,
  /// carried from grow to retrace. Zero for the moves that do not use it.
  double storedR;

  ChainGrowData() : molecule(), atoms(), energies(), RosenbluthWeight(), storedR() {}

  ChainGrowData(const Molecule &molecule, std::vector<Atom> atoms, RunningEnergy energies, double RosenbluthWeight,
                double storedR) noexcept
      : molecule(molecule), atoms(atoms), energies(energies), RosenbluthWeight(RosenbluthWeight), storedR(storedR)
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
