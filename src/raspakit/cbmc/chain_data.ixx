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
  /// Natural logarithm of the Rosenbluth weight, accumulated per growth step so it stays exact where
  /// the raw product underflows. The weight of a long chain is a product of hundreds of per-step
  /// factors of order exp(-beta u); for a polymer it drops below the smallest double (~1e-308) and the
  /// raw 'RosenbluthWeight' flushes to zero, turning every acceptance ratio W_new/W_old into 0/0 (NaN,
  /// silent reject) or x/0 (inf, unconditional accept). Ratio-based acceptances (reinsertion, partial
  /// reinsertion, ideal-gas conformation equilibration) must therefore use
  /// exp(logRosenbluthWeight_new - logRosenbluthWeight_old) instead of the raw quotient.
  double logRosenbluthWeight;

  ChainGrowData() : molecule(), atoms(), energies(), RosenbluthWeight(), storedR(), logRosenbluthWeight() {}

  ChainGrowData(const Molecule &molecule, std::vector<Atom> atoms, RunningEnergy energies, double RosenbluthWeight,
                double storedR) noexcept
      : molecule(molecule),
        atoms(atoms),
        energies(energies),
        RosenbluthWeight(RosenbluthWeight),
        storedR(storedR),
        logRosenbluthWeight(std::log(RosenbluthWeight))
  {
  }

  ChainGrowData(const Molecule &molecule, std::vector<Atom> atoms, RunningEnergy energies, double RosenbluthWeight,
                double storedR, double logRosenbluthWeight) noexcept
      : molecule(molecule),
        atoms(atoms),
        energies(energies),
        RosenbluthWeight(RosenbluthWeight),
        storedR(storedR),
        logRosenbluthWeight(logRosenbluthWeight)
  {
  }
};

export struct ChainRetraceData
{
  RunningEnergy energies;
  double RosenbluthWeight;
  double storedR;
  /// Exact logarithm of the Rosenbluth weight; see ChainGrowData::logRosenbluthWeight.
  double logRosenbluthWeight;

  ChainRetraceData() : energies(), RosenbluthWeight(), storedR(), logRosenbluthWeight() {}

  ChainRetraceData(RunningEnergy energies, double RosenbluthWeight, double storedR) noexcept
      : energies(energies),
        RosenbluthWeight(RosenbluthWeight),
        storedR(storedR),
        logRosenbluthWeight(std::log(RosenbluthWeight))
  {
  }

  ChainRetraceData(RunningEnergy energies, double RosenbluthWeight, double storedR,
                   double logRosenbluthWeight) noexcept
      : energies(energies),
        RosenbluthWeight(RosenbluthWeight),
        storedR(storedR),
        logRosenbluthWeight(logRosenbluthWeight)
  {
  }
};
