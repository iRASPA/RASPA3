module;

export module cbmc_chain_data;

import std;

import atom;
import molecule;
import double3x3;
import double3;
import randomnumbers;
import running_energy;

/// Result of growing a molecule (or the remainder of a molecule) with CBMC or recoil growth.
///
/// The Rosenbluth weight is stored ONLY as its natural logarithm, 'logRosenbluthWeight', accumulated
/// per growth step so it stays exact where the raw product underflows. The weight of a long chain is
/// a product of hundreds of per-step factors of order exp(-beta u); for a polymer it drops below the
/// smallest double (~1e-308) and a raw weight flushes to zero, turning every acceptance ratio
/// W_new/W_old into 0/0 (NaN, silent reject) or x/0 (inf, unconditional accept). Acceptance rules
/// must therefore be written in log space: exp(logRosenbluthWeight_new - logRosenbluthWeight_old),
/// and corrections (dual cut-off, Ewald) are applied with 'multiplyRosenbluthWeight'. The raw weight
/// is available through 'rosenbluthWeight()' for reporting and averaging (Widom) only.
export struct ChainGrowData
{
  Molecule molecule;
  std::vector<Atom> atoms;  ///< All atoms of the grown molecule.
  RunningEnergy energies;
  /// Retained partial Rosenbluth weight for the multiple-first-bead reinsertion scheme (Esselink et
  /// al., 'r' in Eq. 16-18): the Rosenbluth weight minus the Boltzmann factor of the selected trial,
  /// carried from grow to retrace. Zero for the moves that do not use it.
  double storedR;
  /// Natural logarithm of the Rosenbluth weight; the single stored representation (see above).
  double logRosenbluthWeight;

  ChainGrowData() : molecule(), atoms(), energies(), storedR(), logRosenbluthWeight() {}

  ChainGrowData(const Molecule &molecule, std::vector<Atom> atoms, RunningEnergy energies, double logRosenbluthWeight,
                double storedR) noexcept
      : molecule(molecule),
        atoms(std::move(atoms)),
        energies(energies),
        storedR(storedR),
        logRosenbluthWeight(logRosenbluthWeight)
  {
  }

  /// The Rosenbluth weight itself, exp(logRosenbluthWeight). Underflows to zero (or overflows) for
  /// long chains: use for reporting and averaging only, never inside an acceptance ratio.
  [[nodiscard]] double rosenbluthWeight() const noexcept { return std::exp(logRosenbluthWeight); }

  /// Multiplies the Rosenbluth weight by exp(logFactor), e.g. a dual cut-off or Ewald correction
  /// exp(-beta dU) passed as -beta * dU.
  void multiplyRosenbluthWeight(double logFactor) noexcept { logRosenbluthWeight += logFactor; }
};

/// Result of retracing an existing molecule; the weight is stored as its logarithm exactly as in
/// 'ChainGrowData'.
export struct ChainRetraceData
{
  RunningEnergy energies;
  double storedR;
  /// Natural logarithm of the Rosenbluth weight; the single stored representation.
  double logRosenbluthWeight;

  ChainRetraceData() : energies(), storedR(), logRosenbluthWeight() {}

  ChainRetraceData(RunningEnergy energies, double logRosenbluthWeight, double storedR) noexcept
      : energies(energies), storedR(storedR), logRosenbluthWeight(logRosenbluthWeight)
  {
  }

  /// The Rosenbluth weight itself, exp(logRosenbluthWeight); for reporting and averaging only.
  [[nodiscard]] double rosenbluthWeight() const noexcept { return std::exp(logRosenbluthWeight); }

  /// Multiplies the Rosenbluth weight by exp(logFactor).
  void multiplyRosenbluthWeight(double logFactor) noexcept { logRosenbluthWeight += logFactor; }
};
