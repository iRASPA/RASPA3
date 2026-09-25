module;

export module cbmc_results;

import std;

import atom;
import molecule;
import running_energy;

// The result types of the CBMC growth stages: the first bead (FirstBeadData) and the remaining chain
// (ChainGrowData / ChainRetraceData). The entry points in the 'cbmc' module return the combination.

/// Result of placing (or retracing) the first bead of a molecule with one of the first-bead schemes.
/// The first-bead weight is a single-bead average of Boltzmann factors (order one), so it is stored
/// as the plain weight; the chain weight below is stored as a logarithm.
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
