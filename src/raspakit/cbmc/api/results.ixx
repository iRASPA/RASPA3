module;

export module cbmc_results;

import std;

import atom;
import molecule;
import running_energy;

// The vocabulary shared by the entry points of the 'cbmc' module and the growth stages behind them:
// the first-bead scheme enumeration (part of every request) and the result types -- the first bead
// (FirstBeadData, internal to the schemes) and the whole molecule (GrowResult, RetraceResult).
//
// Every weight is stored ONLY as its natural logarithm, accumulated per growth step so it stays exact
// where the raw product underflows. The weight of a long chain is a product of hundreds of per-step
// factors of order exp(-beta u); for a polymer it drops below the smallest double (~1e-308) and a
// raw weight flushes to zero, turning every acceptance ratio W_new/W_old into 0/0 (NaN, silent
// reject) or x/0 (inf, unconditional accept). Acceptance rules must therefore be written in log
// space: exp(logRosenbluthWeight_new - logRosenbluthWeight_old), and corrections (dual cut-off,
// Ewald) are applied with 'multiplyRosenbluthWeight'. The raw weight is available through
// 'rosenbluthWeight()' for reporting and averaging (Widom) only.

export namespace CBMC
{
/// How the first bead of a molecule is placed (grow) or weighted (retrace).
enum class FirstBeadScheme : std::size_t
{
  /// 'numberOfFirstBeadPositions' uniformly random positions in the box; weight = sum of Boltzmann
  /// factors / number of positions. Insertion and deletion.
  MultipleFirstBead = 0,
  /// The multiple-first-bead reinsertion of Esselink et al.: as above against a background without
  /// the molecule itself; the grow retains the partial weight 'GrowResult::firstBeadStoredR', which
  /// the retrace of the old configuration needs ('RetraceRequest::storedR').
  Reinsertion = 1,
  /// A single trial at 'GrowRequest::firstBeadPosition', its Boltzmann factor as weight. Identity
  /// change: the new molecule takes the position of the old one.
  Pinned = 2,
  /// A single trial at 'GrowRequest::firstBeadPosition' with weight one: the caller sampled the
  /// position and accounts for its bias (distance-biased pair and group insertion).
  Fixed = 3,
  /// No first-bead stage: the beads in 'GrowRequest::beadsAlreadyPlaced' keep their positions from
  /// the given molecule and only the remaining beads are (re)grown. Partial reinsertion.
  AlreadyPlaced = 4,
};

/// The name of a scheme, for error messages.
[[nodiscard]] constexpr const char *firstBeadSchemeName(FirstBeadScheme scheme) noexcept
{
  switch (scheme)
  {
    case FirstBeadScheme::MultipleFirstBead:
      return "MultipleFirstBead";
    case FirstBeadScheme::Reinsertion:
      return "Reinsertion";
    case FirstBeadScheme::Pinned:
      return "Pinned";
    case FirstBeadScheme::Fixed:
      return "Fixed";
    case FirstBeadScheme::AlreadyPlaced:
      return "AlreadyPlaced";
  }
  return "?";
}

/// One non-overlapping trial position of the first bead and its external energy (the candidates the
/// multiple-first-bead schemes select among).
struct FirstBeadTrial
{
  Atom position;
  RunningEnergy energy;
};

/// Result of placing (or retracing) the first bead of a molecule with one of the first-bead schemes.
struct FirstBeadData
{
  Atom atom;
  RunningEnergy energies;
  /// Natural logarithm of the first-bead Rosenbluth weight.
  double logRosenbluthWeight;
  /// Retained partial weight of the multiple-first-bead reinsertion scheme (Esselink et al., 'r' in
  /// Eq. 16-18): the Rosenbluth sum minus the Boltzmann factor of the selected trial. Deliberately
  /// linear -- the retrace adds it to a Boltzmann factor, w(o) = exp(-beta u(o)) + r; the first-bead
  /// weight is a single-bead average of Boltzmann factors, order one, so no underflow is possible.
  /// Zero for the other schemes.
  double storedR;

  FirstBeadData() noexcept = delete;
  FirstBeadData(Atom atom, RunningEnergy energies, double logRosenbluthWeight, double storedR) noexcept
      : atom(atom), energies(energies), logRosenbluthWeight(logRosenbluthWeight), storedR(storedR)
  {
  }

  /// The first-bead Rosenbluth weight itself, exp(logRosenbluthWeight); for reporting only.
  [[nodiscard]] double rosenbluthWeight() const noexcept { return std::exp(logRosenbluthWeight); }
};

/// Result of growing a molecule (or the remainder of a molecule) with CBMC or recoil growth.
struct GrowResult
{
  Molecule molecule;
  std::vector<Atom> atoms;  ///< All atoms of the grown molecule.
  RunningEnergy energies;
  /// Natural logarithm of the Rosenbluth weight; the single stored representation (see above).
  double logRosenbluthWeight;
  /// Retained partial first-bead weight of the multiple-first-bead reinsertion scheme
  /// ('FirstBeadScheme::Reinsertion'; Esselink et al., 'r' in Eq. 16-18). The caller hands it to the
  /// retrace of the old configuration through 'RetraceRequest::storedR'. Zero for every other
  /// first-bead scheme.
  double firstBeadStoredR;

  GrowResult() : molecule(), atoms(), energies(), logRosenbluthWeight(), firstBeadStoredR() {}

  GrowResult(const Molecule &molecule, std::vector<Atom> atoms, RunningEnergy energies, double logRosenbluthWeight,
             double firstBeadStoredR = 0.0) noexcept
      : molecule(molecule),
        atoms(std::move(atoms)),
        energies(energies),
        logRosenbluthWeight(logRosenbluthWeight),
        firstBeadStoredR(firstBeadStoredR)
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
/// 'GrowResult'.
struct RetraceResult
{
  RunningEnergy energies;
  /// Natural logarithm of the Rosenbluth weight; the single stored representation.
  double logRosenbluthWeight;

  RetraceResult() : energies(), logRosenbluthWeight() {}

  RetraceResult(RunningEnergy energies, double logRosenbluthWeight) noexcept
      : energies(energies), logRosenbluthWeight(logRosenbluthWeight)
  {
  }

  /// The Rosenbluth weight itself, exp(logRosenbluthWeight); for reporting and averaging only.
  [[nodiscard]] double rosenbluthWeight() const noexcept { return std::exp(logRosenbluthWeight); }

  /// Multiplies the Rosenbluth weight by exp(logFactor).
  void multiplyRosenbluthWeight(double logFactor) noexcept { logRosenbluthWeight += logFactor; }
};
}  // namespace CBMC
