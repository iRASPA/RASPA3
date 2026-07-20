module;

export module cbmc_move_statistics;

import std;

import archive;
import move_statistics;

export struct CBMCMoveStatistics
{
  std::uint64_t versionNumber{3};

  MoveStatistics<double> bondLengthChange{
      .maxChange = 0.3, .lowerLimit = 0.01, .upperLimit = 0.5};
  MoveStatistics<double> bendAngleChange{
      .maxChange = 0.3, .lowerLimit = 0.01, .upperLimit = 0.5};
  MoveStatistics<double> conePositionChange{
      .maxChange = 0.3, .lowerLimit = 0.01, .upperLimit = 0.5};
  // Internal ring-closure Monte-Carlo step sizes: 'ringDisplacementChange' is the maximum per-atom
  // and rigid-fragment translation (Angstrom), 'ringRotationChange' the maximum whole-ring tilt and
  // rigid-fragment rotation angle (radians). Adapted towards the target acceptance by 'optimize'.
  MoveStatistics<double> ringDisplacementChange{
      .maxChange = 0.1, .lowerLimit = 0.01, .upperLimit = 1.0};
  MoveStatistics<double> ringRotationChange{
      .maxChange = 0.15, .lowerLimit = 0.01, .upperLimit = std::numbers::pi};
  // The conformer-hopping crankshaft is deliberately a full-range (uniform +/- pi) rotation -- its job
  // is to cross barriers between ring conformers, so its angle is never adapted; the statistics exist
  // to report the acceptance rate (and to let the attempt probability, a force-field option, be tuned).
  MoveStatistics<double> ringCrankshaftMove{
      .maxChange = std::numbers::pi, .lowerLimit = std::numbers::pi, .upperLimit = std::numbers::pi};

  // Adapts the ring-closure Monte-Carlo step sizes towards their target acceptance ratios from the
  // counters accumulated since the previous call. Invoked from 'System::optimizeMCMoves'. The
  // flexible-bead step sizes (bond-length / bend-angle / cone-position) are deliberately left fixed
  // at their defaults, matching the historical behavior of the flexible-bead sampler, and the ring
  // crankshaft angle is fixed at full range (see above).
  void optimize()
  {
    ringDisplacementChange.optimizeAcceptance();
    ringRotationChange.optimizeAcceptance();
  }

  const std::string writeMCMoveStatistics() const;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const CBMCMoveStatistics& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, CBMCMoveStatistics& p);
};
