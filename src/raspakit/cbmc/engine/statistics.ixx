module;

export module cbmc_statistics;

import std;

import archive;
import move_statistics;

export namespace CBMC
{
/// Acceptance counters and adaptive step sizes of the internal Monte-Carlo samplers of the operator
/// engine (ring closure and rigid tilt), kept per bead of a component ('Component::cbmcMoveStatistics',
/// indexed by the anchor bead of the step). Serialized into the restart file.
struct InternalMoveStatistics
{
  // Version 2 added 'rigidTiltRotationChange' (version-1 restart files are read with its default).
  // Version 3 dropped the bond-length / bend-angle / cone-position entries of the former internal
  // flexible-bead Monte-Carlo, which the exact base sampler replaced; older files still carry them and
  // the reader skips them.
  std::uint64_t versionNumber{3};

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
  // Rigid-body tilt Monte-Carlo: the maximum rotation angle (radians) of the small rigid rotations
  // that relax the junction bends of a hinged rigid fragment. Adapted like the ring step sizes.
  MoveStatistics<double> rigidTiltRotationChange{
      .maxChange = 0.15, .lowerLimit = 0.01, .upperLimit = std::numbers::pi};

  // Adapts the ring-closure and rigid-tilt Monte-Carlo step sizes towards their target acceptance
  // ratios from the counters accumulated since the previous call. Invoked from
  // 'System::optimizeMCMoves'. The ring crankshaft angle is fixed at full range (see above); the
  // flexible-bead base is sampled exactly and has no step size.
  void optimize()
  {
    ringDisplacementChange.optimizeAcceptance();
    ringRotationChange.optimizeAcceptance();
    rigidTiltRotationChange.optimizeAcceptance();
  }

  const std::string writeMCMoveStatistics() const;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const InternalMoveStatistics& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, InternalMoveStatistics& p);
};

/// Diagnostic counters of the recoil-growth chain scheme, kept per component
/// ('Component::recoilGrowthStatistics') and reported with the CBMC statistics. They tell a user tuning
/// the trial count 'k' and recoil length 'l' what the constructed/trial ratio of the move can not: how
/// the failed grows fail, and how many directions were available per step.
///  - A 'dead end' grow exhausted every direction within the recoil length and recoiled all the way
///    back to the start (more directions or a longer recoil helps).
///  - A 'discarded' grow committed to a direction (the chain got 'l' steps past it) and dead-ended
///    later, where recoiling is no longer allowed (a longer recoil length helps).
///  - The mean available directions m_i per step (grow side) measures how crowded the environment is
///    for this k: close to k means almost every direction is open, close to 1 means the weight is
///    dominated by single available directions and has a large variance.
/// Not serialized: the counters restart from zero after a restart, like the CPU timings.
struct RecoilGrowthStatistics
{
  double grows{0.0};
  double completed{0.0};
  double deadEnds{0.0};
  double discarded{0.0};
  double availableDirectionsSum{0.0};
  double growSteps{0.0};

  const std::string writeStatistics(std::size_t numberOfTrialDirections) const;
};
}  // namespace CBMC
