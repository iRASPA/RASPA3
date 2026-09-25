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
}  // namespace CBMC
