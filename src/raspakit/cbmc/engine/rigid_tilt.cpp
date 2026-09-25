module;

module cbmc_rigid_tilt;

import std;

import atom;
import double3;
import double3x3;
import randomnumbers;
import component;
import move_statistics;
import bend_potential;
import intra_molecular_potentials;
import cbmc_util;
import cbmc_constants;
import cbmc_statistics;
import cbmc_grow_step;

namespace
{
namespace Constants = CBMC::Constants;

// Samples the orientation of a rigid-body fragment hinged on the step's anchor. Carries no
// Rosenbluth weight. The body's internal geometry is that of the component's reference atoms,
// exactly; only its orientation about the anchor is sampled:
//
//   uniformOrientation()  -- seed step or a junction without bends: every tilt equally likely.
//   alignToJunctionBend() -- the bend angle previous-current-inner is drawn from its Boltzmann
//                            distribution and the body is aligned so its inner atom lies on that cone;
//   selectRoll()          -- the roll about the cone direction is Rosenbluth-selected on a grid of
//                            the step's bend energies (the other junction bends);
//   relax()               -- a rigid-rotation Metropolis Monte-Carlo on the step's bends, adaptive
//                            step size from the anchor bead's CBMC statistics.
class RigidTiltSampler
{
 public:
  RigidTiltSampler(RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead, double beta,
                   const Component &component, std::vector<Atom> &chainAtoms, const CBMC::GrowStep &step)
      : random_(random),
        numberOfTrialMovesPerOpenBead_(numberOfTrialMovesPerOpenBead),
        beta_(beta),
        component_(component),
        step_(step),
        nextBeads_(step.nextBeads),
        chain_(chainAtoms),
        scratch_(chainAtoms, step),
        anchorReference_(component.atoms[step.currentBead].position),
        anchorPosition_(chainAtoms[step.currentBead].position),
        rotationStats_(component.cbmcMoveStatistics[step.currentBead].rigidTiltRotationChange),
        saved_(step.nextBeads.size())
  {
  }

  std::vector<Atom> sample()
  {
    if (!step_.previousBead.has_value() || step_.intra.bends.empty())
    {
      uniformOrientation();
    }
    else
    {
      const double3 coneDirection = alignToJunctionBend();
      selectRoll(coneDirection);
      relax();
    }

    std::vector<Atom> result(nextBeads_.size());
    for (std::size_t k = 0; k != nextBeads_.size(); ++k) result[k] = chain_[nextBeads_[k]];
    return result;
  }

 private:
  // Body offset of a next bead relative to the anchor, in the reference geometry.
  double3 bodyOffset(std::size_t atom) const { return component_.atoms[atom].position - anchorReference_; }

  void uniformOrientation()
  {
    const double3x3 rotation = random_.randomRotationMatrix();
    for (std::size_t atom : nextBeads_) chain_[atom].position = anchorPosition_ + rotation * bodyOffset(atom);
  }

  // Draws the junction bend angle, aligns the body's inner atom to the resulting cone direction about
  // the previous-current axis, and stores the aligned body offsets. Returns the cone direction (the
  // roll axis).
  double3 alignToJunctionBend()
  {
    const double3 lastBondVector = (chain_[step_.previousBead.value()].position - anchorPosition_).normalized();

    // The inner bead and the junction bend are step constants, looked up when the plan was built.
    const std::size_t inner = step_.rigidTilt.innerBead;
    const double bendAngle = step_.rigidTilt.junctionBend.has_value()
                                 ? step_.rigidTilt.junctionBend->generateBendAngle(random_, beta_)
                                 : Constants::defaultRigidJunctionBendAngle;

    const double3 coneDirection = random_.randomVectorOnCone(lastBondVector, bendAngle);
    const double3x3 alignment = double3x3::computeRotationMatrix(bodyOffset(inner).normalized(), coneDirection);

    alignedOffsets_.resize(nextBeads_.size());
    for (std::size_t k = 0; k != nextBeads_.size(); ++k) alignedOffsets_[k] = alignment * bodyOffset(nextBeads_[k]);
    return coneDirection;
  }

  void placeRoll(const double3 &coneDirection, double rollAngle)
  {
    for (std::size_t k = 0; k != nextBeads_.size(); ++k)
    {
      chain_[nextBeads_[k]].position = anchorPosition_ + coneDirection.rotateAroundAxis(alignedOffsets_[k], rollAngle);
    }
  }

  double bendEnergy() { return step_.intra.calculateBendSmallMCEnergies(chain_); }

  // Rosenbluth selection of the roll about the cone direction on a randomly offset grid.
  void selectRoll(const double3 &coneDirection)
  {
    constexpr std::size_t numberOfRollAngles = Constants::rigidTiltRollGridPoints;
    const double rollOffset = 2.0 * std::numbers::pi * random_.uniform();
    auto rollAngleOf = [&](std::size_t r)
    { return rollOffset + 2.0 * std::numbers::pi * static_cast<double>(r) / numberOfRollAngles; };

    std::vector<double> logRollBoltzmannFactors(numberOfRollAngles);
    for (std::size_t r = 0; r != numberOfRollAngles; ++r)
    {
      placeRoll(coneDirection, rollAngleOf(r));
      logRollBoltzmannFactors[r] = -beta_ * bendEnergy();
    }

    const std::size_t selected = CBMC::selectTrialPosition(random_, logRollBoltzmannFactors);
    placeRoll(coneDirection, rollAngleOf(selected));
    currentEnergy_ = bendEnergy();
  }

  // Adaptive step size: the maximum rotation angle is read from the anchor bead's CBMC statistics and
  // adapted towards the target acceptance ratio between sweeps by 'System::optimizeMCMoves', exactly
  // like the ring-closure step sizes. The tilt carries no Rosenbluth weight, so the step size affects
  // only sampling efficiency, not detailed balance.
  void relax()
  {
    const std::size_t numberOfTrials = 2 * numberOfTrialMovesPerOpenBead_ * nextBeads_.size();
    for (std::size_t trial = 0; trial != numberOfTrials; ++trial) rotationMove();
  }

  // One symmetric rigid rotation of the body about a random axis through the anchor, Metropolis on
  // the step's bend energy.
  void rotationMove()
  {
    rotationStats_.counts += 1.0;
    rotationStats_.totalCounts += 1.0;
    rotationStats_.constructed += 1.0;
    rotationStats_.totalConstructed += 1.0;

    const double3 axis = random_.randomVectorOnUnitSphere();
    const double angle = (2.0 * random_.uniform() - 1.0) * rotationStats_.maxChange;
    for (std::size_t k = 0; k != nextBeads_.size(); ++k)
    {
      saved_[k] = chain_[nextBeads_[k]].position;
      chain_[nextBeads_[k]].position = anchorPosition_ + axis.rotateAroundAxis(saved_[k] - anchorPosition_, angle);
    }
    const double trialEnergy = bendEnergy();
    if (random_.uniform() < std::exp(-beta_ * (trialEnergy - currentEnergy_)))
    {
      currentEnergy_ = trialEnergy;
      rotationStats_.accepted += 1.0;
      rotationStats_.totalAccepted += 1.0;
    }
    else
    {
      for (std::size_t k = 0; k != nextBeads_.size(); ++k) chain_[nextBeads_[k]].position = saved_[k];
    }
  }

  RandomNumber &random_;
  const std::size_t numberOfTrialMovesPerOpenBead_;
  const double beta_;
  const Component &component_;
  const CBMC::GrowStep &step_;
  const std::vector<std::size_t> &nextBeads_;

  /// The chain the step grows in; the body's beads are edited in place and restored by 'scratch_' when
  /// the sampler is destroyed (see the scratch contract in cbmc_operators).
  std::vector<Atom> &chain_;
  const CBMC::ScratchBeads scratch_;
  const double3 anchorReference_;
  const double3 anchorPosition_;
  std::vector<double3> alignedOffsets_{};
  double currentEnergy_{0.0};

  MoveStatistics<double> &rotationStats_;
  std::vector<double3> saved_;
};
}  // namespace

std::vector<Atom> CBMC::generateRigidTilt(RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead,
                                          double beta, const Component &component, std::vector<Atom> &chainAtoms,
                                          const GrowStep &step)
{
  return RigidTiltSampler(random, numberOfTrialMovesPerOpenBead, beta, component, chainAtoms, step).sample();
}
