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
import chiral_center;
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
//                            step size from the anchor bead's CBMC statistics; about the inner bond only
//                            when the junction bend is a FIXED or RIGID constraint.
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
        saved_(step.nextBeads.size()),
        holonomicJunction_(step.rigidTilt.junctionBend.has_value() &&
                           (step.rigidTilt.junctionBend->type == BendType::Fixed ||
                            step.rigidTilt.junctionBend->type == BendType::Rigid))
  {
  }

  // Empty when no orientation with the declared parity of the hinge's chiral centres was found.
  std::optional<std::vector<Atom>> sample()
  {
    const bool uniform = !step_.previousBead.has_value() || step_.rigidTilt.bends.empty();
    bool placed = false;
    for (std::size_t attempt = 0; attempt != Constants::ringParityMaximumReorientations && !placed; ++attempt)
    {
      if (uniform)
      {
        uniformOrientation();
        placed = parityPreserved();
      }
      else
      {
        const double3 coneDirection = alignToJunctionBend();
        placed = selectRoll(coneDirection);
      }
    }
    if (!placed) return std::nullopt;
    if (!uniform) relax();

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

  // Whether every chiral centre the tilt decides has its declared parity (that of the reference geometry). The
  // bonded model is achiral, so without this the body is placed on either side of the junction alike.
  bool parityPreserved() const
  {
    for (const ChiralCenter& center : step_.rigidTilt.chiralCenters)
    {
      // strictly the declared side, as the other parity tests of the growth: a centre exactly in its plane has
      // no parity (the reference geometry of a declared centre is never planar, see the component reader)
      if (!(CBMC::chiralSignedVolume(center.ids, chain_) * CBMC::chiralSignedVolume(center.ids, component_.atoms) >
            0.0))
        return false;
    }
    return true;
  }

  // the spin-invariant bends of the step (see 'GrowStep::RigidTiltData::bends')
  double bendEnergy()
  {
    double energy = 0.0;
    for (const BendPotential& bend : step_.rigidTilt.bends)
    {
      energy += bend.calculateEnergy(chain_[bend.identifiers[0]].position, chain_[bend.identifiers[1]].position,
                                     chain_[bend.identifiers[2]].position, std::nullopt);
    }
    return energy;
  }

  // Rosenbluth selection of the roll about the cone direction on a randomly offset grid.
  // Returns false when no roll on the grid has the declared parity.
  bool selectRoll(const double3 &coneDirection)
  {
    constexpr double deadLogFactor = -1.0e300;
    std::size_t numberOfDeadRolls = 0;
    constexpr std::size_t numberOfRollAngles = Constants::rigidTiltRollGridPoints;
    const double rollOffset = 2.0 * std::numbers::pi * random_.uniform();
    auto rollAngleOf = [&](std::size_t r)
    { return rollOffset + 2.0 * std::numbers::pi * static_cast<double>(r) / numberOfRollAngles; };

    std::vector<double> logRollBoltzmannFactors(numberOfRollAngles);
    for (std::size_t r = 0; r != numberOfRollAngles; ++r)
    {
      placeRoll(coneDirection, rollAngleOf(r));
      logRollBoltzmannFactors[r] = -beta_ * bendEnergy();
      if (!parityPreserved())
      {
        logRollBoltzmannFactors[r] = deadLogFactor;
        ++numberOfDeadRolls;
      }
    }
    if (numberOfDeadRolls == numberOfRollAngles) return false;

    const std::size_t selected = CBMC::selectTrialPosition(random_, logRollBoltzmannFactors).index;
    placeRoll(coneDirection, rollAngleOf(selected));
    currentEnergy_ = bendEnergy();
    return true;
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

    // A FIXED or RIGID junction bend has no energy to hold it (its energy is zero by construction): the body then
    // only rolls about its inner bond, which leaves the angle previous-current-inner where 'alignToJunctionBend'
    // put it. A rotation about a free axis took the body off that cone at no cost, so the constrained angle came
    // out uniform on the sphere. The rotation about the previous-current axis is the spin, selected afterwards.
    const double3 axis = holonomicJunction_
                             ? (chain_[step_.rigidTilt.innerBead].position - anchorPosition_).normalized()
                             : random_.randomVectorOnUnitSphere();
    const double angle = (2.0 * random_.uniform() - 1.0) * rotationStats_.maxChange;
    for (std::size_t k = 0; k != nextBeads_.size(); ++k)
    {
      saved_[k] = chain_[nextBeads_[k]].position;
      chain_[nextBeads_[k]].position = anchorPosition_ + axis.rotateAroundAxis(saved_[k] - anchorPosition_, angle);
    }
    const double trialEnergy = bendEnergy();
    if (parityPreserved() && random_.uniform() < std::exp(-beta_ * (trialEnergy - currentEnergy_)))
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
  /// The junction bend is a constraint (FIXED or RIGID): the relaxation keeps its angle.
  const bool holonomicJunction_;
};
}  // namespace

std::optional<std::vector<Atom>> CBMC::generateRigidTilt(RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead,
                                          double beta, const Component &component, std::vector<Atom> &chainAtoms,
                                          const GrowStep &step)
{
  return RigidTiltSampler(random, numberOfTrialMovesPerOpenBead, beta, component, chainAtoms, step).sample();
}
