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
import cbmc_growth_plan;

namespace Constants = CBMC::Constants;

// Rigid-body tilt: samples the junction-bend tilt of a rigid fragment hinged on the anchor with a
// rigid-rotation Metropolis MC. Carries no Rosenbluth weight. Ported from the former
// 'generateRigidUnitOrientationMonteCarloScheme'.
std::vector<Atom> CBMC::generateRigidTilt(RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead,
                                          double beta, const Component &component,
                                          const std::vector<Atom> &chainAtoms, const GrowStep &step)
{
  const std::optional<std::size_t> previousBead = step.previousBead;
  const std::size_t currentBead = step.currentBead;
  const std::vector<std::size_t> &nextBeads = step.nextBeads;
  const Potentials::IntraMolecularPotentials &intra = step.intra;

  double3 anchor_reference = component.atoms[currentBead].position;
  double3 anchor_position = chainAtoms[currentBead].position;
  std::vector<Atom> chain_atoms(chainAtoms.begin(), chainAtoms.end());

  auto placeWithRotation = [&](const double3x3 &rotation)
  {
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      double3 offset = component.atoms[nextBeads[k]].position - anchor_reference;
      chain_atoms[nextBeads[k]].position = anchor_position + rotation * offset;
    }
  };

  // Seed group or no junction bends: every tilt equally likely, return a uniform orientation.
  if (!previousBead.has_value() || intra.bends.empty())
  {
    placeWithRotation(random.randomRotationMatrix());
    std::vector<Atom> result(nextBeads.size());
    for (std::size_t k = 0; k != nextBeads.size(); ++k) result[k] = chain_atoms[nextBeads[k]];
    return result;
  }

  double3 last_bond_vector = (chainAtoms[previousBead.value()].position - anchor_position).normalized();

  // The inner bead and the junction bend are step constants, looked up when the plan was built.
  const std::size_t inner = step.rigidTilt.innerBead;
  const double bend_angle = step.rigidTilt.junctionBend.has_value()
                                ? step.rigidTilt.junctionBend->generateBendAngle(random, beta)
                                : Constants::defaultRigidJunctionBendAngle;

  double3 target_direction = random.randomVectorOnCone(last_bond_vector, bend_angle);
  double3 body_inner = (component.atoms[inner].position - anchor_reference).normalized();
  double3x3 alignment = double3x3::computeRotationMatrix(body_inner, target_direction);

  std::vector<double3> aligned_offsets(nextBeads.size());
  for (std::size_t k = 0; k != nextBeads.size(); ++k)
  {
    aligned_offsets[k] = alignment * (component.atoms[nextBeads[k]].position - anchor_reference);
  }

  constexpr std::size_t numberOfRollAngles = Constants::rigidTiltRollGridPoints;
  std::vector<double> logRollBoltzmannFactors(numberOfRollAngles);
  double roll_offset = 2.0 * std::numbers::pi * random.uniform();
  for (std::size_t r = 0; r != numberOfRollAngles; ++r)
  {
    double roll_angle = roll_offset + 2.0 * std::numbers::pi * static_cast<double>(r) / numberOfRollAngles;
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      chain_atoms[nextBeads[k]].position =
          anchor_position + target_direction.rotateAroundAxis(aligned_offsets[k], roll_angle);
    }
    logRollBoltzmannFactors[r] = -beta * intra.calculateBendSmallMCEnergies(chain_atoms);
  }

  std::size_t selected_roll = CBMC::selectTrialPosition(random, logRollBoltzmannFactors);
  double roll_angle = roll_offset + 2.0 * std::numbers::pi * static_cast<double>(selected_roll) / numberOfRollAngles;
  for (std::size_t k = 0; k != nextBeads.size(); ++k)
  {
    chain_atoms[nextBeads[k]].position =
        anchor_position + target_direction.rotateAroundAxis(aligned_offsets[k], roll_angle);
  }

  double current_bend_energy = intra.calculateBendSmallMCEnergies(chain_atoms);

  // Adaptive step size: the maximum rotation angle is read from the anchor bead's CBMC statistics
  // and adapted towards the target acceptance ratio between sweeps by 'System::optimizeMCMoves',
  // exactly like the ring-closure step sizes. The tilt carries no Rosenbluth weight, so the step
  // size affects only sampling efficiency, not detailed balance.
  MoveStatistics<double> &rotationStats = component.cbmc_moves_statistics[currentBead].rigidTiltRotationChange;
  const double maximumRotationAngle = rotationStats.maxChange;
  std::size_t number_of_trials = 2 * numberOfTrialMovesPerOpenBead * nextBeads.size();
  std::vector<double3> saved_positions(nextBeads.size());

  for (std::size_t trial = 0; trial != number_of_trials; ++trial)
  {
    rotationStats.counts += 1.0;
    rotationStats.totalCounts += 1.0;
    rotationStats.constructed += 1.0;
    rotationStats.totalConstructed += 1.0;

    double3 axis = random.randomVectorOnUnitSphere();
    double angle = (2.0 * random.uniform() - 1.0) * maximumRotationAngle;
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      saved_positions[k] = chain_atoms[nextBeads[k]].position;
      chain_atoms[nextBeads[k]].position =
          anchor_position + axis.rotateAroundAxis(saved_positions[k] - anchor_position, angle);
    }
    double trial_bend_energy = intra.calculateBendSmallMCEnergies(chain_atoms);
    if (random.uniform() < std::exp(-beta * (trial_bend_energy - current_bend_energy)))
    {
      current_bend_energy = trial_bend_energy;
      rotationStats.accepted += 1.0;
      rotationStats.totalAccepted += 1.0;
    }
    else
    {
      for (std::size_t k = 0; k != nextBeads.size(); ++k) chain_atoms[nextBeads[k]].position = saved_positions[k];
    }
  }

  std::vector<Atom> result(nextBeads.size());
  for (std::size_t k = 0; k != nextBeads.size(); ++k) result[k] = chain_atoms[nextBeads[k]];
  return result;
}
