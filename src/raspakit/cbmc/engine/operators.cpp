module;

module cbmc_operators;

import std;

import randomnumbers;
import forcefield;
import component;
import atom;
import double3;
import bond_potential;
import cbmc_constants;
import cbmc_growth_plan;
import cbmc_torsion_selection;
import cbmc_flexible_base;
import cbmc_rigid_tilt;
import cbmc_ring_closure;

namespace
{
namespace Constants = CBMC::Constants;

// Flexible single-bond seed: a Boltzmann bond length in a uniformly random direction.
Atom placeFlexibleSeedBead(RandomNumber &random, double beta, const std::vector<Atom> &chainAtoms,
                           const CBMC::GrowStep &step)
{
  const BondPotential *bond = step.intra.bonds.empty() ? nullptr : &step.intra.bonds.front();
  double bond_length = bond ? bond->generateBondLength(random, beta) : Constants::defaultBondLength;
  double3 unit_vector = random.randomVectorOnUnitSphere();
  Atom trial_atom = chainAtoms[step.nextBeads[0]];
  trial_atom.position = chainAtoms[step.currentBead].position + bond_length * unit_vector;
  return trial_atom;
}

// Base-conformation dispatch shared by grow / retrace / recoil: ring closure and rigid fragments keep
// their internal Monte-Carlo samplers (no clamp weight), every flexible step draws its base from the
// exact bonded sampler.
CBMC::FlexibleBase sampleBaseConformation(RandomNumber &random, const ForceField &forceField, double beta,
                                          const Component &component, std::vector<Atom> &chainAtoms,
                                          const CBMC::GrowStep &step)
{
  if (step.kind == CBMC::GrowStep::Kind::CloseRing)
  {
    return {CBMC::generateRingConformation(random, forceField, beta, component, chainAtoms, step), 1.0};
  }
  if (step.rigidBody)
  {
    return {CBMC::generateRigidTilt(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, chainAtoms,
                                    step),
            1.0};
  }
  return CBMC::sampleExactFlexibleBase(random, beta, component, chainAtoms, step);
}

// The previous-current axis of a step with a junction (a degenerate axis falls back to z). The one
// definition of the torsion-spin axis for grow, retrace, and recoil: the three must spin about the
// same axis for their torsion weights to be comparable.
double3 junctionAxis(const std::vector<Atom> &chainAtoms, const CBMC::GrowStep &step)
{
  double3 last_bond_vector =
      chainAtoms[step.previousBead.value()].position - chainAtoms[step.currentBead].position;
  if (last_bond_vector.length() < Constants::degenerateAxisLength) last_bond_vector = double3{0.0, 0.0, 1.0};
  return last_bond_vector.normalized();
}
}  // namespace

bool CBMC::stepHandlesUnsampledInternalTerms(const GrowStep &step) { return step.flexibleAttach; }

std::vector<CBMC::StepTrial> CBMC::generateGrowTrials(RandomNumber &random, const ForceField &forceField, double beta,
                                                      const Component &component, std::vector<Atom> &chainAtoms,
                                                      const GrowStep &step, std::size_t numberOfTrialDirections)
{
  std::vector<StepTrial> trials(numberOfTrialDirections);

  // Seed step: no orientational reference exists yet, every orientation is equally likely.
  if (!step.previousBead.has_value())
  {
    if (step.rigidBody)
    {
      // Rigid seed: each direction is an independent uniform orientation of the body about the anchor.
      for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {generateRigidTilt(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, chainAtoms,
                                       step),
                     1.0};
      }
      return trials;
    }
    if (step.kind == GrowStep::Kind::CloseRing)
    {
      // Ring seed: sample one internal conformation and rigidly rotate it for the other directions.
      std::vector<Atom> base = generateRingConformation(random, forceField, beta, component, chainAtoms, step);
      for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {randomlyOrientRing(random, base, chainAtoms[step.currentBead].position), 1.0};
      }
      trials[0] = {std::move(base), 1.0};
      return trials;
    }

    for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
    {
      trials[i] = {{placeFlexibleSeedBead(random, beta, chainAtoms, step)}, 1.0};
    }
    return trials;
  }

  // Attach / ring-closure with a junction: one shared base conformation, one torsion spin per
  // direction about the junction bond. The base's clamp-excess weight is shared by every direction.
  FlexibleBase base = sampleBaseConformation(random, forceField, beta, component, chainAtoms, step);
  const double3 last_bond_vector = junctionAxis(chainAtoms, step);

  for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
  {
    TorsionOrientation torsion =
        selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta, chainAtoms,
                                 base.nextBeadAtoms, step, last_bond_vector, false);
    trials[i] = {std::move(torsion.positions), torsion.rosenbluthWeight * base.clampWeight};
  }
  return trials;
}

std::vector<CBMC::StepTrial> CBMC::generateRetraceTrials(RandomNumber &random, const ForceField &forceField,
                                                         double beta, const Component &component,
                                                         std::vector<Atom> &chainAtoms, const GrowStep &step,
                                                         std::size_t numberOfTrialDirections)
{
  std::vector<StepTrial> trials(numberOfTrialDirections);

  // The old positions of the step's next-beads (trial direction 0).
  std::vector<Atom> old_orientation(step.nextBeads.size());
  for (std::size_t k = 0; k != step.nextBeads.size(); ++k) old_orientation[k] = chainAtoms[step.nextBeads[k]];

  if (!step.previousBead.has_value())
  {
    if (step.rigidBody)
    {
      for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {generateRigidTilt(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, chainAtoms,
                                       step),
                     1.0};
      }
    }
    else if (step.kind == GrowStep::Kind::CloseRing)
    {
      for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {randomlyOrientRing(random, old_orientation, chainAtoms[step.currentBead].position), 1.0};
      }
    }
    else
    {
      for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {{placeFlexibleSeedBead(random, beta, chainAtoms, step)}, 1.0};
      }
    }
    trials[0] = {std::move(old_orientation), 1.0};
    return trials;
  }

  // Attach / ring-closure with a junction. The old orientation is the shared torsion base for every
  // trial direction, pinned as torsion trial 0 of the first trial direction; its clamp-excess weight
  // (the old positions ARE the base) is shared by every direction, mirroring the grow side.
  const double base_clamp_weight = flexibleBaseClampWeight(beta, component, step, chainAtoms);
  const double3 last_bond_vector = junctionAxis(chainAtoms, step);

  for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
  {
    TorsionOrientation torsion = selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta,
                                                          chainAtoms, old_orientation, step, last_bond_vector, i == 0);
    trials[i] = {i == 0 ? old_orientation : std::move(torsion.positions),
                 torsion.rosenbluthWeight * base_clamp_weight};
  }
  return trials;
}

CBMC::StepTrial CBMC::generateRecoilTrial(RandomNumber &random, const ForceField &forceField, double beta,
                                          const Component &component, std::vector<Atom> &contextAtoms,
                                          const GrowStep &step)
{
  if (!step.previousBead.has_value())
  {
    if (step.rigidBody)
    {
      return {generateRigidTilt(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, contextAtoms,
                                step),
              1.0};
    }
    if (step.kind == GrowStep::Kind::CloseRing)
    {
      return {generateRingConformation(random, forceField, beta, component, contextAtoms, step), 1.0};
    }
    return {{placeFlexibleSeedBead(random, beta, contextAtoms, step)}, 1.0};
  }

  // Every flexible step (single bead or branch) draws its base from the exact bonded sampler; ring
  // closure and rigid fragments keep their internal Monte Carlo. The spin about the junction bond is
  // then always torsion-selected -- also when the trial is a feeler bead (see the interface note).
  FlexibleBase base = sampleBaseConformation(random, forceField, beta, component, contextAtoms, step);
  const double3 last_bond_vector = junctionAxis(contextAtoms, step);

  TorsionOrientation torsion =
      selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta, contextAtoms,
                               base.nextBeadAtoms, step, last_bond_vector, false);
  return {std::move(torsion.positions), torsion.rosenbluthWeight * base.clampWeight};
}

double CBMC::oldConfigurationTorsionWeight(RandomNumber &random, const ForceField &forceField, double beta,
                                           const Component &component, std::vector<Atom> &oldAtoms,
                                           const GrowStep &step)
{
  if (!step.previousBead.has_value()) return 1.0;

  std::vector<Atom> old_orientation(step.nextBeads.size());
  for (std::size_t k = 0; k != step.nextBeads.size(); ++k) old_orientation[k] = oldAtoms[step.nextBeads[k]];

  const double base_clamp_weight = flexibleBaseClampWeight(beta, component, step, oldAtoms);
  const double3 last_bond_vector = junctionAxis(oldAtoms, step);

  TorsionOrientation torsion = selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta,
                                                        oldAtoms, old_orientation, step, last_bond_vector, true);
  return torsion.rosenbluthWeight * base_clamp_weight;
}
