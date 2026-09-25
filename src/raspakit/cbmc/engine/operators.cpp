module;

module cbmc_operators;

import std;

import randomnumbers;
import cbmc_grow_context;
import component;
import atom;
import double3;
import bond_potential;
import cbmc_util;
import cbmc_constants;
import cbmc_grow_step;
import cbmc_torsion_selection;
import cbmc_flexible_base;
import cbmc_rigid_tilt;
import cbmc_ring_closure;

// The dispatch is written once, as two primitives that grow, retrace, and recoil growth share:
//
//   seedTrials(pinnedOld)   -- a seed step (no orientational reference): n independent uniformly
//                              oriented placements, or the old one pinned as trial 0 plus n-1 fresh.
//   spinTrials(base, pinOld) -- a step with a junction: one shared base conformation spun about the
//                              previous-current axis, one torsion selection per trial; on the retrace
//                              the old positions ARE the base and are pinned as spin 0 of trial 0.
//
// grow    = seedTrials(none)      | spinTrials(sampleBase(),   pin = false)
// retrace = seedTrials(old)       | spinTrials(oldBase(),      pin = true)
// recoil  = the grow with n = 1;  the old torsion weight of recoil = the retrace with n = 1.
//
// The random-number draw order of every primitive equals that of the former per-operator ladders, so
// seeded trajectories are unchanged.
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

// Base-conformation dispatch of a step with a junction: ring closure and rigid fragments keep their
// internal Monte-Carlo samplers (no clamp weight), every flexible step draws its base from the exact
// bonded sampler.
CBMC::FlexibleBase sampleBase(RandomNumber &random, const CBMC::GrowthSettings &settings, double beta,
                              const Component &component, std::vector<Atom> &chainAtoms, const CBMC::GrowStep &step)
{
  if (step.kind == CBMC::GrowStep::Kind::CloseRing)
  {
    return {CBMC::generateRingConformation(random, settings, beta, component, chainAtoms, step), 1.0};
  }
  if (step.rigidBody)
  {
    return {CBMC::generateRigidTilt(random, settings.numberOfTrialMovesPerOpenBead, beta, component, chainAtoms,
                                    step),
            1.0};
  }
  return CBMC::sampleExactFlexibleBase(random, beta, component, chainAtoms, step);
}

// The old positions as the base of a retrace: they carry the clamp-excess weight a fresh draw of
// exactly these positions would have carried (one for ring and rigid steps).
CBMC::FlexibleBase oldBase(double beta, const Component &component, const std::vector<Atom> &chainAtoms,
                           const CBMC::GrowStep &step)
{
  return {CBMC::stepBeadPositions(chainAtoms, step), CBMC::flexibleBaseClampWeight(beta, component, step, chainAtoms)};
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

// Seed step: no orientational reference exists yet, every orientation is equally likely and every
// torsion weight is one. With 'pinnedOld' the old orientation is trial 0 and only n-1 are generated.
//  - rigid seed: each direction an independent uniform orientation of the body about the anchor;
//  - ring seed: one internal conformation (sampled, or the old one), rigidly rotated for the others;
//  - flexible seed: a Boltzmann bond length in a uniformly random direction per direction.
std::vector<CBMC::StepTrial> seedTrials(RandomNumber &random, const CBMC::GrowthSettings &settings, double beta,
                                        const Component &component, std::vector<Atom> &chainAtoms,
                                        const CBMC::GrowStep &step, std::size_t numberOfTrialDirections,
                                        std::optional<std::vector<Atom>> pinnedOld)
{
  std::vector<CBMC::StepTrial> trials(numberOfTrialDirections);
  const std::size_t first = pinnedOld.has_value() ? 1 : 0;

  if (step.rigidBody)
  {
    for (std::size_t i = first; i != numberOfTrialDirections; ++i)
    {
      trials[i] = {CBMC::generateRigidTilt(random, settings.numberOfTrialMovesPerOpenBead, beta, component,
                                           chainAtoms, step),
                   1.0};
    }
    if (pinnedOld.has_value()) trials[0] = {std::move(pinnedOld.value()), 1.0};
    return trials;
  }

  if (step.kind == CBMC::GrowStep::Kind::CloseRing)
  {
    std::vector<Atom> base = pinnedOld.has_value()
                                 ? std::move(pinnedOld.value())
                                 : CBMC::generateRingConformation(random, settings, beta, component, chainAtoms, step);
    const double3 anchor = chainAtoms[step.currentBead].position;
    for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
    {
      trials[i] = {CBMC::randomlyOrientRing(random, base, anchor), 1.0};
    }
    trials[0] = {std::move(base), 1.0};
    return trials;
  }

  for (std::size_t i = first; i != numberOfTrialDirections; ++i)
  {
    trials[i] = {{placeFlexibleSeedBead(random, beta, chainAtoms, step)}, 1.0};
  }
  if (pinnedOld.has_value()) trials[0] = {std::move(pinnedOld.value()), 1.0};
  return trials;
}

// Step with a junction: one shared base, one torsion spin per direction about the junction bond. The
// base's clamp-excess weight rides on every direction. With 'pinOld' the base positions themselves are
// spin 0 of trial direction 0 (the retrace's old configuration).
std::vector<CBMC::StepTrial> spinTrials(RandomNumber &random, const CBMC::GrowthSettings &settings, double beta,
                                        std::vector<Atom> &chainAtoms, const CBMC::GrowStep &step,
                                        const CBMC::FlexibleBase &base, std::size_t numberOfTrialDirections,
                                        bool pinOld)
{
  std::vector<CBMC::StepTrial> trials(numberOfTrialDirections);
  const double3 axis = junctionAxis(chainAtoms, step);

  for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
  {
    const bool pinned = pinOld && i == 0;
    CBMC::TorsionOrientation torsion = CBMC::selectTorsionOrientation(
        random, settings.numberOfTorsionTrialDirections, beta, chainAtoms, base.nextBeadAtoms, step, axis, pinned);
    trials[i] = {pinned ? base.nextBeadAtoms : std::move(torsion.positions),
                 torsion.rosenbluthWeight * base.clampWeight};
  }
  return trials;
}
}  // namespace

bool CBMC::stepHandlesUnsampledInternalTerms(const GrowStep &step) { return step.flexibleAttach; }

std::vector<CBMC::StepTrial> CBMC::generateGrowTrials(RandomNumber &random, const GrowthSettings &settings, double beta,
                                                      const Component &component, std::vector<Atom> &chainAtoms,
                                                      const GrowStep &step, std::size_t numberOfTrialDirections)
{
  if (!step.previousBead.has_value())
  {
    return seedTrials(random, settings, beta, component, chainAtoms, step, numberOfTrialDirections, std::nullopt);
  }
  const FlexibleBase base = sampleBase(random, settings, beta, component, chainAtoms, step);
  return spinTrials(random, settings, beta, chainAtoms, step, base, numberOfTrialDirections, false);
}

std::vector<CBMC::StepTrial> CBMC::generateRetraceTrials(RandomNumber &random, const GrowthSettings &settings,
                                                         double beta, const Component &component,
                                                         std::vector<Atom> &chainAtoms, const GrowStep &step,
                                                         std::size_t numberOfTrialDirections)
{
  if (!step.previousBead.has_value())
  {
    return seedTrials(random, settings, beta, component, chainAtoms, step, numberOfTrialDirections,
                      CBMC::stepBeadPositions(chainAtoms, step));
  }
  const FlexibleBase base = oldBase(beta, component, chainAtoms, step);
  return spinTrials(random, settings, beta, chainAtoms, step, base, numberOfTrialDirections, true);
}

CBMC::StepTrial CBMC::generateRecoilTrial(RandomNumber &random, const GrowthSettings &settings, double beta,
                                          const Component &component, std::vector<Atom> &contextAtoms,
                                          const GrowStep &step)
{
  return std::move(generateGrowTrials(random, settings, beta, component, contextAtoms, step, 1).front());
}

double CBMC::oldConfigurationTorsionWeight(RandomNumber &random, const GrowthSettings &settings, double beta,
                                           const Component &component, std::vector<Atom> &oldAtoms,
                                           const GrowStep &step)
{
  return generateRetraceTrials(random, settings, beta, component, oldAtoms, step, 1).front().torsionWeight;
}
