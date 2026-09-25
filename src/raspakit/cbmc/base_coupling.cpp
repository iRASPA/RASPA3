module;

module cbmc_base_coupling;

import std;

import randomnumbers;
import atom;
import double3;
import bond_potential;
import bend_potential;
import intra_molecular_potentials;
import cbmc_constants;
import cbmc_growth_plan;

// Minimum of a bend potential over [0, pi]: the envelope offset that makes the sibling-bend
// rejection acceptance exp(-beta (u - u_min)) valid. Analytic for the harmonic forms; a fine grid
// otherwise, lowered by a small slack because a grid minimum can only overestimate the true minimum.
// Evaluated only at setup (once per step signature), so no memo is needed.
static double bendMinimumEnergy(const BendPotential &bend)
{
  if ((bend.type == BendType::Harmonic || bend.type == BendType::CoreShell) && bend.parameters[1] >= 0.0 &&
      bend.parameters[1] <= std::numbers::pi)
  {
    return 0.0;
  }

  constexpr std::size_t numberOfGridPoints = CBMC::Constants::bendMinimumGridPoints;
  const double3 posA{1.0, 0.0, 0.0};
  const double3 posB{0.0, 0.0, 0.0};
  double minimum = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i != numberOfGridPoints; ++i)
  {
    double theta = std::numbers::pi * static_cast<double>(i) / static_cast<double>(numberOfGridPoints - 1);
    double3 posC{std::cos(theta), std::sin(theta), 0.0};
    minimum = std::min(minimum, bend.calculateEnergy(posA, posB, posC, std::nullopt));
  }
  return minimum - CBMC::Constants::bendMinimumSlack;  // the grid can only overestimate the true minimum
}

double CBMC::baseCouplingEnergy(const GrowStep &step, std::span<const Atom> atoms)
{
  double energy = 0.0;
  for (const BendPotential &bend : step.siblingBends)
  {
    energy += bend.calculateEnergy(atoms[bend.identifiers[0]].position, atoms[bend.identifiers[1]].position,
                                   atoms[bend.identifiers[2]].position, std::nullopt);
  }
  if (step.hasBaseCouplingTerms)
  {
    energy += step.baseCouplingTerms.computeInternalEnergiesNotSampledDuringGrowth(atoms).potentialEnergy();
  }
  return energy;
}

CBMC::BaseCouplingConstants CBMC::estimateBaseCouplingConstants(double beta, std::size_t numberOfAtoms,
                                                                const GrowStep &step)
{
  if (!step.flexibleAttach || !step.hasBaseCoupling) return {0.0, 0.0};

  const std::size_t previousBead = step.previousBead.value();
  const std::size_t currentBead = step.currentBead;

  RandomNumber random(Constants::baseCouplingEstimateSeed);  // frozen: one constant per signature
  const double3 axis{0.0, 0.0, 1.0};

  // Evaluation frame: current bead at the origin, previous bead at unit distance along the axis (the
  // classification guarantees no base term depends on the placed previous-current distance).
  std::vector<Atom> atoms(numberOfAtoms);
  atoms[currentBead].position = double3{0.0, 0.0, 0.0};
  atoms[previousBead].position = axis;

  const auto drawIndependentBase = [&]()
  {
    for (std::size_t i = 0; i != step.nextBeads.size(); ++i)
    {
      const std::optional<BondPotential> &bond = step.nextBeadBonds[i];
      const double bondLength =
          bond.has_value() ? bond->generateBondLength(random, beta) : Constants::defaultBondLength;
      const std::optional<BendPotential> &anchor = step.nextBeadAnchorBends[i];
      const double3 direction = anchor.has_value()
                                    ? random.randomVectorOnCone(axis, anchor->generateBendAngle(random, beta))
                                    : random.randomVectorOnUnitSphere();
      atoms[step.nextBeads[i]].position = bondLength * direction;
    }
    return baseCouplingEnergy(step, atoms);
  };

  // Reference energy: rigorous per-bend minima when only sibling bends couple; otherwise the
  // Monte-Carlo minimum with slack (exactness does not depend on it, see BaseCouplingConstants).
  double referenceEnergy = 0.0;
  if (!step.hasBaseCouplingTerms)
  {
    for (const BendPotential &bend : step.siblingBends) referenceEnergy += bendMinimumEnergy(bend);
  }
  else
  {
    double minimum = std::numeric_limits<double>::max();
    for (std::size_t s = 0; s != Constants::baseCouplingMinimumSearchSamples; ++s)
    {
      minimum = std::min(minimum, drawIndependentBase());
    }
    referenceEnergy = minimum - Constants::baseCouplingReferenceSlack;
  }

  double sum = 0.0;
  double sumOfSquares = 0.0;
  double numberOfSamples = 0.0;
  constexpr std::size_t batchSize = Constants::baseCouplingBatchSize;
  for (std::size_t batch = 0; batch != Constants::baseCouplingMaximumBatches; ++batch)
  {
    for (std::size_t s = 0; s != batchSize; ++s)
    {
      const double clampedBoltzmannFactor = std::min(1.0, std::exp(-beta * (drawIndependentBase() - referenceEnergy)));
      sum += clampedBoltzmannFactor;
      sumOfSquares += clampedBoltzmannFactor * clampedBoltzmannFactor;
    }
    numberOfSamples += static_cast<double>(batchSize);
    const double mean = sum / numberOfSamples;
    const double variance = std::max(0.0, sumOfSquares / numberOfSamples - mean * mean);
    if (std::sqrt(variance / numberOfSamples) < Constants::baseCouplingRelativeTolerance * mean) break;
  }

  return {referenceEnergy, std::log(sum / numberOfSamples)};
}

void CBMC::prepareBaseCouplingConstants(double beta, std::size_t numberOfAtoms, std::span<GrowStep> plan,
                                        std::map<std::string, BaseCouplingConstants> &memo)
{
  for (GrowStep &step : plan)
  {
    if (!step.flexibleAttach || !step.hasBaseCoupling)
    {
      step.baseCouplingConstants.reset();
      continue;
    }
    auto it = memo.find(step.baseCouplingSignature);
    if (it == memo.end())
    {
      it = memo.emplace(step.baseCouplingSignature, estimateBaseCouplingConstants(beta, numberOfAtoms, step)).first;
    }
    step.baseCouplingConstants = it->second;
  }
}
