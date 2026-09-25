module;

module cbmc_flexible_base;

import std;

import randomnumbers;
import atom;
import double3;
import component;
import bond_potential;
import bend_potential;
import chiral_center;
import intra_molecular_potentials;
import running_energy;
import cbmc_util;
import cbmc_constants;
import cbmc_growth_plan;

namespace
{
namespace Constants = CBMC::Constants;

// Minimum of a bend potential over [0, pi]: the envelope offset that makes the sibling-bend
// rejection acceptance exp(-beta (u - u_min)) valid. Analytic for the harmonic forms; a fine grid
// otherwise, lowered by a small slack because a grid minimum can only overestimate the true minimum.
// Evaluated only at setup (once per step signature), so no memo is needed.
double bendMinimumEnergy(const BendPotential &bend)
{
  if ((bend.type == BendType::Harmonic || bend.type == BendType::CoreShell) && bend.parameters[1] >= 0.0 &&
      bend.parameters[1] <= std::numbers::pi)
  {
    return 0.0;
  }

  constexpr std::size_t numberOfGridPoints = Constants::bendMinimumGridPoints;
  const double3 posA{1.0, 0.0, 0.0};
  const double3 posB{0.0, 0.0, 0.0};
  double minimum = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i != numberOfGridPoints; ++i)
  {
    double theta = std::numbers::pi * static_cast<double>(i) / static_cast<double>(numberOfGridPoints - 1);
    double3 posC{std::cos(theta), std::sin(theta), 0.0};
    minimum = std::min(minimum, bend.calculateEnergy(posA, posB, posC, std::nullopt));
  }
  return minimum - Constants::bendMinimumSlack;  // the grid can only overestimate the true minimum
}

// The frozen base-coupling constants of a flexible attach step at 'beta'. Normally a plain lookup:
// the component prepared its plans for the system temperature at setup ('Component::prepareGrowthPlans')
// and every cached step carries its constants. The fall-back covers a component that was never
// prepared for this 'beta' (e.g. a unit test driving the operators directly): prepare it now, which
// fills the cached steps in place; a step that is not owned by the component's cache (a copied plan)
// is estimated directly -- the estimate is deterministic per signature, so it equals the memoised one.
CBMC::BaseCouplingConstants baseCouplingConstantsOf(double beta, const Component &component,
                                                    const CBMC::GrowStep &step)
{
  if (!step.hasBaseCoupling) return {0.0, 0.0};
  if (step.baseCouplingConstants.has_value() && component.growthPlanBeta == beta)
  {
    return step.baseCouplingConstants.value();
  }
  component.prepareGrowthPlans(beta);
  if (step.baseCouplingConstants.has_value() && component.growthPlanBeta == beta)
  {
    return step.baseCouplingConstants.value();
  }
  return CBMC::estimateBaseCouplingConstants(beta, component.atoms.size(), step);
}
}  // namespace

// ---------------------------------------------------------------------------------------------------
// The sampler (see the interface documentation).
// ---------------------------------------------------------------------------------------------------
CBMC::FlexibleBase CBMC::sampleExactFlexibleBase(RandomNumber &random, double beta, const Component &component,
                                                 const std::vector<Atom> &moleculeAtoms, const GrowStep &step)
{
  const std::size_t previousBead = step.previousBead.value();
  const std::size_t currentBead = step.currentBead;
  const std::vector<std::size_t> &nextBeads = step.nextBeads;
  const std::size_t numberOfNextBeads = nextBeads.size();

  std::vector<Atom> chain_atoms(moleculeAtoms.begin(), moleculeAtoms.end());

  double3 last_bond_vector = chain_atoms[previousBead].position - chain_atoms[currentBead].position;
  if (last_bond_vector.length() < Constants::degenerateAxisLength) last_bond_vector = double3{0.0, 0.0, 1.0};
  last_bond_vector = last_bond_vector.normalized();

  const bool has_coupling = step.hasBaseCoupling;
  const BaseCouplingConstants coupling_constants =
      has_coupling ? baseCouplingConstantsOf(beta, component, step) : BaseCouplingConstants{0.0, 0.0};

  for (std::size_t attempt = 0; attempt != Constants::baseSamplerMaximumAttempts; ++attempt)
  {
    for (std::size_t i = 0; i != numberOfNextBeads; ++i)
    {
      const std::optional<BondPotential> &bond = step.nextBeadBonds[i];
      const std::optional<BendPotential> &anchor = step.nextBeadAnchorBends[i];
      double bond_length = bond.has_value() ? bond->generateBondLength(random, beta) : Constants::defaultBondLength;
      double3 direction = anchor.has_value()
                              ? random.randomVectorOnCone(last_bond_vector, anchor->generateBendAngle(random, beta))
                              : random.randomVectorOnUnitSphere();
      chain_atoms[nextBeads[i]].position = chain_atoms[currentBead].position + bond_length * direction;
    }

    double clamp_weight = 1.0;
    if (has_coupling)
    {
      const double coupling_energy = baseCouplingEnergy(step, chain_atoms);
      const double boltzmann_excess = std::exp(-beta * (coupling_energy - coupling_constants.referenceEnergy));
      if (random.uniform() > boltzmann_excess) continue;
      clamp_weight = std::max(1.0, boltzmann_excess);
    }

    bool parity_ok = true;
    for (const ChiralCenter &center : step.determinedChiralCenters)
    {
      double reference_sign = center.type == ChiralCenter::Chirality::R_Chiral ? 1.0 : -1.0;
      parity_ok = parity_ok && (reference_sign * chiralSignedVolume(center.ids, chain_atoms) > 0.0);
    }
    if (!parity_ok) continue;

    std::vector<Atom> next_bead_atoms(numberOfNextBeads);
    for (std::size_t i = 0; i != numberOfNextBeads; ++i) next_bead_atoms[i] = chain_atoms[nextBeads[i]];
    return {std::move(next_bead_atoms), clamp_weight};
  }
  throw std::runtime_error(
      "CBMC: the exact base-conformation sampler exceeded its rejection budget (pathologically stiff "
      "base couplings or a nearly unsatisfiable chiral constraint)\n");
}

double CBMC::flexibleBaseClampWeight(double beta, const Component &component, const GrowStep &step,
                                     std::span<const Atom> atoms)
{
  if (!step.flexibleAttach || !step.hasBaseCoupling) return 1.0;

  const BaseCouplingConstants constants = baseCouplingConstantsOf(beta, component, step);
  const double couplingEnergy = baseCouplingEnergy(step, atoms);
  return std::max(1.0, std::exp(-beta * (couplingEnergy - constants.referenceEnergy)));
}

// ---------------------------------------------------------------------------------------------------
// The base-sampler normalization of a plan: for each flexible attach step, the bead densities
// r^2 exp(-beta u_bond) and sin(theta) exp(-beta u_anchor) (uniform azimuth) integrate to the
// product of the potentials' one-dimensional normalizations, the coupling rejection (sibling bends
// plus base-routed unsampled terms) multiplies in <a> e^{-beta u_ref} (see BaseCouplingConstants),
// and each fully determined chiral center restricts the base to a sector of probability exactly one
// half. See the interface documentation for why reptation needs this.
// ---------------------------------------------------------------------------------------------------
double CBMC::logBaseSamplerNormalization(double beta, const Component &component, const std::vector<GrowStep> &plan)
{
  double logNormalization = 0.0;

  for (const GrowStep &step : plan)
  {
    // Rigid-body and ring-closure steps keep internal-MC samplers without a closed-form
    // normalization; reptation enforces congruent end plans for units containing them, so their
    // (equal) factors cancel and are skipped here. Seed steps do not occur in partial plans.
    if (!step.flexibleAttach) continue;

    for (std::size_t i = 0; i != step.nextBeads.size(); ++i)
    {
      const std::optional<BondPotential> &bond = step.nextBeadBonds[i];
      // A bond without potential is placed at the fixed default length: r^2 dr integrates to r^2.
      logNormalization += bond.has_value() ? bond->logBoltzmannVolumeNormalization(beta)
                                           : 2.0 * std::log(Constants::defaultBondLength);

      const std::optional<BendPotential> &anchorBend = step.nextBeadAnchorBends[i];
      logNormalization += anchorBend.has_value() ? anchorBend->logBoltzmannConeNormalization(beta)
                                                 : std::log(4.0 * std::numbers::pi);  // uniform sphere
    }

    if (step.hasBaseCoupling)
    {
      const BaseCouplingConstants constants = baseCouplingConstantsOf(beta, component, step);
      logNormalization += constants.logMeanClampedBoltzmann - beta * constants.referenceEnergy;
    }

    logNormalization += static_cast<double>(step.determinedChiralCenters.size()) * std::log(0.5);
  }

  return logNormalization;
}

// ---------------------------------------------------------------------------------------------------
// The coupling energy and its frozen constants.
// ---------------------------------------------------------------------------------------------------
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
