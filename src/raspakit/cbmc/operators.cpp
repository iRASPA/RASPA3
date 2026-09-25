module;

module cbmc_operators;

import std;

import randomnumbers;
import units;
import forcefield;
import component;
import atom;
import double3;
import double3x3;
import cbmc_util;
import cbmc_growth_plan;
import cbmc_move_statistics;
import move_statistics;
import intra_molecular_potentials;
import chiral_center;
import bond_potential;
import bend_potential;
import torsion_potential;
import urey_bradley_potential;
import inversion_bend_potential;
import out_of_plane_bend_potential;
import bond_bond_potential;
import bond_bend_potential;
import bend_bend_potential;
import bond_torsion_potential;
import bend_torsion_potential;
import running_energy;

namespace
{
// ---------------------------------------------------------------------------------------------------
// Numerical constants of the operator engine, gathered in one place. None of them affects the
// sampled distribution (the samplers are exact or Metropolis-correct for any value); they set
// fall-backs for under-specified topologies, the resolution of one-off numerical estimates, and the
// budgets of rejection loops.
// ---------------------------------------------------------------------------------------------------
namespace Constants
{
/// Bond length (Angstrom) used when a step's bond carries no potential (e.g. a connectivity entry
/// without a matching 'Bonds' term): a fixed C-C single-bond length. Also the value the normalization
/// integrates over in that case.
constexpr double defaultBondLength = 1.54;

/// Junction bend angle used by the rigid tilt when no previous-current-inner bend potential exists:
/// a trigonal 120 degrees.
constexpr double defaultRigidJunctionBendAngle = 120.0 * Units::DegreesToRadians;

/// Squared-length threshold below which a previous-current vector is treated as degenerate (the two
/// beads coincide) and replaced by the z axis.
constexpr double degenerateAxisLength = 1e-8;

/// Rejection budget of the exact flexible base sampler; exceeding it signals a pathologically stiff
/// coupling or an unsatisfiable chiral constraint and throws.
constexpr std::size_t baseSamplerMaximumAttempts = 1'000'000;

/// Random re-orientation attempts of the ring seed geometry to recover the declared parity of a
/// chiral centre before giving up on that seed.
constexpr std::size_t ringParityMaximumReorientations = 1000;

/// Grid over [0, pi] for the numerical minimum of a non-harmonic bend potential, and the slack
/// (Kelvin) subtracted because a grid minimum can only overestimate the true minimum.
constexpr std::size_t bendMinimumGridPoints = 8192;
constexpr double bendMinimumSlack = 0.01;

/// Fixed seed of the one-off Monte-Carlo estimate of the base-coupling constants: the estimate must
/// be one frozen number per step signature so grow and retrace share it.
constexpr std::size_t baseCouplingEstimateSeed = 1806;
/// Samples of the minimum search of the coupling energy, and the slack (Kelvin) below that minimum
/// for the rejection reference (see stepBaseCouplingConstants; exactness does not depend on it).
constexpr std::size_t baseCouplingMinimumSearchSamples = 1uz << 20;
constexpr double baseCouplingReferenceSlack = 1.0;
/// Batching of the clamped-acceptance mean: batches of 'batchSize' samples until the standard error
/// drops below 'relativeTolerance' of the mean, or 'maximumBatches' is reached.
constexpr std::size_t baseCouplingBatchSize = 1uz << 20;
constexpr std::size_t baseCouplingMaximumBatches = 32;
constexpr double baseCouplingRelativeTolerance = 0.005;

/// Roll angles about the junction bond tried when seating a hinged rigid body (a Rosenbluth selection
/// over this grid, randomly offset, seeds the tilt Monte-Carlo).
constexpr std::size_t rigidTiltRollGridPoints = 72;
}  // namespace Constants

// A bend is 'spin-variant' when it involves a placed atom other than the previous or current bead.
// Rotating the next-beads about the previous-current axis (the torsion spin) changes such bends,
// while bends among {previous, current, next-beads} transform rigidly (previous and current lie on
// the axis). Spin-variant bends are excluded from the base conformation (no Rosenbluth weight) and
// weighted exactly in the torsion selection instead.
bool isSpinVariantBend(const BendPotential &bend, std::optional<std::size_t> previousBead, std::size_t currentBead,
                       const std::vector<std::size_t> &nextBeads)
{
  for (std::size_t id : bend.identifiers)
  {
    if (id == currentBead) continue;
    if (previousBead.has_value() && id == previousBead.value()) continue;
    if (std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end()) continue;
    return true;
  }
  return false;
}

// A sibling bend couples two beads grown in the same step (next_i - current - next_j). It is
// spin-invariant (the whole branch rotates rigidly) and is imposed on the base conformation by
// rejection sampling; its coupling integral is part of the base normalization (see
// logBaseSamplerNormalization).
bool isSiblingBend(const BendPotential &bend, std::size_t currentBead, const std::vector<std::size_t> &nextBeads)
{
  if (bend.identifiers[1] != currentBead) return false;
  const bool endA_is_next = std::find(nextBeads.begin(), nextBeads.end(), bend.identifiers[0]) != nextBeads.end();
  const bool endC_is_next = std::find(nextBeads.begin(), nextBeads.end(), bend.identifiers[2]) != nextBeads.end();
  return endA_is_next && endC_is_next;
}

// ---------------------------------------------------------------------------------------------------
// Partition of a flexible attach step's internal terms that the classic growth stages do not sample
// (Urey-Bradley, inversion and out-of-plane bends, improper torsions, and the cross terms). Each such
// term gets one of two exact homes, replacing the former post-selection Boltzmann factor:
//
//  - BASE COUPLING: terms fully determined by the step's own sampled coordinates. Every distance the
//    energy uses must lie within {current} + nextBeads (those pairwise distances are set by this
//    step); the previous bead may participate only in purely direction-based terms centered on the
//    current bead (the direction to the previous bead is the cone axis; the placed previous-current
//    DISTANCE must not enter, or the base normalization would depend on the placed configuration).
//    Such terms are automatically spin-invariant -- the previous bead lies on the spin axis, so the
//    spin preserves every internal coordinate of {previous, current} + nextBeads. They are imposed on
//    the base conformation by rejection, and their coupling average is part of the base
//    normalization (see logBaseSamplerNormalization).
//
//  - SPIN TERMS: every other unsampled term of the step (dependence on placed geometry beyond the
//    axis direction, or on the spin angle itself). They enter the torsion-spin selection energy,
//    steering the spin choice exactly like torsions and spin-variant bends -- Rosenbluth-weighted
//    identically on growth and retrace, with no normalization consequences.
// ---------------------------------------------------------------------------------------------------
struct UnsampledStepTerms
{
  Potentials::IntraMolecularPotentials baseCoupling{};
  Potentials::IntraMolecularPotentials spinTerms{};
};

static bool hasUnsampledTerms(const Potentials::IntraMolecularPotentials &terms)
{
  return !(terms.ureyBradleys.empty() && terms.inversionBends.empty() && terms.outOfPlaneBends.empty() &&
           terms.improperTorsions.empty() && terms.bondBonds.empty() && terms.bondBends.empty() &&
           terms.bondTorsions.empty() && terms.bendBends.empty() && terms.bendTorsions.empty());
}

static UnsampledStepTerms splitUnsampledStepTerms(const CBMC::GrowStep &step)
{
  UnsampledStepTerms split{};
  const std::size_t currentBead = step.currentBead;
  const std::size_t previousBead = step.previousBead.value();

  auto isGrown = [&](std::size_t id)
  { return std::find(step.nextBeads.begin(), step.nextBeads.end(), id) != step.nextBeads.end(); };
  auto inStep = [&](std::size_t id) { return id == currentBead || isGrown(id); };
  auto inStepOrPrevious = [&](std::size_t id) { return inStep(id) || id == previousBead; };

  // Urey-Bradley: a plain distance between its two identifiers.
  for (const UreyBradleyPotential &term : step.intra.ureyBradleys)
  {
    const bool base = inStep(term.identifiers[0]) && inStep(term.identifiers[1]);
    (base ? split.baseCoupling : split.spinTerms).ureyBradleys.push_back(term);
  }

  // Inversion bend: the B-centered functional forms depend only on directions from the central atom
  // B, so the previous bead is admissible there; the plane-ACD forms mix distances among A, C, D.
  for (const InversionBendPotential &term : step.intra.inversionBends)
  {
    const bool centeredOnCurrent =
        term.identifiers[1] == currentBead &&
        (term.type == InversionBendType::Harmonic || term.type == InversionBendType::HarmonicCosine ||
         term.type == InversionBendType::Planar);
    const bool base = std::ranges::all_of(term.identifiers, inStep) ||
                      (centeredOnCurrent && std::ranges::all_of(term.identifiers, inStepOrPrevious));
    (base ? split.baseCoupling : split.spinTerms).inversionBends.push_back(term);
  }

  // Out-of-plane bends evaluate to zero for every implemented form; route them to the spin terms
  // (a zero contribution either way).
  split.spinTerms.outOfPlaneBends = step.intra.outOfPlaneBends;

  // Improper torsion: a dihedral A-B-C-D depends only on its three bond unit vectors, so the
  // previous bead is admissible as a terminal atom adjacent to the current bead (that unit vector is
  // the cone axis itself).
  for (const TorsionPotential &term : step.intra.improperTorsions)
  {
    const auto &ids = term.identifiers;
    const bool base = std::ranges::all_of(ids, inStep) ||
                      (ids[0] == previousBead && ids[1] == currentBead && inStep(ids[2]) && inStep(ids[3])) ||
                      (ids[3] == previousBead && ids[2] == currentBead && inStep(ids[0]) && inStep(ids[1]));
    (base ? split.baseCoupling : split.spinTerms).improperTorsions.push_back(term);
  }

  // Bond-bond: two distances from the central identifier.
  for (const BondBondPotential &term : step.intra.bondBonds)
  {
    const bool base = std::ranges::all_of(term.identifiers, inStep);
    (base ? split.baseCoupling : split.spinTerms).bondBonds.push_back(term);
  }

  // Bond-bend: distances A-B and C-B plus the angle at B (the fourth identifier is unused).
  for (const BondBendPotential &term : step.intra.bondBends)
  {
    const bool base = inStep(term.identifiers[0]) && inStep(term.identifiers[1]) && inStep(term.identifiers[2]);
    (base ? split.baseCoupling : split.spinTerms).bondBends.push_back(term);
  }

  // Bond-torsion and bend-torsion mix distances and dihedrals; base only when fully within the step.
  for (const BondTorsionPotential &term : step.intra.bondTorsions)
  {
    const bool base = std::ranges::all_of(term.identifiers, inStep);
    (base ? split.baseCoupling : split.spinTerms).bondTorsions.push_back(term);
  }
  for (const BendTorsionPotential &term : step.intra.bendTorsions)
  {
    const bool base = std::ranges::all_of(term.identifiers, inStep);
    (base ? split.baseCoupling : split.spinTerms).bendTorsions.push_back(term);
  }

  // Bend-bend: both angles at the central identifier B, directions only.
  for (const BendBendPotential &term : step.intra.bendBends)
  {
    const bool base = std::ranges::all_of(term.identifiers, inStepOrPrevious) &&
                      (term.identifiers[1] == currentBead || std::ranges::all_of(term.identifiers, inStep));
    (base ? split.baseCoupling : split.spinTerms).bendBends.push_back(term);
  }

  return split;
}

// A torsion transforms rigidly (is spin-invariant) when all four atoms belong to the ring body
// ('currentBead' plus the 'nextBeads'); it is spin-variant when it also involves an outside atom.
bool isSpinVariantRingTorsion(const std::array<std::size_t, 4> &identifiers, std::size_t currentBead,
                              const std::vector<std::size_t> &nextBeads)
{
  bool hasRingAtom = false;
  bool hasOutsideAtom = false;
  for (std::size_t id : identifiers)
  {
    bool inRingBody = (id == currentBead) || (std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end());
    bool moves = std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end();
    if (moves) hasRingAtom = true;
    if (!inRingBody) hasOutsideAtom = true;
  }
  return hasRingAtom && hasOutsideAtom;
}

// Potentials restricted to the spin-selected (junction-crossing) terms of a ring-closure step: the
// spin-variant torsions (the torsion selection filters the spin-variant bends itself, so all bends
// are kept). The internal ring torsions are sampled by the conformational MC and must be excluded
// here to avoid double counting.
Potentials::IntraMolecularPotentials ringSpinPotentials(const Potentials::IntraMolecularPotentials &intra,
                                                        std::size_t currentBead,
                                                        const std::vector<std::size_t> &nextBeads)
{
  Potentials::IntraMolecularPotentials spin{};
  spin.bends = intra.bends;
  spin.torsions.reserve(intra.torsions.size());
  for (const TorsionPotential &torsion : intra.torsions)
  {
    if (isSpinVariantRingTorsion(torsion.identifiers, currentBead, nextBeads)) spin.torsions.push_back(torsion);
  }
  return spin;
}

// Signed volume of the tetrahedron of a chiral center; its sign is the center's parity.
double chiralSignedVolume(const std::array<std::size_t, 4> &ids, const std::vector<Atom> &atoms)
{
  double3 p0 = atoms[ids[0]].position;
  double3 d1 = atoms[ids[1]].position - p0;
  double3 d2 = atoms[ids[2]].position - p0;
  double3 d3 = atoms[ids[3]].position - p0;
  return double3::dot(d1, double3::cross(d2, d3));
}

}  // namespace

// ---------------------------------------------------------------------------------------------------
// Coupled-decoupled torsion (spin) step: shared by all operators and both directions.
// ---------------------------------------------------------------------------------------------------
struct TorsionOrientation
{
  std::vector<Atom> positions;
  double rosenbluthWeight;
};

static TorsionOrientation selectTorsionOrientation(RandomNumber &random, std::size_t numberOfTorsionTrials, double beta,
                                                   const std::vector<Atom> &chainAtoms,
                                                   const std::vector<Atom> &baseOrientation, std::size_t previousBead,
                                                   std::size_t currentBead, const std::vector<std::size_t> &nextBeads,
                                                   double3 lastBondVector,
                                                   const Potentials::IntraMolecularPotentials &intra,
                                                   bool pinFirstToBase)
{
  std::vector<std::pair<std::vector<Atom>, double>> torsion_orientations(numberOfTorsionTrials);
  std::vector<Atom> chain_atoms(chainAtoms.begin(), chainAtoms.end());

  // Bends to placed atoms other than the previous bead are not invariant under the spin about the
  // previous-current axis and must enter the selection energy here (Rosenbluth-weighted identically
  // on growth and retrace).
  std::vector<const BendPotential *> spin_variant_bends{};
  for (const BendPotential &bend : intra.bends)
  {
    if (isSpinVariantBend(bend, previousBead, currentBead, nextBeads)) spin_variant_bends.push_back(&bend);
  }

  for (std::size_t j = 0; j != numberOfTorsionTrials; ++j)
  {
    double random_angle = (pinFirstToBase && j == 0) ? 0.0 : (2.0 * random.uniform() - 1.0) * std::numbers::pi;

    std::vector<Atom> rotated_atoms = baseOrientation;
    for (std::size_t k = 0; k != rotated_atoms.size(); ++k)
    {
      rotated_atoms[k].position =
          chainAtoms[currentBead].position +
          lastBondVector.rotateAroundAxis(baseOrientation[k].position - chainAtoms[currentBead].position, random_angle);
    }

    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      chain_atoms[nextBeads[k]] = rotated_atoms[k];
    }

    double torsion_energy = intra.calculateTorsionEnergies(chain_atoms);
    for (const BendPotential *bend : spin_variant_bends)
    {
      torsion_energy += bend->calculateEnergy(chain_atoms[bend->identifiers[0]].position,
                                              chain_atoms[bend->identifiers[1]].position,
                                              chain_atoms[bend->identifiers[2]].position, std::nullopt);
    }
    // The step's spin-routed unsampled terms (see splitUnsampledStepTerms): couplings to placed
    // geometry or to the spin angle, steering the spin choice through the Rosenbluth selection.
    torsion_energy += intra.computeInternalEnergiesNotSampledDuringGrowth(chain_atoms).potentialEnergy();
    torsion_orientations[j] = {rotated_atoms, torsion_energy};
  }

  std::vector<double> logTorsionBoltzmannFactors{};
  logTorsionBoltzmannFactors.reserve(numberOfTorsionTrials);
  std::transform(torsion_orientations.begin(), torsion_orientations.end(),
                 std::back_inserter(logTorsionBoltzmannFactors),
                 [&](const std::pair<std::vector<Atom>, double> &v) { return -beta * std::get<1>(v); });

  double rosenbluth_weight_torsion =
      std::accumulate(logTorsionBoltzmannFactors.begin(), logTorsionBoltzmannFactors.end(), 0.0,
                      [](const double &acc, const double &logFactor) { return acc + std::exp(logFactor); });

  std::size_t selected_torsion = pinFirstToBase ? 0 : CBMC::selectTrialPosition(random, logTorsionBoltzmannFactors);

  return {torsion_orientations[selected_torsion].first,
          rosenbluth_weight_torsion / static_cast<double>(numberOfTorsionTrials)};
}


// Minimum of a bend potential over [0, pi]: the envelope offset that makes the sibling-bend
// rejection acceptance exp(-beta (u - u_min)) valid. Analytic for the harmonic forms; a fine grid
// otherwise, lowered by a small slack because a grid minimum can only overestimate the true minimum
// (memoized per potential).
static double bendMinimumEnergy(const BendPotential &bend)
{
  if ((bend.type == BendType::Harmonic || bend.type == BendType::CoreShell) && bend.parameters[1] >= 0.0 &&
      bend.parameters[1] <= std::numbers::pi)
  {
    return 0.0;
  }

  using ParameterArray = std::remove_cvref_t<decltype(bend.parameters)>;
  using CacheKey = std::pair<BendType, ParameterArray>;
  thread_local std::map<CacheKey, double> cache{};
  CacheKey key{bend.type, bend.parameters};
  auto it = cache.find(key);
  if (it != cache.end()) return it->second;

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
  minimum -= Constants::bendMinimumSlack;  // the grid can only overestimate the true minimum

  cache[key] = minimum;
  return minimum;
}

// The anchor bend of a next bead within a step (previous - current - next, centered on the current
// bead), the potential that shapes its cone direction in the exact flexible base sampler.
static const BendPotential *findAnchorBend(const Potentials::IntraMolecularPotentials &intra,
                                           std::size_t previousBead, std::size_t currentBead, std::size_t nextBead)
{
  for (const BendPotential &bend : intra.bends)
  {
    if (bend.identifiers[1] != currentBead) continue;
    if ((bend.identifiers[0] == previousBead && bend.identifiers[2] == nextBead) ||
        (bend.identifiers[2] == previousBead && bend.identifiers[0] == nextBead))
    {
      return &bend;
    }
  }
  return nullptr;
}

// The sibling-sibling bends of a step (next_i - current - next_j).
static std::vector<const BendPotential *> collectSiblingBends(const Potentials::IntraMolecularPotentials &intra,
                                                              std::size_t currentBead,
                                                              const std::vector<std::size_t> &nextBeads)
{
  std::vector<const BendPotential *> sibling_bends{};
  for (const BendPotential &bend : intra.bends)
  {
    if (isSiblingBend(bend, currentBead, nextBeads)) sibling_bends.push_back(&bend);
  }
  return sibling_bends;
}

// The total coupling energy imposed on a step's base conformation by rejection: the sibling bends
// plus the base-routed unsampled terms (see splitUnsampledStepTerms), evaluated on real positions.
static double baseCouplingEnergy(const std::vector<const BendPotential *> &siblingBends,
                                 const Potentials::IntraMolecularPotentials &baseCouplingTerms,
                                 const std::span<const Atom> atoms)
{
  double energy = 0.0;
  for (const BendPotential *bend : siblingBends)
  {
    energy += bend->calculateEnergy(atoms[bend->identifiers[0]].position, atoms[bend->identifiers[1]].position,
                                    atoms[bend->identifiers[2]].position, std::nullopt);
  }
  if (hasUnsampledTerms(baseCouplingTerms))
  {
    energy += baseCouplingTerms.computeInternalEnergiesNotSampledDuringGrowth(atoms).potentialEnergy();
  }
  return energy;
}

// ---------------------------------------------------------------------------------------------------
// The frozen per-step-signature constants of the base coupling, estimated once by a fixed-seed Monte
// Carlo that mirrors the sampler's own independent per-bead draws (memoized):
//
//  - referenceEnergy u_ref: the offset of the rejection acceptance min(1, e^{-beta (u - u_ref)}).
//    Sibling-bends-only steps use the rigorous per-bend minima; steps with general coupling terms
//    (whose cross terms may be negative) use the Monte-Carlo minimum with slack. u_ref does NOT need
//    to be a rigorous bound: whenever u < u_ref the acceptance clamps at one and the excess
//    e^{-beta (u - u_ref)} > 1 rides the trial's Rosenbluth weight instead (see
//    flexibleBaseClampWeight), so the sampled-density x weight product is exactly Boltzmann for any
//    u_ref -- the constant only tunes efficiency.
//
//  - logMeanClampedBoltzmann log<a>: the log mean CLAMPED acceptance a = min(1, e^{-beta (u -
//    u_ref)}) over the independent base. The step's base normalization is
//    Z_indep x <a> x e^{-beta u_ref} (which reduces to Z_indep x <e^{-beta u}> when no clamping
//    occurs); the frozen relative error (< 0.5%) enters the reptation acceptance as a constant at
//    the same sub-percent level.
// ---------------------------------------------------------------------------------------------------
struct BaseCouplingConstants
{
  double referenceEnergy;
  double logMeanClampedBoltzmann;
};

static BaseCouplingConstants stepBaseCouplingConstants(double beta, const Component &component,
                                                       const CBMC::GrowStep &step,
                                                       const std::vector<const BendPotential *> &siblingBends,
                                                       const Potentials::IntraMolecularPotentials &baseCouplingTerms)
{
  const std::size_t previousBead = step.previousBead.value();
  const std::size_t currentBead = step.currentBead;
  const bool hasGeneralCoupling = hasUnsampledTerms(baseCouplingTerms);

  if (siblingBends.empty() && !hasGeneralCoupling) return {0.0, 0.0};

  // The memo key: the step's per-bead samplers plus every coupling term, with atom identifiers
  // mapped to their step-local roles so congruent steps share one entry.
  const auto roleOf = [&](std::size_t id) -> std::string
  {
    if (id == currentBead) return "c";
    if (id == previousBead) return "p";
    return std::format(
        "n{}", std::distance(step.nextBeads.begin(), std::find(step.nextBeads.begin(), step.nextBeads.end(), id)));
  };

  std::string key = std::format("beta={:.12e}", beta);
  for (std::size_t nextBead : step.nextBeads)
  {
    const std::optional<BondPotential> bond = step.intra.findBondPotential(currentBead, nextBead);
    key += bond.has_value() ? std::format(";bond{}", static_cast<std::size_t>(bond->type)) : ";bond-none";
    if (bond.has_value())
      for (double p : bond->parameters) key += std::format(",{:.12e}", p);
    const BendPotential *anchor = findAnchorBend(step.intra, previousBead, currentBead, nextBead);
    key += anchor != nullptr ? std::format(";anchor{}", static_cast<std::size_t>(anchor->type)) : ";anchor-none";
    if (anchor != nullptr)
      for (double p : anchor->parameters) key += std::format(",{:.12e}", p);
  }
  const auto appendTerm = [&](std::string_view tag, std::size_t type, const auto &identifiers, const auto &parameters)
  {
    key += std::format(";{}{}", tag, type);
    for (std::size_t id : identifiers) key += ":" + roleOf(id);
    for (double p : parameters) key += std::format(",{:.12e}", p);
  };
  for (const BendPotential *bend : siblingBends)
    appendTerm("sib", static_cast<std::size_t>(bend->type), bend->identifiers, bend->parameters);
  for (const auto &t : baseCouplingTerms.ureyBradleys)
    appendTerm("ub", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : baseCouplingTerms.inversionBends)
    appendTerm("inv", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : baseCouplingTerms.improperTorsions)
    appendTerm("imp", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : baseCouplingTerms.bondBonds)
    appendTerm("bb", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : baseCouplingTerms.bondBends)
    appendTerm("bB", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : baseCouplingTerms.bondTorsions)
    appendTerm("bt", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : baseCouplingTerms.bendBends)
    appendTerm("BB", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : baseCouplingTerms.bendTorsions)
    appendTerm("Bt", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);

  thread_local std::map<std::string, BaseCouplingConstants> cache{};
  if (auto it = cache.find(key); it != cache.end()) return it->second;

  RandomNumber random(Constants::baseCouplingEstimateSeed);  // frozen: one constant per signature
  const double3 axis{0.0, 0.0, 1.0};

  // Evaluation frame: current bead at the origin, previous bead at unit distance along the axis (the
  // classification guarantees no base term depends on the placed previous-current distance).
  std::vector<Atom> atoms(component.atoms.begin(), component.atoms.end());
  atoms[currentBead].position = double3{0.0, 0.0, 0.0};
  atoms[previousBead].position = axis;

  const auto drawIndependentBase = [&]()
  {
    for (std::size_t nextBead : step.nextBeads)
    {
      const std::optional<BondPotential> bond = step.intra.findBondPotential(currentBead, nextBead);
      const double bondLength =
          bond.has_value() ? bond->generateBondLength(random, beta) : Constants::defaultBondLength;
      const BendPotential *anchor = findAnchorBend(step.intra, previousBead, currentBead, nextBead);
      const double3 direction = anchor != nullptr
                                    ? random.randomVectorOnCone(axis, anchor->generateBendAngle(random, beta))
                                    : random.randomVectorOnUnitSphere();
      atoms[nextBead].position = bondLength * direction;
    }
    return baseCouplingEnergy(siblingBends, baseCouplingTerms, atoms);
  };

  // Reference energy: rigorous per-bend minima when only sibling bends couple; otherwise the
  // Monte-Carlo minimum with slack (exactness does not depend on it, see above).
  double referenceEnergy = 0.0;
  if (!hasGeneralCoupling)
  {
    for (const BendPotential *bend : siblingBends) referenceEnergy += bendMinimumEnergy(*bend);
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

  const BaseCouplingConstants result{referenceEnergy, std::log(sum / numberOfSamples)};
  cache[key] = result;
  return result;
}

// The clamp-excess weight of a base conformation, max(1, e^{-beta (u - u_ref)}): one wherever the
// rejection acceptance was unclamped, the acceptance excess where the coupling energy fell below the
// reference. It multiplies the step's trial weights on growth (freshly sampled base) and retrace
// (the old positions ARE the base), restoring exactness of the density x weight product for any
// reference constant. Returns one for steps without base coupling.
static double flexibleBaseClampWeight(double beta, const Component &component, const CBMC::GrowStep &step,
                                      const std::span<const Atom> atoms)
{
  if (step.rigidBody || step.kind == CBMC::GrowStep::Kind::CloseRing || !step.previousBead.has_value()) return 1.0;

  const std::vector<const BendPotential *> siblingBends =
      collectSiblingBends(step.intra, step.currentBead, step.nextBeads);
  const Potentials::IntraMolecularPotentials baseCouplingTerms = splitUnsampledStepTerms(step).baseCoupling;
  if (siblingBends.empty() && !hasUnsampledTerms(baseCouplingTerms)) return 1.0;

  const BaseCouplingConstants constants =
      stepBaseCouplingConstants(beta, component, step, siblingBends, baseCouplingTerms);
  const double couplingEnergy = baseCouplingEnergy(siblingBends, baseCouplingTerms, atoms);
  return std::max(1.0, std::exp(-beta * (couplingEnergy - constants.referenceEnergy)));
}

// A sampled flexible base conformation: the positions of the step's next beads, and the clamp-excess
// weight of the accepted draw (one unless the coupling energy fell below the reference constant).
struct FlexibleBase
{
  std::vector<Atom> nextBeadAtoms;
  double clampWeight;
};

// ---------------------------------------------------------------------------------------------------
// Flexible-bead base conformation, sampled exactly: bond lengths from their one-dimensional
// Boltzmann densities, each direction from its anchor-bend density on the cone (uniform azimuth),
// and the step's internal couplings -- sibling-sibling bends and the base-routed unsampled terms
// (see splitUnsampledStepTerms) -- imposed by rejection sampling. Carries no Rosenbluth weight
// (apart from the rare clamp excess, see flexibleBaseClampWeight), so this distribution must be the
// exact bonded Boltzmann of its terms: the former internal Metropolis MC (finite, reservoir-seeded,
// adaptive step sizes) only approximated it, and its deviations depended on the structure of the
// growth step -- harmless for moves that pair grow and retrace on the same growth plan, but a
// systematic bias for reptation, which pairs the grow weight of one chain end's plan against the
// retrace weight of the other's.
//
// The base density is Boltzmann but its NORMALIZATION is plan-dependent (per-bead bond and
// anchor-bend integrals, the coupling average, and a factor one half per determined chiral center);
// logBaseSamplerNormalization below computes it, and reptation corrects its acceptance by the
// normalization ratio of its two plans. Spin-variant bends (to placed atoms other than the previous
// bead), all torsions, and the spin-routed unsampled terms must NOT shape the base: they are
// Rosenbluth-weighted in the torsion-spin stage, identically on growth and retrace. Declared chiral
// centers that are fully determined by this step are enforced by parity rejection: the exact
// conditional distribution within the declared-parity sector, which carries probability exactly one
// half by the reflection symmetry (through planes containing the previous-current axis) of the base
// density -- the coupling terms preserve it, since distances, angles, and the cosine-even dihedral
// forms are reflection-invariant. (The torsion spin afterwards is a proper rotation and preserves
// that parity.)
// ---------------------------------------------------------------------------------------------------
static FlexibleBase sampleExactFlexibleBase(RandomNumber &random, double beta, const Component &component,
                                            const std::vector<Atom> &moleculeAtoms, const CBMC::GrowStep &step)
{
  const std::size_t previousBead = step.previousBead.value();
  const std::size_t currentBead = step.currentBead;
  const std::vector<std::size_t> &nextBeads = step.nextBeads;
  const Potentials::IntraMolecularPotentials &intra = step.intra;

  std::vector<Atom> chain_atoms(moleculeAtoms.begin(), moleculeAtoms.end());

  double3 last_bond_vector = chain_atoms[previousBead].position - chain_atoms[currentBead].position;
  if (last_bond_vector.length() < Constants::degenerateAxisLength) last_bond_vector = double3{0.0, 0.0, 1.0};
  last_bond_vector = last_bond_vector.normalized();

  const std::size_t numberOfNextBeads = nextBeads.size();

  // Per next bead: its bond to the current bead and its anchor bend (previous-current-next, centered
  // on the current bead).
  std::vector<std::optional<BondPotential>> bonds(numberOfNextBeads);
  std::vector<const BendPotential *> anchor_bends(numberOfNextBeads, nullptr);
  for (std::size_t i = 0; i != numberOfNextBeads; ++i)
  {
    bonds[i] = intra.findBondPotential(currentBead, nextBeads[i]);
    for (const BendPotential &bend : intra.bends)
    {
      if (bend.identifiers[1] != currentBead) continue;
      bool matches = (bend.identifiers[0] == previousBead && bend.identifiers[2] == nextBeads[i]) ||
                     (bend.identifiers[2] == previousBead && bend.identifiers[0] == nextBeads[i]);
      if (matches)
      {
        anchor_bends[i] = &bend;
        break;
      }
    }
  }

  // The coupling imposed by rejection: sibling-sibling bends plus the base-routed unsampled terms,
  // with the frozen reference constant shared with the normalization (see stepBaseCouplingConstants).
  const std::vector<const BendPotential *> sibling_bends = collectSiblingBends(intra, currentBead, nextBeads);
  const Potentials::IntraMolecularPotentials base_coupling_terms = splitUnsampledStepTerms(step).baseCoupling;
  const bool has_coupling = !sibling_bends.empty() || hasUnsampledTerms(base_coupling_terms);
  const BaseCouplingConstants coupling_constants =
      has_coupling ? stepBaseCouplingConstants(beta, component, step, sibling_bends, base_coupling_terms)
                   : BaseCouplingConstants{0.0, 0.0};

  // Declared chiral centers fully determined by this step: centered on the current bead, with every
  // neighbor either the previous bead or grown here.
  std::vector<const ChiralCenter *> chiral_centers{};
  for (const ChiralCenter &center : component.intraMolecularPotentials.chiralCenters)
  {
    if (center.ids[0] != currentBead) continue;
    bool determined = true;
    for (std::size_t k = 1; k != 4; ++k)
    {
      std::size_t id = center.ids[k];
      determined =
          determined && (id == previousBead || std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end());
    }
    if (determined) chiral_centers.push_back(&center);
  }

  for (std::size_t attempt = 0; attempt != Constants::baseSamplerMaximumAttempts; ++attempt)
  {
    for (std::size_t i = 0; i != numberOfNextBeads; ++i)
    {
      double bond_length =
          bonds[i]
              .transform([&](const BondPotential &b) { return b.generateBondLength(random, beta); })
              .value_or(Constants::defaultBondLength);
      double3 direction =
          anchor_bends[i] != nullptr
              ? random.randomVectorOnCone(last_bond_vector, anchor_bends[i]->generateBendAngle(random, beta))
              : random.randomVectorOnUnitSphere();
      chain_atoms[nextBeads[i]].position = chain_atoms[currentBead].position + bond_length * direction;
    }

    double clamp_weight = 1.0;
    if (has_coupling)
    {
      const double coupling_energy = baseCouplingEnergy(sibling_bends, base_coupling_terms, chain_atoms);
      const double boltzmann_excess = std::exp(-beta * (coupling_energy - coupling_constants.referenceEnergy));
      if (random.uniform() > boltzmann_excess) continue;
      clamp_weight = std::max(1.0, boltzmann_excess);
    }

    bool parity_ok = true;
    for (const ChiralCenter *center : chiral_centers)
    {
      double reference_sign = center->type == ChiralCenter::Chirality::R_Chiral ? 1.0 : -1.0;
      parity_ok = parity_ok && (reference_sign * chiralSignedVolume(center->ids, chain_atoms) > 0.0);
    }
    if (!parity_ok) continue;

    std::vector<Atom> next_bead_atoms(numberOfNextBeads);
    for (std::size_t i = 0; i != numberOfNextBeads; ++i) next_bead_atoms[i] = chain_atoms[nextBeads[i]];
    return {next_bead_atoms, clamp_weight};
  }
  throw std::runtime_error(
      "CBMC: the exact base-conformation sampler exceeded its rejection budget (pathologically stiff "
      "base couplings or a nearly unsatisfiable chiral constraint)\n");
}

// ---------------------------------------------------------------------------------------------------
// The base-sampler normalization of a plan: for each flexible attach step, the bead densities
// r^2 exp(-beta u_bond) and sin(theta) exp(-beta u_anchor) (uniform azimuth) integrate to the
// product of the potentials' one-dimensional normalizations, the coupling rejection (sibling bends
// plus base-routed unsampled terms) multiplies in <a> e^{-beta u_ref} (see
// stepBaseCouplingConstants), and each fully determined chiral center restricts the base to a
// sector of probability exactly one half. See the interface documentation for why reptation needs
// this.
// ---------------------------------------------------------------------------------------------------
bool CBMC::stepHandlesUnsampledInternalTerms(const GrowStep &step)
{
  return step.kind == GrowStep::Kind::AttachFragment && !step.rigidBody && step.previousBead.has_value();
}

double CBMC::logBaseSamplerNormalization(double beta, const Component &component,
                                         const std::vector<CBMC::GrowStep> &plan)
{
  double logNormalization = 0.0;

  for (const CBMC::GrowStep &step : plan)
  {
    // Rigid-body and ring-closure steps keep internal-MC samplers without a closed-form
    // normalization; reptation enforces congruent end plans for units containing them, so their
    // (equal) factors cancel and are skipped here. Seed steps do not occur in partial plans.
    if (step.rigidBody || step.kind == CBMC::GrowStep::Kind::CloseRing || !step.previousBead.has_value()) continue;

    const std::size_t previousBead = step.previousBead.value();
    const std::size_t currentBead = step.currentBead;

    for (std::size_t nextBead : step.nextBeads)
    {
      const std::optional<BondPotential> bond = step.intra.findBondPotential(currentBead, nextBead);
      // A bond without potential is placed at the fixed default length: r^2 dr integrates to r^2.
      logNormalization += bond.has_value() ? bond->logBoltzmannVolumeNormalization(beta)
                                           : 2.0 * std::log(Constants::defaultBondLength);

      const BendPotential *anchorBend = findAnchorBend(step.intra, previousBead, currentBead, nextBead);
      logNormalization += anchorBend != nullptr ? anchorBend->logBoltzmannConeNormalization(beta)
                                                : std::log(4.0 * std::numbers::pi);  // uniform sphere
    }

    const std::vector<const BendPotential *> siblingBends =
        collectSiblingBends(step.intra, currentBead, step.nextBeads);
    const Potentials::IntraMolecularPotentials baseCouplingTerms = splitUnsampledStepTerms(step).baseCoupling;
    if (!siblingBends.empty() || hasUnsampledTerms(baseCouplingTerms))
    {
      const BaseCouplingConstants constants =
          stepBaseCouplingConstants(beta, component, step, siblingBends, baseCouplingTerms);
      logNormalization += constants.logMeanClampedBoltzmann - beta * constants.referenceEnergy;
    }

    for (const ChiralCenter &center : component.intraMolecularPotentials.chiralCenters)
    {
      if (center.ids[0] != currentBead) continue;
      bool determined = true;
      for (std::size_t k = 1; k != 4; ++k)
      {
        std::size_t id = center.ids[k];
        determined = determined && (id == previousBead ||
                                    std::find(step.nextBeads.begin(), step.nextBeads.end(), id) != step.nextBeads.end());
      }
      if (determined) logNormalization += std::log(0.5);
    }
  }

  return logNormalization;
}

// ---------------------------------------------------------------------------------------------------
// Rigid-body tilt: samples the junction-bend tilt of a rigid fragment hinged on the anchor with a
// rigid-rotation Metropolis MC. Carries no Rosenbluth weight. Ported from the former
// 'generateRigidUnitOrientationMonteCarloScheme'.
// ---------------------------------------------------------------------------------------------------
static std::vector<Atom> generateRigidTilt(RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead, double beta,
                                           const Component &component, const std::vector<Atom> &chainAtoms,
                                           std::optional<std::size_t> previousBead, std::size_t currentBead,
                                           const std::vector<std::size_t> &nextBeads,
                                           const Potentials::IntraMolecularPotentials &intra)
{
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

  std::size_t inner = nextBeads[0];
  for (std::size_t atom : nextBeads)
  {
    if (component.connectivityTable[atom, currentBead])
    {
      inner = atom;
      break;
    }
  }

  double bend_angle = Constants::defaultRigidJunctionBendAngle;
  for (const BendPotential &bend : intra.bends)
  {
    if (bend.identifiers[1] != currentBead) continue;
    if ((bend.identifiers[0] == previousBead.value() && bend.identifiers[2] == inner) ||
        (bend.identifiers[0] == inner && bend.identifiers[2] == previousBead.value()))
    {
      bend_angle = bend.generateBendAngle(random, beta);
      break;
    }
  }

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

// ---------------------------------------------------------------------------------------------------
// Ring-closure conformation: samples the internal conformation of a cyclic cluster (kept closed by
// its closure bonds) plus the junction tilt. Carries no Rosenbluth weight. Ported from the former
// 'generateRingConformationMonteCarloScheme'. A cyclic cluster may pass through rigid-body
// fragments (e.g. a macrocycle of rigid rings linked by flexible bridges); those fragments are
// moved as whole rigid units so their internal geometry stays exact.
// ---------------------------------------------------------------------------------------------------
static bool bendIsSpinVariantRing(const std::array<std::size_t, 3> &identifiers, std::optional<std::size_t> previousBead,
                                  std::size_t currentBead, const std::vector<std::size_t> &nextBeads)
{
  for (std::size_t id : identifiers)
  {
    if (id == currentBead) continue;
    if (previousBead.has_value() && id == previousBead.value()) continue;
    if (std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end()) continue;
    return true;
  }
  return false;
}

static std::vector<Atom> randomlyOrientRing(RandomNumber &random, const std::vector<Atom> &ringAtoms,
                                            double3 anchorPosition)
{
  double3x3 rotation = random.randomRotationMatrix();
  std::vector<Atom> result = ringAtoms;
  for (Atom &atom : result) atom.position = anchorPosition + rotation * (atom.position - anchorPosition);
  return result;
}

static std::vector<Atom> generateRingConformation(RandomNumber &random, const ForceField &forceField, double beta,
                                                  const Component &component, const std::vector<Atom> &chainAtoms,
                                                  std::optional<std::size_t> previousBead, std::size_t currentBead,
                                                  const std::vector<std::size_t> &nextBeads,
                                                  const Potentials::IntraMolecularPotentials &intra)
{
  // Seed geometry: an independent, well-mixed ideal-gas conformation from the reservoir (so different
  // grows start from different ring puckers), or the reference geometry while the reservoir is still
  // being built. A reservoir member keeps the exact internal geometry of any rigid sub-fragment, so
  // rigid parts of the ring stay rigid. The declared chirality reference below stays 'component.atoms'.
  const std::vector<Atom> &seed =
      component.conformationReservoir.empty()
          ? component.atoms
          : component.conformationReservoir[random.uniform_integer(0, component.conformationReservoir.size() - 1)];
  double3 anchor_reference = seed[currentBead].position;
  double3 anchor_position = chainAtoms[currentBead].position;
  std::vector<Atom> chain_atoms(chainAtoms.begin(), chainAtoms.end());

  auto placeWithRotation = [&](const double3x3 &rotation)
  {
    for (std::size_t atom : nextBeads)
    {
      double3 offset = seed[atom].position - anchor_reference;
      chain_atoms[atom].position = anchor_position + rotation * offset;
    }
  };

  // Chiral centers whose four atoms all have known positions during this step (the ring body, the
  // anchor, and the junction's placed neighbor) keep the parity of the reference geometry: the
  // bond/bend/torsion model is achiral (a mirror image has the same energy), so without this guard
  // the internal MC could invert a declared stereocenter, e.g. flip a cis ring fusion to trans.
  std::vector<std::pair<const std::array<std::size_t, 4> *, double>> monitored_centers{};
  {
    auto isKnown = [&](std::size_t id)
    {
      if (id == currentBead) return true;
      if (previousBead.has_value() && id == previousBead.value()) return true;
      return std::find(nextBeads.begin(), nextBeads.end(), id) != nextBeads.end();
    };
    for (const ChiralCenter &center : component.intraMolecularPotentials.chiralCenters)
    {
      if (std::all_of(center.ids.begin(), center.ids.end(), isKnown))
      {
        monitored_centers.push_back({&center.ids, chiralSignedVolume(center.ids, component.atoms)});
      }
    }
  }
  auto parityPreserved = [&]()
  {
    for (const auto &[ids, referenceVolume] : monitored_centers)
    {
      if (chiralSignedVolume(*ids, chain_atoms) * referenceVolume < 0.0) return false;
    }
    return true;
  };

  // A proper rotation of the reference geometry preserves the parity of every chiral center that
  // lies entirely inside the rigidly placed ring body (its signed volume keeps the reference sign),
  // so only a center that also involves the junction's placed neighbour can come out with the wrong
  // parity. Reflecting the whole placed body through a plane that contains the previous-current axis
  // flips exactly those junction-involving centers while preserving all internal distances and the
  // bend angles to the anchor -- the bonded model is achiral, so the reflected seed has identical
  // energy. Because a reflection is improper it also flips any body-internal center, so it is only a
  // valid fix when it restores every monitored parity at once; we therefore try it, keep it only when
  // 'parityPreserved()' then holds, and otherwise fall back to re-rolling random orientations. This
  // replaces an unconditional (up to 1000) random re-roll that could also fail outright, and it makes
  // the dominant single-junction-stereocentre case an O(1), deterministic correction.
  placeWithRotation(random.randomRotationMatrix());
  if (!parityPreserved() && previousBead.has_value())
  {
    double3 axis = (chain_atoms[previousBead.value()].position - anchor_position).normalized();
    double3 normal = double3::perpendicular(axis, random.randomVectorOnUnitSphere());
    std::vector<double3> beforeReflection(nextBeads.size());
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      beforeReflection[k] = chain_atoms[nextBeads[k]].position;
      double3 relative = beforeReflection[k] - anchor_position;
      chain_atoms[nextBeads[k]].position =
          anchor_position + (relative - 2.0 * double3::dot(relative, normal) * normal);
    }
    if (!parityPreserved())
    {
      for (std::size_t k = 0; k != nextBeads.size(); ++k) chain_atoms[nextBeads[k]].position = beforeReflection[k];
    }
  }
  for (std::size_t attempt = 0; attempt != Constants::ringParityMaximumReorientations && !parityPreserved();
       ++attempt)
  {
    placeWithRotation(random.randomRotationMatrix());
  }

  // Move units of the conformational MC: a single-atom (flexible) fragment moves by per-atom
  // displacement, a rigid-body fragment moves as one unit (translation or rotation about its
  // center) so its internal geometry is preserved exactly. The growth plan places a rigid fragment
  // either entirely inside this step or entirely before it, never partially.
  const FragmentGraph &graph = component.fragmentGraph;
  std::vector<std::vector<std::size_t>> moveUnits{};
  {
    std::map<std::size_t, std::size_t> rigidFragmentUnits{};
    for (std::size_t atom : nextBeads)
    {
      std::size_t fragmentIndex = graph.atomFragmentIds[atom];
      if (!graph.fragments[fragmentIndex].isRigidBody())
      {
        moveUnits.push_back({atom});
        continue;
      }
      auto [it, inserted] = rigidFragmentUnits.insert({fragmentIndex, moveUnits.size()});
      if (inserted) moveUnits.push_back({});
      moveUnits[it->second].push_back(atom);
    }
  }

  // Fixed-bond neighbours of each atom. A Fixed bond is a holonomic distance constraint that carries
  // no energy, so a free Cartesian displacement of a ring atom would stretch it with no penalty and be
  // accepted. Any move of a flexible ring atom is therefore built as a rotation about its fixed
  // neighbour(s), which preserves those bond lengths exactly (a rotation about an axis through a point
  // preserves the distance to it). 'intra.bonds' only involves atoms placed in this step, all of which
  // have positions, so every listed endpoint is usable as a pivot.
  std::vector<std::vector<std::size_t>> fixedNeighbors(chain_atoms.size());
  for (const BondPotential &bond : intra.bonds)
  {
    if (bond.type != BondType::Fixed) continue;
    fixedNeighbors[bond.identifiers[0]].push_back(bond.identifiers[1]);
    fixedNeighbors[bond.identifiers[1]].push_back(bond.identifiers[0]);
  }

  // Conformer-hopping (crankshaft) candidates: a flexible (single-atom fragment) ring atom rotated by
  // a large angle about the line through two of its positioned bonded neighbours. This preserves both
  // of those bond lengths exactly while flipping the local pucker (chair <-> twist-boat) and
  // axial <-> equatorial placement -- the barrier crossing that the small adaptive moves almost never
  // make on their own. The axis is chosen to include every Fixed-bond neighbour of the atom (fixed
  // neighbours first), so the crankshaft can never break a Fixed bond; an atom with three or more
  // Fixed bonds is over-constrained (its position is pinned, e.g. a fused-ring junction with fixed
  // bond lengths) and is left to a concerted move, not offered here. The proposal is symmetric
  // (uniform +/- angle), so plain Metropolis on the base energy keeps the sampled distribution exact.
  struct CrankshaftCandidate
  {
    std::size_t atom;
    std::size_t axisA;
    std::size_t axisB;
  };
  std::vector<CrankshaftCandidate> crankshafts{};
  {
    std::vector<bool> positioned(chain_atoms.size(), false);
    positioned[currentBead] = true;
    if (previousBead.has_value()) positioned[previousBead.value()] = true;
    for (std::size_t atom : nextBeads) positioned[atom] = true;
    for (std::size_t atom : nextBeads)
    {
      if (graph.fragments[graph.atomFragmentIds[atom]].isRigidBody()) continue;
      if (fixedNeighbors[atom].size() >= 3) continue;

      // Fixed neighbours first (they must lie on the axis), then any other positioned bonded
      // neighbours; the first two form the crankshaft axis.
      std::vector<std::size_t> axisNeighbors = fixedNeighbors[atom];
      for (std::size_t other = 0; other != chain_atoms.size(); ++other)
      {
        if (other == atom || !positioned[other]) continue;
        if (!component.connectivityTable[atom, other]) continue;
        if (std::find(fixedNeighbors[atom].begin(), fixedNeighbors[atom].end(), other) != fixedNeighbors[atom].end())
          continue;
        axisNeighbors.push_back(other);
      }
      if (axisNeighbors.size() >= 2) crankshafts.push_back({atom, axisNeighbors[0], axisNeighbors[1]});
    }
  }

  Potentials::IntraMolecularPotentials baseIntra{};
  baseIntra.bonds = intra.bonds;
  for (const BendPotential &bend : intra.bends)
  {
    if (!bendIsSpinVariantRing(bend.identifiers, previousBead, currentBead, nextBeads)) baseIntra.bends.push_back(bend);
  }
  for (const TorsionPotential &torsion : intra.torsions)
  {
    if (!isSpinVariantRingTorsion(torsion.identifiers, currentBead, nextBeads)) baseIntra.torsions.push_back(torsion);
  }

  auto baseEnergy = [&](std::vector<Atom> &atoms)
  {
    return baseIntra.calculateBondSmallMCEnergies(atoms) + baseIntra.calculateBendSmallMCEnergies(atoms) +
           baseIntra.calculateTorsionEnergies(atoms);
  };

  double current_energy = baseEnergy(chain_atoms);

  // Adaptive internal-MC step sizes: the maximum displacement and rotation angle are read from the
  // anchor bead's CBMC statistics and adapted towards the target acceptance ratio between sweeps by
  // 'System::optimizeMCMoves'. These moves carry no Rosenbluth weight, so their step size affects
  // only sampling efficiency, not detailed balance.
  MoveStatistics<double> &displacementStats = component.cbmc_moves_statistics[currentBead].ringDisplacementChange;
  MoveStatistics<double> &rotationStats = component.cbmc_moves_statistics[currentBead].ringRotationChange;
  MoveStatistics<double> &crankshaftStats = component.cbmc_moves_statistics[currentBead].ringCrankshaftMove;
  const double maximumDisplacement = displacementStats.maxChange;
  const double maximumRotationAngle = rotationStats.maxChange;
  const bool haveJunction = previousBead.has_value();
  std::size_t number_of_trials = 2 * forceField.numberOfTrialMovesPerOpenBead * nextBeads.size();
  std::vector<double3> saved(nextBeads.size());

  auto acceptOrReject = [&](MoveStatistics<double> &stats, std::size_t unitSize, auto restore)
  {
    stats.counts += 1.0;
    stats.totalCounts += 1.0;
    stats.constructed += 1.0;
    stats.totalConstructed += 1.0;
    double trial_energy = baseEnergy(chain_atoms);
    if (parityPreserved() && random.uniform() < std::exp(-beta * (trial_energy - current_energy)))
    {
      current_energy = trial_energy;
      stats.accepted += 1.0;
      stats.totalAccepted += 1.0;
    }
    else
    {
      restore(unitSize);
    }
  };

  for (std::size_t trial = 0; trial != number_of_trials; ++trial)
  {
    // Conformer-hopping crankshaft: a large-angle rotation of one flexible ring atom about the line
    // through two of its neighbours (see 'crankshafts' above). Attempted a fraction of the time
    // ('CBMCRingCrankshaftProbability') so local relaxation still dominates; it supplies the barrier
    // crossings between ring conformers. Tracked in its own statistics: the angle is deliberately
    // full-range and never adapted, and pooling its acceptances into the adaptive rotation statistics
    // would distort that step-size optimization.
    if (!crankshafts.empty() && random.uniform() < forceField.cbmcRingCrankshaftProbability)
    {
      const CrankshaftCandidate &c = crankshafts[random.uniform_integer(0, crankshafts.size() - 1)];
      double3 pivot = chain_atoms[c.axisA].position;
      double3 axis = (chain_atoms[c.axisB].position - pivot).normalized();
      double angle = (2.0 * random.uniform() - 1.0) * std::numbers::pi;
      double3 savedPosition = chain_atoms[c.atom].position;
      chain_atoms[c.atom].position = pivot + axis.rotateAroundAxis(savedPosition - pivot, angle);
      acceptOrReject(crankshaftStats, 1, [&](std::size_t) { chain_atoms[c.atom].position = savedPosition; });
      continue;
    }

    if (haveJunction && random.uniform() < forceField.cbmcRingTiltProbability)
    {
      double3 axis = random.randomVectorOnUnitSphere();
      double angle = (2.0 * random.uniform() - 1.0) * maximumRotationAngle;
      for (std::size_t k = 0; k != nextBeads.size(); ++k)
      {
        saved[k] = chain_atoms[nextBeads[k]].position;
        chain_atoms[nextBeads[k]].position = anchor_position + axis.rotateAroundAxis(saved[k] - anchor_position, angle);
      }
      acceptOrReject(rotationStats, nextBeads.size(),
                     [&](std::size_t n)
                     {
                       for (std::size_t k = 0; k != n; ++k) chain_atoms[nextBeads[k]].position = saved[k];
                     });
    }
    else
    {
      const std::vector<std::size_t> &unit = moveUnits[random.uniform_integer(0, moveUnits.size() - 1)];
      for (std::size_t k = 0; k != unit.size(); ++k) saved[k] = chain_atoms[unit[k]].position;
      auto restore = [&](std::size_t n)
      { for (std::size_t k = 0; k != n; ++k) chain_atoms[unit[k]].position = saved[k]; };

      if (unit.size() == 1)
      {
        // Local single-atom move, constraint-preserving. With no Fixed bond the atom is displaced
        // freely; with Fixed bonds it is rotated about its fixed neighbour(s) so those exact lengths
        // are kept: one fixed neighbour leaves a full sphere (rotate about a random axis through it),
        // two leave a circle (rotate about the line through both), and three or more pin the atom, so
        // it moves only through the whole-cluster rotation or a concerted move.
        const std::vector<std::size_t> &fixed = fixedNeighbors[unit[0]];
        if (fixed.empty())
        {
          double3 displacement{(2.0 * random.uniform() - 1.0) * maximumDisplacement,
                               (2.0 * random.uniform() - 1.0) * maximumDisplacement,
                               (2.0 * random.uniform() - 1.0) * maximumDisplacement};
          chain_atoms[unit[0]].position += displacement;
          acceptOrReject(displacementStats, unit.size(), restore);
        }
        else if (fixed.size() <= 2)
        {
          double3 pivot = chain_atoms[fixed[0]].position;
          double3 axis = fixed.size() == 1 ? random.randomVectorOnUnitSphere()
                                           : (chain_atoms[fixed[1]].position - pivot).normalized();
          double angle = (2.0 * random.uniform() - 1.0) * maximumRotationAngle;
          chain_atoms[unit[0]].position = pivot + axis.rotateAroundAxis(chain_atoms[unit[0]].position - pivot, angle);
          acceptOrReject(rotationStats, unit.size(), restore);
        }
      }
      else if (random.uniform() < 0.5)
      {
        // Rigid fragment: symmetric whole-unit translation.
        double3 displacement{(2.0 * random.uniform() - 1.0) * maximumDisplacement,
                             (2.0 * random.uniform() - 1.0) * maximumDisplacement,
                             (2.0 * random.uniform() - 1.0) * maximumDisplacement};
        for (std::size_t atom : unit) chain_atoms[atom].position += displacement;
        acceptOrReject(displacementStats, unit.size(), restore);
      }
      else
      {
        // Rigid fragment: symmetric rotation about the unit center.
        double3 center{};
        for (std::size_t atom : unit) center += chain_atoms[atom].position;
        center = center / static_cast<double>(unit.size());
        double3 axis = random.randomVectorOnUnitSphere();
        double angle = (2.0 * random.uniform() - 1.0) * maximumRotationAngle;
        for (std::size_t atom : unit)
        {
          chain_atoms[atom].position = center + axis.rotateAroundAxis(chain_atoms[atom].position - center, angle);
        }
        acceptOrReject(rotationStats, unit.size(), restore);
      }
    }
  }

  std::vector<Atom> result(nextBeads.size());
  for (std::size_t k = 0; k != nextBeads.size(); ++k) result[k] = chain_atoms[nextBeads[k]];
  return result;
}

// ---------------------------------------------------------------------------------------------------
// Base-conformation dispatch shared by grow / retrace / recoil.
// ---------------------------------------------------------------------------------------------------
static FlexibleBase sampleBaseConformation(RandomNumber &random, const ForceField &forceField, double beta,
                                           const Component &component, const std::vector<Atom> &chainAtoms,
                                           const CBMC::GrowStep &step)
{
  if (step.kind == CBMC::GrowStep::Kind::CloseRing)
  {
    return {generateRingConformation(random, forceField, beta, component, chainAtoms, step.previousBead,
                                     step.currentBead, step.nextBeads, step.intra),
            1.0};
  }
  if (step.rigidBody)
  {
    return {generateRigidTilt(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, chainAtoms,
                              step.previousBead, step.currentBead, step.nextBeads, step.intra),
            1.0};
  }
  return sampleExactFlexibleBase(random, beta, component, chainAtoms, step);
}

// The intramolecular potentials used in the torsion (spin) selection: the junction-crossing subset
// for a ring-closure step, the full step potentials otherwise.
// The potentials evaluated by the torsion-spin selection of a step: all torsions and bends (the
// selection filters the spin-variant ones itself), plus -- for flexible attach steps -- the
// spin-routed share of the unsampled terms (see splitUnsampledStepTerms). Ring and rigid-body steps
// carry no unsampled terms here: theirs stay in the caller's post-selection factor.
static Potentials::IntraMolecularPotentials torsionSelectionPotentials(const CBMC::GrowStep &step)
{
  if (step.kind == CBMC::GrowStep::Kind::CloseRing)
  {
    return ringSpinPotentials(step.intra, step.currentBead, step.nextBeads);
  }

  Potentials::IntraMolecularPotentials spin{};
  spin.torsions = step.intra.torsions;
  spin.bends = step.intra.bends;
  if (!step.rigidBody && step.previousBead.has_value())
  {
    UnsampledStepTerms split = splitUnsampledStepTerms(step);
    spin.ureyBradleys = std::move(split.spinTerms.ureyBradleys);
    spin.inversionBends = std::move(split.spinTerms.inversionBends);
    spin.outOfPlaneBends = std::move(split.spinTerms.outOfPlaneBends);
    spin.improperTorsions = std::move(split.spinTerms.improperTorsions);
    spin.bondBonds = std::move(split.spinTerms.bondBonds);
    spin.bondBends = std::move(split.spinTerms.bondBends);
    spin.bondTorsions = std::move(split.spinTerms.bondTorsions);
    spin.bendBends = std::move(split.spinTerms.bendBends);
    spin.bendTorsions = std::move(split.spinTerms.bendTorsions);
  }
  return spin;
}

// ---------------------------------------------------------------------------------------------------
// Public operators.
// ---------------------------------------------------------------------------------------------------
std::vector<CBMC::StepTrial> CBMC::generateGrowTrials(RandomNumber &random, const ForceField &forceField, double beta,
                                                      const Component &component, const std::vector<Atom> &chainAtoms,
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
                                       std::nullopt, step.currentBead, step.nextBeads, step.intra),
                     1.0};
      }
      return trials;
    }
    if (step.kind == GrowStep::Kind::CloseRing)
    {
      // Ring seed: sample one internal conformation and rigidly rotate it for the other directions.
      std::vector<Atom> base =
          generateRingConformation(random, forceField, beta, component, chainAtoms,
                                   std::nullopt, step.currentBead, step.nextBeads, step.intra);
      trials[0] = {base, 1.0};
      for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {randomlyOrientRing(random, base, chainAtoms[step.currentBead].position), 1.0};
      }
      return trials;
    }

    // Flexible single-bond seed: Boltzmann bond length in a uniformly random direction.
    const BondPotential *bond = step.intra.bonds.empty() ? nullptr : &step.intra.bonds.front();
    for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
    {
      double bond_length = bond ? bond->generateBondLength(random, beta) : Constants::defaultBondLength;
      double3 unit_vector = random.randomVectorOnUnitSphere();
      Atom trial_atom = chainAtoms[step.nextBeads[0]];
      trial_atom.position = chainAtoms[step.currentBead].position + bond_length * unit_vector;
      trials[i] = {{trial_atom}, 1.0};
    }
    return trials;
  }

  // Attach / ring-closure with a junction: one shared base conformation, one torsion spin per
  // direction about the junction bond. The base's clamp-excess weight is shared by every direction.
  FlexibleBase base = sampleBaseConformation(random, forceField, beta, component, chainAtoms, step);
  Potentials::IntraMolecularPotentials torsionIntra = torsionSelectionPotentials(step);
  double3 last_bond_vector =
      (chainAtoms[step.previousBead.value()].position - chainAtoms[step.currentBead].position).normalized();

  for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
  {
    TorsionOrientation torsion =
        selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta, chainAtoms,
                                 base.nextBeadAtoms, step.previousBead.value(), step.currentBead, step.nextBeads,
                                 last_bond_vector, torsionIntra, false);
    trials[i] = {torsion.positions, torsion.rosenbluthWeight * base.clampWeight};
  }
  return trials;
}

std::vector<CBMC::StepTrial> CBMC::generateRetraceTrials(RandomNumber &random, const ForceField &forceField,
                                                         double beta, const Component &component,
                                                         const std::vector<Atom> &chainAtoms, const GrowStep &step,
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
      trials[0] = {old_orientation, 1.0};
      for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {generateRigidTilt(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, chainAtoms,
                                       std::nullopt, step.currentBead, step.nextBeads, step.intra),
                     1.0};
      }
      return trials;
    }
    if (step.kind == GrowStep::Kind::CloseRing)
    {
      trials[0] = {old_orientation, 1.0};
      for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        trials[i] = {randomlyOrientRing(random, old_orientation, chainAtoms[step.currentBead].position), 1.0};
      }
      return trials;
    }

    const BondPotential *bond = step.intra.bonds.empty() ? nullptr : &step.intra.bonds.front();
    trials[0] = {old_orientation, 1.0};
    for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
    {
      double bond_length = bond ? bond->generateBondLength(random, beta) : Constants::defaultBondLength;
      double3 unit_vector = random.randomVectorOnUnitSphere();
      Atom trial_atom = chainAtoms[step.nextBeads[0]];
      trial_atom.position = chainAtoms[step.currentBead].position + bond_length * unit_vector;
      trials[i] = {{trial_atom}, 1.0};
    }
    return trials;
  }

  // Attach / ring-closure with a junction. The old orientation is the shared torsion base for every
  // trial direction, pinned as torsion trial 0 of the first trial direction; its clamp-excess weight
  // (the old positions ARE the base) is shared by every direction, mirroring the grow side.
  Potentials::IntraMolecularPotentials torsionIntra = torsionSelectionPotentials(step);
  const double base_clamp_weight = flexibleBaseClampWeight(beta, component, step, chainAtoms);
  double3 last_bond_vector =
      (chainAtoms[step.previousBead.value()].position - chainAtoms[step.currentBead].position).normalized();

  for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
  {
    TorsionOrientation torsion =
        selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta, chainAtoms, old_orientation,
                                 step.previousBead.value(), step.currentBead, step.nextBeads, last_bond_vector,
                                 torsionIntra, i == 0);
    trials[i] = {i == 0 ? old_orientation : torsion.positions, torsion.rosenbluthWeight * base_clamp_weight};
  }
  return trials;
}

CBMC::StepTrial CBMC::generateRecoilTrial(RandomNumber &random, const ForceField &forceField, double beta,
                                          const Component &component, const std::vector<Atom> &contextAtoms,
                                          const GrowStep &step)
{
  if (!step.previousBead.has_value())
  {
    if (step.rigidBody)
    {
      return {generateRigidTilt(random, forceField.numberOfTrialMovesPerOpenBead, beta, component, contextAtoms,
                                std::nullopt, step.currentBead, step.nextBeads, step.intra),
              1.0};
    }
    if (step.kind == GrowStep::Kind::CloseRing)
    {
      return {generateRingConformation(random, forceField, beta, component, contextAtoms,
                                       std::nullopt, step.currentBead, step.nextBeads, step.intra),
              1.0};
    }
    const BondPotential *bond = step.intra.bonds.empty() ? nullptr : &step.intra.bonds.front();
    double bond_length = bond ? bond->generateBondLength(random, beta) : Constants::defaultBondLength;
    double3 unit_vector = random.randomVectorOnUnitSphere();
    Atom trial_atom = contextAtoms[step.nextBeads[0]];
    trial_atom.position = contextAtoms[step.currentBead].position + bond_length * unit_vector;
    return {{trial_atom}, 1.0};
  }

  // Every flexible step (single bead or branch) draws its base from the exact bonded sampler; ring
  // closure and rigid fragments keep their internal Monte Carlo. The spin about the junction bond is
  // then always torsion-selected -- also when the trial is a feeler bead (see the interface note).
  FlexibleBase base = sampleBaseConformation(random, forceField, beta, component, contextAtoms, step);
  double3 last_bond_vector = contextAtoms[step.previousBead.value()].position - contextAtoms[step.currentBead].position;
  if (last_bond_vector.length() < Constants::degenerateAxisLength) last_bond_vector = double3{0.0, 0.0, 1.0};
  last_bond_vector = last_bond_vector.normalized();

  Potentials::IntraMolecularPotentials torsionIntra = torsionSelectionPotentials(step);
  TorsionOrientation torsion =
      selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta, contextAtoms,
                               base.nextBeadAtoms, step.previousBead.value(), step.currentBead, step.nextBeads,
                               last_bond_vector, torsionIntra, false);
  return {torsion.positions, torsion.rosenbluthWeight * base.clampWeight};
}

double CBMC::oldConfigurationTorsionWeight(RandomNumber &random, const ForceField &forceField, double beta,
                                           const Component &component, const std::vector<Atom> &oldAtoms,
                                           const GrowStep &step)
{
  if (!step.previousBead.has_value()) return 1.0;

  std::vector<Atom> old_orientation(step.nextBeads.size());
  for (std::size_t k = 0; k != step.nextBeads.size(); ++k) old_orientation[k] = oldAtoms[step.nextBeads[k]];

  Potentials::IntraMolecularPotentials torsionIntra = torsionSelectionPotentials(step);
  const double base_clamp_weight = flexibleBaseClampWeight(beta, component, step, oldAtoms);
  double3 last_bond_vector =
      (oldAtoms[step.previousBead.value()].position - oldAtoms[step.currentBead].position).normalized();

  TorsionOrientation torsion =
      selectTorsionOrientation(random, forceField.numberOfTorsionTrialDirections, beta, oldAtoms, old_orientation,
                               step.previousBead.value(), step.currentBead, step.nextBeads, last_bond_vector,
                               torsionIntra, true);
  return torsion.rosenbluthWeight * base_clamp_weight;
}
