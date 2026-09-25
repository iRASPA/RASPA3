module;

module cbmc_step_terms;

import std;

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
import connectivity_table;
import fragment;
import fragment_graph;
import cbmc_growth_plan;

namespace
{
// ---------------------------------------------------------------------------------------------------
// Topological classification of a step's internal terms. Everything here depends only on the step
// (which beads are previous / current / grown) and is evaluated once per step when the plan is built.
// ---------------------------------------------------------------------------------------------------

bool contains(const std::vector<std::size_t> &beads, std::size_t id)
{
  return std::find(beads.begin(), beads.end(), id) != beads.end();
}

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
    if (contains(nextBeads, id)) continue;
    return true;
  }
  return false;
}

// A sibling bend couples two beads grown in the same step (next_i - current - next_j). It is
// spin-invariant (the whole branch rotates rigidly) and is imposed on the base conformation by
// rejection sampling; its coupling integral is part of the base normalization (see
// CBMC::logBaseSamplerNormalization).
bool isSiblingBend(const BendPotential &bend, std::size_t currentBead, const std::vector<std::size_t> &nextBeads)
{
  if (bend.identifiers[1] != currentBead) return false;
  return contains(nextBeads, bend.identifiers[0]) && contains(nextBeads, bend.identifiers[2]);
}

// The anchor bend of a next bead within a step (previous - current - next, centered on the current
// bead), the potential that shapes its cone direction in the exact flexible base sampler.
std::optional<BendPotential> findAnchorBend(const Potentials::IntraMolecularPotentials &intra,
                                            std::size_t previousBead, std::size_t currentBead, std::size_t nextBead)
{
  for (const BendPotential &bend : intra.bends)
  {
    if (bend.identifiers[1] != currentBead) continue;
    if ((bend.identifiers[0] == previousBead && bend.identifiers[2] == nextBead) ||
        (bend.identifiers[2] == previousBead && bend.identifiers[0] == nextBead))
    {
      return bend;
    }
  }
  return std::nullopt;
}

bool hasUnsampledTerms(const Potentials::IntraMolecularPotentials &terms)
{
  return !(terms.ureyBradleys.empty() && terms.inversionBends.empty() && terms.outOfPlaneBends.empty() &&
           terms.improperTorsions.empty() && terms.bondBonds.empty() && terms.bondBends.empty() &&
           terms.bondTorsions.empty() && terms.bendBends.empty() && terms.bendTorsions.empty());
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
//    normalization (see CBMC::logBaseSamplerNormalization).
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

UnsampledStepTerms splitUnsampledStepTerms(const CBMC::GrowStep &step)
{
  UnsampledStepTerms split{};
  const std::size_t currentBead = step.currentBead;
  const std::size_t previousBead = step.previousBead.value();

  auto isGrown = [&](std::size_t id) { return contains(step.nextBeads, id); };
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

// A ring torsion transforms rigidly (is spin-invariant) when all four atoms belong to the ring body
// ('currentBead' plus the 'nextBeads'); it is spin-variant when it also involves an outside atom.
bool isSpinVariantRingTorsion(const std::array<std::size_t, 4> &identifiers, std::size_t currentBead,
                              const std::vector<std::size_t> &nextBeads)
{
  bool hasRingAtom = false;
  bool hasOutsideAtom = false;
  for (std::size_t id : identifiers)
  {
    const bool moves = contains(nextBeads, id);
    const bool inRingBody = (id == currentBead) || moves;
    if (moves) hasRingAtom = true;
    if (!inRingBody) hasOutsideAtom = true;
  }
  return hasRingAtom && hasOutsideAtom;
}

// The temperature-independent memo key of a step's base coupling: the per-bead samplers plus every
// coupling term, with atom identifiers mapped to their step-local roles so congruent steps share
// one entry (and thus one frozen estimate, which is what makes their factors cancel in reptation).
std::string baseCouplingSignature(const CBMC::GrowStep &step)
{
  const std::size_t currentBead = step.currentBead;
  const std::size_t previousBead = step.previousBead.value();

  const auto roleOf = [&](std::size_t id) -> std::string
  {
    if (id == currentBead) return "c";
    if (id == previousBead) return "p";
    return std::format(
        "n{}", std::distance(step.nextBeads.begin(), std::find(step.nextBeads.begin(), step.nextBeads.end(), id)));
  };

  std::string key{};
  for (std::size_t i = 0; i != step.nextBeads.size(); ++i)
  {
    const std::optional<BondPotential> &bond = step.nextBeadBonds[i];
    key += bond.has_value() ? std::format(";bond{}", static_cast<std::size_t>(bond->type)) : ";bond-none";
    if (bond.has_value())
      for (double p : bond->parameters) key += std::format(",{:.12e}", p);
    const std::optional<BendPotential> &anchor = step.nextBeadAnchorBends[i];
    key += anchor.has_value() ? std::format(";anchor{}", static_cast<std::size_t>(anchor->type)) : ";anchor-none";
    if (anchor.has_value())
      for (double p : anchor->parameters) key += std::format(",{:.12e}", p);
  }
  const auto appendTerm = [&](std::string_view tag, std::size_t type, const auto &identifiers, const auto &parameters)
  {
    key += std::format(";{}{}", tag, type);
    for (std::size_t id : identifiers) key += ":" + roleOf(id);
    for (double p : parameters) key += std::format(",{:.12e}", p);
  };
  for (const BendPotential &bend : step.siblingBends)
    appendTerm("sib", static_cast<std::size_t>(bend.type), bend.identifiers, bend.parameters);
  const Potentials::IntraMolecularPotentials &terms = step.baseCouplingTerms;
  for (const auto &t : terms.ureyBradleys)
    appendTerm("ub", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : terms.inversionBends)
    appendTerm("inv", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : terms.improperTorsions)
    appendTerm("imp", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : terms.bondBonds) appendTerm("bb", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : terms.bondBends) appendTerm("bB", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : terms.bondTorsions)
    appendTerm("bt", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : terms.bendBends) appendTerm("BB", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  for (const auto &t : terms.bendTorsions)
    appendTerm("Bt", static_cast<std::size_t>(t.type), t.identifiers, t.parameters);
  return key;
}

// ---------------------------------------------------------------------------------------------------
// Step-constant data of the rigid tilt and the ring-closure Monte-Carlo (formerly rebuilt per trial).
// ---------------------------------------------------------------------------------------------------

// The body atom bonded to the anchor and the junction bend previous-current-inner.
void prepareRigidTilt(CBMC::GrowStep &step, const ConnectivityTable &connectivity)
{
  if (!step.rigidBody || !step.previousBead.has_value()) return;

  const std::size_t currentBead = step.currentBead;
  const std::size_t previousBead = step.previousBead.value();

  step.rigidTilt.innerBead = step.nextBeads[0];
  for (std::size_t atom : step.nextBeads)
  {
    if (connectivity[atom, currentBead])
    {
      step.rigidTilt.innerBead = atom;
      break;
    }
  }
  step.rigidTilt.junctionBend = findAnchorBend(step.intra, previousBead, currentBead, step.rigidTilt.innerBead);
}

// Move units, fixed-bond pivots, crankshaft candidates and monitored chiral centres of a ring step.
void prepareRingSampler(CBMC::GrowStep &step, const ConnectivityTable &connectivity, const FragmentGraph &graph,
                        const std::vector<ChiralCenter> &chiralCenters)
{
  if (step.kind != CBMC::GrowStep::Kind::CloseRing) return;

  const std::size_t numberOfBeads = connectivity.numberOfBeads;
  const std::size_t currentBead = step.currentBead;
  const std::vector<std::size_t> &nextBeads = step.nextBeads;
  CBMC::GrowStep::RingSamplerData &ring = step.ring;

  // Move units: a single-atom (flexible) fragment moves by per-atom displacement, a rigid-body fragment
  // moves as one unit so its internal geometry is preserved exactly. The growth plan places a rigid
  // fragment either entirely inside this step or entirely before it, never partially.
  {
    std::map<std::size_t, std::size_t> rigidFragmentUnits{};
    for (std::size_t atom : nextBeads)
    {
      std::size_t fragmentIndex = graph.atomFragmentIds[atom];
      if (!graph.fragments[fragmentIndex].isRigidBody())
      {
        ring.moveUnits.push_back({atom});
        continue;
      }
      auto [it, inserted] = rigidFragmentUnits.insert({fragmentIndex, ring.moveUnits.size()});
      if (inserted) ring.moveUnits.push_back({});
      ring.moveUnits[it->second].push_back(atom);
    }
  }

  // Fixed-bond neighbours of each atom: a Fixed bond is a holonomic distance constraint without energy,
  // so any move of a flexible ring atom is built as a rotation about its fixed neighbour(s).
  ring.fixedNeighbors.assign(numberOfBeads, {});
  for (const BondPotential &bond : step.intra.bonds)
  {
    if (bond.type != BondType::Fixed) continue;
    ring.fixedNeighbors[bond.identifiers[0]].push_back(bond.identifiers[1]);
    ring.fixedNeighbors[bond.identifiers[1]].push_back(bond.identifiers[0]);
  }

  // Crankshaft candidates: a flexible ring atom with two positioned bonded neighbours (fixed neighbours
  // first, so a crankshaft never breaks a Fixed bond; three or more fixed bonds pin the atom).
  {
    std::vector<bool> positioned(numberOfBeads, false);
    positioned[currentBead] = true;
    if (step.previousBead.has_value()) positioned[step.previousBead.value()] = true;
    for (std::size_t atom : nextBeads) positioned[atom] = true;

    for (std::size_t atom : nextBeads)
    {
      if (graph.fragments[graph.atomFragmentIds[atom]].isRigidBody()) continue;
      const std::vector<std::size_t> &fixed = ring.fixedNeighbors[atom];
      if (fixed.size() >= 3) continue;

      std::vector<std::size_t> axisNeighbors = fixed;
      for (std::size_t other = 0; other != numberOfBeads; ++other)
      {
        if (other == atom || !positioned[other]) continue;
        if (!connectivity[atom, other]) continue;
        if (contains(fixed, other)) continue;
        axisNeighbors.push_back(other);
      }
      if (axisNeighbors.size() >= 2) ring.crankshafts.push_back({atom, axisNeighbors[0], axisNeighbors[1]});
    }
  }

  // Chiral centres whose four atoms all have positions during this step.
  {
    auto isKnown = [&](std::size_t id)
    {
      if (id == currentBead) return true;
      if (step.previousBead.has_value() && id == step.previousBead.value()) return true;
      return contains(nextBeads, id);
    };
    for (const ChiralCenter &center : chiralCenters)
    {
      if (std::ranges::all_of(center.ids, isKnown)) ring.monitoredChiralCenters.push_back(center);
    }
  }
}
}  // namespace

// Fills the derived (temperature-independent) sampler data of a step from its topology and filtered
// potentials; see the field documentation of 'GrowStep'.
void CBMC::prepareStep(GrowStep &step, const ConnectivityTable &connectivity, const FragmentGraph &fragmentGraph,
                       const std::vector<ChiralCenter> &chiralCenters)
{
  const std::size_t currentBead = step.currentBead;
  const std::vector<std::size_t> &nextBeads = step.nextBeads;

  step.flexibleAttach =
      step.kind == CBMC::GrowStep::Kind::AttachFragment && !step.rigidBody && step.previousBead.has_value();

  prepareRigidTilt(step, connectivity);
  prepareRingSampler(step, connectivity, fragmentGraph, chiralCenters);

  // Torsion-selection potentials: the junction-crossing torsions of a ring-closure step (the internal
  // ring torsions are sampled by the conformational MC and would be double counted), all torsions
  // otherwise; the spin-variant bends are weighted alongside them for every kind of step.
  if (step.kind == CBMC::GrowStep::Kind::CloseRing)
  {
    // The complement -- bonds and the spin-invariant bends and torsions -- is what the internal
    // conformational MC of the ring samples.
    step.ringInternalPotentials.bonds = step.intra.bonds;
    for (const TorsionPotential &torsion : step.intra.torsions)
    {
      (isSpinVariantRingTorsion(torsion.identifiers, currentBead, nextBeads) ? step.torsionSelectionPotentials
                                                                             : step.ringInternalPotentials)
          .torsions.push_back(torsion);
    }
  }
  else
  {
    step.torsionSelectionPotentials.torsions = step.intra.torsions;
  }
  for (const BendPotential &bend : step.intra.bends)
  {
    if (isSpinVariantBend(bend, step.previousBead, currentBead, nextBeads))
    {
      step.spinVariantBends.push_back(bend);
    }
    else if (step.kind == CBMC::GrowStep::Kind::CloseRing)
    {
      step.ringInternalPotentials.bends.push_back(bend);
    }
  }

  if (!step.previousBead.has_value()) return;
  const std::size_t previousBead = step.previousBead.value();

  step.nextBeadBonds.resize(nextBeads.size());
  step.nextBeadAnchorBends.resize(nextBeads.size());
  for (std::size_t i = 0; i != nextBeads.size(); ++i)
  {
    step.nextBeadBonds[i] = step.intra.findBondPotential(currentBead, nextBeads[i]);
    step.nextBeadAnchorBends[i] = findAnchorBend(step.intra, previousBead, currentBead, nextBeads[i]);
  }

  if (!step.flexibleAttach) return;

  // Flexible attach: the base sampler imposes the sibling bends and the base-routed unsampled terms;
  // the spin-routed unsampled terms join the torsion selection.
  for (const BendPotential &bend : step.intra.bends)
  {
    if (isSiblingBend(bend, currentBead, nextBeads)) step.siblingBends.push_back(bend);
  }

  UnsampledStepTerms split = splitUnsampledStepTerms(step);
  step.baseCouplingTerms = std::move(split.baseCoupling);
  step.hasBaseCouplingTerms = hasUnsampledTerms(step.baseCouplingTerms);
  step.hasBaseCoupling = !step.siblingBends.empty() || step.hasBaseCouplingTerms;

  Potentials::IntraMolecularPotentials &spin = step.torsionSelectionPotentials;
  spin.ureyBradleys = std::move(split.spinTerms.ureyBradleys);
  spin.inversionBends = std::move(split.spinTerms.inversionBends);
  spin.outOfPlaneBends = std::move(split.spinTerms.outOfPlaneBends);
  spin.improperTorsions = std::move(split.spinTerms.improperTorsions);
  spin.bondBonds = std::move(split.spinTerms.bondBonds);
  spin.bondBends = std::move(split.spinTerms.bondBends);
  spin.bondTorsions = std::move(split.spinTerms.bondTorsions);
  spin.bendBends = std::move(split.spinTerms.bendBends);
  spin.bendTorsions = std::move(split.spinTerms.bendTorsions);
  step.torsionSelectionHasUnsampledTerms = hasUnsampledTerms(spin);

  // Declared chiral centres fully determined by this step: centred on the current bead, with every
  // neighbour either the previous bead or grown here.
  for (const ChiralCenter &center : chiralCenters)
  {
    if (center.ids[0] != currentBead) continue;
    bool determined = true;
    for (std::size_t k = 1; k != 4; ++k)
    {
      const std::size_t id = center.ids[k];
      determined = determined && (id == previousBead || contains(nextBeads, id));
    }
    if (determined) step.determinedChiralCenters.push_back(center);
  }

  if (step.hasBaseCoupling) step.baseCouplingSignature = baseCouplingSignature(step);
}
