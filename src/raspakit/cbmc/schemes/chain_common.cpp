module;

module cbmc_chain_common;

import std;

import atom;
import molecule;
import component;
import running_energy;
import intra_molecular_potentials;
import cbmc_util;
import cbmc_results;
import cbmc_grow_context;
import cbmc_grow_step;
import cbmc_operators;
import cbmc_external_energy;
import interactions_cross_link;
import cross_links;

CBMC::TetherSchedule CBMC::scheduleTethers(const GrowContext &context, const Component &component,
                                           const std::vector<GrowStep> &plan)
{
  TetherSchedule schedule{};
  if (context.crossLinkTethers.empty()) return schedule;
  schedule.perStep.resize(plan.size());

  // The step that places each atom; atoms placed by no step are the pre-placed set.
  std::vector<std::optional<std::size_t>> stepOfAtom(component.atoms.size());
  for (std::size_t seg = 0; seg != plan.size(); ++seg)
  {
    for (std::size_t bead : plan[seg].nextBeads) stepOfAtom[bead] = seg;
  }

  for (std::size_t index = 0; index != context.crossLinkTethers.size(); ++index)
  {
    const std::size_t site = context.crossLinkTethers[index].siteAtom;
    std::optional<std::size_t> lastStep = stepOfAtom[site];
    if (site < component.connectivityTable.numberOfBeads)
    {
      for (std::size_t neighbour : component.connectivityTable.findAllNeighbors(site))
      {
        if (stepOfAtom[neighbour].has_value() && (!lastStep.has_value() || stepOfAtom[neighbour] > lastStep))
        {
          lastStep = stepOfAtom[neighbour];
        }
      }
    }
    if (lastStep.has_value()) schedule.perStep[lastStep.value()].tethers.push_back(index);

    // the excluded pairs of the regrown atoms, exempt from the overlap test at the step that places the atom
    for (const CrossLinkTether::Exclusion &exclusion : context.crossLinkTethers[index].exclusions)
    {
      if (!stepOfAtom[exclusion.atom].has_value()) continue;
      const std::size_t seg = stepOfAtom[exclusion.atom].value();
      const auto it = std::ranges::find(plan[seg].nextBeads, exclusion.atom);
      schedule.perStep[seg].overlapExempt.push_back(
          {static_cast<std::size_t>(std::distance(plan[seg].nextBeads.begin(), it)), exclusion.partnerIndex});

      // the exclusion goes with the step that places the atom (once per atom of a tether)
      if (lastStep.has_value() && seg < lastStep.value())
      {
        const std::pair<std::size_t, std::size_t> entry{index, exclusion.atom};
        if (std::ranges::find(schedule.perStep[seg].earlyExclusions, entry) ==
            schedule.perStep[seg].earlyExclusions.end())
        {
          schedule.perStep[seg].earlyExclusions.push_back(entry);
          schedule.perStep[lastStep.value()].exclusionsTakenEarlier.push_back(entry);
        }
      }
    }
  }
  return schedule;
}

std::optional<CBMC::StepTrialEnergy> CBMC::evaluateStepTrial(const GrowContext &context, const Component &component,
                                                             const GrowStep &step, std::vector<Atom> &chainAtoms,
                                                             std::span<const Atom> positions,
                                                             const StepTethers &stepTethers)
{
  std::optional<RunningEnergy> external =
      computeExternalNonOverlappingEnergy(context, component, positions, stepTethers.overlapExempt);
  if (!external.has_value()) return std::nullopt;

  // Only the external-stage share of the step's intramolecular non-bonded pairs; the short-range
  // pairs were weighted in the torsion-spin selection that generated 'positions' (GrowStep::NonBondedData).
  const ScratchBeads scratch(chainAtoms, step);
  placeStepBeads(chainAtoms, step, positions);
  const RunningEnergy intra = step.nonBonded.external.computeInternalIntraVanDerWaalsAndCoulombEnergies(
      context.forceField, context.simulationBox, chainAtoms);

  // The cross-link terms that become evaluable at this step (every atom they need is now in place).
  if (!stepTethers.tethers.empty())
  {
    external.value() += Interactions::computeCrossLinkTetherEnergy(context.forceField, context.simulationBox,
                                                                    component, chainAtoms, context.crossLinkTethers,
                                                                    stepTethers.tethers,
                                                                    stepTethers.exclusionsTakenEarlier);
  }
  if (!stepTethers.earlyExclusions.empty())
  {
    external.value() += Interactions::computeCrossLinkTetherExclusionEnergy(
        context.forceField, context.simulationBox, chainAtoms, context.crossLinkTethers,
        stepTethers.earlyExclusions);
  }

  return StepTrialEnergy{external.value(), intra};
}

double CBMC::unsampledInternalLogFactor(double beta, const GrowStep &step, const std::vector<Atom> &chainAtoms)
{
  if (stepHandlesUnsampledInternalTerms(step)) return 0.0;
  // The terms among the beads of the seed fragment are a constant of its (rigid or single-bead) geometry: they cancel
  // against the single-molecule partition function and are left out, as the conventional swap and the device
  // multi-trial grower leave them out. A rigid component has no ideal-gas Rosenbluth weight to cancel them, and its
  // CBMC insertions were weighted by exp(-beta U) of its constant internal terms.
  //
  // That holds for the terms among the beads of the seed fragment only. In a fixed-endpoint regrowth the placed
  // beads are not one connected set: a seed grown from a fixed bead without a placed neighbour (chain 0-1-2-3
  // with the ends fixed, bead 1 grown from bead 0) completes the terms that reach another fixed bead (a
  // Urey-Bradley 1-3, an improper or a cross term through bead 3). Those depend on where the seed is placed,
  // no later step carries them and the seed sampler does not sample them: they are weighted here.
  if (step.kind == GrowStep::Kind::PlaceSeedFragment)
  {
    const auto insideSeed = [&step](std::size_t bead)
    { return bead == step.currentBead || std::ranges::find(step.nextBeads, bead) != step.nextBeads.end(); };
    const auto reachingOutside = [&insideSeed](const auto &potentials)
    {
      std::remove_cvref_t<decltype(potentials)> selected{};
      for (const auto &potential : potentials)
      {
        if (!std::ranges::all_of(potential.identifiers, insideSeed)) selected.push_back(potential);
      }
      return selected;
    };

    Potentials::IntraMolecularPotentials outside{};
    outside.ureyBradleys = reachingOutside(step.intra.ureyBradleys);
    outside.inversionBends = reachingOutside(step.intra.inversionBends);
    outside.outOfPlaneBends = reachingOutside(step.intra.outOfPlaneBends);
    outside.improperTorsions = reachingOutside(step.intra.improperTorsions);
    outside.bondBonds = reachingOutside(step.intra.bondBonds);
    outside.bondBends = reachingOutside(step.intra.bondBends);
    outside.bondTorsions = reachingOutside(step.intra.bondTorsions);
    outside.bendBends = reachingOutside(step.intra.bendBends);
    outside.bendTorsions = reachingOutside(step.intra.bendTorsions);
    outside.cmaps = reachingOutside(step.intra.cmaps);
    if (!outside.cmaps.empty()) outside.cmapMaps = step.intra.cmapMaps;
    return -beta * outside.computeInternalEnergiesNotSampledDuringGrowth(chainAtoms).potentialEnergy();
  }
  return -beta * step.intra.computeInternalEnergiesNotSampledDuringGrowth(chainAtoms).potentialEnergy();
}

bool CBMC::ChainAccumulator::addGrownStep(double stepLogWeight, const RunningEnergy &stepExternalEnergy,
                                          double minimumRosenbluthFactor)
{
  // compared as logarithms: the exponential of a large step weight overflows
  if (minimumRosenbluthFactor > 0.0 && stepLogWeight < std::log(minimumRosenbluthFactor)) return false;
  logRosenbluthWeight += stepLogWeight;
  externalEnergies += stepExternalEnergy;
  return true;
}

void CBMC::ChainAccumulator::addRetracedStep(double stepLogWeight, const RunningEnergy &stepExternalEnergy)
{
  logRosenbluthWeight += stepLogWeight;
  externalEnergies += stepExternalEnergy;
}

CBMC::GrowResult CBMC::finishGrownChain(const GrowContext &context, const Component &component,
                                        std::vector<Atom> chainAtoms, const ChainAccumulator &accumulated)
{
  const RunningEnergy internalEnergies = component.intraMolecularPotentials.computeInternalEnergies(
      context.forceField, context.simulationBox, chainAtoms);
  const Molecule molecule = component.createMoleculeRecord(chainAtoms);
  return GrowResult(molecule, std::move(chainAtoms), accumulated.externalEnergies + internalEnergies,
                    accumulated.logRosenbluthWeight);
}

CBMC::RetraceResult CBMC::finishRetracedChain(const GrowContext &context, const Component &component,
                                              std::span<const Atom> chainAtoms, const ChainAccumulator &accumulated)
{
  const RunningEnergy internalEnergies = component.intraMolecularPotentials.computeInternalEnergies(
      context.forceField, context.simulationBox, chainAtoms);
  return RetraceResult(accumulated.externalEnergies + internalEnergies, accumulated.logRosenbluthWeight);
}

void CBMC::throwExistingConfigurationOverlaps(std::string_view scheme, const Component &component,
                                              std::size_t stepIndex, const GrowStep &step)
{
  std::string beads{};
  for (std::size_t bead : step.nextBeads) beads += std::format(" {}", bead);
  throw std::runtime_error(std::format(
      "{}: the existing configuration of component '{}' overlaps at growth step {} (bead(s){}), or its geometry "
      "there is outside what the step generates (a declared chiral centre of the other parity, a degenerate "
      "bridge). The retrace of such a molecule has no defined weight. The simulation state is "
      "inconsistent (overlapping molecules in the initial/restart configuration, a molecule inside a blocked "
      "pocket, an inverted stereocentre, or a force field or scaling changed after placement).\n",
      scheme, component.name, stepIndex, beads));
}
