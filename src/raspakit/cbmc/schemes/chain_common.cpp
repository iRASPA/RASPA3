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
    if (lastStep.has_value()) schedule.perStep[lastStep.value()].push_back(index);
  }
  return schedule;
}

std::optional<CBMC::StepTrialEnergy> CBMC::evaluateStepTrial(const GrowContext &context, const Component &component,
                                                             const GrowStep &step, std::vector<Atom> &chainAtoms,
                                                             std::span<const Atom> positions,
                                                             std::span<const std::size_t> tethersOfStep)
{
  std::optional<RunningEnergy> external = computeExternalNonOverlappingEnergy(context, component, positions);
  if (!external.has_value()) return std::nullopt;

  // Only the external-stage share of the step's intramolecular non-bonded pairs; the short-range
  // pairs were weighted in the torsion-spin selection that generated 'positions' (GrowStep::NonBondedData).
  const ScratchBeads scratch(chainAtoms, step);
  placeStepBeads(chainAtoms, step, positions);
  const RunningEnergy intra = step.nonBonded.external.computeInternalIntraVanDerWaalsAndCoulombEnergies(chainAtoms);

  // The cross-link terms that become evaluable at this step (every atom they need is now in place).
  if (!tethersOfStep.empty())
  {
    external.value() += Interactions::computeCrossLinkTetherEnergy(
        context.forceField, context.simulationBox, component, chainAtoms, context.crossLinkTethers, tethersOfStep);
  }

  return StepTrialEnergy{external.value(), intra};
}

double CBMC::unsampledInternalLogFactor(double beta, const GrowStep &step, const std::vector<Atom> &chainAtoms)
{
  if (stepHandlesUnsampledInternalTerms(step)) return 0.0;
  return -beta * step.intra.computeInternalEnergiesNotSampledDuringGrowth(chainAtoms).potentialEnergy();
}

bool CBMC::ChainAccumulator::addGrownStep(double stepLogWeight, const RunningEnergy &stepExternalEnergy,
                                          double minimumRosenbluthFactor)
{
  if (std::exp(stepLogWeight) < minimumRosenbluthFactor) return false;
  logRosenbluthWeight += stepLogWeight;
  externalEnergies += stepExternalEnergy;
  return true;
}

void CBMC::ChainAccumulator::addRetracedStep(double stepLogWeight, const RunningEnergy &stepExternalEnergy)
{
  logRosenbluthWeight += stepLogWeight;
  externalEnergies += stepExternalEnergy;
}

CBMC::GrowResult CBMC::finishGrownChain(const Component &component, std::vector<Atom> chainAtoms,
                                        const ChainAccumulator &accumulated)
{
  const RunningEnergy internalEnergies = component.intraMolecularPotentials.computeInternalEnergies(chainAtoms);
  const Molecule molecule = component.createMoleculeRecord(chainAtoms);
  return GrowResult(molecule, std::move(chainAtoms), accumulated.externalEnergies + internalEnergies,
                    accumulated.logRosenbluthWeight);
}

CBMC::RetraceResult CBMC::finishRetracedChain(const Component &component, std::span<const Atom> chainAtoms,
                                              const ChainAccumulator &accumulated)
{
  const RunningEnergy internalEnergies = component.intraMolecularPotentials.computeInternalEnergies(chainAtoms);
  return RetraceResult(accumulated.externalEnergies + internalEnergies, accumulated.logRosenbluthWeight);
}

void CBMC::throwExistingConfigurationOverlaps(std::string_view scheme, const Component &component,
                                              std::size_t stepIndex, const GrowStep &step)
{
  std::string beads{};
  for (std::size_t bead : step.nextBeads) beads += std::format(" {}", bead);
  throw std::runtime_error(std::format(
      "{}: the existing configuration of component '{}' overlaps at growth step {} (bead(s){}); the retrace of an "
      "overlapping molecule has no defined weight. The simulation state is inconsistent (overlapping molecules in "
      "the initial/restart configuration, a molecule inside a blocked pocket, or a force field or scaling changed "
      "after placement).\n",
      scheme, component.name, stepIndex, beads));
}
