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

std::optional<CBMC::StepTrialEnergy> CBMC::evaluateStepTrial(const GrowContext &context, const Component &component,
                                                             const GrowStep &step, std::vector<Atom> &chainAtoms,
                                                             std::span<const Atom> positions,
                                                             std::optional<std::size_t> skipBackgroundMolecule)
{
  const std::optional<RunningEnergy> external =
      computeExternalNonOverlappingEnergy(context, component, positions, skipBackgroundMolecule);
  if (!external.has_value()) return std::nullopt;

  const ScratchBeads scratch(chainAtoms, step);
  placeStepBeads(chainAtoms, step, positions);
  const RunningEnergy intra = step.intra.computeInternalIntraVanDerWaalsAndCoulombEnergies(chainAtoms);

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
