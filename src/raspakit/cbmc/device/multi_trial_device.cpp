module;

module cbmc_multi_trial_device;

import std;

import atom;
import molecule;
import double3;
import double3x3;
import simd_quatd;
import randomnumbers;
import component;
import forcefield;
import running_energy;
import spatial_decomposition_settings;
import cbmc;
import cbmc_external_energy;
import cbmc_multi_trial;
import cbmc_multi_trial_cpu;
import cbmc_device_environment;

namespace
{
using namespace CBMC;

/// One rigid trial molecule at a uniform random position with a uniform random orientation; 'templateAtoms'
/// supply the attributes (identity, charge, scaling), the component's reference geometry the shape.
void placeRigid(RandomNumber &random, const GrowContext &context, const Component &component,
                std::span<const Atom> templateAtoms, std::span<Atom> out)
{
  const double3 center = context.simulationBox.randomPosition(random);
  if (out.size() == 1)
  {
    out[0] = templateAtoms[0];
    out[0].position = center;
    return;
  }
  const simd_quatd orientation = random.randomSimdQuatd();
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(orientation);
  for (std::size_t i = 0; i < out.size(); ++i)
  {
    out[i] = templateAtoms[i];
    out[i].position = center + rotation * component.atoms[i].position;
  }
}

bool unitScaling(std::span<const Atom> atoms)
{
  return std::ranges::all_of(atoms, [](const Atom &atom)
                             { return atom.scalingVDW == 1.0 && atom.scalingCoulomb == 1.0 && !atom.isFractional; });
}

/// A group of k trials after the device evaluation: its weight, the selected trial, and its long-range energy.
struct GroupSelection
{
  std::size_t selectedSet;  ///< index into the batch
  double logWeight;         ///< log W_i = log( (1/k) sum_j exp(-beta U_ij) )
  double longRangeEnergy;   ///< dU_lr of the selected trial
  double logWeightTilde;    ///< log W~_i = log W_i - beta dU_lr
};

/// The non-overlapping trials of the group of sets [first, first + count) with their log Boltzmann factors.
struct GroupMembers
{
  std::vector<double> logBoltzmann;
  std::vector<std::size_t> sets;
};

GroupMembers collectGroup(double beta, std::size_t first, std::size_t count, std::span<const DeviceTrialEnergy> energies,
                          std::span<const std::uint8_t> invalid)
{
  GroupMembers members;
  members.logBoltzmann.reserve(count);
  members.sets.reserve(count);
  for (std::size_t t = 0; t < count; ++t)
  {
    const std::size_t set = first + t;
    if (set >= energies.size() || invalid[set] || energies[set].overlap) continue;
    members.logBoltzmann.push_back(-beta * energies[set].total());
    members.sets.push_back(set);
  }
  return members;
}

/// log W = log( (1/k) sum_j exp(-beta U_j) ) of a group with at least one member.
double groupLogWeight(const GroupMembers &members, std::size_t k)
{
  return logSumExp(members.logBoltzmann) - std::log(static_cast<double>(k));
}

/// Rosenbluth selection within a fresh group of k trials with the device energies, and the long-range energy of
/// the selected trial; std::nullopt when every trial overlaps.
std::optional<GroupSelection> selectInGroup(RandomNumber &random, double beta, std::size_t first, std::size_t k,
                                            std::span<const DeviceTrialEnergy> energies,
                                            std::span<const std::uint8_t> invalid, std::span<const Atom> batch,
                                            std::size_t atomsPerSet, const LongRangeCorrection &longRange)
{
  const GroupMembers members = collectGroup(beta, first, k, energies, invalid);
  if (members.sets.empty()) return std::nullopt;

  const double logWeight = groupLogWeight(members, k);
  const std::size_t selectedSet = members.sets[selectByLogWeight(random, members.logBoltzmann)];
  const double longRangeEnergy = longRange(batch.subspan(selectedSet * atomsPerSet, atomsPerSet));
  return GroupSelection{.selectedSet = selectedSet,
                        .logWeight = logWeight,
                        .longRangeEnergy = longRangeEnergy,
                        .logWeightTilde = logWeight - beta * longRangeEnergy};
}

std::vector<std::uint8_t> presentTypes(const GrowContext &context, const Component &component)
{
  std::vector<std::uint8_t> present(context.forceField.pseudoAtoms.size(), 0);
  auto mark = [&](std::span<const Atom> atoms)
  {
    for (const Atom &atom : atoms)
    {
      if (atom.type < present.size()) present[atom.type] = 1;
    }
  };
  mark(context.frameworkAtoms);
  mark(context.moleculeAtoms);
  mark(component.atoms);
  return present;
}
}  // namespace

CBMC::DeviceMultiTrialGrower::DeviceMultiTrialGrower(PairDevice device)
    : environment(std::make_unique<DeviceGrowEnvironment>(device))
{
}

PairDevice CBMC::DeviceMultiTrialGrower::device() const { return environment->device(); }

std::string CBMC::DeviceMultiTrialGrower::deviceName() const { return environment->deviceName(); }

bool CBMC::DeviceMultiTrialGrower::supports(const GrowContext &context, const Component &component,
                                            std::string &reason)
{
  if (!component.rigid || component.atoms.empty())
  {
    reason = std::format("component '{}' is not rigid; the device grower places rigid molecules only",
                         component.name);
    return false;
  }
  if (component.atoms.size() > DeviceGrowEnvironment::maximumAtomsPerSet)
  {
    reason = std::format("component '{}' has {} atoms, the device grower holds at most {}", component.name,
                         component.atoms.size(), DeviceGrowEnvironment::maximumAtomsPerSet);
    return false;
  }
  if (context.hasExternalField)
  {
    reason = "the device grower does not evaluate external fields";
    return false;
  }
  if (!context.crossLinkTethers.empty())
  {
    reason = "the device grower does not evaluate cross-link tethers";
    return false;
  }
  if (!unitScaling(component.atoms))
  {
    reason = std::format("component '{}' has scaled (fractional) atoms", component.name);
    return false;
  }
  const std::vector<std::uint8_t> present = presentTypes(context, component);
  return DeviceGrowEnvironment::supportsForceField(context.forceField, present, reason);
}

std::optional<CBMC::MultiTrialGrowResult> CBMC::DeviceMultiTrialGrower::growBatch(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<const Atom> templateAtoms,
    std::size_t numberOfTrials, const LongRangeCorrection &longRange) const
{
  const std::size_t M = std::max<std::size_t>(1, numberOfTrials);
  const std::size_t k = std::max<std::size_t>(1, context.settings.numberOfTrialDirections);
  const std::size_t n = component.atoms.size();
  const double beta = context.beta;

  std::vector<Atom> batch(M * k * n);
  std::vector<std::uint8_t> invalid(M * k, 0);
  for (std::size_t set = 0; set < M * k; ++set)
  {
    const std::span<Atom> atoms(batch.data() + set * n, n);
    placeRigid(random, context, component, templateAtoms, atoms);
    if (insideBlockedPockets(context.framework, component, atoms)) invalid[set] = 1;
  }

  std::vector<DeviceTrialEnergy> energies;
  environment->setForceField(context.forceField, context.simulationBox);
  environment->setBackground(context.frameworkAtoms, context.moleculeAtoms);
  environment->evaluate(batch, n, templateAtoms.front().moleculeId, context.skipBackgroundMolecule, energies);

  std::vector<GroupSelection> groups;
  std::vector<double> logWeightsTilde;
  groups.reserve(M);
  logWeightsTilde.reserve(M);
  for (std::size_t g = 0; g < M; ++g)
  {
    std::optional<GroupSelection> selection =
        selectInGroup(random, beta, g * k, k, energies, invalid, batch, n, longRange);
    if (!selection.has_value()) continue;
    logWeightsTilde.push_back(selection->logWeightTilde);
    groups.push_back(selection.value());
  }
  if (groups.empty()) return std::nullopt;

  const GroupSelection &selected = groups[selectByLogWeight(random, logWeightsTilde)];
  const std::span<const Atom> selectedAtoms(batch.data() + selected.selectedSet * n, n);

  // The host energies of the selected trial, in double at the full cut-offs: the running energies receive these.
  // A trial that overlaps in double has zero weight in the target distribution: the move is rejected (the
  // correction factor would be zero).
  const std::optional<RunningEnergy> hostEnergy =
      computeExternalNonOverlappingEnergy(context.withFullCutOffs(), component, selectedAtoms);
  if (!hostEnergy.has_value()) return std::nullopt;

  const double deviceEnergy = energies[selected.selectedSet].total();
  std::vector<Atom> atoms(selectedAtoms.begin(), selectedAtoms.end());
  GrowResult grow(component.createMoleculeRecord(atoms), std::move(atoms), hostEnergy.value(), selected.logWeight);
  return MultiTrialGrowResult{.selected = std::move(grow),
                              .selectedLongRangeEnergy = selected.longRangeEnergy,
                              .logMeanWeight = logSumExp(logWeightsTilde) - std::log(static_cast<double>(M)),
                              .numberOfConstructed = groups.size(),
                              .logProposalCorrection = -beta * (hostEnergy->potentialEnergy() - deviceEnergy)};
}

CBMC::MultiTrialRetraceResult CBMC::DeviceMultiTrialGrower::retraceBatch(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<const Atom> moleculeAtoms,
    std::span<const Atom> templateAtoms, std::size_t numberOfTrials, const LongRangeCorrection &longRange) const
{
  const std::size_t M = std::max<std::size_t>(1, numberOfTrials);
  const std::size_t k = std::max<std::size_t>(1, context.settings.numberOfTrialDirections);
  const std::size_t n = component.atoms.size();
  const double beta = context.beta;

  // set 0 is the existing molecule; sets 1..k-1 complete its group; the M-1 reference groups follow
  std::vector<Atom> batch(M * k * n);
  std::vector<std::uint8_t> invalid(M * k, 0);
  std::copy(moleculeAtoms.begin(), moleculeAtoms.end(), batch.begin());
  for (std::size_t set = 1; set < M * k; ++set)
  {
    const std::span<Atom> atoms(batch.data() + set * n, n);
    placeRigid(random, context, component, templateAtoms, atoms);
    if (insideBlockedPockets(context.framework, component, atoms)) invalid[set] = 1;
  }

  std::vector<DeviceTrialEnergy> energies;
  environment->setForceField(context.forceField, context.simulationBox);
  environment->setBackground(context.frameworkAtoms, context.moleculeAtoms);
  environment->evaluate(batch, n, moleculeAtoms.front().moleculeId, context.skipBackgroundMolecule, energies);

  // The existing molecule always counts (an overlap flagged by the device for an accepted configuration can only
  // be a borderline rounding case: it is given the overlap criterion as its energy, i.e. a finite tiny weight).
  const DeviceTrialEnergy &oldDevice = energies[0];
  const double oldDeviceEnergy =
      oldDevice.overlap ? std::max(oldDevice.total(), context.forceField.energyOverlapCriteria) : oldDevice.total();
  GroupMembers oldGroup = collectGroup(beta, 1, k - 1, energies, invalid);
  oldGroup.logBoltzmann.insert(oldGroup.logBoltzmann.begin(), -beta * oldDeviceEnergy);
  const double logWeightOld = groupLogWeight(oldGroup, k);
  // no selection on the retrace side: the weight is the group's, the long-range energy the existing molecule's
  const double oldLongRange = longRange(moleculeAtoms);
  const double logWeightTildeOld = logWeightOld - beta * oldLongRange;

  std::vector<double> logWeightsTilde;
  logWeightsTilde.reserve(M);
  logWeightsTilde.push_back(logWeightTildeOld);
  std::size_t constructed = 0;
  for (std::size_t g = 1; g < M; ++g)
  {
    std::optional<GroupSelection> selection =
        selectInGroup(random, beta, g * k, k, energies, invalid, batch, n, longRange);
    if (!selection.has_value()) continue;
    ++constructed;
    logWeightsTilde.push_back(selection->logWeightTilde);
  }

  const std::optional<RunningEnergy> hostEnergy =
      computeExternalNonOverlappingEnergy(context.withFullCutOffs(), component, moleculeAtoms);
  if (!hostEnergy.has_value())
  {
    throw std::runtime_error(std::format(
        "CBMC device grower: the existing configuration of component '{}' (molecule {}) overlaps at the full "
        "cut-offs; the retrace of an overlapping molecule has no defined weight. The simulation state is "
        "inconsistent.\n",
        component.name, moleculeAtoms.front().moleculeId));
  }

  return MultiTrialRetraceResult{.old = RetraceResult(hostEnergy.value(), logWeightOld),
                                 .oldLongRangeEnergy = oldLongRange,
                                 .logMeanWeight = logSumExp(logWeightsTilde) - std::log(static_cast<double>(M)),
                                 .numberOfConstructed = constructed,
                                 .logProposalCorrection = beta * (hostEnergy->potentialEnergy() - oldDeviceEnergy)};
}

std::optional<CBMC::MultiTrialGrowResult> CBMC::DeviceMultiTrialGrower::growNew(
    RandomNumber &random, const GrowContext &context, const Component &component, const NewMoleculeIdentity &identity,
    std::size_t numberOfTrials, const LongRangeCorrection &longRange) const
{
  std::string reason;
  if (identity.scaling != 1.0 || identity.isFractional || !supports(context, component, reason))
  {
    return cpuMultiTrialGrower().growNew(random, context, component, identity, numberOfTrials, longRange);
  }
  std::vector<Atom> templateAtoms = component.atoms;
  for (Atom &atom : templateAtoms)
  {
    atom.moleculeId = static_cast<std::uint32_t>(identity.moleculeId);
    atom.groupId = identity.groupId;
    atom.isFractional = identity.isFractional;
    atom.setScaling(identity.scaling);
  }
  return growBatch(random, context, component, templateAtoms, numberOfTrials, longRange);
}

CBMC::MultiTrialRetraceResult CBMC::DeviceMultiTrialGrower::retraceNew(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<const Atom> moleculeAtoms,
    const NewMoleculeIdentity &referenceIdentity, std::size_t numberOfTrials,
    const LongRangeCorrection &longRange) const
{
  std::string reason;
  if (referenceIdentity.scaling != 1.0 || referenceIdentity.isFractional || !unitScaling(moleculeAtoms) ||
      moleculeAtoms.size() != component.atoms.size() || !supports(context, component, reason))
  {
    return cpuMultiTrialGrower().retraceNew(random, context, component, moleculeAtoms, referenceIdentity,
                                            numberOfTrials, longRange);
  }
  std::vector<Atom> templateAtoms = component.atoms;
  for (Atom &atom : templateAtoms)
  {
    atom.moleculeId = static_cast<std::uint32_t>(referenceIdentity.moleculeId);
    atom.groupId = referenceIdentity.groupId;
    atom.isFractional = referenceIdentity.isFractional;
    atom.setScaling(referenceIdentity.scaling);
  }
  return retraceBatch(random, context, component, moleculeAtoms, templateAtoms, numberOfTrials, longRange);
}

std::optional<CBMC::MultiTrialGrowResult> CBMC::DeviceMultiTrialGrower::regrow(
    RandomNumber &random, const GrowContext &context, const Component &component, const Molecule &molecule,
    std::span<const Atom> moleculeAtoms, std::size_t numberOfTrials, const LongRangeCorrection &longRange) const
{
  std::string reason;
  if (!unitScaling(moleculeAtoms) || moleculeAtoms.size() != component.atoms.size() ||
      !supports(context, component, reason))
  {
    return cpuMultiTrialGrower().regrow(random, context, component, molecule, moleculeAtoms, numberOfTrials,
                                        longRange);
  }
  return growBatch(random, context, component, moleculeAtoms, numberOfTrials, longRange);
}

CBMC::MultiTrialRetraceResult CBMC::DeviceMultiTrialGrower::retraceRegrow(
    RandomNumber &random, const GrowContext &context, const Component &component, const Molecule &molecule,
    std::span<const Atom> moleculeAtoms, double storedR, std::size_t numberOfTrials,
    const LongRangeCorrection &longRange) const
{
  std::string reason;
  if (!unitScaling(moleculeAtoms) || moleculeAtoms.size() != component.atoms.size() ||
      !supports(context, component, reason))
  {
    return cpuMultiTrialGrower().retraceRegrow(random, context, component, molecule, moleculeAtoms, storedR,
                                               numberOfTrials, longRange);
  }
  // the one-stage rigid scheme has no retained first-bead weight: 'storedR' is not used
  return retraceBatch(random, context, component, moleculeAtoms, moleculeAtoms, numberOfTrials, longRange);
}

const CBMC::MultiTrialGrower &CBMC::deviceMultiTrialGrower(PairDevice device)
{
  static std::map<PairDevice, std::unique_ptr<DeviceMultiTrialGrower>> growers;
  auto it = growers.find(device);
  if (it == growers.end())
  {
    it = growers.emplace(device, std::make_unique<DeviceMultiTrialGrower>(device)).first;
  }
  return *it->second;
}
