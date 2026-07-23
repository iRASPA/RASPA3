module;

module mc_moves_lambda_exchange;

import std;

import component;
import atom;
import double3;
import simulationbox;
import randomnumbers;
import system;
import property_lambda_probability_histogram;
import running_energy;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

// The lambda-exchange is a Hamiltonian parallel-tempering move: positions are untouched and only
// the interaction scaling of the pinned fractional molecule(s) changes. The energy-difference
// machinery below mirrors the lambda-change branches of the CFCMC swap moves (swap_cfcmc,
// pair_swap_cfcmc, group_swap_cfcmc): molecules are rescaled in place one at a time, so the
// cross-interaction between two fractional molecules of the same replica is counted exactly once.

// the component whose lambda is pinned for thermodynamic integration (exactly one per replica)
static std::optional<std::size_t> fixedLambdaComponent(const System &system)
{
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (system.components[componentId].fixedLambdaBin.has_value())
    {
      return componentId;
    }
  }
  return std::nullopt;
}

// all fractional molecules coupled to the pinned lambda coordinate: the single GC fractional
// molecule, both molecules of an ion-pair, or all members of a CFCMC group
static std::vector<std::pair<std::size_t, std::size_t>> pinnedFractionalMolecules(System &system,
                                                                                  std::size_t tiComponentId)
{
  std::vector<std::pair<std::size_t, std::size_t>> molecules;
  const Component &component = system.components[tiComponentId];

  switch (component.fixedLambdaCoordinate)
  {
    case Component::FixedLambdaCoordinate::GC:
    {
      molecules.emplace_back(tiComponentId, system.indexOfGCFractionalMoleculesPerComponent_CFCMC(tiComponentId));
      break;
    }
    case Component::FixedLambdaCoordinate::PairSwap:
    {
      const std::size_t partner = component.pairComponentId.value();
      molecules.emplace_back(tiComponentId,
                             system.indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(tiComponentId));
      molecules.emplace_back(partner, system.indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(partner));
      break;
    }
    case Component::FixedLambdaCoordinate::GroupSwap:
    {
      for (std::size_t member = 0; member < system.components.size(); ++member)
      {
        if (system.groupSwapLambdaDriver(member, Move::Types::GroupSwapCFCMC) != tiComponentId)
        {
          continue;
        }
        const std::size_t begin = system.indexOfGroupSwapFractionalMoleculesPerComponent_CFCMC(member);
        for (std::size_t k = 0; k < system.numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC[member]; ++k)
        {
          molecules.emplace_back(member, begin + k);
        }
      }
      break;
    }
  }
  return molecules;
}

// difference in scaled net charge when 'oldAtoms' is replaced by 'newAtoms'
static double scaledChargeDifference(std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms)
{
  double charge = 0.0;
  for (const Atom &atom : newAtoms) charge += atom.scalingCoulomb * atom.charge;
  for (const Atom &atom : oldAtoms) charge -= atom.scalingCoulomb * atom.charge;
  return charge;
}

// per-group summed (unscaled) charge of the dU/dlambda group-tagged atoms
static std::array<double, maximumNumberOfDUDlambdaGroups> groupChargeSum(std::span<const Atom> atoms)
{
  std::array<double, maximumNumberOfDUDlambdaGroups> charge{};
  for (const Atom &atom : atoms)
  {
    if (atom.groupId != 0) charge[atom.groupId - 1] += atom.charge;
  }
  return charge;
}

// external-field + framework + intermolecular energy difference for replacing 'oldAtoms' by 'newAtoms'
static std::optional<RunningEnergy> computeNonEwaldEnergyDifference(System &system, std::span<const Atom> newAtoms,
                                                                    std::span<const Atom> oldAtoms)
{
  std::optional<RunningEnergy> externalFieldDifference = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid, newAtoms,
      oldAtoms);
  if (!externalFieldDifference.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
  if (!frameworkDifference.has_value()) return std::nullopt;

  std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), newAtoms, oldAtoms);
  if (!moleculeDifference.has_value()) return std::nullopt;

  return externalFieldDifference.value() + frameworkDifference.value() + moleculeDifference.value();
}

// Brick-CFCMC-style aggregated tail-correction difference with running effective type counts
static RunningEnergy computeTailEnergyDifference(
    System &system, std::vector<double> &tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> &tailGroupCounts, std::span<const Atom> newAtoms,
    std::span<const Atom> oldAtoms)
{
  RunningEnergy result =
      Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
          system.forceField, system.simulationBox, tailEffectiveCounts, tailGroupCounts, newAtoms, oldAtoms) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
  Interactions::updateEffectiveTypeCounts(tailEffectiveCounts, tailGroupCounts, newAtoms, oldAtoms);
  return result;
}

// trial state of one replica: pinned fractional molecules rescaled in place to a new lambda
struct LambdaRescaleTrial
{
  std::vector<std::span<Atom>> molecules;               // spans into the system's atom storage (rescaled)
  std::vector<std::vector<Atom>> oldAtoms;              // saved copies at the original lambda
  std::vector<double3> electricFieldNeighborDelta;      // polarization field change on the other molecules
  RunningEnergy energyDifference;                       // U(newLambda) - U(oldLambda)
};

static void restoreLambdaRescale(LambdaRescaleTrial &trial)
{
  for (std::size_t i = 0; i < trial.molecules.size(); ++i)
  {
    std::copy(trial.oldAtoms[i].begin(), trial.oldAtoms[i].end(), trial.molecules[i].begin());
  }
}

// Rescale all pinned fractional molecules of 'system' to 'newLambda' and compute the total energy
// difference. On return the atoms are already at the new scaling and system.trialEik holds the new
// Ewald state; the caller either accepts (acceptEwaldMove + tail-count refresh) or restores.
static std::optional<LambdaRescaleTrial> rescaleToLambda(System &system, std::size_t tiComponentId, double newLambda)
{
  LambdaRescaleTrial trial{};

  const std::vector<std::pair<std::size_t, std::size_t>> moleculeIds = pinnedFractionalMolecules(system, tiComponentId);
  trial.molecules.reserve(moleculeIds.size());
  trial.oldAtoms.reserve(moleculeIds.size());
  for (const auto &[componentId, moleculeId] : moleculeIds)
  {
    std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
    trial.molecules.push_back(molecule);
    trial.oldAtoms.emplace_back(molecule.begin(), molecule.end());
  }

  // group-tagged charge of all pinned molecules; per-molecule the 'external' part is the total
  // minus its own contribution (needed for the dU/dlambda-consistent Ewald net-charge correction)
  std::array<double, maximumNumberOfDUDlambdaGroups> totalGroupCharge{};
  for (const std::vector<Atom> &oldMolecule : trial.oldAtoms)
  {
    const std::array<double, maximumNumberOfDUDlambdaGroups> own = groupChargeSum(oldMolecule);
    for (std::size_t g = 0; g < maximumNumberOfDUDlambdaGroups; ++g) totalGroupCharge[g] += own[g];
  }

  std::vector<double> tailEffectiveCounts = system.effectiveNumberOfPseudoAtomsVDW;
  std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> tailGroupCounts =
      system.fractionalPseudoAtomCountsPerGroup;
  double runningNetCharge = system.netCharge;

  for (std::size_t i = 0; i < trial.molecules.size(); ++i)
  {
    std::span<Atom> molecule = trial.molecules[i];
    const std::vector<Atom> &oldMolecule = trial.oldAtoms[i];

    for (Atom &atom : molecule) atom.setScaling(newLambda);

    std::optional<RunningEnergy> nonEwaldDifference = computeNonEwaldEnergyDifference(system, molecule, oldMolecule);
    if (!nonEwaldDifference.has_value())
    {
      restoreLambdaRescale(trial);
      return std::nullopt;
    }
    trial.energyDifference += nonEwaldDifference.value();

    std::array<double, maximumNumberOfDUDlambdaGroups> externalGroupCharge = totalGroupCharge;
    const std::array<double, maximumNumberOfDUDlambdaGroups> ownGroupCharge = groupChargeSum(oldMolecule);
    for (std::size_t g = 0; g < maximumNumberOfDUDlambdaGroups; ++g) externalGroupCharge[g] -= ownGroupCharge[g];

    // chain the Ewald state: the first molecule starts from storedEik, later ones from trialEik
    trial.energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
        (i == 0) ? system.storedEik : system.trialEik, system.trialEik, system.forceField, system.simulationBox,
        molecule, oldMolecule, runningNetCharge, externalGroupCharge);
    runningNetCharge += scaledChargeDifference(molecule, oldMolecule);

    trial.energyDifference +=
        computeTailEnergyDifference(system, tailEffectiveCounts, tailGroupCounts, molecule, oldMolecule);
  }

  // Polarization: rescaling changes the coupling of the fractional molecules to the field they feel
  // (positions unchanged) and the field they produce on the other molecules.
  if (system.forceField.computePolarization)
  {
    for (std::size_t i = 0; i < trial.molecules.size(); ++i)
    {
      const auto &[componentId, moleculeId] = moleculeIds[i];
      std::span<double3> storedField = system.spanElectricFieldOld(componentId, moleculeId);
      trial.energyDifference += Interactions::computePolarizationEnergyDifference(
          system.forceField, storedField, storedField, trial.molecules[i], trial.oldAtoms[i]);
    }

    if (!system.forceField.omitInterPolarization)
    {
      trial.electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
      for (std::size_t i = 0; i < trial.molecules.size(); ++i)
      {
        std::vector<double3> fieldNew(trial.molecules[i].size());
        std::vector<double3> fieldOld(trial.oldAtoms[i].size());
        [[maybe_unused]] std::optional<RunningEnergy> e =
            Interactions::computeInterMolecularPolarizationElectricFieldDifference(
                system.forceField, system.simulationBox, trial.electricFieldNeighborDelta, fieldNew, fieldOld,
                system.spanOfMoleculeAtoms(), trial.molecules[i], trial.oldAtoms[i]);
      }
      trial.energyDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), trial.electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  return trial;
}

static void acceptLambdaRescale(System &system, const LambdaRescaleTrial &trial)
{
  Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

  if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
  {
    std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
    for (std::size_t i = 0; i < storedElectricField.size(); ++i)
    {
      storedElectricField[i] += trial.electricFieldNeighborDelta[i];
    }
  }

  system.computeTailCorrectionCounts();
}

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::LambdaExchange(RandomNumber &random, System &systemA,
                                                                                System &systemB)
{
  const Move::Types move = Move::Types::ParallelTempering;

  const std::optional<std::size_t> tiComponentA = fixedLambdaComponent(systemA);
  const std::optional<std::size_t> tiComponentB = fixedLambdaComponent(systemB);
  if (!tiComponentA.has_value() || !tiComponentB.has_value())
  {
    return std::nullopt;
  }

  Component &componentA = systemA.components[tiComponentA.value()];
  Component &componentB = systemB.components[tiComponentB.value()];

  PropertyLambdaProbabilityHistogram &lambdaA = componentA.fixedLambdaHistogram();
  PropertyLambdaProbabilityHistogram &lambdaB = componentB.fixedLambdaHistogram();

  const std::size_t binA = componentA.fixedLambdaBin.value();
  const std::size_t binB = componentB.fixedLambdaBin.value();
  if (binA == binB)
  {
    return std::nullopt;
  }
  const double lambdaValueA = static_cast<double>(binA) * lambdaA.delta;
  const double lambdaValueB = static_cast<double>(binB) * lambdaB.delta;

  systemA.mc_moves_statistics.addTrial(move);
  systemB.mc_moves_statistics.addTrial(move);

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  // rescale A: lambda_A -> lambda_B
  std::optional<LambdaRescaleTrial> trialA = rescaleToLambda(systemA, tiComponentA.value(), lambdaValueB);
  if (!trialA.has_value())
  {
    return std::nullopt;
  }

  // rescale B: lambda_B -> lambda_A
  std::optional<LambdaRescaleTrial> trialB = rescaleToLambda(systemB, tiComponentB.value(), lambdaValueA);
  if (!trialB.has_value())
  {
    restoreLambdaRescale(trialA.value());
    return std::nullopt;
  }

  systemA.mc_moves_statistics.addConstructed(move);
  systemB.mc_moves_statistics.addConstructed(move);

  const double exponent = -systemA.beta * trialA->energyDifference.potentialEnergy() -
                          systemB.beta * trialB->energyDifference.potentialEnergy();

  if (random.uniform() < std::exp(exponent))
  {
    acceptLambdaRescale(systemA, trialA.value());
    acceptLambdaRescale(systemB, trialB.value());

    // swap the pinned lambda-bins; dU/dlambda samples now accumulate in the exchanged bins
    componentA.fixedLambdaBin = binB;
    lambdaA.setCurrentBin(binB);
    componentB.fixedLambdaBin = binA;
    lambdaB.setCurrentBin(binA);

    systemA.runningEnergies += trialA->energyDifference;
    systemB.runningEnergies += trialB->energyDifference;

    systemA.mc_moves_statistics.addAccepted(move);
    systemB.mc_moves_statistics.addAccepted(move);

    std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
    systemA.mc_moves_cputime[move][Move::Timing::Total] += (time_end - time_begin);
    systemB.mc_moves_cputime[move][Move::Timing::Total] += (time_end - time_begin);

    return std::make_pair(trialA->energyDifference, trialB->energyDifference);
  }

  restoreLambdaRescale(trialA.value());
  restoreLambdaRescale(trialB.value());

  std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
  systemA.mc_moves_cputime[move][Move::Timing::Total] += (time_end - time_begin);
  systemB.mc_moves_cputime[move][Move::Timing::Total] += (time_end - time_begin);

  return std::nullopt;
}
