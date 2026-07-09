module;

module mc_moves_pair_swap_cfcmc_cbmc;

import std;

import component;
import molecule;
import atom;
import double3;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import randomnumbers;
import system;
import property_lambda_probability_histogram;
import running_energy;
import forcefield;
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;
import scaling;

static void accumulateInterMolecularFieldDelta(System& system, std::span<const Atom> newAtoms,
                                               std::span<const Atom> oldAtoms,
                                               std::span<double3> electricFieldNeighborDelta)
{
  std::vector<double3> fieldNew(newAtoms.size());
  std::vector<double3> fieldOld(oldAtoms.size());
  [[maybe_unused]] std::optional<RunningEnergy> energy =
      Interactions::computeInterMolecularPolarizationElectricFieldDifference(
          system.forceField, system.simulationBox, electricFieldNeighborDelta, fieldNew, fieldOld,
          system.spanOfMoleculeAtoms(), newAtoms, oldAtoms);
}

// The pair CB/CFCMC move handles two molecules (one of each linked component) that share a single
// lambda. Energy differences are computed sequentially, one molecule at a time: the first molecule is
// updated in place before the difference of the second molecule is computed, so that the
// cross-interaction between the two molecules of the pair is counted exactly once. The CBMC growth of
// the pair grows the molecule of component A first and the molecule of component B afterwards (with
// the trial molecule of A included in the background); the retrace mirrors this order.

// difference in scaled net charge when 'oldAtoms' is replaced by 'newAtoms'
static double scaledChargeDifference(std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms)
{
  double charge = 0.0;
  for (const Atom& atom : newAtoms) charge += atom.scalingCoulomb * atom.charge;
  for (const Atom& atom : oldAtoms) charge -= atom.scalingCoulomb * atom.charge;
  return charge;
}

// per-group summed (unscaled) charge of the dU/dlambda group-tagged atoms; used as the 'external'
// charge derivative for the net-charge correction when the partner molecule of the pair is updated
static std::array<double, maximumNumberOfDUDlambdaGroups> groupChargeSum(std::span<const Atom> atoms)
{
  std::array<double, maximumNumberOfDUDlambdaGroups> charge{};
  for (const Atom& atom : atoms)
  {
    if (atom.groupId != 0) charge[atom.groupId - 1] += atom.charge;
  }
  return charge;
}

// external-field + framework + intermolecular energy difference for replacing 'oldAtoms' by 'newAtoms'
static std::optional<RunningEnergy> computeNonEwaldEnergyDifference(System& system, std::span<const Atom> newAtoms,
                                                                    std::span<const Atom> oldAtoms,
                                                                    std::span<const Atom> moleculeAtomData)
{
  std::optional<RunningEnergy> externalFieldDifference = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
      newAtoms, oldAtoms);
  if (!externalFieldDifference.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
  if (!frameworkDifference.has_value()) return std::nullopt;

  std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, moleculeAtomData, newAtoms, oldAtoms);
  if (!moleculeDifference.has_value()) return std::nullopt;

  return externalFieldDifference.value() + frameworkDifference.value() + moleculeDifference.value();
}

static RunningEnergy computeTailEnergyDifference(System& system, std::span<const Atom> newAtoms,
                                                 std::span<const Atom> oldAtoms,
                                                 std::span<const Atom> moleculeAtomData)
{
  return Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 moleculeAtomData, newAtoms, oldAtoms) +
         Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                    system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
}

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::pairSwapMove_CFCMC_CBMC(RandomNumber& random,
                                                                                   System& system,
                                                                                   std::size_t selectedComponent)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  const Move::Types move = Move::Types::PairSwapCBCFCMC;
  Component& componentA = system.components[selectedComponent];

  if (!componentA.pairComponentId.has_value())
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
  const std::size_t componentB = componentA.pairComponentId.value();

  // Only the lower-index component performs pair moves to avoid double counting.
  if (componentB >= system.components.size() || selectedComponent >= componentB)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
  Component& componentBRef = system.components[componentB];

  // both fractional molecules of the pair are coupled to the lambda histogram of component A
  PropertyLambdaProbabilityHistogram& lambda = componentA.lambdaPairSwapCB;
  const std::size_t oldBin = lambda.currentBin;
  const double deltaLambda = lambda.delta;
  const double maxChange = componentA.mc_moves_statistics.getMaxChange(move, 2);
  std::make_signed_t<std::size_t> selectedNewBin = lambda.selectNewBin(random, maxChange);

  const std::size_t oldN_A = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  const std::size_t oldN_B = system.numberOfIntegerMoleculesPerComponent[componentB];

  const std::size_t indexFractionalA = system.indexOfFractionalMoleculeForMove(move, selectedComponent);
  const std::size_t indexFractionalB = system.indexOfFractionalMoleculeForMove(move, componentB);

  const double fugacityA = componentA.molFraction * componentA.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double fugacityB = componentBRef.molFraction * componentBRef.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double idealGasA = componentA.idealGasRosenbluthWeight.value_or(1.0);
  const double idealGasB = componentBRef.idealGasRosenbluthWeight.value_or(1.0);

  const double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb = system.forceField.cutOffCoulomb;

  if (selectedNewBin >= std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints))  // Insertion move
  {
    // Steps for insertion Lambda_new = 1 + epsilon
    // ===================================================================
    // (1) Unbiased: the fractional pair with lambda=lambda_old is made integer (lambda=1)
    // (2) Biased: a new fractional pair is grown with CBMC with lambda_new = epsilon

    std::size_t newBin =
        static_cast<std::size_t>(selectedNewBin - std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints));
    double newLambda = deltaLambda * static_cast<double>(newBin);

    componentA.mc_moves_statistics.addTrial(move, 0);

    std::span<Atom> fractionalMoleculeA = system.spanOfMolecule(selectedComponent, indexFractionalA);
    std::span<Atom> fractionalMoleculeB = system.spanOfMolecule(componentB, indexFractionalB);

    std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());

    auto restoreFractionalPair = [&]()
    {
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
    };

    RunningEnergy energyDifference{};
    double runningNetCharge = system.netCharge;

    // (1a) fractional molecule of component A becomes integer
    for (Atom& atom : fractionalMoleculeA) atom.setScalingToInteger();

    if (system.insideBlockedPockets(componentA, fractionalMoleculeA))
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceA =
        computeNonEwaldEnergyDifference(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
    if (!nonEwaldDifferenceA.has_value())
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceA.value();

    time_begin = std::chrono::steady_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMoleculeA, oldFractionalMoleculeA, runningNetCharge,
        groupChargeSum(fractionalMoleculeB));
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);
    runningNetCharge += scaledChargeDifference(fractionalMoleculeA, oldFractionalMoleculeA);

    time_begin = std::chrono::steady_clock::now();
    energyDifference += computeTailEnergyDifference(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);

    // (1b) fractional molecule of component B becomes integer
    for (Atom& atom : fractionalMoleculeB) atom.setScalingToInteger();

    if (system.insideBlockedPockets(componentBRef, fractionalMoleculeB))
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceB =
        computeNonEwaldEnergyDifference(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
    if (!nonEwaldDifferenceB.has_value())
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceB.value();

    time_begin = std::chrono::steady_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB, runningNetCharge,
        groupChargeSum(fractionalMoleculeA));
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);
    runningNetCharge += scaledChargeDifference(fractionalMoleculeB, oldFractionalMoleculeB);

    time_begin = std::chrono::steady_clock::now();
    energyDifference += computeTailEnergyDifference(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);

    // (2a) grow a new fractional molecule of component A with lambda_new
    const std::size_t newMoleculeA = system.numberOfMolecules();
    const std::size_t newMoleculeB = system.numberOfMolecules() + 1;

    time_begin = std::chrono::steady_clock::now();
    std::optional<ChainGrowData> growDataA = CBMC::growMoleculeSwapInsertion(
        random,
        CBMC::GrowContext{system.hasExternalField, system.forceField, system.simulationBox, system.interpolationGrids,
                          system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(),
                          system.spanOfMoleculeAtoms(), system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW,
                          cutOffCoulomb},
        componentA, selectedComponent, componentA.growType, newMoleculeA, newLambda,
        componentA.lambdaPairSwapCB.dUdlambdaGroupId, true);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);

    if (!growDataA)
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    if (system.insideBlockedPockets(componentA, growDataA->atoms))
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // (2b) grow a new fractional molecule of component B with lambda_new; the background includes the
    // trial molecule of A
    std::vector<Atom> moleculeAtomDataWithTrialA(system.spanOfMoleculeAtoms().begin(),
                                                 system.spanOfMoleculeAtoms().end());
    moleculeAtomDataWithTrialA.insert(moleculeAtomDataWithTrialA.end(), growDataA->atoms.begin(),
                                      growDataA->atoms.end());

    time_begin = std::chrono::steady_clock::now();
    std::optional<ChainGrowData> growDataB = CBMC::growMoleculeSwapInsertion(
        random,
        CBMC::GrowContext{system.hasExternalField, system.forceField, system.simulationBox, system.interpolationGrids,
                          system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(),
                          moleculeAtomDataWithTrialA, system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW,
                          cutOffCoulomb},
        componentBRef, componentB, componentBRef.growType, newMoleculeB, newLambda,
        componentA.lambdaPairSwapCB.dUdlambdaGroupId, true);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionNonEwald] += (time_end - time_begin);

    if (!growDataB)
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    if (system.insideBlockedPockets(componentBRef, growDataB->atoms))
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    componentA.mc_moves_statistics.addConstructed(move, 0);

    std::vector<double3> electricFieldNeighborDelta;
    std::vector<double3> grownElectricFieldA(growDataA->atoms.size());
    std::vector<double3> grownElectricFieldB(growDataB->atoms.size());
    RunningEnergy polarizationDifference;
    if (system.forceField.computePolarization)
    {
      if (!system.forceField.omitInterPolarization)
      {
        electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
        accumulateInterMolecularFieldDelta(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                           electricFieldNeighborDelta);
        accumulateInterMolecularFieldDelta(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                           electricFieldNeighborDelta);

        Interactions::computeFrameworkMoleculeElectricFieldDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), grownElectricFieldA, {},
            growDataA->atoms, {});
        Interactions::computeEwaldFourierElectricFieldDifference(
            system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
            system.totalEik, system.forceField, system.simulationBox, grownElectricFieldA, {}, growDataA->atoms, {});

        Interactions::computeFrameworkMoleculeElectricFieldDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), grownElectricFieldB, {},
            growDataB->atoms, {});
        Interactions::computeEwaldFourierElectricFieldDifference(
            system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
            system.totalEik, system.forceField, system.simulationBox, grownElectricFieldB, {}, growDataB->atoms, {});

        accumulateInterMolecularFieldDelta(system, growDataA->atoms, {}, electricFieldNeighborDelta);
        accumulateInterMolecularFieldDelta(system, growDataB->atoms, {}, electricFieldNeighborDelta);

        std::vector<double3> fieldOnAfromB(growDataA->atoms.size());
        std::vector<double3> fieldOnBfromA(growDataB->atoms.size());
        [[maybe_unused]] std::optional<RunningEnergy> mutualEnergy =
            Interactions::computeInterMolecularPolarizationElectricFieldDifference(
                system.forceField, system.simulationBox, fieldOnAfromB, fieldOnBfromA, std::span<double3>{},
                growDataA->atoms, growDataB->atoms, {});
        for (std::size_t i = 0; i < grownElectricFieldA.size(); ++i) grownElectricFieldA[i] += fieldOnAfromB[i];
        for (std::size_t i = 0; i < grownElectricFieldB.size(); ++i) grownElectricFieldB[i] += fieldOnBfromA[i];

        polarizationDifference =
            Interactions::computePolarizationEnergyDifference(system.forceField, grownElectricFieldA, {},
                                                              growDataA->atoms, {}) +
            Interactions::computePolarizationEnergyDifference(system.forceField, grownElectricFieldB, {},
                                                              growDataB->atoms, {}) +
            Interactions::computePolarizationEnergyNeighborDifference(
                system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
                system.spanOfMoleculeAtoms());
      }
      else
      {
        Interactions::computeFrameworkMoleculeElectricFieldDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), grownElectricFieldA, {},
            growDataA->atoms, {});
        Interactions::computeEwaldFourierElectricFieldDifference(
            system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
            system.totalEik, system.forceField, system.simulationBox, grownElectricFieldA, {}, growDataA->atoms, {});

        Interactions::computeFrameworkMoleculeElectricFieldDifference(
            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), grownElectricFieldB, {},
            growDataB->atoms, {});
        Interactions::computeEwaldFourierElectricFieldDifference(
            system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
            system.totalEik, system.forceField, system.simulationBox, grownElectricFieldB, {}, growDataB->atoms, {});

        polarizationDifference =
            Interactions::computePolarizationEnergyDifference(system.forceField, grownElectricFieldA, {},
                                                              growDataA->atoms, {}) +
            Interactions::computePolarizationEnergyDifference(system.forceField, grownElectricFieldB, {},
                                                              growDataB->atoms, {});
      }
    }

    // Ewald and tail-correction energy contributions of the grown pair
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy energyFourierDifferenceA = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, std::span(growDataA->atoms.begin(), growDataA->atoms.end()), {}, runningNetCharge);
    runningNetCharge += scaledChargeDifference(growDataA->atoms, {});
    RunningEnergy energyFourierDifferenceB = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, std::span(growDataB->atoms.begin(), growDataB->atoms.end()), {}, runningNetCharge,
        groupChargeSum(growDataA->atoms));
    RunningEnergy energyFourierDifference = energyFourierDifferenceA + energyFourierDifferenceB;
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionEwald] += (time_end - time_begin);

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy tailEnergyDifferenceGrow =
        computeTailEnergyDifference(system, growDataA->atoms, {}, system.spanOfMoleculeAtoms()) +
        computeTailEnergyDifference(system, growDataB->atoms, {}, moleculeAtomDataWithTrialA);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::InsertionTail] += (time_end - time_begin);

    const double correctionFactorEwald = std::exp(
        -system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifferenceGrow.potentialEnergy() +
                        polarizationDifference.potentialEnergy()));

    // Calculate acceptance probability
    const double preFactor =
        correctionFactorEwald * (system.beta * fugacityA * system.simulationBox.volume / static_cast<double>(1 + oldN_A)) *
        (system.beta * fugacityB * system.simulationBox.volume / static_cast<double>(1 + oldN_B));
    const double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    const double Pacc = preFactor * (growDataA->RosenbluthWeight / idealGasA) *
                        (growDataB->RosenbluthWeight / idealGasB) *
                        std::exp(-system.beta * energyDifference.potentialEnergy() + biasTerm);

    const double biasTransitionMatrix = system.tmmc.biasFactor(oldN_A + 1, oldN_A);

    if (system.tmmc.doTMMC)
    {
      if (oldN_A + 1 > system.tmmc.maxMacrostate)
      {
        restoreFractionalPair();
        return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
      }
    }

    if (random.uniform() < biasTransitionMatrix * Pacc)
    {
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

      lambda.setCurrentBin(newBin);

      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      if (system.forceField.computePolarization)
      {
        system.insertMoleculePolarization(selectedComponent, growDataA->molecule, growDataA->atoms,
                                          grownElectricFieldA);
      }
      else
      {
        system.insertMolecule(selectedComponent, growDataA->molecule, growDataA->atoms);
      }
      std::size_t lastMoleculeIdA = system.numberOfMoleculesPerComponent[selectedComponent] - 1;
      std::span<Atom> lastMoleculeA = system.spanOfMolecule(selectedComponent, lastMoleculeIdA);
      fractionalMoleculeA = system.spanOfMolecule(selectedComponent, indexFractionalA);
      std::swap_ranges(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), lastMoleculeA.begin());
      if (system.forceField.computePolarization)
      {
        std::span<double3> fieldFractionalA = system.spanElectricFieldOld(selectedComponent, indexFractionalA);
        std::span<double3> fieldLastA = system.spanElectricFieldOld(selectedComponent, lastMoleculeIdA);
        std::swap_ranges(fieldFractionalA.begin(), fieldFractionalA.end(), fieldLastA.begin());
      }
      std::swap(system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, indexFractionalA)],
                system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, lastMoleculeIdA)]);

      if (system.forceField.computePolarization)
      {
        system.insertMoleculePolarization(componentB, growDataB->molecule, growDataB->atoms, grownElectricFieldB);
      }
      else
      {
        system.insertMolecule(componentB, growDataB->molecule, growDataB->atoms);
      }
      std::size_t lastMoleculeIdB = system.numberOfMoleculesPerComponent[componentB] - 1;
      std::span<Atom> lastMoleculeB = system.spanOfMolecule(componentB, lastMoleculeIdB);
      fractionalMoleculeB = system.spanOfMolecule(componentB, indexFractionalB);
      std::swap_ranges(fractionalMoleculeB.begin(), fractionalMoleculeB.end(), lastMoleculeB.begin());
      if (system.forceField.computePolarization)
      {
        std::span<double3> fieldFractionalB = system.spanElectricFieldOld(componentB, indexFractionalB);
        std::span<double3> fieldLastB = system.spanElectricFieldOld(componentB, lastMoleculeIdB);
        std::swap_ranges(fieldFractionalB.begin(), fieldFractionalB.end(), fieldLastB.begin());
      }
      std::swap(system.moleculeData[system.moleculeIndexOfComponent(componentB, indexFractionalB)],
                system.moleculeData[system.moleculeIndexOfComponent(componentB, lastMoleculeIdB)]);

      system.updateMoleculeAtomInformation();

      componentA.mc_moves_statistics.addAccepted(move, 0);

      return {energyDifference + growDataA->energies + growDataB->energies + energyFourierDifference +
                  tailEnergyDifferenceGrow + polarizationDifference,
              double3(0.0, 1.0 - Pacc, Pacc)};
    }

    restoreFractionalPair();
    return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
  }
  else if (selectedNewBin < 0)  // Deletion move
  {
    // Steps for deletion Lambda_new = -epsilon
    // ===================================================================
    // (1) Biased: the existing fractional pair is retraced using CBMC with lambda=lambda_old and removed
    // (2) Unbiased: a randomly selected integer molecule of each component becomes the new
    //     fractional pair with lambda_new = 1 - epsilon

    componentA.mc_moves_statistics.addTrial(move, 1);

    if (oldN_A == 0 || oldN_B == 0)
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    std::size_t newBin =
        static_cast<std::size_t>(selectedNewBin + std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints));
    double newLambda = deltaLambda * static_cast<double>(newBin);

    std::size_t selectedMoleculeA = system.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::size_t selectedMoleculeB = system.randomIntegerMoleculeOfComponent(random, componentB);

    std::span<Atom> fractionalMoleculeA = system.spanOfMolecule(selectedComponent, indexFractionalA);
    std::span<Atom> fractionalMoleculeB = system.spanOfMolecule(componentB, indexFractionalB);
    std::span<Atom> newFractionalMoleculeA = system.spanOfMolecule(selectedComponent, selectedMoleculeA);
    std::span<Atom> newFractionalMoleculeB = system.spanOfMolecule(componentB, selectedMoleculeB);

    std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
    std::vector<Atom> oldNewFractionalMoleculeA(newFractionalMoleculeA.begin(), newFractionalMoleculeA.end());
    std::vector<Atom> oldNewFractionalMoleculeB(newFractionalMoleculeB.begin(), newFractionalMoleculeB.end());

    auto restoreMolecules = [&]()
    {
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldNewFractionalMoleculeA.begin(), oldNewFractionalMoleculeA.end(), newFractionalMoleculeA.begin());
      std::copy(oldNewFractionalMoleculeB.begin(), oldNewFractionalMoleculeB.end(), newFractionalMoleculeB.begin());
    };

    // (1) retrace the fractional pair; the growth order of the reverse (insertion) move is A first,
    // then B with A present, so B is retraced with A in the background and A without B
    const std::uint32_t moleculeIdFractionalB = fractionalMoleculeB.front().moleculeId;
    std::vector<Atom> backgroundWithoutFractionalB;
    backgroundWithoutFractionalB.reserve(system.spanOfMoleculeAtoms().size());
    for (const Atom& atom : system.spanOfMoleculeAtoms())
    {
      if (atom.moleculeId == moleculeIdFractionalB) continue;
      backgroundWithoutFractionalB.push_back(atom);
    }

    time_begin = std::chrono::steady_clock::now();
    ChainRetraceData retraceDataB = CBMC::retraceMoleculeSwapDeletion(
        random,
        CBMC::GrowContext{system.hasExternalField, system.forceField, system.simulationBox, system.interpolationGrids,
                          system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(),
                          system.spanOfMoleculeAtoms(), system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW,
                          cutOffCoulomb},
        componentBRef, componentBRef.growType, fractionalMoleculeB);
    ChainRetraceData retraceDataA = CBMC::retraceMoleculeSwapDeletion(
        random,
        CBMC::GrowContext{system.hasExternalField, system.forceField, system.simulationBox, system.interpolationGrids,
                          system.externalFieldInterpolationGrid, system.framework, system.spanOfFrameworkAtoms(),
                          backgroundWithoutFractionalB, system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW,
                          cutOffCoulomb},
        componentA, componentA.growType, fractionalMoleculeA);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);

    double runningNetCharge = system.netCharge;

    // Ewald and tail-correction energy contributions of removing the fractional pair
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy energyFourierDifferenceA = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, {}, fractionalMoleculeA, runningNetCharge, groupChargeSum(fractionalMoleculeB));
    runningNetCharge += scaledChargeDifference({}, fractionalMoleculeA);
    RunningEnergy energyFourierDifferenceB = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, {}, fractionalMoleculeB, runningNetCharge);
    runningNetCharge += scaledChargeDifference({}, fractionalMoleculeB);
    RunningEnergy energyFourierDifference = energyFourierDifferenceA + energyFourierDifferenceB;
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);

    time_begin = std::chrono::steady_clock::now();
    RunningEnergy tailEnergyDifferenceRetrace =
        computeTailEnergyDifference(system, {}, fractionalMoleculeB, system.spanOfMoleculeAtoms()) +
        computeTailEnergyDifference(system, {}, fractionalMoleculeA, backgroundWithoutFractionalB);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);

    // deactivate the fractional pair
    for (Atom& atom : fractionalMoleculeA)
    {
      atom.setScalingOff();
    }
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.setScalingOff();
    }

    RunningEnergy energyDifference{};

    // (2a) the selected integer molecule of component A becomes fractional with lambda_new
    const std::uint8_t groupIdA = componentA.lambdaPairSwapCB.dUdlambdaGroupId;
    for (Atom& atom : newFractionalMoleculeA)
    {
      atom.setScalingToFractional(newLambda, groupIdA);
    }

    if (system.insideBlockedPockets(componentA, newFractionalMoleculeA))
    {
      restoreMolecules();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceNewA =
        computeNonEwaldEnergyDifference(system, newFractionalMoleculeA, oldNewFractionalMoleculeA,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
    if (!nonEwaldDifferenceNewA.has_value())
    {
      restoreMolecules();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceNewA.value();

    time_begin = std::chrono::steady_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, newFractionalMoleculeA, oldNewFractionalMoleculeA, runningNetCharge);
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);
    runningNetCharge += scaledChargeDifference(newFractionalMoleculeA, oldNewFractionalMoleculeA);

    time_begin = std::chrono::steady_clock::now();
    energyDifference += computeTailEnergyDifference(system, newFractionalMoleculeA, oldNewFractionalMoleculeA,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);

    // (2b) the selected integer molecule of component B becomes fractional with lambda_new
    const std::uint8_t groupIdB = componentA.lambdaPairSwapCB.dUdlambdaGroupId;
    for (Atom& atom : newFractionalMoleculeB)
    {
      atom.setScalingToFractional(newLambda, groupIdB);
    }

    if (system.insideBlockedPockets(componentBRef, newFractionalMoleculeB))
    {
      restoreMolecules();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceNewB =
        computeNonEwaldEnergyDifference(system, newFractionalMoleculeB, oldNewFractionalMoleculeB,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::DeletionNonEwald] += (time_end - time_begin);
    if (!nonEwaldDifferenceNewB.has_value())
    {
      restoreMolecules();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceNewB.value();

    time_begin = std::chrono::steady_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, newFractionalMoleculeB, oldNewFractionalMoleculeB, runningNetCharge,
        groupChargeSum(newFractionalMoleculeA));
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::DeletionEwald] += (time_end - time_begin);

    time_begin = std::chrono::steady_clock::now();
    energyDifference += computeTailEnergyDifference(system, newFractionalMoleculeB, oldNewFractionalMoleculeB,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::DeletionTail] += (time_end - time_begin);

    componentA.mc_moves_statistics.addConstructed(move, 1);

    std::vector<double3> electricFieldNeighborDelta;
    RunningEnergy polarizationDifference;
    if (system.forceField.computePolarization)
    {
      std::span<double3> storedFieldFractionalA =
          system.spanElectricFieldOld(selectedComponent, indexFractionalA);
      std::span<double3> storedFieldFractionalB = system.spanElectricFieldOld(componentB, indexFractionalB);
      polarizationDifference =
          Interactions::computePolarizationEnergyDifference(system.forceField, {}, storedFieldFractionalA, {},
                                                            fractionalMoleculeA) +
          Interactions::computePolarizationEnergyDifference(system.forceField, {}, storedFieldFractionalB, {},
                                                            fractionalMoleculeB);

      if (!system.forceField.omitInterPolarization)
      {
        electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
        accumulateInterMolecularFieldDelta(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                           electricFieldNeighborDelta);
        accumulateInterMolecularFieldDelta(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                           electricFieldNeighborDelta);
        accumulateInterMolecularFieldDelta(system, newFractionalMoleculeA, oldNewFractionalMoleculeA,
                                           electricFieldNeighborDelta);
        accumulateInterMolecularFieldDelta(system, newFractionalMoleculeB, oldNewFractionalMoleculeB,
                                           electricFieldNeighborDelta);

        std::span<Atom> allAtoms = system.spanOfMoleculeAtoms();
        std::size_t offsetFractionalA =
            static_cast<std::size_t>(fractionalMoleculeA.data() - allAtoms.data());
        for (std::size_t k = 0; k < fractionalMoleculeA.size(); ++k)
        {
          electricFieldNeighborDelta[offsetFractionalA + k] = double3(0.0, 0.0, 0.0);
        }
        std::size_t offsetFractionalB =
            static_cast<std::size_t>(fractionalMoleculeB.data() - allAtoms.data());
        for (std::size_t k = 0; k < fractionalMoleculeB.size(); ++k)
        {
          electricFieldNeighborDelta[offsetFractionalB + k] = double3(0.0, 0.0, 0.0);
        }

        polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
            system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
            system.spanOfMoleculeAtoms());
      }
    }

    const double correctionFactorEwald = std::exp(
        -system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifferenceRetrace.potentialEnergy() +
                        polarizationDifference.potentialEnergy()));

    // Calculate acceptance probability
    const double preFactor =
        correctionFactorEwald * (static_cast<double>(oldN_A) / (system.beta * fugacityA * system.simulationBox.volume)) *
        (static_cast<double>(oldN_B) / (system.beta * fugacityB * system.simulationBox.volume));
    const double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    const double Pacc = preFactor * (idealGasA / retraceDataA.RosenbluthWeight) *
                        (idealGasB / retraceDataB.RosenbluthWeight) *
                        std::exp(-system.beta * energyDifference.potentialEnergy() + biasTerm);

    const double biasTransitionMatrix = system.tmmc.biasFactor(oldN_A - 1, oldN_A);

    if (system.tmmc.doTMMC)
    {
      if (oldN_A - 1 < system.tmmc.minMacrostate)
      {
        restoreMolecules();
        return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
      }
    }

    if (random.uniform() < biasTransitionMatrix * Pacc)
    {
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      lambda.setCurrentBin(newBin);

      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      // component A: swap the new fractional molecule into the fractional slot and delete the old one
      std::swap_ranges(newFractionalMoleculeA.begin(), newFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      if (system.forceField.computePolarization)
      {
        std::span<double3> storedFieldFractionalA =
            system.spanElectricFieldOld(selectedComponent, indexFractionalA);
        std::span<double3> storedFieldSelectedA =
            system.spanElectricFieldOld(selectedComponent, selectedMoleculeA);
        std::swap_ranges(storedFieldFractionalA.begin(), storedFieldFractionalA.end(),
                         storedFieldSelectedA.begin());
      }
      std::swap(system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMoleculeA)],
                system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, indexFractionalA)]);
      system.deleteMolecule(selectedComponent, selectedMoleculeA, newFractionalMoleculeA);

      // component B: re-acquire the spans (deletion invalidates them) and do the same
      fractionalMoleculeB = system.spanOfMolecule(componentB, indexFractionalB);
      newFractionalMoleculeB = system.spanOfMolecule(componentB, selectedMoleculeB);
      std::swap_ranges(newFractionalMoleculeB.begin(), newFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      if (system.forceField.computePolarization)
      {
        std::span<double3> storedFieldFractionalB = system.spanElectricFieldOld(componentB, indexFractionalB);
        std::span<double3> storedFieldSelectedB =
            system.spanElectricFieldOld(componentB, selectedMoleculeB);
        std::swap_ranges(storedFieldFractionalB.begin(), storedFieldFractionalB.end(),
                         storedFieldSelectedB.begin());
      }
      std::swap(system.moleculeData[system.moleculeIndexOfComponent(componentB, selectedMoleculeB)],
                system.moleculeData[system.moleculeIndexOfComponent(componentB, indexFractionalB)]);
      system.deleteMolecule(componentB, selectedMoleculeB, newFractionalMoleculeB);

      componentA.mc_moves_statistics.addAccepted(move, 1);

      return {energyDifference + energyFourierDifference + tailEnergyDifferenceRetrace - retraceDataA.energies -
                  retraceDataB.energies + polarizationDifference,
              double3(Pacc, 1.0 - Pacc, 0.0)};
    }

    restoreMolecules();
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }
  else  // Lambda-change move
  {
    std::size_t newBin = static_cast<std::size_t>(selectedNewBin);
    double newLambda = deltaLambda * static_cast<double>(newBin);

    componentA.mc_moves_statistics.addTrial(move, 2);

    std::span<Atom> fractionalMoleculeA = system.spanOfMolecule(selectedComponent, indexFractionalA);
    std::span<Atom> fractionalMoleculeB = system.spanOfMolecule(componentB, indexFractionalB);

    std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());

    auto restoreFractionalPair = [&]()
    {
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
    };

    RunningEnergy energyDifference{};
    double runningNetCharge = system.netCharge;

    // rescale the fractional molecule of component A
    for (Atom& atom : fractionalMoleculeA) atom.setScaling(newLambda);

    if (system.insideBlockedPockets(componentA, fractionalMoleculeA))
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceA =
        computeNonEwaldEnergyDifference(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaNonEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::LambdaNonEwald] += (time_end - time_begin);
    if (!nonEwaldDifferenceA.has_value())
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceA.value();

    time_begin = std::chrono::steady_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMoleculeA, oldFractionalMoleculeA, runningNetCharge,
        groupChargeSum(fractionalMoleculeB));
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::LambdaEwald] += (time_end - time_begin);
    runningNetCharge += scaledChargeDifference(fractionalMoleculeA, oldFractionalMoleculeA);

    time_begin = std::chrono::steady_clock::now();
    energyDifference += computeTailEnergyDifference(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::LambdaTail] += (time_end - time_begin);

    // rescale the fractional molecule of component B
    for (Atom& atom : fractionalMoleculeB) atom.setScaling(newLambda);

    if (system.insideBlockedPockets(componentBRef, fractionalMoleculeB))
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    time_begin = std::chrono::steady_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceB =
        computeNonEwaldEnergyDifference(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaNonEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::LambdaNonEwald] += (time_end - time_begin);
    if (!nonEwaldDifferenceB.has_value())
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceB.value();

    time_begin = std::chrono::steady_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB, runningNetCharge,
        groupChargeSum(fractionalMoleculeA));
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::LambdaEwald] += (time_end - time_begin);

    time_begin = std::chrono::steady_clock::now();
    energyDifference += computeTailEnergyDifference(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    componentA.mc_moves_cputime[move][Move::Timing::LambdaTail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::LambdaTail] += (time_end - time_begin);

    componentA.mc_moves_statistics.addConstructed(move, 2);

    std::vector<double3> electricFieldNeighborDelta;
    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
      accumulateInterMolecularFieldDelta(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                         electricFieldNeighborDelta);
      accumulateInterMolecularFieldDelta(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                         electricFieldNeighborDelta);
      energyDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }

    const double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];

    if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + biasTerm))
    {
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      componentA.mc_moves_statistics.addAccepted(move, 2);

      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      lambda.setCurrentBin(newBin);

      return {energyDifference, double3(0.0, 1.0, 0.0)};
    }

    restoreFractionalPair();
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
}
