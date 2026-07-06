module;

module mc_moves_pair_swap_cfcmc;

import std;

import component;
import molecule;
import atom;
import double3;
import simulationbox;
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
import mc_moves_move_types;
import scaling;

// The pair CFCMC move handles two molecules (one of each linked component) that share a single lambda.
// Energy differences are computed sequentially, one molecule at a time: the first molecule is updated
// in place before the difference of the second molecule is computed, so that the cross-interaction
// between the two molecules of the pair is counted exactly once.

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

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::pairSwapMove_CFCMC(RandomNumber& random, System& system,
                                                                              std::size_t selectedComponent)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  const Move::Types move = Move::Types::PairSwapCFCMC;
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
  PropertyLambdaProbabilityHistogram& lambda = componentA.lambdaPairSwap;
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

  if (selectedNewBin >= std::make_signed_t<std::size_t>(lambda.numberOfSamplePoints))  // Insertion move
  {
    // Steps for insertion Lambda_new = 1 + epsilon
    // ===================================================================
    // (1) Unbiased: the fractional pair with lambda=lambda_old is made integer (lambda=1)
    // (2) Unbiased: a new fractional pair is inserted at random positions with lambda_new = epsilon

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

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceA =
        computeNonEwaldEnergyDifference(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);
    if (!nonEwaldDifferenceA.has_value())
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceA.value();

    time_begin = std::chrono::system_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMoleculeA, oldFractionalMoleculeA, runningNetCharge,
        groupChargeSum(fractionalMoleculeB));
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);
    runningNetCharge += scaledChargeDifference(fractionalMoleculeA, oldFractionalMoleculeA);

    time_begin = std::chrono::system_clock::now();
    energyDifference += computeTailEnergyDifference(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);

    // (1b) fractional molecule of component B becomes integer
    for (Atom& atom : fractionalMoleculeB) atom.setScalingToInteger();

    if (system.insideBlockedPockets(componentBRef, fractionalMoleculeB))
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceB =
        computeNonEwaldEnergyDifference(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);
    if (!nonEwaldDifferenceB.has_value())
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceB.value();

    time_begin = std::chrono::system_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB, runningNetCharge,
        groupChargeSum(fractionalMoleculeA));
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);
    runningNetCharge += scaledChargeDifference(fractionalMoleculeB, oldFractionalMoleculeB);

    time_begin = std::chrono::system_clock::now();
    energyDifference += computeTailEnergyDifference(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);

    // (2) insert a new fractional pair with lambda_new at random positions
    std::pair<Molecule, std::vector<Atom>> trialMoleculeA =
        componentA.equilibratedMoleculeRandomInBox(random, selectedComponent, system.simulationBox);
    std::pair<Molecule, std::vector<Atom>> trialMoleculeB =
        componentBRef.equilibratedMoleculeRandomInBox(random, componentB, system.simulationBox);

    const std::size_t upcomingMoleculeIdA = system.numberOfMolecules();
    const std::size_t upcomingMoleculeIdB = system.numberOfMolecules() + 1;

    const std::uint8_t groupIdA = componentA.lambdaPairSwap.dUdlambdaGroupId;
    for (Atom& atom : trialMoleculeA.second)
    {
      atom.moleculeId = static_cast<std::uint32_t>(upcomingMoleculeIdA);
      atom.componentId = static_cast<std::uint8_t>(selectedComponent);
      atom.setScalingToFractional(newLambda, groupIdA);
    }
    const std::uint8_t groupIdB = componentA.lambdaPairSwap.dUdlambdaGroupId;
    for (Atom& atom : trialMoleculeB.second)
    {
      atom.moleculeId = static_cast<std::uint32_t>(upcomingMoleculeIdB);
      atom.componentId = static_cast<std::uint8_t>(componentB);
      atom.setScalingToFractional(newLambda, groupIdB);
    }

    if (system.insideBlockedPockets(componentA, trialMoleculeA.second) ||
        system.insideBlockedPockets(componentBRef, trialMoleculeB.second))
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // (2a) new fractional molecule of component A
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceTrialA =
        computeNonEwaldEnergyDifference(system, trialMoleculeA.second, {}, system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);
    if (!nonEwaldDifferenceTrialA.has_value())
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceTrialA.value();

    time_begin = std::chrono::system_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, trialMoleculeA.second, {}, runningNetCharge);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);
    runningNetCharge += scaledChargeDifference(trialMoleculeA.second, {});

    time_begin = std::chrono::system_clock::now();
    energyDifference += computeTailEnergyDifference(system, trialMoleculeA.second, {}, system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);

    // (2b) new fractional molecule of component B; the background includes the trial molecule of A
    std::vector<Atom> moleculeAtomDataWithTrialA(system.spanOfMoleculeAtoms().begin(),
                                                 system.spanOfMoleculeAtoms().end());
    moleculeAtomDataWithTrialA.insert(moleculeAtomDataWithTrialA.end(), trialMoleculeA.second.begin(),
                                      trialMoleculeA.second.end());

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceTrialB =
        computeNonEwaldEnergyDifference(system, trialMoleculeB.second, {}, moleculeAtomDataWithTrialA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);
    if (!nonEwaldDifferenceTrialB.has_value())
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceTrialB.value();

    time_begin = std::chrono::system_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, trialMoleculeB.second, {}, runningNetCharge, groupChargeSum(trialMoleculeA.second));
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    energyDifference += computeTailEnergyDifference(system, trialMoleculeB.second, {}, moleculeAtomDataWithTrialA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Insertion-Tail"] += (time_end - time_begin);

    componentA.mc_moves_statistics.addConstructed(move, 0);

    // Calculate acceptance probability
    const double preFactor =
        (system.beta * fugacityA * system.simulationBox.volume / static_cast<double>(1 + oldN_A)) *
        (system.beta * fugacityB * system.simulationBox.volume / static_cast<double>(1 + oldN_B));
    const double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    const double Pacc = preFactor * std::exp(-system.beta * energyDifference.potentialEnergy() + biasTerm);

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

      // insert the new fractional molecule of component A and swap it into the fractional slot
      system.insertMolecule(selectedComponent, trialMoleculeA.first, trialMoleculeA.second);
      std::size_t lastMoleculeIdA = system.numberOfMoleculesPerComponent[selectedComponent] - 1;
      std::span<Atom> lastMoleculeA = system.spanOfMolecule(selectedComponent, lastMoleculeIdA);
      fractionalMoleculeA = system.spanOfMolecule(selectedComponent, indexFractionalA);
      std::swap_ranges(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), lastMoleculeA.begin());
      std::swap(system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, indexFractionalA)],
                system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, lastMoleculeIdA)]);

      // insert the new fractional molecule of component B and swap it into the fractional slot
      system.insertMolecule(componentB, trialMoleculeB.first, trialMoleculeB.second);
      std::size_t lastMoleculeIdB = system.numberOfMoleculesPerComponent[componentB] - 1;
      std::span<Atom> lastMoleculeB = system.spanOfMolecule(componentB, lastMoleculeIdB);
      fractionalMoleculeB = system.spanOfMolecule(componentB, indexFractionalB);
      std::swap_ranges(fractionalMoleculeB.begin(), fractionalMoleculeB.end(), lastMoleculeB.begin());
      std::swap(system.moleculeData[system.moleculeIndexOfComponent(componentB, indexFractionalB)],
                system.moleculeData[system.moleculeIndexOfComponent(componentB, lastMoleculeIdB)]);

      system.updateMoleculeAtomInformation();

      componentA.mc_moves_statistics.addAccepted(move, 0);

      return {energyDifference, double3(0.0, 1.0 - Pacc, Pacc)};
    }

    restoreFractionalPair();
    return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
  }
  else if (selectedNewBin < 0)  // Deletion move
  {
    // Steps for deletion Lambda_new = -epsilon
    // ===================================================================
    // (1) Unbiased: the existing fractional pair with lambda=lambda_old is removed (lambda=0)
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

    RunningEnergy energyDifference{};
    double runningNetCharge = system.netCharge;

    // (1a) fractional molecule of component A is switched off (lambda=0)
    for (Atom& atom : fractionalMoleculeA)
    {
      atom.setScalingOff();
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceA =
        computeNonEwaldEnergyDifference(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-NonEwald"] += (time_end - time_begin);
    if (!nonEwaldDifferenceA.has_value())
    {
      restoreMolecules();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceA.value();

    time_begin = std::chrono::system_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMoleculeA, oldFractionalMoleculeA, runningNetCharge,
        groupChargeSum(fractionalMoleculeB));
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);
    runningNetCharge += scaledChargeDifference(fractionalMoleculeA, oldFractionalMoleculeA);

    time_begin = std::chrono::system_clock::now();
    energyDifference += computeTailEnergyDifference(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);

    // (1b) fractional molecule of component B is switched off (lambda=0)
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.setScalingOff();
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceB =
        computeNonEwaldEnergyDifference(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-NonEwald"] += (time_end - time_begin);
    if (!nonEwaldDifferenceB.has_value())
    {
      restoreMolecules();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceB.value();

    time_begin = std::chrono::system_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB, runningNetCharge);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);
    runningNetCharge += scaledChargeDifference(fractionalMoleculeB, oldFractionalMoleculeB);

    time_begin = std::chrono::system_clock::now();
    energyDifference += computeTailEnergyDifference(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);

    // (2a) the selected integer molecule of component A becomes fractional with lambda_new
    const std::uint8_t groupIdA = componentA.lambdaPairSwap.dUdlambdaGroupId;
    for (Atom& atom : newFractionalMoleculeA)
    {
      atom.setScalingToFractional(newLambda, groupIdA);
    }

    if (system.insideBlockedPockets(componentA, newFractionalMoleculeA))
    {
      restoreMolecules();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceNewA =
        computeNonEwaldEnergyDifference(system, newFractionalMoleculeA, oldNewFractionalMoleculeA,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-NonEwald"] += (time_end - time_begin);
    if (!nonEwaldDifferenceNewA.has_value())
    {
      restoreMolecules();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceNewA.value();

    time_begin = std::chrono::system_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, newFractionalMoleculeA, oldNewFractionalMoleculeA, runningNetCharge,
        groupChargeSum(fractionalMoleculeB));
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);
    runningNetCharge += scaledChargeDifference(newFractionalMoleculeA, oldNewFractionalMoleculeA);

    time_begin = std::chrono::system_clock::now();
    energyDifference += computeTailEnergyDifference(system, newFractionalMoleculeA, oldNewFractionalMoleculeA,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);

    // (2b) the selected integer molecule of component B becomes fractional with lambda_new
    const std::uint8_t groupIdB = componentA.lambdaPairSwap.dUdlambdaGroupId;
    for (Atom& atom : newFractionalMoleculeB)
    {
      atom.setScalingToFractional(newLambda, groupIdB);
    }

    if (system.insideBlockedPockets(componentBRef, newFractionalMoleculeB))
    {
      restoreMolecules();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceNewB =
        computeNonEwaldEnergyDifference(system, newFractionalMoleculeB, oldNewFractionalMoleculeB,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-NonEwald"] += (time_end - time_begin);
    if (!nonEwaldDifferenceNewB.has_value())
    {
      restoreMolecules();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceNewB.value();

    time_begin = std::chrono::system_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, newFractionalMoleculeB, oldNewFractionalMoleculeB, runningNetCharge,
        groupChargeSum(newFractionalMoleculeA));
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    energyDifference += computeTailEnergyDifference(system, newFractionalMoleculeB, oldNewFractionalMoleculeB,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Deletion-Tail"] += (time_end - time_begin);

    componentA.mc_moves_statistics.addConstructed(move, 1);

    // Calculate acceptance probability
    const double preFactor =
        (static_cast<double>(oldN_A) / (system.beta * fugacityA * system.simulationBox.volume)) *
        (static_cast<double>(oldN_B) / (system.beta * fugacityB * system.simulationBox.volume));
    const double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
    const double Pacc = preFactor * std::exp(-system.beta * energyDifference.potentialEnergy() + biasTerm);

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

      // component A: swap the new fractional molecule into the fractional slot and delete the old one
      std::swap_ranges(newFractionalMoleculeA.begin(), newFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::swap(system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMoleculeA)],
                system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, indexFractionalA)]);
      system.deleteMolecule(selectedComponent, selectedMoleculeA, newFractionalMoleculeA);

      // component B: re-acquire the spans (deletion invalidates them) and do the same
      fractionalMoleculeB = system.spanOfMolecule(componentB, indexFractionalB);
      newFractionalMoleculeB = system.spanOfMolecule(componentB, selectedMoleculeB);
      std::swap_ranges(newFractionalMoleculeB.begin(), newFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::swap(system.moleculeData[system.moleculeIndexOfComponent(componentB, selectedMoleculeB)],
                system.moleculeData[system.moleculeIndexOfComponent(componentB, indexFractionalB)]);
      system.deleteMolecule(componentB, selectedMoleculeB, newFractionalMoleculeB);

      componentA.mc_moves_statistics.addAccepted(move, 1);

      return {energyDifference, double3(Pacc, 1.0 - Pacc, 0.0)};
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

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceA =
        computeNonEwaldEnergyDifference(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Lambda-NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Lambda-NonEwald"] += (time_end - time_begin);
    if (!nonEwaldDifferenceA.has_value())
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceA.value();

    time_begin = std::chrono::system_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMoleculeA, oldFractionalMoleculeA, runningNetCharge,
        groupChargeSum(fractionalMoleculeB));
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Lambda-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Lambda-Ewald"] += (time_end - time_begin);
    runningNetCharge += scaledChargeDifference(fractionalMoleculeA, oldFractionalMoleculeA);

    time_begin = std::chrono::system_clock::now();
    energyDifference += computeTailEnergyDifference(system, fractionalMoleculeA, oldFractionalMoleculeA,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Lambda-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Lambda-Tail"] += (time_end - time_begin);

    // rescale the fractional molecule of component B
    for (Atom& atom : fractionalMoleculeB) atom.setScaling(newLambda);

    if (system.insideBlockedPockets(componentBRef, fractionalMoleculeB))
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> nonEwaldDifferenceB =
        computeNonEwaldEnergyDifference(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                        system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Lambda-NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Lambda-NonEwald"] += (time_end - time_begin);
    if (!nonEwaldDifferenceB.has_value())
    {
      restoreFractionalPair();
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    energyDifference += nonEwaldDifferenceB.value();

    time_begin = std::chrono::system_clock::now();
    energyDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.totalEik, system.forceField,
        system.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB, runningNetCharge,
        groupChargeSum(fractionalMoleculeA));
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Lambda-Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Lambda-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    energyDifference += computeTailEnergyDifference(system, fractionalMoleculeB, oldFractionalMoleculeB,
                                                    system.spanOfMoleculeAtoms());
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["Lambda-Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Lambda-Tail"] += (time_end - time_begin);

    componentA.mc_moves_statistics.addConstructed(move, 2);

    const double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];

    if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + biasTerm))
    {
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      componentA.mc_moves_statistics.addAccepted(move, 2);

      lambda.setCurrentBin(newBin);

      return {energyDifference, double3(0.0, 1.0, 0.0)};
    }

    restoreFractionalPair();
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
}
