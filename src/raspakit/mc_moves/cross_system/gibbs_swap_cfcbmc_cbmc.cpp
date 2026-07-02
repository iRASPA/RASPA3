module;

module mc_moves_gibbs_swap_cfcbmc_cbmc;

import std;

import randomnumbers;
import running_energy;
import system;
import molecule;
import atom;
import cbmc;
import cbmc_chain_data;
import energy_factor;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import running_energy;
import forcefield;
import component;
import simulationbox;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;

// All systems have a fractional molecule, only one of these is 'active', the others are switched off with 'lambda=0'.
// Serial CFCMC: at most one simulation box has containsTheFractionalMolecule==true; systemA must be that box.
// Integer insertions and deletions use CBMC grow and retrace.

using GibbsSwapFractionalSnapshot = std::vector<std::pair<std::pair<std::size_t, std::size_t>, std::vector<Atom>>>;

GibbsSwapFractionalSnapshot saveGibbsSwapFractionalMolecules(const System& system)
{
  GibbsSwapFractionalSnapshot snapshot;
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }
    const std::size_t index =
        system.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, componentId);
    std::span<const Atom> fractionalMolecule = system.spanOfMolecule(componentId, index);
    snapshot.emplace_back(std::make_pair(componentId, index),
                          std::vector<Atom>(fractionalMolecule.begin(), fractionalMolecule.end()));
  }
  return snapshot;
}

std::vector<Atom> flattenGibbsSwapFractionalSnapshot(const GibbsSwapFractionalSnapshot& snapshot)
{
  std::vector<Atom> atoms;
  for (const auto& [indices, moleculeAtoms] : snapshot)
  {
    (void)indices;
    atoms.insert(atoms.end(), moleculeAtoms.begin(), moleculeAtoms.end());
  }
  return atoms;
}

void restoreGibbsSwapFractionalMolecules(System& system, const GibbsSwapFractionalSnapshot& snapshot)
{
  for (const auto& [indices, moleculeAtoms] : snapshot)
  {
    const auto& [componentId, index] = indices;
    std::span<Atom> fractionalMolecule = system.spanOfMolecule(componentId, index);
    std::copy(moleculeAtoms.begin(), moleculeAtoms.end(), fractionalMolecule.begin());
  }
}

std::vector<Atom> buildSerialFlagSwapTrialFractionalAtoms(const System& system, std::size_t selectedComponent,
                                                          bool systemBecomesActive, double selectedActiveLambda,
                                                          const System& lambdaSourceSystem)
{
  std::vector<Atom> atoms;
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }

    const std::size_t indexFractionalMolecule =
        system.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, componentId);
    std::span<const Atom> fractionalMolecule = system.spanOfMolecule(componentId, indexFractionalMolecule);
    for (Atom atom : fractionalMolecule)
    {
      atom.isFractional = true;
      if (systemBecomesActive)
      {
        if (componentId == selectedComponent)
        {
          atom.setScaling(selectedActiveLambda);
          atom.groupId = static_cast<std::uint8_t>(system.components[componentId].lambdaGC.computeDUdlambda);
        }
        else
        {
          const double lambda = lambdaSourceSystem.components[componentId].lambdaGC.lambdaValue();
          atom.setScaling(lambda);
          atom.groupId = static_cast<std::uint8_t>(system.components[componentId].lambdaGC.computeDUdlambda);
        }
      }
      else
      {
        // inactive fractional molecules must not contribute dUdlambda
        atom.setScalingFullyOff();
        atom.groupId = std::uint8_t{0};
      }
      atoms.push_back(atom);
    }
  }
  return atoms;
}

std::optional<RunningEnergy> computeGibbsSwapFractionalMoleculesEnergyDifference(
    System& system, std::span<const Atom> newAtomsCombined, std::span<const Atom> oldAtomsCombined) noexcept
{
  std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), newAtomsCombined, oldAtomsCombined);
  if (!frameworkDifference.has_value())
  {
    return std::nullopt;
  }

  std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), newAtomsCombined, oldAtomsCombined);
  if (!moleculeDifference.has_value())
  {
    return std::nullopt;
  }

  RunningEnergy ewaldDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newAtomsCombined, oldAtomsCombined, system.netCharge);

  RunningEnergy tailDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), newAtomsCombined,
                                                              oldAtomsCombined) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newAtomsCombined,
                                                                 oldAtomsCombined);

  return frameworkDifference.value() + moleculeDifference.value() + ewaldDifference + tailDifference;
}

std::optional<RunningEnergy> computeSerialFlagSwapFractionalEnergyDifference(
    System& system, const GibbsSwapFractionalSnapshot& beforeSnapshot, std::size_t selectedComponent,
    bool systemBecomesActive, double selectedActiveLambda, const System& lambdaSourceSystem) noexcept
{
  const std::vector<Atom> oldAtomsCombined = flattenGibbsSwapFractionalSnapshot(beforeSnapshot);
  const std::vector<Atom> newAtomsCombined = buildSerialFlagSwapTrialFractionalAtoms(
      system, selectedComponent, systemBecomesActive, selectedActiveLambda, lambdaSourceSystem);
  return computeGibbsSwapFractionalMoleculesEnergyDifference(system, newAtomsCombined, oldAtomsCombined);
}

std::optional<RunningEnergy> computeSerialFlagSwapOtherComponentsEnergyDifference(
    System& system, const GibbsSwapFractionalSnapshot& beforeSnapshot, std::size_t selectedComponent,
    bool systemBecomesActive, const System& lambdaSourceSystem) noexcept
{
  std::vector<Atom> oldAtomsCombined;
  std::vector<Atom> newAtomsCombined;

  for (const auto& [indices, moleculeAtoms] : beforeSnapshot)
  {
    const auto& [componentId, index] = indices;
    (void)index;
    if (componentId == selectedComponent)
    {
      continue;
    }
    oldAtomsCombined.insert(oldAtomsCombined.end(), moleculeAtoms.begin(), moleculeAtoms.end());

    std::span<const Atom> fractionalMolecule = system.spanOfMolecule(componentId, index);
    for (Atom atom : fractionalMolecule)
    {
      atom.isFractional = true;
      if (systemBecomesActive)
      {
        const double lambda = lambdaSourceSystem.components[componentId].lambdaGC.lambdaValue();
        atom.setScaling(lambda);
        atom.groupId = static_cast<std::uint8_t>(system.components[componentId].lambdaGC.computeDUdlambda);
      }
      else
      {
        // inactive fractional molecules must not contribute dUdlambda
        atom.setScalingFullyOff();
        atom.groupId = std::uint8_t{0};
      }
      newAtomsCombined.push_back(atom);
    }
  }

  if (oldAtomsCombined.empty())
  {
    return RunningEnergy{};
  }

  return computeGibbsSwapFractionalMoleculesEnergyDifference(system, newAtomsCombined, oldAtomsCombined);
}

void applySerialGibbsFlagSwapScaling(System& system, std::size_t selectedComponent, bool systemBecomesActive,
                                     double selectedActiveLambda, const System& lambdaSourceSystem)
{
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }

    const std::size_t indexFractionalMolecule =
        system.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, componentId);
    std::span<Atom> fractionalMolecule = system.spanOfMolecule(componentId, indexFractionalMolecule);
    for (Atom& atom : fractionalMolecule)
    {
      atom.isFractional = true;
      if (systemBecomesActive)
      {
        if (componentId == selectedComponent)
        {
          atom.setScaling(selectedActiveLambda);
          atom.groupId = static_cast<std::uint8_t>(system.components[componentId].lambdaGC.computeDUdlambda);
        }
        else
        {
          const double lambda = lambdaSourceSystem.components[componentId].lambdaGC.lambdaValue();
          atom.setScaling(lambda);
          atom.groupId = static_cast<std::uint8_t>(system.components[componentId].lambdaGC.computeDUdlambda);
        }
      }
      else
      {
        // inactive fractional molecules must not contribute dUdlambda
        atom.setScalingFullyOff();
        atom.groupId = std::uint8_t{0};
      }
    }
  }
}

void applyOtherComponentFlagSwapScaling(System& systemBecomingInactive, System& systemBecomingActive,
                                        std::size_t selectedComponent, double selectedActiveLambda)
{
  applySerialGibbsFlagSwapScaling(systemBecomingInactive, selectedComponent, false, selectedActiveLambda,
                                  systemBecomingInactive);
  applySerialGibbsFlagSwapScaling(systemBecomingActive, selectedComponent, true, selectedActiveLambda,
                                  systemBecomingInactive);
}

void syncOtherComponentLambdaBins(System& systemBecomingInactive, System& systemBecomingActive,
                                  std::size_t selectedComponent)
{
  for (std::size_t componentId = 0; componentId < systemBecomingInactive.components.size(); ++componentId)
  {
    if (componentId == selectedComponent ||
        systemBecomingInactive.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }

    systemBecomingActive.components[componentId].lambdaGC.setCurrentBin(
        systemBecomingInactive.components[componentId].lambdaGC.currentBin);
  }
}

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsSwapMove_CBCFCMC(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent,
    [[maybe_unused]] std::size_t& fractionalMoleculeSystem)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::GibbsSwapCBCFCMC;
  Component& componentA = systemA.components[selectedComponent];
  Component& componentB = systemB.components[selectedComponent];

  PropertyLambdaProbabilityHistogram& lambdaA = componentA.lambdaGC;
  PropertyLambdaProbabilityHistogram& lambdaB = componentB.lambdaGC;
  std::size_t oldBin = lambdaA.currentBin;
  double deltaLambda = lambdaA.delta;
  double oldLambda = componentA.lambdaGC.lambdaValue();

  double maxChange = componentA.mc_moves_statistics.getMaxChange(move, 2);
  std::make_signed_t<std::size_t> selectedNewBin = lambdaA.selectNewBin(random, maxChange);

  double switchValue = random.uniform();

  double cutOffFrameworkVDWA = systemA.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDWA = systemA.forceField.cutOffMoleculeVDW;
  double cutOffCoulombA = systemA.forceField.cutOffCoulomb;
  double cutOffFrameworkVDWB = systemB.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDWB = systemB.forceField.cutOffMoleculeVDW;
  double cutOffCoulombB = systemB.forceField.cutOffCoulomb;
  Component::GrowType growType = componentA.growType;

  if (!systemA.containsTheFractionalMolecule || systemB.containsTheFractionalMolecule)
  {
    return std::nullopt;
  }

  if (switchValue < 0.25)
  {
    // Lambda interchange: fractional molecule exchange between boxes with CBMC grow (A) and retrace (B).

    componentA.mc_moves_statistics.addTrial(move, 0);

    if (systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

    std::size_t indexFractionalMoleculeA =
        systemA.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, selectedComponent);
    std::size_t indexFractionalMoleculeB =
        systemB.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);

    const std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    const std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
    const GibbsSwapFractionalSnapshot snapshotA = saveGibbsSwapFractionalMolecules(systemA);
    const GibbsSwapFractionalSnapshot snapshotB = saveGibbsSwapFractionalMolecules(systemB);

    // System A: copy inactive fractional molecule from B into the active fractional slot
    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
    for (Atom& atom : fractionalMoleculeA)
    {
      atom.moleculeId = static_cast<std::uint32_t>(
          systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
    }

    RunningEnergy intraFractionalDifferenceA =
        componentA.intraMolecularPotentials.computeInternalEnergies(fractionalMoleculeA) -
        componentA.intraMolecularPotentials.computeInternalEnergies(oldFractionalMoleculeA);

    RunningEnergy energyDifferenceA = intraFractionalDifferenceA;

    // Combined flag-swap energy for all fractional molecules in A (selected: new position + off; others: off),
    // computed and applied to memory before the grow so the grow sees the post-flag-swap background.
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> flagSwapFractionalDifferenceA = computeSerialFlagSwapFractionalEnergyDifference(
        systemA, snapshotA, selectedComponent, false, oldLambda, systemA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    if (!flagSwapFractionalDifferenceA.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }
    energyDifferenceA += flagSwapFractionalDifferenceA.value();

    applySerialGibbsFlagSwapScaling(systemA, selectedComponent, false, oldLambda, systemA);

    // System A: CBMC grow a new integer molecule
    // Trial global molecule id (unique, not yet in the system); insertMolecule assigns the final id.
    std::size_t newMoleculeIndex = systemA.numberOfMolecules();
    time_begin = std::chrono::system_clock::now();
    std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
        random, componentA, selectedComponent, systemA.hasExternalField, systemA.forceField, systemA.simulationBox,
        systemA.interpolationGrids, systemA.externalFieldInterpolationGrid, systemA.framework,
        systemA.spanOfFrameworkAtoms(), systemA.spanOfMoleculeAtoms(), systemA.beta, growType, cutOffFrameworkVDWA,
        cutOffMoleculeVDWA, cutOffCoulombA, newMoleculeIndex, 1.0, false, false);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);

    if (!growData)
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }

    if (systemA.insideBlockedPockets(componentA, growData->atom))
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifferenceGrowA = Interactions::energyDifferenceEwaldFourier(
        systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.totalEik, systemA.totalEik,
        systemA.forceField, systemA.simulationBox, std::span(growData->atom.begin(), growData->atom.end()), {},
        systemA.netCharge);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceGrowA =
        Interactions::computeInterMolecularTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                                systemA.spanOfMoleculeAtoms(),
                                                                std::span(growData->atom.begin(), growData->atom.end()),
                                                                {}) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(
            systemA.forceField, systemA.simulationBox, systemA.spanOfFrameworkAtoms(),
            std::span(growData->atom.begin(), growData->atom.end()), {});
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);

    energyDifferenceA += growData->energies + EwaldFourierDifferenceGrowA + tailEnergyDifferenceGrowA;

    double correctionFactorEwaldGrowA = std::exp(
        -systemA.beta * (EwaldFourierDifferenceGrowA.potentialEnergy() + tailEnergyDifferenceGrowA.potentialEnergy()));

    // System B: CBMC retrace a randomly selected integer molecule
    std::size_t indexSelectedIntegerMoleculeB = systemB.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> selectedIntegerMoleculeB = systemB.spanOfMolecule(selectedComponent, indexSelectedIntegerMoleculeB);
    std::vector<Atom> oldSelectedIntegerMoleculeB(selectedIntegerMoleculeB.begin(), selectedIntegerMoleculeB.end());

    time_begin = std::chrono::system_clock::now();
    ChainRetraceData retraceData = CBMC::retraceMoleculeSwapDeletion(
        random, componentB, systemB.hasExternalField, systemB.forceField, systemB.simulationBox,
        systemB.interpolationGrids, systemB.externalFieldInterpolationGrid, systemB.framework,
        systemB.spanOfFrameworkAtoms(), systemB.spanOfMoleculeAtoms(), systemB.beta, growType, cutOffFrameworkVDWB,
        cutOffMoleculeVDWB, cutOffCoulombB, selectedIntegerMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifferenceRetraceB = Interactions::energyDifferenceEwaldFourier(
        systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.totalEik,
        systemB.forceField, systemB.simulationBox, {}, selectedIntegerMoleculeB, systemB.netCharge);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceRetraceB =
        Interactions::computeInterMolecularTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                                systemB.spanOfMoleculeAtoms(), {},
                                                                selectedIntegerMoleculeB) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                                   systemB.spanOfFrameworkAtoms(), {},
                                                                   selectedIntegerMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);

    double correctionFactorEwaldRetraceB = std::exp(
        systemB.beta *
        (EwaldFourierDifferenceRetraceB.potentialEnergy() + tailEnergyDifferenceRetraceB.potentialEnergy()));

    for (Atom& atom : selectedIntegerMoleculeB)
    {
      atom.setScalingFullyOff();
      atom.groupId = std::uint8_t{0};
      atom.position = systemB.simulationBox.randomPosition(random);
    }

    std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.moleculeId = static_cast<std::uint16_t>(
          systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB));
      atom.setScaling(oldLambda);
      atom.groupId = static_cast<std::uint8_t>(componentB.lambdaGC.computeDUdlambda);
      atom.isFractional = true;
    }

    RunningEnergy intraFractionalDifferenceB =
        componentB.intraMolecularPotentials.computeInternalEnergies(fractionalMoleculeB) -
        componentB.intraMolecularPotentials.computeInternalEnergies(oldFractionalMoleculeB);

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceB2 = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
        systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    if (!frameworkDifferenceB2.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceB2 = Interactions::computeInterMolecularEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.spanOfMoleculeAtoms(), fractionalMoleculeB,
        oldFractionalMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    if (!moleculeDifferenceB2.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifferenceB2 = Interactions::energyDifferenceEwaldFourier(
        systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.totalEik, systemB.totalEik,
        systemB.forceField, systemB.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB, systemB.netCharge);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceB2 =
        Interactions::computeInterMolecularTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                                systemB.spanOfMoleculeAtoms(), fractionalMoleculeB,
                                                                oldFractionalMoleculeB) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                                   systemB.spanOfFrameworkAtoms(), fractionalMoleculeB,
                                                                   oldFractionalMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);

    std::optional<RunningEnergy> flagSwapOtherComponentsDifferenceB =
        computeSerialFlagSwapOtherComponentsEnergyDifference(systemB, snapshotB, selectedComponent, true, systemA);
    if (!flagSwapOtherComponentsDifferenceB.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy energyDifferenceB =
        -(retraceData.energies - EwaldFourierDifferenceRetraceB - tailEnergyDifferenceRetraceB) +
        frameworkDifferenceB2.value() + moleculeDifferenceB2.value() + EwaldFourierDifferenceB2 +
        tailEnergyDifferenceB2 + intraFractionalDifferenceB + flagSwapOtherComponentsDifferenceB.value();

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];
    double preFactor = static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent]));
    double idealGasRosenbluthWeight = componentA.idealGasRosenbluthWeight.value_or(1.0);

    componentA.mc_moves_statistics.addConstructed(move, 0);

    if (random.uniform() <
        preFactor * (growData->RosenbluthWeight / idealGasRosenbluthWeight) / retraceData.RosenbluthWeight *
            correctionFactorEwaldGrowA * correctionFactorEwaldRetraceB *
            std::exp(-systemA.beta * (energyDifferenceA.potentialEnergy() + energyDifferenceB.potentialEnergy()) +
                     biasTerm))
    {
      componentA.mc_moves_statistics.addAccepted(move, 0);

      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());

      syncOtherComponentLambdaBins(systemA, systemB, selectedComponent);
      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
      for (Atom& atom : fractionalMoleculeA)
      {
        atom.moleculeId = static_cast<std::uint16_t>(
            systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
        atom.isFractional = true;
      }

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);
      systemA.insertMolecule(selectedComponent, growData->molecule, growData->atom);
      systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)] =
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)];

      std::swap(
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)],
          systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexSelectedIntegerMoleculeB)]);
      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);
      systemB.deleteMolecule(selectedComponent, indexSelectedIntegerMoleculeB, selectedIntegerMoleculeB);

      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
      for (Atom& atom : fractionalMoleculeB)
      {
        atom.moleculeId = static_cast<std::uint16_t>(
            systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB));
        atom.isFractional = true;
      }

      const double activeLambda = componentB.lambdaGC.lambdaValue();
      applyOtherComponentFlagSwapScaling(systemA, systemB, selectedComponent, activeLambda);

      systemA.updateMoleculeAtomInformation();
      systemB.updateMoleculeAtomInformation();

      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
    std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), selectedIntegerMoleculeB.begin());

    return std::nullopt;
  }
  else if (switchValue < 0.5)
  {
    // Shuffle: move active fractional molecule from A to B using CBMC grow in B.

    componentA.mc_moves_statistics.addTrial(move, 1);

    std::size_t indexFractionalMoleculeA =
        systemA.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, selectedComponent);
    std::size_t indexFractionalMoleculeB =
        systemB.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);

    std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
    const GibbsSwapFractionalSnapshot snapshotA = saveGibbsSwapFractionalMolecules(systemA);
    const GibbsSwapFractionalSnapshot snapshotB = saveGibbsSwapFractionalMolecules(systemB);

    std::swap_ranges(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), fractionalMoleculeB.begin());

    // swap_ranges copies moleculeId across systems; restore the ids valid in each system so that
    // self-exclusion in the span-based energy differences and the CBMC grow works correctly.
    // System A becomes inactive: its fractional must not contribute dUdlambda (groupId 0).
    for (Atom& atom : fractionalMoleculeA)
    {
      atom.isFractional = true;
      atom.groupId = std::uint8_t{0};
      atom.moleculeId = static_cast<std::uint32_t>(
          systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA));
    }
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.isFractional = true;
      atom.groupId = static_cast<std::uint8_t>(componentB.lambdaGC.computeDUdlambda);
      atom.moleculeId = static_cast<std::uint32_t>(
          systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB));
    }

    RunningEnergy intraEnergyDifferenceA =
        componentA.intraMolecularPotentials.computeInternalEnergies(fractionalMoleculeA) -
        componentA.intraMolecularPotentials.computeInternalEnergies(oldFractionalMoleculeA);

    // Combined flag-swap energy for all fractional molecules in A (selected: new position + off; others: off),
    // computed in one batch to count the selected-other fractional cross-terms exactly once, and applied to
    // memory so the remainder of the move sees the post-flag-swap background.
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> flagSwapFractionalDifferenceA = computeSerialFlagSwapFractionalEnergyDifference(
        systemA, snapshotA, selectedComponent, false, oldLambda, systemA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);
    if (!flagSwapFractionalDifferenceA.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy energyDifferenceA = intraEnergyDifferenceA + flagSwapFractionalDifferenceA.value();

    applySerialGibbsFlagSwapScaling(systemA, selectedComponent, false, oldLambda, systemA);

    time_begin = std::chrono::system_clock::now();
    const std::size_t globalFractionalMoleculeIndexB =
        systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB);
    std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
        random, componentB, selectedComponent, systemB.hasExternalField, systemB.forceField, systemB.simulationBox,
        systemB.interpolationGrids, systemB.externalFieldInterpolationGrid, systemB.framework,
        systemB.spanOfFrameworkAtoms(), systemB.spanOfMoleculeAtoms(), systemB.beta, growType, cutOffFrameworkVDWB,
        cutOffMoleculeVDWB, cutOffCoulombB, globalFractionalMoleculeIndexB, oldLambda,
        componentB.lambdaGC.computeDUdlambda, true);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);

    if (!growData)
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    if (systemB.insideBlockedPockets(componentB, growData->atom))
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    std::copy(growData->atom.begin(), growData->atom.end(), fractionalMoleculeB.begin());
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.moleculeId = static_cast<std::uint32_t>(globalFractionalMoleculeIndexB);
    }

    RunningEnergy intraEnergyDifference =
        componentB.intraMolecularPotentials.computeInternalEnergies(fractionalMoleculeB) -
        componentB.intraMolecularPotentials.computeInternalEnergies(oldFractionalMoleculeB);

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceB = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
        systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);
    if (!frameworkDifferenceB.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceB = Interactions::computeInterMolecularEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.spanOfMoleculeAtoms(), fractionalMoleculeB,
        oldFractionalMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);
    if (!moleculeDifferenceB.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifferenceB = Interactions::energyDifferenceEwaldFourier(
        systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.totalEik,
        systemB.forceField, systemB.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB, systemB.netCharge);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-Ewald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceB =
        Interactions::computeInterMolecularTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                                systemB.spanOfMoleculeAtoms(), fractionalMoleculeB,
                                                                oldFractionalMoleculeB) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                                   systemB.spanOfFrameworkAtoms(), fractionalMoleculeB,
                                                                   oldFractionalMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-Tail"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-Tail"] += (time_end - time_begin);

    std::optional<RunningEnergy> flagSwapOtherComponentsDifferenceB =
        computeSerialFlagSwapOtherComponentsEnergyDifference(systemB, snapshotB, selectedComponent, true, systemA);
    if (!flagSwapOtherComponentsDifferenceB.has_value())
    {
      restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    RunningEnergy energyDifferenceB = frameworkDifferenceB.value() + moleculeDifferenceB.value() +
                                      EwaldEnergyDifferenceB + tailEnergyDifferenceB + intraEnergyDifference +
                                      flagSwapOtherComponentsDifferenceB.value();

    componentA.mc_moves_statistics.addConstructed(move, 1);

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];
    double preFactor = systemB.simulationBox.volume / systemA.simulationBox.volume;
    double idealGasRosenbluthWeight = componentB.idealGasRosenbluthWeight.value_or(1.0);

    if (random.uniform() < preFactor * (growData->RosenbluthWeight / idealGasRosenbluthWeight) *
                                std::exp(-systemA.beta * (energyDifferenceA.potentialEnergy() +
                                                           energyDifferenceB.potentialEnergy()) +
                                         biasTerm))
    {
      componentA.mc_moves_statistics.addAccepted(move, 1);

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);
      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);

      std::swap(systemA.moleculeData[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)],
                systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)]);
      systemB.moleculeData[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)] =
          growData->molecule;
      std::copy(growData->atom.begin(), growData->atom.end(), fractionalMoleculeB.begin());

      syncOtherComponentLambdaBins(systemA, systemB, selectedComponent);
      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

      const double activeLambda = componentB.lambdaGC.lambdaValue();
      applyOtherComponentFlagSwapScaling(systemA, systemB, selectedComponent, activeLambda);

      systemA.updateMoleculeAtomInformation();
      systemB.updateMoleculeAtomInformation();

      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    restoreGibbsSwapFractionalMolecules(systemA, snapshotA);
    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());

    return std::nullopt;
  }
  else  // lambda move
  {
    componentA.mc_moves_statistics.addTrial(move, 2);

    if (selectedNewBin < 0) return std::nullopt;
    if (selectedNewBin >= std::make_signed_t<std::size_t>(lambdaA.numberOfSamplePoints)) return std::nullopt;

    std::size_t newBin = static_cast<std::size_t>(selectedNewBin);
    double newLambda = deltaLambda * static_cast<double>(newBin);

    std::size_t indexFractionalMoleculeA =
        systemA.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCBCFCMC, selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);

    std::vector<Atom> trialPositions(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::transform(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), trialPositions.begin(),
                   [&](Atom a)
                   {
                     a.setScaling(newLambda);
                     return a;
                   });

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkEnergyDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemA.forceField, systemA.simulationBox, systemA.interpolationGrids, systemA.framework,
        systemA.spanOfFrameworkAtoms(), trialPositions, fractionalMoleculeA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaChange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaChange-NonEwald"] += (time_end - time_begin);

    if (!frameworkEnergyDifference.has_value()) return std::nullopt;

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> interEnergyDifference = Interactions::computeInterMolecularEnergyDifference(
        systemA.forceField, systemA.simulationBox, systemA.spanOfMoleculeAtoms(), trialPositions, fractionalMoleculeA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaChange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaChange-NonEwald"] += (time_end - time_begin);

    if (!interEnergyDifference.has_value()) return std::nullopt;

    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifference = Interactions::energyDifferenceEwaldFourier(
        systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.totalEik,
        systemA.forceField, systemA.simulationBox, trialPositions, fractionalMoleculeA, systemA.netCharge);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaChange-Ewald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaChange-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifference = Interactions::computeInterMolecularTailEnergyDifference(
                                             systemA.forceField, systemA.simulationBox, systemA.spanOfMoleculeAtoms(),
                                             trialPositions, fractionalMoleculeA) +
                                         Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                             systemA.forceField, systemA.simulationBox, systemA.spanOfFrameworkAtoms(),
                                             trialPositions, fractionalMoleculeA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaChange-Tail"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaChange-Tail"] += (time_end - time_begin);

    RunningEnergy energyDifference = frameworkEnergyDifference.value() + interEnergyDifference.value() +
                                     EwaldFourierDifference + tailEnergyDifference;

    componentA.mc_moves_statistics.addConstructed(move, 2);

    double biasTerm = lambdaA.biasFactor[newBin] - lambdaA.biasFactor[oldBin];

    if (random.uniform() < std::exp(-systemA.beta * energyDifference.potentialEnergy() + biasTerm))
    {
      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);

      componentA.mc_moves_statistics.addAccepted(move, 2);

      std::copy(trialPositions.begin(), trialPositions.end(), fractionalMoleculeA.begin());

      componentA.lambdaGC.setCurrentBin(newBin);

      return std::make_pair(energyDifference, RunningEnergy());
    };

    return std::nullopt;
  }

  return std::nullopt;
}
