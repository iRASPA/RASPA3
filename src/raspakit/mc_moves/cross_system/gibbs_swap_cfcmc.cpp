module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#endif

module mc_moves_gibbs_swap_cfcmc;

#ifndef USE_LEGACY_HEADERS
import <optional>;
import <span>;
import <chrono>;
import <vector>;
import <cmath>;
import <tuple>;
import <algorithm>;
import <utility>;
import <type_traits>;
#endif

import randomnumbers;
import running_energy;
import system;
import molecule;
import atom;
import cbmc;
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
// Implementation advantage: the number of fractional molecules per system remains constant.

// systemA contains the fractional molecule
std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsSwapMove_CFCMC(
    RandomNumber& random, System& systemA, System& systemB, size_t selectedComponent,
    [[maybe_unused]] size_t& fractionalMoleculeSystem)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::GibbsSwapCFCMC;
  Component& componentA = systemA.components[selectedComponent];
  Component& componentB = systemB.components[selectedComponent];

  PropertyLambdaProbabilityHistogram& lambdaA = componentA.lambdaGC;
  PropertyLambdaProbabilityHistogram& lambdaB = componentB.lambdaGC;
  size_t oldBin = lambdaA.currentBin;
  double deltaLambda = lambdaA.delta;
  double oldLambda = componentA.lambdaGC.lambdaValue();

  double maxChange = componentA.mc_moves_statistics.getMaxChange(move, 2);
  std::make_signed_t<std::size_t> selectedNewBin = lambdaA.selectNewBin(random, maxChange);

  double switchValue = random.uniform();

  // if (selectedNewBin >= std::make_signed_t<std::size_t>(lambdaA.numberOfBins))
  if (switchValue < 0.25)
  {
    // Swap move:
    // Changing the fractional molecule into a whole molecule, keeping its position fixed
    // Changing a randomly selected molecule in the other simulation box into a fractional molecule (at same lambda)

    componentA.mc_moves_statistics.addTrial(move, 0);

    if (systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

    size_t indexFractionalMoleculeA = systemA.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
    size_t indexFractionalMoleculeB = systemB.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);

    // assert(fractionalMoleculeA.front().groupId == uint8_t{ 1 });
    // assert(fractionalMoleculeB.front().groupId == uint8_t{ 1 });

    // make copy of old fractional molecule for reference and restoring
    const std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    const std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
    std::vector<Atom> oldFractionalMoleculeB2(fractionalMoleculeB.begin(), fractionalMoleculeB.end());

    // System A: Changing the fractional molecule into a whole molecule, keeping its position fixed
    //=============================================================================================

    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
    for (Atom& atom : fractionalMoleculeA)
    {
      atom.moleculeId = static_cast<uint32_t>(indexFractionalMoleculeA);
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceA = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemA.forceField, systemA.simulationBox, systemA.interpolationGrids, systemA.framework,
        systemA.spanOfFrameworkAtoms(), fractionalMoleculeA, oldFractionalMoleculeA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);

    if (!frameworkDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceA = Interactions::computeInterMolecularEnergyDifference(
        systemA.forceField, systemA.simulationBox, systemA.spanOfMoleculeAtoms(), fractionalMoleculeA,
        oldFractionalMoleculeA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);

    if (!moleculeDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifferenceA = Interactions::energyDifferenceEwaldFourier(
        systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.totalEik,
        systemA.forceField, systemA.simulationBox, fractionalMoleculeA, oldFractionalMoleculeA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceA = Interactions::computeInterMolecularTailEnergyDifference(
                                              systemA.forceField, systemA.simulationBox, systemA.spanOfMoleculeAtoms(),
                                              fractionalMoleculeA, oldFractionalMoleculeA) +
                                          Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                              systemA.forceField, systemA.simulationBox, systemA.spanOfFrameworkAtoms(),
                                              fractionalMoleculeA, oldFractionalMoleculeA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);

    // step 2

    std::vector<Atom> newMolecule(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end());
    for (Atom& atom : newMolecule)
    {
      atom.setScalingToInteger();
      atom.moleculeId = static_cast<uint32_t>(systemA.numberOfMoleculesPerComponent[selectedComponent]);
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceA2 = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemA.forceField, systemA.simulationBox, systemA.interpolationGrids, systemA.framework,
        systemA.spanOfFrameworkAtoms(), newMolecule, {});
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);

    if (!frameworkDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceA2 = Interactions::computeInterMolecularEnergyDifference(
        systemA.forceField, systemA.simulationBox, systemA.spanOfMoleculeAtoms(), newMolecule, {});
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);

    if (!moleculeDifferenceA2.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifferenceA2 = Interactions::energyDifferenceEwaldFourier(
        systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.totalEik, systemA.totalEik,
        systemA.forceField, systemA.simulationBox, newMolecule, {});
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceA2 =
        Interactions::computeInterMolecularTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                                systemA.spanOfMoleculeAtoms(), newMolecule, {}) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                                   systemA.spanOfFrameworkAtoms(), newMolecule, {});
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);

    RunningEnergy energyDifferenceA = frameworkDifferenceA.value() + moleculeDifferenceA.value() +
                                      EwaldFourierDifferenceA + tailEnergyDifferenceA + frameworkDifferenceA2.value() +
                                      moleculeDifferenceA2.value() + EwaldFourierDifferenceA2 + tailEnergyDifferenceA2;

    // System B: Changing a randomly selected molecule in the other simulation box into a fractional molecule
    // (at the same lambda)
    //================================================================================================================

    size_t indexSelectedIntegerMoleculeB = systemB.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> selectedIntegerMoleculeB = systemB.spanOfMolecule(selectedComponent, indexSelectedIntegerMoleculeB);

    // make copy of selected molecule for reference and restoring
    std::vector<Atom> oldSelectedIntegerMoleculeB(selectedIntegerMoleculeB.begin(), selectedIntegerMoleculeB.end());
    std::vector<Atom> oldSelectedIntegerMoleculeB2(selectedIntegerMoleculeB.begin(), selectedIntegerMoleculeB.end());

    for (Atom& atom : selectedIntegerMoleculeB)
    {
      atom.setScalingFullyOff();
      atom.groupId = uint8_t{0};
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceB = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
        systemB.spanOfFrameworkAtoms(), selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    if (!frameworkDifferenceB.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceB = Interactions::computeInterMolecularEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.spanOfMoleculeAtoms(), selectedIntegerMoleculeB,
        oldSelectedIntegerMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    if (!moleculeDifferenceB.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifferenceB = Interactions::energyDifferenceEwaldFourier(
        systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.totalEik,
        systemB.forceField, systemB.simulationBox, selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceB = Interactions::computeInterMolecularTailEnergyDifference(
                                              systemB.forceField, systemB.simulationBox, systemB.spanOfMoleculeAtoms(),
                                              selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB) +
                                          Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                              systemB.forceField, systemB.simulationBox, systemB.spanOfFrameworkAtoms(),
                                              selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-Tail"] += (time_end - time_begin);

    std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
    for (Atom& atom : fractionalMoleculeB)
    {
      atom.moleculeId = static_cast<uint16_t>(indexFractionalMoleculeB);
      atom.setScaling(oldLambda);
      atom.groupId = uint8_t{1};
    }

    for (Atom& atom : selectedIntegerMoleculeB)
    {
      atom.setScalingFullyOff();
      atom.groupId = uint8_t{0};
      atom.position = systemA.simulationBox.randomPosition(random);
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceB2 = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
        systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaInterchange-NonEwald"] += (time_end - time_begin);
    if (!frameworkDifferenceB2.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
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
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifferenceB2 = Interactions::energyDifferenceEwaldFourier(
        systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.totalEik, systemB.totalEik,
        systemB.forceField, systemB.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB);
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

    RunningEnergy energyDifferenceB = frameworkDifferenceB.value() + moleculeDifferenceB.value() +
                                      EwaldFourierDifferenceB + tailEnergyDifferenceB + frameworkDifferenceB2.value() +
                                      moleculeDifferenceB2.value() + EwaldFourierDifferenceB2 + tailEnergyDifferenceB2;

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];

    double preFactor = static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent]));

    componentA.mc_moves_statistics.addConstructed(move, 0);

    // apply acceptance/rejection rule
    if (random.uniform() < preFactor * std::exp(-systemA.beta * (energyDifferenceA.potentialEnergy() +
                                                                 energyDifferenceB.potentialEnergy()) +
                                                biasTerm))
    {
      componentA.mc_moves_statistics.addAccepted(move, 0);

      // restore
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(),
                selectedIntegerMoleculeB.begin());

      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

      // copy the fractional molecule from B to A
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
      for (Atom& atom : fractionalMoleculeA)
      {
        atom.moleculeId = static_cast<uint16_t>(indexFractionalMoleculeA);
      }

      // make old fractional molecule integer
      std::vector<Atom> addedMolecule = std::vector<Atom>(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end());
      for (Atom& atom : addedMolecule)
      {
        atom.setScalingToInteger();
      }
      systemA.insertMolecule(
          selectedComponent,
          systemA.moleculePositions[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)],
          addedMolecule);
      systemA.moleculePositions[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)] =
          systemB.moleculePositions[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)];

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);

      std::swap(
          systemB.moleculePositions[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)],
          systemB
              .moleculePositions[systemB.moleculeIndexOfComponent(selectedComponent, indexSelectedIntegerMoleculeB)]);
      systemB.deleteMolecule(selectedComponent, indexSelectedIntegerMoleculeB, selectedIntegerMoleculeB);

      std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
      for (Atom& atom : fractionalMoleculeB)
      {
        atom.moleculeId = static_cast<uint16_t>(indexFractionalMoleculeB);
        atom.setScaling(oldLambda);
        atom.groupId = uint8_t{1};
      }

      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);

      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
    std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), selectedIntegerMoleculeB.begin());

    return std::nullopt;
  }
  // else if (selectedNewBin < 0)
  else if (switchValue < 0.5)
  {
    // Move fractional molecule to the other box (from A to a random position in B)

    componentA.mc_moves_statistics.addTrial(move, 1);

    size_t indexFractionalMoleculeA = systemA.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
    size_t indexFractionalMoleculeB = systemB.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);

    // make copy of old fractional molecule for reference and restoring
    std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());

    // swap the active and the inactive fractional molecule
    std::swap_ranges(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), fractionalMoleculeB.begin());

    std::pair<Molecule, std::vector<Atom>> trialMolecule =
        componentB.equilibratedMoleculeRandomInBox(random, systemB.simulationBox);

    // std::copy(trialMolecule.second.begin(), trialMolecule.second.end(), fractionalMoleculeB.begin());
    std::transform(fractionalMoleculeB.begin(), fractionalMoleculeB.end(), trialMolecule.second.begin(),
                   fractionalMoleculeB.begin(),
                   [](const Atom& a, const Atom& b)
                   {
                     return Atom(b.position, a.charge, a.scalingVDW, a.scalingCoulomb, a.moleculeId, a.type,
                                 a.componentId, a.groupId);
                   });

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceA = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemA.forceField, systemA.simulationBox, systemA.interpolationGrids, systemA.framework,
        systemA.spanOfFrameworkAtoms(), fractionalMoleculeA, oldFractionalMoleculeA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);

    if (!frameworkDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceA = Interactions::computeInterMolecularEnergyDifference(
        systemA.forceField, systemA.simulationBox, systemA.spanOfMoleculeAtoms(), fractionalMoleculeA,
        oldFractionalMoleculeA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);

    if (!moleculeDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifferenceA = Interactions::energyDifferenceEwaldFourier(
        systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.totalEik,
        systemA.forceField, systemA.simulationBox, fractionalMoleculeA, oldFractionalMoleculeA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-Ewald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceA = Interactions::computeInterMolecularTailEnergyDifference(
                                              systemA.forceField, systemA.simulationBox, systemA.spanOfMoleculeAtoms(),
                                              fractionalMoleculeA, oldFractionalMoleculeA) +
                                          Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                              systemA.forceField, systemA.simulationBox, systemA.spanOfFrameworkAtoms(),
                                              fractionalMoleculeA, oldFractionalMoleculeA);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-Tail"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-Tail"] += (time_end - time_begin);

    RunningEnergy energyDifferenceA =
        frameworkDifferenceA.value() + moleculeDifferenceA.value() + EwaldEnergyDifferenceA + tailEnergyDifferenceA;

    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceB = Interactions::computeFrameworkMoleculeEnergyDifference(
        systemB.forceField, systemB.simulationBox, systemB.interpolationGrids, systemB.framework,
        systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-NonEwald"] += (time_end - time_begin);

    if (!frameworkDifferenceB.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
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
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    time_begin = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifferenceB = Interactions::energyDifferenceEwaldFourier(
        systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.totalEik,
        systemB.forceField, systemB.simulationBox, fractionalMoleculeB, oldFractionalMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-Ewald"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-Ewald"] += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceB = Interactions::computeInterMolecularTailEnergyDifference(
                                              systemB.forceField, systemB.simulationBox, systemB.spanOfMoleculeAtoms(),
                                              fractionalMoleculeB, oldFractionalMoleculeB) +
                                          Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                              systemB.forceField, systemB.simulationBox, systemB.spanOfFrameworkAtoms(),
                                              fractionalMoleculeB, oldFractionalMoleculeB);
    time_end = std::chrono::system_clock::now();
    componentA.mc_moves_cputime[move]["LambdaShuffle-Tail"] += (time_end - time_begin);
    systemA.mc_moves_cputime[move]["LambdaShuffle-Tail"] += (time_end - time_begin);

    RunningEnergy energyDifferenceB =
        frameworkDifferenceB.value() + moleculeDifferenceB.value() + EwaldEnergyDifferenceB + tailEnergyDifferenceB;

    componentA.mc_moves_statistics.addConstructed(move, 1);

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];

    double preFactor = systemB.simulationBox.volume / systemA.simulationBox.volume;

    // apply acceptance/rejection rule
    if (random.uniform() <
        preFactor *
            exp(-systemA.beta * (energyDifferenceA.potentialEnergy() + energyDifferenceB.potentialEnergy()) + biasTerm))
    {
      componentA.mc_moves_statistics.addAccepted(move, 1);

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);
      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);

      std::swap(
          systemA.moleculePositions[systemA.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeA)],
          systemB.moleculePositions[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)]);
      systemB.moleculePositions[systemB.moleculeIndexOfComponent(selectedComponent, indexFractionalMoleculeB)] =
          trialMolecule.first;

      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(componentA.lambdaGC.currentBin, componentB.lambdaGC.currentBin);

      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    // reject, set fractional molecule back to old state
    std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());

    return std::nullopt;
  }
  else  // lambda move
  {
    componentA.mc_moves_statistics.addTrial(move, 2);

    if (selectedNewBin < 0) return std::nullopt;
    if (selectedNewBin >= std::make_signed_t<std::size_t>(lambdaA.numberOfSamplePoints)) return std::nullopt;

    size_t newBin = static_cast<size_t>(selectedNewBin);
    double newLambda = deltaLambda * static_cast<double>(newBin);

    size_t indexFractionalMoleculeA = systemA.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
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
        systemA.forceField, systemA.simulationBox, trialPositions, fractionalMoleculeA);
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

    // apply acceptance/rejection rule
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
