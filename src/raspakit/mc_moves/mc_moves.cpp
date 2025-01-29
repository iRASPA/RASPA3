module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <exception>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <print>
#include <source_location>
#include <span>
#include <tuple>
#include <vector>
#endif

module mc_moves;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <array>;
import <tuple>;
import <optional>;
import <span>;
import <optional>;
import <tuple>;
import <algorithm>;
import <chrono>;
import <cmath>;
import <iostream>;
import <iomanip>;
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <print>;
#endif

import archive;
import double3;
import double3x3;
import simd_quatd;
import randomnumbers;

import component;
import atom;
import molecule;
import simulationbox;
import cbmc;
import system;
import energy_status;
import energy_status_inter;
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_cputime;
import transition_matrix;
import mc_moves_translation;
import mc_moves_random_translation;
import mc_moves_rotation;
import mc_moves_random_rotation;
import mc_moves_reinsertion;
import mc_moves_insertion;
import mc_moves_deletion;
import mc_moves_insertion_cbmc;
import mc_moves_deletion_cbmc;
import mc_moves_swap_cfcmc;
import mc_moves_swap_cfcmc_cbmc;
import mc_moves_gibbs_swap_cbmc;
import mc_moves_volume;
import mc_moves_gibbs_volume;
import mc_moves_identity_change;
import mc_moves_swap_cfcmc;
import mc_moves_swap_cfcmc_cbmc;
import mc_moves_gibbs_swap_cbmc;
import mc_moves_gibbs_swap_cfcmc;
import mc_moves_reaction;
import mc_moves_reaction_cfcmc_cbmc;
import mc_moves_widom;
import mc_moves_parallel_tempering_swap;
import mc_moves_hybridmc;

void MC_Moves::performRandomMove(RandomNumber &random, System &selectedSystem, System &selectedSecondSystem,
                                 size_t selectedComponent, size_t &fractionalMoleculeSystem)
{
  // pick move type from probabilities object
  MoveTypes moveType = selectedSystem.components[selectedComponent].mc_moves_probabilities.sample(random);

  // save old number of molecules for reference
  size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  switch (moveType)
  {
    case MoveTypes::Translation:
    {
      // select molecule and only move if there are actually molecules
      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculePositions[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::translationMove(random, selectedSystem, selectedComponent, selectedMolecule,
                                      selectedSystem.components, molecule, molecule_atoms);

        // accept if energy difference is not 0
        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }

      break;
    }
    case MoveTypes::RandomTranslation:
    {
      // select molecule and only move if there are actually molecules
      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculePositions[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference = MC_Moves::randomTranslationMove(
            random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);

        // accept if energy difference is not 0
        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      break;
    }
    case MoveTypes::Rotation:
    {
      // select molecule and only move if there are actually molecules
      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculePositions[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference = MC_Moves::rotationMove(
            random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);

        // accept if energy difference is not 0
        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }

      break;
    }
    case MoveTypes::RandomRotation:
    {
      // select molecule and only move if there are actually molecules
      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculePositions[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference = MC_Moves::randomRotationMove(
            random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);

        // accept if energy difference is not 0
        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }

      break;
    }
    case MoveTypes::VolumeChange:
    {
      // perform move
      std::optional<RunningEnergy> energy = MC_Moves::volumeMove(random, selectedSystem);

      // accept if energy difference is not 0
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value();
      }
      break;
    }
    case MoveTypes::ReinsertionCBMC:
    {
      // select molecule and only move if there are actually molecules
      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculePositions[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference = MC_Moves::reinsertionMove(
            random, selectedSystem, selectedComponent, selectedMolecule, molecule, molecule_atoms);

        // accept if energy difference is not 0
        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }

        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      break;
    }
    case MoveTypes::IdentityChangeCBMC:
    {
      break;
    }
    case MoveTypes::Swap:
    {
      if (random.uniform() < 0.5)
      {
        // perform move
        const auto [energyDifference, Pacc] = MC_Moves::insertionMove(random, selectedSystem, selectedComponent);

        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      else
      {
        size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

        // perform move
        const auto [energyDifference, Pacc] =
            MC_Moves::deletionMove(random, selectedSystem, selectedComponent, selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      break;
    }
    case MoveTypes::SwapCBMC:
    {
      if (random.uniform() < 0.5)
      {
        const auto [energyDifference, Pacc] = MC_Moves::insertionMoveCBMC(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      else
      {
        size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

        const auto [energyDifference, Pacc] =
            MC_Moves::deletionMoveCBMC(random, selectedSystem, selectedComponent, selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      break;
    }
    case MoveTypes::SwapCFCMC:
    {
      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC(random, selectedSystem, selectedComponent, selectedMolecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      break;
    }
    case MoveTypes::SwapCBCFCMC:
    {
      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, selectedMolecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      break;
    }
    case MoveTypes::GibbsVolume:
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
          MC_Moves::GibbsVolumeMove(random, selectedSystem, selectedSecondSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value().first;
        selectedSecondSystem.runningEnergies = energy.value().second;
      }
      break;
    }
    case MoveTypes::GibbsSwapCBMC:
    {
      if (random.uniform() < 0.5)
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            MC_Moves::GibbsSwapMove_CBMC(random, selectedSystem, selectedSecondSystem, selectedComponent);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }
      }
      else
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            MC_Moves::GibbsSwapMove_CBMC(random, selectedSecondSystem, selectedSystem, selectedComponent);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }
      }
      break;
    }
    case MoveTypes::GibbsSwapCFCMC:
    {
      if (selectedSystem.containsTheFractionalMolecule)
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = MC_Moves::GibbsSwapMove_CFCMC(
            random, selectedSystem, selectedSecondSystem, fractionalMoleculeSystem, selectedComponent);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }
      }
      else if (selectedSecondSystem.containsTheFractionalMolecule)
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = MC_Moves::GibbsSwapMove_CFCMC(
            random, selectedSecondSystem, selectedSystem, fractionalMoleculeSystem, selectedComponent);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }
      }
      break;
    }
    case MoveTypes::Widom:
    {
      break;
    }
    case MoveTypes::WidomCFCMC:
    {
      break;
    }
    case MoveTypes::WidomCBCFCMC:
    {
      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, 0, true, true);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      break;
    }
    case MoveTypes::ParallelTempering:
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
          MC_Moves::ParallelTemperingSwap(random, selectedSystem, selectedSecondSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value().first;
        selectedSecondSystem.runningEnergies = energy.value().second;
      }
      break;
    }
    case MoveTypes::HybridMC:
    {
      std::optional<RunningEnergy> energy = MC_Moves::hybridMCMove(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value();
      }
      break;
    }
    case MoveTypes::Count:
    {
      throw std::runtime_error("Move count called, invalid sampling of move probabilities");
      break;
    }
    default:
    {
      throw std::runtime_error("No move called, invalid sampling of move probabilities");
      break;
    }
  }
}

void MC_Moves::performRandomMoveProduction(RandomNumber &random, System &selectedSystem, System &selectedSecondSystem,
                                           size_t selectedComponent, size_t &fractionalMoleculeSystem,
                                           size_t currentBlock)
{
  // pick move type from probabilities object
  MoveTypes moveType = selectedSystem.components[selectedComponent].mc_moves_probabilities.sample(random);

  // save old number of molecules for reference
  size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  switch (moveType)
  {
    case MoveTypes::Translation:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::Translation);

      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::translationMove(random, selectedSystem, selectedComponent, selectedMolecule,
                                      selectedSystem.components, molecule, molecule_atoms);
        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.translationMove += (t2 - t1);
      selectedSystem.mc_moves_cputime.translationMove += (t2 - t1);
      break;
    }
    case MoveTypes::RandomTranslation:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::RandomTranslation);

      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
        std::optional<RunningEnergy> energyDifference = MC_Moves::randomTranslationMove(
            random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);
        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.randomTranslationMove += (t2 - t1);
      selectedSystem.mc_moves_cputime.randomTranslationMove += (t2 - t1);
      break;
    }
    case MoveTypes::Rotation:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::Rotation);

      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
        std::optional<RunningEnergy> energyDifference = MC_Moves::rotationMove(
            random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);
        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.rotationMove += (t2 - t1);
      selectedSystem.mc_moves_cputime.rotationMove += (t2 - t1);
      break;
    }
    case MoveTypes::RandomRotation:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::RandomRotation);

      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
        std::optional<RunningEnergy> energyDifference = MC_Moves::randomRotationMove(
            random, selectedSystem, selectedComponent, selectedSystem.components, molecule, molecule_atoms);
        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.randomRotationMove += (t2 - t1);
      selectedSystem.mc_moves_cputime.randomRotationMove += (t2 - t1);
      break;
    }
    case MoveTypes::VolumeChange:
    {
      selectedSystem.mc_moves_statistics.addAllCounts(MoveTypes::VolumeChange);

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<RunningEnergy> energy = MC_Moves::volumeMove(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value();
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      
      selectedSystem.mc_moves_cputime.volumeMove += (t2 - t1);
      break;
    }
    case MoveTypes::ReinsertionCBMC:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::ReinsertionCBMC);

      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculePositions[molecule_index];
        std::optional<RunningEnergy> energyDifference = MC_Moves::reinsertionMove(
            random, selectedSystem, selectedComponent, selectedMolecule, molecule, molecule_atoms);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.reinsertionMoveCBMC += (t2 - t1);
      selectedSystem.mc_moves_cputime.reinsertionMoveCBMC += (t2 - t1);
      break;
    }
    case MoveTypes::IdentityChangeCBMC:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::IdentityChangeCBMC);
      break;
    }
    case MoveTypes::Swap:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::Swap);
      if (random.uniform() < 0.5)
      {
        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

        const auto [energyDifference, Pacc] = MC_Moves::insertionMove(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

        selectedSystem.components[selectedComponent].mc_moves_cputime.swapInsertionMove += (t2 - t1);
        selectedSystem.mc_moves_cputime.swapInsertionMove += (t2 - t1);
      }
      else
      {
        size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);
        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

        const auto [energyDifference, Pacc] =
            MC_Moves::deletionMove(random, selectedSystem, selectedComponent, selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

        selectedSystem.components[selectedComponent].mc_moves_cputime.swapDeletionMove += (t2 - t1);
        selectedSystem.mc_moves_cputime.swapDeletionMove += (t2 - t1);
      }
      break;
    }
    case MoveTypes::SwapCBMC:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::SwapCBMC);
      if (random.uniform() < 0.5)
      {
        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

        const auto [energyDifference, Pacc] = MC_Moves::insertionMoveCBMC(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

        selectedSystem.components[selectedComponent].mc_moves_cputime.swapInsertionMoveCBMC += (t2 - t1);
        selectedSystem.mc_moves_cputime.swapInsertionMoveCBMC += (t2 - t1);
      }
      else
      {
        size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);
        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

        const auto [energyDifference, Pacc] =
            MC_Moves::deletionMoveCBMC(random, selectedSystem, selectedComponent, selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

        selectedSystem.components[selectedComponent].mc_moves_cputime.swapDeletionMoveCBMC += (t2 - t1);
        selectedSystem.mc_moves_cputime.swapDeletionMoveCBMC += (t2 - t1);
      }
      break;
    }
    case MoveTypes::SwapCFCMC:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::SwapCFCMC);

      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC(random, selectedSystem, selectedComponent, selectedMolecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);

      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.swapLambdaMoveCFCMC += (t2 - t1);
      selectedSystem.mc_moves_cputime.swapLambdaMoveCFCMC += (t2 - t1);
      break;
    }
    case MoveTypes::SwapCBCFCMC:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::SwapCBCFCMC);

      size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, selectedMolecule);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);

      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.swapLambdaMoveCBCFCMC += (t2 - t1);
      selectedSystem.mc_moves_cputime.swapLambdaMoveCBCFCMC += (t2 - t1);
      break;
    }
    case MoveTypes::GibbsVolume:
    {
      selectedSystem.mc_moves_statistics.addAllCounts(MoveTypes::GibbsVolume);

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
          MC_Moves::GibbsVolumeMove(random, selectedSystem, selectedSecondSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value().first;
        selectedSecondSystem.runningEnergies = energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
      selectedSystem.mc_moves_cputime.GibbsVolumeMove += (t2 - t1);
      break;
    }
    case MoveTypes::GibbsSwapCBMC:
    {
      if (random.uniform() < 0.5)
      {
        selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::GibbsSwapCBMC);

        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            MC_Moves::GibbsSwapMove_CBMC(random, selectedSystem, selectedSecondSystem, selectedComponent);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }
        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

        selectedSystem.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMC += (t2 - t1);
        selectedSystem.mc_moves_cputime.GibbsSwapMoveCBMC += (t2 - t1);
      }
      else
      {
        selectedSecondSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::GibbsSwapCBMC);

        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            MC_Moves::GibbsSwapMove_CBMC(random, selectedSecondSystem, selectedSystem, selectedComponent);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }
        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

        selectedSecondSystem.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMC += (t2 - t1);
        selectedSecondSystem.mc_moves_cputime.GibbsSwapMoveCBMC += (t2 - t1);
      }
      break;
    }
    case MoveTypes::GibbsSwapCFCMC:
    {
      if (selectedSystem.containsTheFractionalMolecule)
      {
        selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::GibbsSwapCFCMC);

        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = MC_Moves::GibbsSwapMove_CFCMC(
            random, selectedSystem, selectedSecondSystem, selectedComponent, fractionalMoleculeSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }
        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

        selectedSystem.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaMoveCFCMC += (t2 - t1);
        selectedSystem.mc_moves_cputime.GibbsSwapLambdaMoveCFCMC += (t2 - t1);
      }
      else if (selectedSecondSystem.containsTheFractionalMolecule)
      {
        selectedSecondSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::GibbsSwapCFCMC);
        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = MC_Moves::GibbsSwapMove_CFCMC(
            random, selectedSecondSystem, selectedSystem, selectedComponent, fractionalMoleculeSystem);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }
        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

        selectedSecondSystem.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaMoveCFCMC += (t2 - t1);
        selectedSecondSystem.mc_moves_cputime.GibbsSwapLambdaMoveCFCMC += (t2 - t1);
      }
      break;
    }
    case MoveTypes::Widom:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::Widom);

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      std::pair<double, double> value = MC_Moves::WidomMove(random, selectedSystem, selectedComponent);

      selectedSystem.components[selectedComponent].averageRosenbluthWeights.addWidomSample(
          currentBlock, value.first, 2.0 * value.second, selectedSystem.weight());

      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.WidomMoveCBMC += (t2 - t1);
      selectedSystem.mc_moves_cputime.WidomMoveCBMC += (t2 - t1);
      break;
    }
    case MoveTypes::WidomCFCMC:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::WidomCFCMC);

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC(random, selectedSystem, selectedComponent, 0, true, true);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.WidomMoveCFCMC += (t2 - t1);
      selectedSystem.mc_moves_cputime.WidomMoveCFCMC += (t2 - t1);
      break;
    }
    case MoveTypes::WidomCBCFCMC:
    {
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(MoveTypes::WidomCBCFCMC);

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, 0, true, true);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }

      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.components[selectedComponent].mc_moves_cputime.WidomMoveCBCFCMC += (t2 - t1);
      selectedSystem.mc_moves_cputime.WidomMoveCBCFCMC += (t2 - t1);
      break;
    }
    case MoveTypes::ParallelTempering:
    {
      selectedSystem.mc_moves_statistics.addAllCounts(MoveTypes::ParallelTempering);

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
          MC_Moves::ParallelTemperingSwap(random, selectedSystem, selectedSecondSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value().first;
        selectedSecondSystem.runningEnergies = energy.value().second;
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.mc_moves_cputime.ParallelTemperingSwap += (t2 - t1);
      break;
    }
    case MoveTypes::HybridMC:
    {
      selectedSystem.mc_moves_statistics.addAllCounts(MoveTypes::HybridMC);

      std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
      std::optional<RunningEnergy> energy = MC_Moves::hybridMCMove(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value();
      }
      std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

      selectedSystem.mc_moves_cputime.hybridMC += (t2 - t1);
      break;
    }
    case MoveTypes::Count:
    {
      throw std::runtime_error("Move count called, invalid sampling of move probabilities");
      break;
    }
    default:
    {
      throw std::runtime_error("No move called, invalid sampling of move probabilities");
      break;
    }
  }
}
