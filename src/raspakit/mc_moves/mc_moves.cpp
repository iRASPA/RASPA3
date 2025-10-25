module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
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

#ifdef USE_STD_IMPORT
import std;
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
import mc_moves_partial_reinsertion;
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

void MC_Moves::performRandomMoveInitialization(RandomNumber &random, System &selectedSystem,
                                               System &selectedSecondSystem, std::size_t selectedComponent,
                                               [[maybe_unused]] std::size_t &fractionalMoleculeSystem)
{
  // pick move type from probabilities object
  MoveTypes moveType = selectedSystem.components[selectedComponent].mc_moves_probabilities.sample(random);

  // save old number of molecules for reference
  std::size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  switch (moveType)
  {
    case MoveTypes::Translation:
    {
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

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
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::randomTranslationMove(random, selectedSystem, selectedComponent, selectedMolecule,
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
    case MoveTypes::Rotation:
    {
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::rotationMove(random, selectedSystem, selectedComponent, selectedMolecule,
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
    case MoveTypes::RandomRotation:
    {
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::randomRotationMove(random, selectedSystem, selectedComponent, selectedMolecule,
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
    case MoveTypes::VolumeChange:
    {
      // perform move
      std::optional<RunningEnergy> energy = MC_Moves::volumeMove(random, selectedSystem);

      // accept if energy difference is not 0
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value();
      }
      break;
    }
    case MoveTypes::ReinsertionCBMC:
    {
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

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
    case MoveTypes::PartialReinsertionCBMC:
    {
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference = MC_Moves::partialReinsertionMove(
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
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference = MC_Moves::identityChangeMove(
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
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

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
    case MoveTypes::SwapCFCMC:
    case MoveTypes::SwapCBCFCMC:
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
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

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
    case MoveTypes::GibbsVolume:
    {
      std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
          MC_Moves::GibbsVolumeMove(random, selectedSystem, selectedSecondSystem);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
      break;
    }
    case MoveTypes::GibbsSwapCBMC:
    case MoveTypes::GibbsSwapCFCMC:
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

void MC_Moves::performRandomMoveEquilibration(RandomNumber &random, System &selectedSystem,
                                              System &selectedSecondSystem, std::size_t selectedComponent,
                                              std::size_t &fractionalMoleculeSystem)
{
  // pick move type from probabilities object
  MoveTypes moveType = selectedSystem.components[selectedComponent].mc_moves_probabilities.sample(random);

  // save old number of molecules for reference
  std::size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  switch (moveType)
  {
    case MoveTypes::Translation:
    {
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

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
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::randomTranslationMove(random, selectedSystem, selectedComponent, selectedMolecule,
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
    case MoveTypes::Rotation:
    {
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::rotationMove(random, selectedSystem, selectedComponent, selectedMolecule,
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
    case MoveTypes::RandomRotation:
    {
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::randomRotationMove(random, selectedSystem, selectedComponent, selectedMolecule,
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
    case MoveTypes::VolumeChange:
    {
      // perform move
      std::optional<RunningEnergy> energy = MC_Moves::volumeMove(random, selectedSystem);

      // accept if energy difference is not 0
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value();
      }
      break;
    }
    case MoveTypes::ReinsertionCBMC:
    {
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

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
    case MoveTypes::PartialReinsertionCBMC:
    {
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference = MC_Moves::partialReinsertionMove(
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
      // select molecule and only move if there are actually molecules
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        // load molecule atoms
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];

        // perform move
        std::optional<RunningEnergy> energyDifference = MC_Moves::identityChangeMove(
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
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

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
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

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
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

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
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

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
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
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
                                           std::size_t selectedComponent, std::size_t &fractionalMoleculeSystem,
                                           std::size_t currentBlock)
{
  // pick move type from probabilities object
  MoveTypes moveType = selectedSystem.components[selectedComponent].mc_moves_probabilities.sample(random);

  // save old number of molecules for reference
  std::size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  switch (moveType)
  {
    case MoveTypes::Translation:
    {
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::translationMove(random, selectedSystem, selectedComponent, selectedMolecule,
                                      selectedSystem.components, molecule, molecule_atoms);
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
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::randomTranslationMove(random, selectedSystem, selectedComponent, selectedMolecule,
                                            selectedSystem.components, molecule, molecule_atoms);
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
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::rotationMove(random, selectedSystem, selectedComponent, selectedMolecule,
                                   selectedSystem.components, molecule, molecule_atoms);
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
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::randomRotationMove(random, selectedSystem, selectedComponent, selectedMolecule,
                                         selectedSystem.components, molecule, molecule_atoms);
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
      std::optional<RunningEnergy> energy = MC_Moves::volumeMove(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value();
      }
      break;
    }
    case MoveTypes::ReinsertionCBMC:
    {
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];
        std::optional<RunningEnergy> energyDifference = MC_Moves::reinsertionMove(
            random, selectedSystem, selectedComponent, selectedMolecule, molecule, molecule_atoms);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      break;
    }
    case MoveTypes::PartialReinsertionCBMC:
    {
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];
        std::optional<RunningEnergy> energyDifference = MC_Moves::partialReinsertionMove(
            random, selectedSystem, selectedComponent, selectedMolecule, molecule, molecule_atoms);

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
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

      if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
      {
        std::span<Atom> molecule_atoms = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
        std::size_t molecule_index = selectedSystem.moleculeIndexOfComponent(selectedComponent, selectedMolecule);
        Molecule &molecule = selectedSystem.moleculeData[molecule_index];
        std::optional<RunningEnergy> energyDifference = MC_Moves::identityChangeMove(
            random, selectedSystem, selectedComponent, selectedMolecule, molecule, molecule_atoms);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      break;
    }
    case MoveTypes::Swap:
    {
      if (random.uniform() < 0.5)
      {
        const auto [energyDifference, Pacc] = MC_Moves::insertionMove(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);

        // extra time keeping for insertion / deletion split
        std::chrono::system_clock::time_point t3 = std::chrono::system_clock::now();
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType]["Insertion-Total"] += (t3 - t1);
        selectedSystem.mc_moves_cputime[moveType]["Insertion-Total"] += (t3 - t1);
      }
      else
      {
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);
        const auto [energyDifference, Pacc] =
            MC_Moves::deletionMove(random, selectedSystem, selectedComponent, selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);

        // extra time keeping for insertion / deletion split
        std::chrono::system_clock::time_point t3 = std::chrono::system_clock::now();
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType]["Deletion-Total"] += (t3 - t1);
        selectedSystem.mc_moves_cputime[moveType]["Deletion-Total"] += (t3 - t1);
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

        // extra time keeping for insertion / deletion split
        std::chrono::system_clock::time_point t3 = std::chrono::system_clock::now();
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType]["Insertion-Total"] += (t3 - t1);
        selectedSystem.mc_moves_cputime[moveType]["Insertion-Total"] += (t3 - t1);
      }
      else
      {
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

        const auto [energyDifference, Pacc] =
            MC_Moves::deletionMoveCBMC(random, selectedSystem, selectedComponent, selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);

        // extra time keeping for insertion / deletion split
        std::chrono::system_clock::time_point t3 = std::chrono::system_clock::now();
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType]["Deletion-Total"] += (t3 - t1);
        selectedSystem.mc_moves_cputime[moveType]["Deletion-Total"] += (t3 - t1);
      }
      break;
    }
    case MoveTypes::SwapCFCMC:
    {
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);

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
      std::size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(random, selectedComponent);
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
        selectedSystem.runningEnergies += energy.value().first;
        selectedSecondSystem.runningEnergies += energy.value().second;
      }
      break;
    }
    case MoveTypes::GibbsSwapCBMC:
    {
      if (random.uniform() < 0.5)
      {
        selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);

        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            MC_Moves::GibbsSwapMove_CBMC(random, selectedSystem, selectedSecondSystem, selectedComponent);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType]["Total"] += (t2 - t1);
        selectedSystem.mc_moves_cputime[moveType]["Total"] += (t2 - t1);
      }
      else
      {
        selectedSecondSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);

        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            MC_Moves::GibbsSwapMove_CBMC(random, selectedSecondSystem, selectedSystem, selectedComponent);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }
        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        selectedSecondSystem.components[selectedComponent].mc_moves_cputime[moveType]["Total"] += (t2 - t1);
        selectedSecondSystem.mc_moves_cputime[moveType]["Total"] += (t2 - t1);
      }
      break;
    }
    case MoveTypes::GibbsSwapCFCMC:
    {
      if (selectedSystem.containsTheFractionalMolecule)
      {
        selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);

        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = MC_Moves::GibbsSwapMove_CFCMC(
            random, selectedSystem, selectedSecondSystem, selectedComponent, fractionalMoleculeSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType]["Total"] += (t2 - t1);
        selectedSystem.mc_moves_cputime[moveType]["Total"] += (t2 - t1);
      }
      else if (selectedSecondSystem.containsTheFractionalMolecule)
      {
        selectedSecondSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy = MC_Moves::GibbsSwapMove_CFCMC(
            random, selectedSecondSystem, selectedSystem, selectedComponent, fractionalMoleculeSystem);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        selectedSecondSystem.components[selectedComponent].mc_moves_cputime[moveType]["Total"] += (t2 - t1);
        selectedSecondSystem.mc_moves_cputime[moveType]["Total"] += (t2 - t1);
      }
      break;
    }
    case MoveTypes::Widom:
    {
      double value = MC_Moves::WidomMove(random, selectedSystem, selectedComponent);

      selectedSystem.components[selectedComponent].averageRosenbluthWeights.addWidomSample(currentBlock, value,
                                                                                           selectedSystem.weight());
      break;
    }
    case MoveTypes::WidomCFCMC:
    {
      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC(random, selectedSystem, selectedComponent, 0, true, true);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
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
  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

  // Use lookup to skip adding count to system stats for cross system
  if (moveType != MoveTypes::GibbsSwapCBMC && moveType != MoveTypes::GibbsSwapCFCMC)
  {
    selectedSystem.mc_moves_cputime[moveType]["Total"] += (t2 - t1);
    if (componentMoves.count(moveType))
    {
      selectedSystem.components[selectedComponent].mc_moves_cputime[moveType]["Total"] += (t2 - t1);
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);
    }
    else
    {
      selectedSystem.mc_moves_statistics.addAllCounts(moveType);
    }
  }
}
