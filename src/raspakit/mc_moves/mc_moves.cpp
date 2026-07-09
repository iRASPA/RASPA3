module;

module mc_moves;

import std;

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
import mc_moves_force_biased_translation;
import mc_moves_force_biased_translation_all;
import mc_moves_rotation;
import mc_moves_random_rotation;
import mc_moves_reinsertion;
import mc_moves_partial_reinsertion;
import mc_moves_insertion;
import mc_moves_deletion;
import mc_moves_insertion_cbmc;
import mc_moves_deletion_cbmc;
import mc_moves_pair_insertion_cbmc;
import mc_moves_pair_deletion_cbmc;
import mc_moves_pair_swap_cfcmc;
import mc_moves_pair_swap_cfcmc_cbmc;
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
import mc_moves_gibbs_swap_cbcfcmc;
import mc_moves_gibbs_conventional_cfcmc;
import mc_moves_gibbs_conventional_cbcfcmc;
import mc_moves_gibbs_identity_change;
import mc_moves_reaction;
import mc_moves_reaction_conventional_cfcmc;
import mc_moves_reaction_conventional_cbcfcmc;
import mc_moves_reaction_cfcmc;
import mc_moves_reaction_cbcfcmc;
import mc_moves_widom;
import mc_moves_parallel_tempering_swap;
import mc_moves_hybridmc;

void MC_Moves::performRandomMoveInitialization(RandomNumber &random, System &selectedSystem,
                                               System &selectedSecondSystem, std::size_t selectedComponent,
                                               [[maybe_unused]] std::size_t &fractionalMoleculeSystem)
{
  // pick move type from probabilities object
  Move::Types moveType = selectedSystem.components[selectedComponent].mc_moves_probabilities.sample(random);

  // save old number of molecules for reference
  std::size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  switch (moveType)
  {
    case Move::Types::Translation:
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
    case Move::Types::RandomTranslation:
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
    case Move::Types::ForceBiasTranslation:
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
            MC_Moves::forceBiasTranslationMove(random, selectedSystem, selectedComponent, selectedMolecule,
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
    case Move::Types::Rotation:
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
    case Move::Types::RandomRotation:
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
    case Move::Types::VolumeChange:
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
    case Move::Types::AnisotropicVolumeChange:
    {
      std::optional<RunningEnergy> energy = MC_Moves::anisotropicVolumeMove(random, selectedSystem);

      if (energy)
      {
        selectedSystem.runningEnergies += energy.value();
      }
      break;
    }
    case Move::Types::ReinsertionCBMC:
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
    case Move::Types::PartialReinsertionCBMC:
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
    case Move::Types::IdentityChangeCBMC:
    {
      if (!selectedSystem.components[selectedComponent].identityChanges.empty())
      {
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::identityChangeMove(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }

        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      break;
    }
    case Move::Types::Swap:
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
    case Move::Types::SwapCBMC:
    case Move::Types::SwapCFCMC:
    case Move::Types::SwapCBCFCMC:
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
    case Move::Types::PairSwapCBMC:
    case Move::Types::PairSwapCFCMC:
    case Move::Types::PairSwapCBCFCMC:
    {
      if (random.uniform() < 0.5)
      {
        const auto [energyDifference, Pacc] =
            MC_Moves::pairInsertionMoveCBMC(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      else
      {
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

        const auto [energyDifference, Pacc] = MC_Moves::pairDeletionMoveCBMC(random, selectedSystem, selectedComponent,
                                                                             selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      break;
    }
    case Move::Types::PairSwap:
    {
      if (random.uniform() < 0.5)
      {
        const auto [energyDifference, Pacc] =
            MC_Moves::pairInsertionMove(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      else
      {
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

        const auto [energyDifference, Pacc] = MC_Moves::pairDeletionMove(random, selectedSystem, selectedComponent,
                                                                         selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      break;
    }
    case Move::Types::GibbsVolume:
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
    case Move::Types::GibbsSwapCBMC:
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
    case Move::Types::GibbsSwapCFCMC:
    case Move::Types::GibbsSwapCBCFCMC:
    {
      if (selectedSystem.containsTheFractionalMolecule)
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsSwapCFCMC
                ? MC_Moves::GibbsSwapMove_CFCMC(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                fractionalMoleculeSystem)
                : MC_Moves::GibbsSwapMove_CBCFCMC(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                   fractionalMoleculeSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }
      }
      else if (selectedSecondSystem.containsTheFractionalMolecule)
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsSwapCFCMC
                ? MC_Moves::GibbsSwapMove_CFCMC(random, selectedSecondSystem, selectedSystem, selectedComponent,
                                                fractionalMoleculeSystem)
                : MC_Moves::GibbsSwapMove_CBCFCMC(random, selectedSecondSystem, selectedSystem, selectedComponent,
                                                   fractionalMoleculeSystem);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }
      }
      break;
    }
    case Move::Types::GibbsConventionalCFCMC:
    case Move::Types::GibbsConventionalCBCFCMC:
    {
      if (random.uniform() < 0.5)
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsConventionalCBCFCMC
                ? MC_Moves::GibbsConventionalCBCFCMCMove(random, selectedSystem, selectedSecondSystem,
                                                             selectedComponent)
                : MC_Moves::GibbsConventionalCFCMCMove(random, selectedSystem, selectedSecondSystem,
                                                       selectedComponent);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }
      }
      else
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsConventionalCBCFCMC
                ? MC_Moves::GibbsConventionalCBCFCMCMove(random, selectedSecondSystem, selectedSystem,
                                                             selectedComponent)
                : MC_Moves::GibbsConventionalCFCMCMove(random, selectedSecondSystem, selectedSystem,
                                                       selectedComponent);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }
      }
      break;
    }
    case Move::Types::GibbsIdentityChangeCBMC:
    {
      if (!selectedSystem.components[selectedComponent].gibbsIdentityChanges.empty())
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            MC_Moves::GibbsIdentityChangeMove_CBMC(random, selectedSystem, selectedSecondSystem, selectedComponent);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }
      }
      break;
    }
    case Move::Types::Widom:
    {
      break;
    }
    case Move::Types::WidomCFCMC:
    {
      break;
    }
    case Move::Types::WidomCBCFCMC:
    {
      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, 0, true, true);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      break;
    }
    case Move::Types::ParallelTempering:
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
    case Move::Types::ForceBiasTranslationAll:
    {
      std::optional<RunningEnergy> energy = MC_Moves::forceBiasTranslationMoveAll(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value();
      }
      break;
    }
    case Move::Types::HybridMC:
    {
      std::optional<RunningEnergy> energy = MC_Moves::hybridMCMove(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value();
      }
      break;
    }
    case Move::Types::ReactionCBMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_CBMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionConventionalCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_ConventionalCFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionConventionalCBCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_ConventionalCBCFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_CFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionCBCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_CBCFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::Count:
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
  Move::Types moveType = selectedSystem.components[selectedComponent].mc_moves_probabilities.sample(random);

  // save old number of molecules for reference
  std::size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  switch (moveType)
  {
    case Move::Types::Translation:
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
    case Move::Types::RandomTranslation:
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
    case Move::Types::ForceBiasTranslation:
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
            MC_Moves::forceBiasTranslationMove(random, selectedSystem, selectedComponent, selectedMolecule,
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
    case Move::Types::Rotation:
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
    case Move::Types::RandomRotation:
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
    case Move::Types::VolumeChange:
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
    case Move::Types::AnisotropicVolumeChange:
    {
      std::optional<RunningEnergy> energy = MC_Moves::anisotropicVolumeMove(random, selectedSystem);

      if (energy)
      {
        selectedSystem.runningEnergies += energy.value();
      }
      break;
    }
    case Move::Types::ReinsertionCBMC:
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
    case Move::Types::PartialReinsertionCBMC:
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
    case Move::Types::IdentityChangeCBMC:
    {
      if (!selectedSystem.components[selectedComponent].identityChanges.empty())
      {
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::identityChangeMove(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }

        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      break;
    }
    case Move::Types::Swap:
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
    case Move::Types::SwapCBMC:
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
    case Move::Types::SwapCFCMC:
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
    case Move::Types::SwapCBCFCMC:
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
    case Move::Types::PairSwapCFCMC:
    {
      const auto [energyDifference, Pacc] =
          MC_Moves::pairSwapMove_CFCMC(random, selectedSystem, selectedComponent);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      break;
    }
    case Move::Types::PairSwapCBCFCMC:
    {
      const auto [energyDifference, Pacc] =
          MC_Moves::pairSwapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      break;
    }
    case Move::Types::PairSwapCBMC:
    {
      if (random.uniform() < 0.5)
      {
        const auto [energyDifference, Pacc] =
            MC_Moves::pairInsertionMoveCBMC(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      else
      {
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

        const auto [energyDifference, Pacc] = MC_Moves::pairDeletionMoveCBMC(random, selectedSystem, selectedComponent,
                                                                             selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      break;
    }
    case Move::Types::PairSwap:
    {
      if (random.uniform() < 0.5)
      {
        const auto [energyDifference, Pacc] =
            MC_Moves::pairInsertionMove(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      else
      {
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

        const auto [energyDifference, Pacc] = MC_Moves::pairDeletionMove(random, selectedSystem, selectedComponent,
                                                                         selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      break;
    }
    case Move::Types::GibbsVolume:
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
    case Move::Types::GibbsSwapCBMC:
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
    case Move::Types::GibbsSwapCFCMC:
    case Move::Types::GibbsSwapCBCFCMC:
    {
      if (selectedSystem.containsTheFractionalMolecule)
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsSwapCFCMC
                ? MC_Moves::GibbsSwapMove_CFCMC(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                fractionalMoleculeSystem)
                : MC_Moves::GibbsSwapMove_CBCFCMC(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                   fractionalMoleculeSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }
      }
      else if (selectedSecondSystem.containsTheFractionalMolecule)
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsSwapCFCMC
                ? MC_Moves::GibbsSwapMove_CFCMC(random, selectedSecondSystem, selectedSystem, selectedComponent,
                                                fractionalMoleculeSystem)
                : MC_Moves::GibbsSwapMove_CBCFCMC(random, selectedSecondSystem, selectedSystem, selectedComponent,
                                                   fractionalMoleculeSystem);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }
      }
      break;
    }
    case Move::Types::GibbsConventionalCFCMC:
    case Move::Types::GibbsConventionalCBCFCMC:
    {
      if (random.uniform() < 0.5)
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsConventionalCBCFCMC
                ? MC_Moves::GibbsConventionalCBCFCMCMove(random, selectedSystem, selectedSecondSystem,
                                                             selectedComponent)
                : MC_Moves::GibbsConventionalCFCMCMove(random, selectedSystem, selectedSecondSystem,
                                                       selectedComponent);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }
      }
      else
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsConventionalCBCFCMC
                ? MC_Moves::GibbsConventionalCBCFCMCMove(random, selectedSecondSystem, selectedSystem,
                                                             selectedComponent)
                : MC_Moves::GibbsConventionalCFCMCMove(random, selectedSecondSystem, selectedSystem,
                                                       selectedComponent);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }
      }
      break;
    }
    case Move::Types::GibbsIdentityChangeCBMC:
    {
      if (!selectedSystem.components[selectedComponent].gibbsIdentityChanges.empty())
      {
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            MC_Moves::GibbsIdentityChangeMove_CBMC(random, selectedSystem, selectedSecondSystem, selectedComponent);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }
      }
      break;
    }
    case Move::Types::Widom:
    {
      break;
    }
    case Move::Types::WidomCFCMC:
    {
      break;
    }
    case Move::Types::WidomCBCFCMC:
    {
      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, 0, true, true);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      break;
    }
    case Move::Types::ParallelTempering:
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
    case Move::Types::ForceBiasTranslationAll:
    {
      std::optional<RunningEnergy> energy = MC_Moves::forceBiasTranslationMoveAll(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value();
      }
      break;
    }
    case Move::Types::HybridMC:
    {
      std::optional<RunningEnergy> energy = MC_Moves::hybridMCMove(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value();
      }
      break;
    }
    case Move::Types::ReactionCBMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_CBMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionConventionalCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_ConventionalCFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionConventionalCBCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_ConventionalCBCFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_CFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionCBCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_CBCFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::Count:
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
  Move::Types moveType = selectedSystem.components[selectedComponent].mc_moves_probabilities.sample(random);

  // save old number of molecules for reference
  std::size_t oldN = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];

  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  switch (moveType)
  {
    case Move::Types::Translation:
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
    case Move::Types::RandomTranslation:
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
    case Move::Types::ForceBiasTranslation:
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
            MC_Moves::forceBiasTranslationMove(random, selectedSystem, selectedComponent, selectedMolecule,
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
    case Move::Types::Rotation:
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
    case Move::Types::RandomRotation:
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
    case Move::Types::VolumeChange:
    {
      std::optional<RunningEnergy> energy = MC_Moves::volumeMove(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value();
      }
      break;
    }
    case Move::Types::AnisotropicVolumeChange:
    {
      std::optional<RunningEnergy> energy = MC_Moves::anisotropicVolumeMove(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies += energy.value();
      }
      break;
    }
    case Move::Types::ReinsertionCBMC:
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
    case Move::Types::PartialReinsertionCBMC:
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
    case Move::Types::IdentityChangeCBMC:
    {
      if (!selectedSystem.components[selectedComponent].identityChanges.empty())
      {
        std::optional<RunningEnergy> energyDifference =
            MC_Moves::identityChangeMove(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldN);
      }
      break;
    }
    case Move::Types::Swap:
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
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::InsertionTotal] += (t3 - t1);
        selectedSystem.mc_moves_cputime[moveType][Move::Timing::InsertionTotal] += (t3 - t1);
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
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::DeletionTotal] += (t3 - t1);
        selectedSystem.mc_moves_cputime[moveType][Move::Timing::DeletionTotal] += (t3 - t1);
      }
      break;
    }
    case Move::Types::SwapCBMC:
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
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::InsertionTotal] += (t3 - t1);
        selectedSystem.mc_moves_cputime[moveType][Move::Timing::InsertionTotal] += (t3 - t1);
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
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::DeletionTotal] += (t3 - t1);
        selectedSystem.mc_moves_cputime[moveType][Move::Timing::DeletionTotal] += (t3 - t1);
      }
      break;
    }
    case Move::Types::SwapCFCMC:
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
    case Move::Types::SwapCBCFCMC:
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
    case Move::Types::PairSwapCFCMC:
    {
      const auto [energyDifference, Pacc] =
          MC_Moves::pairSwapMove_CFCMC(random, selectedSystem, selectedComponent);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      break;
    }
    case Move::Types::PairSwapCBCFCMC:
    {
      const auto [energyDifference, Pacc] =
          MC_Moves::pairSwapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      break;
    }
    case Move::Types::PairSwapCBMC:
    {
      if (random.uniform() < 0.5)
      {
        const auto [energyDifference, Pacc] =
            MC_Moves::pairInsertionMoveCBMC(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      else
      {
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

        const auto [energyDifference, Pacc] = MC_Moves::pairDeletionMoveCBMC(random, selectedSystem, selectedComponent,
                                                                             selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      break;
    }
    case Move::Types::PairSwap:
    {
      if (random.uniform() < 0.5)
      {
        const auto [energyDifference, Pacc] =
            MC_Moves::pairInsertionMove(random, selectedSystem, selectedComponent);

        if (energyDifference)
        {
          selectedSystem.runningEnergies += energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      else
      {
        std::size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(random, selectedComponent);

        const auto [energyDifference, Pacc] = MC_Moves::pairDeletionMove(random, selectedSystem, selectedComponent,
                                                                         selectedMolecule);

        if (energyDifference)
        {
          selectedSystem.runningEnergies -= energyDifference.value();
        }
        selectedSystem.tmmc.updateMatrix(Pacc, oldN);
      }
      break;
    }
    case Move::Types::GibbsVolume:
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
    case Move::Types::GibbsSwapCBMC:
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
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
        selectedSystem.mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
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
        selectedSecondSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
        selectedSecondSystem.mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
      }
      break;
    }
    case Move::Types::GibbsSwapCFCMC:
    case Move::Types::GibbsSwapCBCFCMC:
    {
      if (selectedSystem.containsTheFractionalMolecule)
      {
        selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);

        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsSwapCFCMC
                ? MC_Moves::GibbsSwapMove_CFCMC(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                fractionalMoleculeSystem)
                : MC_Moves::GibbsSwapMove_CBCFCMC(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                   fractionalMoleculeSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
        selectedSystem.mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
      }
      else if (selectedSecondSystem.containsTheFractionalMolecule)
      {
        selectedSecondSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);
        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsSwapCFCMC
                ? MC_Moves::GibbsSwapMove_CFCMC(random, selectedSecondSystem, selectedSystem, selectedComponent,
                                                fractionalMoleculeSystem)
                : MC_Moves::GibbsSwapMove_CBCFCMC(random, selectedSecondSystem, selectedSystem, selectedComponent,
                                                   fractionalMoleculeSystem);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        selectedSecondSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
        selectedSecondSystem.mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
      }
      break;
    }
    case Move::Types::GibbsConventionalCFCMC:
    case Move::Types::GibbsConventionalCBCFCMC:
    {
      if (random.uniform() < 0.5)
      {
        selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);

        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsConventionalCBCFCMC
                ? MC_Moves::GibbsConventionalCBCFCMCMove(random, selectedSystem, selectedSecondSystem,
                                                             selectedComponent)
                : MC_Moves::GibbsConventionalCFCMCMove(random, selectedSystem, selectedSecondSystem,
                                                       selectedComponent);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
        selectedSystem.mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
      }
      else
      {
        selectedSecondSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);

        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            moveType == Move::Types::GibbsConventionalCBCFCMC
                ? MC_Moves::GibbsConventionalCBCFCMCMove(random, selectedSecondSystem, selectedSystem,
                                                             selectedComponent)
                : MC_Moves::GibbsConventionalCFCMCMove(random, selectedSecondSystem, selectedSystem,
                                                       selectedComponent);
        if (energy)
        {
          selectedSecondSystem.runningEnergies += energy.value().first;
          selectedSystem.runningEnergies += energy.value().second;
        }

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        selectedSecondSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
        selectedSecondSystem.mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
      }
      break;
    }
    case Move::Types::GibbsIdentityChangeCBMC:
    {
      if (!selectedSystem.components[selectedComponent].gibbsIdentityChanges.empty())
      {
        selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);
        selectedSecondSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);

        std::optional<std::pair<RunningEnergy, RunningEnergy>> energy =
            MC_Moves::GibbsIdentityChangeMove_CBMC(random, selectedSystem, selectedSecondSystem, selectedComponent);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value().first;
          selectedSecondSystem.runningEnergies += energy.value().second;
        }

        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        selectedSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
        selectedSystem.mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
        selectedSecondSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
        selectedSecondSystem.mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
      }
      break;
    }
    case Move::Types::Widom:
    {
      double value = MC_Moves::WidomMove(random, selectedSystem, selectedComponent);

      std::size_t N = selectedSystem.numberOfIntegerMoleculesPerComponent[selectedComponent];
      double V = selectedSystem.simulationBox.volume;

      selectedSystem.components[selectedComponent].averageRosenbluthWeights.addWidomSample(currentBlock, 
           value, N, V, selectedSystem.weight());
      selectedSystem.components[selectedComponent].averageGibbsRosenbluthWeights.addWidomSample(currentBlock, 
           value, N, V, selectedSystem.weight());
      break;
    }
    case Move::Types::WidomCFCMC:
    {
      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC(random, selectedSystem, selectedComponent, 0, true, true);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      break;
    }
    case Move::Types::WidomCBCFCMC:
    {
      const auto [energyDifference, Pacc] =
          MC_Moves::swapMove_CFCMC_CBMC(random, selectedSystem, selectedComponent, 0, true, true);
      if (energyDifference)
      {
        selectedSystem.runningEnergies += energyDifference.value();
      }
      break;
    }
    case Move::Types::ParallelTempering:
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
    case Move::Types::ForceBiasTranslationAll:
    {
      std::optional<RunningEnergy> energy = MC_Moves::forceBiasTranslationMoveAll(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value();
      }
      break;
    }
    case Move::Types::HybridMC:
    {
      std::optional<RunningEnergy> energy = MC_Moves::hybridMCMove(random, selectedSystem);
      if (energy)
      {
        selectedSystem.runningEnergies = energy.value();
      }
      break;
    }
    case Move::Types::ReactionCBMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_CBMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionConventionalCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_ConventionalCFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionConventionalCBCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_ConventionalCBCFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_CFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::ReactionCBCFCMC:
    {
      if (!selectedSystem.reactions.list.empty())
      {
        std::optional<RunningEnergy> energy = MC_Moves::reactionMove_CBCFCMC(random, selectedSystem);
        if (energy)
        {
          selectedSystem.runningEnergies += energy.value();
        }
      }
      break;
    }
    case Move::Types::Count:
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
  if (moveType != Move::Types::GibbsSwapCBMC && moveType != Move::Types::GibbsSwapCFCMC &&
      moveType != Move::Types::GibbsSwapCBCFCMC && moveType != Move::Types::GibbsConventionalCFCMC &&
      moveType != Move::Types::GibbsConventionalCBCFCMC &&
      moveType != Move::Types::GibbsIdentityChangeCBMC)
  {
    selectedSystem.mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
    if (componentMoves.count(moveType))
    {
      selectedSystem.components[selectedComponent].mc_moves_cputime[moveType][Move::Timing::Total] += (t2 - t1);
      selectedSystem.components[selectedComponent].mc_moves_statistics.addAllCounts(moveType);
    }
    else
    {
      selectedSystem.mc_moves_statistics.addAllCounts(moveType);
    }
  }
}
