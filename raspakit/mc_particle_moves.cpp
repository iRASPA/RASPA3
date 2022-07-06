module;

module mc_particle_moves;

import component;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_status;
import energy_status_inter;
import lambda;
import property_widom;
import averages;

import <complex>;
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


void MC_Particle_Moves::performRandomMove(System& selectedSystem, size_t selectedComponent)
{
    double randomNumber = RandomNumber::Uniform();

    if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityTranslationMove)
    {
	    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);

		if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
		{
			std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
			std::optional<EnergyStatus> energyDifference = translationMove(selectedSystem, selectedComponent, molecule);
			if (energyDifference)
			{
				selectedSystem.runningEnergies += *energyDifference;
			}
		}
	}
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityRotationMove)
	{
	    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);

		if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
		{
			std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
			std::optional<EnergyStatus> energyDifference = rotationMove(selectedSystem, selectedComponent, molecule);
			if (energyDifference)
			{
				selectedSystem.runningEnergies += *energyDifference;
			}
		}
	}
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityVolumeMove)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityReinsertionMove_CBMC)
	{
	    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);

		if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
		{
			std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
			std::optional<EnergyStatus> energyDifference = reinsertionMove(selectedSystem, selectedComponent, selectedMolecule, molecule);

			if (energyDifference)
			{
				selectedSystem.runningEnergies += *energyDifference;
			}
		}
	}
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityIdentityChangeMove_CBMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilitySwapMove_CBMC)
	{
		if (RandomNumber::Uniform() < 0.5)
		{
			std::optional<EnergyStatus> energyDifference = insertionMove(selectedSystem, selectedComponent);

			if (energyDifference)
			{
				selectedSystem.runningEnergies += *energyDifference;
			}
		}
		else
		{
	        size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(selectedComponent);

			std::optional<EnergyStatus> energyDifference = deletionMove(selectedSystem, selectedComponent, selectedMolecule);

			if (energyDifference)
			{
				selectedSystem.runningEnergies -= *energyDifference;
			}
		}
		
	}
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilitySwapMove_CFCMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilitySwapMove_CFCMC_CBMC)
	{
	    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);

		std::optional<EnergyStatus> energyDifference = swapMove_CFCMC_CBMC(selectedSystem, selectedComponent, selectedMolecule);
		if (energyDifference)
		{
			selectedSystem.runningEnergies += *energyDifference;
		}

	}
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityGibbsVolumeMove)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityGibbsSwapMove_CBMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityGibbsSwapMove_CFCMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityWidomMove)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityWidomMove_CFCMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityWidomMove_CFCMC_CBMC)
    {
      std::optional<EnergyStatus> energyDifference = swapMove_CFCMC_CBMC(selectedSystem, selectedComponent, 0, true, true);
      if (energyDifference)
	  {
		selectedSystem.runningEnergies += *energyDifference;
      }
    }
}

void MC_Particle_Moves::performRandomMoveProduction(System& selectedSystem, size_t selectedComponent, size_t currentBlock)
{
    double randomNumber = RandomNumber::Uniform();

    if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityTranslationMove)
    {
	    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);
		std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

		if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
		{
			std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
			std::optional<EnergyStatus> energyDifference = translationMove(selectedSystem, selectedComponent, molecule);
			if (energyDifference)
			{
				selectedSystem.runningEnergies += *energyDifference;
			}
		}
		std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
		selectedSystem.components[selectedComponent].cpuTime_TranslationMove += (t2 - t1);
	}
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityRotationMove)
	{
	    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);
		std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

		if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
		{
			std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
			std::optional<EnergyStatus> energyDifference = rotationMove(selectedSystem, selectedComponent, molecule);
			if (energyDifference)
			{
				selectedSystem.runningEnergies += *energyDifference;
			}
		}
		std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
		selectedSystem.components[selectedComponent].cpuTime_RotationMove += (t2 - t1);
	}
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityVolumeMove)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityReinsertionMove_CBMC)
	{
	    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);
		std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

		if (selectedSystem.numberOfMoleculesPerComponent[selectedComponent] > 0)
		{
			std::span<Atom> molecule = selectedSystem.spanOfMolecule(selectedComponent, selectedMolecule);
			std::optional<EnergyStatus> energyDifference = reinsertionMove(selectedSystem, selectedComponent, selectedMolecule, molecule);

			if (energyDifference)
			{
				selectedSystem.runningEnergies += *energyDifference;
			}
		}
		std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
		selectedSystem.components[selectedComponent].cpuTime_ReinsertionMove_CBMC += (t2 - t1);
	}
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityIdentityChangeMove_CBMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilitySwapMove_CBMC)
	{
		if (RandomNumber::Uniform() < 0.5)
		{
			std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

			std::optional<EnergyStatus> energyDifference = insertionMove(selectedSystem, selectedComponent);

			if (energyDifference)
			{
				selectedSystem.runningEnergies += *energyDifference;
			}
			std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
			selectedSystem.components[selectedComponent].cpuTime_SwapInsertionMove_CBMC += (t2 - t1);
		}
		else
		{
	        size_t selectedMolecule = selectedSystem.randomIntegerMoleculeOfComponent(selectedComponent);
			std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

			std::optional<EnergyStatus> energyDifference = deletionMove(selectedSystem, selectedComponent, selectedMolecule);

			if (energyDifference)
			{
				selectedSystem.runningEnergies -= *energyDifference;
			}
			
			std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
			selectedSystem.components[selectedComponent].cpuTime_SwapDeletionMove_CBMC += (t2 - t1);
		}
		
	}
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilitySwapMove_CFCMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilitySwapMove_CFCMC_CBMC)
	{
	    size_t selectedMolecule = selectedSystem.randomMoleculeOfComponent(selectedComponent);
		std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

		std::optional<EnergyStatus> energyDifference = swapMove_CFCMC_CBMC(selectedSystem, selectedComponent, selectedMolecule);
		if (energyDifference)
		{
			selectedSystem.runningEnergies += *energyDifference;
		}

		std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
		selectedSystem.components[selectedComponent].cpuTime_SwapMove_CFCMC_CBMC += (t2 - t1);
	}
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityGibbsVolumeMove)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityGibbsSwapMove_CBMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityGibbsSwapMove_CFCMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityWidomMove)
    {
	  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

	  std::optional<double> RosenbluthWeight = WidomMove(selectedSystem, selectedComponent);
      selectedSystem.components[selectedComponent].averageRosenbluthWeights.addWidomSample(currentBlock, RosenbluthWeight.value_or(0.0), selectedSystem.weight());

	  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
	  selectedSystem.components[selectedComponent].cpuTime_WidomMove_CBMC += (t2 - t1);
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityWidomMove_CFCMC)
    {
    }
    else if (randomNumber < selectedSystem.components[selectedComponent].accumulatedProbabilityWidomMove_CFCMC_CBMC)
    {
	  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

      std::optional<EnergyStatus> energyDifference = swapMove_CFCMC_CBMC(selectedSystem, selectedComponent, 0, true, true);
      if (energyDifference)
	  {
		selectedSystem.runningEnergies += *energyDifference;
      }

	  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
	  selectedSystem.components[selectedComponent].cpuTime_WidomMove_CFCMC_CBMC += (t2 - t1);
    }
}

