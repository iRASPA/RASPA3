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


void MC_Particle_Moves::performRandomMove(System& selectedSystem, size_t selectedComponent, bool production, size_t currentBlock)
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
      if(production)
      {
	    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();

	    std::optional<double> RosenbluthWeight = WidomMove(selectedSystem, selectedComponent);
        selectedSystem.components[selectedComponent].averageRosenbluthWeights.addWidomSample(currentBlock, RosenbluthWeight.value_or(0.0), selectedSystem.weight());

	    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
	    selectedSystem.components[selectedComponent].cpuTime_WidomMove += (t2 - t1);
      }
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

/*
std::optional<EnergyStatus> MC_Particle_Moves::translationMove(System & system, size_t selectedComponent, std::span<Atom> molecule)
{
	double3 displacement{};
	double3 maxDisplacement = system.components[selectedComponent].statistics_TranslationMove.maxChange;
	size_t selectedDirection = size_t(3.0 * RandomNumber::Uniform());
	displacement[selectedDirection] = maxDisplacement[selectedDirection] * 2.0 * (RandomNumber::Uniform() - 0.5);
	system.components[selectedComponent].statistics_TranslationMove.counts[selectedDirection] += 1;

	std::vector<Atom> trialPositions(molecule.size());
	std::transform(molecule.begin(), molecule.end(), trialPositions.begin(),
		[&](Atom a) { a.position += displacement; return a; });
    std::span<Atom> newMolecule{trialPositions.begin(), trialPositions.end()};

    std::optional<EnergyStatus> frameworkMolecule = system.computeFrameworkMoleculeEnergyDifference(newMolecule, molecule);
	if (!frameworkMolecule.has_value()) return std::nullopt;

    std::optional<EnergyStatus> interMolecule = system.computeInterMolecularEnergyDifference(newMolecule, molecule);
	if (!interMolecule.has_value()) return std::nullopt;

    EnergyStatus ewaldFourierEnergy = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, molecule);
    EnergyStatus energyDifference = frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy;

	system.components[selectedComponent].statistics_TranslationMove.constructed[selectedDirection] += 1;

    if (RandomNumber::Uniform() < std::exp(-system.simulationBox.Beta * energyDifference.totalEnergy))
	{
		system.components[selectedComponent].statistics_TranslationMove.accepted[selectedDirection] += 1;

        system.acceptEwaldMove();
		std::copy(trialPositions.cbegin(), trialPositions.cend(), molecule.begin());

		return energyDifference;
	};
	return std::nullopt;
}

std::optional<EnergyStatus> MC_Particle_Moves::rotationMove(System& system, size_t selectedComponent, std::span<Atom> molecule)
{
	double3 angle{};
	std::array<double3,3> axes{double3(1.0,0.0,0.0), double3(0.0,1.0,0.0) ,double3(0.0,0.0,1.0) };
	double3 maxAngle = system.components[selectedComponent].statistics_RotationMove.maxChange;
	size_t selectedDirection = size_t(3.0 * RandomNumber::Uniform());
	angle[selectedDirection] = maxAngle[selectedDirection] * 2.0 * (RandomNumber::Uniform() - 0.5);
	system.components[selectedComponent].statistics_RotationMove.counts[selectedDirection] += 1;

	size_t startingBead = system.components[selectedComponent].startingBead;
	std::vector<Atom> trialPositions(molecule.size());
    double rotationAngle = angle[selectedDirection];
    double3 rotationAxis = double3(axes[selectedDirection]);
	double3x3 rotationMatrix = double3x3(simd_quatd::fromAxisAngle(rotationAngle, rotationAxis));
	std::transform(molecule.begin(), molecule.end(), trialPositions.begin(),
			[&](Atom a) { a.position = rotationMatrix * (a.position - molecule[startingBead].position) 
                          + molecule[startingBead].position; return a; });
    std::span<Atom> newMolecule{trialPositions.begin(), trialPositions.end()};

    std::optional<EnergyStatus> frameworkMolecule = system.computeFrameworkMoleculeEnergyDifference(newMolecule, molecule);
    if (!frameworkMolecule.has_value()) return std::nullopt;

    std::optional<EnergyStatus> interMolecule = system.computeInterMolecularEnergyDifference(newMolecule, molecule);
    if (!interMolecule.has_value()) return std::nullopt;

    EnergyStatus ewaldFourierEnergy = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, molecule);
    EnergyStatus energyDifference = frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy;

	system.components[selectedComponent].statistics_RotationMove.constructed[selectedDirection] += 1;

	if (RandomNumber::Uniform() < std::exp(-system.simulationBox.Beta * energyDifference.totalEnergy))
	{
		system.components[selectedComponent].statistics_RotationMove.accepted[selectedDirection] += 1;

        system.acceptEwaldMove();
		std::copy(trialPositions.cbegin(), trialPositions.cend(), molecule.begin());

		return energyDifference;
	};
	return std::nullopt;
}

std::optional<EnergyStatus> MC_Particle_Moves::reinsertionMove(System& system, size_t selectedComponent, size_t selectedMolecule, std::span<Atom> molecule)
{
	system.components[selectedComponent].statistics_ReinsertionMove_CBMC.counts += 1;

	if (system.numberOfMoleculesPerComponent[selectedComponent] > 0)
	{
		std::optional<ChainData> growData = system.growMoleculeReinsertion(selectedComponent, selectedMolecule, molecule);

		if (!growData) return std::nullopt;

        std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

		system.components[selectedComponent].statistics_ReinsertionMove_CBMC.constructed += 1;

		ChainData retraceData = system.retraceMoleculeReinsertion(selectedComponent, selectedMolecule, molecule, growData->storedR);
        EnergyStatus energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, molecule);
        double correctionFactorFourier = std::exp(-system.simulationBox.Beta * energyFourierDifference.totalEnergy);

		if (RandomNumber::Uniform() < correctionFactorFourier * growData->RosenbluthWeight / retraceData.RosenbluthWeight)
		{
			system.components[selectedComponent].statistics_ReinsertionMove_CBMC.accepted += 1;

            system.acceptEwaldMove();
			std::copy(newMolecule.begin(), newMolecule.end(), molecule.begin());

			return (growData->energies - retraceData.energies) + energyFourierDifference;
		};
	}

	return std::nullopt;
}

std::optional<EnergyStatus> MC_Particle_Moves::insertionMove(System& system, size_t selectedComponent)
{
	size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
	system.components[selectedComponent].statistics_SwapInsertionMove_CBMC.counts += 1;
	
	std::optional<ChainData> growData = system.growMoleculeSwapInsertion(selectedComponent, selectedMolecule, 1.0);
    std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());
	if (!growData) return std::nullopt;

	system.components[selectedComponent].statistics_SwapInsertionMove_CBMC.constructed += 1;

    EnergyStatus energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, {});
    EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWAddEnergy(selectedComponent) - 
                                        system.computeTailCorrectionVDWOldEnergy();
    double correctionFactorEwald = std::exp(-system.simulationBox.Beta * (energyFourierDifference.totalEnergy + tailEnergyDifference.totalEnergy));


	double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);
	double preFactor = correctionFactorEwald * system.simulationBox.Beta * system.components[selectedComponent].molFraction * 
                       system.simulationBox.pressure * system.simulationBox.volume /
		               double(1 + system.numberOfMoleculesPerComponent[selectedComponent]);

	if (RandomNumber::Uniform() < preFactor * growData->RosenbluthWeight / idealGasRosenbluthWeight)
	{
		system.components[selectedComponent].statistics_SwapInsertionMove_CBMC.accepted += 1;

        system.acceptEwaldMove();
		system.insertMolecule(selectedComponent, growData->atom);

		// Debug
		//assert(system.checkMoleculeIds());

		return growData->energies + energyFourierDifference + tailEnergyDifference;
	};
	
	return std::nullopt;
}

std::optional<EnergyStatus> MC_Particle_Moves::deletionMove(System& system, size_t selectedComponent, size_t selectedMolecule)
{
	system.components[selectedComponent].statistics_SwapDeletionMove_CBMC.counts += 1;
	
	if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
	{
		system.components[selectedComponent].statistics_SwapDeletionMove_CBMC.constructed += 1;

		std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

		ChainData retraceData = system.retraceMoleculeSwapDeletion(selectedComponent, selectedMolecule, molecule, 1.0, 0.0);

        EnergyStatus energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, {}, molecule);
        EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWRemoveEnergy(selectedComponent) - 
                                            system.computeTailCorrectionVDWOldEnergy();
        double correctionFactorEwald = std::exp(-system.simulationBox.Beta * (energyFourierDifference.totalEnergy + tailEnergyDifference.totalEnergy));

		double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);
		double preFactor = correctionFactorEwald * double(system.numberOfMoleculesPerComponent[selectedComponent]) /
			               (system.simulationBox.Beta * system.components[selectedComponent].molFraction * 
                            system.simulationBox.pressure * system.simulationBox.volume);
		if (RandomNumber::Uniform() < preFactor * idealGasRosenbluthWeight / retraceData.RosenbluthWeight)
		{
			system.components[selectedComponent].statistics_SwapDeletionMove_CBMC.accepted += 1;

            system.acceptEwaldMove();
			system.deleteMolecule(selectedComponent, selectedMolecule, molecule);

			// Debug
			//assert(system.checkMoleculeIds());

			return retraceData.energies - energyFourierDifference - tailEnergyDifference;
		};
	}

	return std::nullopt;
}



std::optional<EnergyStatus> MC_Particle_Moves::swapMove_CFCMC_CBMC(System& system, size_t selectedComponent, size_t selectedMolecule,
         bool insertionDisabled, bool deletionDisabled)
{
	Lambda& lambda = system.components[selectedComponent].lambda;
	size_t oldBin = lambda.currentBin;
	double deltaLambda = lambda.delta;
	double oldLambda = deltaLambda * static_cast<double>(oldBin);
	std::make_signed_t<std::size_t> selectedNewBin = lambda.selectNewBin();
	
	if (selectedNewBin >= std::make_signed_t<std::size_t>(lambda.numberOfBins)) // Insertion move
	{
        if(insertionDisabled)
        {
          return std::nullopt;
        }

		// Steps for insertion Lambda_new = 1 + epsilon
		// ===================================================================
		// (1) Unbiased: the fractional molecule with lambda=lambda_old is made integer (lambda=1) and deltaU is computed
		// (2) Biased: a new fractional molecule is grown with lambda_new = epsilon

		size_t newBin = static_cast<size_t>(selectedNewBin - std::make_signed_t<std::size_t>(lambda.numberOfBins));
		double newLambda = deltaLambda * static_cast<double>(newBin);

		system.components[selectedComponent].statistics_SwapMove_CFCMC_CBMC.counts[0] += 1;

		std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, 0);

        std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());
        std::span<Atom> spanOfOldFractionalMolecule{oldFractionalMolecule.begin(), oldFractionalMolecule.end()};

		// fractional particle becomes integer (lambda=1.0)
		for (Atom& atom : fractionalMolecule) { atom.setScaling(1.0); }		

        std::optional<EnergyStatus> frameworkDifference = system.computeFrameworkMoleculeEnergyDifference(fractionalMolecule, spanOfOldFractionalMolecule);
        if(!frameworkDifference.has_value())
        {
		  for (Atom& atom : fractionalMolecule) { atom.setScaling(oldLambda); }
          return std::nullopt;
        }

        std::optional<EnergyStatus> moleculeDifference = system.computeInterMolecularEnergyDifference(fractionalMolecule, spanOfOldFractionalMolecule);
        if(!moleculeDifference.has_value())
        {
		  for (Atom& atom : fractionalMolecule) { atom.setScaling(oldLambda); }
          return std::nullopt;
        }

        EnergyStatus EwaldEnergyDifference = system.energyDifferenceEwaldFourier(system.storedEik, fractionalMolecule, spanOfOldFractionalMolecule);
        EnergyStatus energyDifference = frameworkDifference.value() + moleculeDifference.value() + EwaldEnergyDifference;

		size_t newMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
		std::optional<ChainData> growData = system.growMoleculeSwapInsertion(selectedComponent, newMolecule, newLambda);

		if (!growData)
		{
			for (Atom& atom : fractionalMolecule) { atom.setScaling(oldLambda); }
			return std::nullopt;
		}

		system.components[selectedComponent].statistics_SwapMove_CFCMC_CBMC.constructed[0] += 1;

        EnergyStatus energyFourierDifference = system.energyDifferenceEwaldFourier(system.totalEik, std::span(growData->atom.begin(), growData->atom.end()), {});
        EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWAddEnergy(selectedComponent) - 
                                            system.computeTailCorrectionVDWOldEnergy();
        double correctionFactorEwald = std::exp(-system.simulationBox.Beta * (energyFourierDifference.totalEnergy + tailEnergyDifference.totalEnergy));

		double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);
		double preFactor = correctionFactorEwald * system.simulationBox.Beta * system.components[selectedComponent].molFraction * system.simulationBox.pressure * system.simulationBox.volume /
			static_cast<double>(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
		double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
		if (RandomNumber::Uniform() < preFactor * (growData->RosenbluthWeight / idealGasRosenbluthWeight) * exp(-system.simulationBox.Beta * energyDifference.totalEnergy + biasTerm))
		{
            system.acceptEwaldMove();

			// Note: inserting invalidates iterators and spans (the vector could reallocate memory)
			system.insertMolecule(selectedComponent, growData->atom);

			system.components[selectedComponent].lambda.setCurrentBin(newBin);
			
			// swap first and last molecule (selectedMolecule) so that molecule 0 is always the fractional molecule 
			size_t lastMoleculeId = system.numberOfMoleculesPerComponent[selectedComponent] - 1;
			std::span<Atom> lastMolecule = system.spanOfMolecule(selectedComponent, lastMoleculeId);
			fractionalMolecule = system.spanOfMolecule(selectedComponent, 0);
			std::swap_ranges(fractionalMolecule.begin(), fractionalMolecule.end(), lastMolecule.begin());
			for (Atom& atom : fractionalMolecule) { atom.moleculeId = 0; }
			for (Atom& atom : lastMolecule) { atom.moleculeId = static_cast<short>(lastMoleculeId); }
			
			system.components[selectedComponent].statistics_SwapMove_CFCMC_CBMC.accepted[0] += 1;

			return growData->energies + energyDifference + energyFourierDifference + tailEnergyDifference;
		};

		// Restore old lamba
		for (Atom& atom : fractionalMolecule) { atom.setScaling(oldLambda); }
		
		return std::nullopt;
	}
	else if (selectedNewBin < 0) // Deletion move
	{
        if(deletionDisabled)
        {
          return std::nullopt;
        }
		// Steps for deletion Lambda_new = -epsilon
		// ===================================================================
		// (1) Biased: the existing fractional molecule is retraced using CBMC with lambda=lambda_o, fractional molecule is removed.
		// (2) Unbiased: A new fractional molecule is chosen with lambda_new = 1 - epsilon, deltaU is computed.

		if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
		{
			selectedMolecule = system.randomIntegerMoleculeOfComponent(selectedComponent);

			system.components[selectedComponent].statistics_SwapMove_CFCMC_CBMC.counts[1] += 1;

			std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, 0);
			std::span<Atom> newFractionalMolecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

		    // (1) Biased: the existing fractional molecule is retraced using CBMC with lambda=lambda_o, fractional molecule is removed.
			ChainData retraceData = system.retraceMoleculeSwapDeletion(selectedComponent, 0, fractionalMolecule, oldLambda, 0.0);			
            EnergyStatus energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, {}, fractionalMolecule);
            double correctionFactorEwald = std::exp(-system.simulationBox.Beta * energyFourierDifference.totalEnergy);

			for (Atom &atom : fractionalMolecule) { atom.setScaling(0.0); }

            std::vector<Atom> testFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());
            std::span<Atom> spanOfTestFractionalMolecule{testFractionalMolecule.begin(), testFractionalMolecule.end()};

		    // (2) Unbiased: A new fractional molecule is chosen with lambda_new = 1 - epsilon, deltaU is computed.
			size_t newBin = static_cast<size_t>(selectedNewBin + std::make_signed_t<std::size_t>(lambda.numberOfBins));
			double newLambda = deltaLambda * static_cast<double>(newBin);
			for (Atom& atom : newFractionalMolecule) { atom.setScaling(newLambda); }

            std::optional<EnergyStatus> frameworkDifference = system.computeFrameworkMoleculeEnergyDifference(newFractionalMolecule, spanOfTestFractionalMolecule);
            if(!frameworkDifference.has_value())
            {
			  for (Atom& atom : fractionalMolecule) { atom.setScaling(oldLambda); }
			  for (Atom& atom : newFractionalMolecule) { atom.setScaling(1.0); }
              return std::nullopt;
            }

            std::optional<EnergyStatus> moleculeDifference = system.computeInterMolecularEnergyDifference(newFractionalMolecule, spanOfTestFractionalMolecule);
            if(!moleculeDifference.has_value())
            {
			  for (Atom& atom : fractionalMolecule) { atom.setScaling(oldLambda); }
			  for (Atom& atom : newFractionalMolecule) { atom.setScaling(1.0); }
              return std::nullopt;
            }

            EnergyStatus EwaldEnergyDifference = system.energyDifferenceEwaldFourier(system.totalEik, newFractionalMolecule, spanOfTestFractionalMolecule);
            EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWRemoveEnergy(selectedComponent) - 
                                                system.computeTailCorrectionVDWOldEnergy();
            EnergyStatus energyDifference = frameworkDifference.value() + moleculeDifference.value() + EwaldEnergyDifference + tailEnergyDifference;

			system.components[selectedComponent].statistics_SwapMove_CFCMC_CBMC.constructed[1] += 1;

			double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);
			double preFactor = correctionFactorEwald * double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
				(system.simulationBox.Beta * system.components[selectedComponent].molFraction * system.simulationBox.pressure * system.simulationBox.volume);
			double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
			if (RandomNumber::Uniform() < preFactor * (idealGasRosenbluthWeight / retraceData.RosenbluthWeight) * exp(-system.simulationBox.Beta * energyDifference.totalEnergy + biasTerm))
			{
                system.acceptEwaldMove();
				system.components[selectedComponent].lambda.setCurrentBin(newBin);
			
				// Swap first and last molecule (selectedMolecule) so that molecule 0 is always the fractional molecule 
				std::swap_ranges(newFractionalMolecule.begin(), newFractionalMolecule.end(), fractionalMolecule.begin());
				for (Atom& atom : fractionalMolecule) { atom.moleculeId = 0; }
			
				system.deleteMolecule(selectedComponent, selectedMolecule, newFractionalMolecule);
			
				system.components[selectedComponent].statistics_SwapMove_CFCMC_CBMC.accepted[1] += 1;

			
				return energyDifference + energyFourierDifference - retraceData.energies;
			};

			// Restore old lamba
			for (Atom& atom : fractionalMolecule) { atom.setScaling(oldLambda); }
			for (Atom& atom : newFractionalMolecule) { atom.setScaling(1.0); }
		
			return std::nullopt;

			
		}
		return std::nullopt;
	}
	else // Lambda-move
	{
        size_t newBin = static_cast<size_t>(selectedNewBin);
	    double newLambda = deltaLambda * static_cast<double>(newBin);

		system.components[selectedComponent].statistics_SwapMove_CFCMC_CBMC.counts[2] += 1;

		std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, 0);

		std::vector<Atom> trialPositions(molecule.size());
		std::transform(molecule.begin(), molecule.end(), trialPositions.begin(),
			[&](Atom a) { a.setScaling(newLambda); return a; });
        std::span<Atom> newMolecule{trialPositions.begin(), trialPositions.end()};

        std::optional<EnergyStatus> frameworkEnergyDifference = system.computeFrameworkMoleculeEnergyDifference(newMolecule, molecule);
        if(!frameworkEnergyDifference.has_value()) return std::nullopt;

        std::optional<EnergyStatus> interEnergyDifference = system.computeInterMolecularEnergyDifference(newMolecule, molecule);
        if(!interEnergyDifference.has_value()) return std::nullopt;

        EnergyStatus EwaldFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, molecule);

		EnergyStatus energyDifference = frameworkEnergyDifference.value() + interEnergyDifference.value() + EwaldFourierDifference;

		system.components[selectedComponent].statistics_SwapMove_CFCMC_CBMC.constructed[2] += 1;


		double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
		if (RandomNumber::Uniform() < std::exp(-system.simulationBox.Beta * energyDifference.totalEnergy + biasTerm))
		{
            system.acceptEwaldMove();
			system.components[selectedComponent].statistics_SwapMove_CFCMC_CBMC.accepted[2] += 1;
			
			std::span<Atom> span = system.spanOfMolecule(selectedComponent, 0);
			std::copy(trialPositions.cbegin(), trialPositions.cend(), span.begin());

			system.components[selectedComponent].lambda.setCurrentBin(newBin);

			return energyDifference;
		};
		return std::nullopt;
	}
}

std::optional<double> MC_Particle_Moves::WidomMove(System& system, size_t selectedComponent)
{
	size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
	system.components[selectedComponent].statistics_WidomMove_CBMC.counts += 1;
	
	std::optional<ChainData> growData = system.growMoleculeSwapInsertion(selectedComponent, selectedMolecule, 1.0);
    std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());
	if (!growData) return std::nullopt;

	system.components[selectedComponent].statistics_WidomMove_CBMC.constructed += 1;

    EnergyStatus energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, {});
    EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWAddEnergy(selectedComponent) - 
                                        system.computeTailCorrectionVDWOldEnergy();
    double correctionFactorEwald = std::exp(-system.simulationBox.Beta * (energyFourierDifference.totalEnergy + tailEnergyDifference.totalEnergy));


	double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);

	return  correctionFactorEwald * growData->RosenbluthWeight / idealGasRosenbluthWeight;
}
*/
