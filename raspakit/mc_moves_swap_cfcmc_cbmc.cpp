module;

module mc_moves;

import component;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_factor;
import energy_status;
import energy_status_inter;
import lambda;
import property_widom;
import averages;
import running_energy;
import forcefield;
import move_statistics;
import mc_moves_probabilities_particles;

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


std::optional<RunningEnergy> MC_Moves::swapMove_CFCMC_CBMC(System& system, size_t selectedComponent, size_t selectedMolecule,
         bool insertionDisabled, bool deletionDisabled)
{
    Lambda& lambda = system.components[selectedComponent].lambda;
    size_t oldBin = lambda.currentBin;
    double deltaLambda = lambda.delta;
    double oldLambda = system.components[selectedComponent].lambda.lambdaValue();
    std::make_signed_t<std::size_t> selectedNewBin = lambda.selectNewBin();

    //return std::nullopt;
    if (selectedNewBin >= std::make_signed_t<std::size_t>(lambda.numberOfBins)) // Insertion move
    {
       //return std::nullopt;
        if(insertionDisabled)
        {
          return std::nullopt;
        }

        // Steps for insertion Lambda_new = 1 + epsilon
        // ===================================================================
        // (1) Unbiased: the fractional molecule with lambda=lambda_old is made integer (lambda=1) and deltaU is computed
        // (2) Biased: a new fractional molecule is grown with lambda_new = epsilon

        double cutOffVDW = system.forceField.cutOffVDW;
        double cutOffCoulomb = system.forceField.cutOffCoulomb;

        size_t newBin = static_cast<size_t>(selectedNewBin - std::make_signed_t<std::size_t>(lambda.numberOfBins));
        double newLambda = deltaLambda * static_cast<double>(newBin);

        system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC_CBMC.counts[0] += 1;

        std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, 0);

        // make copy of old fractional molecule for reference and restoring
        std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());

        // fractional particle becomes integer (lambda=1.0)
        for (Atom& atom : fractionalMolecule) 
        { 
          atom.setScalingToInteger();
        }        

        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
        std::optional<RunningEnergy> frameworkDifference = system.computeFrameworkMoleculeEnergyDifference(fractionalMolecule, oldFractionalMolecule);
        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald += (t2 - t1);

        if(!frameworkDifference.has_value())
        {
          // reject, set fractional molecule back to old state
          std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
          return std::nullopt;
        }

        std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
        std::optional<RunningEnergy> moleculeDifference = system.computeInterMolecularEnergyDifference(fractionalMolecule, oldFractionalMolecule);
        std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
        system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald += (u2 - u1);

        if(!moleculeDifference.has_value())
        {
          // reject, set fractional molecule back to old state
          std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
          //for (Atom& atom : fractionalMolecule) 
          //{ 
          //  atom.setScaling(oldLambda); 
          //  atom.groupId = std::byte{ 1 };
          //}
          return std::nullopt;
        }


        std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
        RunningEnergy EwaldEnergyDifference = system.energyDifferenceEwaldFourier(system.storedEik, fractionalMolecule, oldFractionalMolecule);
        std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
        system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald += (v2 - v1);

        RunningEnergy energyDifference = frameworkDifference.value() + moleculeDifference.value() + EwaldEnergyDifference;

        // grow molecule with newLambda
        size_t newMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
        std::chrono::system_clock::time_point w1 = std::chrono::system_clock::now();
        std::vector<Atom> newatoms = system.components[selectedComponent].newAtoms(newLambda, system.numberOfMoleculesPerComponent[selectedComponent]);
        for (Atom& atom : newatoms)
        {
          atom.setScaling(newLambda);
          atom.groupId = std::byte{ 1 };
        }
        std::optional<ChainData> growData = system.growMoleculeSwapInsertion(cutOffVDW, cutOffCoulomb, selectedComponent, newMolecule, newLambda, newatoms);
        std::chrono::system_clock::time_point w2 = std::chrono::system_clock::now();
        system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_NonEwald += (w2 - w1);

        if (!growData)
        {
          // reject, set fractional molecule back to old state
          std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
          //for (Atom& atom : fractionalMolecule) 
          //{ 
          //  atom.setScaling(oldLambda); 
          //  atom.groupId = std::byte{ 1 };
          //}
          return std::nullopt;
        }

        system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC_CBMC.constructed[0] += 1;

        std::chrono::system_clock::time_point y1 = std::chrono::system_clock::now();
        for (Atom& atom : growData->atom)
        {
          atom.groupId = std::byte{ 1 };

        }
        RunningEnergy energyFourierDifference = system.energyDifferenceEwaldFourier(system.totalEik, std::span(growData->atom.begin(), growData->atom.end()), {});
        std::chrono::system_clock::time_point y2 = std::chrono::system_clock::now();
        system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_Ewald += (y2 - y1);

        //EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWAddEnergy(selectedComponent) - 
        //                                    system.computeTailCorrectionVDWOldEnergy();
        RunningEnergy tailEnergyDifference;
        double correctionFactorEwald = std::exp(-system.beta * (energyFourierDifference.total() + tailEnergyDifference.total()));

        double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);
        double preFactor = correctionFactorEwald * system.beta * system.components[selectedComponent].molFraction * system.pressure * system.simulationBox.volume /
            static_cast<double>(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
        double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
        if (RandomNumber::Uniform() < preFactor * (growData->RosenbluthWeight / idealGasRosenbluthWeight) * exp(-system.beta * energyDifference.total() + biasTerm))
        {
            system.acceptEwaldMove();

            system.components[selectedComponent].lambda.setCurrentBin(newBin);

            // Note: inserting invalidates iterators and spans (the vector could reallocate memory)
            system.insertMolecule(selectedComponent, growData->atom);

            
            
            // swap first and last molecule (selectedMolecule) so that molecule 0 is always the fractional molecule 
            size_t lastMoleculeId = system.numberOfMoleculesPerComponent[selectedComponent] - 1;
            std::span<Atom> lastMolecule = system.spanOfMolecule(selectedComponent, lastMoleculeId);
            fractionalMolecule = system.spanOfMolecule(selectedComponent, 0);
            //for (Atom& atom : fractionalMolecule)
            //{
            //  atom.moleculeId = static_cast<short>(lastMoleculeId);
            //}

            std::swap_ranges(fractionalMolecule.begin(), fractionalMolecule.end(), lastMolecule.begin());
            
           
            for (Atom& atom : fractionalMolecule)
            {
              //atom.setScaling(newLambda);
              atom.groupId = std::byte{ 1 };
            }
            for (Atom& atom : lastMolecule)
            {
              atom.setScalingToInteger();
              atom.groupId = std::byte{ 0 };
            }

            
            system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC_CBMC.accepted[0] += 1;

            return energyDifference + growData->energies + energyFourierDifference + tailEnergyDifference;
        };

        // Restore old lamba
        std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
        //for (Atom& atom : fractionalMolecule) 
        //{ 
        //   atom.setScaling(oldLambda); 
        //   atom.groupId = std::byte{ 1 };
        //}
        
        return std::nullopt;
    }
    else if (selectedNewBin < 0) // Deletion move
    {
       //return std::nullopt;
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
          double cutOffVDW = system.forceField.cutOffVDW;
          double cutOffCoulomb = system.forceField.cutOffCoulomb;

          selectedMolecule = system.randomIntegerMoleculeOfComponent(selectedComponent);
          
          system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC_CBMC.counts[1] += 1;
          
          std::span<Atom> fractionalMolecule = system.spanOfMolecule(selectedComponent, 0);
          std::span<Atom> newFractionalMolecule = system.spanOfMolecule(selectedComponent, selectedMolecule);
          
          // make copy of old fractional molecule for reference and restoring
          std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());
          std::vector<Atom> oldNewFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());
          
          // (1) Biased: the existing fractional molecule is retraced using CBMC with lambda=lambda_o, fractional molecule is removed.
          std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
          ChainData retraceData = system.retraceMoleculeSwapDeletion(cutOffVDW, cutOffCoulomb, selectedComponent, 0, fractionalMolecule, oldLambda, 0.0);            
          std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
          system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_NonEwald += (t2 - t1);
          
          std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
          RunningEnergy energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, {}, fractionalMolecule);
          std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
          system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_Ewald += (u2 - u1);
          
          double correctionFactorEwald = std::exp(-system.beta * energyFourierDifference.total());
          
          for (Atom &atom : fractionalMolecule) 
          { 
            atom.setScalingFullyOff();
            atom.groupId = std::byte{ 0 };
          }
          
          std::vector<Atom> savedFractionalMolecule(newFractionalMolecule.begin(), newFractionalMolecule.end());
          
          // (2) Unbiased: A new fractional molecule is chosen with lambda_new = 1 - epsilon, deltaU is computed.
          size_t newBin = static_cast<size_t>(selectedNewBin + std::make_signed_t<std::size_t>(lambda.numberOfBins));
          double newLambda = deltaLambda * static_cast<double>(newBin);
          for (Atom& atom : newFractionalMolecule) 
          { 
            atom.setScaling(newLambda); 
            atom.groupId = std::byte{ 1 };
          }
          
          std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
          std::optional<RunningEnergy> frameworkDifference = system.computeFrameworkMoleculeEnergyDifference(newFractionalMolecule, savedFractionalMolecule);
          std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
          system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapDeletionMove_CFCMC_CBMC_NonEwald += (v2 - v1);
          if(!frameworkDifference.has_value())
          {
            std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
            std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());
           
            return std::nullopt;
          }
          
          std::chrono::system_clock::time_point w1 = std::chrono::system_clock::now();
          std::optional<RunningEnergy> moleculeDifference = system.computeInterMolecularEnergyDifference(newFractionalMolecule, savedFractionalMolecule);
          std::chrono::system_clock::time_point w2 = std::chrono::system_clock::now();
          system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapDeletionMove_CFCMC_CBMC_NonEwald += (w2 - w1);
          if(!moleculeDifference.has_value())
          {
            std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
            std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());
            
            return std::nullopt;
          }
          
          std::chrono::system_clock::time_point y1 = std::chrono::system_clock::now();
          RunningEnergy EwaldEnergyDifference = system.energyDifferenceEwaldFourier(system.totalEik, newFractionalMolecule, savedFractionalMolecule);
          std::chrono::system_clock::time_point y2 = std::chrono::system_clock::now();
          system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapDeletionMove_CFCMC_CBMC_Ewald += (y2 - y1);
          
          //EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWRemoveEnergy(selectedComponent) - 
          //                                    system.computeTailCorrectionVDWOldEnergy();
          RunningEnergy tailEnergyDifference;
          RunningEnergy energyDifference = frameworkDifference.value() + moleculeDifference.value() + EwaldEnergyDifference + tailEnergyDifference;
          
          system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC_CBMC.constructed[1] += 1;
          
          double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);
          double preFactor = correctionFactorEwald * double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
              (system.beta * system.components[selectedComponent].molFraction * system.pressure * system.simulationBox.volume);
          double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
          if (RandomNumber::Uniform() < preFactor * (idealGasRosenbluthWeight / retraceData.RosenbluthWeight) * exp(-system.beta * energyDifference.total() + biasTerm))
          {
              system.acceptEwaldMove();
              system.components[selectedComponent].lambda.setCurrentBin(newBin);
          
              // Swap first and last molecule (selectedMolecule) so that molecule 0 is always the fractional molecule 
              std::swap_ranges(newFractionalMolecule.begin(), newFractionalMolecule.end(), fractionalMolecule.begin());
              for (Atom& atom : fractionalMolecule) 
              { 
                atom.moleculeId = 0; 
                atom.setScaling(newLambda);
                atom.groupId = std::byte{ 1 };
              }
          
              system.deleteMolecule(selectedComponent, selectedMolecule, newFractionalMolecule);
          
              system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC_CBMC.accepted[1] += 1;
          
          
              return energyDifference + energyFourierDifference - retraceData.energies;
          };

          // Restore old lamba
          std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
          std::copy(oldNewFractionalMolecule.begin(), oldNewFractionalMolecule.end(), newFractionalMolecule.begin());
          
          return std::nullopt;            
        }
        return std::nullopt;
    }
    else // Lambda-move
    {
      
        size_t newBin = static_cast<size_t>(selectedNewBin);
        double newLambda = deltaLambda * static_cast<double>(newBin);

        system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC_CBMC.counts[2] += 1;

        std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, 0);

        std::vector<Atom> trialPositions(molecule.begin(), molecule.end());
        std::transform(molecule.begin(), molecule.end(), trialPositions.begin(),
            [&](Atom a) { a.setScaling(newLambda); return a; });

        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
        std::optional<RunningEnergy> frameworkEnergyDifference = system.computeFrameworkMoleculeEnergyDifference(trialPositions, molecule);
        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        system.components[selectedComponent].mc_moves_probabilities.cpuTime_WidomMove_CFCMC_CBMC_NonEwald += (t2 - t1);
        system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapLambdaMove_CFCMC_CBMC_NonEwald += (t2 - t1);

        if(!frameworkEnergyDifference.has_value()) return std::nullopt;

        std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
        std::optional<RunningEnergy> interEnergyDifference = system.computeInterMolecularEnergyDifference(trialPositions, molecule);
        std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
        system.components[selectedComponent].mc_moves_probabilities.cpuTime_WidomMove_CFCMC_CBMC_NonEwald += (u2 - u1);
        system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapLambdaMove_CFCMC_CBMC_NonEwald += (u2 - u1);

        if(!interEnergyDifference.has_value()) return std::nullopt;

        std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
        RunningEnergy EwaldFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, trialPositions, molecule);
        std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
        system.components[selectedComponent].mc_moves_probabilities.cpuTime_WidomMove_CFCMC_CBMC_Ewald += (v2 - v1);
        system.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapLambdaMove_CFCMC_CBMC_Ewald += (v2 - v1);

        RunningEnergy energyDifference = frameworkEnergyDifference.value() + interEnergyDifference.value() + EwaldFourierDifference;

        system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC_CBMC.constructed[2] += 1;


        double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
        if (RandomNumber::Uniform() < std::exp(-system.beta * energyDifference.total() + biasTerm))
        {
            system.acceptEwaldMove();
            system.components[selectedComponent].mc_moves_probabilities.statistics_SwapMove_CFCMC_CBMC.accepted[2] += 1;
            
            std::copy(trialPositions.begin(), trialPositions.end(), molecule.begin());

            system.components[selectedComponent].lambda.setCurrentBin(newBin);

            return energyDifference;
        };
        return std::nullopt;
    }
}

