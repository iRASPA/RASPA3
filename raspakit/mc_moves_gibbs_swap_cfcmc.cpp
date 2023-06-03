module;

module mc_moves;

import <optional>;
import <span>;
import <chrono>;
import <vector>;
import <cmath>;
import <tuple>;
import <algorithm>;
import <utility>;
import <assert.h>;

import randomnumbers;
import running_energy;
import system;
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
import move_statistics;
import component;
import mc_moves_probabilities_particles;
import simulationbox;

// All systems have a fractional molecule, but only one of these is 'active', the others are switched off with 'lambda=0'.
// Implementation advantage: the number of fractional molecules per system remains constant.

// systemA contains the fractional molecule
std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsSwapMove_CFCMC(System& systemA, System& systemB, size_t selectedComponent, size_t &fractionalMoleculeSystem)
{
  PropertyLambdaProbabilityHistogram& lambdaA = systemA.components[selectedComponent].lambdaGC;
  PropertyLambdaProbabilityHistogram& lambdaB = systemB.components[selectedComponent].lambdaGC;
  size_t oldBin = lambdaA.currentBin;
  double deltaLambda = lambdaA.delta;
  double oldLambda = systemA.components[selectedComponent].lambdaGC.lambdaValue();
  std::make_signed_t<std::size_t> selectedNewBin = lambdaA.selectNewBin();

  assert(systemA.containsTheFractionalMolecule == true);
  assert(systemB.containsTheFractionalMolecule == false);

  assert(systemB.components[selectedComponent].lambdaGC.currentBin == 0);

  double switchValue = RandomNumber::Uniform();
  
  //if (selectedNewBin >= std::make_signed_t<std::size_t>(lambdaA.numberOfBins)) 
  if(switchValue < 0.25)
  {
    //return  std::nullopt;
     // Swap move:
     // Changing the fractional molecule into a whole molecule, keeping its position fixed
     // Changing a randomly selected molecule in the other simulation box into a fractional molecule (at the same lambda)

     if (systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

     size_t indexFractionalMoleculeA = systemA.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
     size_t indexFractionalMoleculeB = systemB.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
     std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
     std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);

     // make copy of old fractional molecule for reference and restoring
     std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
     std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());

     systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CFCMC.counts[0] += 1;

     // System A: Changing the fractional molecule into a whole molecule, keeping its position fixed
     //=============================================================================================

     // fractional particle becomes integer (lambda=1.0)
     for (Atom& atom : fractionalMoleculeA)
     {
       //atom.setScalingFullyOff();
       //atom.setScalingToInteger();
       atom.setScalingFullyOn();
     }

     

     std::chrono::system_clock::time_point t1A = std::chrono::system_clock::now();
     std::optional<RunningEnergy> frameworkDifferenceA = systemA.computeFrameworkMoleculeEnergyDifference(fractionalMoleculeA, oldFractionalMoleculeA);
     std::chrono::system_clock::time_point t2A = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald += (t2A - t1A);

     if (!frameworkDifferenceA.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       return std::nullopt;
     }

     std::chrono::system_clock::time_point u1A = std::chrono::system_clock::now();
     std::optional<RunningEnergy> moleculeDifferenceA = systemA.computeInterMolecularEnergyDifference(fractionalMoleculeA, oldFractionalMoleculeA);
     std::chrono::system_clock::time_point u2A = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald += (u2A - u1A);

     if (!moleculeDifferenceA.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       return std::nullopt;
     }

     std::vector<Atom> check(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end());
     for (Atom& atom : check)
     {
       atom.setScalingToInteger();
     }
     std::optional<RunningEnergy> moleculeDifferenceA1 = systemA.computeInterMolecularEnergyDifference(check, fractionalMoleculeA);

     


     std::chrono::system_clock::time_point v1A = std::chrono::system_clock::now();
     RunningEnergy EwaldEnergyDifferenceA = systemA.energyDifferenceEwaldFourier(systemA.storedEik, fractionalMoleculeA, oldFractionalMoleculeA);
     std::chrono::system_clock::time_point v2A = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald += (v2A - v1A);

     RunningEnergy energyDifferenceA = frameworkDifferenceA.value() + moleculeDifferenceA.value() + EwaldEnergyDifferenceA;


     // System B: Changing a randomly selected molecule in the other simulation box into a fractional molecule (at the same lambda)
     //============================================================================================================================

     size_t indexSelectedIntegerMoleculeB = systemB.randomIntegerMoleculeOfComponent(selectedComponent);
     std::span<Atom> selectedIntegerMoleculeB = systemB.spanOfMolecule(selectedComponent, indexSelectedIntegerMoleculeB);

     // make copy of selected molecule for reference and restoring
     std::vector<Atom> oldSelectedIntegerMoleculeB(selectedIntegerMoleculeB.begin(), selectedIntegerMoleculeB.end());

     for (Atom& atom : selectedIntegerMoleculeB)
     {
       atom.setScaling(oldLambda);
       atom.groupId = std::byte{ 0 };
     }

     std::chrono::system_clock::time_point t1B = std::chrono::system_clock::now();
     std::optional<RunningEnergy> frameworkDifferenceB = systemB.computeFrameworkMoleculeEnergyDifference(selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);
     std::chrono::system_clock::time_point t2B = std::chrono::system_clock::now();
     systemB.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald += (t2B - t1B);

     if (!frameworkDifferenceB.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), selectedIntegerMoleculeB.begin());
       return std::nullopt;
     }

     std::chrono::system_clock::time_point u1B = std::chrono::system_clock::now();
     std::optional<RunningEnergy> moleculeDifferenceB = systemB.computeInterMolecularEnergyDifference(selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);
     std::chrono::system_clock::time_point u2B = std::chrono::system_clock::now();
     systemB.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald += (u2B - u1B);

     if (!moleculeDifferenceB.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), selectedIntegerMoleculeB.begin());
       return std::nullopt;
     }


     std::chrono::system_clock::time_point v1B = std::chrono::system_clock::now();
     RunningEnergy EwaldEnergyDifferenceB = systemB.energyDifferenceEwaldFourier(systemB.storedEik, selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);
     std::chrono::system_clock::time_point v2B = std::chrono::system_clock::now();
     systemB.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald += (v2B - v1B);

     RunningEnergy energyDifferenceB = frameworkDifferenceB.value() + moleculeDifferenceB.value() + EwaldEnergyDifferenceB;

     double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];

     // TODO: check whether it is per component or not.
     //double preFactor = static_cast<double>(systemB.numberOfIntegerMolecules() / (1.0 + static_cast<double>(systemA.numberOfIntegerMolecules())));
     double preFactor = static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] / (1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent])));

     systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CFCMC.constructed[0] += 1;

     if (RandomNumber::Uniform() < preFactor * std::exp(-systemA.beta * (energyDifferenceA.total() + energyDifferenceB.total()) + biasTerm))
     {
       systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CFCMC.accepted[0] += 1;
       
       

       // restore
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
       std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), selectedIntegerMoleculeB.begin());

       RunningEnergy oldEnergyBefore = systemA.runningEnergies;
       RunningEnergy oldEnergy = 0.5 * systemA.computeInterMolecularEnergy(systemA.forceField.cutOffVDW, systemA.forceField.cutOffCoulomb, systemA.atomPositions, -1).value();
       
       systemA.acceptEwaldMove();
       systemB.acceptEwaldMove();

       // swap the active and the inactive fractional molecule
       std::swap_ranges(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), fractionalMoleculeB.begin());
       std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
       std::swap(systemA.components[selectedComponent].lambdaGC.currentBin, systemB.components[selectedComponent].lambdaGC.currentBin);

       std::copy(selectedIntegerMoleculeB.begin(), selectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
       for (Atom& atom : fractionalMoleculeB)
       {
         atom.setScaling(oldLambda);
         atom.groupId = std::byte{ 1 };
       }

       for (Atom& atom : oldFractionalMoleculeA)
       {
         atom.setScalingToInteger();
       }

       systemA.insertMolecule(selectedComponent, oldFractionalMoleculeA);

      
       systemB.deleteMolecule(selectedComponent, indexSelectedIntegerMoleculeB, selectedIntegerMoleculeB);      

       for (Atom& atom : fractionalMoleculeA)
       {
         atom.groupId = std::byte{ 1  };
       }
       for (Atom& atom : fractionalMoleculeB)
       {
         atom.groupId = std::byte{ 1 };
       }

       //RunningEnergy oldEnergyBefore1 = systemA.runningEnergies;
       //RunningEnergy difftest = energyDifferenceA + moleculeDifferenceA1.value();
       //
       //systemA.computeTotalEnergies();
       //RunningEnergy newEnergyAfter = systemA.runningEnergies;

       RunningEnergy newEnergy = 0.5 * systemA.computeInterMolecularEnergy(systemA.forceField.cutOffVDW, systemA.forceField.cutOffCoulomb, systemA.atomPositions, -1).value();

       RunningEnergy diff = newEnergy - oldEnergy;

       return std::make_pair(energyDifferenceA, energyDifferenceB);
     }

     std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
     std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
     std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), selectedIntegerMoleculeB.begin());

     return std::nullopt;
  }
  //else if (selectedNewBin < 0) // Move fractional molecule to the other box
  else  if (switchValue < 0.5)
  {
    
   // return std::nullopt;

    systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CFCMC.counts[1] += 1;

    size_t indexFractionalMoleculeA = systemA.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
    size_t indexFractionalMoleculeB = systemB.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);
    
    // make copy of old fractional molecule for reference and restoring
    std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
 
    // swap the active and the inactive fractional molecule
    std::swap_ranges(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), fractionalMoleculeB.begin());
 
    std::chrono::system_clock::time_point t1A = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceA = systemA.computeFrameworkMoleculeEnergyDifference(fractionalMoleculeA, oldFractionalMoleculeA);
    std::chrono::system_clock::time_point t2A = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald += (t2A - t1A);
 
    if (!frameworkDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }
 
    std::chrono::system_clock::time_point u1A = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceA = systemA.computeInterMolecularEnergyDifference(fractionalMoleculeA, oldFractionalMoleculeA);
    std::chrono::system_clock::time_point u2A = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald += (u2A - u1A);
 
    if (!moleculeDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }
 
 
    std::chrono::system_clock::time_point v1A = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifferenceA = systemA.energyDifferenceEwaldFourier(systemA.storedEik, fractionalMoleculeA, oldFractionalMoleculeA);
    std::chrono::system_clock::time_point v2A = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald += (v2A - v1A);
 
    RunningEnergy energyDifferenceA = frameworkDifferenceA.value() + moleculeDifferenceA.value() + EwaldEnergyDifferenceA;


    std::chrono::system_clock::time_point t1B = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceB = systemB.computeFrameworkMoleculeEnergyDifference(fractionalMoleculeB, oldFractionalMoleculeB);
    std::chrono::system_clock::time_point t2B = std::chrono::system_clock::now();
    systemB.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald += (t2B - t1B);

    if (!frameworkDifferenceB.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    std::chrono::system_clock::time_point u1B = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceB = systemB.computeInterMolecularEnergyDifference(fractionalMoleculeB, oldFractionalMoleculeB);
    std::chrono::system_clock::time_point u2B = std::chrono::system_clock::now();
    systemB.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald += (u2B - u1B);

    if (!moleculeDifferenceB.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }


    std::chrono::system_clock::time_point v1B = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifferenceB = systemB.energyDifferenceEwaldFourier(systemB.storedEik, fractionalMoleculeB, oldFractionalMoleculeB);
    std::chrono::system_clock::time_point v2B = std::chrono::system_clock::now();
    systemB.components[selectedComponent].mc_moves_probabilities.cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald += (v2B - v1B);

    RunningEnergy energyDifferenceB = frameworkDifferenceB.value() + moleculeDifferenceB.value() + EwaldEnergyDifferenceB;

    systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CFCMC.constructed[1] += 1;

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];

    double preFactor = systemB.simulationBox.volume / systemA.simulationBox.volume;

    if (RandomNumber::Uniform() < preFactor * exp(-systemA.beta * (energyDifferenceA.total() + energyDifferenceB.total()) + biasTerm))
    {
      systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CFCMC.accepted[1] += 1;

      systemA.acceptEwaldMove();
      systemB.acceptEwaldMove();
    
      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(systemA.components[selectedComponent].lambdaGC.currentBin, systemB.components[selectedComponent].lambdaGC.currentBin);
     
      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    // reject, set fractional molecule back to old state
    std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());

    return std::nullopt;
  }
  else  // lambda move
  {
    //return  std::nullopt;

    if (selectedNewBin < 0) return  std::nullopt;
    if (selectedNewBin >= std::make_signed_t<std::size_t>(lambdaA.numberOfBins)) return  std::nullopt;

    size_t newBin = static_cast<size_t>(selectedNewBin);
    double newLambda = deltaLambda * static_cast<double>(newBin);

    
    
    systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CFCMC.counts[2] += 1;
    
    size_t indexFractionalMoleculeA = systemA.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    
    std::vector<Atom> trialPositions(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::transform(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), trialPositions.begin(),
      [&](Atom a) { a.setScaling(newLambda); return a; });
    
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkEnergyDifference = systemA.computeFrameworkMoleculeEnergyDifference(trialPositions, fractionalMoleculeA);
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_WidomMove_CFCMC_NonEwald += (t2 - t1);
    systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_GibbsSwapLambdaMove_CFCMC_NonEwald += (t2 - t1);
    
    if (!frameworkEnergyDifference.has_value()) return std::nullopt;
    
    std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> interEnergyDifference = systemA.computeInterMolecularEnergyDifference(trialPositions, fractionalMoleculeA);
    std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_WidomMove_CFCMC_NonEwald += (u2 - u1);
    systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_GibbsSwapLambdaMove_CFCMC_NonEwald += (u2 - u1);
    
    if (!interEnergyDifference.has_value()) return std::nullopt;
    
    std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifference = systemA.energyDifferenceEwaldFourier(systemA.storedEik, trialPositions, fractionalMoleculeA);
    std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_WidomMove_CFCMC_Ewald += (v2 - v1);
    systemA.components[selectedComponent].mc_moves_probabilities.cpuTime_GibbsSwapLambdaMove_CFCMC_Ewald += (v2 - v1);
    
    RunningEnergy energyDifference = frameworkEnergyDifference.value() + interEnergyDifference.value() + EwaldFourierDifference;
    
    systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CFCMC.constructed[2] += 1;
    
    
    double biasTerm = lambdaA.biasFactor[newBin] - lambdaA.biasFactor[oldBin];
    if (RandomNumber::Uniform() < std::exp(-systemA.beta * energyDifference.total() + biasTerm))
    {
      systemA.acceptEwaldMove();
      systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CFCMC.accepted[2] += 1;
    
      std::copy(trialPositions.begin(), trialPositions.end(), fractionalMoleculeA.begin());
    
      systemA.components[selectedComponent].lambdaGC.setCurrentBin(newBin);
    
      return std::make_pair(energyDifference, RunningEnergy());
    };

    return std::nullopt;
  }
  
  return std::nullopt;
}