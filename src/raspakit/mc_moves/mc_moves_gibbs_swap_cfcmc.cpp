module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
#include <chrono>
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <utility>
#include <type_traits>
#endif

#if defined(_WIN32)
  import <cassert>;
#else
  #include <assert.h>
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
import move_statistics;
import component;
import mc_moves_probabilities_particles;
import simulationbox;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;


// All systems have a fractional molecule, only one of these is 'active', the others are switched off with 'lambda=0'.
// Implementation advantage: the number of fractional molecules per system remains constant.

// systemA contains the fractional molecule
std::optional<std::pair<RunningEnergy, RunningEnergy>> 
MC_Moves::GibbsSwapMove_CFCMC(RandomNumber &random, System& systemA, System& systemB, size_t selectedComponent, 
                              [[maybe_unused]] size_t &fractionalMoleculeSystem)
{
  PropertyLambdaProbabilityHistogram& lambdaA = systemA.components[selectedComponent].lambdaGC;
  PropertyLambdaProbabilityHistogram& lambdaB = systemB.components[selectedComponent].lambdaGC;
  size_t oldBin = lambdaA.currentBin;
  double deltaLambda = lambdaA.delta;
  double oldLambda = systemA.components[selectedComponent].lambdaGC.lambdaValue();

  double maxChange = 
    systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.maxChange[2];
  std::make_signed_t<std::size_t> selectedNewBin = lambdaA.selectNewBin(random, maxChange);

  assert(systemA.containsTheFractionalMolecule == true);
  assert(systemB.containsTheFractionalMolecule == false);

  assert(systemB.components[selectedComponent].lambdaGC.currentBin == 0);

  //assert(systemA.runningEnergies == systemA.computeTotalEnergies());
  //assert(systemB.runningEnergies == systemB.computeTotalEnergies());

  double switchValue = random.uniform();
  
  //if (selectedNewBin >= std::make_signed_t<std::size_t>(lambdaA.numberOfBins)) 
  if(switchValue < 0.25)
  {
     // Swap move:
     // Changing the fractional molecule into a whole molecule, keeping its position fixed
     // Changing a randomly selected molecule in the other simulation box into a fractional molecule (at same lambda)

     if (systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

     size_t indexFractionalMoleculeA = systemA.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
     size_t indexFractionalMoleculeB = systemB.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
     std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
     std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);

     //assert(fractionalMoleculeA.front().groupId == uint8_t{ 1 });
     //assert(fractionalMoleculeB.front().groupId == uint8_t{ 1 });

     // make copy of old fractional molecule for reference and restoring
     const std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
     const std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
     std::vector<Atom> oldFractionalMoleculeB2(fractionalMoleculeB.begin(), fractionalMoleculeB.end());

     systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.counts[0] += 1;
     systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.totalCounts[0] += 1;


     // System A: Changing the fractional molecule into a whole molecule, keeping its position fixed
     //=============================================================================================

     std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
     for (Atom& atom : fractionalMoleculeA)
     {
       atom.moleculeId = static_cast<uint32_t>(indexFractionalMoleculeA);
     }

     std::chrono::system_clock::time_point t1A = std::chrono::system_clock::now();
     std::optional<RunningEnergy> frameworkDifferenceA = 
       Interactions::computeFrameworkMoleculeEnergyDifference(systemA.forceField, systemA.simulationBox,                 
                                          systemA.spanOfFrameworkAtoms(), fractionalMoleculeA, oldFractionalMoleculeA);
     std::chrono::system_clock::time_point t2A = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (t2A - t1A);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (t2A - t1A);

     if (!frameworkDifferenceA.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       return std::nullopt;
     }

     std::chrono::system_clock::time_point u1A = std::chrono::system_clock::now();
     std::optional<RunningEnergy> moleculeDifferenceA = 
       Interactions::computeInterMolecularEnergyDifference(systemA.forceField, systemA.simulationBox,                      
                                         systemA.spanOfMoleculeAtoms(), fractionalMoleculeA, oldFractionalMoleculeA);
     std::chrono::system_clock::time_point u2A = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (u2A - u1A);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (u2A - u1A);

     if (!moleculeDifferenceA.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       return std::nullopt;
     }

     std::chrono::system_clock::time_point v1A = std::chrono::system_clock::now();
     RunningEnergy EwaldFourierDifferenceA = 
       Interactions::energyDifferenceEwaldFourier(systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy,
                                                  systemA.storedEik, systemA.totalEik,
                                                  systemA.forceField, systemA.simulationBox,
                                                  fractionalMoleculeA, oldFractionalMoleculeA);
     std::chrono::system_clock::time_point v2A = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCEwald += (v2A - v1A);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCEwald += (v2A - v1A);

     std::chrono::system_clock::time_point m1A = std::chrono::system_clock::now();
     RunningEnergy tailEnergyDifferenceA = 
       Interactions::computeInterMolecularTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                         systemA.spanOfMoleculeAtoms(), fractionalMoleculeA, oldFractionalMoleculeA) +
       Interactions::computeFrameworkMoleculeTailEnergyDifference(systemA.forceField, systemA.simulationBox,                
                                         systemA.spanOfFrameworkAtoms(), fractionalMoleculeA, oldFractionalMoleculeA);
     std::chrono::system_clock::time_point m2A = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCTail += (m2A - m1A);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCTail += (m2A - m1A);



     // step 2

     std::vector<Atom> newMolecule(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end());
     for (Atom& atom : newMolecule)
     {
       atom.setScalingToInteger();
       atom.moleculeId = static_cast<uint32_t>(systemA.numberOfMoleculesPerComponent[selectedComponent]);
     }

     std::chrono::system_clock::time_point w1A = std::chrono::system_clock::now();
     std::optional<RunningEnergy> frameworkDifferenceA2 = 
       Interactions::computeFrameworkMoleculeEnergyDifference(systemA.forceField, systemA.simulationBox,            
                                                              systemA.spanOfFrameworkAtoms(), newMolecule, {});
     std::chrono::system_clock::time_point w2A = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (w2A - w1A);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (w2A - w1A);

     if (!frameworkDifferenceA.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       return std::nullopt;
     }


     std::chrono::system_clock::time_point x1A = std::chrono::system_clock::now();
     std::optional<RunningEnergy> moleculeDifferenceA2 = 
       Interactions::computeInterMolecularEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                           systemA.spanOfMoleculeAtoms(), newMolecule, {});
     std::chrono::system_clock::time_point x2A = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (x2A - x1A);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (x2A - x1A);

     if (!moleculeDifferenceA2.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       return std::nullopt;
     }


     std::chrono::system_clock::time_point y1A = std::chrono::system_clock::now();
     RunningEnergy EwaldFourierDifferenceA2 = 
       Interactions::energyDifferenceEwaldFourier(systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy,
                                                  systemA.totalEik, systemA.totalEik,
                                                  systemA.forceField, systemA.simulationBox,
                                                  newMolecule, {});
     std::chrono::system_clock::time_point y2A = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCEwald += (y2A - y1A);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCEwald += (y2A - y1A);

     std::chrono::system_clock::time_point z1A = std::chrono::system_clock::now();
     RunningEnergy tailEnergyDifferenceA2 = 
       Interactions::computeInterMolecularTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                         systemA.spanOfMoleculeAtoms(), newMolecule, {}) +
       Interactions::computeFrameworkMoleculeTailEnergyDifference(systemA.forceField, systemA.simulationBox,            
                                         systemA.spanOfFrameworkAtoms(), newMolecule, {});
     std::chrono::system_clock::time_point z2A = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCTail += (z2A - z1A);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCTail += (z2A - z1A);

     RunningEnergy energyDifferenceA = 
       frameworkDifferenceA.value() + moleculeDifferenceA.value() + EwaldFourierDifferenceA + tailEnergyDifferenceA + 
       frameworkDifferenceA2.value() + moleculeDifferenceA2.value() + EwaldFourierDifferenceA2 + tailEnergyDifferenceA2;

     // System B: Changing a randomly selected molecule in the other simulation box into a fractional molecule 
     // (at the same lambda)
     //================================================================================================================

     size_t indexSelectedIntegerMoleculeB = systemB.randomIntegerMoleculeOfComponent(random, selectedComponent);
     std::span<Atom> selectedIntegerMoleculeB = 
       systemB.spanOfMolecule(selectedComponent, indexSelectedIntegerMoleculeB);

     // make copy of selected molecule for reference and restoring
     std::vector<Atom> oldSelectedIntegerMoleculeB(selectedIntegerMoleculeB.begin(), selectedIntegerMoleculeB.end());
     std::vector<Atom> oldSelectedIntegerMoleculeB2(selectedIntegerMoleculeB.begin(), selectedIntegerMoleculeB.end());

     for (Atom& atom : selectedIntegerMoleculeB)
     {
       atom.setScalingFullyOff();
       atom.groupId = uint8_t{ 0 };
     }
     
     std::chrono::system_clock::time_point t1B = std::chrono::system_clock::now();
     std::optional<RunningEnergy> frameworkDifferenceB = 
       Interactions::computeFrameworkMoleculeEnergyDifference(systemB.forceField, systemB.simulationBox,            
                                systemB.spanOfFrameworkAtoms(), selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);
     std::chrono::system_clock::time_point t2B = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (t2B - t1B);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (t2B - t1B);
     if (!frameworkDifferenceB.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
       std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), 
                 selectedIntegerMoleculeB.begin());
       return std::nullopt;
     }

     std::chrono::system_clock::time_point u1B = std::chrono::system_clock::now();
     std::optional<RunningEnergy> moleculeDifferenceB = 
       Interactions::computeInterMolecularEnergyDifference(systemB.forceField, systemB.simulationBox,
                                 systemB.spanOfMoleculeAtoms(), selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);
     std::chrono::system_clock::time_point u2B = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (u2B - u1B);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (u2B - u1B);
     if (!moleculeDifferenceB.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
       std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), 
                 selectedIntegerMoleculeB.begin());
       return std::nullopt;
     }

     std::chrono::system_clock::time_point v1B = std::chrono::system_clock::now();
     RunningEnergy EwaldFourierDifferenceB = 
       Interactions::energyDifferenceEwaldFourier(systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy,
                                                  systemB.storedEik, systemB.totalEik,
                                                  systemB.forceField, systemB.simulationBox,
                                                  selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);
     std::chrono::system_clock::time_point v2B = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCEwald += (v2B - v1B);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCEwald += (v2B - v1B);

     std::chrono::system_clock::time_point m1B = std::chrono::system_clock::now();
     RunningEnergy tailEnergyDifferenceB = 
       Interactions::computeInterMolecularTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                systemB.spanOfMoleculeAtoms(), selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB) +
       Interactions::computeFrameworkMoleculeTailEnergyDifference(systemB.forceField, systemB.simulationBox,            
                                systemB.spanOfFrameworkAtoms(), selectedIntegerMoleculeB, oldSelectedIntegerMoleculeB);
     std::chrono::system_clock::time_point m2B = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCTail += (m2B - m1B);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCTail += (m2B - m1B);
    

     std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
     for (Atom& atom : fractionalMoleculeB)
     {
       atom.moleculeId = static_cast<uint16_t>(indexFractionalMoleculeB);
       atom.setScaling(oldLambda);
       atom.groupId = uint8_t{ 1 };
     }

     for (Atom& atom : selectedIntegerMoleculeB)
     {
       atom.setScalingFullyOff();
       atom.groupId = uint8_t{ 0 };
       atom.position = systemA.simulationBox.randomPosition(random);
     }

     std::chrono::system_clock::time_point w1B = std::chrono::system_clock::now();
     std::optional<RunningEnergy> frameworkDifferenceB2 = 
       Interactions::computeFrameworkMoleculeEnergyDifference(systemB.forceField, systemB.simulationBox,            
                                         systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
     std::chrono::system_clock::time_point w2B = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (w2B - w1B);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (w2B - w1B);
     if (!frameworkDifferenceB2.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
       std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), 
                 selectedIntegerMoleculeB.begin());
       return std::nullopt;
     }

     std::chrono::system_clock::time_point x1B = std::chrono::system_clock::now();
     std::optional<RunningEnergy> moleculeDifferenceB2 = 
       Interactions::computeInterMolecularEnergyDifference(systemB.forceField, systemB.simulationBox,
                                         systemB.spanOfMoleculeAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
     std::chrono::system_clock::time_point x2B = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (x2B - x1B);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCNonEwald += (x2B - x1B);
     if (!moleculeDifferenceB2.has_value())
     {
       // reject, set fractional molecule back to old state
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
       std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), 
                 selectedIntegerMoleculeB.begin());
       return std::nullopt;
     }

     std::chrono::system_clock::time_point y1B = std::chrono::system_clock::now();
     RunningEnergy EwaldFourierDifferenceB2 = 
       Interactions::energyDifferenceEwaldFourier(systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy,
                                                  systemB.totalEik, systemB.totalEik,
                                                  systemB.forceField, systemB.simulationBox,
                                                  fractionalMoleculeB, oldFractionalMoleculeB);
     std::chrono::system_clock::time_point y2B = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCEwald += (y2B - y1B);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCEwald += (y2B - y1B);

     std::chrono::system_clock::time_point z1B = std::chrono::system_clock::now();
     RunningEnergy tailEnergyDifferenceB2 = 
       Interactions::computeInterMolecularTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                         systemB.spanOfMoleculeAtoms(), fractionalMoleculeB, oldFractionalMoleculeB) +
       Interactions::computeFrameworkMoleculeTailEnergyDifference(systemB.forceField, systemB.simulationBox,            
                                         systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
     std::chrono::system_clock::time_point z2B = std::chrono::system_clock::now();
     systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCTail += (z2B - z1B);
     systemA.mc_moves_cputime.GibbsSwapLambdaInterChangeMoveCFCMCTail += (z2B - z1B);


     RunningEnergy energyDifferenceB = 
       frameworkDifferenceB.value() + moleculeDifferenceB.value() + EwaldFourierDifferenceB + tailEnergyDifferenceB +
       frameworkDifferenceB2.value() + moleculeDifferenceB2.value() + EwaldFourierDifferenceB2 + tailEnergyDifferenceB2;


     double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];

     double preFactor = static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent]) / 
                        (1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent]));

     systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.constructed[0] += 1;
     systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.totalConstructed[0] += 1;

     // apply acceptance/rejection rule
     if (random.uniform() < preFactor * std::exp(-systemA.beta * 
                                          (energyDifferenceA.total() + energyDifferenceB.total()) + biasTerm))
     {
       systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.accepted[0] += 1;
       systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.totalAccepted[0] += 1;
       
       // restore
       std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
       std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
       std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), selectedIntegerMoleculeB.begin());

       std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
       std::swap(systemA.components[selectedComponent].lambdaGC.currentBin, 
                 systemB.components[selectedComponent].lambdaGC.currentBin);


       // copy the fractional molecule from B to A
       std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeA.begin());
       for (Atom& atom : fractionalMoleculeA)
       {
         atom.moleculeId = static_cast<uint16_t>(indexFractionalMoleculeA);
       }

       // make old fractional molecule integer
       std::vector<Atom> addedMolecule = 
         std::vector<Atom>(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end());
       for (Atom& atom : addedMolecule)
       {
         atom.setScalingToInteger();
       }
       systemA.insertMolecule(selectedComponent, Molecule(), addedMolecule);
       
       Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);


       systemB.deleteMolecule(selectedComponent, indexSelectedIntegerMoleculeB, selectedIntegerMoleculeB);
       
       std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), fractionalMoleculeB.begin());
       for (Atom& atom : fractionalMoleculeB)
       {
         atom.moleculeId = static_cast<uint16_t>(indexFractionalMoleculeB);
         atom.setScaling(oldLambda);
         atom.groupId = uint8_t{ 1 };
       }

       Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);

       return std::make_pair(energyDifferenceA, energyDifferenceB);
     }

     std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
     std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
     std::copy(oldSelectedIntegerMoleculeB.begin(), oldSelectedIntegerMoleculeB.end(), 
               selectedIntegerMoleculeB.begin());

     return std::nullopt;
  }
  //else if (selectedNewBin < 0)
  else  if (switchValue < 0.5)
  {
    // Move fractional molecule to the other box (from A to a random position in B)

    systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.counts[1] += 1;
    systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.totalCounts[1] += 1;

    size_t indexFractionalMoleculeA = systemA.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
    size_t indexFractionalMoleculeB = systemB.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    std::span<Atom> fractionalMoleculeB = systemB.spanOfMolecule(selectedComponent, indexFractionalMoleculeB);
    
    // make copy of old fractional molecule for reference and restoring
    std::vector<Atom> oldFractionalMoleculeA(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::vector<Atom> oldFractionalMoleculeB(fractionalMoleculeB.begin(), fractionalMoleculeB.end());
 
    // swap the active and the inactive fractional molecule
    std::swap_ranges(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), fractionalMoleculeB.begin());


    //systemA.components[selectedComponent].
    std::vector<Atom> newatoms = systemB.randomConfiguration(random, selectedComponent, fractionalMoleculeB);
    std::copy(newatoms.begin(), newatoms.end(), fractionalMoleculeB.begin());

    std::chrono::system_clock::time_point t1A = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceA = 
      Interactions::computeFrameworkMoleculeEnergyDifference(systemA.forceField, systemA.simulationBox,            
                                         systemA.spanOfFrameworkAtoms(), fractionalMoleculeA, oldFractionalMoleculeA);
    std::chrono::system_clock::time_point t2A = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCNonEwald += (t2A - t1A);
    systemA.mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCNonEwald += (t2A - t1A);
 
    if (!frameworkDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }
 
    std::chrono::system_clock::time_point u1A = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceA = 
      Interactions::computeInterMolecularEnergyDifference(systemA.forceField, systemA.simulationBox,
                                         systemA.spanOfMoleculeAtoms(), fractionalMoleculeA, oldFractionalMoleculeA);
    std::chrono::system_clock::time_point u2A = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCNonEwald += (u2A - u1A);
    systemA.mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCNonEwald += (u2A - u1A);
 
    if (!moleculeDifferenceA.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }
 
 
    std::chrono::system_clock::time_point v1A = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifferenceA = 
      Interactions::energyDifferenceEwaldFourier(systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy,
                                                 systemA.storedEik, systemA.totalEik,
                                                 systemA.forceField, systemA.simulationBox,
                                                 fractionalMoleculeA, oldFractionalMoleculeA);
    std::chrono::system_clock::time_point v2A = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCEwald += (v2A - v1A);
    systemA.mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCEwald += (v2A - v1A);

    std::chrono::system_clock::time_point w1A = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceA = 
      Interactions::computeInterMolecularTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                         systemA.spanOfMoleculeAtoms(), fractionalMoleculeA, oldFractionalMoleculeA) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(systemA.forceField, systemA.simulationBox,            
                                         systemA.spanOfFrameworkAtoms(), fractionalMoleculeA, oldFractionalMoleculeA);
    std::chrono::system_clock::time_point w2A = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCTail += (w2A - w1A);
    systemA.mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCTail += (w2A - w1A);
 
    RunningEnergy energyDifferenceA = frameworkDifferenceA.value() + moleculeDifferenceA.value() + 
                                      EwaldEnergyDifferenceA + tailEnergyDifferenceA;


    std::chrono::system_clock::time_point t1B = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkDifferenceB = 
      Interactions::computeFrameworkMoleculeEnergyDifference(systemB.forceField, systemB.simulationBox,            
                                         systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
    std::chrono::system_clock::time_point t2B = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCNonEwald += (t2B - t1B);
    systemA.mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCNonEwald += (t2B - t1B);



    if (!frameworkDifferenceB.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }

    std::chrono::system_clock::time_point u1B = std::chrono::system_clock::now();
    std::optional<RunningEnergy> moleculeDifferenceB = 
      Interactions::computeInterMolecularEnergyDifference(systemB.forceField, systemB.simulationBox,
                                         systemB.spanOfMoleculeAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
    std::chrono::system_clock::time_point u2B = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCNonEwald += (u2B - u1B);
    systemA.mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCNonEwald += (u2B - u1B);

    if (!moleculeDifferenceB.has_value())
    {
      // reject, set fractional molecule back to old state
      std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
      std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());
      return std::nullopt;
    }


    std::chrono::system_clock::time_point v1B = std::chrono::system_clock::now();
    RunningEnergy EwaldEnergyDifferenceB = 
      Interactions::energyDifferenceEwaldFourier(systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy,
                                                 systemB.storedEik, systemB.totalEik,
                                                 systemB.forceField, systemB.simulationBox,
                                                 fractionalMoleculeB, oldFractionalMoleculeB);
    std::chrono::system_clock::time_point v2B = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCEwald += (v2B - v1B);
    systemA.mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCEwald += (v2B - v1B);

    std::chrono::system_clock::time_point w1B = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceB = 
      Interactions::computeInterMolecularTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                         systemB.spanOfMoleculeAtoms(), fractionalMoleculeB, oldFractionalMoleculeB) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                         systemB.spanOfFrameworkAtoms(), fractionalMoleculeB, oldFractionalMoleculeB);
    std::chrono::system_clock::time_point w2B = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCTail += (w2B - w1B);
    systemA.mc_moves_cputime.GibbsSwapLambdaShuffleMoveCFCMCTail += (w2B - w1B);


    RunningEnergy energyDifferenceB = frameworkDifferenceB.value() + moleculeDifferenceB.value() + 
                                      EwaldEnergyDifferenceB + tailEnergyDifferenceB;

    systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.constructed[1] += 1;
    systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.totalConstructed[1] += 1;

    double biasTerm = lambdaB.biasFactor[oldBin] - lambdaA.biasFactor[oldBin];

    double preFactor = systemB.simulationBox.volume / systemA.simulationBox.volume;

    // apply acceptance/rejection rule
    if (random.uniform() < preFactor * exp(-systemA.beta * (energyDifferenceA.total() + energyDifferenceB.total()) 
                                           + biasTerm))
    {
      systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.accepted[1] += 1;
      systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.totalAccepted[1] += 1;

      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);
      Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);
    
      std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
      std::swap(systemA.components[selectedComponent].lambdaGC.currentBin, 
                systemB.components[selectedComponent].lambdaGC.currentBin);
     
      return std::make_pair(energyDifferenceA, energyDifferenceB);
    }

    // reject, set fractional molecule back to old state
    std::copy(oldFractionalMoleculeA.begin(), oldFractionalMoleculeA.end(), fractionalMoleculeA.begin());
    std::copy(oldFractionalMoleculeB.begin(), oldFractionalMoleculeB.end(), fractionalMoleculeB.begin());

    return std::nullopt;
  }
  else  // lambda move
  {
    if (selectedNewBin < 0) return  std::nullopt;
    if (selectedNewBin >= std::make_signed_t<std::size_t>(lambdaA.numberOfBins)) return  std::nullopt;

    size_t newBin = static_cast<size_t>(selectedNewBin);
    double newLambda = deltaLambda * static_cast<double>(newBin);

    
    
    systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.counts[2] += 1;
    systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.totalCounts[2] += 1;
    
    size_t indexFractionalMoleculeA = systemA.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);
    std::span<Atom> fractionalMoleculeA = systemA.spanOfMolecule(selectedComponent, indexFractionalMoleculeA);
    
    std::vector<Atom> trialPositions(fractionalMoleculeA.begin(), fractionalMoleculeA.end());
    std::transform(fractionalMoleculeA.begin(), fractionalMoleculeA.end(), trialPositions.begin(),
      [&](Atom a) { a.setScaling(newLambda); return a; });
    
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkEnergyDifference = 
      Interactions::computeFrameworkMoleculeEnergyDifference(systemA.forceField, systemA.simulationBox,            
                                         systemA.spanOfFrameworkAtoms(), trialPositions, fractionalMoleculeA);
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaChangeMoveCFCMCNonEwald += (t2 - t1);
    systemA.mc_moves_cputime.GibbsSwapLambdaChangeMoveCFCMCNonEwald += (t2 - t1);
    
    if (!frameworkEnergyDifference.has_value()) return std::nullopt;
    
    std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> interEnergyDifference = 
      Interactions::computeInterMolecularEnergyDifference(systemA.forceField, systemA.simulationBox,
                                         systemA.spanOfMoleculeAtoms(), trialPositions, fractionalMoleculeA);
    std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaChangeMoveCFCMCNonEwald += (u2 - u1);
    systemA.mc_moves_cputime.GibbsSwapLambdaChangeMoveCFCMCNonEwald += (u2 - u1);
    
    if (!interEnergyDifference.has_value()) return std::nullopt;
    
    std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
    RunningEnergy EwaldFourierDifference = 
      Interactions::energyDifferenceEwaldFourier(systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy,
                                                 systemA.storedEik, systemA.totalEik,
                                                 systemA.forceField, systemA.simulationBox,
                                                 trialPositions, fractionalMoleculeA);
    std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaChangeMoveCFCMCEwald += (v2 - v1);
    systemA.mc_moves_cputime.GibbsSwapLambdaChangeMoveCFCMCEwald += (v2 - v1);

    std::chrono::system_clock::time_point w1 = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifference = 
      Interactions::computeInterMolecularTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                         systemA.spanOfMoleculeAtoms(), trialPositions, fractionalMoleculeA) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(systemA.forceField, systemA.simulationBox,            
                                         systemA.spanOfFrameworkAtoms(), trialPositions, fractionalMoleculeA);
    std::chrono::system_clock::time_point w2 = std::chrono::system_clock::now();
    systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapLambdaChangeMoveCFCMCTail += (w2 - w1);
    systemA.mc_moves_cputime.GibbsSwapLambdaChangeMoveCFCMCTail += (w2 - w1);

    
    RunningEnergy energyDifference = frameworkEnergyDifference.value() + interEnergyDifference.value() + 
                                     EwaldFourierDifference + tailEnergyDifference;
    
    systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.constructed[2] += 1;
    systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.totalConstructed[2] += 1;
    
    double biasTerm = lambdaA.biasFactor[newBin] - lambdaA.biasFactor[oldBin];

    // apply acceptance/rejection rule
    if (random.uniform() < std::exp(-systemA.beta * energyDifference.total() + biasTerm))
    {
      Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);

      systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.accepted[2] += 1;
      systemA.components[selectedComponent].mc_moves_statistics.GibbsSwapMove_CFCMC.totalAccepted[2] += 1;
    
      std::copy(trialPositions.begin(), trialPositions.end(), fractionalMoleculeA.begin());
    
      systemA.components[selectedComponent].lambdaGC.setCurrentBin(newBin);
    
      return std::make_pair(energyDifference, RunningEnergy());
    };

    return std::nullopt;
  }
  
  return std::nullopt;
}
