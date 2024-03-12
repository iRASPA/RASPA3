module;

module mc_moves_gibbs_volume;

import <complex>;
import <vector>;
import <array>;
import <tuple>;
import <optional>;
import <span>;
import <optional>;
import <tuple>;
import <algorithm>;
import <numeric>;
import <chrono>;
import <cmath>;
import <iostream>;
import <iomanip>;

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
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import forcefield;
import move_statistics;
import mc_moves_probabilities_system;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;


std::optional<std::pair<RunningEnergy, RunningEnergy>> 
MC_Moves::GibbsVolumeMove(RandomNumber &random, System &systemA, System &systemB)
{
  std::chrono::system_clock::time_point time_begin, time_end;

  systemA.mc_moves_statistics.GibbsVolumeMove.counts += 1;
  systemA.mc_moves_statistics.GibbsVolumeMove.totalCounts += 1;
  systemB.mc_moves_statistics.GibbsVolumeMove.counts += 1;
  systemB.mc_moves_statistics.GibbsVolumeMove.totalCounts += 1;

  // determine New box-volumes leaving the total volume constant
  double oldVolumeA = systemA.simulationBox.volume;
  double maxVolumeChangeA = systemA.mc_moves_statistics.GibbsVolumeMove.maxChange;
  double oldVolumeB = systemB.simulationBox.volume;
  double totalVolume = oldVolumeA + oldVolumeB;
  double expdv = std::exp(std::log(oldVolumeA/oldVolumeB) + maxVolumeChangeA * (2.0 * random.uniform() - 1.0));
  double newVolumeA = expdv * totalVolume / (1.0 + expdv);
  double newVolumeB = totalVolume - newVolumeA;

  RunningEnergy oldTotalEnergyA = systemA.runningEnergies;
  double numberOfMoleculesA = static_cast<double>(std::reduce(systemA.numberOfIntegerMoleculesPerComponent.begin(), 
                                                              systemA.numberOfIntegerMoleculesPerComponent.end()));
  double scaleA = std::pow(newVolumeA/oldVolumeA, 1.0/3.0);
  SimulationBox newBoxA = systemA.simulationBox.scaled(scaleA);
  std::vector<Atom> newPositionsA = systemA.scaledCenterOfMassPositions(scaleA);

  RunningEnergy newTotalEnergyA;
  time_begin = std::chrono::system_clock::now();
  Interactions::computeInterMolecularEnergy(systemA.forceField, newBoxA, newPositionsA, newTotalEnergyA);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime.GibbsVolumeMoveNonEwald += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  Interactions::computeInterMolecularTailEnergy(systemA.forceField, newBoxA, newPositionsA, newTotalEnergyA);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime.GibbsVolumeMoveTail += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  Interactions::computeEwaldFourierEnergy(systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy,
                                          systemA.fixedFrameworkStoredEik, systemA.totalEik,
                                          systemA.forceField, newBoxA,
                                          systemA.components, systemA.numberOfMoleculesPerComponent,
                                          newPositionsA, newPositionsA, newTotalEnergyA);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime.GibbsVolumeMoveEwald += (time_end - time_begin);

  systemA.mc_moves_statistics.GibbsVolumeMove.constructed += 1;
  systemA.mc_moves_statistics.GibbsVolumeMove.totalConstructed += 1;


  RunningEnergy oldTotalEnergyB = systemB.runningEnergies;
  double numberOfMoleculesB = static_cast<double>(std::reduce(systemB.numberOfIntegerMoleculesPerComponent.begin(), 
                        systemB.numberOfIntegerMoleculesPerComponent.end()));
  double scaleB = std::pow(newVolumeB/oldVolumeB, 1.0/3.0);
  SimulationBox newBoxB = systemB.simulationBox.scaled(scaleB);
  std::vector<Atom> newPositionsB = systemB.scaledCenterOfMassPositions(scaleB);


  RunningEnergy newTotalEnergyB;
  time_begin = std::chrono::system_clock::now();
  Interactions::computeInterMolecularEnergy(systemB.forceField, newBoxB, newPositionsB, newTotalEnergyB);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime.GibbsVolumeMoveNonEwald += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  Interactions::computeInterMolecularTailEnergy(systemB.forceField, newBoxB, newPositionsB, newTotalEnergyB);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime.GibbsVolumeMoveTail += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  Interactions::computeEwaldFourierEnergy(systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy,
                                          systemB.fixedFrameworkStoredEik, systemB.totalEik,
                                          systemB.forceField, newBoxB,
                                          systemB.components, systemB.numberOfMoleculesPerComponent,
                                          newPositionsB, newPositionsB, newTotalEnergyB);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime.GibbsVolumeMoveEwald += (time_end - time_begin);

  systemB.mc_moves_statistics.GibbsVolumeMove.constructed += 1;
  systemB.mc_moves_statistics.GibbsVolumeMove.totalConstructed += 1;

  double deltaU = (newTotalEnergyA.total() - oldTotalEnergyA.total()) + 
                  (newTotalEnergyB.total() - oldTotalEnergyB.total());

  // apply acceptance/rejection rule
  if(random.uniform() < std::exp(-systemA.beta * deltaU +
           (static_cast<double>(numberOfMoleculesA + 1.0) * std::log(newVolumeA/oldVolumeA))+
           (static_cast<double>(numberOfMoleculesB + 1.0) * std::log(newVolumeB/oldVolumeB)) ))
  {
    systemA.mc_moves_statistics.GibbsVolumeMove.accepted += 1;
    systemA.mc_moves_statistics.GibbsVolumeMove.totalAccepted += 1;

    systemA.simulationBox = newBoxA;
    std::copy(newPositionsA.begin(), newPositionsA.end(), systemA.atomPositions.begin());
    Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);

    systemB.mc_moves_statistics.GibbsVolumeMove.accepted += 1;
    systemB.mc_moves_statistics.GibbsVolumeMove.totalAccepted += 1;

    systemB.simulationBox = newBoxB;
    std::copy(newPositionsB.begin(), newPositionsB.end(), systemB.atomPositions.begin());
    Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);

    return std::make_pair(newTotalEnergyA, newTotalEnergyB);
  }

  return std::nullopt;
}
