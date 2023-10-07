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
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import forcefield;
import move_statistics;
import mc_moves_probabilities_system;

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

// mc_moves_gibbs_volume.cpp

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsVolumeMove(System &systemA, System &systemB) const
{
  systemA.mc_moves_probabilities.statistics_GibbsVolumeMove.counts += 1;
  systemA.mc_moves_probabilities.statistics_GibbsVolumeMove.totalCounts += 1;
  systemB.mc_moves_probabilities.statistics_GibbsVolumeMove.counts += 1;
  systemB.mc_moves_probabilities.statistics_GibbsVolumeMove.totalCounts += 1;

  // determine New box-volumes leaving the total volume constant
  double oldVolumeA = systemA.simulationBox.volume;
  double maxVolumeChangeA = systemA.mc_moves_probabilities.statistics_GibbsVolumeMove.maxChange;
  double oldVolumeB = systemB.simulationBox.volume;
  double totalVolume = oldVolumeA + oldVolumeB;
  double expdv = std::exp(std::log(oldVolumeA/oldVolumeB) + maxVolumeChangeA * (2.0 * RandomNumber::Uniform() - 1.0));
  double newVolumeA = expdv * totalVolume / (1.0 + expdv);
  double newVolumeB = totalVolume - newVolumeA;

  RunningEnergy oldTotalEnergyA = systemA.runningEnergies;
  double numberOfMoleculesA = static_cast<double>(std::reduce(systemA.numberOfIntegerMoleculesPerComponent.begin(), 
                        systemA.numberOfIntegerMoleculesPerComponent.end()));
  double scaleA = std::pow(newVolumeA/oldVolumeA, 1.0/3.0);
  SimulationBox newBoxA = systemA.simulationBox.scaled(scaleA);
  std::vector<Atom> newPositionsA = systemA.scaledCenterOfMassPositions(scaleA);

  RunningEnergy newTotalEnergyA;
  std::chrono::system_clock::time_point t1A = std::chrono::system_clock::now();
  systemA.computeInterMolecularEnergy(newBoxA, newPositionsA, newTotalEnergyA);
  std::chrono::system_clock::time_point t2A = std::chrono::system_clock::now();
  systemA.mc_moves_cputime.GibbsVolumeMoveNonEwald += (t2A - t1A);

  std::chrono::system_clock::time_point t3A = std::chrono::system_clock::now();
  systemA.computeEwaldFourierEnergy(newBoxA, newPositionsA, newTotalEnergyA);
  std::chrono::system_clock::time_point t4A = std::chrono::system_clock::now();
  systemA.mc_moves_cputime.GibbsVolumeMoveEwald += (t4A - t3A);

  systemA.mc_moves_probabilities.statistics_GibbsVolumeMove.constructed += 1;
  systemA.mc_moves_probabilities.statistics_GibbsVolumeMove.totalConstructed += 1;


  RunningEnergy oldTotalEnergyB = systemB.runningEnergies;
  double numberOfMoleculesB = static_cast<double>(std::reduce(systemB.numberOfIntegerMoleculesPerComponent.begin(), 
                        systemB.numberOfIntegerMoleculesPerComponent.end()));
  double scaleB = std::pow(newVolumeB/oldVolumeB, 1.0/3.0);
  SimulationBox newBoxB = systemB.simulationBox.scaled(scaleB);
  std::vector<Atom> newPositionsB = systemB.scaledCenterOfMassPositions(scaleB);


  RunningEnergy newTotalEnergyB;
  std::chrono::system_clock::time_point t1B = std::chrono::system_clock::now();
  systemB.computeInterMolecularEnergy(newBoxB, newPositionsB, newTotalEnergyB);
  std::chrono::system_clock::time_point t2B = std::chrono::system_clock::now();
  systemA.mc_moves_cputime.GibbsVolumeMoveNonEwald += (t2B - t1B);

  std::chrono::system_clock::time_point t3B = std::chrono::system_clock::now();
  systemB.computeEwaldFourierEnergy(newBoxB, newPositionsB, newTotalEnergyB);
  std::chrono::system_clock::time_point t4B = std::chrono::system_clock::now();
  systemA.mc_moves_cputime.GibbsVolumeMoveEwald += (t4B - t3B);

  systemB.mc_moves_probabilities.statistics_GibbsVolumeMove.constructed += 1;
  systemB.mc_moves_probabilities.statistics_GibbsVolumeMove.totalConstructed += 1;

  double deltaU = (newTotalEnergyA.total() - oldTotalEnergyA.total()) + (newTotalEnergyB.total() - oldTotalEnergyB.total());

  if(RandomNumber::Uniform() < std::exp(-systemA.beta * deltaU +
           (static_cast<double>(numberOfMoleculesA + 1.0) * std::log(newVolumeA/oldVolumeA))+
           (static_cast<double>(numberOfMoleculesB + 1.0) * std::log(newVolumeB/oldVolumeB)) ))
  {
    systemA.mc_moves_probabilities.statistics_GibbsVolumeMove.accepted += 1;
    systemA.mc_moves_probabilities.statistics_GibbsVolumeMove.totalAccepted += 1;

    systemA.simulationBox = newBoxA;
    std::copy(newPositionsA.begin(), newPositionsA.end(), systemA.atomPositions.begin());
    systemA.acceptEwaldMove();

    systemB.mc_moves_probabilities.statistics_GibbsVolumeMove.accepted += 1;
    systemB.mc_moves_probabilities.statistics_GibbsVolumeMove.totalAccepted += 1;

    systemB.simulationBox = newBoxB;
    std::copy(newPositionsB.begin(), newPositionsB.end(), systemB.atomPositions.begin());
    systemB.acceptEwaldMove();

    return std::make_pair(newTotalEnergyA, newTotalEnergyB);
  }

  return std::nullopt;
}
