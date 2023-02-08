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
import lambda;
import property_widom;
import averages;
import forcefield;

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


std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsVolumeMove(System &systemA, System &systemB) const
{
  systemA.statistics_VolumeMove.counts += 1;
  systemB.statistics_VolumeMove.counts += 1;

  // determine New box-volumes leaving the total volume constant
  double oldVolumeA = systemA.simulationBox.volume;
  double maxVolumeChangeA = systemA.statistics_GibbsVolumeMove.maxChange;
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
  systemA.computeInterMolecularEnergy(newBoxA, std::span(newPositionsA.begin(), newPositionsA.end()), newTotalEnergyA);
  std::chrono::system_clock::time_point t2A = std::chrono::system_clock::now();
  systemA.cpuTime_GibbsVolumeMove_NonEwald += (t2A - t1A);

  systemA.statistics_GibbsVolumeMove.constructed += 1;


  RunningEnergy oldTotalEnergyB = systemB.runningEnergies;
  double numberOfMoleculesB = static_cast<double>(std::reduce(systemB.numberOfIntegerMoleculesPerComponent.begin(), 
                        systemB.numberOfIntegerMoleculesPerComponent.end()));
  double scaleB = std::pow(newVolumeB/oldVolumeB, 1.0/3.0);
  SimulationBox newBoxB = systemB.simulationBox.scaled(scaleB);
  std::vector<Atom> newPositionsB = systemB.scaledCenterOfMassPositions(scaleB);
  RunningEnergy newTotalEnergyB;
  std::chrono::system_clock::time_point t1B = std::chrono::system_clock::now();
  systemB.computeInterMolecularEnergy(newBoxB, std::span(newPositionsB.begin(), newPositionsB.end()), newTotalEnergyB);
  std::chrono::system_clock::time_point t2B = std::chrono::system_clock::now();
  systemB.cpuTime_GibbsVolumeMove_NonEwald += (t2B - t1B);

  systemB.statistics_GibbsVolumeMove.constructed += 1;

  if(RandomNumber::Uniform() < std::exp(-systemA.Beta * (
           ((newTotalEnergyA.total() - oldTotalEnergyA.total()) + (newTotalEnergyB.total() - oldTotalEnergyB.total())) +
           ((numberOfMoleculesA + 1.0) * std::log(newVolumeA/oldVolumeA))+
           ((numberOfMoleculesB + 1.0) * std::log(newVolumeB/oldVolumeB)) )))
  {
    systemA.statistics_GibbsVolumeMove.accepted += 1;

    systemA.simulationBox = newBoxA;
    std::copy(newPositionsA.begin(), newPositionsA.end(), systemA.atomPositions.begin());

    systemB.statistics_GibbsVolumeMove.accepted += 1;

    systemB.simulationBox = newBoxB;
    std::copy(newPositionsB.begin(), newPositionsB.end(), systemB.atomPositions.begin());

    return std::make_pair(newTotalEnergyA, newTotalEnergyB);
  }

  return std::nullopt;
}
