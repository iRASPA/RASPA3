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


std::optional<RunningEnergy> MC_Moves::volumeMove([[maybe_unused]] System &system) const
{
  system.statistics_VolumeMove.counts += 1;

  RunningEnergy oldTotalEnergy = system.runningEnergies;
  double numberOfMolecules = static_cast<double>(std::reduce(system.numberOfIntegerMoleculesPerComponent.begin(), 
                        system.numberOfIntegerMoleculesPerComponent.end()));
  double oldVolume = system.simulationBox.volume;
  double maxVolumeChange = system.statistics_VolumeMove.maxChange;
  double newVolume = std::exp(std::log(oldVolume) + maxVolumeChange * (2.0 * RandomNumber::Uniform() - 1.0));
  double scale = std::pow(newVolume/oldVolume, 1.0/3.0);

  SimulationBox newBox = system.simulationBox.scaled(scale);
  std::vector<Atom> newPositions = system.scaledCenterOfMassPositions(scale);

  RunningEnergy newTotalEnergy;
  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  system.computeInterMolecularEnergy(newBox, std::span(newPositions.begin(), newPositions.end()), newTotalEnergy);
  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
  system.cpuTime_VolumeMove_NonEwald += (t2 - t1);

  system.statistics_VolumeMove.constructed += 1;

  if(RandomNumber::Uniform() < std::exp((numberOfMolecules + 1.0) * std::log(newVolume/oldVolume)
        - (system.pressure * (newVolume - oldVolume)+ (newTotalEnergy.total() - oldTotalEnergy.total())) * system.Beta))
  {
    system.statistics_VolumeMove.accepted += 1;

    system.simulationBox = newBox;
    std::copy(newPositions.begin(), newPositions.end(), system.atomPositions.begin());
    return newTotalEnergy;
  }

  return std::nullopt;
}