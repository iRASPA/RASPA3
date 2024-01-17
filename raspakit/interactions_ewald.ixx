export module interactions_ewald;

import <span>;
import <optional>;
import <tuple>;
import <complex>;
import <vector>;
             
import double3;
import double3x3;                                                                                                       
import atom;                                                                                                            
import running_energy;                                                                                                  
import energy_status;                                                                                                   
import simulationbox;                                                                                                   
import force_factor;                                                                                                    
import forcefield;                                                                                                      
import component;


export struct Ewald
{
  std::vector<std::complex<double>> eik_xy;
  std::vector<std::complex<double>> eik_x;
  std::vector<std::complex<double>> eik_y;
  std::vector<std::complex<double>> eik_z;
  std::vector<std::pair<std::complex<double>, std::complex<double>>> storedEik;
  std::vector<std::pair<std::complex<double>, std::complex<double>>> fixedFrameworkStoredEik;
  std::vector<std::pair<std::complex<double>, std::complex<double>>> totalEik;
  double CoulombicFourierEnergySingleIon{ 0.0 };

  void registerEwaldFourierEnergySingleIon(const ForceField &forceField, const SimulationBox &simulationBox,
                                           double3 position, double charge);

  void computeEwaldFourierRigidEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                      std::span<const Atom> frameworkAtoms, RunningEnergy& energyStatus);

  void computeEwaldFourierEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                 std::span<const Atom> atoms, RunningEnergy &energyStatus);


  ForceFactor computeEwaldFourierGradient(const ForceField &forceField, const SimulationBox &simulationBox,
                                               std::span<Atom> atomPositions);

  RunningEnergy energyDifferenceEwaldFourier(const ForceField &forceField, const SimulationBox &simulationBox,
                                             std::vector<std::pair<std::complex<double>, 
                                             std::complex<double>>> &storedWavevectors, 
                                             std::span<const Atom> newatoms, std::span<const Atom> oldatoms);

  void acceptEwaldMove(const ForceField &forceField);
};
