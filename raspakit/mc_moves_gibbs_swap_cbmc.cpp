module;

module mc_moves;

import <optional>;
import <span>;
import <chrono>;
import <vector>;
import <cmath>;
import <tuple>;

import randomnumbers;
import running_energy;
import system;
import atom;
import cbmc;
import cbmc_chain_data;
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


std::optional<std::pair<RunningEnergy, RunningEnergy>> 
MC_Moves::GibbsSwapMove_CBMC(RandomNumber &random, System& systemA, System& systemB, size_t selectedComponent)
{
  if (systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

  size_t newMoleculeIndex = systemA.numberOfMoleculesPerComponent[selectedComponent];
  systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.counts += 1;
  systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.totalCounts += 1;
  systemB.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.counts += 1;
  systemB.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.totalCounts += 1;

  double cutOffVDW = systemA.forceField.cutOffVDW;
  double cutOffCoulomb = systemA.forceField.cutOffCoulomb;
  Component::GrowType growType = systemA.components[selectedComponent].growType;

  std::vector<Atom> atoms = systemA.components[selectedComponent].recenteredCopy(1.0, 
                                                          systemA.numberOfMoleculesPerComponent[selectedComponent]);
  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  std::optional<ChainData> growData = 
    CBMC::growMoleculeSwapInsertion(random, systemA.components, systemA.forceField, systemA.simulationBox, 
                                    systemA.spanOfFrameworkAtoms(), systemA.spanOfMoleculeAtoms(), systemA.beta,
                                    growType, cutOffVDW, cutOffCoulomb, selectedComponent, newMoleculeIndex, 1.0, 
                                    atoms, systemA.numberOfTrialDirections);
  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCNonEwald += (t2 - t1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCNonEwald += (t2 - t1);
  if (!growData) return std::nullopt;

  std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

  systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.constructed += 1;
  systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.totalConstructed += 1;

  std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifferenceA = systemA.energyDifferenceEwaldFourier(systemA.storedEik, newMolecule, {});
  std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCEwald += (u2 - u1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCEwald += (u2 - u1);

  std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifferenceA =
      systemA.computeInterMolecularTailEnergyDifference(newMolecule, {}) +
      systemA.computeFrameworkMoleculeTailEnergyDifference(newMolecule, {});
  std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCTail += (v2 - v1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCTail += (v2 - v1);

  double correctionFactorEwaldA = 
    std::exp(-systemA.beta * (energyFourierDifferenceA.total() + tailEnergyDifferenceA.total()));
  
  size_t selectedMolecule = systemB.randomMoleculeOfComponent(random, selectedComponent);
  std::span<Atom> molecule = systemB.spanOfMolecule(selectedComponent, selectedMolecule);

  std::chrono::system_clock::time_point tB1 = std::chrono::system_clock::now();
  ChainData retraceData = 
    CBMC::retraceMoleculeSwapDeletion(random, systemB.components, systemB.forceField, systemB.simulationBox, 
                                      systemB.spanOfFrameworkAtoms(), systemB.spanOfMoleculeAtoms(), systemB.beta,
                                      cutOffVDW, cutOffCoulomb, selectedComponent, selectedMolecule, molecule, 
                                      1.0, 0.0, systemB.numberOfTrialDirections);
  std::chrono::system_clock::time_point tB2 = std::chrono::system_clock::now();
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCNonEwald += (tB2 - tB1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCNonEwald += (tB2 - tB1);

  std::chrono::system_clock::time_point uB1 = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifferenceB = systemB.energyDifferenceEwaldFourier(systemB.storedEik, {}, molecule);
  std::chrono::system_clock::time_point uB2 = std::chrono::system_clock::now();
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCEwald += (uB2 - uB1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCEwald += (uB2 - uB1);

  std::chrono::system_clock::time_point w1 = std::chrono::system_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifferenceB =
      systemB.computeInterMolecularTailEnergyDifference({}, molecule) +
      systemB.computeFrameworkMoleculeTailEnergyDifference({}, molecule);
  std::chrono::system_clock::time_point w2 = std::chrono::system_clock::now();
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCTail += (w2 - w1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCTail += (w2 - w1);

  systemB.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.constructed += 1;
  systemB.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.totalConstructed += 1;

  double correctionFactorEwaldB = 
    std::exp(systemB.beta * (energyFourierDifferenceB.total() + tailEnergyDifferenceB.total()));


  if (random.uniform() < (correctionFactorEwaldA * growData->RosenbluthWeight * 
                          static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent]) * 
                          systemA.simulationBox.volume) / (correctionFactorEwaldB * retraceData.RosenbluthWeight *
                          (1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent])) 
                          * systemB.simulationBox.volume))
  {
    systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.accepted += 1;
    systemA.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.totalAccepted += 1;

    systemA.acceptEwaldMove();
    systemA.insertMolecule(selectedComponent, growData->atom);

    // Debug
    //assert(system.checkMoleculeIds());

    systemB.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.accepted += 1;
    systemB.components[selectedComponent].mc_moves_probabilities.statistics_GibbsSwapMove_CBMC.totalAccepted += 1;

    systemB.acceptEwaldMove();
    systemB.deleteMolecule(selectedComponent, selectedMolecule, molecule);


    return std::make_pair(growData->energies + energyFourierDifferenceA + tailEnergyDifferenceA,
                          -(retraceData.energies - energyFourierDifferenceB - tailEnergyDifferenceB));
  };

  return std::nullopt;
}
