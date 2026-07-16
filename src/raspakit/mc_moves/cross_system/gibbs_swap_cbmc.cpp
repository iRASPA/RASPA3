module;

module mc_moves_gibbs_swap_cbmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;
import atom;
import cbmc;
import cbmc_chain_data;
import cbmc_interactions;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import running_energy;
import forcefield;
import component;
import simulationbox;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;

namespace
{
class DualTMMCTrial
{
 public:
  DualTMMCTrial(System& systemA, System& systemB, std::size_t selectedComponent)
      : systemA_(systemA),
        systemB_(systemB),
        oldNA_(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent]),
        oldNB_(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent])
  {
  }

  ~DualTMMCTrial()
  {
    if (!recorded_)
    {
      systemA_.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldNA_);
      systemB_.tmmc.updateMatrix(double3(0.0, 1.0, 0.0), oldNB_);
    }
  }

  bool transferIsInBounds() const
  {
    return (!systemA_.tmmc.doTMMC || !systemA_.tmmc.rejectOutOfBound || oldNA_ < systemA_.tmmc.maxMacrostate) &&
           (!systemB_.tmmc.doTMMC || !systemB_.tmmc.rejectOutOfBound ||
            (oldNB_ > 0 && oldNB_ - 1 >= systemB_.tmmc.minMacrostate));
  }

  double biasFactor() const
  {
    return systemA_.tmmc.biasFactor(oldNA_ + 1, oldNA_) * systemB_.tmmc.biasFactor(oldNB_ - 1, oldNB_);
  }

  void recordTransfer(double physicalAcceptance)
  {
    systemA_.tmmc.updateMatrix(double3(0.0, 1.0 - physicalAcceptance, physicalAcceptance), oldNA_);
    systemB_.tmmc.updateMatrix(double3(physicalAcceptance, 1.0 - physicalAcceptance, 0.0), oldNB_);
    recorded_ = true;
  }

 private:
  System& systemA_;
  System& systemB_;
  std::size_t oldNA_;
  std::size_t oldNB_;
  bool recorded_{false};
};
}  // namespace

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsSwapMove_CBMC(RandomNumber& random,
                                                                                    System& systemA, System& systemB,
                                                                                    std::size_t selectedComponent)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::GibbsSwapCBMC;
  Component& componentA = systemA.components[selectedComponent];
  Component& componentB = systemB.components[selectedComponent];
  DualTMMCTrial tmmcTrial(systemA, systemB, selectedComponent);

  // Return if there are no molecules of the selected component in system B
  if (systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

  // Index for the new molecule in system A
  std::size_t newMoleculeIndex = systemA.numberOfMolecules();

  // Update move counts statistics for both systems
  componentA.mc_moves_statistics.addTrial(move);
  componentB.mc_moves_statistics.addTrial(move);

  // Retrieve cutoff distances (dual-cutoff aware) and grow type from system A
  double cutOffFrameworkVDWA =
      systemA.forceField.useDualCutOff ? systemA.forceField.dualCutOff : systemA.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDWA =
      systemA.forceField.useDualCutOff ? systemA.forceField.dualCutOff : systemA.forceField.cutOffMoleculeVDW;
  double cutOffCoulombA =
      systemA.forceField.useDualCutOff ? systemA.forceField.dualCutOff : systemA.forceField.cutOffCoulomb;
  Component::GrowType growType = componentA.growType;

  const CBMC::GrowContext growContext{systemA.hasExternalField, systemA.forceField, systemA.simulationBox,
                                      systemA.interpolationGrids, systemA.externalFieldInterpolationGrid,
                                      systemA.framework, systemA.spanOfFrameworkAtoms(), systemA.spanOfMoleculeAtoms(),
                                      systemA.beta, cutOffFrameworkVDWA, cutOffMoleculeVDWA, cutOffCoulombA};

  // Attempt to grow a new molecule in system A using CBMC insertion
  time_begin = std::chrono::steady_clock::now();
  std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
      random, growContext, componentA, selectedComponent, growType, newMoleculeIndex, 1.0, false, false);
  time_end = std::chrono::steady_clock::now();

  // Update CPU time statistics for CBMC insertion (non-Ewald part)
  componentA.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  systemA.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  if (!growData) return std::nullopt;  // Insertion failed, return

  if (systemA.forceField.useDualCutOff)
  {
    // Dual cut-off scheme: correct the grown configuration from the inner cut-off to the full
    // cut-offs, so that Rosenbluth weight and energies behave as if grown at the full cut-offs.
    std::optional<RunningEnergy> correctionNew =
        CBMC::computeDualCutOffCorrection(growContext, componentA, growData->atoms);
    if (!correctionNew.has_value()) return std::nullopt;

    growData->energies += correctionNew.value();
    growData->RosenbluthWeight *= std::exp(-systemA.beta * correctionNew->potentialEnergy());
  }

  // Get new molecule atoms
  std::span<const Atom> newMolecule = std::span(growData->atoms.begin(), growData->atoms.end());

  // Update statistics for successfully constructed molecules in system A
  componentA.mc_moves_statistics.addConstructed(move);

  // Compute Ewald Fourier energy difference for system A
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifferenceA = Interactions::energyDifferenceEwaldFourier(
      systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.trialEik,
      systemA.forceField, systemA.simulationBox, newMolecule, {}, systemA.netCharge);
  time_end = std::chrono::steady_clock::now();

  // Update CPU time statistics for Ewald Fourier computation
  componentA.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  systemA.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Compute tail energy difference for system A
  time_begin = std::chrono::steady_clock::now();
  [[maybe_unused]] RunningEnergy tailEnergyDifferenceA =
      Interactions::computeInterMolecularTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                              systemA.spanOfMoleculeAtoms(), newMolecule, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                                 systemA.spanOfFrameworkAtoms(), newMolecule, {});
  time_end = std::chrono::steady_clock::now();

  // Update CPU time statistics for tail energy computation
  componentA.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);
  systemA.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  // Compute correction factor for Ewald energies in system A
  double correctionFactorEwaldA =
      std::exp(-systemA.beta * (energyFourierDifferenceA.potentialEnergy() + tailEnergyDifferenceA.potentialEnergy()));

  // Select a random molecule of the selected component from system B
  std::size_t selectedMolecule = systemB.randomIntegerMoleculeOfComponent(random, selectedComponent);
  std::span<Atom> molecule = systemB.spanOfMolecule(selectedComponent, selectedMolecule);

  // Retrieve cutoff distances (dual-cutoff aware) from system B
  double cutOffFrameworkVDWB =
      systemB.forceField.useDualCutOff ? systemB.forceField.dualCutOff : systemB.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDWB =
      systemB.forceField.useDualCutOff ? systemB.forceField.dualCutOff : systemB.forceField.cutOffMoleculeVDW;
  double cutOffCoulombB =
      systemB.forceField.useDualCutOff ? systemB.forceField.dualCutOff : systemB.forceField.cutOffCoulomb;

  const CBMC::GrowContext retraceContext{systemB.hasExternalField, systemB.forceField, systemB.simulationBox,
                                         systemB.interpolationGrids, systemB.externalFieldInterpolationGrid,
                                         systemB.framework, systemB.spanOfFrameworkAtoms(),
                                         systemB.spanOfMoleculeAtoms(), systemB.beta, cutOffFrameworkVDWB,
                                         cutOffMoleculeVDWB, cutOffCoulombB};

  // Retrace the selected molecule in system B for deletion using CBMC
  time_begin = std::chrono::steady_clock::now();
  ChainRetraceData retraceData = CBMC::retraceMoleculeSwapDeletion(random, retraceContext, componentB, growType,
                                                                   molecule);
  time_end = std::chrono::steady_clock::now();

  // Update CPU time statistics for CBMC deletion (non-Ewald part)
  componentA.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  systemA.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  if (systemB.forceField.useDualCutOff)
  {
    // Dual cut-off scheme: correct the retraced configuration from the inner cut-off to the full
    // cut-offs, so that Rosenbluth weight and energies behave as if retraced at the full cut-offs.
    std::vector<Atom> oldMolecule(molecule.begin(), molecule.end());
    std::optional<RunningEnergy> correctionOld =
        CBMC::computeDualCutOffCorrection(retraceContext, componentB, oldMolecule);
    if (!correctionOld.has_value()) return std::nullopt;

    retraceData.energies += correctionOld.value();
    retraceData.RosenbluthWeight *= std::exp(-systemB.beta * correctionOld->potentialEnergy());
  }

  // Compute Ewald Fourier energy difference for system B
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifferenceB = Interactions::energyDifferenceEwaldFourier(
      systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.trialEik,
      systemB.forceField, systemB.simulationBox, {}, molecule, systemB.netCharge);
  time_end = std::chrono::steady_clock::now();

  // Update CPU time statistics for Ewald Fourier computation
  componentA.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  systemA.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Compute tail energy difference for system B
  time_begin = std::chrono::steady_clock::now();
  [[maybe_unused]] RunningEnergy tailEnergyDifferenceB =
      Interactions::computeInterMolecularTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                              systemB.spanOfMoleculeAtoms(), {}, molecule) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                                 systemB.spanOfFrameworkAtoms(), {}, molecule);
  time_end = std::chrono::steady_clock::now();

  // Update CPU time statistics for tail energy computation
  componentA.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);
  systemA.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  // Update statistics for retraced molecules in system B
  componentB.mc_moves_statistics.addConstructed(move);

  // Compute correction factor for Ewald energies in system B
  double correctionFactorEwaldB =
      std::exp(systemB.beta * (energyFourierDifferenceB.potentialEnergy() + tailEnergyDifferenceB.potentialEnergy()));

  const double physicalAcceptance =
      (correctionFactorEwaldA * growData->RosenbluthWeight *
       static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent]) *
       systemA.simulationBox.volume) /
      (correctionFactorEwaldB * retraceData.RosenbluthWeight *
       (1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent])) *
       systemB.simulationBox.volume);

  // Bounds must be checked before evaluating a TMMC bias, whose lookup assumes
  // both macrostates are representable.
  if (!tmmcTrial.transferIsInBounds())
  {
    tmmcTrial.recordTransfer(physicalAcceptance);
    return std::nullopt;
  }

  const double tmmcBias = tmmcTrial.biasFactor();
  tmmcTrial.recordTransfer(physicalAcceptance);

  // Apply Metropolis acceptance criterion.  The collection matrix above sees
  // only the unbiased physical probability.
  if (random.uniform() < tmmcBias * physicalAcceptance)
  {
    // Update accepted move statistics for system A
    componentA.mc_moves_statistics.addAccepted(move);

    // Accept Ewald updates and insert the new molecule into system A
    Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.trialEik);
    systemA.insertMolecule(selectedComponent, growData->molecule, growData->atoms);

    // Update accepted move statistics for system B
    componentB.mc_moves_statistics.addAccepted(move);

    // Accept Ewald updates and delete the selected molecule from system B
    Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.trialEik);
    systemB.deleteMolecule(selectedComponent, selectedMolecule, molecule);

    // Return the energy differences for both systems
    return std::make_pair(growData->energies + energyFourierDifferenceA + tailEnergyDifferenceA,
                          -(retraceData.energies - energyFourierDifferenceB - tailEnergyDifferenceB));
  };

  return std::nullopt;
}
