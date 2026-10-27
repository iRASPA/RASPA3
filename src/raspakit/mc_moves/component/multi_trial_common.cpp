module;

module mc_moves_multi_trial_common;

import std;

import atom;
import running_energy;
import system;
import forcefield;
import interactions_ewald;
import interactions_intermolecular;
import interactions_framework_molecule;
import cbmc_multi_trial;
import cbmc_multi_trial_cpu;
import cbmc_multi_trial_device;
import spatial_decomposition_settings;

const CBMC::MultiTrialGrower& MC_Moves::multiTrialGrower(const System& system)
{
  // The device grower proposes rigid molecules on the system's CBMC device and forwards everything else to the
  // CPU grower itself; a system without a device uses the CPU grower directly.
  if (system.cbmcDevice != PairDevice::CPU) return CBMC::deviceMultiTrialGrower(system.cbmcDevice);
  return CBMC::cpuMultiTrialGrower();
}

MC_Moves::LongRangeDifference MC_Moves::longRangeDifference(System& system, std::span<const Atom> inserted,
                                                            std::span<const Atom> removed, bool includeTail)
{
  LongRangeDifference difference;
  difference.fourier = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, system.components, inserted, removed, system.netCharge);
  if (includeTail)
  {
    difference.tail = Interactions::computeInterMolecularTailEnergyDifference(
                          system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), inserted, removed) +
                      Interactions::computeFrameworkMoleculeTailEnergyDifference(
                          system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), inserted, removed);
  }
  return difference;
}

CBMC::LongRangeCorrection MC_Moves::longRangeCorrection(System& system, std::span<const Atom> removed,
                                                        bool includeTail, double baseline)
{
  return [&system, removed, includeTail, baseline](std::span<const Atom> trial) -> double
  {
    if (!removed.empty() && trial.data() == removed.data() && trial.size() == removed.size()) return -baseline;
    return longRangeDifference(system, trial, removed, includeTail).potentialEnergy() - baseline;
  };
}

void MC_Moves::checkMultiTrialSupported(const System& system, std::string_view moveName)
{
  if (system.forceField.computePolarization)
  {
    throw std::runtime_error(
        std::format("[{}]: polarization is not supported by the multiple-try moves\n", moveName));
  }
}
