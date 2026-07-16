module;

module mc_moves_widom;

import std;

import component;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import cbmc_interactions;
import randomnumbers;
import system;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import running_energy;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

double MC_Moves::WidomMove(RandomNumber& random, System& system, std::size_t selectedComponent)
{
  // Set trial moleculeId to something that does not overlap with the current molecules
  std::size_t selectedMolecule = system.numberOfMolecules();

  Move::Types move = Move::Types::Widom;
  Component& component = system.components[selectedComponent];
  std::chrono::steady_clock::time_point t1, t2;

  // Update move statistics for Widom insertion move.
  component.mc_moves_statistics.addTrial(move);

  // Determine cutoff distances based on whether dual cutoff is used.
  double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;
  Component::GrowType growType = component.growType;

  const CBMC::GrowContext growContext{system.hasExternalField, system.forceField, system.simulationBox,
                                      system.interpolationGrids, system.externalFieldInterpolationGrid,
                                      system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
                                      system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb};

  // Attempt to grow a new molecule using Configurational Bias Monte Carlo (CBMC) insertion.
  t1 = std::chrono::steady_clock::now();
  std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
      random, growContext, component, selectedComponent, growType, selectedMolecule, 1.0, false, false);
  t2 = std::chrono::steady_clock::now();

  component.mc_moves_cputime[move][Move::Timing::NonEwald] += (t2 - t1);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (t2 - t1);

  // If molecule growth failed, terminate the move.
  if (!growData) return 0.0;

  if (system.forceField.useDualCutOff)
  {
    // Dual cut-off scheme: correct the Widom Rosenbluth weight from the inner cut-off to the full
    // cut-offs.
    std::optional<RunningEnergy> correctionNew =
        CBMC::computeDualCutOffCorrection(growContext, component, growData->atoms);
    if (!correctionNew.has_value()) return 0.0;

    growData->energies += correctionNew.value();
    growData->RosenbluthWeight *= std::exp(-system.beta * correctionNew->potentialEnergy());
  }

  [[maybe_unused]] std::span<const Atom> newMolecule = std::span(growData->atoms.begin(), growData->atoms.end());

  // Check if the new molecule is inside blocked pockets; if so, abort the move.
  if (system.insideBlockedPockets(component, newMolecule))
  {
    return 0.0;
  }

  // Update statistics for successfully constructed molecules.
  component.mc_moves_statistics.addConstructed(move);

  // Compute the energy difference in Ewald Fourier space due to the new molecule.
  t1 = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, newMolecule, {}, system.netCharge);
  t2 = std::chrono::steady_clock::now();

  component.mc_moves_cputime[move][Move::Timing::Ewald] += (t2 - t1);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (t2 - t1);

  // Compute the tail corrections for the energy due to the new molecule.
  t1 = std::chrono::steady_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), newMolecule, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMolecule, {});
  t2 = std::chrono::steady_clock::now();

  component.mc_moves_cputime[move][Move::Timing::Tail] += (t2 - t1);
  system.mc_moves_cputime[move][Move::Timing::Tail] += (t2 - t1);

  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    std::vector<double3> newElectricField(newMolecule.size());
    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), newElectricField, {},
                                                                  growData->atoms, {});

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.trialEik, system.forceField, system.simulationBox, newElectricField, {}, growData->atoms, {});

    if (!system.forceField.omitInterPolarization)
    {
      std::vector<double3> electricFieldNeighborDelta(system.spanOfMoleculeAtoms().size(),
                                                      double3(0.0, 0.0, 0.0));
      [[maybe_unused]] std::optional<RunningEnergy> interPolarizationEnergy =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, electricFieldNeighborDelta, newElectricField,
              std::span<double3>{}, system.spanOfMoleculeAtoms(), growData->atoms, {});

      polarizationDifference = Interactions::computePolarizationEnergyDifference(system.forceField, newElectricField,
                                                                               {}, growData->atoms, {});
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
    else
    {
      polarizationDifference = Interactions::computePolarizationEnergyDifference(system.forceField, newElectricField,
                                                                               {}, growData->atoms, {});
    }
  }

  // Compute the correction factor from Ewald, tail and polarization energy differences.
  double correctionFactorEwald = std::exp(-system.beta * (energyFourierDifference.potentialEnergy() +
                                                          tailEnergyDifference.potentialEnergy() +
                                                          polarizationDifference.potentialEnergy()));

  double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);

  return correctionFactorEwald * growData->RosenbluthWeight / idealGasRosenbluthWeight;
}
