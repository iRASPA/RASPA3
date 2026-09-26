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
import mc_moves_cputime;

double MC_Moves::WidomMove(RandomNumber& random, System& system, std::size_t selectedComponent)
{
  // Set trial moleculeId to something that does not overlap with the current molecules
  std::size_t selectedMolecule = system.numberOfMolecules();

  Move::Types move = Move::Types::Widom;
  Component& component = system.components[selectedComponent];

  // Update move statistics for Widom insertion move.
  component.mc_moves_statistics.addTrial(move);

  // Widom sampling averages the Rosenbluth weight itself, so the chain must be grown with
  // configurational bias regardless of the production scheme: the recoil-growth weight is a valid
  // factor in an acceptance ratio, not the Rosenbluth weight whose average is the excess chemical
  // potential (see 'CBMC::ChainScheme').
  const CBMC::GrowContext growContext =
      system.makeGrowContext().withChainScheme(CBMC::ChainScheme::ConfigurationalBias);

  // Attempt to grow a new molecule using Configurational Bias Monte Carlo (CBMC) insertion.
  std::optional<CBMC::GrowResult> growData =
      timed(system, component, move, Move::Timing::NonEwald,
            [&]
            {
              return CBMC::growNewMolecule(random, growContext, component,
                                           {.componentId = selectedComponent, .moleculeId = selectedMolecule});
            });

  // If molecule growth failed, terminate the move.
  if (!growData) return 0.0;

  std::span<const Atom> newMolecule = std::span(growData->atoms.begin(), growData->atoms.end());

  // Update statistics for successfully constructed molecules.
  component.mc_moves_statistics.addConstructed(move);

  // Compute the energy difference in Ewald Fourier space due to the new molecule.
  RunningEnergy energyFourierDifference =
      timed(system, component, move, Move::Timing::Ewald,
            [&]
            {
              return Interactions::energyDifferenceEwaldFourier(
                  system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik,
                  system.forceField, system.simulationBox, newMolecule, {}, system.netCharge);
            });

  // Compute the tail corrections for the energy due to the new molecule.
  RunningEnergy tailEnergyDifference =
      timed(system, component, move, Move::Timing::Tail,
            [&]
            {
              return Interactions::computeInterMolecularTailEnergyDifference(
                         system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), newMolecule, {}) +
                     Interactions::computeFrameworkMoleculeTailEnergyDifference(
                         system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), newMolecule, {});
            });

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

  // The Rosenbluth weight enters through its exact logarithm: the raw weight of a long chain underflows
  // to zero even when the normalized sample W/W_ideal is of order one.
  return std::exp(std::log(correctionFactorEwald) + growData->logRosenbluthWeight -
                  std::log(idealGasRosenbluthWeight));
}
