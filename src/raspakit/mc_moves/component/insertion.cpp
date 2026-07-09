module;

module mc_moves_insertion;

import std;

import component;
import molecule;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import randomnumbers;
import system;
import energy_factor;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import running_energy;
import forcefield;
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::insertionMove(RandomNumber& random, System& system,
                                                                         std::size_t selectedComponent)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::Swap;
  Component& component = system.components[selectedComponent];

  // Set trial moleculeId to something that does not overlap with the current molecules
  // The 'insertMolecule' routine will place it after the last molecule of the component
  std::size_t selectedMolecule = system.numberOfMolecules();

  // Update swap insertion move counts.
  component.mc_moves_statistics.addTrial(move, 0);

  // Generate a trial molecule with a random position inside the simulation box.
  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      component.equilibratedMoleculeRandomInBox(random, selectedComponent, system.simulationBox);

  std::vector<double3> electricFieldMoleculeNew(trialMolecule.second.size());
  std::vector<double3> electricFieldNeighborDelta;

  // Check if the trial molecule is inside blocked pockets; reject if true.
  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  // Assign molecule ID, component ID, group ID, and set scaling factors for each atom.
  std::for_each(std::begin(trialMolecule.second), std::end(trialMolecule.second),
                [selectedComponent, selectedMolecule](Atom& atom)
                {
                  atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
                  atom.componentId = static_cast<std::uint8_t>(selectedComponent);
                  atom.groupId = static_cast<std::uint8_t>(0);
                  atom.setScaling(1.0);
                });

  // Update constructed counts for swap insertion moves.
  component.mc_moves_statistics.addConstructed(move, 0);

  // compute external field energy contribution
  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, 
      system.externalFieldInterpolationGrid, trialMolecule.second, {});
  if (!externalFieldMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  // compute framework-molecule energy contribution
  std::optional<RunningEnergy> frameworkMolecule;
  if (system.forceField.computePolarization)
  {
    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), electricFieldMoleculeNew, {}, trialMolecule.second, {});
  }
  else
  {
    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialMolecule.second, {});
  }
  if (!frameworkMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  // compute molecule-molecule energy contribution (and, for molecule-molecule polarization, the electric field on
  // the inserted molecule as well as the change of the field on every existing molecule)
  std::optional<RunningEnergy> interMolecule;
  if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
  {
    electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
    interMolecule = Interactions::computeInterMolecularPolarizationElectricFieldDifference(
        system.forceField, system.simulationBox, electricFieldNeighborDelta, electricFieldMoleculeNew,
        std::span<double3>{}, system.spanOfMoleculeAtoms(), trialMolecule.second, {});
  }
  else
  {
    interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, {});
  }
  if (!interMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  // Compute Ewald Fourier energy difference and update CPU time statistics.
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifference;
  if (system.forceField.computePolarization)
  {
    energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, electricFieldMoleculeNew, {}, trialMolecule.second,
        {}, system.netCharge);
  }
  else
  {
    energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, trialMolecule.second, {}, system.netCharge);
  }
  time_end = std::chrono::steady_clock::now();

  component.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Compute tail energy difference and update CPU time statistics.
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), trialMolecule.second, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(
          system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), trialMolecule.second, {});
  time_end = std::chrono::steady_clock::now();

  component.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    // Polarization energy of the inserted molecule
    polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, electricFieldMoleculeNew, {}, trialMolecule.second, {});

    // Polarization energy change of all existing molecules whose field changes due to the insertion
    if (!system.forceField.omitInterPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  // get the total difference in energy
  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                   energyFourierDifference + tailEnergyDifference + polarizationDifference;

  double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
  double preFactor = system.beta * fugacity * system.simulationBox.volume /
                     double(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
  double Pacc = preFactor * std::exp(-system.beta * energyDifference.potentialEnergy());
  std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

  // Calculate acceptance probability and bias from the transition matrix.
  if (system.tmmc.doTMMC)
  {
    std::size_t newN = oldN + 1;
    if (newN > system.tmmc.maxMacrostate)
    {
      return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
    }
  }
  // Check if the new macrostate exceeds the maximum allowed; reject if true.

  // apply acceptance/rejection rule
  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    component.mc_moves_statistics.addAccepted(move, 0);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

    // Apply the field changes on the existing molecules before the new molecule (and its field) is appended.
    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
      for (std::size_t i = 0; i < storedElectricField.size(); ++i)
      {
        storedElectricField[i] += electricFieldNeighborDelta[i];
      }
    }

    system.insertMoleculePolarization(selectedComponent, trialMolecule.first, trialMolecule.second,
                                      electricFieldMoleculeNew);

    return {energyDifference, double3(0.0, 1.0 - Pacc, Pacc)};
  };

  return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}
