module;

module mc_moves_tethered_proton_hop;

import std;

import component;
import molecule;
import atom;
import double3;
import simulationbox;
import randomnumbers;
import system;
import running_energy;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::tetheredProtonHopMove(RandomNumber& random, System& system,
                                                             std::size_t selectedComponent,
                                                             std::size_t selectedMolecule)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  const Move::Types move = Move::Types::TetheredProtonHop;
  Component& component = system.components[selectedComponent];

  component.mc_moves_statistics.addTrial(move, 0);

  // The candidate group of the selected proton. With one proton per tether site (one integer
  // molecule per group) the molecule index selects the group directly.
  if (selectedMolecule >= component.tetheredProtonHopSiteGroups.size()) return std::nullopt;
  const std::vector<double3>& fractionalSites = component.tetheredProtonHopSiteGroups[selectedMolecule];
  if (fractionalSites.size() < 2) return std::nullopt;

  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule& molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  const std::size_t startingBead = component.startingBead;
  const double3 currentPosition = molecule_atoms[startingBead].position;

  // Candidate positions in Cartesian coordinates (sites are stored as fractional coordinates in the
  // simulation box).
  std::vector<double3> cartesianSites(fractionalSites.size());
  for (std::size_t i = 0; i < fractionalSites.size(); ++i)
  {
    cartesianSites[i] = system.simulationBox.cell * fractionalSites[i];
  }

  // The site the proton currently occupies is the nearest candidate (minimum image).
  std::size_t currentSite = 0;
  double nearestDistanceSquared = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < cartesianSites.size(); ++i)
  {
    const double3 dr = system.simulationBox.applyPeriodicBoundaryConditions(cartesianSites[i] - currentPosition);
    const double distanceSquared = dr.length_squared();
    if (distanceSquared < nearestDistanceSquared)
    {
      nearestDistanceSquared = distanceSquared;
      currentSite = i;
    }
  }

  // Propose a different site of the same group, drawn uniformly among the remaining candidates.
  std::size_t targetSite = static_cast<std::size_t>(random.uniform_integer(0, cartesianSites.size() - 2));
  if (targetSite >= currentSite) ++targetSite;

  // The displacement carries the proton from its current position to the chosen candidate site,
  // wrapped to the minimum image so that the jump crosses periodic boundaries correctly.
  const double3 displacement =
      system.simulationBox.applyPeriodicBoundaryConditions(cartesianSites[targetSite] - currentPosition);

  std::pair<Molecule, std::vector<Atom>> trialMolecule = component.translate(molecule, molecule_atoms, displacement);

  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return std::nullopt;
  }

  std::vector<double3> electricFieldMoleculeNew(molecule_atoms.size());
  std::vector<double3> electricFieldMoleculeOld(molecule_atoms.size());

  // External-field contribution.
  time_begin = std::chrono::steady_clock::now();
  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
      trialMolecule.second, molecule_atoms);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::ExternalFieldMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::ExternalFieldMolecule] += (time_end - time_begin);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // Framework-molecule contribution.
  time_begin = std::chrono::steady_clock::now();
  std::optional<RunningEnergy> frameworkMolecule;
  if (system.forceField.computePolarization)
  {
    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second,
        molecule_atoms);
  }
  else
  {
    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialMolecule.second, molecule_atoms);
  }
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::FrameworkMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::FrameworkMolecule] += (time_end - time_begin);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // Molecule-molecule contribution.
  time_begin = std::chrono::steady_clock::now();
  std::vector<double3> electricFieldNeighborDelta;
  std::optional<RunningEnergy> interMolecule;
  if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
  {
    electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
    interMolecule = Interactions::computeInterMolecularPolarizationElectricFieldDifference(
        system.forceField, system.simulationBox, electricFieldNeighborDelta, electricFieldMoleculeNew,
        electricFieldMoleculeOld, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
  }
  else
  {
    interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
  }
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  if (!interMolecule.has_value()) return std::nullopt;

  // Ewald reciprocal contribution.
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy ewaldFourierEnergy;
  if (system.forceField.computePolarization)
  {
    ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.trialEik, system.forceField, system.simulationBox, electricFieldMoleculeNew, electricFieldMoleculeOld,
        trialMolecule.second, molecule_atoms);
  }
  else
  {
    ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
        system.simulationBox, trialMolecule.second, molecule_atoms);
  }
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second, molecule_atoms);

    if (!system.forceField.omitInterPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                   ewaldFourierEnergy + polarizationDifference;

  component.mc_moves_statistics.addConstructed(move, 0);

  // Symmetric proposal: plain Metropolis acceptance on the energy difference.
  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    component.mc_moves_statistics.addAccepted(move, 0);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), molecule_atoms.begin());
    molecule = trialMolecule.first;

    if (system.forceField.computePolarization)
    {
      if (!system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }
      std::span<double3> electricFieldMolecule = system.spanElectricFieldOld(selectedComponent, selectedMolecule);
      std::copy(electricFieldMoleculeNew.begin(), electricFieldMoleculeNew.end(), electricFieldMolecule.begin());
    }

    return energyDifference;
  }

  return std::nullopt;
}
