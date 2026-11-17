module;

module mc_moves_tethered_proton_hop;

import std;

import mc_moves_acceptance;

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
import mc_moves_cputime;

std::optional<RunningEnergy> MC_Moves::tetheredProtonHopMove(RandomNumber& random, System& system,
                                                             std::size_t selectedComponent,
                                                             std::size_t selectedMolecule)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  const Move::Types move = Move::Types::TetheredProtonHop;
  Component& component = system.components[selectedComponent];

  component.mc_moves_statistics.addTrial(move, 0);

  // Fractional molecules (a Widom or reaction slot of the component) come first among the molecules of a
  // component and are not tethered protons.
  if (selectedMolecule < system.numberOfFractionalMoleculesPerComponent[selectedComponent]) return std::nullopt;
  if (component.tetheredProtonHopSiteGroups.empty()) return std::nullopt;
  if (system.crossLinks.moleculeIsLinked(selectedComponent, selectedMolecule)) return std::nullopt;

  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule& molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  const std::size_t startingBead = component.startingBead;
  const double3 currentPosition = molecule_atoms[startingBead].position;

  // The candidate group of the selected proton is the group that holds the site nearest to it (minimum image).
  // The molecule index does not tell: the fractional molecules of the component shift it, and a deletion moves
  // the last molecule into the freed place, so the index named the group of another Al, whose site-to-site
  // vector then carried the proton into the framework. A site listed in two groups goes to the first.
  struct SiteOfGroup
  {
    std::size_t group;
    std::size_t site;
  };
  const auto nearestSiteOfAllGroups = [&](const double3& position)
  {
    SiteOfGroup nearest{0, 0};
    double nearestDistanceSquared = std::numeric_limits<double>::max();
    for (std::size_t g = 0; g < component.tetheredProtonHopSiteGroups.size(); ++g)
    {
      const std::vector<double3>& sites = component.tetheredProtonHopSiteGroups[g];
      for (std::size_t i = 0; i < sites.size(); ++i)
      {
        const double3 dr =
            system.simulationBox.applyPeriodicBoundaryConditions(system.simulationBox.cell * sites[i] - position);
        const double distanceSquared = dr.length_squared();
        if (distanceSquared < nearestDistanceSquared)
        {
          nearestDistanceSquared = distanceSquared;
          nearest = {g, i};
        }
      }
    }
    return nearest;
  };
  const SiteOfGroup current = nearestSiteOfAllGroups(currentPosition);
  const std::vector<double3>& fractionalSites = component.tetheredProtonHopSiteGroups[current.group];
  if (fractionalSites.size() < 2) return std::nullopt;

  // Candidate positions in Cartesian coordinates (sites are stored as fractional coordinates in the
  // simulation box).
  std::vector<double3> cartesianSites(fractionalSites.size());
  for (std::size_t i = 0; i < fractionalSites.size(); ++i)
  {
    cartesianSites[i] = system.simulationBox.cell * fractionalSites[i];
  }
  const std::size_t currentSite = current.site;

  // Propose a different site of the same group, drawn uniformly among the remaining candidates.
  std::size_t targetSite = static_cast<std::size_t>(random.uniform_integer(0, cartesianSites.size() - 2));
  if (targetSite >= currentSite) ++targetSite;

  // The displacement is the site-to-site vector (minimum image), not the vector from the current position
  // to the target site: the proton keeps its thermal offset from the site centre. Combined with
  // translation moves the proton is generally off-centre, and snapping it onto the exact centre would make
  // the hop irreversible (the reverse hop could never regenerate the off-centre position), breaking
  // detailed balance. With the site-to-site vector the reverse hop (target -> current) applies exactly the
  // opposite displacement and restores the old position, so the proposal is symmetric.
  const double3 displacement =
      system.simulationBox.applyPeriodicBoundaryConditions(cartesianSites[targetSite] - cartesianSites[currentSite]);

  std::pair<Molecule, std::vector<Atom>> trialMolecule = component.translate(molecule, molecule_atoms, displacement);

  // The reverse hop must classify the trial position as belonging to the target site; when the offset is so
  // large that the trial position is nearer to another site (the Voronoi cells of the sites differ in
  // shape), the reverse proposal would never be generated and the move is rejected.
  const SiteOfGroup reverse = nearestSiteOfAllGroups(trialMolecule.second[startingBead].position);
  if (reverse.group != current.group || reverse.site != targetSite)
  {
    return std::nullopt;
  }

  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return std::nullopt;
  }

  std::vector<double3> electricFieldMoleculeNew(molecule_atoms.size());
  std::vector<double3> electricFieldMoleculeOld(molecule_atoms.size());

  // External-field contribution.
  std::optional<RunningEnergy> externalFieldMolecule =
      timed(system, component, move, Move::Timing::ExternalFieldMolecule,
            [&]
            {
              return Interactions::computeExternalFieldEnergyDifference(
                  system.hasExternalField, system.forceField, system.simulationBox,
                  system.externalFieldInterpolationGrid, trialMolecule.second, molecule_atoms);
            });
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // Framework-molecule contribution.
  std::optional<RunningEnergy> frameworkMolecule =
      timed(system, component, move, Move::Timing::FrameworkMolecule,
            [&]() -> std::optional<RunningEnergy>
            {
              if (system.forceField.computePolarization)
              {
                return Interactions::computeFrameworkMoleculeEnergyDifference(
                    system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
                    system.spanOfFrameworkAtoms(), electricFieldMoleculeNew, electricFieldMoleculeOld,
                    trialMolecule.second, molecule_atoms);
              }
              return Interactions::computeFrameworkMoleculeEnergyDifference(
                  system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
                  system.frameworkCellList(), system.spanOfFrameworkAtoms(), trialMolecule.second, molecule_atoms);
            });
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
        system.forceField, system.simulationBox, system.cellList(), system.spanOfMoleculeAtoms(), trialMolecule.second,
        molecule_atoms);
  }
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  if (!interMolecule.has_value()) return std::nullopt;

  // Ewald reciprocal contribution.
  RunningEnergy ewaldFourierEnergy =
      timed(system, component, move, Move::Timing::Ewald,
            [&]() -> RunningEnergy
            {
              if (system.forceField.computePolarization)
              {
                return Interactions::energyDifferenceEwaldFourier(
                    system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik,
                    system.storedEik, system.trialEik, system.forceField, system.simulationBox, system.components,
                    electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second, molecule_atoms);
              }
              return Interactions::energyDifferenceEwaldFourier(
                  system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik,
                  system.forceField, system.simulationBox, system.components, trialMolecule.second, molecule_atoms);
            });

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
  if (random.uniform() < MC_Moves::boundedExp(-system.beta * energyDifference.potentialEnergy()))
  {
    component.mc_moves_statistics.addAccepted(move, 0);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), molecule_atoms.begin());
    molecule = trialMolecule.first;
    system.cellListAtomsMoved(molecule_atoms);

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
