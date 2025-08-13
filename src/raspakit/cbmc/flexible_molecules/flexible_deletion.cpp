module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <type_traits>
#include <vector>
#include <numbers>
#endif

module cbmc_flexible_deletion;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import randomnumbers;
import component;
import molecule;
import atom;
import double3;
import simd_quatd;
import double3x3;
import simulationbox;
import energy_status;
import forcefield;
import energy_factor;
import energy_status;
import running_energy;
import framework;
import component;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_interactions;
import cbmc_generate_trialorientations_mc;
import cbmc_multiple_first_bead;
import interpolation_energy_grid;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;


[[nodiscard]] ChainData CBMC::retraceFlexibleMoleculeSwapDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forcefield, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, [[maybe_unused]] std::size_t selectedComponent, [[maybe_unused]] std::size_t selectedMolecule,
    std::span<Atom> molecule_atoms, double scaling, std::size_t numberOfTrialDirections) noexcept
{
  std::size_t startingBead = components[selectedComponent].startingBead;

  const FirstBeadData firstBeadData = CBMC::retraceMultipleFirstBeadSwapDeletion(
      random, component, hasExternalField, forcefield, simulationBox, interpolationGrids, framework, frameworkAtomData,
      moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, molecule_atoms[startingBead], scaling,
      numberOfTrialDirections);

  if (molecule_atoms.size() == 1)
  {
    return ChainData(
        Molecule(double3(), simd_quatd(), component.totalMass, component.componentId, component.definedAtoms.size()),
        std::vector<Atom>(molecule_atoms.begin(), molecule_atoms.end()), firstBeadData.energies, firstBeadData.RosenbluthWeight,
        0.0);
  }

  const ChainData chainData =
      retraceFlexibleMoleculeChainDeletion(random, component, hasExternalField, forcefield, simulationBox, interpolationGrids, framework,
                        frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
                        startingBead, scaling, molecule_atoms, numberOfTrialDirections);

  return ChainData(
      Molecule(double3(), simd_quatd(), component.totalMass, component.componentId, component.definedAtoms.size()),
      std::vector<Atom>(molecule_atoms.begin(), molecule_atoms.end()), firstBeadData.energies + chainData.energies,
      firstBeadData.RosenbluthWeight * chainData.RosenbluthWeight, 0.0);
}

[[nodiscard]] ChainData retraceFlexibleMoleculeChainDeletion(RandomNumber &random, const Component &component, bool hasExternalField,
                                                             const ForceField &forceField, const SimulationBox &simulationBox,
                                                             const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
                                                             const std::optional<Framework> &framework,
                                                             std::span<const Atom> frameworkAtomData, std::span<const Atom> moleculeAtomData,
                                                             double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
                                                             double cutOffCoulomb, std::size_t startingBead,
                                                             [[maybe_unused]] double scaling, std::span<Atom> molecule_atoms,
                                                             std::size_t numberOfTrialDirections) noexcept
{
  std::size_t numberOfBeads = component.connectivityTable.numberOfBeads;
  std::vector<std::vector<Atom>> trialPositions(numberOfTrialDirections);

  std::vector<std::size_t> beads_already_placed{component.startingBead};
  const Atom first_bead = molecule_atoms[component.startingBead];

  double chain_rosen_bluth_weight = 1.0;
  RunningEnergy chain_external_energies{};

  do
  {
    auto [previous_bead, current_bead, nextBeads] = component.connectivityTable.nextBeads(beads_already_placed);

    Potentials::IntraMolecularPotentials intraMolecularPotentials =
        component.intraMolecularPotentials.filteredInteractions(numberOfBeads, beads_already_placed, nextBeads);

    if (!previous_bead.has_value())
    {
      // Case: growing a single bond with no previous beads
      //       for example: dimer, or starting in the middle of a linear chain

      const BondPotential bond = intraMolecularPotentials.bonds.front();

      trialPositions[0] = { molecule_atoms[nextBeads[0]] };
      for(std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        double bond_length = bond.generateBondLength(random, beta);
        double3 unit_vector = random.randomVectorOnUnitSphere();

        Atom trial_atom = molecule_atoms[nextBeads[0]];
        trial_atom.position = first_bead.position + bond_length * unit_vector;

        trialPositions[i] = { trial_atom };
      }

    }
    else
    {
      // Case: growing a single or multiple bonds with a previous bead present
 
      std::vector<Atom> trial_orientations(nextBeads.size());
      for(std::size_t i = 0; i < nextBeads.size(); ++i)
      {
        trial_orientations[i] = molecule_atoms[nextBeads[i]];
      }

      // Rotate around to obtain the other trial-orientations
      double3 last_bond_vector =
          (molecule_atoms[previous_bead.value()].position - molecule_atoms[current_bead].position).normalized();

      trialPositions[0] = trial_orientations;
      for (std::size_t i = 1; i != numberOfTrialDirections; ++i)
      {
        double random_angle = (2.0 * random.uniform() - 1.0) * std::numbers::pi;

        std::vector<Atom> rotated_atoms = trial_orientations;
        for(Atom &rotated_atom : rotated_atoms)
        {
          rotated_atom.position = molecule_atoms[current_bead].position +
               last_bond_vector.rotateAroundAxis(rotated_atom.position - molecule_atoms[current_bead].position, random_angle);
        }
        trialPositions[i] = rotated_atoms;
      }
    }

    // Compute the external-energies for the next-beads
    const std::vector<std::pair<std::vector<Atom>, RunningEnergy>> externalEnergies =
        CBMC::computeExternalNonOverlappingEnergies(
            component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtomData,
            moleculeAtomData, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions, -1);

    std::vector<double> logBoltmannFactors{};
    std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltmannFactors),
                   [&](const std::pair<std::vector<Atom>, RunningEnergy> &v)
                   { return -beta * std::get<1>(v).potentialEnergy(); });

    std::size_t selected = 0;

    double rosen_bluth_weight = std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                                [](const double &acc, const double &logBoltmannFactor)
                                                { return acc + std::exp(logBoltmannFactor); });

    chain_rosen_bluth_weight *= rosen_bluth_weight / static_cast<double>(numberOfTrialDirections);

    const std::vector<Atom> &selected_trial_atoms = externalEnergies[selected].first;
    chain_external_energies += externalEnergies[selected].second;

    // Add 'nextBeads' to 'beads_already_placed'
    beads_already_placed.insert(beads_already_placed.end(), nextBeads.begin(), nextBeads.end());

  } while (beads_already_placed.size() < component.connectivityTable.numberOfBeads);

  // Recompute all the internal interactions
  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(molecule_atoms);

  return ChainData({}, std::vector<Atom>(molecule_atoms.begin(), molecule_atoms.end()), chain_external_energies + internal_energies, chain_rosen_bluth_weight, 0.0);
}
