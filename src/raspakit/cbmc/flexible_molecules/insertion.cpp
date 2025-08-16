module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <type_traits>
#include <vector>
#endif

module cbmc_flexible_insertion;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import randomnumbers;
import units;
import component;
import molecule;
import atom;
import double3;
import simd_quatd;
import double3x3;
import simulationbox;
import energy_status;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_multiple_first_bead;
import cbmc_interactions;
import cbmc_generate_trialorientations_mc;
import forcefield;
import energy_factor;
import running_energy;
import framework;
import component;
import interpolation_energy_grid;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;


[[nodiscard]] std::optional<ChainGrowData> CBMC::growFlexibleMoleculeChainInsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t startingBead, std::vector<Atom> molecule_atoms,
    std::size_t numberOfTrialDirections)
{
  std::size_t numberOfBeads = component.connectivityTable.numberOfBeads;
  std::vector<std::vector<Atom>> trialPositions(numberOfTrialDirections);
  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());

  std::vector<std::size_t> beads_already_placed{startingBead};
  const Atom first_bead = molecule_atoms[startingBead];

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

      if (intraMolecularPotentials.bonds.size() != 1)
      {
        throw std::runtime_error(
            std::format("[CBMC]: multiple bonds detected in 'growFlexibleMoleculeChainInsertion'\n"));
      }

      const BondPotential bond = intraMolecularPotentials.bonds.front();

      for(std::size_t i = 0; i != numberOfTrialDirections; ++i)
      {
        double bond_length = bond.generateBondLength(random, beta);
        double3 unit_vector = random.randomVectorOnUnitSphere();

        Atom trial_atom = chain_atoms[nextBeads[0]];
        trial_atom.position = first_bead.position + bond_length * unit_vector;

        trialPositions[i] = { trial_atom };
      }
    }
    else
    {
      // Case: growing a single or multiple bonds with a previous bead present
      
      std::vector<Atom> trial_orientations = generateTrialOrientationsMonteCarloScheme(
          random, beta, component, chain_atoms, previous_bead.value(), current_bead, nextBeads, numberOfTrialDirections,
          intraMolecularPotentials);

      // Rotate around the last bond-vector to obtain the other trial-orientations
      double3 last_bond_vector =
          (chain_atoms[previous_bead.value()].position - chain_atoms[current_bead].position).normalized();

      for (std::size_t i = 0; i != numberOfTrialDirections; ++i)
      {
        double random_angle = (2.0 * random.uniform() - 1.0) * std::numbers::pi;

        std::vector<Atom> rotated_atoms = trial_orientations;
        for(Atom & rotated_atom : rotated_atoms)
        {
          rotated_atom.position = chain_atoms[current_bead].position +
              last_bond_vector.rotateAroundAxis(rotated_atom.position - chain_atoms[current_bead].position, random_angle);
        }

        trialPositions[i] = rotated_atoms;
      }
    }

    // Compute the external-energies for the next-beads
    const std::vector<std::pair<std::vector<Atom>, RunningEnergy>> externalEnergies =
        CBMC::computeExternalNonOverlappingEnergies(
            component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtomData,
            moleculeAtomData, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions, -1);

    if (externalEnergies.empty()) return std::nullopt;

    std::vector<double> logBoltmannFactors{};
    std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltmannFactors),
                   [&](const std::pair<std::vector<Atom>, RunningEnergy> &v)
                   { return -beta * std::get<1>(v).potentialEnergy(); });

    std::size_t selected = CBMC::selectTrialPosition(random, logBoltmannFactors);

    double rosen_bluth_weight = std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                                [](const double &acc, const double &logBoltmannFactor)
                                                { return acc + std::exp(logBoltmannFactor); });

    if (rosen_bluth_weight < forceField.minimumRosenbluthFactor) return std::nullopt;

    chain_rosen_bluth_weight *= rosen_bluth_weight / static_cast<double>(numberOfTrialDirections);

    const std::vector<Atom> &selected_trial_atoms = externalEnergies[selected].first;
    chain_external_energies += externalEnergies[selected].second;

    // Add 'nextBeads' to 'beads_already_placed'
    beads_already_placed.insert(beads_already_placed.end(), nextBeads.begin(), nextBeads.end());

    // Add the selected atoms
    for (std::size_t i = 0; i != nextBeads.size(); ++i)
    {
      std::size_t index = nextBeads[i];
      chain_atoms[index] = selected_trial_atoms[i];
    }

  } while (beads_already_placed.size() < component.connectivityTable.numberOfBeads);

  // Recompute all the internal interactions
  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(chain_atoms);

  return ChainGrowData({}, chain_atoms, chain_external_energies + internal_energies, chain_rosen_bluth_weight, 0.0);
}
