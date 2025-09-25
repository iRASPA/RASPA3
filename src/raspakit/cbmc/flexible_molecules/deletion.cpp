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

[[nodiscard]] ChainRetraceData CBMC::retraceFlexibleMoleculeChainDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::span<Atom> molecule_atoms, const std::vector<std::size_t> beadsAlreadyPlaced) noexcept
{
  std::size_t numberOfBeads = component.connectivityTable.numberOfBeads;
  std::vector<std::vector<Atom>> trialPositions(forceField.numberOfTrialDirections);

  std::vector<std::size_t> beads_already_placed(beadsAlreadyPlaced.begin(), beadsAlreadyPlaced.end());

  double chain_rosen_bluth_weight = 1.0;
  std::vector<double> RosenBluthWeightTorsion(forceField.numberOfTrialDirections, 1.0);
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

      trialPositions[0] = {molecule_atoms[nextBeads[0]]};
      for (std::size_t i = 1; i != forceField.numberOfTrialDirections; ++i)
      {
        double bond_length = bond.generateBondLength(random, beta);
        double3 unit_vector = random.randomVectorOnUnitSphere();

        Atom trial_atom = molecule_atoms[nextBeads[0]];
        trial_atom.position = molecule_atoms[current_bead].position + bond_length * unit_vector;

        trialPositions[i] = {trial_atom};
      }
    }
    else
    {
      // Case: growing a single or multiple bonds with a previous bead present

      // Get old positions
      std::vector<Atom> trial_orientations(nextBeads.size());
      for (std::size_t k = 0; k < nextBeads.size(); ++k)
      {
        std::size_t index = nextBeads[k];
        trial_orientations[k] = molecule_atoms[index];
      }

      // Rotate around to obtain the other trial-orientations
      double3 last_bond_vector =
          (molecule_atoms[previous_bead.value()].position - molecule_atoms[current_bead].position).normalized();

      for (std::size_t i = 0; i != forceField.numberOfTrialDirections; ++i)
      {
        std::vector<std::pair<std::vector<Atom>, double>> torsion_orientations(
            forceField.numberOfTorsionTrialDirections);

        for (std::size_t j = 0; j != forceField.numberOfTorsionTrialDirections; ++j)
        {
          std::vector<Atom> rotated_atoms = trial_orientations;
          double random_angle = (i == 0 && j == 0) ? 0.0 : (2.0 * random.uniform() - 1.0) * std::numbers::pi;
          for (std::size_t k = 0; k < rotated_atoms.size(); ++k)
          {
            rotated_atoms[k].position =
                molecule_atoms[current_bead].position +
                last_bond_vector.rotateAroundAxis(
                    trial_orientations[k].position - molecule_atoms[current_bead].position, random_angle);
          }

          std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());
          for (std::size_t k = 0; k != nextBeads.size(); ++k)
          {
            std::size_t index = nextBeads[k];
            chain_atoms[index] = rotated_atoms[k];
          }

          double torsion_energy = intraMolecularPotentials.calculateTorsionEnergies(chain_atoms);

          torsion_orientations[j] = {rotated_atoms, torsion_energy};
        }

        std::vector<double> logTorsionBoltzmannFactors{};
        logTorsionBoltzmannFactors.reserve(forceField.numberOfTorsionTrialDirections);
        std::transform(torsion_orientations.begin(), torsion_orientations.end(),
                       std::back_inserter(logTorsionBoltzmannFactors),
                       [&](const std::pair<std::vector<Atom>, double> &v) { return -beta * std::get<1>(v); });

        double rosen_bluth_weight_torsion =
            std::accumulate(logTorsionBoltzmannFactors.begin(), logTorsionBoltzmannFactors.end(), 0.0,
                            [](const double &acc, const double &logTorsionBoltzmannFactor)
                            { return acc + std::exp(logTorsionBoltzmannFactor); });

        std::size_t selected_torsion = (i == 0) ? 0 : CBMC::selectTrialPosition(random, logTorsionBoltzmannFactors);

        auto &[selected_positions, selected_energy] = torsion_orientations[selected_torsion];

        trialPositions[i] = (i == 0) ? trial_orientations : selected_positions;

        RosenBluthWeightTorsion[i] =
            rosen_bluth_weight_torsion / static_cast<double>(forceField.numberOfTorsionTrialDirections);
      }
    }

    // Compute the external-energies for the next-beads
    std::vector<std::tuple<std::vector<Atom>, RunningEnergy, double>> externalEnergies =
        CBMC::computeExternalNonOverlappingEnergies(component, hasExternalField, forceField, simulationBox,
                                                    interpolationGrids, framework, frameworkAtomData, moleculeAtomData,
                                                    cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
                                                    trialPositions, RosenBluthWeightTorsion, -1);

    // add van der Waals
    std::vector<std::tuple<std::vector<Atom>, RunningEnergy, double>> totalExternalEnergies = externalEnergies;
    for (auto &[external_positions, external_energies, external_torsion] : totalExternalEnergies)
    {
      RunningEnergy recomputed_internal_energies =
          intraMolecularPotentials.computeInternalIntraVanDerWaalsEnergies(molecule_atoms);
      external_energies += recomputed_internal_energies;
    }

    std::vector<double> logBoltzmannFactors{};
    logBoltzmannFactors.reserve(forceField.numberOfTrialDirections);
    std::transform(totalExternalEnergies.begin(), totalExternalEnergies.end(), std::back_inserter(logBoltzmannFactors),
                   [&](const std::tuple<std::vector<Atom>, RunningEnergy, double> &v)
                   { return -beta * std::get<1>(v).potentialEnergy(); });

    double rosen_bluth_weight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                                [](const double &acc, const double &logBoltzmannFactor)
                                                { return acc + std::exp(logBoltzmannFactor); });

    std::size_t selected = 0;

    auto &[selected_atom_positions, selected_energies, selected_torsion_rosenbluth_factor] = externalEnergies[selected];

    chain_rosen_bluth_weight *= selected_torsion_rosenbluth_factor * rosen_bluth_weight /
                                static_cast<double>(forceField.numberOfTrialDirections);

    chain_external_energies += selected_energies;

    // Add 'nextBeads' to 'beads_already_placed'
    beads_already_placed.insert(beads_already_placed.end(), nextBeads.begin(), nextBeads.end());

  } while (beads_already_placed.size() < component.connectivityTable.numberOfBeads);

  // Recompute all the internal interactions
  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(molecule_atoms);

  return ChainRetraceData(chain_external_energies + internal_energies, chain_rosen_bluth_weight, 0.0);
}
