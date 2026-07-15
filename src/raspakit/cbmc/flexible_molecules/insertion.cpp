module;

module cbmc_flexible_insertion;

import std;

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
import cbmc_growth_context;
import cbmc_generate_trialorientations_mc;
import forcefield;
import running_energy;
import framework;
import component;
import interpolation_energy_grid;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;

[[nodiscard]] std::optional<ChainGrowData> CBMC::growFlexibleMoleculeChainInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> beadsAlreadyPlaced, std::make_signed_t<std::size_t> skipBackgroundMolecule)
{
  const ForceField &forceField = context.forceField;
  double beta = context.beta;

  std::size_t numberOfBeads = component.connectivityTable.numberOfBeads;
  std::vector<std::vector<Atom>> trialPositions(forceField.numberOfTrialDirections);
  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());

  std::vector<std::size_t> beads_already_placed(beadsAlreadyPlaced.begin(), beadsAlreadyPlaced.end());

  double chain_rosenbluth_weight = 1.0;
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

      if (intraMolecularPotentials.bonds.size() != 1)
      {
        throw std::runtime_error(
            std::format("[CBMC]: multiple bonds detected in 'growFlexibleMoleculeChainInsertion'\n"));
      }

      const BondPotential bond = intraMolecularPotentials.bonds.front();

      for (std::size_t i = 0; i != forceField.numberOfTrialDirections; ++i)
      {
        double bond_length = bond.generateBondLength(random, beta);
        double3 unit_vector = random.randomVectorOnUnitSphere();

        Atom trial_atom = chain_atoms[nextBeads[0]];
        trial_atom.position = chain_atoms[current_bead].position + bond_length * unit_vector;

        trialPositions[i] = {trial_atom};
      }
    }
    else
    {
      // Case: growing a single or multiple bonds with a previous bead present

      std::vector<Atom> trial_orientations =
          generateTrialOrientationsMonteCarloScheme(random, forceField.numberOfTrialMovesPerOpenBead, beta, 
                                                    component, chain_atoms, previous_bead.value(),
                                                    current_bead, nextBeads, intraMolecularPotentials);

      // Rotate around the last bond-vector to obtain the other trial-orientations
      double3 last_bond_vector =
          (chain_atoms[previous_bead.value()].position - chain_atoms[current_bead].position).normalized();

      for (std::size_t i = 0; i != forceField.numberOfTrialDirections; ++i)
      {
        CBMC::TorsionOrientation torsion = CBMC::selectTorsionOrientation(
            random, forceField.numberOfTorsionTrialDirections, beta, chain_atoms, trial_orientations,
            previous_bead.value(), current_bead, nextBeads, last_bond_vector, intraMolecularPotentials, false);

        RosenBluthWeightTorsion[i] = torsion.rosenbluthWeight;
        trialPositions[i] = torsion.positions;
      }
    }

    // Compute the external-energies for the next-beads
    std::vector<CBMC::ChainTrialTorsion> externalEnergies = CBMC::computeExternalNonOverlappingEnergies(
        context, component, trialPositions, RosenBluthWeightTorsion, -1, skipBackgroundMolecule);

    if (externalEnergies.empty()) return std::nullopt;

    // add the intramolecular van der Waals and Coulomb energy for the bead selection
    std::vector<CBMC::ChainTrialTorsion> totalExternalEnergies = externalEnergies;
    for (auto &[external_positions, external_energy, external_torsion] : totalExternalEnergies)
    {
      for (std::size_t k = 0; k != external_positions.size(); ++k)
      {
        std::size_t index = nextBeads[k];
        chain_atoms[index] = external_positions[k];
      }

      external_energy += intraMolecularPotentials.computeInternalIntraVanDerWaalsAndCoulombEnergies(chain_atoms);
    }

    // Select based on the external energy plus the intramolecular non-bonded energy
    std::vector<double> logBoltzmannFactors{};
    logBoltzmannFactors.reserve(forceField.numberOfTrialDirections);
    std::transform(totalExternalEnergies.begin(), totalExternalEnergies.end(), std::back_inserter(logBoltzmannFactors),
                   [&](const CBMC::ChainTrialTorsion &v) { return -beta * v.energy.potentialEnergy(); });

    double rosenbluth_weight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                                [](const double &acc, const double &logBoltzmannFactor)
                                                { return acc + std::exp(logBoltzmannFactor); });

    std::size_t selected = CBMC::selectTrialPosition(random, logBoltzmannFactors);

    const CBMC::ChainTrialTorsion &selectedTrial = externalEnergies[selected];

    chain_rosenbluth_weight *= selectedTrial.torsionWeight * rosenbluth_weight /
                                static_cast<double>(forceField.numberOfTrialDirections);

    if (chain_rosenbluth_weight < forceField.minimumRosenbluthFactor) return std::nullopt;

    chain_external_energies += selectedTrial.energy;

    // Add 'nextBeads' to 'beads_already_placed'
    beads_already_placed.insert(beads_already_placed.end(), nextBeads.begin(), nextBeads.end());

    // Add the selected atoms
    for (std::size_t i = 0; i != nextBeads.size(); ++i)
    {
      std::size_t index = nextBeads[i];
      chain_atoms[index] = selectedTrial.positions[i];
    }

  } while (beads_already_placed.size() < component.connectivityTable.numberOfBeads);

  // Recompute all the internal interactions (including the cross-terms) for the returned energy
  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(chain_atoms);

  // Only bond, bend, and torsion are taken into account during the growth (intra van der Waals and
  // Coulomb enter through the selection of the beads). Correct the Rosenbluth weight with the
  // Boltzmann factor of the remaining internal interactions (Urey-Bradley, inversion-bend,
  // out-of-plane-bend, improper torsion, and the cross-terms).
  RunningEnergy unsampled_internal_energies =
      component.intraMolecularPotentials.computeInternalEnergiesNotSampledDuringGrowth(chain_atoms);
  chain_rosenbluth_weight *= std::exp(-beta * unsampled_internal_energies.potentialEnergy());

  // Copy this configuration so that it can be used as a starting point
  component.grownAtoms = chain_atoms;

  // Build a valid molecule record (center-of-mass, mass) for the grown flexible chain
  double3 com{};
  for (std::size_t i = 0; i != chain_atoms.size(); ++i)
  {
    com += component.definedAtoms[i].second * chain_atoms[i].position;
  }
  com = com / component.totalMass;
  Molecule molecule = Molecule(com, simd_quatd(0.0, 0.0, 0.0, 1.0), component.totalMass,
                               static_cast<std::size_t>(chain_atoms.front().componentId), chain_atoms.size());

  return ChainGrowData(molecule, chain_atoms, chain_external_energies + internal_energies, chain_rosenbluth_weight,
                       0.0);
}
