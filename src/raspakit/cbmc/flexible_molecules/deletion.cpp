module;

module cbmc_flexible_deletion;

import std;

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
import energy_status;
import running_energy;
import framework;
import component;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_interactions;
import cbmc_growth_context;
import cbmc_generate_trialorientations_mc;
import cbmc_multiple_first_bead;
import interpolation_energy_grid;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import cbmc_segments;

[[nodiscard]] ChainRetraceData CBMC::retraceFlexibleMoleculeChainDeletion(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> beadsAlreadyPlaced) noexcept
{
  const ForceField &forceField = context.forceField;
  double beta = context.beta;

  std::vector<std::vector<Atom>> trialPositions(forceField.numberOfTrialDirections);

  double chain_rosenbluth_weight = 1.0;
  std::vector<double> RosenBluthWeightTorsion(forceField.numberOfTrialDirections, 1.0);
  RunningEnergy chain_external_energies{};

  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());

  // Same deterministic, group-aware growth sequence as the insertion so that grow and retrace are
  // exactly reversible.
  std::vector<CBMC::GrowSegment> segments = CBMC::buildGrowSegments(component, beadsAlreadyPlaced);

  for (const CBMC::GrowSegment &segment : segments)
  {
    const std::optional<std::size_t> previous_bead = segment.previousBead;
    const std::size_t current_bead = segment.currentBead;
    const std::vector<std::size_t> &nextBeads = segment.nextBeads;
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials = segment.intra;

    std::fill(RosenBluthWeightTorsion.begin(), RosenBluthWeightTorsion.end(), 1.0);

    // Old positions of the beads grown in this segment.
    std::vector<Atom> trial_orientations(nextBeads.size());
    for (std::size_t k = 0; k < nextBeads.size(); ++k)
    {
      trial_orientations[k] = molecule_atoms[nextBeads[k]];
    }

    if (segment.rigidUnit)
    {
      // Case: rigid group retraced as a single rigid body hinged on the anchor, mirroring the
      // insertion's coupled-decoupled scheme: the bend MC contributes no weight, and the torsion
      // Rosenbluth weight is computed with the old configuration pinned as torsion trial 0 of the
      // first trial direction.
      if (previous_bead.has_value())
      {
        double3 last_bond_vector =
            (molecule_atoms[previous_bead.value()].position - molecule_atoms[current_bead].position).normalized();

        // As in the flexible retrace, the old orientation is the shared torsion base for every trial
        // direction (the insertion shares one freshly MC-sampled base among its directions).
        for (std::size_t i = 0; i != forceField.numberOfTrialDirections; ++i)
        {
          CBMC::TorsionOrientation torsion = CBMC::selectTorsionOrientation(
              random, forceField.numberOfTorsionTrialDirections, beta, chain_atoms, trial_orientations,
              previous_bead.value(), current_bead, nextBeads, last_bond_vector, intraMolecularPotentials, i == 0);

          trialPositions[i] = (i == 0) ? trial_orientations : torsion.positions;
          RosenBluthWeightTorsion[i] = torsion.rosenbluthWeight;
        }
      }
      else
      {
        // Seed group: no junction terms, every orientation equally likely, all torsion weights 1.
        trialPositions[0] = trial_orientations;
        for (std::size_t i = 1; i != forceField.numberOfTrialDirections; ++i)
        {
          trialPositions[i] = CBMC::generateRigidUnitOrientationMonteCarloScheme(
              random, forceField.numberOfTrialMovesPerOpenBead, beta, component, chain_atoms, previous_bead,
              current_bead, nextBeads, intraMolecularPotentials);
        }
      }
    }
    else if (!previous_bead.has_value())
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

      // Rotate around to obtain the other trial-orientations
      double3 last_bond_vector =
          (molecule_atoms[previous_bead.value()].position - molecule_atoms[current_bead].position).normalized();

      for (std::size_t i = 0; i != forceField.numberOfTrialDirections; ++i)
      {
        // The first trial direction is the old configuration itself (its first torsion angle is 0
        // and it is the selected one); the others are random torsion rotations.
        CBMC::TorsionOrientation torsion = CBMC::selectTorsionOrientation(
            random, forceField.numberOfTorsionTrialDirections, beta, chain_atoms, trial_orientations,
            previous_bead.value(), current_bead, nextBeads, last_bond_vector, intraMolecularPotentials, i == 0);

        trialPositions[i] = (i == 0) ? trial_orientations : torsion.positions;
        RosenBluthWeightTorsion[i] = torsion.rosenbluthWeight;
      }
    }

    // Compute the external-energies for the next-beads
    std::vector<CBMC::ChainTrialTorsion> externalEnergies =
        CBMC::computeExternalNonOverlappingEnergies(context, component, trialPositions, RosenBluthWeightTorsion, -1);

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

    std::vector<double> logBoltzmannFactors{};
    logBoltzmannFactors.reserve(forceField.numberOfTrialDirections);
    std::transform(totalExternalEnergies.begin(), totalExternalEnergies.end(), std::back_inserter(logBoltzmannFactors),
                   [&](const CBMC::ChainTrialTorsion &v) { return -beta * v.energy.potentialEnergy(); });

    double rosenbluth_weight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                                [](const double &acc, const double &logBoltzmannFactor)
                                                { return acc + std::exp(logBoltzmannFactor); });

    // The old configuration is always the first trial direction of the retrace
    const CBMC::ChainTrialTorsion &selectedTrial = externalEnergies.front();

    chain_rosenbluth_weight *= selectedTrial.torsionWeight * rosenbluth_weight /
                                static_cast<double>(forceField.numberOfTrialDirections);

    chain_external_energies += selectedTrial.energy;

    // Restore the old positions of this segment before growing the next one.
    for (std::size_t k = 0; k != nextBeads.size(); ++k)
    {
      chain_atoms[nextBeads[k]] = molecule_atoms[nextBeads[k]];
    }
  }

  // Recompute all the internal interactions (including the terms not sampled during growth,
  // so that the returned energies contain the cross-terms)
  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(molecule_atoms);

  // Only bond, bend, and torsion are taken into account during the retracing (intra van der Waals
  // and Coulomb enter through the selection of the beads). Correct the Rosenbluth weight with
  // the Boltzmann factor of the remaining internal interactions (Urey-Bradley, inversion-bend,
  // out-of-plane-bend, improper torsion, and the cross-terms).
  RunningEnergy unsampled_internal_energies =
      component.intraMolecularPotentials.computeInternalEnergiesNotSampledDuringGrowth(molecule_atoms);
  chain_rosenbluth_weight *= std::exp(-beta * unsampled_internal_energies.potentialEnergy());

  return ChainRetraceData(chain_external_energies + internal_energies, chain_rosenbluth_weight, 0.0);
}
