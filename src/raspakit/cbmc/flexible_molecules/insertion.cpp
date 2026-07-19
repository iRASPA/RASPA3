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
import cbmc_interactions;
import cbmc_growth_context;
import forcefield;
import running_energy;
import framework;
import interpolation_energy_grid;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import cbmc_growth_plan;
import cbmc_operators;

[[nodiscard]] std::optional<ChainGrowData> CBMC::growFlexibleMoleculeChainInsertion(
    RandomNumber &random, const GrowContext &context, Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> beadsAlreadyPlaced, std::make_signed_t<std::size_t> skipBackgroundMolecule)
{
  const ForceField &forceField = context.forceField;
  double beta = context.beta;

  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());

  double chain_rosenbluth_weight = 1.0;
  RunningEnergy chain_external_energies{};

  // Deterministic growth plan over the fragment graph (flexible beads, hinged rigid bodies, and
  // ring-closure of cyclic clusters), shared with the retrace so grow and retrace are reversible.
  std::vector<CBMC::GrowStep> plan = CBMC::buildGrowthPlan(component, beadsAlreadyPlaced);

  for (const CBMC::GrowStep &step : plan)
  {
    const std::vector<std::size_t> &nextBeads = step.nextBeads;
    const Potentials::IntraMolecularPotentials &intra = step.intra;

    // Generate all trial directions for this step (the operator engine handles the seed / attach /
    // ring-closure cases, the rigid-body tilt, and the coupled-decoupled torsion selection).
    std::vector<CBMC::StepTrial> stepTrials =
        CBMC::generateGrowTrials(random, forceField, beta, component, chain_atoms, step, forceField.numberOfTrialDirections);

    std::vector<std::vector<Atom>> trialPositions(forceField.numberOfTrialDirections);
    std::vector<double> RosenBluthWeightTorsion(forceField.numberOfTrialDirections, 1.0);
    for (std::size_t i = 0; i != forceField.numberOfTrialDirections; ++i)
    {
      trialPositions[i] = stepTrials[i].positions;
      RosenBluthWeightTorsion[i] = stepTrials[i].torsionWeight;
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
        chain_atoms[nextBeads[k]] = external_positions[k];
      }
      external_energy += intra.computeInternalIntraVanDerWaalsAndCoulombEnergies(chain_atoms);
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

    // Fold the internal terms whose atoms are all placed by this step, but which are not sampled by
    // the growth scheme (Urey-Bradley, inversion/out-of-plane bend, improper torsion, cross-terms),
    // into this step's Rosenbluth weight instead of a single end-of-growth correction.
    for (std::size_t i = 0; i != nextBeads.size(); ++i)
    {
      chain_atoms[nextBeads[i]] = selectedTrial.positions[i];
    }
    RunningEnergy stepUnsampled = intra.computeInternalEnergiesNotSampledDuringGrowth(chain_atoms);
    chain_rosenbluth_weight *= std::exp(-beta * stepUnsampled.potentialEnergy());
    if (chain_rosenbluth_weight < forceField.minimumRosenbluthFactor) return std::nullopt;
  }

  // Recompute all the internal interactions (including the cross-terms) for the returned energy.
  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(chain_atoms);

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
