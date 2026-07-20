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
import running_energy;
import framework;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_interactions;
import cbmc_growth_context;
import interpolation_energy_grid;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import cbmc_growth_plan;
import cbmc_operators;

[[nodiscard]] ChainRetraceData CBMC::retraceFlexibleMoleculeChainDeletion(
    RandomNumber &random, const GrowContext &context, const Component &component, std::span<Atom> molecule_atoms,
    const std::vector<std::size_t> beadsAlreadyPlaced) noexcept
{
  const ForceField &forceField = context.forceField;
  double beta = context.beta;

  double chain_rosenbluth_weight = 1.0;
  RunningEnergy chain_external_energies{};

  std::vector<Atom> chain_atoms(molecule_atoms.begin(), molecule_atoms.end());

  // Same deterministic growth plan as the insertion so grow and retrace are exactly reversible.
  const std::vector<CBMC::GrowStep> &plan = component.growthPlan(beadsAlreadyPlaced);

  for (const CBMC::GrowStep &step : plan)
  {
    const std::vector<std::size_t> &nextBeads = step.nextBeads;
    const Potentials::IntraMolecularPotentials &intra = step.intra;

    // The old positions of this step's beads must be present in 'chain_atoms' for the retrace.
    for (std::size_t k = 0; k != nextBeads.size(); ++k) chain_atoms[nextBeads[k]] = molecule_atoms[nextBeads[k]];

    std::vector<CBMC::StepTrial> stepTrials = CBMC::generateRetraceTrials(
        random, forceField, beta, component, chain_atoms, step, forceField.numberOfTrialDirections);

    std::vector<std::vector<Atom>> trialPositions(forceField.numberOfTrialDirections);
    std::vector<double> RosenBluthWeightTorsion(forceField.numberOfTrialDirections, 1.0);
    for (std::size_t i = 0; i != forceField.numberOfTrialDirections; ++i)
    {
      trialPositions[i] = stepTrials[i].positions;
      RosenBluthWeightTorsion[i] = stepTrials[i].torsionWeight;
    }

    std::vector<CBMC::ChainTrialTorsion> externalEnergies =
        CBMC::computeExternalNonOverlappingEnergies(context, component, trialPositions, RosenBluthWeightTorsion, -1);

    std::vector<CBMC::ChainTrialTorsion> totalExternalEnergies = externalEnergies;
    for (auto &[external_positions, external_energy, external_torsion] : totalExternalEnergies)
    {
      for (std::size_t k = 0; k != external_positions.size(); ++k)
      {
        chain_atoms[nextBeads[k]] = external_positions[k];
      }
      external_energy += intra.computeInternalIntraVanDerWaalsAndCoulombEnergies(chain_atoms);
    }

    std::vector<double> logBoltzmannFactors{};
    logBoltzmannFactors.reserve(forceField.numberOfTrialDirections);
    std::transform(totalExternalEnergies.begin(), totalExternalEnergies.end(), std::back_inserter(logBoltzmannFactors),
                   [&](const CBMC::ChainTrialTorsion &v) { return -beta * v.energy.potentialEnergy(); });

    double rosenbluth_weight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                                [](const double &acc, const double &logBoltzmannFactor)
                                                { return acc + std::exp(logBoltzmannFactor); });

    // The old configuration is always the first trial direction of the retrace.
    const CBMC::ChainTrialTorsion &selectedTrial = externalEnergies.front();

    chain_rosenbluth_weight *= selectedTrial.torsionWeight * rosenbluth_weight /
                                static_cast<double>(forceField.numberOfTrialDirections);

    chain_external_energies += selectedTrial.energy;

    // Fold this step's not-sampled internal terms into the weight (mirrors the insertion).
    for (std::size_t k = 0; k != nextBeads.size(); ++k) chain_atoms[nextBeads[k]] = molecule_atoms[nextBeads[k]];
    RunningEnergy stepUnsampled = intra.computeInternalEnergiesNotSampledDuringGrowth(chain_atoms);
    chain_rosenbluth_weight *= std::exp(-beta * stepUnsampled.potentialEnergy());
  }

  RunningEnergy internal_energies = component.intraMolecularPotentials.computeInternalEnergies(molecule_atoms);

  return ChainRetraceData(chain_external_energies + internal_energies, chain_rosenbluth_weight, 0.0);
}
