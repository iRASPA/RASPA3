module;

module cbmc_rigid_deletion;

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
import cbmc_multiple_first_bead;
import interpolation_energy_grid;

[[nodiscard]] ChainRetraceData CBMC::retraceRigidMoleculeChainDeletion(RandomNumber &random, const GrowContext &context,
                                                                       const Component &component,
                                                                       std::span<Atom> molecule_atoms) noexcept
{
  std::vector<Atom> trialPosition = std::vector<Atom>(molecule_atoms.begin(), molecule_atoms.end());
  std::vector<std::vector<Atom>> trialPositions = {trialPosition};

  for (std::size_t i = 1; i < context.forceField.numberOfTrialDirections; ++i)
  {
    trialPositions.push_back(CBMC::rotateRandomlyAround(random, trialPosition, component.startingBead));
  };

  const std::vector<ChainTrial> externalEnergies = CBMC::computeExternalNonOverlappingEnergies(
      context, component, trialPositions, std::make_signed_t<std::size_t>(component.startingBead));

  std::vector<double> logBoltzmannFactors{};
  std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltzmannFactors),
                 [&](const ChainTrial &v) { return -context.beta * v.energy.potentialEnergy(); });

  double RosenbluthWeight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                            [](const double &acc, const double &logBoltzmannFactor)
                                            { return acc + std::exp(logBoltzmannFactor); });

  return ChainRetraceData(externalEnergies[0].energy,
                          RosenbluthWeight / double(context.forceField.numberOfTrialDirections), 0.0);
}
