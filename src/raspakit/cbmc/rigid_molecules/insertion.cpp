module;

module cbmc_rigid_insertion;

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
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_multiple_first_bead;
import cbmc_interactions;
import cbmc_growth_context;
import forcefield;
import running_energy;
import framework;
import component;
import interpolation_energy_grid;

[[nodiscard]] std::optional<ChainGrowData> CBMC::growRigidMoleculeChainInsertion(
    RandomNumber &random, const GrowContext &context, const Component &component, std::size_t selectedComponent,
    std::span<Atom> molecule_atoms, std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<std::pair<Molecule, std::vector<Atom>>> trialPositions(context.forceField.numberOfTrialDirections);

  std::size_t starting_bead = component.startingBead;

  // randomly rotated configurations around the starting bead
  for (std::size_t i = 0; i < context.forceField.numberOfTrialDirections; ++i)
  {
    simd_quatd orientation = random.randomSimdQuatd();
    std::vector<Atom> randomlyRotatedAtoms = CBMC::rotateRandomlyAround(orientation, molecule_atoms, starting_bead);

    double3 com = component.computeCenterOfMass(randomlyRotatedAtoms);

    trialPositions[i] = {
        Molecule(com, orientation, component.totalMass, selectedComponent, component.definedAtoms.size()),
        randomlyRotatedAtoms};
  };

  const std::vector<MoleculeTrial> externalEnergies = CBMC::computeExternalNonOverlappingEnergies(
      context, component, trialPositions, std::make_signed_t<std::size_t>(component.startingBead),
      skipBackgroundMolecule);
  if (externalEnergies.empty()) return std::nullopt;

  std::vector<double> logBoltzmannFactors{};
  std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltzmannFactors),
                 [&](const MoleculeTrial &v) { return -context.beta * v.energy.potentialEnergy(); });

  std::size_t selected = CBMC::selectTrialPosition(random, logBoltzmannFactors);

  const MoleculeTrial &selectedTrial = externalEnergies[selected];

  double RosenbluthWeight = std::accumulate(logBoltzmannFactors.begin(), logBoltzmannFactors.end(), 0.0,
                                            [](const double &acc, const double &logBoltzmannFactor)
                                            { return acc + std::exp(logBoltzmannFactor); });

  if (RosenbluthWeight < context.forceField.minimumRosenbluthFactor) return std::nullopt;

  return ChainGrowData(selectedTrial.molecule, selectedTrial.positions, selectedTrial.energy,
                       RosenbluthWeight / double(context.forceField.numberOfTrialDirections), 0.0);
}
