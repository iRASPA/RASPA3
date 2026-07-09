module;

module mc_moves_force_biased_translation;

import std;

import component;
import molecule;
import atom;
import double3;
import randomnumbers;
import system;
import running_energy;
import integrators_update;
import interactions_ewald;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::forceBiasTranslationMove(RandomNumber &random, System &system,
                                                               std::size_t selectedComponent,
                                                               [[maybe_unused]] std::size_t selectedMolecule,
                                                               const std::vector<Component> &components,
                                                               Molecule &molecule, std::span<Atom> molecule_atoms)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::ForceBiasTranslation;
  Component &component = system.components[selectedComponent];

  // Register the trial. The step size 'sigma' is the standard deviation of the Gaussian part of the
  // trial displacement; the drift coefficient follows the smart-MC relation b = beta * sigma^2 / 2.
  component.mc_moves_statistics.addTrial(move);

  double sigma = component.mc_moves_statistics.getMaxChange(move);
  double b = 0.5 * system.beta * sigma * sigma;

  // Keep the true (maintained) running energy so that tail/polarization baselines are preserved: only the
  // gradient-based potential difference is added to it on acceptance.
  RunningEnergy savedRunningEnergy = system.runningEnergies;

  // Force on the selected molecule in the current configuration.
  time_begin = std::chrono::system_clock::now();
  system.precomputeTotalGradients();
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                        system.spanOfMoleculeDynamics(), system.components);
  RunningEnergy oldEnergy = system.runningEnergies;
  double3 forceOld = -molecule.gradient;

  // Biased trial displacement of the center of mass.
  double3 xi(random.Gaussian(), random.Gaussian(), random.Gaussian());
  double3 displacement = b * forceOld + sigma * xi;

  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      components[selectedComponent].translate(molecule, molecule_atoms, displacement);

  // Reject when the trial position lies inside a blocked pocket. The configuration was not modified yet.
  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    system.runningEnergies = savedRunningEnergy;
    return std::nullopt;
  }

  // Apply the trial configuration to the live system so the force at the new position can be evaluated.
  std::vector<Atom> savedAtoms(molecule_atoms.begin(), molecule_atoms.end());
  Molecule savedMolecule = molecule;
  std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), molecule_atoms.begin());
  molecule = trialMolecule.first;

  // Force on the selected molecule in the trial configuration.
  system.precomputeTotalGradients();
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                        system.spanOfMoleculeDynamics(), system.components);
  RunningEnergy newEnergy = system.runningEnergies;
  double3 forceNew = -molecule.gradient;
  time_end = std::chrono::system_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);

  RunningEnergy energyDifference = newEnergy - oldEnergy;

  component.mc_moves_statistics.addConstructed(move);

  // Metropolis-Hastings acceptance with the asymmetric-proposal correction. The forward proposal draws
  // Delta r ~ N(b * F_old, sigma^2); the reverse proposal (n -> o) draws -Delta r ~ N(b * F_new, sigma^2).
  double3 forwardDeviation = displacement - b * forceOld;
  double3 reverseDeviation = displacement + b * forceNew;
  double logBias = (forwardDeviation.length_squared() - reverseDeviation.length_squared()) / (2.0 * sigma * sigma);

  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + logBias))
  {
    component.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

    // Restore the maintained running energy; the caller adds the gradient-based potential difference.
    system.runningEnergies = savedRunningEnergy;
    return energyDifference;
  }

  // Reject: restore the original configuration and running energy.
  std::copy(savedAtoms.cbegin(), savedAtoms.cend(), molecule_atoms.begin());
  molecule = savedMolecule;
  system.runningEnergies = savedRunningEnergy;
  return std::nullopt;
}
