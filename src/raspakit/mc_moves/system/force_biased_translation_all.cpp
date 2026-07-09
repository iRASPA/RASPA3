module;

module mc_moves_force_biased_translation_all;

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

std::optional<RunningEnergy> MC_Moves::forceBiasTranslationMoveAll(RandomNumber &random, System &system)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::ForceBiasTranslationAll;

  system.mc_moves_statistics.addTrial(move);

  std::size_t numberOfMolecules = system.moleculeData.size();
  if (numberOfMolecules == 0)
  {
    return std::nullopt;
  }

  // The step size 'sigma' is the standard deviation of the Gaussian part of each per-molecule
  // displacement; the drift coefficient follows the smart-MC relation b = beta * sigma^2 / 2.
  double sigma = system.mc_moves_statistics.getMaxChange(move);
  double b = 0.5 * system.beta * sigma * sigma;

  // Preserve the maintained running energy so that tail/polarization baselines are kept: only the
  // gradient-based potential difference is applied to it on acceptance.
  RunningEnergy savedRunningEnergy = system.runningEnergies;

  // Forces on every molecule in the current configuration.
  time_begin = std::chrono::steady_clock::now();
  system.precomputeTotalGradients();
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                        system.spanOfMoleculeDynamics(), system.components);
  RunningEnergy oldEnergy = system.runningEnergies;

  // Biased trial displacement of each molecule's center of mass.
  std::vector<double3> forceOld(numberOfMolecules);
  std::vector<double3> displacement(numberOfMolecules);
  for (std::size_t i = 0; i != numberOfMolecules; ++i)
  {
    forceOld[i] = -system.moleculeData[i].gradient;
    double3 xi(random.Gaussian(), random.Gaussian(), random.Gaussian());
    displacement[i] = b * forceOld[i] + sigma * xi;
  }

  // Save the current configuration and apply the collective trial displacement to the live system.
  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::vector<Atom> savedAtoms(atomData.begin(), atomData.end());
  std::vector<Molecule> savedMolecules(system.moleculeData);

  std::size_t index = 0;
  for (std::size_t i = 0; i != numberOfMolecules; ++i)
  {
    Molecule &molecule = system.moleculeData[i];
    molecule.centerOfMassPosition += displacement[i];
    std::span<Atom> span(&atomData[index], molecule.numberOfAtoms);
    for (Atom &atom : span)
    {
      atom.position += displacement[i];
    }
    index += molecule.numberOfAtoms;
  }

  // Forces on every molecule in the trial configuration.
  system.precomputeTotalGradients();
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                        system.spanOfMoleculeDynamics(), system.components);
  RunningEnergy newEnergy = system.runningEnergies;
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);

  system.mc_moves_statistics.addConstructed(move);

  // Metropolis-Hastings acceptance with the asymmetric-proposal correction summed over all molecules.
  double logBias = 0.0;
  for (std::size_t i = 0; i != numberOfMolecules; ++i)
  {
    double3 forceNew = -system.moleculeData[i].gradient;
    double3 forwardDeviation = displacement[i] - b * forceOld[i];
    double3 reverseDeviation = displacement[i] + b * forceNew;
    logBias += forwardDeviation.length_squared() - reverseDeviation.length_squared();
  }
  logBias /= (2.0 * sigma * sigma);

  RunningEnergy energyDifference = newEnergy - oldEnergy;

  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + logBias))
  {
    system.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

    // Add only the gradient-based potential difference to the maintained running energy.
    system.runningEnergies = savedRunningEnergy + energyDifference;
    return system.runningEnergies;
  }

  // Reject: restore the original configuration and running energy.
  std::copy(savedAtoms.cbegin(), savedAtoms.cend(), atomData.begin());
  system.moleculeData = savedMolecules;
  system.runningEnergies = savedRunningEnergy;
  return std::nullopt;
}
