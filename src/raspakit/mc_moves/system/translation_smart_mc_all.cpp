module;

module mc_moves_translation_smart_mc_all;

import std;

import component;
import molecule;
import atom;
import atom_dynamics;
import double3;
import randomnumbers;
import system;
import running_energy;
import ewald_structure_factor;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::translationSmartMCMoveAll(RandomNumber &random, System &system)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::TranslationSmartMCAll;

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

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const AtomDynamics> dynamics = system.spanOfMoleculeDynamics();

  // The net force on a molecule is the sum of its atomic forces whatever the molecule type (rigid,
  // flexible or semi-flexible); 'Molecule::gradient' is not maintained for semi-flexible molecules
  // (their rigid groups carry the gradients), so the sum is taken here directly.
  const auto netForces = [&]()
  {
    std::vector<double3> forces(numberOfMolecules);
    std::size_t index = 0;
    for (std::size_t i = 0; i != numberOfMolecules; ++i)
    {
      double3 gradient{};
      for (std::size_t k = 0; k != system.moleculeData[i].numberOfAtoms; ++k)
      {
        gradient += dynamics[index + k].gradient;
      }
      forces[i] = -gradient;
      index += system.moleculeData[i].numberOfAtoms;
    }
    return forces;
  };

  // The gradient routines cover the pair, Ewald and intra-molecular terms only; the acceptance uses the
  // exact total energies (external field and polarization included) of the old and the new configuration.
  // 'computeTotalEnergies' also refreshes the stored Ewald structure factors and the electric field for the
  // configuration it evaluates, so those are saved for the rejection path.
  time_begin = std::chrono::steady_clock::now();
  const RunningEnergy savedRunningEnergy = system.runningEnergies;
  const RunningEnergy oldEnergy = system.computeTotalEnergies();
  const EwaldStructureFactor savedStoredEik(system.storedEik);
  std::vector<double3> savedElectricField;
  if (system.forceField.computePolarization)
  {
    const std::span<const double3> electricField = system.spanOfMoleculeElectricField();
    savedElectricField.assign(electricField.begin(), electricField.end());
  }

  // Forces on every molecule in the current configuration (the group states of the semi-flexible molecules
  // rebuilt first: the rigid-group torques are taken about their centres of mass).
  system.initializeGroupData();
  system.precomputeTotalGradients();
  const std::vector<double3> forceOld = netForces();
  // The gradients of the old configuration, the framework's and the group states (rebuilt above from the positions,
  // which the Monte Carlo moves do not keep them in step with) go back on rejection: a refused trial leaves the
  // system as it was, gradients and all.
  const std::vector<AtomDynamics> savedDynamics(dynamics.begin(), dynamics.end());
  const std::span<const AtomDynamics> frameworkDynamics = system.spanOfFrameworkDynamics();
  const std::vector<AtomDynamics> savedFrameworkDynamics(frameworkDynamics.begin(), frameworkDynamics.end());
  const std::vector<GroupState> savedGroupData(system.groupData);

  // Biased trial displacement of each molecule's center of mass.
  std::vector<double3> displacement(numberOfMolecules);
  for (std::size_t i = 0; i != numberOfMolecules; ++i)
  {
    double3 xi(random.Gaussian(), random.Gaussian(), random.Gaussian());
    displacement[i] = b * forceOld[i] + sigma * xi;
  }

  // Save the current configuration and apply the collective trial displacement to the live system.
  std::vector<Atom> savedAtoms(atomData.begin(), atomData.end());
  std::vector<Molecule> savedMolecules(system.moleculeData);

  const auto restore = [&]()
  {
    std::copy(savedAtoms.cbegin(), savedAtoms.cend(), atomData.begin());
    system.moleculeData = savedMolecules;
    system.storedEik = savedStoredEik;
    if (system.forceField.computePolarization)
    {
      std::copy(savedElectricField.cbegin(), savedElectricField.cend(), system.spanOfMoleculeElectricField().begin());
    }
    system.runningEnergies = savedRunningEnergy;
    std::copy(savedDynamics.cbegin(), savedDynamics.cend(), system.spanOfMoleculeDynamics().begin());
    std::copy(savedFrameworkDynamics.cbegin(), savedFrameworkDynamics.cend(), system.spanOfFrameworkDynamics().begin());
    system.groupData = savedGroupData;
    // the molecule cell list may have been rebuilt by 'computeTotalEnergies' for the rejected trial positions
    system.invalidateCellList();
  };

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

  // A molecule displaced into a blocked pocket rejects the collective move.
  index = 0;
  for (std::size_t i = 0; i != numberOfMolecules; ++i)
  {
    const Molecule &molecule = system.moleculeData[i];
    if (system.insideBlockedPockets(system.components[molecule.componentId],
                                    std::span<const Atom>(&atomData[index], molecule.numberOfAtoms)))
    {
      restore();
      time_end = std::chrono::steady_clock::now();
      system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
      return std::nullopt;
    }
    index += molecule.numberOfAtoms;
  }

  // Forces on every molecule in the trial configuration (group states from the trial positions), and its exact
  // energy.
  system.initializeGroupData();
  system.precomputeTotalGradients();
  const std::vector<double3> forceNew = netForces();
  const RunningEnergy newEnergy = system.computeTotalEnergies();
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);

  system.mc_moves_statistics.addConstructed(move);

  // Metropolis-Hastings acceptance with the asymmetric-proposal correction summed over all molecules.
  double logBias = 0.0;
  for (std::size_t i = 0; i != numberOfMolecules; ++i)
  {
    double3 forwardDeviation = displacement[i] - b * forceOld[i];
    double3 reverseDeviation = displacement[i] + b * forceNew[i];
    logBias += forwardDeviation.length_squared() - reverseDeviation.length_squared();
  }
  logBias /= (2.0 * sigma * sigma);

  const RunningEnergy energyDifference = newEnergy - oldEnergy;

  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + logBias))
  {
    system.mc_moves_statistics.addAccepted(move);

    // The stored Ewald structure factors and the electric field already describe the new configuration.
    system.runningEnergies = savedRunningEnergy + energyDifference;
    return system.runningEnergies;
  }

  // Reject: restore the original configuration, structure factors, electric field and running energy.
  restore();
  return std::nullopt;
}
