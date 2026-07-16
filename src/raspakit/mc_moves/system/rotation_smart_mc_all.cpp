module;

module mc_moves_rotation_smart_mc_all;

import std;

import component;
import molecule;
import atom;
import atom_dynamics;
import double3;
import simd_quatd;
import randomnumbers;
import system;
import running_energy;
import integrators_update;
import interactions_ewald;
import mc_moves_move_types;

static double3 labFrameTorque(const Component &component, const Molecule &molecule, std::span<const Atom> atoms,
                              std::span<const AtomDynamics> dynamics)
{
  if (molecule.numberOfAtoms < 2 || !component.rigid)
  {
    return double3{};
  }

  double3 torque{};
  const double3 com = molecule.centerOfMassPosition;
  for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i)
  {
    torque += double3::cross(atoms[i].position - com, -dynamics[i].gradient);
  }
  return torque;
}

static simd_quatd quaternionFromAngularDisplacement(double3 deltaPhi)
{
  const double angle = deltaPhi.length();
  if (angle < 1.0e-14)
  {
    return simd_quatd(0.0, 0.0, 0.0, 1.0);
  }
  return simd_quatd::fromAxisAngle(angle, deltaPhi * (1.0 / angle));
}

std::optional<RunningEnergy> MC_Moves::rotationSmartMCMoveAll(RandomNumber &random, System &system)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::RotationSmartMCAll;

  system.mc_moves_statistics.addTrial(move);

  std::size_t numberOfMolecules = system.moleculeData.size();
  if (numberOfMolecules == 0)
  {
    return std::nullopt;
  }

  // Angular step size 'sigma' (radians); drift coefficient b = beta * sigma^2 / 2.
  double sigma = system.mc_moves_statistics.getMaxChange(move);
  double b = 0.5 * system.beta * sigma * sigma;

  RunningEnergy savedRunningEnergy = system.runningEnergies;

  time_begin = std::chrono::steady_clock::now();
  system.precomputeTotalGradients();
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                        system.spanOfMoleculeDynamics(), system.components);
  RunningEnergy oldEnergy = system.runningEnergies;

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();

  std::vector<double3> torqueOld(numberOfMolecules);
  std::vector<double3> angularDisplacement(numberOfMolecules);
  std::size_t index = 0;
  for (std::size_t i = 0; i != numberOfMolecules; ++i)
  {
    Molecule &molecule = system.moleculeData[i];
    const Component &component = system.components[molecule.componentId];
    std::span<const Atom> atoms(&atomData[index], molecule.numberOfAtoms);
    std::span<const AtomDynamics> moleculeDynamics(&dynamics[index], molecule.numberOfAtoms);

    // Only rigid multi-atomic molecules participate; others keep Delta phi = 0 for detailed balance.
    if (molecule.numberOfAtoms >= 2 && component.rigid)
    {
      torqueOld[i] = labFrameTorque(component, molecule, atoms, moleculeDynamics);
      double3 xi(random.Gaussian(), random.Gaussian(), random.Gaussian());
      angularDisplacement[i] = b * torqueOld[i] + sigma * xi;
    }
    else
    {
      torqueOld[i] = double3{};
      angularDisplacement[i] = double3{};
    }
    index += molecule.numberOfAtoms;
  }

  std::vector<Atom> savedAtoms(atomData.begin(), atomData.end());
  std::vector<Molecule> savedMolecules(system.moleculeData);

  index = 0;
  for (std::size_t i = 0; i != numberOfMolecules; ++i)
  {
    Molecule &molecule = system.moleculeData[i];
    const Component &component = system.components[molecule.componentId];
    std::span<Atom> atoms(&atomData[index], molecule.numberOfAtoms);

    if (molecule.numberOfAtoms >= 2 && component.rigid)
    {
      simd_quatd rotationQuaternion = quaternionFromAngularDisplacement(angularDisplacement[i]);
      auto trial = component.rotate(molecule, atoms, rotationQuaternion);
      molecule = trial.first;
      std::copy(trial.second.cbegin(), trial.second.cend(), atoms.begin());
    }
    index += molecule.numberOfAtoms;
  }

  system.precomputeTotalGradients();
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                        system.spanOfMoleculeDynamics(), system.components);
  RunningEnergy newEnergy = system.runningEnergies;
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);

  system.mc_moves_statistics.addConstructed(move);

  double logBias = 0.0;
  index = 0;
  for (std::size_t i = 0; i != numberOfMolecules; ++i)
  {
    Molecule &molecule = system.moleculeData[i];
    const Component &component = system.components[molecule.componentId];
    std::span<const Atom> atoms(&atomData[index], molecule.numberOfAtoms);
    std::span<const AtomDynamics> moleculeDynamics(&dynamics[index], molecule.numberOfAtoms);

    double3 torqueNew{};
    if (molecule.numberOfAtoms >= 2 && component.rigid)
    {
      torqueNew = labFrameTorque(component, molecule, atoms, moleculeDynamics);
    }

    double3 forwardDeviation = angularDisplacement[i] - b * torqueOld[i];
    double3 reverseDeviation = angularDisplacement[i] + b * torqueNew;
    logBias += forwardDeviation.length_squared() - reverseDeviation.length_squared();
    index += molecule.numberOfAtoms;
  }
  logBias /= (2.0 * sigma * sigma);

  RunningEnergy energyDifference = newEnergy - oldEnergy;

  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + logBias))
  {
    system.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    system.runningEnergies = savedRunningEnergy + energyDifference;
    return system.runningEnergies;
  }

  std::copy(savedAtoms.cbegin(), savedAtoms.cend(), atomData.begin());
  system.moleculeData = savedMolecules;
  system.runningEnergies = savedRunningEnergy;
  return std::nullopt;
}
