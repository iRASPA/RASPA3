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
import ewald_structure_factor;
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

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const AtomDynamics> dynamics = system.spanOfMoleculeDynamics();

  // Only rigid multi-atomic molecules participate; the others keep Delta phi = 0 and a zero torque.
  const auto torques = [&]()
  {
    std::vector<double3> result(numberOfMolecules);
    std::size_t index = 0;
    for (std::size_t i = 0; i != numberOfMolecules; ++i)
    {
      const Molecule &molecule = system.moleculeData[i];
      const Component &component = system.components[molecule.componentId];
      result[i] = labFrameTorque(component, molecule, std::span<const Atom>(&atomData[index], molecule.numberOfAtoms),
                                 std::span<const AtomDynamics>(&dynamics[index], molecule.numberOfAtoms));
      index += molecule.numberOfAtoms;
    }
    return result;
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

  // Torques on every rigid molecule in the current configuration.
  system.precomputeTotalGradients();
  const std::vector<double3> torqueOld = torques();

  std::vector<double3> angularDisplacement(numberOfMolecules);
  for (std::size_t i = 0; i != numberOfMolecules; ++i)
  {
    const Molecule &molecule = system.moleculeData[i];
    const Component &component = system.components[molecule.componentId];
    if (molecule.numberOfAtoms >= 2 && component.rigid)
    {
      double3 xi(random.Gaussian(), random.Gaussian(), random.Gaussian());
      angularDisplacement[i] = b * torqueOld[i] + sigma * xi;
    }
  }

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
  };

  std::size_t index = 0;
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

  // A molecule rotated into a blocked pocket rejects the collective move.
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

  // Torques in the trial configuration, and its exact energy.
  system.precomputeTotalGradients();
  const std::vector<double3> torqueNew = torques();
  const RunningEnergy newEnergy = system.computeTotalEnergies();
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);

  system.mc_moves_statistics.addConstructed(move);

  double logBias = 0.0;
  for (std::size_t i = 0; i != numberOfMolecules; ++i)
  {
    double3 forwardDeviation = angularDisplacement[i] - b * torqueOld[i];
    double3 reverseDeviation = angularDisplacement[i] + b * torqueNew[i];
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
