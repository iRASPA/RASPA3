module;

module mc_moves_volume;

import std;

import component;
import atom;
import molecule;
import int3;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_status;
import energy_status_inter;
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import forcefield;
import ewald_structure_factor;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import mc_moves_move_types;
import mc_moves_cputime;

MC_Moves::SavedVolumeState MC_Moves::saveVolumeState(const System &system)
{
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  SavedVolumeState saved{.box = system.simulationBox,
                         .molecules = system.moleculeData,
                         .atoms = std::vector<Atom>(atoms.begin(), atoms.end()),
                         .storedEik = system.storedEik,
                         .electricField = {},
                         .cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW,
                         .cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW,
                         .cutOffCoulomb = system.forceField.cutOffCoulomb,
                         .ewaldAlpha = system.forceField.EwaldAlpha,
                         .numberOfWaveVectors = system.forceField.numberOfWaveVectors,
                         .reciprocalIntegerCutOffSquared = system.forceField.reciprocalIntegerCutOffSquared};
  if (system.forceField.computePolarization)
  {
    const std::span<const double3> electricField = system.spanOfMoleculeElectricField();
    saved.electricField.assign(electricField.begin(), electricField.end());
  }
  return saved;
}

RunningEnergy MC_Moves::applyTrialVolume(System &system, const SimulationBox &newBox)
{
  // Scale molecules by their centre of mass: only the centre-of-mass fractional position is kept fixed while
  // each molecule is translated rigidly. This preserves the internal geometry of both rigid and flexible
  // molecules. The framework atoms are untouched.
  std::pair<std::vector<Molecule>, std::vector<Atom>> newPositions =
      system.scaledCenterOfMassPositions(system.simulationBox, newBox);

  system.simulationBox = newBox;
  std::copy(newPositions.first.begin(), newPositions.first.end(), system.moleculeData.begin());
  std::copy(newPositions.second.begin(), newPositions.second.end(), system.spanOfMoleculeAtoms().begin());

  system.forceField.initializeAutomaticCutOff(newBox);
  system.forceField.initializeEwaldParameters(newBox);

  return system.computeTotalEnergies();
}

void MC_Moves::restoreVolumeState(System &system, const SavedVolumeState &saved)
{
  system.simulationBox = saved.box;
  system.moleculeData = saved.molecules;
  std::copy(saved.atoms.begin(), saved.atoms.end(), system.spanOfMoleculeAtoms().begin());
  system.storedEik = saved.storedEik;
  if (system.forceField.computePolarization)
  {
    std::copy(saved.electricField.begin(), saved.electricField.end(), system.spanOfMoleculeElectricField().begin());
  }

  system.forceField.cutOffFrameworkVDW = saved.cutOffFrameworkVDW;
  system.forceField.cutOffMoleculeVDW = saved.cutOffMoleculeVDW;
  system.forceField.cutOffCoulomb = saved.cutOffCoulomb;
  system.forceField.EwaldAlpha = saved.ewaldAlpha;
  system.forceField.numberOfWaveVectors = saved.numberOfWaveVectors;
  system.forceField.reciprocalIntegerCutOffSquared = saved.reciprocalIntegerCutOffSquared;
  system.forceField.updateEwaldRealSpaceTable();
}

bool MC_Moves::volumeTrialOverlaps(const System &system, const RunningEnergy &trialEnergy)
{
  return trialEnergy.moleculeMoleculeVDW + trialEnergy.moleculeMoleculeCharge > system.forceField.energyOverlapCriteria;
}

std::optional<RunningEnergy> MC_Moves::volumeMove(RandomNumber &random, System &system)
{
  Move::Types move = Move::Types::VolumeChange;

  // Update volume move counts
  system.mc_moves_statistics.addTrial(move);

  RunningEnergy oldTotalEnergy = system.runningEnergies;
  // Calculate the total number of molecules
  double numberOfMolecules = static_cast<double>(std::accumulate(system.numberOfIntegerMoleculesPerComponent.begin(),
                                                                 system.numberOfIntegerMoleculesPerComponent.end(), 0));
  double oldVolume = system.simulationBox.volume;
  double maxVolumeChange = system.mc_moves_statistics.getMaxChange(move);

  // Propose a new volume change
  double newVolume = std::exp(std::log(oldVolume) + maxVolumeChange * (2.0 * random.uniform() - 1.0));

  // Compute scaling factor for box dimensions
  double scale = std::pow(newVolume / oldVolume, 1.0 / 3.0);

  SimulationBox newBox = system.simulationBox.scaled(scale);

  // The trial box is applied to the live system and its exact total energy computed; every non-accepting exit
  // restores the saved state.
  const SavedVolumeState saved = saveVolumeState(system);
  RunningEnergy newTotalEnergy =
      timed(system, move, Move::Timing::NonEwald, [&] { return applyTrialVolume(system, newBox); });

  if (volumeTrialOverlaps(system, newTotalEnergy))
  {
    restoreVolumeState(system, saved);
    return std::nullopt;
  }

  // Update constructed move counts
  system.mc_moves_statistics.addConstructed(move);

  // Apply acceptance/rejection rule
  if (random.uniform() < std::exp((numberOfMolecules + 1.0) * std::log(newVolume / oldVolume) -
                                  (system.pressure * (newVolume - oldVolume) +
                                   (newTotalEnergy.potentialEnergy() - oldTotalEnergy.potentialEnergy())) *
                                      system.beta))
  {
    // Move accepted: the system already holds the new box, positions, structure factors and electric field.
    system.mc_moves_statistics.addAccepted(move);
    return newTotalEnergy - oldTotalEnergy;
  }

  restoreVolumeState(system, saved);
  return std::nullopt;
}

std::optional<RunningEnergy> MC_Moves::anisotropicVolumeMove(RandomNumber &random, System &system)
{
  Move::Types move = Move::Types::AnisotropicVolumeChange;

  system.mc_moves_statistics.addTrial(move);

  RunningEnergy oldTotalEnergy = system.runningEnergies;
  const double numberOfMolecules = static_cast<double>(std::accumulate(
      system.numberOfIntegerMoleculesPerComponent.begin(), system.numberOfIntegerMoleculesPerComponent.end(), 0));
  const double oldVolume = system.simulationBox.volume;

  const double3 maxVolumeChange(system.mc_moves_statistics.getMaxChange(move, 0),
                                system.mc_moves_statistics.getMaxChange(move, 1),
                                system.mc_moves_statistics.getMaxChange(move, 2));

  const double3 scale(std::exp(maxVolumeChange.x * (2.0 * random.uniform() - 1.0)),
                      std::exp(maxVolumeChange.y * (2.0 * random.uniform() - 1.0)),
                      std::exp(maxVolumeChange.z * (2.0 * random.uniform() - 1.0)));

  SimulationBox newBox = system.simulationBox.scaled(scale);

  const SavedVolumeState saved = saveVolumeState(system);
  RunningEnergy newTotalEnergy =
      timed(system, move, Move::Timing::NonEwald, [&] { return applyTrialVolume(system, newBox); });

  if (volumeTrialOverlaps(system, newTotalEnergy))
  {
    restoreVolumeState(system, saved);
    return std::nullopt;
  }

  system.mc_moves_statistics.addConstructed(move);

  const double volumeRatio = scale.x * scale.y * scale.z;
  const double pressureWork =
      oldVolume * (system.pressureTensorDiagonal.x * (scale.x - 1.0) + system.pressureTensorDiagonal.y * (scale.y - 1.0) +
                   system.pressureTensorDiagonal.z * (scale.z - 1.0));

  if (random.uniform() < std::exp((numberOfMolecules + 1.0) * std::log(volumeRatio) -
                                  (pressureWork + (newTotalEnergy.potentialEnergy() - oldTotalEnergy.potentialEnergy())) *
                                      system.beta))
  {
    system.mc_moves_statistics.addAccepted(move);
    return newTotalEnergy - oldTotalEnergy;
  }

  restoreVolumeState(system, saved);
  return std::nullopt;
}
