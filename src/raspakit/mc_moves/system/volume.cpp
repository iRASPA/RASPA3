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
import energy_factor;
import energy_status;
import energy_status_inter;
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::volumeMove(RandomNumber &random, System &system)
{
  std::chrono::system_clock::time_point time_begin, time_end;
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
  std::pair<std::vector<Molecule>, std::vector<Atom>> newPositions = system.scaledCenterOfMassPositions(scale);

  double cutOffFrameworkVDW_stored = system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW_stored = system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb_stored = system.forceField.cutOffCoulomb;
  double ewald_alpha_stored = system.forceField.EwaldAlpha;
  int3 ewald_k_stored = system.forceField.numberOfWaveVectors;

  system.forceField.initializeAutomaticCutOff(newBox);

  time_begin = std::chrono::system_clock::now();
  // Compute new intermolecular energy
  RunningEnergy newTotalInterEnergy =
      Interactions::computeInterMolecularEnergy(system.forceField, newBox, newPositions.second);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  // Compute new tail corrections
  RunningEnergy newTotalTailEnergy =
      Interactions::computeInterMolecularTailEnergy(system.forceField, newBox, newPositions.second);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  // Compute new Ewald Fourier energy
  RunningEnergy newTotalEwaldEnergy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.totalEik,
      system.forceField, newBox, system.components, system.numberOfMoleculesPerComponent, newPositions.second,
      system.netChargeFramework);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Sum up all energy contributions
  RunningEnergy newTotalEnergy = newTotalInterEnergy + newTotalTailEnergy + newTotalEwaldEnergy;

  // The intra-molecular energies have not changed by the com-scaling
  newTotalEnergy.bond = oldTotalEnergy.bond;
  newTotalEnergy.ureyBradley = oldTotalEnergy.ureyBradley;
  newTotalEnergy.bend = oldTotalEnergy.bend;
  newTotalEnergy.inversionBend = oldTotalEnergy.inversionBend;
  newTotalEnergy.outOfPlaneBend = oldTotalEnergy.outOfPlaneBend;
  newTotalEnergy.torsion = oldTotalEnergy.torsion;
  newTotalEnergy.improperTorsion = oldTotalEnergy.improperTorsion;
  newTotalEnergy.bondBond = oldTotalEnergy.bondBond;
  newTotalEnergy.bondBend = oldTotalEnergy.bondBend;
  newTotalEnergy.bondTorsion = oldTotalEnergy.bondTorsion;
  newTotalEnergy.bendBend = oldTotalEnergy.bendBend;
  newTotalEnergy.bendTorsion = oldTotalEnergy.bendTorsion;

  // Update constructed move counts
  system.mc_moves_statistics.addConstructed(move);

  // Apply acceptance/rejection rule
  if (random.uniform() < std::exp((numberOfMolecules + 1.0) * std::log(newVolume / oldVolume) -
                                  (system.pressure * (newVolume - oldVolume) +
                                   (newTotalEnergy.potentialEnergy() - oldTotalEnergy.potentialEnergy())) *
                                      system.beta))
  {
    // Move accepted: update system state
    system.mc_moves_statistics.addAccepted(move);

    system.simulationBox = newBox;
    std::copy(newPositions.first.begin(), newPositions.first.end(), system.moleculeData.begin());
    std::copy(newPositions.second.begin(), newPositions.second.end(), system.atomData.begin());

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

    return newTotalEnergy - oldTotalEnergy;
  }

  system.forceField.cutOffFrameworkVDW = cutOffFrameworkVDW_stored;
  system.forceField.cutOffMoleculeVDW = cutOffMoleculeVDW_stored;
  system.forceField.cutOffCoulomb = cutOffCoulomb_stored;
  system.forceField.EwaldAlpha = ewald_alpha_stored;
  system.forceField.numberOfWaveVectors = ewald_k_stored;

  return std::nullopt;
}

std::optional<RunningEnergy> MC_Moves::anisotropicVolumeMove(RandomNumber& random, System& system)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::AnisotropicVolumeChange;

  system.mc_moves_statistics.addTrial(move);

  RunningEnergy oldTotalEnergy = system.runningEnergies;
  const double numberOfMolecules = static_cast<double>(std::accumulate(
      system.numberOfIntegerMoleculesPerComponent.begin(), system.numberOfIntegerMoleculesPerComponent.end(), 0));
  const SimulationBox& oldBox = system.simulationBox;
  const double oldVolume = oldBox.volume;

  const double3 maxVolumeChange(
      system.mc_moves_statistics.getMaxChange(move, 0),
      system.mc_moves_statistics.getMaxChange(move, 1),
      system.mc_moves_statistics.getMaxChange(move, 2));

  const double3 scale(std::exp(maxVolumeChange.x * (2.0 * random.uniform() - 1.0)),
                      std::exp(maxVolumeChange.y * (2.0 * random.uniform() - 1.0)),
                      std::exp(maxVolumeChange.z * (2.0 * random.uniform() - 1.0)));

  SimulationBox newBox = oldBox.scaled(scale);
  std::pair<std::vector<Molecule>, std::vector<Atom>> newPositions{system.moleculeData, system.atomData};
  for (Molecule& molecule : newPositions.first)
  {
    std::span<Atom> moleculeAtoms = {&newPositions.second[molecule.atomIndex], molecule.numberOfAtoms};

    for (Atom& atom : moleculeAtoms)
    {
      const double3 fractional = oldBox.inverseCell * atom.position;
      atom.position = newBox.cell * fractional;
    }

    double totalMass = 0.0;
    double3 newCom(0.0, 0.0, 0.0);
    for (const Atom& atom : moleculeAtoms)
    {
      const double mass = system.forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].mass;
      newCom += mass * atom.position;
      totalMass += mass;
    }
    molecule.centerOfMassPosition = newCom / totalMass;
  }

  const double cutOffFrameworkVDW_stored = system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW_stored = system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb_stored = system.forceField.cutOffCoulomb;
  const double ewald_alpha_stored = system.forceField.EwaldAlpha;
  const int3 ewald_k_stored = system.forceField.numberOfWaveVectors;

  system.forceField.initializeAutomaticCutOff(newBox);

  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalInterEnergy =
      Interactions::computeInterMolecularEnergy(system.forceField, newBox, newPositions.second);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalTailEnergy =
      Interactions::computeInterMolecularTailEnergy(system.forceField, newBox, newPositions.second);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalEwaldEnergy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.totalEik,
      system.forceField, newBox, system.components, system.numberOfMoleculesPerComponent, newPositions.second,
      system.netChargeFramework);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  RunningEnergy newTotalEnergy = newTotalInterEnergy + newTotalTailEnergy + newTotalEwaldEnergy;

  newTotalEnergy.bond = oldTotalEnergy.bond;
  newTotalEnergy.ureyBradley = oldTotalEnergy.ureyBradley;
  newTotalEnergy.bend = oldTotalEnergy.bend;
  newTotalEnergy.inversionBend = oldTotalEnergy.inversionBend;
  newTotalEnergy.outOfPlaneBend = oldTotalEnergy.outOfPlaneBend;
  newTotalEnergy.torsion = oldTotalEnergy.torsion;
  newTotalEnergy.improperTorsion = oldTotalEnergy.improperTorsion;
  newTotalEnergy.bondBond = oldTotalEnergy.bondBond;
  newTotalEnergy.bondBend = oldTotalEnergy.bondBend;
  newTotalEnergy.bondTorsion = oldTotalEnergy.bondTorsion;
  newTotalEnergy.bendBend = oldTotalEnergy.bendBend;
  newTotalEnergy.bendTorsion = oldTotalEnergy.bendTorsion;

  system.mc_moves_statistics.addConstructed(move);

  const double volumeRatio = scale.x * scale.y * scale.z;
  const double pressureWork =
      oldVolume * (system.pressureTensorDiagonal.x * (scale.x - 1.0) +
                   system.pressureTensorDiagonal.y * (scale.y - 1.0) +
                   system.pressureTensorDiagonal.z * (scale.z - 1.0));

  if (random.uniform() < std::exp((numberOfMolecules + 1.0) * std::log(volumeRatio) -
                                  (pressureWork + (newTotalEnergy.potentialEnergy() - oldTotalEnergy.potentialEnergy())) *
                                      system.beta))
  {
    system.mc_moves_statistics.addAccepted(move);

    system.simulationBox = newBox;
    std::copy(newPositions.first.begin(), newPositions.first.end(), system.moleculeData.begin());
    std::copy(newPositions.second.begin(), newPositions.second.end(), system.atomData.begin());

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

    return newTotalEnergy - oldTotalEnergy;
  }

  system.forceField.cutOffFrameworkVDW = cutOffFrameworkVDW_stored;
  system.forceField.cutOffMoleculeVDW = cutOffMoleculeVDW_stored;
  system.forceField.cutOffCoulomb = cutOffCoulomb_stored;
  system.forceField.EwaldAlpha = ewald_alpha_stored;
  system.forceField.numberOfWaveVectors = ewald_k_stored;

  return std::nullopt;
}
