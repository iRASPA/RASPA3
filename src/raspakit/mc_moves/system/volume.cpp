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

MC_Moves::SavedVolumeState MC_Moves::saveVolumeState(const System& system)
{
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  SavedVolumeState saved{.box = system.simulationBox,
                         .molecules = system.moleculeData,
                         .atoms = std::vector<Atom>(atoms.begin(), atoms.end()),
                         .storedEik = system.storedEik,
                         .fixedFrameworkStoredEik = system.fixedFrameworkStoredEik,
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

RunningEnergy MC_Moves::applyTrialVolume(System& system, const SimulationBox& newBox)
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
  // The stored structure factor of the rigid framework belongs to the old box: its reciprocal vectors, Ewald
  // alpha and wave-vector count all change with the cell (the atoms stay where they are). The total energy
  // only computes it when it is absent, so it is rebuilt here for the trial box.
  if (!system.spanOfRigidFrameworkAtoms().empty()) system.precomputeTotalRigidEnergy();

  return system.computeTotalEnergies();
}

void MC_Moves::restoreVolumeState(System& system, const SavedVolumeState& saved)
{
  system.simulationBox = saved.box;
  system.moleculeData = saved.molecules;
  std::copy(saved.atoms.begin(), saved.atoms.end(), system.spanOfMoleculeAtoms().begin());
  system.storedEik = saved.storedEik;
  system.fixedFrameworkStoredEik = saved.fixedFrameworkStoredEik;
  if (system.forceField.computePolarization)
  {
    std::copy(saved.electricField.begin(), saved.electricField.end(), system.spanOfMoleculeElectricField().begin());
  }

  const bool vdwCutOffChanged = system.forceField.cutOffFrameworkVDW != saved.cutOffFrameworkVDW ||
                                system.forceField.cutOffMoleculeVDW != saved.cutOffMoleculeVDW;
  system.forceField.cutOffFrameworkVDW = saved.cutOffFrameworkVDW;
  system.forceField.cutOffMoleculeVDW = saved.cutOffMoleculeVDW;
  // the shifts and tail corrections follow an automatic VDW cut-off (see ForceField::initializeAutomaticCutOff)
  if (vdwCutOffChanged) system.forceField.preComputeCutOffDependentParameters();
  system.forceField.cutOffCoulomb = saved.cutOffCoulomb;
  system.forceField.EwaldAlpha = saved.ewaldAlpha;
  system.forceField.numberOfWaveVectors = saved.numberOfWaveVectors;
  system.forceField.reciprocalIntegerCutOffSquared = saved.reciprocalIntegerCutOffSquared;
  system.forceField.updateEwaldRealSpaceTable();

  // 'applyTrialVolume' rebuilt the molecule cell list (via 'computeTotalEnergies') for the trial box and
  // positions; it must not survive the restore (the box comparison in 'isCurrent' would catch a changed box,
  // but a restored cut-off or an identical box would not).
  system.invalidateCellList();
}

bool MC_Moves::volumeTrialOverlaps(const System& system, const RunningEnergy& trialEnergy)
{
  return trialEnergy.moleculeMoleculeVDW + trialEnergy.moleculeMoleculeCharge > system.forceField.energyOverlapCriteria;
}

std::optional<RunningEnergy> MC_Moves::volumeMove(RandomNumber& random, System& system)
{
  Move::Types move = Move::Types::VolumeChange;

  // Update volume move counts
  system.mc_moves_statistics.addTrial(move);

  RunningEnergy oldTotalEnergy = system.runningEnergies;
  // Every centre of mass the trial scales, including reserved fractional slots (the integer loading
  // is a different N: it is the observable, not the Jacobian exponent).
  const double numberOfMolecules = numberOfScaledMolecules(system);
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

double MC_Moves::anisotropicPressureWork(const System& system, const SimulationBox& oldBox, const SimulationBox& newBox)
{
  const double3& p = system.pressureTensorDiagonal;
  const double hydrostaticPressure = (p.x + p.y + p.z) / 3.0;

  // hydrostatic part: the exact p dV of the NpT ensemble
  const double hydrostaticWork = hydrostaticPressure * (newBox.volume - oldBox.volume);

  // deviatoric part: a constant tension -V_0 (p_aa - p) / L_a^0 on each pair of faces (linear strain with
  // respect to the reference cell); it vanishes identically for a hydrostatic pressure tensor
  const SimulationBox& reference = system.referenceSimulationBox;
  const double deviatoricWork =
      reference.volume * ((p.x - hydrostaticPressure) * (newBox.lengthA - oldBox.lengthA) / reference.lengthA +
                          (p.y - hydrostaticPressure) * (newBox.lengthB - oldBox.lengthB) / reference.lengthB +
                          (p.z - hydrostaticPressure) * (newBox.lengthC - oldBox.lengthC) / reference.lengthC);

  return hydrostaticWork + deviatoricWork;
}

std::optional<RunningEnergy> MC_Moves::anisotropicVolumeMove(RandomNumber& random, System& system)
{
  Move::Types move = Move::Types::AnisotropicVolumeChange;

  system.mc_moves_statistics.addTrial(move);

  RunningEnergy oldTotalEnergy = system.runningEnergies;
  const double numberOfMolecules = numberOfScaledMolecules(system);

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

  // The box lengths are sampled uniformly in ln L_a, so the Jacobian is (V_n/V_o)^{N+1} as for the isotropic
  // move. The pressure work is the exact difference of the constant-stress enthalpy between the two boxes.
  const double volumeRatio = scale.x * scale.y * scale.z;
  const double pressureWork = anisotropicPressureWork(system, saved.box, newBox);

  if (random.uniform() <
      std::exp((numberOfMolecules + 1.0) * std::log(volumeRatio) -
               (pressureWork + (newTotalEnergy.potentialEnergy() - oldTotalEnergy.potentialEnergy())) * system.beta))
  {
    system.mc_moves_statistics.addAccepted(move);
    return newTotalEnergy - oldTotalEnergy;
  }

  restoreVolumeState(system, saved);
  return std::nullopt;
}
