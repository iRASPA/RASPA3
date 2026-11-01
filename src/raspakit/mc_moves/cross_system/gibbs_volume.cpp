module;

module mc_moves_gibbs_volume;

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
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;
import mc_moves_cputime;
import mc_moves_volume;

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsVolumeMove(
    RandomNumber& random, System& systemA, System& systemB, [[maybe_unused]] std::size_t selectedComponent)
{
  Move::Types move = Move::Types::GibbsVolume;

  // Volume moves do not change either particle-number macrostate and leave the
  // TMMC collection matrices untouched.
  systemA.mc_moves_statistics.addTrial(move);
  systemB.mc_moves_statistics.addTrial(move);

  // determine New box-volumes leaving the total volume constant
  double oldVolumeA = systemA.simulationBox.volume;
  double maxVolumeChangeA = systemA.mc_moves_statistics.getMaxChange(move);
  double oldVolumeB = systemB.simulationBox.volume;
  double totalVolume = oldVolumeA + oldVolumeB;
  double expdv = std::exp(std::log(oldVolumeA / oldVolumeB) + maxVolumeChangeA * (2.0 * random.uniform() - 1.0));
  double newVolumeA = expdv * totalVolume / (1.0 + expdv);
  double newVolumeB = totalVolume - newVolumeA;

  // Scale systemA to the new volume: the trial box is applied to the live system and its exact total energy
  // computed (see 'applyTrialVolume'); every non-accepting exit restores the saved state.
  RunningEnergy oldTotalEnergyA = systemA.runningEnergies;
  const double numberOfMoleculesA = numberOfScaledMolecules(systemA);
  double scaleA = std::pow(newVolumeA / oldVolumeA, 1.0 / 3.0);
  SimulationBox newBoxA = systemA.simulationBox.scaled(scaleA);

  const SavedVolumeState savedA = saveVolumeState(systemA);
  RunningEnergy newTotalEnergyA =
      timed(systemA, move, Move::Timing::NonEwald, [&] { return applyTrialVolume(systemA, newBoxA); });

  if (volumeTrialOverlaps(systemA, newTotalEnergyA))
  {
    restoreVolumeState(systemA, savedA);
    return std::nullopt;
  }

  systemA.mc_moves_statistics.addConstructed(move);

  // Scale systemB to the new volume
  RunningEnergy oldTotalEnergyB = systemB.runningEnergies;
  const double numberOfMoleculesB = numberOfScaledMolecules(systemB);
  double scaleB = std::pow(newVolumeB / oldVolumeB, 1.0 / 3.0);
  SimulationBox newBoxB = systemB.simulationBox.scaled(scaleB);

  const SavedVolumeState savedB = saveVolumeState(systemB);
  RunningEnergy newTotalEnergyB =
      timed(systemA, move, Move::Timing::NonEwald, [&] { return applyTrialVolume(systemB, newBoxB); });

  if (volumeTrialOverlaps(systemB, newTotalEnergyB))
  {
    restoreVolumeState(systemA, savedA);
    restoreVolumeState(systemB, savedB);
    return std::nullopt;
  }

  systemB.mc_moves_statistics.addConstructed(move);

  // Calculate total energy change deltaU
  double deltaU = (newTotalEnergyA.potentialEnergy() - oldTotalEnergyA.potentialEnergy()) +
                  (newTotalEnergyB.potentialEnergy() - oldTotalEnergyB.potentialEnergy());

  // apply acceptance/rejection rule
  if (random.uniform() <
      std::exp(-systemA.beta * deltaU + ((numberOfMoleculesA + 1.0) * std::log(newVolumeA / oldVolumeA)) +
               ((numberOfMoleculesB + 1.0) * std::log(newVolumeB / oldVolumeB))))
  {
    // Accept the move: both systems already hold their new boxes, positions, structure factors and electric fields.
    systemA.mc_moves_statistics.addAccepted(move);
    systemB.mc_moves_statistics.addAccepted(move);

    return std::make_pair(newTotalEnergyA - oldTotalEnergyA, newTotalEnergyB - oldTotalEnergyB);
  }

  restoreVolumeState(systemA, savedA);
  restoreVolumeState(systemB, savedB);

  return std::nullopt;
}
