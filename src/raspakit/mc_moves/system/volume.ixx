module;

export module mc_moves_volume;

import std;

import int3;
import double3;
import randomnumbers;
import running_energy;
import atom;
import molecule;
import simulationbox;
import ewald_structure_factor;
import system;

export namespace MC_Moves
{
/**
 * \brief The state of a system a volume move changes, saved before the trial box is applied and restored
 * when the move is rejected (see 'applyTrialVolume').
 */
struct SavedVolumeState
{
  SimulationBox box;
  std::vector<Molecule> molecules;
  std::vector<Atom> atoms;
  EwaldStructureFactor storedEik;
  std::vector<double3> electricField;  ///< Only when the force field computes polarization.
  double cutOffFrameworkVDW;
  double cutOffMoleculeVDW;
  double cutOffCoulomb;
  double ewaldAlpha;
  int3 numberOfWaveVectors;
  std::size_t reciprocalIntegerCutOffSquared;
};

SavedVolumeState saveVolumeState(const System &system);

/**
 * \brief Applies a trial box to the system and returns the exact total energy of the resulting configuration.
 *
 * The molecules are translated rigidly so that their centres of mass keep their fractional positions
 * ('System::scaledCenterOfMassPositions'); the automatic cut-offs and the Ewald parameters are re-derived for
 * the new box ('initializeAutomaticCutOff' followed by 'initializeEwaldParameters', as every driver does at
 * stage entry: with a user-fixed Coulomb cut-off the number of wave vectors still follows the box). The energy
 * is a full recomputation ('System::computeTotalEnergies'), so it carries every term the running energy
 * carries -- inter-molecular, tail, Ewald, intra-molecular, cross-links, polymer links, external field and
 * polarization -- and the stored Ewald structure factors and electric field describe the trial configuration
 * afterwards. A rejected move restores all of that with 'restoreVolumeState'.
 */
RunningEnergy applyTrialVolume(System &system, const SimulationBox &newBox);

void restoreVolumeState(System &system, const SavedVolumeState &saved);

/// Whether the inter-molecular part of a trial energy marks an overlap (the move is then rejected outright).
bool volumeTrialOverlaps(const System &system, const RunningEnergy &trialEnergy);

/**
 * \brief Performs a volume move in a Monte Carlo simulation.
 *
 * Attempts to change the system volume by a random amount, scaling the simulation box and molecule positions
 * accordingly. Computes the new energies and decides whether to accept the move based on the Metropolis criterion.
 *
 * \param random A random number generator.
 * \param system The simulation system to modify.
 * \return The new total energy if the move is accepted; std::nullopt otherwise.
 */
std::optional<RunningEnergy> volumeMove(RandomNumber &random, System &system);

/**
 * \brief Performs an anisotropic volume move with independent pressures along a, b, and c.
 *
 * Scales each lattice vector length independently and applies the corresponding diagonal pressure component
 * in the acceptance criterion.
 */
std::optional<RunningEnergy> anisotropicVolumeMove(RandomNumber &random, System &system);
}  // namespace MC_Moves
