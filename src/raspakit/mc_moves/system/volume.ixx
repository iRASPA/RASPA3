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

SavedVolumeState saveVolumeState(const System& system);

/// Centres of mass a volume move scales: every record in 'moleculeData', including reserved
/// fractional slots. The configuration-space Jacobian is (V_new / V_old)^n with this n; loadings
/// and TMMC keep using the integer count.
[[nodiscard]] inline double numberOfScaledMolecules(const System& system)
{
  return static_cast<double>(system.moleculeData.size());
}

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
RunningEnergy applyTrialVolume(System& system, const SimulationBox& newBox);

void restoreVolumeState(System& system, const SavedVolumeState& saved);

/// Whether the inter-molecular part of a trial energy marks an overlap (the move is then rejected outright).
bool volumeTrialOverlaps(const System& system, const RunningEnergy& trialEnergy);

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
std::optional<RunningEnergy> volumeMove(RandomNumber& random, System& system);

/**
 * \brief The work of the imposed pressure tensor for a change of the box lengths in the anisotropic volume move.
 *
 * The move samples the constant-stress enthalpy (Parrinello and Rahman, J. Appl. Phys. 52, 7182 (1981))
 *
 *   H(L) = U + p V + V_0 sum_a (p_aa - p) (L_a / L_a^0 - 1)
 *
 * with p = (p_xx + p_yy + p_zz) / 3 the hydrostatic part of the pressure tensor and the deviatoric part acting
 * as a constant tension on each pair of box faces, measured relative to the reference cell h_0
 * ('System::referenceSimulationBox', volume V_0 and lengths L_a^0). The returned value is the difference of the
 * pressure terms of H between the new and the old box, i.e. an exact difference of a state function, so that the
 * acceptance rule obeys detailed balance. For a hydrostatic pressure tensor it reduces to p (V_new - V_old).
 *
 * The earlier first-order form V_old sum_a p_aa (s_a - 1) is not a state function: the forward and reverse
 * works do not cancel (even for a hydrostatic tensor it differs from p dV at second order) and the resulting
 * chain samples no equilibrium ensemble.
 */
double anisotropicPressureWork(const System& system, const SimulationBox& oldBox, const SimulationBox& newBox);

/**
 * \brief Performs an anisotropic volume move with independent pressures along a, b, and c.
 *
 * Scales each lattice vector length independently (the angles are kept) and accepts with the constant-stress
 * enthalpy of 'anisotropicPressureWork'.
 */
std::optional<RunningEnergy> anisotropicVolumeMove(RandomNumber& random, System& system);
}  // namespace MC_Moves
