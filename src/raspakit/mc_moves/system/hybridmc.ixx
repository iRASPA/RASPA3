module;

export module mc_moves_hybridmc;

import std;

import running_energy;
import randomnumbers;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a hybrid Monte Carlo move on the system.
 *
 * Short NVE Velocity-Verlet trajectory used as a collective Monte Carlo proposal.
 * Supports rigid and flexible molecules in a rigid framework, a flexible framework,
 * or no framework. Flexible-framework atom positions and velocities are included in
 * the trial state and restored only on acceptance.
 *
 * Without a framework (and the other restrictions of ForceEngine::supports) the trajectory is integrated with the
 * spatial-decomposition force engine of the MD driver ('System::hybridMCUseForceEngine', input
 * 'HybridMCUseForceEngine'): the Verlet-list pair kernels on the worker threads (or a device), particle-mesh
 * Ewald, and the molecule-range integrator passes. The acceptance uses the exact energies: velocity Verlet is
 * reversible and volume-preserving for any position-dependent force, so the surrogate forces make a valid
 * proposal and the Metropolis test on the exact Hamiltonian keeps the sampled distribution exact. The returned
 * energies are the exact energies of the new configuration (System::computeTotalEnergies).
 *
 * The same holds for the exact integrator path: its gradient routines carry no external-field or polarization
 * forces, so in such systems the trajectory is a surrogate proposal, and the acceptance is on the exact energies
 * (U_old and U_new from System::computeTotalEnergies, external field and polarization included). The centre-of-mass
 * velocity is removed only when the total momentum is conserved along the trajectory (no rigid framework, no
 * external field).
 *
 * \param random A reference to the random number generator.
 * \param system A reference to the system on which the move is performed.
 * \return An optional RunningEnergy object containing the system's energy if the move is accepted, or std::nullopt if
 * rejected.
 */
std::optional<RunningEnergy> hybridMCMove(RandomNumber& random, System& system);

/// Status and timings of the force engine used by the hybrid MC move of the system (empty when none is used).
std::string hybridMCForceEngineReport(const System& system);
}  // namespace MC_Moves
