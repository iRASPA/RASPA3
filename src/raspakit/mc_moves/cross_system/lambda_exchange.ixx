module;

export module mc_moves_lambda_exchange;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Hamiltonian replica-exchange of the fixed lambda values of two neighboring systems.
 *
 * Both systems are replicas of the same physical system, each holding fractional molecule(s)
 * pinned at a fixed lambda (thermodynamic integration at constant lambda). The move proposes to
 * exchange the two lambda values: the fractional molecule(s) of system A are rescaled from
 * lambda_A to lambda_B and those of system B from lambda_B to lambda_A, while all positions stay
 * unchanged. Only the interactions of the fractional molecule(s) change, so the energy differences
 *
 *     dU_A = U_A(lambda_B) - U_A(lambda_A),   dU_B = U_B(lambda_A) - U_B(lambda_B)
 *
 * are cheap single-molecule-type difference computations. The exchange is accepted with
 *
 *     acc = min(1, exp(-beta_A dU_A - beta_B dU_B))
 *
 * which for equal temperatures reduces to the standard Hamiltonian parallel-tempering rule.
 * On acceptance the pinned lambda-bins of the two systems are swapped and both running energies
 * are updated with the (group-resolved, dU/dlambda-exact) energy differences.
 *
 * \param random   Random number generator used for the acceptance test.
 * \param systemA  First replica (owns lambda_A).
 * \param systemB  Neighboring replica (owns lambda_B).
 * \return         The pair (dU_A, dU_B) if the exchange was accepted; std::nullopt otherwise.
 */
std::optional<std::pair<RunningEnergy, RunningEnergy>> LambdaExchange(RandomNumber &random, System &systemA,
                                                                      System &systemB);
}  // namespace MC_Moves
