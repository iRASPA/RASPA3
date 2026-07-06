module;

export module janaf;

import std;

// Computation of ideal-gas molecular partition functions from the NIST-JANAF Thermochemical
// Tables (M.W. Chase, "NIST-JANAF Thermochemical Tables", Fourth Edition, J. Phys. Chem. Ref.
// Data, Monograph 9, 1998; https://janaf.nist.gov).
//
// The method follows appendix A.2 of:
//   A. Rahbari, "Advanced Monte Carlo Simulations of Chemical Reaction Equilibrium and
//   Thermodynamic Properties of Mixtures", PhD thesis, Delft University of Technology, 2020.
//
// The JANAF tables list the Gibbs energy function -[G°(T) - H°(298.15 K)]/T and the enthalpy
// increment H°(0 K) - H°(298.15 K). Shifting the reference enthalpy to T = 0 K (Eq. A82):
//
//   -[G°(T) - H°(0 K)]/T = -[G°(T) - H°(298.15 K)]/T + [H°(0 K) - H°(298.15 K)]/T
//                        = R ln[ (q₀(V,T)/V) (k_B T / P°) ]
//
// gives the temperature-dependent part of the molecular partition function q₀(V,T)/V, in which
// the ground state (vibrational and electronic) of the molecule is taken as the zero of energy.
// For use in reaction-ensemble Monte Carlo, the zero of energy of all components must be
// consistent; this is achieved by referencing each molecule to its dissociated atoms in their
// electronic ground state (Eq. A74):
//
//   q(V,T)/V = (q₀(V,T)/V) exp[D₀ / (k_B T)]
//
// with the atomization energy D₀ obtained from the JANAF enthalpies of formation at 0 K
// (Eq. A78):
//
//   D₀ = Σᵢ yᵢ ΔfH°ᵢ(0 K) - ΔfH°(0 K)
//
// in which yᵢ is the number of atoms of kind i in the molecule.

export namespace JANAF
{
/**
 * \brief Checks whether JANAF data is available for a molecule.
 *
 * \param moleculeName Name of the molecule; the chemical formula (e.g. "N2", "CO2") or a
 *                     common name (e.g. "nitrogen", "carbon dioxide"). Matching is
 *                     case-insensitive and ignores spaces, hyphens, and underscores.
 *
 * \return True if the molecule is present in the embedded JANAF data.
 */
bool contains(std::string_view moleculeName);

/**
 * \brief Returns the list of canonical names of all molecules with embedded JANAF data.
 */
std::vector<std::string_view> availableMolecules();

/**
 * \brief Natural logarithm of the ideal-gas molecular partition function per unit volume,
 *        ln(q(V,T)/V), in units of [Å⁻³].
 *
 * The zero of energy is the dissociated atoms in their electronic ground state, i.e. the
 * factor exp[D₀/(k_B T)] is included (thesis Eqs. A74, A78, and A82). This is the quantity
 * required for the reaction-ensemble acceptance rules ('LnPartitionFunction' of a component).
 *
 * \param moleculeName Name of the molecule (see contains()).
 * \param temperature  Temperature in Kelvin; must lie within the tabulated range.
 *
 * \return ln(q(V,T)/V) with q/V in [Å⁻³].
 *
 * \throws std::runtime_error If the molecule is unknown or the temperature is outside the
 *                            tabulated range.
 */
double logPartitionFunction(std::string_view moleculeName, double temperature);

/**
 * \brief Natural logarithm of the temperature-dependent part of the partition function,
 *        ln(q₀(V,T)/V), in units of [Å⁻³].
 *
 * The zero of energy is the ground state (vibrational and electronic) of the molecule itself
 * (thesis Eqs. A76 and A82); the factor exp[D₀/(k_B T)] is NOT included. Note that these
 * values are only mutually consistent between different molecules when combined with the
 * ground-state energies (e.g. via atomizationEnergy()).
 *
 * \param moleculeName Name of the molecule (see contains()).
 * \param temperature  Temperature in Kelvin; must lie within the tabulated range.
 *
 * \return ln(q₀(V,T)/V) with q₀/V in [Å⁻³].
 *
 * \throws std::runtime_error If the molecule is unknown or the temperature is outside the
 *                            tabulated range.
 */
double logPartitionFunctionGroundStateReference(std::string_view moleculeName, double temperature);

/**
 * \brief Atomization energy D₀ of the molecule at 0 K in [kJ/mol].
 *
 * Computed from the JANAF enthalpies of formation at 0 K of the molecule and its constituent
 * atoms (thesis Eq. A78): D₀ = Σᵢ yᵢ ΔfH°ᵢ(0 K) - ΔfH°(0 K).
 *
 * \param moleculeName Name of the molecule (see contains()).
 *
 * \return D₀ in [kJ/mol].
 *
 * \throws std::runtime_error If the molecule is unknown.
 */
double atomizationEnergy(std::string_view moleculeName);
}  // namespace JANAF
