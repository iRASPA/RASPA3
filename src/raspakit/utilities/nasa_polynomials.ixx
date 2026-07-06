module;

export module nasa_polynomials;

import std;

// Computation of ideal-gas molecular partition functions from 7-term NASA polynomials, using
// the "Extended Third Millennium Ideal Gas and Condensed Phase Thermochemical Database for
// Combustion with Updates from Active Thermochemical Tables" by E. Goos, A. Burcat, and
// B. Ruscic (https://respecth.elte.hu/burcat.php). All neutral gas-phase species of the
// database composed of the elements H, C, N, O, F, Si, P, S, Cl, Br, I, B, He, Ne, Ar, Kr,
// and Xe are embedded (about 1900 species), typically valid for 200 K to 6000 K.
//
// The method follows appendix A.2 of:
//   A. Rahbari, "Advanced Monte Carlo Simulations of Chemical Reaction Equilibrium and
//   Thermodynamic Properties of Mixtures", PhD thesis, Delft University of Technology, 2020.
//
// The NASA polynomials give Cp°(T)/R, H°(T)/RT, and S°(T)/R of the ideal gas, with the
// enthalpy referenced such that H°(298.15 K) = ΔfH°(298.15 K) (elements in their reference
// states at 298.15 K have zero enthalpy). The standard chemical potential follows from
// G°(T) = H°(T) - T S°(T), and is related to the molecular partition function by (Eq. A79):
//
//   -G°(T)/T = R ln[ (q(V,T)/V) (k_B T / P°) ]
//
// To reference the partition function of each molecule to its dissociated atoms in their
// electronic ground state (the convention used for reaction-ensemble Monte Carlo, consistent
// with the janaf module), a fixed per-atom constant is added, obtained from the JANAF values
// of ΔfH°(0 K) of the gaseous atoms and H°(298.15 K) - H°(0 K) of the reference elements:
//
//   R ln[ (q/V)(k_B T/P°) ] = -G°(T)/T + (1/T) Σᵢ yᵢ [ ΔfH°ᵢ,atom(0 K) - (H°(298.15 K) - H°(0 K))ᵢ,element / mᵢ ]
//
// in which yᵢ is the number of atoms of element i in the molecule and mᵢ the number of atoms
// in the reference-state molecule of element i (e.g. m = 2 for H2). These constants cancel
// for balanced reactions, so this convention is interchangeable per species with the janaf
// module.

export namespace NASAPolynomials
{
/**
 * \brief Checks whether NASA-polynomial data is available for a species.
 *
 * \param speciesName Name as listed in the Burcat database (e.g. "C3H6 propylene",
 *                    "C4H8,tr2-butene") or a common name/formula (e.g. "propene",
 *                    "trans-2-butene", "CO2"). Matching is case-insensitive and ignores
 *                    all non-alphanumeric characters.
 *
 * \return True if the species is present in the embedded database.
 */
bool contains(std::string_view speciesName);

/**
 * \brief Returns the names of all species in the embedded database.
 */
std::vector<std::string_view> availableSpecies();

/**
 * \brief Natural logarithm of the ideal-gas molecular partition function per unit volume,
 *        ln(q(V,T)/V), in units of [Å⁻³].
 *
 * The zero of energy is the dissociated atoms in their electronic ground state, consistent
 * with JANAF::logPartitionFunction(). This is the quantity required for the reaction-ensemble
 * acceptance rules ('LnPartitionFunction' of a component).
 *
 * \param speciesName Name of the species (see contains()).
 * \param temperature Temperature in Kelvin; must lie within the valid range of the polynomial.
 *
 * \return ln(q(V,T)/V) with q/V in [Å⁻³].
 *
 * \throws std::runtime_error If the species is unknown or the temperature is outside the
 *                            valid range.
 */
double logPartitionFunction(std::string_view speciesName, double temperature);

/**
 * \brief Standard molar Gibbs energy G°(T) = H°(T) - T S°(T) in [kJ/mol] at P° = 1 bar.
 *
 * The enthalpy reference is the NASA convention: elements in their reference states at
 * 298.15 K have zero enthalpy, so that H°(298.15 K) = ΔfH°(298.15 K).
 *
 * \param speciesName Name of the species (see contains()).
 * \param temperature Temperature in Kelvin; must lie within the valid range of the polynomial.
 *
 * \return G°(T) in [kJ/mol].
 *
 * \throws std::runtime_error If the species is unknown or the temperature is outside the
 *                            valid range.
 */
double standardGibbsEnergy(std::string_view speciesName, double temperature);
}  // namespace NASAPolynomials
