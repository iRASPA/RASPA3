module;

export module partition_function;

import std;

// Unified lookup of ideal-gas molecular partition functions ln(q(V,T)/V) for the
// reaction-ensemble acceptance rules ('LnPartitionFunction' of a component).
//
// The Burcat/NASA-polynomial database (module nasa_polynomials, ~1900 gas-phase species,
// 200-6000 K) is the primary source. The NIST-JANAF tables (module janaf, small molecules,
// 100-3000 K) are used as a fallback when a species is missing from the primary source or the
// temperature lies outside the range of its polynomials.
//
// Both sources use the same zero of energy (the dissociated atoms in their electronic ground
// state), so values from the two sources are mutually consistent up to the small differences
// between the underlying databases (see the unit test nasa_polynomials.consistent_with_janaf_tables).

export namespace PartitionFunction
{
/**
 * \brief Checks whether partition-function data is available for a species in any source.
 *
 * \param speciesName Name of the species; the chemical formula (e.g. "CO2"), a common name
 *                    (e.g. "carbon dioxide", "propene"), or a Burcat database name. Matching
 *                    is case-insensitive and ignores spaces, hyphens, and underscores.
 *
 * \return True if the species is present in the NASA-polynomial or JANAF data.
 */
bool contains(std::string_view speciesName);

/**
 * \brief Natural logarithm of the ideal-gas molecular partition function per unit volume,
 *        ln(q(V,T)/V), in units of [Å⁻³].
 *
 * The NASA-polynomial (Burcat) database is the primary source; the JANAF tables are the
 * fallback. The zero of energy is the dissociated atoms in their electronic ground state.
 *
 * \param speciesName Name of the species (see contains()).
 * \param temperature Temperature in Kelvin.
 *
 * \return ln(q(V,T)/V) with q/V in [Å⁻³].
 *
 * \throws std::runtime_error If the species is unknown to all sources, or the temperature is
 *                            outside the valid range of every source containing the species.
 */
double logPartitionFunction(std::string_view speciesName, double temperature);
}  // namespace PartitionFunction
