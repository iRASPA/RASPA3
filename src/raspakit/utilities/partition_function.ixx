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
// Within a single simulation system, one database must be used for all components; do not mix
// JANAF and NASA partition-function lookups in the same reaction ensemble.

export namespace PartitionFunction
{
/**
 * \brief Thermochemical database used for automatic 'LnPartitionFunction' lookup.
 */
enum class Source
{
  NASA,   ///< NASA polynomials (Burcat) as primary source, JANAF tables as fallback.
  JANAF,  ///< NIST-JANAF tables only.
};

/**
 * \brief Checks whether partition-function data is available for a species.
 *
 * \param speciesName Name of the species; the chemical formula (e.g. "CO2"), a common name
 *                    (e.g. "carbon dioxide", "propene"), or a Burcat database name. Matching
 *                    is case-insensitive and ignores spaces, hyphens, and underscores.
 * \param source Database to query. For \c Source::NASA, either database counts as available.
 *
 * \return True if the species is present in the selected database (or its fallback).
 */
bool contains(std::string_view speciesName, Source source = Source::NASA);

/**
 * \brief Natural logarithm of the ideal-gas molecular partition function per unit volume,
 *        ln(q(V,T)/V), in units of [Å⁻³].
 *
 * For \c Source::NASA, the Burcat/NASA-polynomial database is the primary source and the
 * JANAF tables are the fallback. For \c Source::JANAF, only the JANAF tables are used.
 * The zero of energy is the dissociated atoms in their electronic ground state.
 *
 * \param speciesName Name of the species (see contains()).
 * \param temperature Temperature in Kelvin.
 * \param source Database to use.
 *
 * \return ln(q(V,T)/V) with q/V in [Å⁻³].
 *
 * \throws std::runtime_error If the species is unknown to the selected database, or the
 *                            temperature is outside its valid range.
 */
double logPartitionFunction(std::string_view speciesName, double temperature, Source source = Source::NASA);
}  // namespace PartitionFunction
