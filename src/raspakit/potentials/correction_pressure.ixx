module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#endif

export module potential_correction_pressure;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import vdwparameters;
import double4;

export namespace Potentials
{
/**
 * \brief Calculates the potential correction for van der Waals (VDW) interactions.
 *
 * This function computes the potential correction for VDW interactions between two
 * atom types based on the provided force field parameters. It determines the type
 * of VDW potential and calculates the correction accordingly.
 *
 * \param forcefield The force field containing interaction parameters.
 * \param typeA The type identifier for the first atom.
 * \param typeB The type identifier for the second atom.
 *
 * \return The calculated VDW potential correction.
 */
inline double potentialCorrectionPressure(VDWParameters::Type potentialType, double4 &parameters, double cutOffVDW, const std::size_t& typeA,
                                          const std::size_t& typeB)
{
  switch (potentialType)
  {
    case VDWParameters::Type::LennardJones:
    {
      double arg1 = parameters.x;
      double arg2 = parameters.y;
      double term3 = (arg2 / cutOffVDW) * (arg2 / cutOffVDW) * (arg2 / cutOffVDW);
      double term6 = term3 * term3;
      return (8.0 / 3.0) * arg1 * arg2 * arg2 * arg2 * ((2.0 / 3.0) * term6 * term3 - term3);
    }
    default:
      return 0.0;
  }
};
}  // namespace Potentials
