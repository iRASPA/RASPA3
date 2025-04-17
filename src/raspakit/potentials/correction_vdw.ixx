module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <memory>
#endif

export module potential_correction_vdw;

#ifndef USE_LEGACY_HEADERS
import <memory>;
#endif

import vdwparameters;
import forcefield;
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
inline double potentialCorrectionVDW(const ForceField& forceField, const size_t& typeA, const size_t& typeB)
{
  VDWParameters::Type potentialType = forceField(typeA, typeB).type;

  double cutOffVDW = forceField.cutOffVDW(typeA, typeB);

  switch (potentialType)
  {
    case VDWParameters::Type::LennardJones:
    {
      double arg1 = forceField(typeA, typeB).parameters.x;
      double arg2 = forceField(typeA, typeB).parameters.y;
      double term3 = (arg2 / cutOffVDW) * (arg2 / cutOffVDW) * (arg2 / cutOffVDW);
      double term6 = term3 * term3;
      return (4.0 / 3.0) * arg1 * arg2 * arg2 * arg2 * ((1.0 / 3.0) * term6 * term3 - term3);
    }
    default:
      return 0.0;
  }
};
}  // namespace Potentials
