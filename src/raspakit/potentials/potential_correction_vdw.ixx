module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#endif

export module potential_correction_vdw;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
#endif

import vdwparameters;
import forcefield;
import double4;

export inline double potentialCorrectionVDW(const ForceField& forcefield, const size_t& typeA, const size_t& typeB)
{
  VDWParameters::Type potentialType = forcefield(typeA, typeB).type;

  double cutOffVDW = forcefield.cutOffVDW;

  switch (potentialType)
  {
    case VDWParameters::Type::LennardJones:
    {
      double arg1 = forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y;
      double term3 = (arg2 / cutOffVDW) * (arg2 / cutOffVDW) * (arg2 / cutOffVDW);
      double term6 = term3 * term3;
      return (4.0 / 3.0) * arg1 * arg2 * arg2 * arg2 * ((1.0 / 3.0) * term6 * term3 - term3);
    }
    default:
      return 0.0;
  }
};
