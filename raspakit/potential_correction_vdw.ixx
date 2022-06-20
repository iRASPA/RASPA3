export module potential_correction_vdw;

import forcefield;

import double4;
import <cmath>;

export inline double potentialCorrectionVDW(const ForceField& forcefield, const size_t& typeA, const size_t& typeB)
{
	VDWParameters::Type potentialType = forcefield(typeA, typeB).type;

    double cutOffVDW = forcefield.cutOff;

	switch (potentialType)
	{
	case VDWParameters::Type::LennardJones:
	{
		double arg1 = forcefield(typeA, typeB).parameters.x;
		double arg2 = forcefield(typeA, typeB).parameters.y;
        double term1 = (arg2/cutOffVDW) * (arg2/cutOffVDW) * (arg2/cutOffVDW);
        double term2 = term1 * term1 * term1;
        return (4.0/3.0) * arg1 * arg2 * arg2 *arg2 * ((1.0 / 3.0) * term2 - term1);
	}
    default:
      return 0.0;
	}
};
