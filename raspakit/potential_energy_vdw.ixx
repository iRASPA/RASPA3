export module potential_energy_vdw;

import forcefield;

import double4;
import <cmath>;
import <iostream>;

export struct EnergyFactor
{
    double energy;
    double dUdlambda;

    EnergyFactor(double energy, double dUdlambda):
        energy(energy),
        dUdlambda(dUdlambda) {}
};

export inline EnergyFactor potentialVDWEnergy(const ForceField& forcefield, const double& scaling, const double& rr, const size_t& typeA, const size_t& typeB)
{
	VDWParameters::Type potentialType = forcefield(typeA, typeB).type;

	switch (potentialType)
	{
	case VDWParameters::Type::LennardJones:
	{
		double arg1 = 4.0 * forcefield(typeA, typeB).parameters.x;
		double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
		double arg3 = forcefield(typeA, typeB).shift;
		double temp = (rr / arg2);
		double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
        double rri6 = rri3 * rri3;
        double term = arg1 * (rri3 * (rri3 - 1.0)) - arg3;
		return EnergyFactor(scaling * term, scaling < 1.0 ? term + scaling * (1.0 - scaling) * arg1 * (rri6 * (2.0 * rri3 - 1.0)) : 0.0);
        //double alpha = 0.5;
        //double b = 1.0;
        //double c = 6.0;
        //double Y = 0.5*pow(1.0-scaling, b) + pow((sqrt(rr)/forcefield(typeA, typeB).parameters.y), c);
        //return EnergyFactor(arg1 * scaling * (pow(1.0/Y, 12.0/c) - pow(1.0 / Y,6.0/c)), 
        //    arg1 * scaling * pow(1.0 / Y, 6.0/c) * (pow(1.0 / Y, 6.0/c)-1.0+((6.0* scaling * b * alpha)/(c*Y)) * pow(1.0-scaling, b - 1.0)*(2.0*pow(1.0/Y, 6.0/c)-1.0)));
	}
	case VDWParameters::Type::BuckingHam:
	{
		double arg1 = forcefield(typeA, typeB).parameters.x;
		double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
		double arg3 = forcefield(typeA, typeB).shift;
		double temp = (rr / arg2);
		double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
		return EnergyFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3),0.0);
	}
	case VDWParameters::Type::Morse:
	{
		double arg1 = forcefield(typeA, typeB).parameters.x;
		double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
		double arg3 = forcefield(typeA, typeB).shift;
		double temp = (rr / arg2);
		double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
		return EnergyFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0);
	}
	case VDWParameters::Type::FeynmannHibbs:
	{
		double arg1 = forcefield(typeA, typeB).parameters.x;
		double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
		double arg3 = forcefield(typeA, typeB).shift;
		double temp = (rr / arg2);
		double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
		return EnergyFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0);
	}
	case VDWParameters::Type::MM3:
	{
		double arg1 = forcefield(typeA, typeB).parameters.x;
		double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
		double arg3 = forcefield(typeA, typeB).shift;
		double temp = (rr / arg2);
		double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
		return EnergyFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0);
	}
	case VDWParameters::Type::BornHugginsMeyer:
	{
		double arg1 = forcefield(typeA, typeB).parameters.x;
		double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
		double arg3 = forcefield(typeA, typeB).shift;
		double temp = (rr / arg2);
		double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
		return EnergyFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0);
	}
	}
};
