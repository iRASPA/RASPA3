export module potential_gradient_vdw;

import forcefield;
import force_factor;

import double4;
import <cmath>;
import <iostream>;

export inline ForceFactor potentialVDWGradient(const ForceField& forcefield, const double& scaling, const double& rr, const size_t& typeA, const size_t& typeB)
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
        double temp3 = temp * temp * temp;
        double rri3 = 1.0 / (temp3 + 0.5 * (1.0 - scaling) * (1.0 - scaling));
        double rri6 = rri3 * rri3;
        double term = arg1 * (rri3 * (rri3 - 1.0)) - arg3;
        double dlambda_term = scaling * arg1 * (rri6 * (2.0 * rri3 - 1.0));
        return ForceFactor(scaling * term, 
                            //6.0 * temp3 * dlambda_term,
                            12.0 * arg1 * (rri3 * (rri3 - 0.5))/rr,
                            scaling < 1.0 ? term + (1.0 - scaling) * dlambda_term : 0.0);
    }
    case VDWParameters::Type::BuckingHam:
    {
        double arg1 = forcefield(typeA, typeB).parameters.x;
        double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
        double arg3 = forcefield(typeA, typeB).shift;
        double temp = (rr / arg2);
        double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
        return ForceFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0, 0.0);
    }
    case VDWParameters::Type::Morse:
    {
        double arg1 = forcefield(typeA, typeB).parameters.x;
        double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
        double arg3 = forcefield(typeA, typeB).shift;
        double temp = (rr / arg2);
        double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
        return ForceFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0, 0.0);
    }
    case VDWParameters::Type::FeynmannHibbs:
    {
        double arg1 = forcefield(typeA, typeB).parameters.x;
        double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
        double arg3 = forcefield(typeA, typeB).shift;
        double temp = (rr / arg2);
        double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
        return ForceFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0, 0.0);
    }
    case VDWParameters::Type::MM3:
    {
        double arg1 = forcefield(typeA, typeB).parameters.x;
        double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
        double arg3 = forcefield(typeA, typeB).shift;
        double temp = (rr / arg2);
        double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
        return ForceFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0, 0.0);
    }
    case VDWParameters::Type::BornHugginsMeyer:
    {
        double arg1 = forcefield(typeA, typeB).parameters.x;
        double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
        double arg3 = forcefield(typeA, typeB).shift;
        double temp = (rr / arg2);
        double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
        return ForceFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0, 0.0);
    }
    }
};
