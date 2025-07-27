module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <iostream>
#endif

export module potential_energy_vdw;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import vdwparameters;
import forcefield;
import energy_factor;

import double4;

export namespace Potentials
{
/**
 * \brief Calculates the van der Waals potential energy between two atoms.
 *
 * This function computes the van der Waals (VDW) potential energy based on the specified
 * force field parameters, group identifiers, scaling factors, distance between atoms,
 * and their respective types. It supports multiple VDW potential types including
 * Lennard-Jones, Bucking-Ham, Morse, Feynmann-Hibbs, MM3, and Born-Huggins-Meyer.
 *
 * The scaling is linear: it first activates Lennard-Jones interactions from 0 to 0.5,
 * then activates electrostatic interactions from 0.5 to 1.0.
 *
 * \param forcefield The force field parameters used for the calculation.
 * \param groupIdA Boolean indicating if the first atom is part of a specific group.
 * \param groupIdB Boolean indicating if the second atom is part of a specific group.
 * \param scalingA Scaling factor for the first atom's interactions.
 * \param scalingB Scaling factor for the second atom's interactions.
 * \param rr The distance between the two atoms.
 * \param typeA The type identifier for the first atom.
 * \param typeB The type identifier for the second atom.
 * \return An EnergyFactor object containing the calculated potential energy and lambda derivative.
 */
[[clang::always_inline]] inline EnergyFactor potentialVDWEnergy(const ForceField& forcefield, const bool& groupIdA,
                                                                const bool& groupIdB, const double& scalingA,
                                                                const double& scalingB, const double& rr,
                                                                const std::size_t& typeA, const std::size_t& typeB)
{
  VDWParameters::Type potentialType = forcefield(typeA, typeB).type;

  double scaling = scalingA * scalingB;
  switch (potentialType)
  {
    [[likely]] case VDWParameters::Type::LennardJones:
    {
      double arg1 = 4.0 * forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
      double arg3 = forcefield(typeA, typeB).shift;
      double temp = (rr / arg2);
      double temp3 = temp * temp * temp;
      double inv_scaling = 1.0 - scaling;
      double rri3 = 1.0 / (temp3 + 0.5 * inv_scaling * inv_scaling);
      double rri6 = rri3 * rri3;
      double term = arg1 * (rri3 * (rri3 - 1.0)) - arg3;
      double dlambda_term = arg1 * scaling * inv_scaling * (2.0 * rri6 * rri3 - rri6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::Morse:
    {
      double wellDepth = forcefield(typeA, typeB).parameters.x;
      double stiffness = forcefield(typeA, typeB).parameters.y;
      double equilibriumDistance = forcefield(typeA, typeB).parameters.z;
      double r = std::sqrt(rr);
      double scaledDistance = -stiffness * (r - equilibriumDistance);
      double expTerm = std::exp(scaledDistance);
      double energy = wellDepth * (1 - expTerm) * (1 - expTerm);
      return EnergyFactor(energy, 0.0);
    }
    case VDWParameters::Type::RepulsiveHarmonic:
    {
      double r = std::sqrt(rr);
      double arg1 = forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y;
      double temp = (1.0 - r / arg2);
      return EnergyFactor(r >= arg2 ? 0.0 : 0.5 * arg1 * temp * temp, 0.0);
    }
    default:
      return EnergyFactor(0.0, 0.0);
  }

  return EnergyFactor(0.0, 0.0);
};
}  // namespace Potentials
