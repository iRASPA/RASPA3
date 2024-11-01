module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <span>
#include <type_traits>
#include <vector>
#endif

module interactions_polarization;

#ifndef USE_LEGACY_HEADERS
import <complex>;
import <span>;
import <numbers>;
import <cmath>;
import <vector>;
import <iostream>;
import <algorithm>;
import <type_traits>;
#endif

import int3;
import double3;
import double3x3;
import atom;
import simulationbox;
import energy_status;
import energy_status_inter;
import units;
import energy_factor;
import gradient_factor;
import running_energy;
import framework;
import component;
import forcefield;

RunningEnergy Interactions::computePolarizationEnergyDifference(const ForceField &forceField,
                                                                std::span<double3> electricFieldNew,
                                                                std::span<double3> electricField,
                                                                std::span<Atom> moleculeAtomPositions)
{
  RunningEnergy energy;

  for (size_t i = 0; i < moleculeAtomPositions.size(); ++i)
  {
    size_t type = moleculeAtomPositions[i].type;
    double polarizability = forceField.pseudoAtoms[type].polarizability / Units::CoulombicConversionFactor;
    energy.polarization -= 0.5 * polarizability * double3::dot(electricFieldNew[i], electricFieldNew[i]);
  }

  for (size_t i = 0; i < moleculeAtomPositions.size(); ++i)
  {
    size_t type = moleculeAtomPositions[i].type;
    double polarizability = forceField.pseudoAtoms[type].polarizability / Units::CoulombicConversionFactor;
    energy.polarization += 0.5 * polarizability * double3::dot(electricField[i], electricField[i]);
  }

  return energy;
}
