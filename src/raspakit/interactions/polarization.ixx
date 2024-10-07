module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module interactions_polarization;

#ifndef USE_LEGACY_HEADERS
import <span>;
import <optional>;
import <tuple>;
import <complex>;
import <vector>;
#endif

import double3;
import double3x3;
import atom;
import running_energy;
import energy_status;
import simulationbox;
import force_factor;
import forcefield;
import framework;
import component;

export namespace Interactions
{
RunningEnergy computePolarizationEnergyDifference(const ForceField &forceField, std::span<double3> electricField,
                                                  std::span<double3> electricFieldNew,
                                                  std::span<Atom> moleculeAtomPositions);
}
