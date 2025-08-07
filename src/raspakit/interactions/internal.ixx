module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module interactions_internal;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import double3x3;
import atom;
import molecule;
import running_energy;
import energy_status;
import simulationbox;
import energy_factor;
import internal_potentials;
import bond_potential;
import gradient_factor;
import forcefield;
import component;

export namespace Interactions
{
  RunningEnergy computeIntraMolecularEnergy(const Potentials::InternalPotentials &internalPotentials,
                                            std::span<const Molecule> moleculePositions,
                                            std::span<const Atom> moleculeAtoms) noexcept;
};  // namespace Interactions
