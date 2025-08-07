module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <future>
#include <iostream>
#include <numbers>
#include <optional>
#include <span>
#include <thread>
#include <vector>
#endif

module interactions_internal;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import energy_status;
import potential_energy_vdw;
import potential_gradient_vdw;
import potential_energy_coulomb;
import potential_gradient_coulomb;
import potential_correction_vdw;
import potential_electrostatics;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import molecule;
import energy_factor;
import gradient_factor;
import energy_status_inter;
import running_energy;
import internal_potentials;
import bond_potential;
import component;
import units;
import threadpool;


RunningEnergy Interactions::computeIntraMolecularEnergy(const Potentials::InternalPotentials &internalPotentials,
                                                        std::span<const Molecule> moleculePositions,
                                                        std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energy{};

  for(const Molecule &molecule : moleculePositions)
  {
    std::span<const Atom> atom_molecule_span = {&moleculeAtoms[molecule.atomIndex], molecule.numberOfAtoms};
    energy += internalPotentials.computeInternalEnergies(atom_molecule_span);
  }

  return energy;
}

