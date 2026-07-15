module;

export module interactions_polarization;

import std;

import double3;
import double3x3;
import atom;
import running_energy;
import energy_status;
import simulationbox;
import forcefield;
import framework;
import component;

export namespace Interactions
{
RunningEnergy computePolarizationEnergyDifference(const ForceField &forceField, std::span<double3> electricField,
                                                  std::span<double3> electricFieldNew,
                                                  std::span<Atom> moleculeAtomPositions);

RunningEnergy computePolarizationEnergyDifference(const ForceField &forceField, std::span<double3> electricField,
                                                  std::span<double3> electricFieldNew,
                                                  std::span<Atom> moleculeAtomPositionsNew,
                                                  std::span<Atom> moleculeAtomPositionsOld);

/**
 * \brief Polarization energy change of the "neighbor" atoms whose field changed because another molecule moved.
 *
 * Given the stored (pre-move) electric field \p electricField and the field change \p electricFieldDelta of every
 * atom in the system, returns the change in polarization energy summed over all atoms:
 * \f$-\tfrac12\sum_i \alpha_i\,(|\mathbf{E}_i+\Delta\mathbf{E}_i|^2 - |\mathbf{E}_i|^2)\f$. Atoms whose field is
 * unchanged (\f$\Delta\mathbf{E}_i = 0\f$) contribute nothing, so the moved molecule's own atoms may be included
 * harmlessly. This makes molecule-molecule polarization an O(N) incremental update rather than an O(N^2) rebuild.
 */
RunningEnergy computePolarizationEnergyNeighborDifference(const ForceField &forceField,
                                                          std::span<const double3> electricField,
                                                          std::span<const double3> electricFieldDelta,
                                                          std::span<const Atom> moleculeAtoms);
}  // namespace Interactions
