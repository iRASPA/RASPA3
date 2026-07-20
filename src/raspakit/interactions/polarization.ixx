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

/**
 * \brief Polarization contribution to dU/dlambda of a dUdlambda group (thermodynamic integration).
 *
 * The polarization energy is \f$U = -\tfrac12 \sum_i s_i \alpha_i |\mathbf{E}_i|^2\f$ with \f$s_i\f$ the
 * atom's Coulomb scaling. The lambda of group \p groupId enters twice: through the coupling \f$s_i\f$ of the
 * group's own atoms, and through the field the group's atoms produce at atoms of other molecules (which is
 * proportional to the source's scaling). Returns
 * \f$-\tfrac12 \sum_{i \in g} \alpha_i |\mathbf{E}_i|^2 - \sum_j s_j \alpha_j (\mathbf{E}_j \cdot \mathbf{f}_j)\f$
 * where \f$\mathbf{f}_j\f$ is the unscaled real-space field of the group's atoms at \p j (the reciprocal part
 * of the stored field comes from the fixed framework only and is lambda-independent). The caller multiplies by
 * \f$d s_\mathrm{coul}/d\lambda\f$. Evaluated on demand from the stored field; it is not maintained
 * incrementally in RunningEnergy like the VDW/charge/Ewald parts.
 */
double computePolarizationDUdlambda(const ForceField &forceField, const SimulationBox &simulationBox,
                                    std::span<const Atom> moleculeAtoms, std::span<const double3> electricField,
                                    std::size_t groupId);
}  // namespace Interactions
