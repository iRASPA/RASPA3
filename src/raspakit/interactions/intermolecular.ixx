module;

export module interactions_intermolecular;

import std;

import double3;
import double3x3;
import atom;
import atom_dynamics;
import running_energy;
import energy_status;
import simulationbox;
import forcefield;
import component;
import interactions_pair_kernel;

export namespace Interactions
{
/**
 * \brief Computes the inter-molecular energy between atoms.
 *
 * Calculates the van der Waals and Coulombic energy contributions between all pairs of atoms
 * in \p moleculeAtoms, excluding interactions within the same molecule.
 *
 * \param forceField The force field parameters used for the energy calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms for which to compute inter-molecular energies.
 * \return The total inter-molecular energy contributions.
 */
RunningEnergy computeInterMolecularEnergy(const ForceField& forceField, const SimulationBox& simulationBox,
                                          std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the tail correction for inter-molecular van der Waals energy.
 *
 * Calculates the long-range correction to the van der Waals energy between all pairs of atoms
 * in \p moleculeAtoms, excluding interactions within the same molecule.
 *
 * \param forceField The force field parameters used for the energy calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms for which to compute the tail energy correction.
 * \return The tail energy contributions to the total energy.
 */
RunningEnergy computeInterMolecularTailEnergy(const ForceField& forceField, const SimulationBox& simulationBox,
                                              std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Reference (per-atom, O(N^2)) inter-molecular van der Waals tail correction.
 *
 * Retained as the oracle for validating the Brick-CFCMC-style aggregated accounting; not used on the hot path.
 */
RunningEnergy computeInterMolecularTailEnergyReference(const ForceField& forceField, const SimulationBox& simulationBox,
                                                       std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Brick-CFCMC-style aggregated inter-molecular van der Waals tail correction.
 *
 * Computes the tail energy and its per-group dU/dlambda from effective (fractionally-weighted) pseudo-atom-type
 * counts, U_tail = (2 pi / V) sum_a sum_b eff_a eff_b C_ab, in O(nType^2). Mathematically identical to the
 * per-atom sum because the tail constant C_ab depends only on the atom types.
 *
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box (provides the volume).
 * \param effectiveTypeCounts Per-type sum of scalingVDW over all molecule atoms.
 * \param groupCounts Per dU/dlambda group, the number of molecule atoms of each type in that group.
 * \return The tail energy and per-group dU/dlambda contributions.
 */
[[nodiscard]] RunningEnergy computeInterMolecularTailEnergyAggregated(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const double> effectiveTypeCounts,
    const std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& groupCounts) noexcept;

/**
 * \brief Aggregated inter-molecular tail-correction difference for replacing \p oldatoms with \p newatoms.
 *
 * Given the effective type counts of the system in its current ("old") configuration, returns
 * Aggregated(counts - oldatoms + newatoms) - Aggregated(counts). Correct-by-construction parity (energy and
 * dU/dlambda) with computeInterMolecularTailEnergyDifference.
 */
[[nodiscard]] RunningEnergy computeInterMolecularTailEnergyDifferenceAggregated(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const double> effectiveTypeCounts,
    const std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& groupCounts, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms) noexcept;

/**
 * \brief Advances local effective type counts by (+newatoms, -oldatoms).
 *
 * Used by CFCMC moves to thread the intermediate effective counts across sequential sub-steps.
 */
void updateEffectiveTypeCounts(std::vector<double>& effectiveTypeCounts,
                               std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& groupCounts,
                               std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept;

/**
 * \brief Computes the difference in inter-molecular energy due to atom changes.
 *
 * Calculates the energy difference resulting from replacing \p oldatoms with \p newatoms
 * in the system, considering interactions with \p moleculeAtoms. Excludes interactions within
 * the same molecule. Returns std::nullopt if an overlap is detected based on the force field's overlap criteria.
 *
 * \param forceField The force field parameters used for the energy calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of existing atoms in the system.
 * \param newatoms A span of new atoms to be added to the system.
 * \param oldatoms A span of atoms to be removed from the system.
 * \return The energy difference due to the atom changes, or std::nullopt if an overlap occurs.
 */
[[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergyDifference(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> moleculeAtoms,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept;

/**
 * \brief Computes the difference in inter-molecular tail energy due to atom changes.
 *
 * Calculates the change in tail correction to the van der Waals energy resulting from replacing
 * \p oldatoms with \p newatoms in the system, considering interactions with \p moleculeAtoms.
 * Excludes interactions within the same molecule.
 *
 * \param forceField The force field parameters used for the energy calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of existing atoms in the system.
 * \param newatoms A span of new atoms to be added to the system.
 * \param oldatoms A span of atoms to be removed from the system.
 * \return The change in tail energy contributions due to the atom changes.
 */
[[nodiscard]] RunningEnergy computeInterMolecularTailEnergyDifference(const ForceField& forceField,
                                                                      const SimulationBox& simulationBox,
                                                                      std::span<const Atom> moleculeAtoms,
                                                                      std::span<const Atom> newatoms,
                                                                      std::span<const Atom> oldatoms) noexcept;

/**
 * \brief Computes inter-molecular tail energy from pseudo-atom type counts.
 */
[[nodiscard]] RunningEnergy computeInterMolecularTailEnergyFromTypeCounts(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const std::size_t> typeCounts) noexcept;

/**
 * \brief Tail correction difference when replacing one molecule type with another (identity change).
 *
 * Equivalent to RASPA2 TailMolecularEnergyDifferenceAddRemove.
 */
[[nodiscard]] RunningEnergy computeInterMolecularTailEnergyDifferenceAddRemove(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const std::size_t> currentTypeCounts,
    const Component& componentToAdd, const Component& componentToRemove) noexcept;

/**
 * \brief Tail correction difference for a reaction ensemble move (RASPA2 TailMolecularEnergyDifferenceRXMX).
 */
[[nodiscard]] RunningEnergy computeInterMolecularTailEnergyDifferenceReaction(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const std::size_t> currentTypeCounts,
    const std::vector<std::size_t>& reactantStoichiometry, const std::vector<std::size_t>& productStoichiometry,
    const std::vector<Component>& components, bool forward) noexcept;

/**
 * \brief Computes the inter-molecular forces and energy.
 *
 * Calculates the van der Waals and Coulombic forces between all pairs of atoms
 * in \p moleculeAtoms, excluding interactions within the same molecule. Updates the
 * gradient (force) field of each Atom accordingly.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms for which to compute inter-molecular forces and energies.
 * \return The total inter-molecular energy contributions.
 */
RunningEnergy computeInterMolecularGradient(const ForceField& forceField, const SimulationBox& simulationBox,
                                            std::span<const Atom> moleculeAtoms,
                                            std::span<AtomDynamics> moleculeDynamics) noexcept;

/**
 * \brief Computes the inter-molecular gradient (force) acting on a single, selected molecule.
 *
 * Unlike \ref computeInterMolecularGradient, which evaluates the force on every molecule (an O(N^2) loop
 * over all pairs), this variant only accumulates the gradient on the atoms of one selected molecule due to
 * the atoms of all other molecules. The cost is O(n_selected * N), matching the cost of a single-molecule
 * energy-difference evaluation. It is intended for force-biased single-molecule moves.
 *
 * Interactions between the selected molecule and itself are skipped through the molecule-id comparison, so
 * \p selectedAtoms may be a trial configuration of a molecule that is still present (at other positions) in
 * \p moleculeAtoms.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of all molecule atoms currently in the system.
 * \param selectedAtoms A span of the atoms of the selected molecule (current or trial positions).
 * \param selectedDynamics A span (size = selectedAtoms.size()) into which the per-atom gradient is accumulated.
 */
void computeInterMolecularGradientMolecule(const ForceField& forceField, const SimulationBox& simulationBox,
                                           std::span<const Atom> moleculeAtoms, std::span<const Atom> selectedAtoms,
                                           std::span<AtomDynamics> selectedDynamics) noexcept;

/**
 * \brief Computes inter-molecular energy, forces, and strain derivative tensor.
 *
 * Calculates the van der Waals and Coulombic energies and forces between all pairs of atoms
 * in \p moleculeAtoms, excluding interactions within the same molecule. Updates the gradient
 * (force) field of each Atom and computes the strain derivative tensor for stress calculations.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param components The list of components in the system.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms for which to compute energies, forces, and strain derivatives.
 * \param polarizationGather Optional accumulator: when non-null the Coulomb pair loop additionally
 *        gathers the polarization field and its strain response (see PolarizationFieldStrain); the
 *        caller must only pass it when inter-molecular polarization is active.
 * \return A pair containing the energy status and the strain derivative tensor.
 */
std::pair<EnergyStatus, double3x3> computeInterMolecularEnergyStrainDerivative(
    const ForceField& forceField, const std::vector<Component>& components, const SimulationBox& simulationBox,
    std::span<const Atom> moleculeAtoms, std::span<AtomDynamics> moleculeDynamics,
    const PolarizationFieldStrain* polarizationGather = nullptr) noexcept;

/**
 * \brief Computes the inter-molecular electric potential for each atom.
 *
 * Calculates the electric potential at each atom in \p moleculeAtoms due to the charges of
 * other atoms, excluding interactions within the same molecule. The results are stored in
 * \p electricPotentialMolecules.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param box The simulation box containing the atoms.
 * \param electricPotentialMolecules A span to store the computed electric potentials.
 * \param moleculeAtoms A span of atoms for which to compute the electric potentials.
 */
void computeInterMolecularElectrostaticPotential(const ForceField& forceField, const SimulationBox& box,
                                                 std::span<double> electricPotentialMolecules,
                                                 std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the inter-molecular electric field for each atom.
 *
 * Calculates the electric field at each atom in \p moleculeAtoms due to the charges of
 * other atoms, excluding interactions within the same molecule. The results are stored in
 * \p electricFieldMolecules. Also computes inter-molecular energy contributions.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param box The simulation box containing the atoms.
 * \param electricFieldMolecules A span to store the computed electric fields.
 * \param moleculeAtoms A span of atoms for which to compute the electric fields.
 * \return The total inter-molecular energy contributions.
 */
RunningEnergy computeInterMolecularElectricField(const ForceField& forceField, const SimulationBox& box,
                                                 std::span<double3> electricFieldMolecules,
                                                 std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the difference in inter-molecular electric field due to atom changes.
 *
 * Calculates the change in electric field at each atom in \p moleculeAtoms resulting from replacing
 * \p oldatoms with \p newatoms in the system, considering interactions with \p moleculeAtoms.
 * Excludes interactions within the same molecule. Updates \p electricFieldMolecules and \p electricFieldMolecule
 * accordingly. Returns std::nullopt if an overlap is detected based on the force field's overlap criteria.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param box The simulation box containing the atoms.
 * \param electricFieldMolecules A span to store the updated electric fields of existing atoms.
 * \param electricFieldMolecule A span to store the updated electric fields of new atoms.
 * \param moleculeAtoms A span of existing atoms in the system.
 * \param newatoms A span of new atoms to be added to the system.
 * \param oldatoms A span of atoms to be removed from the system.
 * \return The energy difference due to the atom changes, or std::nullopt if an overlap occurs.
 */
std::optional<RunningEnergy> computeInterMolecularElectricFieldDifference(
    const ForceField& forceField, const SimulationBox& box, std::span<double3> electricFieldMolecules,
    std::span<double3> electricFieldMolecule, std::span<const Atom> moleculeAtoms, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms) noexcept;

/**
 * \brief Computes the inter-molecular energy and electric-field difference for a single-molecule move.
 *
 * Variant of \ref computeInterMolecularElectricFieldDifference tailored to polarization in translation/rotation
 * moves. Besides the van der Waals and Coulombic energy difference (with overlap detection), it fills three
 * separate electric-field buffers using bare Coulomb gradients between *different* molecules:
 *   - \p electricFieldMoleculeNew : field on the moved molecule at its trial positions (accumulated, size = moved
 *     molecule),
 *   - \p electricFieldMoleculeOld : field on the moved molecule at its current positions (accumulated, size = moved
 *     molecule),
 *   - \p electricFieldNeighborDelta : the change (new - old) of the field on every other atom in the system caused
 *     by moving the molecule (accumulated, size = all molecule atoms; entries of the moved molecule stay zero).
 *
 * Keeping the moved-molecule new and old contributions separate allows the polarization energy of the moved
 * molecule to be evaluated directly, while the neighbor delta is combined with the stored field to obtain the
 * neighbors' polarization-energy change incrementally.
 *
 * \return The inter-molecular energy difference, or std::nullopt if an overlap is detected.
 */
[[nodiscard]] std::optional<RunningEnergy> computeInterMolecularPolarizationElectricFieldDifference(
    const ForceField& forceField, const SimulationBox& box, std::span<double3> electricFieldNeighborDelta,
    std::span<double3> electricFieldMoleculeNew, std::span<double3> electricFieldMoleculeOld,
    std::span<const Atom> moleculeAtoms, std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept;
};  // namespace Interactions
