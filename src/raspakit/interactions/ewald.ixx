module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module interactions_ewald;

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
import gradient_factor;
import forcefield;
import framework;
import component;

export namespace Interactions
{
/**
 * \brief Computes the Ewald Fourier energy contribution for a single ion.
 *
 * Calculates the Fourier-space part of the Ewald summation for a single ion,
 * given its position and charge. This function is useful for computing
 * corrections due to net charges in the system or for inserting ions.
 *
 * References:
 * - An Exact Ewald Summation Method in Theory and Practice, S. Stenberg and B. Stenqvist,
 *   J. Phys. Chem. A 2020, 124, 3943−3946; https://doi.org/10.1021/acs.jpca.0c01684
 * - Removal of pressure and free energy artifacts in charged periodic systems via net charge corrections
 *   to the Ewald potential, Stephen Bogusz, Thomas E. Cheatham III, and Bernard R. Brooks,
 *   J. Chem. Phys. 108, 7070 (1998); https://doi.org/10.1063/1.476320
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param position The position of the ion.
 * \param charge The charge of the ion.
 * \return The Fourier energy contribution for the ion.
 */
double computeEwaldFourierEnergySingleIon(std::vector<std::complex<double>> &eik_x,
                                          std::vector<std::complex<double>> &eik_y,
                                          std::vector<std::complex<double>> &eik_z,
                                          std::vector<std::complex<double>> &eik_xy, const ForceField &forceField,
                                          const SimulationBox &simulationBox, double3 position, double charge);

/**
 * \brief Precomputes the Ewald Fourier terms for a rigid framework.
 *
 * Calculates and stores the Fourier-space components for a rigid framework,
 * which can be reused in subsequent Ewald summation calculations to improve efficiency.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Storage for the precomputed Fourier components.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param frameworkAtoms The atoms of the rigid framework.
 */
void precomputeEwaldFourierRigid(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms);

/**
 * \brief Computes the Ewald Fourier energy for the system.
 *
 * Calculates the Fourier-space part of the Ewald summation for the entire system,
 * including interactions between molecules and between molecules and the rigid framework.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param storedEik Storage for the total Fourier components.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param components The molecular components in the system.
 * \param numberOfMoleculesPerComponent Number of molecules per component.
 * \param moleculeAtoms Positions and properties of the molecules' atoms.
 * \return The running energy containing the Ewald Fourier energy contributions.
 */
RunningEnergy computeEwaldFourierEnergy(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms);

/**
 * \brief Computes the energy difference due to atom position changes in the Ewald Fourier summation.
 *
 * Calculates the change in Fourier-space Ewald energy when atoms are moved from old positions to new positions.
 * Useful for Monte Carlo moves or molecular dynamics steps.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param storedEik Previously stored Fourier components of the system.
 * \param totalEik Updated Fourier components after the move.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param newatoms The new positions and properties of the atoms.
 * \param oldatoms The old positions and properties of the atoms.
 * \return The running energy containing the Ewald Fourier energy difference.
 */
RunningEnergy energyDifferenceEwaldFourier(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &totalEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<const Atom> newatoms, std::span<const Atom> oldatoms);

/**
 * \brief Computes the Ewald Fourier energy and its gradient (forces) on atoms.
 *
 * Calculates the Fourier-space part of the Ewald summation and computes the forces acting on each atom.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param components The molecular components in the system.
 * \param numberOfMoleculesPerComponent Number of molecules per component.
 * \param atomPositions Positions and properties of the atoms (updated with computed gradients).
 * \return The running energy containing the Ewald Fourier energy contributions.
 */
RunningEnergy computeEwaldFourierGradient(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &totalEik,
    const std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    const ForceField &forceField, const SimulationBox &simulationBox, const std::vector<Component> &components,
    const std::vector<size_t> &numberOfMoleculesPerComponent, std::span<Atom> atomPositions);

/**
 * \brief Computes the Ewald Fourier energy and its strain derivative.
 *
 * Calculates the Fourier-space part of the Ewald summation, including the strain derivative, which is necessary
 * for pressure and stress calculations in simulations involving variable cell shapes.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param storedEik Stored Fourier components of the system (unused).
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param frameworkComponents The framework components in the system.
 * \param components The molecular components in the system.
 * \param numberOfMoleculesPerComponent Number of molecules per component.
 * \param atomPositions Positions and properties of the atoms (updated with computed gradients).
 * \param UIon The energy of a single ion in the Ewald summation.
 * \param netChargeFramework Net charge of the framework.
 * \param netChargePerComponent Net charges per component.
 * \return A pair containing the energy status and the strain derivative tensor.
 */
std::pair<EnergyStatus, double3x3> computeEwaldFourierEnergyStrainDerivative(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::optional<Framework> &framework,
    const std::vector<Component> &components, const std::vector<size_t> &numberOfMoleculesPerComponent,
    std::span<Atom> atomPositions, double UIon, double netChargeFramework,
    std::vector<double> netChargePerComponent) noexcept;

/**
 * \brief Accepts a move by updating the stored Ewald Fourier components.
 *
 * Updates the stored Fourier components after a move is accepted, ensuring that future energy calculations
 * are based on the new configuration.
 *
 * \param forceField The force field parameters.
 * \param storedEik The stored Fourier components to be updated.
 * \param totalEik The new Fourier components after the move.
 */
void acceptEwaldMove(const ForceField &forceField,
                     std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik,
                     std::vector<std::pair<std::complex<double>, std::complex<double>>> &totalEik);

/**
 * \brief Computes the electric potential at molecule positions using the Ewald Fourier summation.
 *
 * Calculates the electrostatic potential at each atom's position due to the Fourier-space Ewald summation,
 * including contributions from the rigid framework.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param electricPotentialMolecules Output array to store the computed electric potentials for each molecule.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param components The molecular components in the system.
 * \param numberOfMoleculesPerComponent Number of molecules per component.
 * \param moleculeAtomPositions Positions and properties of the molecules' atoms.
 */
void computeEwaldFourierElectricPotential(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    std::span<double> electricPotentialMolecules, const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<Component> &components, const std::vector<size_t> &numberOfMoleculesPerComponent,
    std::span<const Atom> moleculeAtomPositions);

/**
 * \brief Computes the electric field at molecule positions using the Ewald Fourier summation.
 *
 * Calculates the electrostatic field at each atom's position due to the Fourier-space Ewald summation,
 * including contributions from the rigid framework.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param storedEik Storage for the total Fourier components.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param electricFieldMolecules Output array to store the computed electric fields for each molecule.
 * \param components The molecular components in the system.
 * \param numberOfMoleculesPerComponent Number of molecules per component.
 * \param atomPositions Positions and properties of the atoms (updated with computed fields).
 * \return The running energy containing the Ewald Fourier energy contributions.
 */
RunningEnergy computeEwaldFourierElectricField(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<double3> electricFieldMolecules,
    const std::vector<Component> &components, const std::vector<size_t> &numberOfMoleculesPerComponent,
    std::span<Atom> atomPositions);

/**
 * \brief Computes the difference in electric field due to atom position changes in the Ewald Fourier summation.
 *
 * Calculates the change in electric field when atoms are moved from old positions to new positions,
 * using the Fourier-space Ewald summation.
 *
 * \param eik_x Preallocated vector to temporarily store exponential terms along x-axis.
 * \param eik_y Preallocated vector to temporarily store exponential terms along y-axis.
 * \param eik_z Preallocated vector to temporarily store exponential terms along z-axis.
 * \param eik_xy Preallocated vector to temporarily store exponential terms along xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Fourier components of the rigid framework.
 * \param storedEik Previously stored Fourier components of the system.
 * \param totalEik Updated Fourier components after the move.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box parameters.
 * \param electricField Output array to store the computed electric fields difference.
 * \param newatoms The new positions and properties of the atoms.
 * \param oldatoms The old positions and properties of the atoms.
 * \return The running energy containing the Ewald Fourier energy difference.
 */
RunningEnergy eletricFieldDifferenceEwaldFourier(
    std::vector<std::complex<double>> &eik_x, std::vector<std::complex<double>> &eik_y,
    std::vector<std::complex<double>> &eik_z, std::vector<std::complex<double>> &eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &fixedFrameworkStoredEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>> &totalEik, const ForceField &forceField,
    const SimulationBox &simulationBox, std::span<double3> electricField, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms);
}  // namespace Interactions
