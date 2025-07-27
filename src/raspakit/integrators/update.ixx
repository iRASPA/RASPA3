module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>
#endif

export module integrators_update;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import molecule;
import atom;
import component;
import running_energy;
import simulationbox;
import forcefield;
import randomnumbers;
import interpolation_energy_grid;

export namespace Integrators
{
/**
 * \brief Scales the velocities and orientation momenta of molecules.
 *
 * \param moleculePositions Span of molecules whose velocities are to be scaled.
 * \param scaling Pair of scaling factors for velocity and orientation momentum.
 */
void scaleVelocities(std::span<Molecule> moleculePositions, std::pair<double, double> scaling);

/**
 * \brief Removes total system velocity.
 *
 * \param moleculePositions Span of molecules whose velocities are to be scaled.
 */
void removeCenterOfMassVelocityDrift(std::span<Molecule> moleculePositions);
/**
 * \brief Updates the positions of molecules based on their velocities.
 *
 * \param moleculePositions Span of molecules whose positions are to be updated.
 * \param dt Time step for the update.
 */
void updatePositions(std::span<Molecule> moleculePositions, double dt);

/**
 * \brief Updates the velocities and orientation momenta of molecules based on gradients.
 *
 * \param moleculePositions Span of molecules whose velocities are to be updated.
 * \param dt Time step for the update.
 */
void updateVelocities(std::span<Molecule> moleculePositions, double dt);

/**
 * \brief Initializes the velocities according to the Boltzmann distribution
 *
 * \param random A random number generator
 * \param moleculePositions Span of molecules to initialize velocity for.
 * \param components Components to get the inertiaVector
 * \param temperature Temperature to set velocities to.
 */
void initializeVelocities(RandomNumber& random, std::span<Molecule> moleculePositions,
                          const std::vector<Component> components, double temperature);

/**
 * \brief Converts molecule positions and orientations into Cartesian atom positions.
 *
 * \param moleculePositions Span of molecules.
 * \param moleculeAtomPositions Span to store the resulting atom positions.
 * \param components Vector of component definitions.
 */
void createCartesianPositions(std::span<const Molecule> moleculePositions, std::span<Atom> moleculeAtomPositions,
                              std::vector<Component> components);

/**
 * \brief Performs second-order NoSquish free rotor integration.
 *
 * \param moleculePositions Span of molecules to be updated.
 * \param components Vector of component definitions.
 * \param dt Time step for the integration.
 */
void noSquishFreeRotorOrderTwo(std::span<Molecule> moleculePositions, const std::vector<Component> components,
                               double dt);

/**
 * \brief Updates center of mass velocities and orientation momenta of molecules.
 *
 * \param moleculePositions Span of molecules to be updated.
 * \param moleculeAtomPositions Span of atom positions corresponding to the molecules.
 * \param components Vector of component definitions.
 */
void updateCenterOfMassAndQuaternionVelocities(std::span<Molecule> moleculePositions,
                                               std::span<Atom> moleculeAtomPositions,
                                               std::vector<Component> components);

/**
 * \brief Updates gradients of center of mass and orientation for molecules.
 *
 * \param moleculePositions Span of molecules to be updated.
 * \param moleculeAtomPositions Span of atom positions corresponding to the molecules.
 * \param components Vector of component definitions.
 */
void updateCenterOfMassAndQuaternionGradients(std::span<Molecule> moleculePositions,
                                              std::span<Atom> moleculeAtomPositions, std::vector<Component> components);

/**
 * \brief Updates gradients and computes energies due to interactions.
 *
 * \param moleculeAtomPositions Span of molecule atom positions.
 * \param frameworkAtomPositions Span of framework atom positions.
 * \param forceField Force field parameters.
 * \param simulationBox Simulation box parameters.
 * \param components Vector of component definitions.
 * \param eik_x Vector of complex exponentials in x-direction.
 * \param eik_y Vector of complex exponentials in y-direction.
 * \param eik_z Vector of complex exponentials in z-direction.
 * \param eik_xy Vector of complex exponentials in xy-plane.
 * \param fixedFrameworkStoredEik Stored complex exponentials for the fixed framework.
 * \param numberOfMoleculesPerComponent Vector of molecule counts per component.
 * \return The total running energy computed from interactions.
 */
RunningEnergy updateGradients(
    std::span<Atom> moleculeAtomPositions, std::span<Atom> frameworkAtomPositions, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<Component> components,
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>>& totalEik,
    const std::vector<std::pair<std::complex<double>, std::complex<double>>>& fixedFrameworkStoredEik,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::vector<std::size_t> numberOfMoleculesPerComponent);
}  // namespace Integrators
