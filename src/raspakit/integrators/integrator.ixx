module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>
#endif

export module integrators;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import molecule;
import atom;
import component;
import running_energy;
import thermostat;
import integrators_compute;
import integrators_update;
import simulationbox;
import forcefield;
import interpolation_energy_grid;

// integrators.ixx

export namespace Integrators
{
/**
 * \brief Performs the Velocity Verlet integration step for molecular dynamics simulation.
 *
 * Advances the positions and velocities of molecules and atoms using the Velocity Verlet algorithm.
 * It handles both translational and rotational motion, applies thermostat if present, computes energies,
 * and updates gradients.
 *
 * \param moleculePositions The span of molecule positions and orientations.
 * \param moleculeAtomPositions The span of atom positions within molecules.
 * \param components The list of component types in the simulation.
 * \param dt The time step size.
 * \param thermostat Optional thermostat for temperature control.
 * \param frameworkAtomPositions The positions of framework atoms.
 * \param forceField The force field parameters used for computing interactions.
 * \param simulationBox The simulation box defining periodic boundaries.
 * \param eik_x Preallocated complex exponentials for Ewald summation in x-direction.
 * \param eik_y Preallocated complex exponentials for Ewald summation in y-direction.
 * \param eik_z Preallocated complex exponentials for Ewald summation in z-direction.
 * \param eik_xy Preallocated complex exponentials for Ewald summation in xy-plane.
 * \param fixedFrameworkStoredEik Precomputed Ewald sums for the fixed framework.
 * \param numberOfMoleculesPerComponent The number of molecules for each component type.
 *
 * \return The updated running energies after the integration step.
 */
RunningEnergy velocityVerlet(
    std::span<Molecule> moleculePositions, std::span<Atom> moleculeAtomPositions,
    const std::vector<Component> components, double dt, std::optional<Thermostat>& thermostat,
    std::span<Atom> frameworkAtomPositions, const ForceField& forceField, const SimulationBox& simulationBox,
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>>& totalEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>>& fixedFrameworkStoredEik,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::vector<std::size_t> numberOfMoleculesPerComponent);
}  // namespace Integrators
