module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <span>
#include <vector>
#endif

export module integrators_update;

#ifndef USE_LEGACY_HEADERS
import <span>;
import <vector>;
import <complex>;
#endif

import molecule;
import atom;
import component;
import running_energy;
import simulationbox;
import forcefield;

export namespace Integrators
{
void scaleVelocities(std::span<Molecule> moleculePositions, std::pair<double, double> scaling);
void updatePositions(std::span<Molecule> moleculePositions, double dt);
void updateVelocities(std::span<Molecule> moleculePositions, double dt);
void createCartesianPositions(std::span<const Molecule> moleculePositions, std::span<Atom> moleculeAtomPositions,
                              std::vector<Component> components);
void noSquishFreeRotorOrderTwo(std::span<Molecule> moleculePositions, const std::vector<Component> components,
                               double dt);
void noSquishRotate(std::span<Molecule> moleculePositions, const std::vector<Component> components, size_t k,
                    double dt);
void updateCenterOfMassAndQuaternionVelocities(std::span<Molecule> moleculePositions,
                                               std::span<Atom> moleculeAtomPositions,
                                               std::vector<Component> components);
void updateCenterOfMassAndQuaternionGradients(std::span<Molecule> moleculePositions,
                                              std::span<Atom> moleculeAtomPositions, std::vector<Component> components);
RunningEnergy updateGradients(
    std::span<Atom> moleculeAtomPositions, std::span<Atom> frameworkAtomPositions, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<Component> components,
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>>& fixedFrameworkStoredEik,
    const std::vector<size_t> numberOfMoleculesPerComponent);
}  // namespace Integrators