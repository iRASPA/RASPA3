module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <optional>
#include <span>
#endif

export module integrators;

#ifndef USE_LEGACY_HEADERS
import <span>;
import <optional>;
import <complex>;
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

export namespace Integrators
{
RunningEnergy velocityVerlet(
    std::span<Molecule> moleculePositions, std::span<Atom> moleculeAtomPositions,
    const std::vector<Component> components, double dt, std::optional<Thermostat>& thermostat,
    std::span<Atom> frameworkAtomPositions, const ForceField& forceField, const SimulationBox& simulationBox,
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>>& fixedFrameworkStoredEik,
    const std::vector<size_t> numberOfMoleculesPerComponent);
}  // namespace Integrators