module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#endif

export module lammps_io;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import component;
import atom;
import simulationbox;
import forcefield;
import framework;

export namespace IO
{
void ReadLAMMPSDataFile();
std::string WriteLAMMPSDataFile(std::span<const Component> components, std::span<const Atom> atomPositions,
                                const SimulationBox simulationBox, const ForceField forceField,
                                std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                                std::optional<Framework> framework);
}  // namespace IO