module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
#include <string>
#include <vector>
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
import molecule;

export namespace IO
{
void ReadLAMMPSDataFile();

/**
 * \brief Writes output for a LAMMPS data file.
 *
 * Takes system information and writes it in a LAMMPS data file format, such that a simulation can be restarted on a
 * LAMMPS engine. Assumes LAMMPS units real unit system (Angstrom, kcal/mol).
 *
 * \note Currently does not support writing of bonded information.
 *
 * \param components system component information
 * \param atomData holds all information on all atoms in the system.
 * \param simulationBox system simulation box (not unit cell).
 * \param forceField contains parameters for the pair interactions.
 * \param numberOfIntegerMoleculesPerComponent amount of molecules per component, necessary for accounting.
 * \param framework meta info on the framework.
 */
std::string WriteLAMMPSDataFile(std::span<const Component> components, std::span<const Atom> atomData,
                                std::span<const Molecule> moleculeData, const SimulationBox simulationBox,
                                const ForceField forceField,
                                std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                                std::optional<Framework> framework);
}  // namespace IO
