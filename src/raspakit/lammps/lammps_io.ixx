module;

export module lammps_io;

import std;

import component;
import atom;
import atom_dynamics;
import simulationbox;
import forcefield;
import framework;
import molecule;
import lammps_topology;

/**
 * \brief LAMMPS export: data file, companion input script, and table files.
 *
 * Everything is derived from a LAMMPS::Topology (see lammps_topology). The data file ('units real',
 * 'atom_style full') carries the box, masses, deduplicated coefficient sections, atoms with image flags,
 * velocities, topology lists and, when rigid fragments exist, a custom 'Fragments' section for
 * 'fix property/atom'. The input script carries what a data file cannot: styles (hybrid where needed),
 * pair_style/pair_coeff, kspace, special_bonds, class2 cross coefficients, table references,
 * fix shake / fix rigid/small / frozen framework atoms, and a 'run 0' so the exported energy can be
 * compared with RASPA's.
 */
export namespace LAMMPS
{
struct ExportFiles
{
  std::string data{};
  std::string input{};
  std::string table{};     ///< empty when no potential had to be tabulated
  std::string pairList{};  ///< empty unless 1-4 pairs need 'pair_style list'
  std::vector<std::string> warnings{};
};

std::string writeDataFile(const Topology &topology);
std::string writeInputScript(const Topology &topology);
std::string writeTableFile(const Topology &topology);
std::string writePairListFile(const Topology &topology);

ExportFiles exportSystem(std::span<const Component> components, std::span<const Atom> atomData,
                         std::span<const AtomDynamics> atomDynamics, std::span<const Molecule> moleculeData,
                         const SimulationBox &simulationBox, const ForceField &forceField,
                         std::span<const std::size_t> numberOfIntegerMoleculesPerComponent,
                         const std::optional<Framework> &framework, ExportOptions options = {});
}  // namespace LAMMPS

export namespace IO
{
/**
 * \brief Writes the LAMMPS data file of a system (kept for the WriteLammpsData property).
 *
 * Equivalent to LAMMPS::exportSystem(...).data; the companion input script is written by the property
 * next to it.
 */
std::string WriteLAMMPSDataFile(std::span<const Component> components, std::span<const Atom> atomData,
                                std::span<const AtomDynamics> atomDynamics, std::span<const Molecule> moleculeData,
                                const SimulationBox simulationBox, const ForceField forceField,
                                std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                                std::optional<Framework> framework);
}  // namespace IO
