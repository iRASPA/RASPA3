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

export namespace IO
{
void ReadLAMMPSDataFile();

/**
 * \brief Writes output for a LAMMPS data file.
 *
 * Takes system information and writes it in a LAMMPS data file format ('units real': Angstrom, kcal/mol, fs;
 * 'atom_style full'), such that a simulation can be continued on a LAMMPS engine.
 *
 * Bonded terms are mapped exactly onto stock LAMMPS styles where one exists (RASPA's (1/2)k conventions become
 * LAMMPS's K without the 1/2; every cosine-series torsion is expanded into 'dihedral_style nharmonic'); terms
 * without a LAMMPS equivalent are written with the 'zero' style and flagged in a comment, so the type numbering of
 * the Bonds/Angles/Dihedrals sections always matches the coefficient sections. When a class needs more than one
 * LAMMPS style the lines carry the style name for 'hybrid'. Pair interactions are listed as 'PairIJ Coeffs' for
 * every i <= j, so no mixing-rule assumption is needed on the LAMMPS side.
 *
 * Settings a data file cannot carry (styles, cut-offs, kspace, special_bonds derived from the components'
 * 1-4 scaling) are written as a comment block at the top, ready to paste into the input script.
 *
 * \param components system component information
 * \param atomData holds all information on all atoms in the system.
 * \param simulationBox system simulation box (not unit cell).
 * \param forceField contains parameters for the pair interactions.
 * \param numberOfIntegerMoleculesPerComponent amount of molecules per component, necessary for accounting.
 * \param framework meta info on the framework.
 */
std::string WriteLAMMPSDataFile(std::span<const Component> components, std::span<const Atom> atomData,
                                std::span<const AtomDynamics> atomDynamics, std::span<const Molecule> moleculeData,
                                const SimulationBox simulationBox, const ForceField forceField,
                                std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                                std::optional<Framework> framework);
}  // namespace IO
