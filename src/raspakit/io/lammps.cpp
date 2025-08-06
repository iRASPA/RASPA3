module;

#ifdef USE_LEGACY_HEADERS
#include <math>
#include <print>
#include <string>
#include <vector>
#endif

module lammps_io;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import component;
import atom;
import simulationbox;
import forcefield;
import pseudo_atom;
import units;
import framework;

std::string IO::WriteLAMMPSDataFile(std::span<const Component> components, std::span<const Atom> atomPositions,
                                    const SimulationBox simulationBox, const ForceField forceField,
                                    std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                                    std::optional<Framework> framework)
{
  std::ostringstream out;

  std::print(out, "LAMMPS Description\n\n");
  std::print(out, "{} atoms\n", atomPositions.size());

  std::print(out, "0 bonds\n");
  std::print(out, "0 angles\n");
  std::print(out, "0 dihedrals\n\n");

  std::print(out, "{} atom types\n", forceField.numberOfPseudoAtoms);
  std::print(out, "0 bond types\n");
  std::print(out, "0 angle types\n");
  std::print(out, "0 dihedral types\n");
  std::print(out, "0 improper types\n\n");

  double3 lengths = simulationBox.lengths();
  double3 angles = simulationBox.angles();

  double xy = lengths.y * std::cos(angles.z);
  double xz = lengths.z * std::cos(angles.y);
  double yz = (lengths.y * lengths.z * std::cos(angles.x) - xy * xz) / (lengths.y * std::sin(angles.z));

  std::print(out, "{} {} xlo xhi\n", 0.0, lengths.x);
  std::print(out, "{} {} ylo yhi\n", 0.0, lengths.y);
  std::print(out, "{} {} zlo zhi\n", 0.0, lengths.z);
  std::print(out, "{} {} {} xy xz yz\n", xy, xz, yz);

  std::print(out, "\nMasses\n\n");
  std::size_t idx = 1;  // lammps works with 1-indexing
  for (const PseudoAtom& pseudoAtom : forceField.pseudoAtoms)
  {
    std::print(out, "  {} {}\n", idx++, pseudoAtom.mass);
  }
  std::print(out, "\nPair Coeffs\n\n");
  idx = 1;  // lammps works with 1-indexing
  for (std::size_t i = 0; i < forceField.pseudoAtoms.size(); ++i)
  {
    std::print(out, "  {} {} {}\n", i + 1, forceField(i, i).parameters.x * Units::EnergyToKCalPerMol,
               forceField(i, i).parameters.y);
  }

  std::vector<std::size_t> molAtomOffset(numberOfIntegerMoleculesPerComponent.size());
  molAtomOffset[0] = 0;
  for (std::size_t i = 1; i < components.size(); ++i)
  {
    molAtomOffset[i] = molAtomOffset[i - 1] + numberOfIntegerMoleculesPerComponent[i - 1];
  }

  std::print(out, "\nAtoms\n\n");
  idx = 1;  // lammps works with 1-indexing
  std::size_t numberOfFrameworkAtoms = (framework.has_value()) ? framework->atoms.size() : 0uz;
  for (const Atom& atom : atomPositions)
  {
    std::size_t fwOffset = (framework.has_value() && idx < numberOfFrameworkAtoms + 1) ? 0 : 1;
    std::print(out, "  {} {} {} {} {} {} {}\n", idx++, atom.moleculeId + molAtomOffset[atom.componentId] + fwOffset + 1,
               atom.type + 1, atom.charge, atom.position.x, atom.position.y, atom.position.z);
  }

  std::print(out, "\nVelocities\n\n");
  idx = 1;  // lammps works with 1-indexing
  for (const Atom& atom : atomPositions)
  {
    std::print(out, "  {} {} {} {}\n", idx++, atom.velocity.x, atom.velocity.y, atom.velocity.z);
  }

  //   std::print(out, "\nBonds\n");
  //   idx = 1; // lammps works with 1-indexing
  //   std::size_t bondIndex = 1;
  //   for (const Component& component : components)
  // {
  //     for (std::size_t k = 0; k < numberOfIntegerMoleculesPerComponent[component.componentId]; ++k)
  //     {
  //         for (auto& bondPotential : component.internalPotentials.bonds)
  //         {
  //             std::size_t A = bondPotential.identifiers[0];
  //             std::size_t B = bondPotential.identifiers[1];
  //             std::print(out, "  {} {} {} {}\n", bondIndex, bondIndex++, A + idx , B + idx);
  //         }
  //         idx += component.definedAtoms.size();
  //     }
  // }

  return out.str();
}