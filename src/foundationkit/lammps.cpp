module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <format>
#include <iostream>
#include <ostream>
#include <print>
#include <span>
#include <sstream>
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
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;

std::string IO::WriteLAMMPSDataFile(std::span<const Component> components, std::span<const Atom> atomData,
                                    std::span<const Molecule> moleculeData, const SimulationBox simulationBox,
                                    const ForceField forceField,
                                    std::vector<std::size_t> numberOfIntegerMoleculesPerComponent,
                                    std::optional<Framework> framework)
{
  std::ostringstream out;

  std::print(out, "LAMMPS Description\n\n");
  std::print(out, "{} atoms\n", atomData.size());

  std::size_t numberOfBonds = 0uz;
  std::size_t numberOfBends = 0uz;
  std::size_t numberOfTorsions = 0uz;
  std::size_t numberOfBondPotentials = 0uz;
  std::size_t numberOfBendPotentials = 0uz;
  std::size_t numberOfTorsionPotentials = 0uz;
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    numberOfBonds += components[i].intraMolecularPotentials.bonds.size() * numberOfIntegerMoleculesPerComponent[i];
    numberOfBends += components[i].intraMolecularPotentials.bends.size() * numberOfIntegerMoleculesPerComponent[i];
    numberOfTorsions +=
        components[i].intraMolecularPotentials.torsions.size() * numberOfIntegerMoleculesPerComponent[i];
    numberOfBondPotentials += components[i].intraMolecularPotentials.bonds.size();
    numberOfBendPotentials += components[i].intraMolecularPotentials.bends.size();
    numberOfTorsionPotentials += components[i].intraMolecularPotentials.torsions.size();
  }

  std::print(out, "{} bonds\n", numberOfBonds);
  std::print(out, "{} angles\n", numberOfBends);
  std::print(out, "{} dihedrals\n\n", numberOfTorsions);

  std::print(out, "{} atom types\n", forceField.numberOfPseudoAtoms);
  std::print(out, "{} bond types\n", numberOfBondPotentials);
  std::print(out, "{} angle types\n", numberOfBendPotentials);
  std::print(out, "{} dihedral types\n", numberOfTorsionPotentials);
  std::print(out, "0 improper types\n\n");

  double3 lengths = simulationBox.lengths();
  double3 angles = simulationBox.angles();

  double xy = lengths.y * std::cos(angles.z);
  double xz = lengths.z * std::cos(angles.y);
  double yz = (lengths.y * lengths.z * std::cos(angles.x) - xy * xz) / (lengths.y * std::sin(angles.z));

  std::print(out, "{} {} xlo xhi\n", 0.0, lengths.x);
  std::print(out, "{} {} ylo yhi\n", 0.0, lengths.y);
  std::print(out, "{} {} zlo zhi\n", 0.0, lengths.z);
  if ((xy > 1e-10) || (xz > 1e-10) || (yz > 1e-10))
  {
    std::print(out, "{} {} {} xy xz yz\n", xy, xz, yz);
  }

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

  if (numberOfBondPotentials)
  {
    std::print(out, "\nBond Coeffs\n\n");
    idx = 1;
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      for (auto& bond : components[i].intraMolecularPotentials.bonds)
      {
        switch (bond.type)
        {
          case BondType::Harmonic:
          {
            std::print(out, "{} {} {}\n", idx++, bond.parameters[0] * Units::EnergyToKCalPerMol, bond.parameters[1]);
            break;
          }
          default:
          {
            break;
          }
        }
      }
    }
  }

  if (numberOfBendPotentials)
  {
    std::print(out, "\nAngle Coeffs\n\n");
    idx = 1;
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      for (auto& bend : components[i].intraMolecularPotentials.bends)
      {
        switch (bend.type)
        {
          case BendType::Harmonic:
          {
            std::print(out, "{} {} {}\n", idx++, bend.parameters[0] * Units::EnergyToKCalPerMol,
                       bend.parameters[1] * Units::RadiansToDegrees);
            break;
          }
          default:
          {
            break;
          }
        }
        // now only supporting Harmonic
      }
    }
  }

  if (numberOfTorsionPotentials)
  {
    std::print(out, "\nDihedral Coeffs\n\n");
    idx = 1;
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      for (auto& torsion : components[i].intraMolecularPotentials.torsions)
      {
        switch (torsion.type)
        {
          case TorsionType::Harmonic:
          {
            std::print(out, "{} {} {}\n", idx++, torsion.parameters[0] * Units::EnergyToKCalPerMol,
                       torsion.parameters[1] * Units::RadiansToDegrees);
            break;
          }
          default:
          {
            break;
          }
        }
        // now only supporting Harmonic
      }
    }
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
  for (const Atom& atom : atomData)
  {
    std::size_t fwOffset = (framework.has_value() && idx < numberOfFrameworkAtoms + 1) ? 0 : 1;
    std::print(out, "  {} {} {} {} {} {} {}\n", idx++, atom.moleculeId + molAtomOffset[atom.componentId] + fwOffset + 1,
               atom.type + 1, atom.charge, atom.position.x, atom.position.y, atom.position.z);
  }

  std::print(out, "\nVelocities\n\n");
  idx = 1;  // lammps works with 1-indexing
  for (const Atom& atom : atomData)
  {
    if (components[atom.componentId].growType == Component::GrowType::Flexible)
    {
      std::print(out, "  {} {} {} {}\n", idx++, atom.velocity.x, atom.velocity.y, atom.velocity.z);
    }
    else
    {
      double3 molVelocity = moleculeData[atom.moleculeId].velocity;
      std::print(out, "  {} {} {} {}\n", idx++, molVelocity.x, molVelocity.y, molVelocity.z);
    }
  }

  if (numberOfBonds)
  {
    std::print(out, "\nBonds\n\n");
    idx = 1;  // lammps works with 1-indexing
    std::size_t bondIndex = 1;
    std::size_t bondCoeffCount = 1;
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      for (std::size_t k = 0; k < numberOfIntegerMoleculesPerComponent[i]; ++k)
      {
        for (std::size_t b = 0; b < components[i].intraMolecularPotentials.bonds.size(); b++)
        {
          auto& bond = components[i].intraMolecularPotentials.bonds[b];
          std::print(out, "  {} {} {} {}\n", bondIndex++, bondCoeffCount + b, bond.identifiers[0] + idx,
                     bond.identifiers[1] + idx);
        }
        idx += components[i].definedAtoms.size();
      }
      bondCoeffCount += components[i].intraMolecularPotentials.bonds.size();
    }
  }

  if (numberOfBends)
  {
    std::print(out, "\nAngles\n\n");
    idx = 1;  // lammps works with 1-indexing
    std::size_t bendIndex = 1;
    std::size_t bendCoeffCount = 1;
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      for (std::size_t k = 0; k < numberOfIntegerMoleculesPerComponent[i]; ++k)
      {
        for (std::size_t b = 0; b < components[i].intraMolecularPotentials.bends.size(); b++)
        {
          auto& bend = components[i].intraMolecularPotentials.bends[b];
          std::print(out, "  {} {} {} {} {}\n", bendIndex++, bendCoeffCount + b, bend.identifiers[0] + idx,
                     bend.identifiers[1] + idx, bend.identifiers[2] + idx);
        }
        idx += components[i].definedAtoms.size();
      }
      bendCoeffCount += components[i].intraMolecularPotentials.bends.size();
    }
  }

  if (numberOfTorsions)
  {
    std::print(out, "\nDihedrals\n\n");
    idx = 1;  // lammps works with 1-indexing
    std::size_t torsionIndex = 1;
    std::size_t torsionCoeffCount = 1;
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      for (std::size_t k = 0; k < numberOfIntegerMoleculesPerComponent[i]; ++k)
      {
        for (std::size_t t = 0; t < components[i].intraMolecularPotentials.torsions.size(); t++)
        {
          auto& torsion = components[i].intraMolecularPotentials.torsions[t];
          std::print(out, "  {} {} {} {} {} {}\n", torsionIndex++, torsionCoeffCount + t, torsion.identifiers[0] + idx,
                     torsion.identifiers[1] + idx, torsion.identifiers[2] + idx, torsion.identifiers[3] + idx);
        }
        idx += components[i].definedAtoms.size();
      }
      torsionCoeffCount += components[i].intraMolecularPotentials.torsions.size();
    }
  }

  return out.str();
}
