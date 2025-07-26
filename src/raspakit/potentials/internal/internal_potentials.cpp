module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <map>
#include <print>
#include <source_location>
#include <vector>
#endif

module internal_potentials;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import chiral_center;
import bond_potential;
import urey_bradley_potential;
import bend_potential;
import inversion_bend_potential;
import out_of_plane_bend_potential;
import torsion_potential;
import bond_bond_potential;
import bond_bend_potential;
import bond_torsion_potential;
import bend_bend_potential;
import bend_torsion_potential;
import van_der_waals_potential;
import coulomb_potential;

// compute the internal interactions that affect laying out the 'beadsToBePlaced'.
Potentials::InternalPotentials Potentials::InternalPotentials::filteredInteractions(size_t numberOfBeads, 
    const std::vector<size_t> &beadsAlreadyPlaced, const std::vector<size_t> &beadsToBePlaced)
{
  Potentials::InternalPotentials filteredPotentials{};
  std::vector<bool> boolToBePlaced(numberOfBeads);
  std::vector<bool> boolAlreadyPlacedToBePlaced(numberOfBeads);

  for(size_t const bead : beadsAlreadyPlaced)
  {
    boolAlreadyPlacedToBePlaced[bead] = true;
  }

  for(size_t const bead : beadsToBePlaced)
  {
    boolToBePlaced[bead] = true;
    boolAlreadyPlacedToBePlaced[bead] = true;
  }

  filteredPotentials.bonds.reserve(bonds.size());
  for(const BondPotential &bond : bonds)
  {
    size_t A = bond.identifiers[0];
    size_t B = bond.identifiers[1];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B]))
    {
      filteredPotentials.bonds.push_back(bond);
    }
  }

  filteredPotentials.ureyBradleys.reserve(ureyBradleys.size());
  for(const UreyBradleyPotential &ureyBradley : ureyBradleys)
  {
    size_t A = ureyBradley.identifiers[0];
    size_t B = ureyBradley.identifiers[1];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B]))
    {
      filteredPotentials.ureyBradleys.push_back(ureyBradley);
    }
  }

  filteredPotentials.bends.reserve(bends.size());
  for(const BendPotential &bend : bends)
  {
    size_t A = bend.identifiers[0];
    size_t B = bend.identifiers[1];
    size_t C = bend.identifiers[2];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
       boolAlreadyPlacedToBePlaced[C] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B] ||
       boolToBePlaced[C]))
    {
      filteredPotentials.bends.push_back(bend);
    }
  }

  filteredPotentials.inversionBends.reserve(inversionBends.size());
  for(const InversionBendPotential &inversionBend : inversionBends)
  {
    size_t A = inversionBend.identifiers[0];
    size_t B = inversionBend.identifiers[1];
    size_t C = inversionBend.identifiers[2];
    size_t D = inversionBend.identifiers[3];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
       boolAlreadyPlacedToBePlaced[C] &&
       boolAlreadyPlacedToBePlaced[D] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B] ||
       boolToBePlaced[C] ||
       boolToBePlaced[D]))
    {
      filteredPotentials.inversionBends.push_back(inversionBend);
    }
  }

  filteredPotentials.outOfPlaneBends.reserve(outOfPlaneBends.size());
  for(const OutOfPlaneBendPotential &outOfPlaneBend : outOfPlaneBends)
  {
    size_t A = outOfPlaneBend.identifiers[0];
    size_t B = outOfPlaneBend.identifiers[1];
    size_t C = outOfPlaneBend.identifiers[2];
    size_t D = outOfPlaneBend.identifiers[3];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
       boolAlreadyPlacedToBePlaced[C] &&
       boolAlreadyPlacedToBePlaced[D] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B] ||
       boolToBePlaced[C] ||
       boolToBePlaced[D]))
    {
      filteredPotentials.outOfPlaneBends.push_back(outOfPlaneBend);
    }
  }

  filteredPotentials.torsions.reserve(torsions.size());
  for(const TorsionPotential &torsion : torsions)
  {
    size_t A = torsion.identifiers[0];
    size_t B = torsion.identifiers[1];
    size_t C = torsion.identifiers[2];
    size_t D = torsion.identifiers[3];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
       boolAlreadyPlacedToBePlaced[C] &&
       boolAlreadyPlacedToBePlaced[D] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B] ||
       boolToBePlaced[C] ||
       boolToBePlaced[D]))
    {
      filteredPotentials.torsions.push_back(torsion);
    }
  }

  filteredPotentials.improperTorsions.reserve(improperTorsions.size());
  for(const TorsionPotential &improperTorsion : improperTorsions)
  {
    size_t A = improperTorsion.identifiers[0];
    size_t B = improperTorsion.identifiers[1];
    size_t C = improperTorsion.identifiers[2];
    size_t D = improperTorsion.identifiers[3];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
       boolAlreadyPlacedToBePlaced[C] &&
       boolAlreadyPlacedToBePlaced[D] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B] ||
       boolToBePlaced[C] ||
       boolToBePlaced[D]))
    {
      filteredPotentials.improperTorsions.push_back(improperTorsion);
    }
  }

  filteredPotentials.bondBonds.reserve(bondBonds.size());
  for(const BondBondPotential &bondBond : bondBonds)
  {
    size_t A = bondBond.identifiers[0];
    size_t B = bondBond.identifiers[1];
    size_t C = bondBond.identifiers[2];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
       boolAlreadyPlacedToBePlaced[C] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B] ||
       boolToBePlaced[C]))
    {
      filteredPotentials.bondBonds.push_back(bondBond);
    }
  }

  filteredPotentials.bondBends.reserve(bondBends.size());
  for(const BondBendPotential &bondBend : bondBends)
  {
    size_t A = bondBend.identifiers[0];
    size_t B = bondBend.identifiers[1];
    size_t C = bondBend.identifiers[2];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
       boolAlreadyPlacedToBePlaced[C] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B] ||
       boolToBePlaced[C]))
    {
      filteredPotentials.bondBends.push_back(bondBend);
    }
  }

  filteredPotentials.bondTorsions.reserve(bondTorsions.size());
  for(const BondTorsionPotential &bondTorsion : bondTorsions)
  {
    size_t A = bondTorsion.identifiers[0];
    size_t B = bondTorsion.identifiers[1];
    size_t C = bondTorsion.identifiers[2];
    size_t D = bondTorsion.identifiers[3];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
       boolAlreadyPlacedToBePlaced[C] &&
       boolAlreadyPlacedToBePlaced[D] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B] ||
       boolToBePlaced[C] ||
       boolToBePlaced[D]))
    {
      filteredPotentials.bondTorsions.push_back(bondTorsion);
    }
  }

  filteredPotentials.bendBends.reserve(bendBends.size());
  for(const BendBendPotential &bendBend : bendBends)
  {
    size_t A = bendBend.identifiers[0];
    size_t B = bendBend.identifiers[1];
    size_t C = bendBend.identifiers[2];
    size_t D = bendBend.identifiers[3];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
       boolAlreadyPlacedToBePlaced[C] &&
       boolAlreadyPlacedToBePlaced[D] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B] ||
       boolToBePlaced[C] ||
       boolToBePlaced[D]))
    {
      filteredPotentials.bendBends.push_back(bendBend);
    }
  }

  filteredPotentials.bendTorsions.reserve(bendTorsions.size());
  for(const BendTorsionPotential &bendTorsion : bendTorsions)
  {
    size_t A = bendTorsion.identifiers[0];
    size_t B = bendTorsion.identifiers[1];
    size_t C = bendTorsion.identifiers[2];
    size_t D = bendTorsion.identifiers[3];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
       boolAlreadyPlacedToBePlaced[C] &&
       boolAlreadyPlacedToBePlaced[D] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B] ||
       boolToBePlaced[C] ||
       boolToBePlaced[D]))
    {
      filteredPotentials.bendTorsions.push_back(bendTorsion);
    }
  }

  filteredPotentials.vanDerWaals.reserve(vanDerWaals.size());
  for(const VanDerWaalsPotential &vanDerWaal : vanDerWaals)
  {
    size_t A = vanDerWaal.identifiers[0];
    size_t B = vanDerWaal.identifiers[1];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B]))
    {
      filteredPotentials.vanDerWaals.push_back(vanDerWaal);
    }
  }

  filteredPotentials.coulombs.reserve(coulombs.size());
  for(size_t i = 0; i != coulombs.size(); ++i)
  for(const CoulombPotential &coulomb : coulombs)
  {
    size_t A = coulomb.identifiers[0];
    size_t B = coulomb.identifiers[1];
    if(boolAlreadyPlacedToBePlaced[A] &&
       boolAlreadyPlacedToBePlaced[B] &&
      (boolToBePlaced[A] ||
       boolToBePlaced[B]))
    {
      filteredPotentials.coulombs.push_back(coulomb);
    }
  }

  return filteredPotentials;
}

Archive<std::ofstream> &Potentials::operator<<(Archive<std::ofstream> &archive, const Potentials::InternalPotentials &p)
{
  archive << p.versionNumber;

  archive << p.bonds;
  archive << p.ureyBradleys;
  archive << p.bends;
  archive << p.inversionBends;
  archive << p.outOfPlaneBends;
  archive << p.torsions;
  archive << p.improperTorsions;
  archive << p.bondBonds;
  archive << p.bondBends;
  archive << p.bondTorsions;
  archive << p.bendBends;
  archive << p.bendTorsions;
  archive << p.vanDerWaals;
  archive << p.coulombs;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &Potentials::operator>>(Archive<std::ifstream> &archive, Potentials::InternalPotentials &p)
{
  uint64_t versionNumber;
  archive >> versionNumber;

  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'InternalPotentials' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> p.bonds;
  archive >> p.ureyBradleys;
  archive >> p.bends;
  archive >> p.inversionBends;
  archive >> p.outOfPlaneBends;
  archive >> p.torsions;
  archive >> p.improperTorsions;
  archive >> p.bondBonds;
  archive >> p.bondBends;
  archive >> p.bondTorsions;
  archive >> p.bendBends;
  archive >> p.bendTorsions;
  archive >> p.vanDerWaals;
  archive >> p.coulombs;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Potentials::InternalPotentials: Error in binary restart\n"));
  }
#endif

  return archive;
}
