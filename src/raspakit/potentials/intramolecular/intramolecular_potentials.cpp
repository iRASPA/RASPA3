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
#include <optional>
#include <ostream>
#include <print>
#include <source_location>
#include <span>
#include <sstream>
#include <vector>
#endif

module intra_molecular_potentials;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import atom;
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
import running_energy;

std::optional<BondPotential> Potentials::IntraMolecularPotentials::findBondPotential(std::size_t A, std::size_t B) const
{
  auto iterator = std::find_if(bonds.begin(), bonds.end(),
                               [A, B](const BondPotential &bond)
                               {
                                 return ((bond.identifiers[0] == A && bond.identifiers[1] == B) ||
                                         (bond.identifiers[0] == B && bond.identifiers[1] == A));
                               });
  if (iterator != bonds.end()) return *iterator;
  return std::nullopt;
}

double Potentials::IntraMolecularPotentials::calculateBondSmallMCEnergies(const std::span<Atom> atoms) const
{
  double energy{};

  for (const BondPotential &bond : bonds)
  {
    std::size_t A = bond.identifiers[0];
    std::size_t B = bond.identifiers[1];
    energy += bond.calculateEnergy(atoms[A].position, atoms[B].position);
  }

  return energy;
}

double Potentials::IntraMolecularPotentials::calculateBendSmallMCEnergies(const std::span<Atom> atoms) const
{
  double energy{};

  // for(const UreyBradleyPotential &ureyBradley : ureyBradleys)
  //{
  //   std::size_t A = ureyBradley.identifiers[0];
  //   std::size_t B = ureyBradley.identifiers[1];
  //   energies.ureyBradley += ureyBradley.calculateEnergy(atoms[A].position, atoms[B].position);
  // }

  for (const BendPotential &bend : bends)
  {
    std::size_t A = bend.identifiers[0];
    std::size_t B = bend.identifiers[1];
    std::size_t C = bend.identifiers[2];
    energy += bend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, std::nullopt);
  }

  // for(const InversionBendPotential &inversionBend : inversionBends)
  //{
  //   std::size_t A = inversionBend.identifiers[0];
  //   std::size_t B = inversionBend.identifiers[1];
  //   std::size_t C = inversionBend.identifiers[2];
  //   std::size_t D = inversionBend.identifiers[3];
  //   energies.inversionBend += inversionBend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position,
  //   atoms[D].position);
  // }

  // for(const OutOfPlaneBendPotential &outOfPlaneBend : outOfPlaneBends)
  //{
  //   std::size_t A = outOfPlaneBend.identifiers[0];
  //   std::size_t B = outOfPlaneBend.identifiers[1];
  //   std::size_t C = outOfPlaneBend.identifiers[2];
  //   std::size_t D = outOfPlaneBend.identifiers[3];
  //   energies.outOfPlaneBend += outOfPlaneBend.calculateEnergy(atoms[A].position, atoms[B].position,
  //   atoms[C].position, atoms[D].position);
  // }

  return energy;
}

double Potentials::IntraMolecularPotentials::calculateTorsionEnergies(const std::span<Atom> atoms) const
{
  double energy{};

  for (const TorsionPotential &torsion : torsions)
  {
    std::size_t A = torsion.identifiers[0];
    std::size_t B = torsion.identifiers[1];
    std::size_t C = torsion.identifiers[2];
    std::size_t D = torsion.identifiers[3];
    energy += torsion.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  return energy;
}

double Potentials::IntraMolecularPotentials::calculateVanDerWaalsEnergies(const std::span<Atom> atoms) const
{
  double energy{};

  for (const VanDerWaalsPotential &vanDerWaal : vanDerWaals)
  {
    std::size_t A = vanDerWaal.identifiers[0];
    std::size_t B = vanDerWaal.identifiers[1];
    energy += vanDerWaal.calculateEnergy(atoms[A].position, atoms[B].position);
  }

  return energy;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalEnergies(const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const BondPotential &bond : bonds)
  {
    std::size_t A = bond.identifiers[0];
    std::size_t B = bond.identifiers[1];
    energies.bond += bond.calculateEnergy(atoms[A].position, atoms[B].position);
  }

  for (const UreyBradleyPotential &ureyBradley : ureyBradleys)
  {
    std::size_t A = ureyBradley.identifiers[0];
    std::size_t B = ureyBradley.identifiers[1];
    energies.ureyBradley += ureyBradley.calculateEnergy(atoms[A].position, atoms[B].position);
  }

  for (const BendPotential &bend : bends)
  {
    std::size_t A = bend.identifiers[0];
    std::size_t B = bend.identifiers[1];
    std::size_t C = bend.identifiers[2];
    energies.bend += bend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, std::nullopt);
  }

  for (const InversionBendPotential &inversionBend : inversionBends)
  {
    std::size_t A = inversionBend.identifiers[0];
    std::size_t B = inversionBend.identifiers[1];
    std::size_t C = inversionBend.identifiers[2];
    std::size_t D = inversionBend.identifiers[3];
    energies.inversionBend +=
        inversionBend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  for (const OutOfPlaneBendPotential &outOfPlaneBend : outOfPlaneBends)
  {
    std::size_t A = outOfPlaneBend.identifiers[0];
    std::size_t B = outOfPlaneBend.identifiers[1];
    std::size_t C = outOfPlaneBend.identifiers[2];
    std::size_t D = outOfPlaneBend.identifiers[3];
    energies.outOfPlaneBend +=
        outOfPlaneBend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  for (const TorsionPotential &torsion : torsions)
  {
    std::size_t A = torsion.identifiers[0];
    std::size_t B = torsion.identifiers[1];
    std::size_t C = torsion.identifiers[2];
    std::size_t D = torsion.identifiers[3];
    energies.torsion +=
        torsion.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  for (const TorsionPotential &improperTorsion : improperTorsions)
  {
    std::size_t A = improperTorsion.identifiers[0];
    std::size_t B = improperTorsion.identifiers[1];
    std::size_t C = improperTorsion.identifiers[2];
    std::size_t D = improperTorsion.identifiers[3];
    energies.improperTorsion +=
        improperTorsion.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  for (const BondBondPotential &bondBond : bondBonds)
  {
    std::size_t A = bondBond.identifiers[0];
    std::size_t B = bondBond.identifiers[1];
    std::size_t C = bondBond.identifiers[2];
    energies.bondBond += bondBond.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position);
  }

  for (const BondBendPotential &bondBend : bondBends)
  {
    std::size_t A = bondBend.identifiers[0];
    std::size_t B = bondBend.identifiers[1];
    std::size_t C = bondBend.identifiers[2];
    std::size_t D = bondBend.identifiers[3];
    energies.bondBend +=
        bondBend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  for (const BondTorsionPotential &bondTorsion : bondTorsions)
  {
    std::size_t A = bondTorsion.identifiers[0];
    std::size_t B = bondTorsion.identifiers[1];
    std::size_t C = bondTorsion.identifiers[2];
    std::size_t D = bondTorsion.identifiers[3];
    energies.bondTorsion +=
        bondTorsion.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  for (const BendBendPotential &bendBend : bendBends)
  {
    std::size_t A = bendBend.identifiers[0];
    std::size_t B = bendBend.identifiers[1];
    std::size_t C = bendBend.identifiers[2];
    std::size_t D = bendBend.identifiers[3];
    energies.bendBend +=
        bendBend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  for (const BendTorsionPotential &bendTorsion : bendTorsions)
  {
    std::size_t A = bendTorsion.identifiers[0];
    std::size_t B = bendTorsion.identifiers[1];
    std::size_t C = bendTorsion.identifiers[2];
    std::size_t D = bendTorsion.identifiers[3];
    energies.bendTorsion +=
        bendTorsion.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  for (const VanDerWaalsPotential &vanDerWaal : vanDerWaals)
  {
    std::size_t A = vanDerWaal.identifiers[0];
    std::size_t B = vanDerWaal.identifiers[1];
    energies.intraVDW += vanDerWaal.calculateEnergy(atoms[A].position, atoms[B].position);
  }

  for (const CoulombPotential &coulomb : coulombs)
  {
    std::size_t A = coulomb.identifiers[0];
    std::size_t B = coulomb.identifiers[1];
    energies.intraCoul += coulomb.calculateEnergy(atoms[A].position, atoms[B].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalBondEnergies(const std::span<const Atom> atoms) const
{
  RunningEnergy energies;

  for (const BondPotential &bond : bonds)
  {
    std::size_t A = bond.identifiers[0];
    std::size_t B = bond.identifiers[1];
    energies.bond += bond.calculateEnergy(atoms[A].position, atoms[B].position);
  }
  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalUreyBradleyEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const UreyBradleyPotential &ureyBradley : ureyBradleys)
  {
    std::size_t A = ureyBradley.identifiers[0];
    std::size_t B = ureyBradley.identifiers[1];
    energies.ureyBradley += ureyBradley.calculateEnergy(atoms[A].position, atoms[B].position);
  }
  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalBendEnergies(const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const BendPotential &bend : bends)
  {
    std::size_t A = bend.identifiers[0];
    std::size_t B = bend.identifiers[1];
    std::size_t C = bend.identifiers[2];
    energies.bend += bend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, std::nullopt);
  }
  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalInversionBendEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const InversionBendPotential &inversionBend : inversionBends)
  {
    std::size_t A = inversionBend.identifiers[0];
    std::size_t B = inversionBend.identifiers[1];
    std::size_t C = inversionBend.identifiers[2];
    std::size_t D = inversionBend.identifiers[3];
    energies.inversionBend +=
        inversionBend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalOutOfPlaneBendEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const OutOfPlaneBendPotential &outOfPlaneBend : outOfPlaneBends)
  {
    std::size_t A = outOfPlaneBend.identifiers[0];
    std::size_t B = outOfPlaneBend.identifiers[1];
    std::size_t C = outOfPlaneBend.identifiers[2];
    std::size_t D = outOfPlaneBend.identifiers[3];
    energies.outOfPlaneBend +=
        outOfPlaneBend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalTorsionEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const TorsionPotential &torsion : torsions)
  {
    std::size_t A = torsion.identifiers[0];
    std::size_t B = torsion.identifiers[1];
    std::size_t C = torsion.identifiers[2];
    std::size_t D = torsion.identifiers[3];
    energies.torsion +=
        torsion.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalImproperTorsionEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const TorsionPotential &improperTorsion : improperTorsions)
  {
    std::size_t A = improperTorsion.identifiers[0];
    std::size_t B = improperTorsion.identifiers[1];
    std::size_t C = improperTorsion.identifiers[2];
    std::size_t D = improperTorsion.identifiers[3];
    energies.improperTorsion +=
        improperTorsion.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalBondBondEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const BondBondPotential &bondBond : bondBonds)
  {
    std::size_t A = bondBond.identifiers[0];
    std::size_t B = bondBond.identifiers[1];
    std::size_t C = bondBond.identifiers[2];
    energies.bondBond += bondBond.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalBondBendEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const BondBendPotential &bondBend : bondBends)
  {
    std::size_t A = bondBend.identifiers[0];
    std::size_t B = bondBend.identifiers[1];
    std::size_t C = bondBend.identifiers[2];
    std::size_t D = bondBend.identifiers[3];
    energies.bondBend +=
        bondBend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalBondTorsionEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const BondTorsionPotential &bondTorsion : bondTorsions)
  {
    std::size_t A = bondTorsion.identifiers[0];
    std::size_t B = bondTorsion.identifiers[1];
    std::size_t C = bondTorsion.identifiers[2];
    std::size_t D = bondTorsion.identifiers[3];
    energies.bondTorsion +=
        bondTorsion.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalBendBendEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const BendBendPotential &bendBend : bendBends)
  {
    std::size_t A = bendBend.identifiers[0];
    std::size_t B = bendBend.identifiers[1];
    std::size_t C = bendBend.identifiers[2];
    std::size_t D = bendBend.identifiers[3];
    energies.bendBend +=
        bendBend.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalBendTorsionEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const BendTorsionPotential &bendTorsion : bendTorsions)
  {
    std::size_t A = bendTorsion.identifiers[0];
    std::size_t B = bendTorsion.identifiers[1];
    std::size_t C = bendTorsion.identifiers[2];
    std::size_t D = bendTorsion.identifiers[3];
    energies.bendTorsion +=
        bendTorsion.calculateEnergy(atoms[A].position, atoms[B].position, atoms[C].position, atoms[D].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalIntraVanDerWaalsEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies{};

  for (const VanDerWaalsPotential &vanDerWaal : vanDerWaals)
  {
    std::size_t A = vanDerWaal.identifiers[0];
    std::size_t B = vanDerWaal.identifiers[1];
    energies.intraVDW += vanDerWaal.calculateEnergy(atoms[A].position, atoms[B].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalIntraCoulombEnergies(
    const std::span<const Atom> atoms) const
{
  RunningEnergy energies;

  for (const CoulombPotential &coulomb : coulombs)
  {
    std::size_t A = coulomb.identifiers[0];
    std::size_t B = coulomb.identifiers[1];
    energies.intraCoul += coulomb.calculateEnergy(atoms[A].position, atoms[B].position);
  }

  return energies;
}

RunningEnergy Potentials::IntraMolecularPotentials::computeInternalGradient(
    [[maybe_unused]] const std::span<Atom> atoms) const
{
  RunningEnergy energies{};

  return energies;
}

// compute the internal interactions that affect laying out the 'beadsToBePlaced'.
Potentials::IntraMolecularPotentials Potentials::IntraMolecularPotentials::filteredInteractions(
    std::size_t numberOfBeads, const std::span<std::size_t> beadsAlreadyPlaced,
    const std::span<std::size_t> beadsToBePlaced) const
{
  Potentials::IntraMolecularPotentials filteredPotentials{};
  std::vector<bool> boolToBePlaced(numberOfBeads);
  std::vector<bool> boolAlreadyPlacedToBePlaced(numberOfBeads);

  for (std::size_t const bead : beadsAlreadyPlaced)
  {
    boolAlreadyPlacedToBePlaced[bead] = true;
  }

  for (std::size_t const bead : beadsToBePlaced)
  {
    boolToBePlaced[bead] = true;
    boolAlreadyPlacedToBePlaced[bead] = true;
  }

  filteredPotentials.bonds.reserve(bonds.size());

  /*
  // guarantee that bond[i] corresponds to beadsToBePlaced[i]
  for (std::size_t i = 0; i != beadsToBePlaced.size(); ++i)
  {
    std::size_t A = currentBead;
    std::size_t B = beadsToBePlaced[i];
    const BondPotential & bond = bonds[i];
    if ((A == bond.identifiers[0] && B == bond.identifiers[1]) || (B == bond.identifiers[0] && A ==
  bond.identifiers[1]))
    {
      filteredPotentials.bonds.push_back(bond);
    }
  }
  */

  filteredPotentials.bonds.reserve(bonds.size());
  for (const BondPotential &bond : bonds)
  {
    std::size_t A = bond.identifiers[0];
    std::size_t B = bond.identifiers[1];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && (boolToBePlaced[A] || boolToBePlaced[B]))
    {
      filteredPotentials.bonds.push_back(bond);
    }
  }

  filteredPotentials.ureyBradleys.reserve(ureyBradleys.size());
  for (const UreyBradleyPotential &ureyBradley : ureyBradleys)
  {
    std::size_t A = ureyBradley.identifiers[0];
    std::size_t B = ureyBradley.identifiers[1];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && (boolToBePlaced[A] || boolToBePlaced[B]))
    {
      filteredPotentials.ureyBradleys.push_back(ureyBradley);
    }
  }

  filteredPotentials.bends.reserve(bends.size());
  for (const BendPotential &bend : bends)
  {
    std::size_t A = bend.identifiers[0];
    std::size_t B = bend.identifiers[1];
    std::size_t C = bend.identifiers[2];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && boolAlreadyPlacedToBePlaced[C] &&
        (boolToBePlaced[A] || boolToBePlaced[B] || boolToBePlaced[C]))
    {
      filteredPotentials.bends.push_back(bend);
    }
  }

  filteredPotentials.inversionBends.reserve(inversionBends.size());
  for (const InversionBendPotential &inversionBend : inversionBends)
  {
    std::size_t A = inversionBend.identifiers[0];
    std::size_t B = inversionBend.identifiers[1];
    std::size_t C = inversionBend.identifiers[2];
    std::size_t D = inversionBend.identifiers[3];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && boolAlreadyPlacedToBePlaced[C] &&
        boolAlreadyPlacedToBePlaced[D] &&
        (boolToBePlaced[A] || boolToBePlaced[B] || boolToBePlaced[C] || boolToBePlaced[D]))
    {
      filteredPotentials.inversionBends.push_back(inversionBend);
    }
  }

  filteredPotentials.outOfPlaneBends.reserve(outOfPlaneBends.size());
  for (const OutOfPlaneBendPotential &outOfPlaneBend : outOfPlaneBends)
  {
    std::size_t A = outOfPlaneBend.identifiers[0];
    std::size_t B = outOfPlaneBend.identifiers[1];
    std::size_t C = outOfPlaneBend.identifiers[2];
    std::size_t D = outOfPlaneBend.identifiers[3];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && boolAlreadyPlacedToBePlaced[C] &&
        boolAlreadyPlacedToBePlaced[D] &&
        (boolToBePlaced[A] || boolToBePlaced[B] || boolToBePlaced[C] || boolToBePlaced[D]))
    {
      filteredPotentials.outOfPlaneBends.push_back(outOfPlaneBend);
    }
  }

  filteredPotentials.torsions.reserve(torsions.size());
  for (const TorsionPotential &torsion : torsions)
  {
    std::size_t A = torsion.identifiers[0];
    std::size_t B = torsion.identifiers[1];
    std::size_t C = torsion.identifiers[2];
    std::size_t D = torsion.identifiers[3];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && boolAlreadyPlacedToBePlaced[C] &&
        boolAlreadyPlacedToBePlaced[D] &&
        (boolToBePlaced[A] || boolToBePlaced[B] || boolToBePlaced[C] || boolToBePlaced[D]))
    {
      filteredPotentials.torsions.push_back(torsion);
    }
  }

  filteredPotentials.improperTorsions.reserve(improperTorsions.size());
  for (const TorsionPotential &improperTorsion : improperTorsions)
  {
    std::size_t A = improperTorsion.identifiers[0];
    std::size_t B = improperTorsion.identifiers[1];
    std::size_t C = improperTorsion.identifiers[2];
    std::size_t D = improperTorsion.identifiers[3];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && boolAlreadyPlacedToBePlaced[C] &&
        boolAlreadyPlacedToBePlaced[D] &&
        (boolToBePlaced[A] || boolToBePlaced[B] || boolToBePlaced[C] || boolToBePlaced[D]))
    {
      filteredPotentials.improperTorsions.push_back(improperTorsion);
    }
  }

  filteredPotentials.bondBonds.reserve(bondBonds.size());
  for (const BondBondPotential &bondBond : bondBonds)
  {
    std::size_t A = bondBond.identifiers[0];
    std::size_t B = bondBond.identifiers[1];
    std::size_t C = bondBond.identifiers[2];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && boolAlreadyPlacedToBePlaced[C] &&
        (boolToBePlaced[A] || boolToBePlaced[B] || boolToBePlaced[C]))
    {
      filteredPotentials.bondBonds.push_back(bondBond);
    }
  }

  filteredPotentials.bondBends.reserve(bondBends.size());
  for (const BondBendPotential &bondBend : bondBends)
  {
    std::size_t A = bondBend.identifiers[0];
    std::size_t B = bondBend.identifiers[1];
    std::size_t C = bondBend.identifiers[2];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && boolAlreadyPlacedToBePlaced[C] &&
        (boolToBePlaced[A] || boolToBePlaced[B] || boolToBePlaced[C]))
    {
      filteredPotentials.bondBends.push_back(bondBend);
    }
  }

  filteredPotentials.bondTorsions.reserve(bondTorsions.size());
  for (const BondTorsionPotential &bondTorsion : bondTorsions)
  {
    std::size_t A = bondTorsion.identifiers[0];
    std::size_t B = bondTorsion.identifiers[1];
    std::size_t C = bondTorsion.identifiers[2];
    std::size_t D = bondTorsion.identifiers[3];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && boolAlreadyPlacedToBePlaced[C] &&
        boolAlreadyPlacedToBePlaced[D] &&
        (boolToBePlaced[A] || boolToBePlaced[B] || boolToBePlaced[C] || boolToBePlaced[D]))
    {
      filteredPotentials.bondTorsions.push_back(bondTorsion);
    }
  }

  filteredPotentials.bendBends.reserve(bendBends.size());
  for (const BendBendPotential &bendBend : bendBends)
  {
    std::size_t A = bendBend.identifiers[0];
    std::size_t B = bendBend.identifiers[1];
    std::size_t C = bendBend.identifiers[2];
    std::size_t D = bendBend.identifiers[3];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && boolAlreadyPlacedToBePlaced[C] &&
        boolAlreadyPlacedToBePlaced[D] &&
        (boolToBePlaced[A] || boolToBePlaced[B] || boolToBePlaced[C] || boolToBePlaced[D]))
    {
      filteredPotentials.bendBends.push_back(bendBend);
    }
  }

  filteredPotentials.bendTorsions.reserve(bendTorsions.size());
  for (const BendTorsionPotential &bendTorsion : bendTorsions)
  {
    std::size_t A = bendTorsion.identifiers[0];
    std::size_t B = bendTorsion.identifiers[1];
    std::size_t C = bendTorsion.identifiers[2];
    std::size_t D = bendTorsion.identifiers[3];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && boolAlreadyPlacedToBePlaced[C] &&
        boolAlreadyPlacedToBePlaced[D] &&
        (boolToBePlaced[A] || boolToBePlaced[B] || boolToBePlaced[C] || boolToBePlaced[D]))
    {
      filteredPotentials.bendTorsions.push_back(bendTorsion);
    }
  }

  filteredPotentials.vanDerWaals.reserve(vanDerWaals.size());
  for (const VanDerWaalsPotential &vanDerWaal : vanDerWaals)
  {
    std::size_t A = vanDerWaal.identifiers[0];
    std::size_t B = vanDerWaal.identifiers[1];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && (boolToBePlaced[A] || boolToBePlaced[B]))
    {
      filteredPotentials.vanDerWaals.push_back(vanDerWaal);
    }
  }

  filteredPotentials.coulombs.reserve(coulombs.size());
  for (const CoulombPotential &coulomb : coulombs)
  {
    std::size_t A = coulomb.identifiers[0];
    std::size_t B = coulomb.identifiers[1];
    if (boolAlreadyPlacedToBePlaced[A] && boolAlreadyPlacedToBePlaced[B] && (boolToBePlaced[A] || boolToBePlaced[B]))
    {
      filteredPotentials.coulombs.push_back(coulomb);
    }
  }

  return filteredPotentials;
}

std::string Potentials::IntraMolecularPotentials::printStatus() const
{
  std::ostringstream stream;

  // std::print(stream, "bonds: {}\n", bonds.size());
  // for(const BondPotential &bond : bonds)
  //{
  //   std::print("  {}\n", bond);
  // }
  // std::print(stream, "bends: {}\n", bends.size());
  // for(const BendPotential &bend : bends)
  //{
  //   std::print("  {}\n", bend);
  // }
  // std::print(stream, "torsions: {}\n", torsions.size());
  // for(const TorsionPotential &torsion : torsions)
  //{
  //   std::print("  {}\n", torsion);
  // }

  std::print(stream, "\n\n");

  return stream.str();
}

Archive<std::ofstream> &Potentials::operator<<(Archive<std::ofstream> &archive,
                                               const Potentials::IntraMolecularPotentials &p)
{
  archive << p.versionNumber;

  archive << p.chiralCenters;
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
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &Potentials::operator>>(Archive<std::ifstream> &archive, Potentials::IntraMolecularPotentials &p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;

  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'IntraMolecularPotentials' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> p.chiralCenters;
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
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Potentials::IntraMolecularPotentials: Error in binary restart\n"));
  }
#endif

  return archive;
}
