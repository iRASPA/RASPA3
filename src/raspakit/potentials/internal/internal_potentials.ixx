module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#include <vector>
#endif

export module internal_potentials;

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

export namespace Potentials
{
struct InternalPotentials
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::vector<ChiralCenter> chiralCenters{};               ///< List of chiral centers in the component.
  std::vector<BondPotential> bonds{};                      ///< List of bond potentials.
  std::vector<UreyBradleyPotential> ureyBradleys{};        ///< List of Urey-Bradley potentials.
  std::vector<BendPotential> bends{};                      ///< List of bend potentials.
  std::vector<InversionBendPotential> inversionBends{};    ///< List of inversion-bend potentials.
  std::vector<OutOfPlaneBendPotential> outOfPlaneBends{};  ///< List of out-of-plane-bend potentials.
  std::vector<TorsionPotential> torsions{};                ///< List of torsion potentials.
  std::vector<TorsionPotential> improperTorsions{};        ///< List of improper torsion potentials.
  std::vector<BondBondPotential> bondBonds{};              ///< List of bond-bond potentials.
  std::vector<BondBendPotential> bondBends{};              ///< List of bond-bend potentials.
  std::vector<BondTorsionPotential> bondTorsions{};        ///< List of bond-torsion potentials.
  std::vector<BendBendPotential> bendBends{};              ///< List of bend-bend potentials.
  std::vector<BendTorsionPotential> bendTorsions{};        ///< List of bend-torsion potentials.
  std::vector<VanDerWaalsPotential> vanDerWaals{};         ///< List of vanDerWaals potentials.
  std::vector<CoulombPotential> coulombs{};                ///< List of Coulomb potentials.
  
  RunningEnergy computeInternalEnergies(const std::vector<Atom> &atoms) const;

  Potentials::InternalPotentials filteredInteractions(std::size_t numberOfBeads,
                                                      const std::vector<std::size_t>& beadsAlreadyPlaced,
                                                      const std::vector<std::size_t>& beadsToBePlaced) const;


  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Potentials::InternalPotentials& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Potentials::InternalPotentials& p);
};
}  // namespace Potentials
