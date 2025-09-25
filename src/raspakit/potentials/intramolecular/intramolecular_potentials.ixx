module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <format>
#include <fstream>
#include <optional>
#include <print>
#include <span>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>
#endif

export module intra_molecular_potentials;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import double3x3;
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
struct IntraMolecularPotentials
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

  std::optional<BondPotential> findBondPotential(std::size_t A, std::size_t B) const;

  double calculateBondSmallMCEnergies(const std::span<Atom> atoms) const;
  double calculateBendSmallMCEnergies(const std::span<Atom> atoms) const;

  double calculateTorsionEnergies(const std::span<Atom> atoms) const;

  double calculateVanDerWaalsEnergies(const std::span<Atom> atoms) const;

  RunningEnergy computeInternalEnergies(const std::span<const Atom> atoms) const;

  RunningEnergy computeInternalBondEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalUreyBradleyEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBendEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalInversionBendEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalOutOfPlaneBendEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalTorsionEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalImproperTorsionEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBondBondEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBondBendEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBondTorsionEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBendBendEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalBendTorsionEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalIntraVanDerWaalsEnergies(const std::span<const Atom> atoms) const;
  RunningEnergy computeInternalIntraCoulombEnergies(const std::span<const Atom> atoms) const;

  RunningEnergy computeInternalGradient(const std::span<Atom> atoms) const;

  Potentials::IntraMolecularPotentials filteredInteractions(std::size_t numberOfBeads,
                                                            const std::span<std::size_t> beadsAlreadyPlaced,
                                                            const std::span<std::size_t> beadsToBePlaced) const;

  std::string printStatus() const;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive,
                                            const Potentials::IntraMolecularPotentials& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Potentials::IntraMolecularPotentials& p);
};
}  // namespace Potentials

/*
export template <>
struct std::formatter<Potentials::IntraMolecularPotentials>: std::formatter<string_view>
{
  auto format(const Potentials::IntraMolecularPotentials& p, std::format_context& ctx) const
  {
    std::string temp{};
    std::format_to(std::back_inserter(temp), "(bonds: {}-{})", p.bonds);
    return std::formatter<string_view>::format(temp, ctx);
  }
};
*/
