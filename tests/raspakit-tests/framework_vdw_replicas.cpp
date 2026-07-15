#include <gtest/gtest.h>

import std;

import units;
import int3;
import double3;
import double4;
import atom;
import forcefield;
import simulationbox;
import framework;
import van_der_waals_potential;
import intra_molecular_potentials;
import potential_pair_vdw;
import potential_pair_derivatives;
import interactions_internal;

namespace
{
// Seed the framework van der Waals pair list exactly as the construction path does: one entry per
// unordered atom pair at the minimum image. This makes the list non-empty (so the replica regeneration
// engages) and provides the minimum-image baseline for cells large enough not to need replicas.
void seedMinimumImageVanDerWaals(Framework& framework, const ForceField& forceField, const SimulationBox& box)
{
  framework.intraMolecularPotentials.vanDerWaals.clear();
  framework.intraMolecularImageShifts.vanDerWaals.clear();
  const std::vector<Atom>& atoms = framework.atoms;
  for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
  {
    for (std::size_t j = i + 1; j < atoms.size(); ++j)
    {
      const std::size_t typeA = static_cast<std::size_t>(atoms[i].type);
      const std::size_t typeB = static_cast<std::size_t>(atoms[j].type);
      const double4 parameters = forceField(typeA, typeB).parameters;
      framework.intraMolecularPotentials.vanDerWaals.emplace_back(
          std::array<std::size_t, 2>{i, j}, VanDerWaalsType::LennardJones,
          std::vector<double>{parameters.x * Units::EnergyToKelvin, parameters.y, parameters.z, parameters.w}, 1.0);

      const double3 raw = atoms[i].position - atoms[j].position;
      const double3 wrapped = box.applyPeriodicBoundaryConditions(raw);
      const double3 fractional = box.inverseCell * (raw - wrapped);
      const int3 shift(static_cast<int>(std::llround(fractional.x)), static_cast<int>(std::llround(fractional.y)),
                       static_cast<int>(std::llround(fractional.z)));
      framework.intraMolecularImageShifts.vanDerWaals.push_back({int3{}, shift});
    }
  }
}

// Independent reference: the full symmetric lattice sum E = 1/2 sum_{i,j,n, (i,n)!=(j,0)} U(r_i - r_j - L n)
// over every periodic image within the cutoff. Enumerated over all shells (positive and negative) with a
// 1/2 prefactor, i.e. a different code path from the implementation's (i<j, all images) + (self, positive
// half-shell) enumeration, so agreement validates the double-counting and self-image bookkeeping.
double bruteForceLatticeVanDerWaals(const std::vector<Atom>& atoms, const ForceField& forceField,
                                    const SimulationBox& box, double cutOff, int shells)
{
  const double cutOffSquared = cutOff * cutOff;
  double energy = 0.0;
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    for (std::size_t j = 0; j < atoms.size(); ++j)
    {
      const std::size_t typeA = static_cast<std::size_t>(atoms[i].type);
      const std::size_t typeB = static_cast<std::size_t>(atoms[j].type);
      for (int a = -shells; a <= shells; ++a)
      {
        for (int b = -shells; b <= shells; ++b)
        {
          for (int c = -shells; c <= shells; ++c)
          {
            if (i == j && a == 0 && b == 0 && c == 0) continue;
            const double3 dr = atoms[i].position - (atoms[j].position + box.cell * double3(static_cast<double>(a),
                                                                                           static_cast<double>(b),
                                                                                           static_cast<double>(c)));
            const double rr = double3::dot(dr, dr);
            if (rr < cutOffSquared)
            {
              energy += 0.5 * Potentials::potentialVDW<0>(forceField, 1.0, 1.0, rr, typeA, typeB).energy;
            }
          }
        }
      }
    }
  }
  return energy;
}

Framework makeSmallFramework(const ForceField& forceField, const SimulationBox& box)
{
  const std::vector<Atom> atoms{Atom({0.10, 0.15, 0.20}, 0.0, 1.0, 0, 0, 0, 0, true),
                                Atom({0.60, 0.35, 0.80}, 0.0, 1.0, 0, 0, 0, 0, true),
                                Atom({0.30, 0.70, 0.45}, 0.0, 1.0, 0, 0, 0, 0, true)};
  Framework framework(forceField, "replica-test", box, 1, atoms, atoms, int3(1, 1, 1));
  framework.rigid = false;
  return framework;
}
}  // namespace

// With a fixed cutoff larger than half the cell, the regenerated replica list must reproduce the exact
// periodic lattice van der Waals energy computed by an independent brute-force sum over image cells.
TEST(framework_vdw_replicas, replica_energy_matches_brute_force_lattice_sum)
{
  const double cutOff = 10.0;
  ForceField forceField({{"C", false, 12.0, 0.0, 0.0, 6, false}}, {{100.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, cutOff, cutOff, cutOff, false, false, false);
  const SimulationBox box(8.0, 8.0, 8.0);  // width 8 < 2 * cutOff, so a single minimum image is insufficient

  Framework framework = makeSmallFramework(forceField, box);
  seedMinimumImageVanDerWaals(framework, forceField, box);
  framework.regenerateVanDerWaalsImageList(forceField, box, cutOff);

  const double replicaEnergy =
      Interactions::computeFrameworkIntraMolecularEnergy(forceField, framework, box, framework.atoms).intraVDW;
  const double referenceEnergy = bruteForceLatticeVanDerWaals(framework.atoms, forceField, box, cutOff, 4);

  EXPECT_NEAR(replicaEnergy, referenceEnergy, 1.0e-6 * (1.0 + std::abs(referenceEnergy)));
}

// When the cell is at least twice the cutoff the minimum image is exact and regeneration must be a no-op,
// leaving the construction-time list (and hence the energy) untouched.
TEST(framework_vdw_replicas, large_cell_is_left_as_single_minimum_image)
{
  const double cutOff = 10.0;
  ForceField forceField({{"C", false, 12.0, 0.0, 0.0, 6, false}}, {{100.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, cutOff, cutOff, cutOff, false, false, false);
  const SimulationBox box(25.0, 25.0, 25.0);  // width 25 >= 2 * cutOff

  Framework framework = makeSmallFramework(forceField, box);
  seedMinimumImageVanDerWaals(framework, forceField, box);
  const std::size_t entriesBefore = framework.intraMolecularPotentials.vanDerWaals.size();

  framework.regenerateVanDerWaalsImageList(forceField, box, cutOff);

  EXPECT_EQ(framework.intraMolecularPotentials.vanDerWaals.size(), entriesBefore);

  const double replicaEnergy =
      Interactions::computeFrameworkIntraMolecularEnergy(forceField, framework, box, framework.atoms).intraVDW;
  const double referenceEnergy = bruteForceLatticeVanDerWaals(framework.atoms, forceField, box, cutOff, 2);
  EXPECT_NEAR(replicaEnergy, referenceEnergy, 1.0e-6 * (1.0 + std::abs(referenceEnergy)));
}
