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
import vdwparameters;
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
          std::array<std::size_t, 2>{i, j}, VDWParameters::Type::LennardJones,
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

// A small cell with a Fixed group {0, 1}, a Rigid five-atom chain {2..6} and a Flexible atom {7}. Fixed atoms have no energy
// among themselves in any image, and neither has a Rigid group in its own geometry (the minimum-image pairs). A Rigid
// atom does interact with the replicas of its group: that distance changes when the group rotates.
TEST(framework_vdw_replicas, replicas_respect_fixed_and_rigid_groups)
{
  const double cutOff = 13.0;
  ForceField forceField({{"C", true, 12.0, 0.0, 0.0, 6, true}}, {{100.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, cutOff, cutOff, cutOff, false, false, false);
  const SimulationBox box(12.0, 12.0, 12.0);
  std::vector<Atom> atoms;
  // chains with ~1.45 A spacing, the bond detection connects the atoms of each group
  // (the chain ends are a nonbonded pair of the group, closer than half the cell)
  for (const double3& position : {double3(1.0, 1.0, 1.0), double3(2.4, 1.3, 1.0), double3(3.0, 6.0, 6.0),
                                  double3(4.4, 6.3, 6.0), double3(5.8, 5.8, 6.0), double3(7.2, 6.3, 6.0),
                                  double3(8.6, 6.0, 6.0), double3(1.5, 10.0, 10.0)})
  {
    atoms.push_back(Atom(position / 12.0, 0.0, 1.0, 0, 0, 0, 0, true));
  }
  Framework framework(forceField, "replica-groups", box, 1, atoms, atoms, int3(1, 1, 1));

  const std::filesystem::path definition = std::filesystem::temp_directory_path() / "raspa3-replica-groups.json";
  std::ofstream(definition) << R"({"Type": "Flexible",
                                   "Groups": [{"Type": "Fixed", "Atoms": [0, 1]},
                                              {"Type": "Rigid", "Atoms": [2, 3, 4, 5, 6]},
                                              {"Type": "Flexible", "Atoms": [7]}],
                                   "Intra14VanDerWaalsScalingValue": 0.0,
                                   "Intra14ChargeChargeScalingValue": 0.0})";
  framework.readFrameworkDefinition(forceField, definition.string());
  std::filesystem::remove(definition);
  ASSERT_TRUE(framework.isMixed());

  framework.regenerateVanDerWaalsImageList(forceField, box, cutOff);

  const auto isFixed = [](std::size_t atom) { return atom < 2; };
  const auto isRigid = [](std::size_t atom) { return atom >= 2 && atom < 7; };
  const auto minimumImageShift = [&](std::size_t i, std::size_t j)
  {
    const double3 raw = framework.atoms[i].position - framework.atoms[j].position;
    const double3 fractional = box.inverseCell * (raw - box.applyPeriodicBoundaryConditions(raw));
    return int3(static_cast<int>(std::llround(fractional.x)), static_cast<int>(std::llround(fractional.y)),
                static_cast<int>(std::llround(fractional.z)));
  };

  std::vector<std::size_t> rigidReplicaTerms;
  std::size_t fixedRigidTerms{0};
  const auto& terms = framework.intraMolecularPotentials.vanDerWaals;
  for (std::size_t index = 0; index < terms.size(); ++index)
  {
    const std::size_t i = terms[index].identifiers[0];
    const std::size_t j = terms[index].identifiers[1];
    const int3 shift = framework.intraMolecularImageShifts.vanDerWaals[index][1];
    EXPECT_FALSE(isFixed(i) && isFixed(j)) << "pair " << i << " " << j;
    if (i != j && isRigid(i) && isRigid(j))
    {
      const int3 minimum = minimumImageShift(i, j);
      EXPECT_FALSE(shift.x == minimum.x && shift.y == minimum.y && shift.z == minimum.z) << "pair " << i << " " << j;
    }
    if (isRigid(i) && isRigid(j) && i != j) rigidReplicaTerms.push_back(index);
    if ((isFixed(i) && isRigid(j)) || (isRigid(i) && isFixed(j))) ++fixedRigidTerms;
  }
  EXPECT_GT(fixedRigidTerms, 0uz);
  ASSERT_FALSE(rigidReplicaTerms.empty());

  // those replica terms are real: their energy depends on the orientation of the group
  const auto replicaEnergy = [&](const std::array<double3, 5>& positions)
  {
    double energy{0.0};
    for (const std::size_t index : rigidReplicaTerms)
    {
      const std::size_t i = terms[index].identifiers[0];
      const std::size_t j = terms[index].identifiers[1];
      const int3 shift = framework.intraMolecularImageShifts.vanDerWaals[index][1];
      const double3 dr = positions[i - 2] - (positions[j - 2] + box.cell * double3(static_cast<double>(shift.x),
                                                                                  static_cast<double>(shift.y),
                                                                                  static_cast<double>(shift.z)));
      energy += Potentials::potentialVDW<0>(forceField, 1.0, 1.0, double3::dot(dr, dr), 0, 0).energy;
    }
    return energy;
  };
  std::array<double3, 5> positions{};
  double3 center{};
  for (std::size_t k = 0; k < 5; ++k)
  {
    positions[k] = framework.atoms[k + 2].position;
    center = center + positions[k] / 5.0;
  }
  std::array<double3, 5> rotated{};
  const double cosine = std::cos(0.5);
  const double sine = std::sin(0.5);
  for (std::size_t k = 0; k < 5; ++k)
  {
    const double3 d = positions[k] - center;
    rotated[k] = center + double3(cosine * d.x - sine * d.y, sine * d.x + cosine * d.y, d.z);
  }
  EXPECT_GT(std::abs(replicaEnergy(rotated) - replicaEnergy(positions)), 1e-3);
}

// One atom per small cell: the generated home-cell pair list is empty (no pairs), but the atom interacts with its own
// periodic images within the cutoff. The regeneration returned on the empty list and dropped them all; it now returns
// only for a framework whose list was never generated (no connectivity table either).
TEST(framework_vdw_replicas, empty_generated_list_still_gets_the_periodic_images)
{
  const double cutOff = 10.0;
  ForceField forceField({{"C", false, 12.0, 0.0, 0.0, 6, false}}, {{100.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, cutOff, cutOff, cutOff, false, false, false);
  const SimulationBox box(4.0, 4.5, 5.0);
  const std::vector<Atom> atoms{Atom({0.1, 0.2, 0.3}, 0.0, 1.0, 0, 0, 0, 0, true)};
  Framework framework(forceField, "single-atom", box, 1, atoms, atoms, int3(1, 1, 1));
  framework.rigid = false;
  framework.generateIntraMolecularPotentials(forceField);
  ASSERT_TRUE(framework.intraMolecularPotentials.vanDerWaals.empty());

  framework.regenerateVanDerWaalsImageList(forceField, box, cutOff);
  const double replicaEnergy =
      Interactions::computeFrameworkIntraMolecularEnergy(forceField, framework, box, framework.atoms).intraVDW;
  const double referenceEnergy = bruteForceLatticeVanDerWaals(framework.atoms, forceField, box, cutOff, 4);
  EXPECT_LT(referenceEnergy, -1.0);
  EXPECT_NEAR(replicaEnergy, referenceEnergy, 1.0e-6 * (1.0 + std::abs(referenceEnergy)));
}

// A framework 1-4 pair takes the 1-4 parameters of the force field (as the 1-4 pairs of a molecule do), also with a
// 1-4 scaling of one (as in CHARMM); before, it used the regular table, and the 1-4 pair was recognized by a scaling
// other than one only.
TEST(framework_vdw_replicas, framework_one_four_pairs_use_the_one_four_parameters)
{
  ForceField forceField({{"C", false, 12.0, 0.0, 0.0, 6, false}}, {{100.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, false, false);
  forceField.setPair14SelfInteraction(0, VDWParameters(40.0, 3.4));
  forceField.preComputeCutOffDependentParameters();

  // a trans-planar butane-like chain, far from its periodic images: one 1-4 pair (0, 3), all others 1-2 or 1-3
  const SimulationBox box(30.0, 30.0, 30.0);
  const std::array<double3, 4> positions{double3(10.0, 10.0, 10.0), double3(11.26, 10.89, 10.0),
                                         double3(12.52, 10.0, 10.0), double3(13.78, 10.89, 10.0)};
  std::vector<Atom> atoms;
  for (const double3& position : positions) atoms.emplace_back(box.inverseCell * position, 0.0, 1.0, 0, 0, 0, 0, true);

  for (const double scaling : {1.0, 0.5})
  {
    Framework framework(forceField, "chain", box, 1, atoms, atoms, int3(1, 1, 1));
    framework.rigid = false;
    framework.intra14VanDerWaalsScaling = scaling;
    framework.generateIntraMolecularPotentials(forceField);
    ASSERT_EQ(framework.intraMolecularPotentials.vanDerWaals.size(), 1uz);

    const double3 dr = framework.atoms[0].position - framework.atoms[3].position;
    const double expected =
        scaling * Potentials::potentialVDW<0>(forceField.pair14(0, 0), 1.0, 1.0, double3::dot(dr, dr)).energy;
    const double energy =
        Interactions::computeFrameworkIntraMolecularEnergy(forceField, framework, box, framework.atoms).intraVDW;
    EXPECT_NEAR(energy, expected, 1e-10 * std::abs(expected)) << "scaling " << scaling;
    EXPECT_GT(std::abs(expected -
                       scaling * Potentials::potentialVDW<0>(forceField, 1.0, 1.0, double3::dot(dr, dr), 0, 0).energy),
              1.0);
  }
}
