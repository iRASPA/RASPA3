#include <gtest/gtest.h>

import std;

import int3;
import double3;
import atom;
import component;
import framework;
import connectivity_table;
import intra_molecular_potentials;
import forcefield;
import system;
import simulationbox;
import interactions_ewald;
import normal_modes;
import phonon_dynamical_matrix;

namespace
{
// An eight-atom charged carbon chain host partitioned into Fixed{0} | Rigid{1, 2, 3, 4, 5} | Flexible{6, 7}, plus a
// charged single-site adsorbate. The rigid group holds 1-4 and 1-5 pairs (1-4, 2-5 and 1-5), which have no
// real-space Coulomb term, so their reciprocal-space energy has to be excluded like that of a rigid molecule.
System makeMixedFrameworkWithLongRigidGroup()
{
  ForceField forceField({{"C", true, 12.0, 0.0, 0.0, 6, true}, {"M", false, 16.0, 0.0, 0.0, 8, false}},
                        {{80.0, 3.0}, {60.0, 2.5}}, ForceField::MixingRule::Lorentz_Berthelot, 5.0, 5.0, 5.0, true,
                        false, true);
  const SimulationBox box(14.0, 14.0, 14.0);
  const std::uint16_t typeC = static_cast<std::uint16_t>(*forceField.findPseudoAtom("C"));
  const std::uint16_t typeM = static_cast<std::uint16_t>(*forceField.findPseudoAtom("M"));

  const std::array<double, 8> charges{0.3, -0.2, 0.25, -0.15, 0.2, -0.3, -0.2, -0.3};
  std::vector<Atom> atoms{};
  for (std::size_t i = 0; i < charges.size(); ++i)
  {
    // a zigzag in y so the rigid group is nonlinear, ~1.43-1.49 A between consecutive atoms
    const double3 position(3.0 + 1.4 * static_cast<double>(i), (i % 2 == 0) ? 7.0 : 7.3, 7.0);
    atoms.push_back(Atom(position / 14.0, charges[i], 1.0, 0, typeC, 0, 0, true));
  }
  Framework framework(forceField, "mixed-chain", box, 1, atoms, atoms, int3(1, 1, 1));

  const std::filesystem::path definition =
      std::filesystem::temp_directory_path() / "raspa3-mixed-framework-rigid-group-exclusions.json";
  {
    std::ofstream stream(definition);
    stream << R"({"Type": "Flexible",
                  "Groups": [
                    {"Type": "Fixed", "Atoms": [0]},
                    {"Type": "Rigid", "Atoms": [1, 2, 3, 4, 5]},
                    {"Type": "Flexible", "Atoms": [6, 7]}
                  ],
                  "Intra14VanDerWaalsScalingValue": 0.5,
                  "Intra14ChargeChargeScalingValue": 0.5,
                  "Bonds": [[["C", "C"], "HARMONIC", [1200.0, 1.45]]],
                  "Bends": [[["C", "C", "C"], "HARMONIC", [250.0, 160.0]]],
                  "Torsions": [[["C", "C", "C", "C"], "TRAPPE", [10.0, 20.0, 30.0, 40.0]]]})";
  }
  framework.readFrameworkDefinition(forceField, definition.string());
  std::filesystem::remove(definition);

  Component molecule(forceField, "single-site", 16.0, 1.0e6, 0.011,
                     {Atom({0.0, 0.0, 0.0}, 0.4, 1.0, 0, typeM, 0, false, false)}, ConnectivityTable(1),
                     Potentials::IntraMolecularPotentials{}, 0, 0);
  molecule.rigid = false;

  System system(forceField, box, false, 300.0, 1.0e4, 1.0, {framework}, {molecule}, {}, {1}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(12.4, 9.8, 7.6);
  return system;
}
}  // namespace

TEST(MIXED_FRAMEWORK_RIGID_GROUP_EXCLUSIONS, every_pair_inside_a_rigid_group_is_ewald_excluded)
{
  const System system = makeMixedFrameworkWithLongRigidGroup();
  ASSERT_TRUE(system.framework->isMixed());
  const std::set<std::array<std::size_t, 2>> excluded = Interactions::frameworkExcludedChargePairs(*system.framework);
  for (std::size_t i = 1; i <= 5; ++i)
  {
    for (std::size_t j = i + 1; j <= 5; ++j)
    {
      EXPECT_TRUE(excluded.contains({i, j})) << "rigid-group pair " << i << "-" << j;
    }
  }
}

TEST(MIXED_FRAMEWORK_RIGID_GROUP_EXCLUSIONS, charged_phonon_gamma_matches_normal_modes)
{
  // The phonon route flattens the mixed framework to independent Cartesian sites before computing the
  // force constants; the rigid-group pairs must keep their Ewald exclusion there, as in the generalized Hessian.
  System system = makeMixedFrameworkWithLongRigidGroup();
  const NormalModesResult modes = computeNormalModes(system);
  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1uz);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), modes.numberOfModes);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);
  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-6 * scale + 1e-10) << "mode " << mode;
  }
}
