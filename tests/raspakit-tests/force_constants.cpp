#include <gtest/gtest.h>

import std;

import double3;
import int3;
import atom;
import atom_dynamics;
import forcefield;
import component;
import framework;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import pseudo_atom;
import system;
import simulationbox;
import generalized_hessian;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import interactions_hessian_intermolecular;
import phonon_force_constants;

namespace
{
Atom fractionalAtom(double3 position, std::size_t type)
{
  return Atom(position, 0.0, 1.0, 0, static_cast<std::uint16_t>(type), 0, 0, true);
}

// Force field with four uncharged Lennard-Jones pseudo atoms; short cutoffs keep the minimum image
// unambiguous inside the test box.
ForceField makeFlexibleFrameworkForceField(double cutoff)
{
  return ForceField({{"C1", true, 12.0, 0.0, 0.0, 6, true},
                     {"C2", true, 12.0, 0.0, 0.0, 6, true},
                     {"C3", true, 12.0, 0.0, 0.0, 6, true},
                     {"H1", true, 1.0, 0.0, 0.0, 1, true}},
                    {{80.0, 3.0}, {70.0, 3.1}, {60.0, 3.2}, {40.0, 2.6}},
                    ForceField::MixingRule::Lorentz_Berthelot, cutoff, cutoff, cutoff, false, false, false);
}

class TemporaryFile
{
 public:
  TemporaryFile(std::string name, std::string_view contents)
      : path(std::filesystem::temp_directory_path() / std::move(name))
  {
    std::ofstream stream(path);
    stream << contents;
  }
  ~TemporaryFile()
  {
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
  }
  std::filesystem::path path;
};
}  // namespace

TEST(force_constants, image_resolved_fold_matches_intermolecular_hessian)
{
  // Flexible single-bead adsorbates with only VDW: the folded image-resolved force
  // constants must reproduce the analytic Cartesian Hessian block for block. Use a
  // small box and a pair straddling the periodic boundary so nonzero lattice images
  // are exercised.
  ForceField forceField = ForceField::makeZeoliteForceField(6.0, false, false, false);
  Component methane = Component::makeMethane(forceField, 0);
  methane.rigid = false;

  System system =
      System(forceField, SimulationBox(12.0, 12.0, 12.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {4}, 5);

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  ASSERT_EQ(atoms.size(), 4u);
  atoms[0].position = double3(1.0, 6.0, 6.0);
  atoms[1].position = double3(11.0, 6.0, 6.0);  // interacts with atom 0 across the x boundary
  atoms[2].position = double3(6.0, 6.0, 6.0);
  atoms[3].position = double3(6.0, 6.0, 9.5);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 12u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  Interactions::computeInterMolecularHessian(system.forceField, system.simulationBox, system.moleculeData,
                                             system.components, system.spanOfMoleculeAtoms(), layout, hessian,
                                             dynamics);

  const ForceConstants forceConstants = computeRealSpaceForceConstants(system);
  EXPECT_EQ(forceConstants.numberOfAtoms(), atoms.size());
  // At least the home image plus one wrapped image must be present.
  EXPECT_GT(forceConstants.numberOfImages(), 1u);

  const std::vector<double> folded = forceConstants.foldToGamma();
  ASSERT_EQ(folded.size(), layout.numDofs() * layout.numDofs());

  double maxAbs = 0.0;
  for (std::size_t index = 0; index < folded.size(); ++index)
  {
    maxAbs = std::max(maxAbs, std::abs(hessian.positionPosition()[index]));
  }
  ASSERT_GT(maxAbs, 0.0);

  for (std::size_t row = 0; row < layout.numDofs(); ++row)
  {
    for (std::size_t column = 0; column < layout.numDofs(); ++column)
    {
      const double reference = hessian.positionPosition()[row * layout.numDofs() + column];
      const double summed = folded[row * layout.numDofs() + column];
      EXPECT_NEAR(summed, reference, 1e-9 * std::max(1.0, maxAbs)) << "row=" << row << " column=" << column;
    }
  }
}

TEST(force_constants, home_image_holds_self_blocks)
{
  ForceField forceField = ForceField::makeZeoliteForceField(6.0, false, false, false);
  Component methane = Component::makeMethane(forceField, 0);
  methane.rigid = false;

  System system =
      System(forceField, SimulationBox(12.0, 12.0, 12.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  atoms[0].position = double3(1.0, 6.0, 6.0);
  atoms[1].position = double3(11.0, 6.0, 6.0);

  const ForceConstants forceConstants = computeRealSpaceForceConstants(system);
  const ForceConstants::BlockMap& blocks = forceConstants.blocks();

  // The single interacting pair wraps by one cell in -x for the (0->1) block.
  ASSERT_TRUE(blocks.contains(int3(0, 0, 0)));
  const auto image = std::ranges::find_if(blocks, [](const auto& entry) { return !(entry.first == int3(0, 0, 0)); });
  ASSERT_NE(image, blocks.end());
  EXPECT_EQ(image->first.y, 0);
  EXPECT_EQ(image->first.z, 0);
  EXPECT_EQ(std::abs(image->first.x), 1);
}

TEST(force_constants, image_resolved_fold_matches_flexible_framework_hessian)
{
  // Flexible framework (a periodic chain with bonds, bends and torsions crossing the cell boundary)
  // plus a flexible single-site adsorbate. The folded image-resolved force constants must reproduce
  // the full analytic position-position Hessian assembled by evaluateDerivatives: framework
  // intramolecular terms, framework-molecule pairs and (trivially) the intermolecular block.
  const double cutoff = 5.0;
  const SimulationBox box(12.0, 12.0, 12.0);
  ForceField forceField = makeFlexibleFrameworkForceField(cutoff);

  const std::size_t c1 = *forceField.findPseudoAtom("C1");
  const std::size_t c2 = *forceField.findPseudoAtom("C2");
  const std::size_t c3 = *forceField.findPseudoAtom("C3");
  const std::size_t h1 = *forceField.findPseudoAtom("H1");

  std::vector<Atom> frameworkAtoms{fractionalAtom({9.5 / 12.0, 4.8 / 12.0, 5.0 / 12.0}, c1),
                                   fractionalAtom({10.8 / 12.0, 5.1 / 12.0, 5.0 / 12.0}, c2),
                                   fractionalAtom({0.1 / 12.0, 4.9 / 12.0, 5.2 / 12.0}, c3),
                                   fractionalAtom({1.4 / 12.0, 5.2 / 12.0, 5.1 / 12.0}, h1)};
  Framework framework(forceField, "periodic-chain", box, 1, frameworkAtoms, frameworkAtoms, {1, 1, 1});
  TemporaryFile definition("raspa3-force-constants-chain.json",
                           R"({"Type": "Flexible",
          "Intra14VanDerWaalsScalingValue": 0.0,
          "Intra14ChargeChargeScalingValue": 0.0,
          "Bonds": [
            [["C1", "C2"], "HARMONIC", [200.0, 1.3]],
            [["C2", "C3"], "HARMONIC", [200.0, 1.3]],
            [["C3", "H1"], "HARMONIC", [200.0, 1.3]]
          ],
          "Bends": [
            [["C1", "C2", "C3"], "HARMONIC", [100.0, 150.0]],
            [["C2", "C3", "H1"], "HARMONIC", [100.0, 150.0]]
          ],
          "Torsions": [
            [["C1", "C2", "C3", "H1"], "TRAPPE", [10.0, 20.0, 30.0, 40.0]]
          ]})");
  framework.readFrameworkDefinition(forceField, definition.path.string());
  ASSERT_FALSE(framework.rigid);

  Component molecule(forceField, "single-site", 16.0, 1.0e6, 0.011,
                     {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(c1), 0, false, false)},
                     ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  molecule.rigid = false;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {framework}, {molecule}, {}, {1}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(10.8, 6.3, 5.0);  // ~1.2 A from framework atom 1

  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, system.spanOfFrameworkAtoms().size());
  ASSERT_EQ(layout.numDofs(), 15u);  // 4 flexible framework atoms + 1 flexible molecule atom

  GeneralizedHessian reference(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = true};
  DerivativeResults results{.gradient = gradient, .hessian = reference};
  evaluateDerivatives(system, layout, capabilities, results);

  const ForceConstants forceConstants = computeRealSpaceForceConstants(system);
  EXPECT_EQ(forceConstants.numberOfAtoms(), 5u);
  EXPECT_GT(forceConstants.numberOfImages(), 1u);

  const std::vector<double> folded = forceConstants.foldToGamma();
  ASSERT_EQ(folded.size(), layout.numDofs() * layout.numDofs());

  double maxAbs = 0.0;
  for (const double value : reference.positionPosition()) maxAbs = std::max(maxAbs, std::abs(value));
  ASSERT_GT(maxAbs, 0.0);

  for (std::size_t row = 0; row < layout.numDofs(); ++row)
  {
    for (std::size_t column = 0; column < layout.numDofs(); ++column)
    {
      const double expected = reference.positionPosition()[row * layout.numDofs() + column];
      const double summed = folded[row * layout.numDofs() + column];
      EXPECT_NEAR(summed, expected, 1e-8 * std::max(1.0, maxAbs)) << "row=" << row << " column=" << column;
    }
  }
}

TEST(force_constants, image_resolved_fold_matches_molecule_intramolecular_hessian)
{
  // Two flexible bonded dimers: the folded image-resolved force constants must reproduce the analytic
  // Hessian including the intramolecular harmonic bonds (home cell) and the intermolecular VDW pairs
  // (including a pair wrapped across the periodic boundary).
  ForceField forceField({{"C", false, 12.0, 0.0, 0.0, 6, false}}, {{120.0, 3.4}},
                        ForceField::MixingRule::Lorentz_Berthelot, 6.0, 6.0, 6.0, false, false, false);
  const SimulationBox box(12.0, 12.0, 12.0);

  Potentials::IntraMolecularPotentials intra{};
  intra.bonds.emplace_back(std::array<std::size_t, 2>{0, 1}, BondType::Harmonic, std::vector<double>{1000.0, 1.2});
  const std::size_t carbon = *forceField.findPseudoAtom("C");
  ConnectivityTable connectivity(2);
  connectivity[0, 1] = true;
  connectivity[1, 0] = true;
  Component dimer(forceField, "dimer", 24.0, 1.0e6, 0.1,
                  {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(carbon), 0, false, false),
                   Atom({1.2, 0.0, 0.0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(carbon), 0, false, false)},
                  connectivity, intra, 0, 0);
  dimer.rigid = false;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {}, {dimer}, {}, {2}, 5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  ASSERT_EQ(atoms.size(), 4u);
  atoms[0].position = double3(1.0, 6.0, 6.0);
  atoms[1].position = double3(2.2, 6.0, 6.0);
  atoms[2].position = double3(11.0, 6.0, 6.0);  // atom 2 interacts with atom 0 across the x boundary
  atoms[3].position = double3(9.8, 6.0, 6.0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 12u);

  GeneralizedHessian reference(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = true};
  DerivativeResults results{.gradient = gradient, .hessian = reference};
  evaluateDerivatives(system, layout, capabilities, results);

  const ForceConstants forceConstants = computeRealSpaceForceConstants(system);
  const std::vector<double> folded = forceConstants.foldToGamma();
  ASSERT_EQ(folded.size(), layout.numDofs() * layout.numDofs());

  double maxAbs = 0.0;
  for (const double value : reference.positionPosition()) maxAbs = std::max(maxAbs, std::abs(value));
  ASSERT_GT(maxAbs, 0.0);

  for (std::size_t row = 0; row < layout.numDofs(); ++row)
  {
    for (std::size_t column = 0; column < layout.numDofs(); ++column)
    {
      const double expected = reference.positionPosition()[row * layout.numDofs() + column];
      const double summed = folded[row * layout.numDofs() + column];
      EXPECT_NEAR(summed, expected, 1e-8 * std::max(1.0, maxAbs)) << "row=" << row << " column=" << column;
    }
  }
}

TEST(force_constants, rejects_rigid_molecules)
{
  ForceField forceField = ForceField::makeZeoliteForceField(6.0, false, false, false);
  Component methane = Component::makeMethane(forceField, 0);  // rigid by default

  System system =
      System(forceField, SimulationBox(12.0, 12.0, 12.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5);

  EXPECT_THROW(computeRealSpaceForceConstants(system), std::runtime_error);
}
