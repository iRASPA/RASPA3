#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "irmof_fixtures.hpp"

import std;

import archive;
import atom;
import atom_dynamics;
import cif_reader;
import double3;
import forcefield;
import framework;
import generalized_hessian;
import int3;
import interactions_hessian_intramolecular;
import interactions_internal;
import minimization_dof_layout;
import pseudo_atom;
import running_energy;
import simulationbox;
import van_der_waals_potential;
import vdwparameters;

namespace
{

ForceField makeFlexibleIrmofForceField()
{
  TemporaryDirectory dir;
  dir.write("force_field.json", irmof_fixtures::kFlexibleIrmofForceFieldJson);
  return ForceField::readForceField(dir.path().string(), "force_field.json").value();
}

ForceField makeTestForceField(double cutoff = 12.0)
{
  return ForceField({{"C1", true, 12.0, 0.0, 0.0, 6, true},
                     {"C2", true, 12.0, 0.0, 0.0, 6, true},
                     {"C3", true, 12.0, 0.0, 0.0, 6, true},
                     {"H1", true, 1.0, 0.0, 0.0, 1, true}},
                    {{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}}, ForceField::MixingRule::Lorentz_Berthelot, cutoff,
                    cutoff, cutoff, false, false, false);
}

Atom fractionalAtom(double3 position, std::size_t type)
{
  return Atom(position, 0.0, 1.0, 0, static_cast<std::uint16_t>(type), 0, 0, true);
}

Framework makePeriodicPairFramework(const ForceField& forceField)
{
  std::vector<Atom> atoms{fractionalAtom({0.01, 0.5, 0.5}, *forceField.findPseudoAtom("C1")),
                          fractionalAtom({0.99, 0.5, 0.5}, *forceField.findPseudoAtom("C2"))};
  return Framework(forceField, "periodic-pair", SimulationBox(10.0, 10.0, 10.0), 1, atoms, atoms, {1, 1, 1});
}

template <std::size_t N>
std::array<std::string, N> canonicalTypes(const std::array<std::size_t, N>& indices, const Framework& framework,
                                          const ForceField& forceField)
{
  std::array<std::string, N> forward;
  std::array<std::string, N> reverse;
  for (std::size_t i = 0; i < N; ++i)
  {
    forward[i] = forceField.pseudoAtoms[framework.atoms[indices[i]].type].name;
    reverse[N - i - 1] = forward[i];
  }
  return std::min(forward, reverse);
}
}  // namespace

TEST(framework_intramolecular_reader, periodic_connectivity_and_reversed_type_matching)
{
  ForceField forceField = makeTestForceField();
  Framework framework = makePeriodicPairFramework(forceField);
  TemporaryFile definition("raspa3-framework-periodic-definition.json",
                           R"({"Type": "Flexible",
                               "Bonds": [[["C2", "C1"], "HARMONIC", [100.0, 1.0]]]})");

  framework.readFrameworkDefinition(forceField, definition.path().string());

  ASSERT_EQ(framework.connectivityTable.findAllBonds().size(), 1);
  ASSERT_EQ(framework.intraMolecularPotentials.bonds.size(), 1);
  EXPECT_FALSE(framework.rigid);
  EXPECT_EQ(framework.intraMolecularPotentials.bonds[0].identifiers, (std::array<std::size_t, 2>{0, 1}));
  ASSERT_EQ(framework.intraMolecularImageShifts.bonds.size(), 1);
  EXPECT_EQ(framework.intraMolecularImageShifts.bonds[0][0], int3(0, 0, 0));
  EXPECT_EQ(framework.intraMolecularImageShifts.bonds[0][1], int3(-1, 0, 0));

  const std::string status = framework.printStatus(forceField);
  EXPECT_NE(status.find("framework type: Flexible\n    number of bonds: 1"), std::string::npos);
  EXPECT_NE(status.find("C1 - C2 : HARMONIC p_0/k_B=100 [K/Å^2], p_1=1 [Å] (1 found)"), std::string::npos);
  EXPECT_NE(status.find("number of bend potentials"), std::string::npos);
  EXPECT_NE(status.find("number of torsion potentials"), std::string::npos);
  EXPECT_NE(status.find("number of intra-framework Van der Waals potentials"), std::string::npos);
  EXPECT_NE(status.find("number of intra-framework Coulomb potentials"), std::string::npos);
}

TEST(framework_intramolecular_reader, rigid_framework_skips_internal_derivatives_and_dofs)
{
  ForceField forceField = makeTestForceField();
  Framework framework = makePeriodicPairFramework(forceField);
  TemporaryFile definition("raspa3-framework-rigid-definition.json",
                           R"({"Type": "Rigid",
                               "Bonds": [[["C1", "C2"], "HARMONIC", [100.0, 1.0]]]})");
  framework.readFrameworkDefinition(forceField, definition.path().string());
  ASSERT_TRUE(framework.rigid);
  ASSERT_FALSE(framework.intraMolecularPotentials.bonds.empty());
  EXPECT_NE(framework.printStatus(forceField).find("framework type: Rigid\n    number of bonds: 1"), std::string::npos);

  const std::size_t numberOfFlexibleFrameworkAtoms = framework.rigid ? 0 : framework.atoms.size();
  const MinimizationDofLayout layout = buildMinimizationDofLayout({}, {}, numberOfFlexibleFrameworkAtoms);
  EXPECT_EQ(layout.numDofs(), 0);
  EXPECT_FALSE(layout.frameworkAtomDof(0, MinimizationDofAxis::X).has_value());

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<AtomDynamics> dynamics(framework.atoms.size());
  const RunningEnergy energy = Interactions::computeFrameworkIntraMolecularHessian(
      forceField, framework, SimulationBox(10.0, 10.0, 10.0), framework.atoms, layout, hessian, dynamics);
  EXPECT_DOUBLE_EQ(energy.potentialEnergy(), 0.0);
  EXPECT_TRUE(std::ranges::all_of(dynamics,
                                  [](const AtomDynamics& value) { return value.gradient == double3(0.0, 0.0, 0.0); }));
}

TEST(framework_intramolecular_reader, periodic_bond_gradient_and_hessian_match_finite_difference)
{
  ForceField forceField = makeTestForceField();
  Framework framework = makePeriodicPairFramework(forceField);
  TemporaryFile definition("raspa3-framework-periodic-derivatives.json",
                           R"({"Type": "Flexible",
                               "Bonds": [[["C1", "C2"], "HARMONIC", [100.0, 1.0]]]})");
  framework.readFrameworkDefinition(forceField, definition.path().string());

  const SimulationBox box(10.0, 10.0, 10.0);
  std::vector<Atom> atoms = framework.atoms;
  std::vector<AtomDynamics> dynamics(atoms.size());
  const MinimizationDofLayout layout = buildMinimizationDofLayout({}, {}, atoms.size());
  GeneralizedHessian hessian(layout.numDofs(), 0);
  const RunningEnergy analytic =
      Interactions::computeFrameworkIntraMolecularHessian(forceField, framework, box, atoms, layout, hessian, dynamics);
  EXPECT_DOUBLE_EQ(analytic.bond,
                   Interactions::computeFrameworkIntraMolecularEnergy(forceField, framework, box, atoms).bond);

  const auto energy = [&]()
  { return Interactions::computeFrameworkIntraMolecularEnergy(forceField, framework, box, atoms).bond; };
  constexpr double delta = 1.0e-5;
  constexpr double gradientTolerance = 1.0e-5;
  constexpr double hessianTolerance = 2.0e-2;
  const double referenceEnergy = energy();

  for (std::size_t row = 0; row < layout.numDofs(); ++row)
  {
    const std::size_t atomA = row / 3;
    const std::size_t axisA = row % 3;
    (&atoms[atomA].position.x)[axisA] += delta;
    const double energyPlus = energy();
    (&atoms[atomA].position.x)[axisA] -= 2.0 * delta;
    const double energyMinus = energy();
    (&atoms[atomA].position.x)[axisA] += delta;
    EXPECT_NEAR((&dynamics[atomA].gradient.x)[axisA], (energyPlus - energyMinus) / (2.0 * delta), gradientTolerance);

    for (std::size_t column = 0; column < layout.numDofs(); ++column)
    {
      const std::size_t atomB = column / 3;
      const std::size_t axisB = column % 3;
      double numerical{};
      if (row == column)
      {
        numerical = (energyPlus - 2.0 * referenceEnergy + energyMinus) / (delta * delta);
      }
      else
      {
        (&atoms[atomA].position.x)[axisA] += delta;
        (&atoms[atomB].position.x)[axisB] += delta;
        const double ePP = energy();
        (&atoms[atomB].position.x)[axisB] -= 2.0 * delta;
        const double ePM = energy();
        (&atoms[atomA].position.x)[axisA] -= 2.0 * delta;
        const double eMM = energy();
        (&atoms[atomB].position.x)[axisB] += 2.0 * delta;
        const double eMP = energy();
        (&atoms[atomA].position.x)[axisA] += delta;
        (&atoms[atomB].position.x)[axisB] -= delta;
        numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
      }
      EXPECT_NEAR(hessian(row, column), numerical, hessianTolerance);
    }
  }
}

TEST(framework_intramolecular_reader, framework_vdw_respects_framework_cutoff)
{
  ForceField forceField = makeTestForceField(0.5);
  Framework framework = makePeriodicPairFramework(forceField);
  framework.rigid = false;
  framework.intraMolecularPotentials.vanDerWaals = {
      VanDerWaalsPotential({0, 1}, VanDerWaalsType::LennardJones, {1.0, 1.0}, 1.0)};
  framework.intraMolecularImageShifts.vanDerWaals = {{{int3{}, int3{}}}};

  const SimulationBox box(10.0, 10.0, 10.0);
  std::vector<Atom> atoms = framework.atoms;
  atoms[0].position = double3(1.0, 1.0, 1.0);
  atoms[1].position = double3(1.6, 1.0, 1.0);
  EXPECT_DOUBLE_EQ(Interactions::computeFrameworkIntraMolecularEnergy(forceField, framework, box, atoms).intraVDW, 0.0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout({}, {}, atoms.size());
  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<AtomDynamics> dynamics(atoms.size());
  const RunningEnergy outside =
      Interactions::computeFrameworkIntraMolecularHessian(forceField, framework, box, atoms, layout, hessian, dynamics);
  EXPECT_DOUBLE_EQ(outside.intraVDW, 0.0);
  EXPECT_TRUE(std::ranges::all_of(hessian.positionPosition(), [](double value) { return value == 0.0; }));
  EXPECT_TRUE(std::ranges::all_of(dynamics,
                                  [](const AtomDynamics& value) { return value.gradient == double3(0.0, 0.0, 0.0); }));

  atoms[1].position.x = 1.4;
  EXPECT_NE(Interactions::computeFrameworkIntraMolecularEnergy(forceField, framework, box, atoms).intraVDW, 0.0);
  GeneralizedHessian insideHessian(layout.numDofs(), 0);
  std::vector<AtomDynamics> insideDynamics(atoms.size());
  Interactions::computeFrameworkIntraMolecularHessian(forceField, framework, box, atoms, layout, insideHessian,
                                                      insideDynamics);
  EXPECT_NEAR(insideHessian.strainGradient().m11,
              (atoms[0].position.x - atoms[1].position.x) * insideDynamics[0].gradient.x, 1.0e-12);
}

TEST(framework_intramolecular_reader, periodic_bend_and_torsion_hessian_match_finite_difference)
{
  ForceField forceField = makeTestForceField();
  const SimulationBox box(12.0, 12.0, 12.0);
  std::vector<Atom> fractionalAtoms{
      fractionalAtom({9.5 / 12.0, 4.8 / 12.0, 5.0 / 12.0}, *forceField.findPseudoAtom("C1")),
      fractionalAtom({10.8 / 12.0, 5.1 / 12.0, 5.0 / 12.0}, *forceField.findPseudoAtom("C2")),
      fractionalAtom({0.1 / 12.0, 4.9 / 12.0, 5.2 / 12.0}, *forceField.findPseudoAtom("C3")),
      fractionalAtom({1.4 / 12.0, 5.2 / 12.0, 5.1 / 12.0}, *forceField.findPseudoAtom("H1"))};
  Framework framework(forceField, "periodic-chain", box, 1, fractionalAtoms, fractionalAtoms, {1, 1, 1});
  TemporaryFile definition("raspa3-framework-periodic-chain.json",
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
  framework.readFrameworkDefinition(forceField, definition.path().string());
  ASSERT_EQ(framework.intraMolecularImageShifts.bends.size(), 2);
  ASSERT_EQ(framework.intraMolecularImageShifts.torsions.size(), 1);
  EXPECT_EQ(framework.intraMolecularImageShifts.torsions[0][2], int3(1, 0, 0));
  EXPECT_EQ(framework.intraMolecularImageShifts.torsions[0][3], int3(1, 0, 0));

  std::vector<Atom> atoms = framework.atoms;
  const MinimizationDofLayout layout = buildMinimizationDofLayout({}, {}, atoms.size());
  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<AtomDynamics> referenceDynamics(atoms.size());
  Interactions::computeFrameworkIntraMolecularHessian(forceField, framework, box, atoms, layout, hessian,
                                                      referenceDynamics);

  constexpr double delta = 2.0e-5;
  constexpr double tolerance = 5.0e-2;
  for (std::size_t column = 0; column < layout.numDofs(); ++column)
  {
    const std::size_t atom = column / 3;
    const std::size_t axis = column % 3;
    (&atoms[atom].position.x)[axis] += delta;
    GeneralizedHessian plusHessian(layout.numDofs(), 0);
    std::vector<AtomDynamics> plusDynamics(atoms.size());
    Interactions::computeFrameworkIntraMolecularHessian(forceField, framework, box, atoms, layout, plusHessian,
                                                        plusDynamics);

    (&atoms[atom].position.x)[axis] -= 2.0 * delta;
    GeneralizedHessian minusHessian(layout.numDofs(), 0);
    std::vector<AtomDynamics> minusDynamics(atoms.size());
    Interactions::computeFrameworkIntraMolecularHessian(forceField, framework, box, atoms, layout, minusHessian,
                                                        minusDynamics);
    (&atoms[atom].position.x)[axis] += delta;

    for (std::size_t row = 0; row < layout.numDofs(); ++row)
    {
      const double numerical =
          ((&plusDynamics[row / 3].gradient.x)[row % 3] - (&minusDynamics[row / 3].gradient.x)[row % 3]) /
          (2.0 * delta);
      EXPECT_NEAR(hessian(row, column), numerical, tolerance) << "row=" << row << " column=" << column;
    }
  }
}

TEST(framework_intramolecular_reader, converted_irmof_definition_expands_all_templates)
{
  TemporaryDirectory fixtureDir;
  fixtureDir.write("force_field.json", irmof_fixtures::kFlexibleIrmofForceFieldJson);
  fixtureDir.write("framework.json", irmof_fixtures::kFlexibleIrmofFrameworkJson);

  ForceField forceField = ForceField::readForceField(fixtureDir.path().string(), "force_field.json").value();

  const auto cif = CIFReader::readCIFString(std::string(irmof_fixtures::kIrmof1Cif), forceField,
                                            CIFReader::UseChargesFrom::PseudoAtoms);
  ASSERT_TRUE(cif.has_value());
  auto [box, hallNumber, definedAtoms, fractionalAtoms] = *cif;
  Framework framework(forceField, "IRMOF-1", box, hallNumber, definedAtoms, fractionalAtoms, {1, 1, 1});

  framework.readFrameworkDefinition(forceField, (fixtureDir.path() / "framework.json").string());
  forceField.initializeEwaldParameters(framework.simulationBox.scaled(framework.numberOfUnitCells));

  const MinimizationDofLayout diagnosticLayout = buildMinimizationDofLayout({}, {}, framework.atoms.size());
  GeneralizedHessian diagnosticHessian(diagnosticLayout.numDofs(), 0);
  std::vector<AtomDynamics> diagnosticDynamics(framework.atoms.size());
  const RunningEnergy diagnosticEnergy = Interactions::computeFrameworkIntraMolecularHessian(
      forceField, framework, framework.simulationBox.scaled(framework.numberOfUnitCells), framework.atoms,
      diagnosticLayout, diagnosticHessian, diagnosticDynamics);
  EXPECT_NEAR(diagnosticEnergy.bond, 126391.76547486844, 1.0e-6);
  EXPECT_NEAR(diagnosticEnergy.bend, 11632.676845790833, 1.0e-6);
  EXPECT_NEAR(diagnosticEnergy.torsion, 0.0, 1.0e-6);
  EXPECT_NEAR(diagnosticEnergy.improperTorsion, 0.0, 1.0e-6);
  EXPECT_NEAR(diagnosticEnergy.intraVDW, 651417.2173727077, 1.0e-6);
  EXPECT_NEAR(diagnosticEnergy.intraCoul, -2218955.250157769, 1.0e-6);

  const auto& potentials = framework.intraMolecularPotentials;
  EXPECT_FALSE(potentials.bonds.empty());
  EXPECT_FALSE(potentials.bends.empty());
  EXPECT_FALSE(potentials.torsions.empty());
  EXPECT_FALSE(potentials.improperTorsions.empty());
  EXPECT_EQ(potentials.bonds.size(), 384);
  EXPECT_EQ(potentials.bends.size(), 576);
  EXPECT_EQ(potentials.torsions.size(), 768);
  EXPECT_EQ(potentials.improperTorsions.size(), 192);
  EXPECT_FALSE(potentials.vanDerWaals.empty());
  EXPECT_FALSE(potentials.coulombs.empty());
  const std::size_t expectedNonBondedCount = 88716;
  EXPECT_EQ(potentials.vanDerWaals.size(), expectedNonBondedCount);
  EXPECT_EQ(potentials.coulombs.size(), expectedNonBondedCount);
  EXPECT_TRUE(
      std::ranges::all_of(potentials.bonds, [](const auto& potential) { return potential.parameters[0] > 0.0; }));
  EXPECT_TRUE(
      std::ranges::all_of(potentials.bends, [](const auto& potential) { return potential.parameters[0] > 0.0; }));
  EXPECT_TRUE(
      std::ranges::all_of(potentials.torsions, [](const auto& potential) { return potential.parameters[2] > 0.0; }));
  EXPECT_TRUE(std::ranges::all_of(potentials.improperTorsions,
                                  [](const auto& potential) { return potential.parameters[2] > 0.0; }));

  std::set<std::array<std::string, 2>> bondTemplates;
  std::set<std::array<std::string, 3>> bendTemplates;
  std::set<std::array<std::string, 4>> torsionTemplates;
  std::set<std::array<std::string, 4>> improperTemplates;
  for (const auto& potential : potentials.bonds)
    bondTemplates.insert(canonicalTypes(potential.identifiers, framework, forceField));
  for (const auto& potential : potentials.bends)
    bendTemplates.insert(canonicalTypes(potential.identifiers, framework, forceField));
  for (const auto& potential : potentials.torsions)
    torsionTemplates.insert(canonicalTypes(potential.identifiers, framework, forceField));
  for (const auto& potential : potentials.improperTorsions)
  {
    std::array<std::string, 3> outer{forceField.pseudoAtoms[framework.atoms[potential.identifiers[0]].type].name,
                                     forceField.pseudoAtoms[framework.atoms[potential.identifiers[2]].type].name,
                                     forceField.pseudoAtoms[framework.atoms[potential.identifiers[3]].type].name};
    std::ranges::sort(outer);
    improperTemplates.insert(
        {outer[0], forceField.pseudoAtoms[framework.atoms[potential.identifiers[1]].type].name, outer[1], outer[2]});
  }

  EXPECT_EQ(bondTemplates.size(), 5);
  EXPECT_EQ(bendTemplates.size(), 7);
  EXPECT_EQ(torsionTemplates.size(), 8);
  EXPECT_EQ(improperTemplates.size(), 3);
}

TEST(framework_intramolecular_reader, applies_independent_14_scaling_factors)
{
  ForceField forceField = makeFlexibleIrmofForceField();

  const std::size_t carbon = *forceField.findPseudoAtom("C3");
  const std::size_t hydrogen = *forceField.findPseudoAtom("H1");
  std::vector<Atom> atoms{fractionalAtom({0.100, 0.5, 0.5}, carbon), fractionalAtom({0.165, 0.5, 0.5}, hydrogen),
                          fractionalAtom({0.230, 0.5, 0.5}, carbon), fractionalAtom({0.295, 0.5, 0.5}, hydrogen)};
  Framework framework(forceField, "scaled-chain", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms, {1, 1, 1});
  TemporaryFile definition("raspa3-framework-14-scaling.json",
                           R"({"Intra14VanDerWaalsScalingValue": 0.5,
          "Intra14ChargeChargeScalingValue": 0.75,
          "Bonds": [[["C3", "H1"], "HARMONIC", [100.0, 1.3]]]})");

  framework.readFrameworkDefinition(forceField, definition.path().string());

  ASSERT_EQ(framework.connectivityTable.findAllTorsions().size(), 1);
  ASSERT_EQ(framework.intraMolecularPotentials.vanDerWaals.size(), 1);
  ASSERT_EQ(framework.intraMolecularPotentials.coulombs.size(), 1);
  const auto scaledVDW = std::ranges::find_if(framework.intraMolecularPotentials.vanDerWaals, [](const auto& potential)
                                               { return potential.identifiers == std::array<std::size_t, 2>{0, 3}; });
  const auto scaledCoulomb = std::ranges::find_if(framework.intraMolecularPotentials.coulombs, [](const auto& potential)
                                                   { return potential.identifiers == std::array<std::size_t, 2>{0, 3}; });
  ASSERT_NE(scaledVDW, framework.intraMolecularPotentials.vanDerWaals.end());
  ASSERT_NE(scaledCoulomb, framework.intraMolecularPotentials.coulombs.end());
  EXPECT_DOUBLE_EQ(scaledVDW->scaling, 0.5);
  EXPECT_DOUBLE_EQ(scaledCoulomb->scaling, 0.75);
}

TEST(framework_intramolecular_reader, controls_12_and_13_exclusions_independently)
{
  ForceField forceField = makeFlexibleIrmofForceField();

  const std::size_t carbon = *forceField.findPseudoAtom("C3");
  const std::size_t hydrogen = *forceField.findPseudoAtom("H1");
  std::vector<Atom> atoms{fractionalAtom({0.100, 0.5, 0.5}, carbon), fractionalAtom({0.165, 0.5, 0.5}, hydrogen),
                          fractionalAtom({0.230, 0.5, 0.5}, carbon), fractionalAtom({0.295, 0.5, 0.5}, hydrogen)};

  Framework include12(forceField, "include-12", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms, {1, 1, 1});
  TemporaryFile include12Definition("raspa3-framework-include-12.json",
                                    R"({"ExcludeIntra12Interactions": false,
          "Bonds": [[["C3", "H1"], "HARMONIC", [100.0, 1.3]]]})");
  include12.readFrameworkDefinition(forceField, include12Definition.path().string());
  EXPECT_EQ(include12.intraMolecularPotentials.vanDerWaals.size(), 3);
  EXPECT_EQ(include12.intraMolecularPotentials.coulombs.size(), 3);

  Framework excludeBonds(forceField, "exclude-bonds", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms, {1, 1, 1});
  TemporaryFile excludeBondsDefinition("raspa3-framework-exclude-bonds.json",
                                       R"({"ExcludeIntra12Interactions": false,
          "ExcludeIntraBondInteractions": true,
          "Bonds": [[["C3", "H1"], "HARMONIC", [100.0, 1.3]]]})");
  excludeBonds.readFrameworkDefinition(forceField, excludeBondsDefinition.path().string());
  EXPECT_TRUE(excludeBonds.intraMolecularPotentials.vanDerWaals.empty());
  EXPECT_TRUE(excludeBonds.intraMolecularPotentials.coulombs.empty());

  Framework include13(forceField, "include-13", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms, {1, 1, 1});
  TemporaryFile include13Definition("raspa3-framework-include-13.json",
                                    R"({"ExcludeIntra13Interactions": false,
          "Bonds": [[["C3", "H1"], "HARMONIC", [100.0, 1.3]]]})");
  include13.readFrameworkDefinition(forceField, include13Definition.path().string());
  EXPECT_EQ(include13.intraMolecularPotentials.vanDerWaals.size(), 2);
  EXPECT_EQ(include13.intraMolecularPotentials.coulombs.size(), 2);

  Framework excludeBends(forceField, "exclude-bends", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms, {1, 1, 1});
  TemporaryFile excludeBendsDefinition(
      "raspa3-framework-exclude-bends.json",
      R"({"ExcludeIntra13Interactions": false,
          "ExcludeIntraBendInteractions": true,
          "Bonds": [[["C3", "H1"], "HARMONIC", [100.0, 1.3]]],
          "Bends": [[["C3", "H1", "C3"], "HARMONIC", [100.0, 120.0]],
                    [["H1", "C3", "H1"], "HARMONIC", [100.0, 120.0]]]})");
  excludeBends.readFrameworkDefinition(forceField, excludeBendsDefinition.path().string());
  EXPECT_TRUE(excludeBends.intraMolecularPotentials.vanDerWaals.empty());
  EXPECT_TRUE(excludeBends.intraMolecularPotentials.coulombs.empty());
}

TEST(framework_intramolecular_reader, rejects_unknown_ambiguous_and_malformed_definitions)
{
  ForceField forceField = makeTestForceField();
  Framework framework = makePeriodicPairFramework(forceField);

  TemporaryFile unknown("raspa3-framework-unknown-definition.json",
                        R"({"Bonds": [[["missing", "C1"], "HARMONIC", [100.0, 1.0]]]})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, unknown.path().string()), std::runtime_error);

  TemporaryFile ambiguous("raspa3-framework-ambiguous-definition.json",
                          R"({"Bonds": [[["C1", "C2"], "HARMONIC", [100.0, 1.0]],
                    [["C2", "C1"], "HARMONIC", [100.0, 1.0]]]})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, ambiguous.path().string()), std::runtime_error);

  TemporaryFile malformed("raspa3-framework-malformed-definition.json",
                          R"({"Bonds": [[["C1", "C2"], "HARMONIC", [100.0]]]})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, malformed.path().string()), std::runtime_error);

  TemporaryFile invalidScaling("raspa3-framework-invalid-scaling.json", R"({"Intra14VanDerWaalsScalingValue": -0.5})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, invalidScaling.path().string()), std::runtime_error);

  TemporaryFile invalidExclusion("raspa3-framework-invalid-exclusion.json", R"({"ExcludeIntra12Interactions": 1})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, invalidExclusion.path().string()), std::runtime_error);

  TemporaryFile invalidType("raspa3-framework-invalid-type.json", R"({"Type": "deformable"})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, invalidType.path().string()), std::runtime_error);
}

TEST(framework_intramolecular_reader, mixed_groups_reorder_fixed_first_and_skip_internal_bonds)
{
  ForceField forceField = makeTestForceField();
  // Linear chain with ~1.5 Å bonds (within covalent-radius detection in a 20 Å box).
  std::vector<Atom> atoms{fractionalAtom({0.40, 0.5, 0.5}, *forceField.findPseudoAtom("C1")),
                          fractionalAtom({0.475, 0.5, 0.5}, *forceField.findPseudoAtom("C2")),
                          fractionalAtom({0.55, 0.5, 0.5}, *forceField.findPseudoAtom("C3"))};
  Framework framework(forceField, "mixed-chain", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms, {1, 1, 1});

  TemporaryFile definition("raspa3-framework-mixed-groups.json",
                           R"({"Type": "Flexible",
                               "Groups": [
                                 {"Type": "Fixed", "Atoms": [0]},
                                 {"Type": "Flexible", "Atoms": [1, 2]}
                               ],
                               "Bonds": [
                                 [["C1", "C2"], "HARMONIC", [100.0, 1.5]],
                                 [["C2", "C3"], "HARMONIC", [100.0, 1.5]]
                               ]})");
  framework.readFrameworkDefinition(forceField, definition.path().string());

  ASSERT_TRUE(framework.isMixed());
  EXPECT_FALSE(framework.rigid);
  EXPECT_EQ(framework.numberOfFixedAtoms(), 1);
  EXPECT_EQ(framework.numberOfMobileAtoms(), 2);
  EXPECT_TRUE(framework.isFixedAtom(0));
  EXPECT_TRUE(framework.isFlexibleAtom(1));
  EXPECT_TRUE(framework.isFlexibleAtom(2));

  // After Fixed-first reorder: Fixed–Flexible and Flexible–Flexible bonds are both kept.
  ASSERT_EQ(framework.intraMolecularPotentials.bonds.size(), 2);

  const MinimizationDofLayout layout =
      buildMinimizationDofLayout({}, {}, framework.atoms.size(), 0, &framework);
  EXPECT_FALSE(layout.frameworkAtomDof(0, MinimizationDofAxis::X).has_value());
  EXPECT_TRUE(layout.frameworkAtomDof(1, MinimizationDofAxis::X).has_value());
  EXPECT_TRUE(layout.frameworkAtomDof(2, MinimizationDofAxis::X).has_value());
  EXPECT_EQ(layout.numDofs(), 6);
}

TEST(framework_intramolecular_reader, mixed_rigid_group_skips_internal_bond)
{
  ForceField forceField = makeTestForceField();
  std::vector<Atom> atoms{fractionalAtom({0.40, 0.5, 0.5}, *forceField.findPseudoAtom("C1")),
                          fractionalAtom({0.475, 0.5, 0.5}, *forceField.findPseudoAtom("C2")),
                          fractionalAtom({0.55, 0.5, 0.5}, *forceField.findPseudoAtom("C3"))};
  Framework framework(forceField, "mixed-rigid", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms, {1, 1, 1});

  TemporaryFile definition("raspa3-framework-mixed-rigid.json",
                           R"({"Type": "Flexible",
                               "Groups": [
                                 {"Type": "Rigid", "Atoms": [0, 1]},
                                 {"Type": "Flexible", "Atoms": [2]}
                               ],
                               "Bonds": [
                                 [["C1", "C2"], "HARMONIC", [100.0, 1.5]],
                                 [["C2", "C3"], "HARMONIC", [100.0, 1.5]]
                               ]})");
  framework.readFrameworkDefinition(forceField, definition.path().string());

  ASSERT_TRUE(framework.isMixed());
  EXPECT_EQ(framework.numberOfRigidGroups(), 1);
  EXPECT_EQ(framework.numberOfFixedAtoms(), 0);
  // Internal C1–C2 bond inside the Rigid group is skipped; only the bridging C2–C3 remains.
  ASSERT_EQ(framework.intraMolecularPotentials.bonds.size(), 1);

  const MinimizationDofLayout layout =
      buildMinimizationDofLayout({}, {}, framework.atoms.size(), 0, &framework);
  // Rigid group: 6 DOF, flexible atom: 3 DOF
  EXPECT_EQ(layout.numDofs(), 9);
  EXPECT_TRUE(layout.frameworkAtomRigidComDof(0).has_value());
  EXPECT_FALSE(layout.frameworkAtomDof(0, MinimizationDofAxis::X).has_value());
  EXPECT_TRUE(layout.frameworkAtomDof(2, MinimizationDofAxis::X).has_value());
}

TEST(framework_intramolecular_reader, groups_rejected_for_rigid_type_and_cycle)
{
  ForceField forceField = makeTestForceField();
  Framework framework = makePeriodicPairFramework(forceField);

  TemporaryFile rigidGroups("raspa3-framework-rigid-groups.json",
                            R"({"Type": "Rigid",
                                "Groups": [{"Type": "Fixed", "Atoms": [0, 1]}]})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, rigidGroups.path().string()), std::runtime_error);

  Framework flexible = makePeriodicPairFramework(forceField);
  TemporaryFile cycle("raspa3-framework-cycle-group.json",
                      R"({"Type": "Flexible",
                          "Groups": [{"Type": "Cycle", "Atoms": [0, 1]}]})");
  EXPECT_THROW(flexible.readFrameworkDefinition(forceField, cycle.path().string()), std::runtime_error);
}

TEST(framework_intramolecular_reader, archive_round_trip_preserves_topology_and_potentials)
{
  ForceField forceField = makeTestForceField();
  Framework framework = makePeriodicPairFramework(forceField);
  TemporaryFile definition("raspa3-framework-archive-definition.json",
                           R"({"Type": "Flexible",
                               "Bonds": [[["C1", "C2"], "HARMONIC", [100.0, 1.0]]]})");
  framework.readFrameworkDefinition(forceField, definition.path().string());

  const std::filesystem::path archivePath =
      std::filesystem::temp_directory_path() / "raspa3-framework-intramolecular-reader.bin";
  {
    std::ofstream stream(archivePath, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << framework;
  }

  Framework restored;
  {
    std::ifstream stream(archivePath, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(archivePath);

  EXPECT_EQ(restored.connectivityTable.findAllBonds(), framework.connectivityTable.findAllBonds());
  EXPECT_FALSE(restored.rigid);
  ASSERT_EQ(restored.intraMolecularPotentials.bonds.size(), 1);
  EXPECT_EQ(restored.intraMolecularPotentials.bonds[0].identifiers,
            framework.intraMolecularPotentials.bonds[0].identifiers);
  EXPECT_EQ(restored.intraMolecularPotentials.bonds[0].type, framework.intraMolecularPotentials.bonds[0].type);
  EXPECT_EQ(restored.intraMolecularImageShifts.bonds, framework.intraMolecularImageShifts.bonds);
}
