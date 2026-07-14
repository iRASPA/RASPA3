#include <gtest/gtest.h>

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
import vdwparameters;

namespace
{
std::filesystem::path repositoryRoot()
{
  return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
}

std::string readLocalFile(const std::filesystem::path& path)
{
  std::ifstream stream(path);
  if (!stream) throw std::runtime_error(std::format("Could not read '{}'", path.string()));
  return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
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

ForceField makeTestForceField()
{
  return ForceField({{"C1", true, 12.0, 0.0, 0.0, 6, true},
                     {"C2", true, 12.0, 0.0, 0.0, 6, true},
                     {"C3", true, 12.0, 0.0, 0.0, 6, true},
                     {"H1", true, 1.0, 0.0, 0.0, 1, true}},
                    {{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0,
                    12.0, 12.0, false, false, false);
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

  framework.readFrameworkDefinition(forceField, definition.path.string());

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
  framework.readFrameworkDefinition(forceField, definition.path.string());
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
      framework, SimulationBox(10.0, 10.0, 10.0), framework.atoms, layout, hessian, dynamics);
  EXPECT_DOUBLE_EQ(energy.potentialEnergy(), 0.0);
  EXPECT_TRUE(std::ranges::all_of(dynamics, [](const AtomDynamics& value) {
    return value.gradient == double3(0.0, 0.0, 0.0);
  }));
}

TEST(framework_intramolecular_reader, periodic_bond_gradient_and_hessian_match_finite_difference)
{
  ForceField forceField = makeTestForceField();
  Framework framework = makePeriodicPairFramework(forceField);
  TemporaryFile definition("raspa3-framework-periodic-derivatives.json",
                           R"({"Type": "Flexible",
                               "Bonds": [[["C1", "C2"], "HARMONIC", [100.0, 1.0]]]})");
  framework.readFrameworkDefinition(forceField, definition.path.string());

  const SimulationBox box(10.0, 10.0, 10.0);
  std::vector<Atom> atoms = framework.atoms;
  std::vector<AtomDynamics> dynamics(atoms.size());
  const MinimizationDofLayout layout = buildMinimizationDofLayout({}, {}, atoms.size());
  GeneralizedHessian hessian(layout.numDofs(), 0);
  const RunningEnergy analytic =
      Interactions::computeFrameworkIntraMolecularHessian(framework, box, atoms, layout, hessian, dynamics);
  EXPECT_DOUBLE_EQ(analytic.bond, Interactions::computeFrameworkIntraMolecularEnergy(framework, box, atoms).bond);

  const auto energy = [&]() { return Interactions::computeFrameworkIntraMolecularEnergy(framework, box, atoms).bond; };
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
  framework.readFrameworkDefinition(forceField, definition.path.string());
  ASSERT_EQ(framework.intraMolecularImageShifts.bends.size(), 2);
  ASSERT_EQ(framework.intraMolecularImageShifts.torsions.size(), 1);
  EXPECT_EQ(framework.intraMolecularImageShifts.torsions[0][2], int3(1, 0, 0));
  EXPECT_EQ(framework.intraMolecularImageShifts.torsions[0][3], int3(1, 0, 0));

  std::vector<Atom> atoms = framework.atoms;
  const MinimizationDofLayout layout = buildMinimizationDofLayout({}, {}, atoms.size());
  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<AtomDynamics> referenceDynamics(atoms.size());
  Interactions::computeFrameworkIntraMolecularHessian(framework, box, atoms, layout, hessian, referenceDynamics);

  constexpr double delta = 2.0e-5;
  constexpr double tolerance = 5.0e-2;
  for (std::size_t column = 0; column < layout.numDofs(); ++column)
  {
    const std::size_t atom = column / 3;
    const std::size_t axis = column % 3;
    (&atoms[atom].position.x)[axis] += delta;
    GeneralizedHessian plusHessian(layout.numDofs(), 0);
    std::vector<AtomDynamics> plusDynamics(atoms.size());
    Interactions::computeFrameworkIntraMolecularHessian(framework, box, atoms, layout, plusHessian, plusDynamics);

    (&atoms[atom].position.x)[axis] -= 2.0 * delta;
    GeneralizedHessian minusHessian(layout.numDofs(), 0);
    std::vector<AtomDynamics> minusDynamics(atoms.size());
    Interactions::computeFrameworkIntraMolecularHessian(framework, box, atoms, layout, minusHessian, minusDynamics);
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
  const std::filesystem::path root = repositoryRoot();
  const std::filesystem::path forceFieldDirectory =
      root / "examples/non_basic/1_mc_adsorption_binary_mixture_co2_ch4_in_irmof_1";
  ForceField forceField = ForceField::readForceField(forceFieldDirectory.string(), "force_field.json").value();

  const std::filesystem::path cifPath = root / "examples/basic/19_minimization_co2_in_irmof_1/IRMOF-1.cif";
  const auto cif = CIFReader::readCIFString(readLocalFile(cifPath), forceField, CIFReader::UseChargesFrom::PseudoAtoms);
  ASSERT_TRUE(cif.has_value());
  auto [box, hallNumber, definedAtoms, fractionalAtoms] = *cif;
  Framework framework(forceField, "IRMOF-1", box, hallNumber, definedAtoms, fractionalAtoms, {1, 1, 1});

  framework.readFrameworkDefinition(forceField,
                                    (root / "data/frameworks/Dubbeldam2007FlexibleIRMOF-1/framework.json").string());

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
  std::set<std::array<std::size_t, 2>> pairs12And13;
  for (const std::array<std::size_t, 2>& bond : framework.connectivityTable.findAllBonds())
  {
    pairs12And13.insert({std::min(bond[0], bond[1]), std::max(bond[0], bond[1])});
  }
  for (const std::array<std::size_t, 3>& bend : framework.connectivityTable.findAllBends())
  {
    pairs12And13.insert({std::min(bend[0], bend[2]), std::max(bend[0], bend[2])});
  }
  std::set<std::array<std::size_t, 2>> pairs14;
  for (const std::array<std::size_t, 4>& torsion : framework.connectivityTable.findAllTorsions())
  {
    const std::array<std::size_t, 2> pair{std::min(torsion[0], torsion[3]), std::max(torsion[0], torsion[3])};
    if (!pairs12And13.contains(pair)) pairs14.insert(pair);
  }
  const std::size_t expectedNonBondedCount = framework.connectivityTable.findAllVanDerWaals().size() + pairs14.size();
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
  const std::filesystem::path root = repositoryRoot();
  ForceField forceField =
      ForceField::readForceField((root / "examples/non_basic/11_framework_intramolecular_reader_irmof_1").string(),
                                 "force_field.json")
          .value();

  const std::size_t carbon = *forceField.findPseudoAtom("C3");
  const std::size_t hydrogen = *forceField.findPseudoAtom("H1");
  std::vector<Atom> atoms{fractionalAtom({0.100, 0.5, 0.5}, carbon), fractionalAtom({0.165, 0.5, 0.5}, hydrogen),
                          fractionalAtom({0.230, 0.5, 0.5}, carbon), fractionalAtom({0.295, 0.5, 0.5}, hydrogen)};
  Framework framework(forceField, "scaled-chain", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms, {1, 1, 1});
  TemporaryFile definition("raspa3-framework-14-scaling.json",
                           R"({"Intra14VanDerWaalsScalingValue": 0.5,
          "Intra14ChargeChargeScalingValue": 0.75,
          "Bonds": [[["C3", "H1"], "HARMONIC", [100.0, 1.3]]]})");

  framework.readFrameworkDefinition(forceField, definition.path.string());

  ASSERT_EQ(framework.connectivityTable.findAllTorsions().size(), 1);
  ASSERT_EQ(framework.intraMolecularPotentials.vanDerWaals.size(), 1);
  ASSERT_EQ(framework.intraMolecularPotentials.coulombs.size(), 1);
  EXPECT_DOUBLE_EQ(framework.intraMolecularPotentials.vanDerWaals[0].scaling, 0.5);
  EXPECT_DOUBLE_EQ(framework.intraMolecularPotentials.coulombs[0].scaling, 0.75);
}

TEST(framework_intramolecular_reader, controls_12_and_13_exclusions_independently)
{
  const std::filesystem::path root = repositoryRoot();
  ForceField forceField =
      ForceField::readForceField((root / "examples/non_basic/11_framework_intramolecular_reader_irmof_1").string(),
                                 "force_field.json")
          .value();

  const std::size_t carbon = *forceField.findPseudoAtom("C3");
  const std::size_t hydrogen = *forceField.findPseudoAtom("H1");
  std::vector<Atom> atoms{fractionalAtom({0.100, 0.5, 0.5}, carbon), fractionalAtom({0.165, 0.5, 0.5}, hydrogen),
                          fractionalAtom({0.230, 0.5, 0.5}, carbon), fractionalAtom({0.295, 0.5, 0.5}, hydrogen)};

  Framework include12(forceField, "include-12", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms, {1, 1, 1});
  TemporaryFile include12Definition("raspa3-framework-include-12.json",
                                    R"({"ExcludeIntra12Interactions": false,
          "Bonds": [[["C3", "H1"], "HARMONIC", [100.0, 1.3]]]})");
  include12.readFrameworkDefinition(forceField, include12Definition.path.string());
  EXPECT_EQ(include12.intraMolecularPotentials.vanDerWaals.size(), 3);
  EXPECT_EQ(include12.intraMolecularPotentials.coulombs.size(), 3);

  Framework include13(forceField, "include-13", SimulationBox(20.0, 20.0, 20.0), 1, atoms, atoms, {1, 1, 1});
  TemporaryFile include13Definition("raspa3-framework-include-13.json",
                                    R"({"ExcludeIntra13Interactions": false,
          "Bonds": [[["C3", "H1"], "HARMONIC", [100.0, 1.3]]]})");
  include13.readFrameworkDefinition(forceField, include13Definition.path.string());
  EXPECT_EQ(include13.intraMolecularPotentials.vanDerWaals.size(), 2);
  EXPECT_EQ(include13.intraMolecularPotentials.coulombs.size(), 2);
}

TEST(framework_intramolecular_reader, rejects_unknown_ambiguous_and_malformed_definitions)
{
  ForceField forceField = makeTestForceField();
  Framework framework = makePeriodicPairFramework(forceField);

  TemporaryFile unknown("raspa3-framework-unknown-definition.json",
                        R"({"Bonds": [[["missing", "C1"], "HARMONIC", [100.0, 1.0]]]})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, unknown.path.string()), std::runtime_error);

  TemporaryFile ambiguous("raspa3-framework-ambiguous-definition.json",
                          R"({"Bonds": [[["C1", "C2"], "HARMONIC", [100.0, 1.0]],
                    [["C2", "C1"], "HARMONIC", [100.0, 1.0]]]})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, ambiguous.path.string()), std::runtime_error);

  TemporaryFile malformed("raspa3-framework-malformed-definition.json",
                          R"({"Bonds": [[["C1", "C2"], "HARMONIC", [100.0]]]})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, malformed.path.string()), std::runtime_error);

  TemporaryFile invalidScaling("raspa3-framework-invalid-scaling.json", R"({"Intra14VanDerWaalsScalingValue": -0.5})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, invalidScaling.path.string()), std::runtime_error);

  TemporaryFile invalidExclusion("raspa3-framework-invalid-exclusion.json", R"({"ExcludeIntra12Interactions": 1})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, invalidExclusion.path.string()), std::runtime_error);

  TemporaryFile invalidType("raspa3-framework-invalid-type.json", R"({"Type": "deformable"})");
  EXPECT_THROW(framework.readFrameworkDefinition(forceField, invalidType.path.string()), std::runtime_error);
}

TEST(framework_intramolecular_reader, archive_round_trip_preserves_topology_and_potentials)
{
  ForceField forceField = makeTestForceField();
  Framework framework = makePeriodicPairFramework(forceField);
  TemporaryFile definition("raspa3-framework-archive-definition.json",
                           R"({"Type": "Flexible",
                               "Bonds": [[["C1", "C2"], "HARMONIC", [100.0, 1.0]]]})");
  framework.readFrameworkDefinition(forceField, definition.path.string());

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
