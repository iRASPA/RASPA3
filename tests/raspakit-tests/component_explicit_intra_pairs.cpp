#include <gtest/gtest.h>

import std;

import atom;
import component;
import forcefield;
import interactions_internal;
import json;
import molecule;
import running_energy;

namespace
{
ForceField makeForceField()
{
  return ForceField({{"C", false, 12.0, 0.5, 0.0, 6, true}}, {{120.0, 3.4}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, false, true);
}

nlohmann::json chainDefinition()
{
  return nlohmann::json::parse(R"({
    "PseudoAtoms": [
      ["C", [0.0, 0.0, 0.0]],
      ["C", [1.5, 0.0, 0.0]],
      ["C", [3.0, 0.0, 0.0]],
      ["C", [4.5, 0.0, 0.0]],
      ["C", [6.0, 0.0, 0.0]]
    ],
    "Connectivity": [[0, 1], [1, 2], [2, 3], [3, 4]],
    "Bonds": [
      [["C", "C"], "FIXED", [1.5]]
    ],
    "Bends": [
      [["C", "C", "C"], "HARMONIC", [100.0, 180.0]]
    ],
    "Torsions": [
      [["C", "C", "C", "C"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]]
    ]
  })");
}

class TemporaryComponentFile
{
 public:
  explicit TemporaryComponentFile(const nlohmann::json& definition)
  {
    static std::size_t counter{};
    path = std::filesystem::temp_directory_path() /
           std::format("raspa3-explicit-intra-pairs-{}.json", counter++);
    std::ofstream stream(path);
    stream << definition.dump(2);
  }

  ~TemporaryComponentFile()
  {
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
  }

  std::string stem() const
  {
    std::filesystem::path result = path;
    result.replace_extension();
    return result.string();
  }

 private:
  std::filesystem::path path;
};

Component readComponent(const ForceField& forceField, const nlohmann::json& definition)
{
  TemporaryComponentFile file(definition);
  return Component(Component::Type::Adsorbate, 0, forceField, "explicit-intra-test", file.stem(), 5, 21);
}
}  // namespace

TEST(component_explicit_intra_pairs, absent_empty_and_explicit_pair_lists)
{
  const ForceField forceField = makeForceField();

  const Component automatic = readComponent(forceField, chainDefinition());
  ASSERT_EQ(automatic.intraMolecularPotentials.vanDerWaals.size(), 1uz);
  ASSERT_EQ(automatic.intraMolecularPotentials.coulombs.size(), 1uz);
  EXPECT_EQ(automatic.intraMolecularPotentials.vanDerWaals.front().identifiers,
            (std::array<std::size_t, 2>{0, 4}));
  EXPECT_EQ(automatic.intraMolecularPotentials.coulombs.front().identifiers,
            (std::array<std::size_t, 2>{0, 4}));

  nlohmann::json disabledDefinition = chainDefinition();
  disabledDefinition["IntraVanDerWaalsPairs"] = nlohmann::json::array();
  disabledDefinition["IntraCoulombPairs"] = nlohmann::json::array();
  const Component disabled = readComponent(forceField, disabledDefinition);
  EXPECT_TRUE(disabled.intraMolecularPotentials.vanDerWaals.empty());
  EXPECT_TRUE(disabled.intraMolecularPotentials.coulombs.empty());

  nlohmann::json explicitDefinition = chainDefinition();
  explicitDefinition["IntraVanDerWaalsPairs"] = {{3, 0}, {4, 1}};
  explicitDefinition["IntraCoulombPairs"] = {{3, 0}, {4, 1}};
  const Component explicitPairs = readComponent(forceField, explicitDefinition);
  ASSERT_EQ(explicitPairs.intraMolecularPotentials.vanDerWaals.size(), 2uz);
  ASSERT_EQ(explicitPairs.intraMolecularPotentials.coulombs.size(), 2uz);
  EXPECT_EQ(explicitPairs.intraMolecularPotentials.vanDerWaals[0].identifiers,
            (std::array<std::size_t, 2>{0, 3}));
  EXPECT_EQ(explicitPairs.intraMolecularPotentials.vanDerWaals[1].identifiers,
            (std::array<std::size_t, 2>{1, 4}));
  EXPECT_EQ(explicitPairs.intraMolecularPotentials.coulombs[0].identifiers,
            (std::array<std::size_t, 2>{0, 3}));
  EXPECT_EQ(explicitPairs.intraMolecularPotentials.coulombs[1].identifiers,
            (std::array<std::size_t, 2>{1, 4}));

  std::vector<Molecule> molecules{explicitPairs.createMoleculeRecord(explicitPairs.atoms)};
  molecules.front().atomIndex = 0;
  const RunningEnergy energy = Interactions::computeIntraMolecularEnergy(
      explicitPairs.intraMolecularPotentials, molecules, explicitPairs.atoms);
  const double expectedVDW =
      explicitPairs.intraMolecularPotentials.vanDerWaals[0].calculateEnergy(
          explicitPairs.atoms[0].position, explicitPairs.atoms[3].position) +
      explicitPairs.intraMolecularPotentials.vanDerWaals[1].calculateEnergy(
          explicitPairs.atoms[1].position, explicitPairs.atoms[4].position);
  const double expectedCoulomb =
      explicitPairs.intraMolecularPotentials.coulombs[0].calculateEnergy(
          explicitPairs.atoms[0].position, explicitPairs.atoms[3].position) +
      explicitPairs.intraMolecularPotentials.coulombs[1].calculateEnergy(
          explicitPairs.atoms[1].position, explicitPairs.atoms[4].position);
  EXPECT_DOUBLE_EQ(energy.intraVDW, expectedVDW);
  EXPECT_DOUBLE_EQ(energy.intraCoul, expectedCoulomb);
}

TEST(component_explicit_intra_pairs, rejects_invalid_pairs)
{
  const ForceField forceField = makeForceField();
  const std::array<nlohmann::json, 4> invalidLists{
      nlohmann::json::array({nlohmann::json::array({0})}),
      nlohmann::json::array({nlohmann::json::array({"0", 4})}),
      nlohmann::json::array({nlohmann::json::array({0, 5})}),
      nlohmann::json::array({nlohmann::json::array({0, 4}), nlohmann::json::array({4, 0})})};

  for (const nlohmann::json& invalidList : invalidLists)
  {
    nlohmann::json definition = chainDefinition();
    definition["IntraVanDerWaalsPairs"] = invalidList;
    EXPECT_THROW(readComponent(forceField, definition), std::runtime_error) << invalidList.dump();
  }
}

TEST(component_explicit_intra_pairs, rejects_ambiguous_14_scaling)
{
  const ForceField forceField = makeForceField();

  nlohmann::json vdwDefinition = chainDefinition();
  vdwDefinition["IntraVanDerWaalsPairs"] = nlohmann::json::array();
  vdwDefinition["Intra14VanDerWaalsScalingValue"] = 0.5;
  EXPECT_THROW(readComponent(forceField, vdwDefinition), std::runtime_error);

  nlohmann::json coulombDefinition = chainDefinition();
  coulombDefinition["IntraCoulombPairs"] = nlohmann::json::array();
  coulombDefinition["Intra14ChargeChargeScalingValue"] = 0.5;
  EXPECT_THROW(readComponent(forceField, coulombDefinition), std::runtime_error);
}
