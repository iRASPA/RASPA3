#include <gtest/gtest.h>

import std;

import atom;
import component;
import forcefield;
import json;

namespace
{
// one Lennard-Jones type 'C' with a default charge of 0.5
ForceField makeForceField()
{
  return ForceField({{"C", false, 12.0, 0.5, 0.0, 6, true}}, {{120.0, 3.4}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, false, true);
}

class TemporaryComponentFile
{
 public:
  explicit TemporaryComponentFile(const nlohmann::json& definition)
  {
    static std::size_t counter{};
    path = std::filesystem::temp_directory_path() / std::format("raspa3-per-atom-charges-{}.json", counter++);
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
  return Component(Component::Type::Adsorbate, 0, forceField, "per-atom-charge-test", file.stem(), 5, 21);
}
}  // namespace

TEST(component_per_atom_charges, pseudo_atom_charge_is_the_default)
{
  const ForceField forceField = makeForceField();
  const Component component = readComponent(forceField, nlohmann::json::parse(R"({
    "PseudoAtoms": [
      ["C", [0.0, 0.0, 0.0]],
      ["C", [1.5, 0.0, 0.0]]
    ]
  })"));

  ASSERT_EQ(component.atoms.size(), 2uz);
  EXPECT_DOUBLE_EQ(component.atoms[0].charge, 0.5);
  EXPECT_DOUBLE_EQ(component.atoms[1].charge, 0.5);
  EXPECT_DOUBLE_EQ(component.netCharge, 1.0);
}

TEST(component_per_atom_charges, third_element_overrides_the_charge_per_atom)
{
  const ForceField forceField = makeForceField();
  const Component component = readComponent(forceField, nlohmann::json::parse(R"({
    "PseudoAtoms": [
      ["C", [0.0, 0.0, 0.0], -0.1825],
      ["C", [1.5, 0.0, 0.0]],
      ["C", [3.0, 0.0, 0.0], 0.0337]
    ],
    "Connectivity": [[0, 1], [1, 2]],
    "Bonds": [[["C", "C"], "FIXED", [1.5]]],
    "Bends": [[["C", "C", "C"], "HARMONIC", [100.0, 180.0]]]
  })"));

  ASSERT_EQ(component.atoms.size(), 3uz);
  EXPECT_DOUBLE_EQ(component.atoms[0].charge, -0.1825);
  EXPECT_DOUBLE_EQ(component.atoms[1].charge, 0.5);
  EXPECT_DOUBLE_EQ(component.atoms[2].charge, 0.0337);
  EXPECT_NEAR(component.netCharge, -0.1825 + 0.5 + 0.0337, 1e-12);
  // all atoms share the pseudo-atom type (mass and Lennard-Jones parameters)
  EXPECT_EQ(component.atoms[0].type, component.atoms[2].type);
  EXPECT_DOUBLE_EQ(component.totalMass, 36.0);
}

TEST(component_per_atom_charges, rejects_non_numeric_charge_and_wrong_arity)
{
  const ForceField forceField = makeForceField();

  EXPECT_THROW(readComponent(forceField, nlohmann::json::parse(R"({
    "PseudoAtoms": [["C", [0.0, 0.0, 0.0], "q"]]
  })")),
               std::runtime_error);

  EXPECT_THROW(readComponent(forceField, nlohmann::json::parse(R"({
    "PseudoAtoms": [["C", [0.0, 0.0, 0.0], 0.1, 0.2]]
  })")),
               std::runtime_error);
}
