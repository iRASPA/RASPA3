#include <gtest/gtest.h>

import std;

import double4;
import units;
import pseudo_atom;
import vdwparameters;
import forcefield;
import potential_pair_vdw;
import int3;
import atom;
import simulationbox;
import framework;
import component;
import system;

namespace
{

// A force field of one pseudo-atom per self interaction, mixed with the given rule.
ForceField makeForceField(const std::vector<VDWParameters>& selfInteractions,
                          ForceField::MixingRule mixingRule = ForceField::MixingRule::Lorentz_Berthelot)
{
  std::vector<PseudoAtom> pseudoAtoms;
  for (std::size_t i = 0; i != selfInteractions.size(); ++i)
  {
    pseudoAtoms.emplace_back(std::format("X{}", i), false, 12.0, 0.0, 0.0, 6, false);
  }
  return ForceField(pseudoAtoms, selfInteractions, mixingRule, 12.0, 12.0, 12.0, false, false, false);
}

// the scaled (CFCMC, thermodynamic integration, scaled 1-4) pair evaluation of every derivative order is finite;
// NaN is detected through its printed form, as the project flags fold std::isnan away
void expectFiniteScaledPair(const ForceField& forceField, std::size_t typeA, std::size_t typeB)
{
  for (double scaling : {0.0, 0.3, 0.999})
  {
    const auto pair = Potentials::potentialVDW<2>(forceField, scaling, 1.0, 9.0, typeA, typeB);
    for (double value : {pair.energy, pair.dUdlambda, pair.firstDerivativeFactor, pair.secondDerivativeFactor})
    {
      EXPECT_EQ(std::to_string(value).find("nan"), std::string::npos)
          << "types " << typeA << "-" << typeB << " at scaling " << scaling;
      EXPECT_EQ(std::to_string(value).find("inf"), std::string::npos)
          << "types " << typeA << "-" << typeB << " at scaling " << scaling;
    }
  }
}

}  // namespace

// A cross pair without a mixing rule (two different forms, or forms such as Buckingham and Morse) has no
// interaction, not a Lennard-Jones pair with sigma = 0 whose soft-core force is 0 * inf under scaling.
TEST(FORCEFIELD_MIXING, unmixed_cross_pairs_have_no_interaction_and_finite_scaled_forces)
{
  const ForceField forceField = makeForceField({
      VDWParameters(VDWParameters::Type::BuckingHam, {1.0e5, 3.5, 1.0e3}),
      VDWParameters(VDWParameters::Type::Morse, {100.0, 1.5, 3.0}),
      VDWParameters(VDWParameters::Type::Mie, {1.0e6, 12.0, 1.0e3, 6.0}),
      VDWParameters(100.0, 3.5),
  });

  for (std::size_t i = 0; i != 4; ++i)
  {
    for (std::size_t j = 0; j != 4; ++j)
    {
      if (i == j) continue;
      EXPECT_EQ(forceField(i, j).type, VDWParameters::Type::None) << i << "-" << j;
      expectFiniteScaledPair(forceField, i, j);
    }
  }
}

// The epsilon-sigma forms other than plain Lennard-Jones mix like it when both pseudo-atoms use the same form.
TEST(FORCEFIELD_MIXING, epsilon_sigma_forms_mix_with_the_mixing_rule)
{
  const std::array types{VDWParameters::Type::LennardJonesShiftedForce, VDWParameters::Type::WeeksChandlerAndersen,
                         VDWParameters::Type::MM3, VDWParameters::Type::FeynmannHibbs,
                         VDWParameters::Type::CFFEpsilonSigma};
  for (VDWParameters::Type type : types)
  {
    const ForceField forceField =
        makeForceField({VDWParameters(type, {100.0, 3.0, 1.0}), VDWParameters(type, {400.0, 4.0, 4.0})});
    const VDWParameters& cross = forceField(0, 1);
    SCOPED_TRACE(VDWParameters::nameOfType(type));
    EXPECT_EQ(cross.type, type);
    EXPECT_EQ(forceField(1, 0).type, type);
    EXPECT_NEAR(cross.parameters.x, 200.0 * Units::KelvinToEnergy, 1e-12 * Units::KelvinToEnergy);
    EXPECT_DOUBLE_EQ(cross.parameters.y, 3.5);
    if (type == VDWParameters::Type::FeynmannHibbs)
    {
      // mu_ij = m_i m_j / (m_i + m_j) with m = 2 mu_ii: 2 * 1 * 4 / (1 + 4)
      EXPECT_DOUBLE_EQ(cross.parameters.z, 1.6);
    }
    expectFiniteScaledPair(forceField, 0, 1);
  }

  const ForceField jorgensen =
      makeForceField({VDWParameters(VDWParameters::Type::WeeksChandlerAndersen, {100.0, 3.0}),
                      VDWParameters(VDWParameters::Type::WeeksChandlerAndersen, {400.0, 12.0})},
                     ForceField::MixingRule::Jorgensen);
  EXPECT_EQ(jorgensen(0, 1).type, VDWParameters::Type::WeeksChandlerAndersen);
  EXPECT_DOUBLE_EQ(jorgensen(0, 1).parameters.y, 6.0);
}

// A pseudo-atom without self interaction (the pair table entries nobody set) has no interaction with anything.
TEST(FORCEFIELD_MIXING, pseudo_atom_without_self_interaction_has_no_interaction)
{
  std::vector<PseudoAtom> pseudoAtoms{{"A", false, 12.0, 0.0, 0.0, 6, false}, {"B", false, 12.0, 0.0, 0.0, 6, false}};
  const ForceField forceField(pseudoAtoms, {VDWParameters(100.0, 3.5)}, ForceField::MixingRule::Lorentz_Berthelot, 12.0,
                              12.0, 12.0, false, false, false);
  EXPECT_EQ(forceField(0, 0).type, VDWParameters::Type::LennardJones);
  EXPECT_EQ(forceField(0, 1).type, VDWParameters::Type::None);
  EXPECT_EQ(forceField(1, 1).type, VDWParameters::Type::None);
  expectFiniteScaledPair(forceField, 0, 1);
  expectFiniteScaledPair(forceField, 1, 1);
}

// The Feynman-Hibbs prefactor is hbar^2 / (24 u k_B) = 2.0212 K A^2 per (amu K): hbar^2 is h^2 / (2 pi)^2, and
// dividing h^2 by 2 pi alone made it 2 pi times too large.
TEST(FORCEFIELD_MIXING, feynman_hibbs_prefactor_is_hbar_squared_over_24_u_kb)
{
  const double hbar = Units::PlanckConstant / (2.0 * std::numbers::pi);
  const double expected =
      hbar * hbar / (24.0 * Units::AtomicMassUnit * Units::BoltzmannConstant * Units::Angstrom * Units::Angstrom);
  EXPECT_NEAR(Units::FeymannHibbsConversionFactor, expected, 1e-12 * expected);
  EXPECT_NEAR(Units::FeymannHibbsConversionFactor, 2.0212, 1e-4);
}

// The shift of a pair with a framework atom is that at the framework cut-off, where the framework-molecule energy is
// truncated, also for a framework type that the force field did not flag 'framework' (the default): the pair constants
// were evaluated at the molecule cut-off, and the shifted potential jumped at the framework cut-off.
TEST(FORCEFIELD_MIXING, framework_pairs_take_the_framework_cut_off_whatever_the_framework_flag)
{
  ForceField forceField({{"Si", false, 28.0855, 0.0, 0.0, 14, false}, {"CH4", false, 16.04246, 0.0, 0.0, 6, false}},
                        {VDWParameters(22.0, 2.30), VDWParameters(158.5, 3.72)},
                        ForceField::MixingRule::Lorentz_Berthelot, 12.0, 8.0, 12.0, true, false, false);
  const Framework framework(forceField, "host", SimulationBox(25.0, 25.0, 25.0), 1,
                            {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, int3(1, 1, 1));
  const Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                          {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false)}, {}, {}, 5, 21);
  System system(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {framework}, {methane}, {}, {0}, 5);

  ForceField reference = forceField;
  reference.pseudoAtoms[0].framework = true;
  reference.preComputeCutOffDependentParameters();
  EXPECT_TRUE(system.forceField.pseudoAtoms[0].framework);
  EXPECT_FALSE(system.forceField.pseudoAtoms[1].framework);
  EXPECT_DOUBLE_EQ(system.forceField(0, 1).shift, reference(0, 1).shift);
  EXPECT_DOUBLE_EQ(system.forceField(1, 1).shift, forceField(1, 1).shift);
  EXPECT_NE(system.forceField(0, 1).shift, forceField(0, 1).shift);

  // one type in both the framework and a molecule can not have its pairs at both cut-offs
  const Component silicon(forceField, "silicon", 190.564, 4599200, 0.01142,
                          {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);
  EXPECT_THROW(System(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {framework}, {silicon}, {}, {0}, 5),
               std::runtime_error);
}
