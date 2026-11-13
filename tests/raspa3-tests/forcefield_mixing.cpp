#include <gtest/gtest.h>

import std;

import double4;
import units;
import pseudo_atom;
import vdwparameters;
import forcefield;
import potential_pair_vdw;

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
