#include <gtest/gtest.h>

import std;

import atom;
import forcefield;
import component;
import units;
import json;
import property_enthalpy;
import mc_moves_move_types;
import mc_moves_probabilities;

// The enthalpy of adsorption from the particle-number and energy fluctuations: a block whose covariance is
// singular has no enthalpy and is skipped rather than read off a failed inverse (the build uses fast math, so the
// data carries an explicit flag instead of a not-a-number); the total over the components
// is weighted by the adsorbed-phase mole fractions of the block and reported in Kelvin once.

namespace
{

constexpr double temperature = 300.0;

Component makeComponent(const ForceField& forceField, std::string name)
{
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Swap, 1.0);
  Component component(forceField, std::move(name), 190.564, 4599200, 0.01142,
                      {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities, std::nullopt,
                      false);
  component.molFraction = 0.5;  // the bulk gas: not the weight of the total
  component.idealGasEnergy = 0.0;
  return component;
}

// one sample: the counts of the two components and the total energy
void addSample(PropertyEnthalpy& property, std::size_t block, std::size_t nA, std::size_t nB, double energy)
{
  property.addSample(block, EnthalpyOfAdsorptionTerms({0, 1}, {nA, nB}, energy, temperature), 1.0);
}

// a block in which both counts fluctuate (the counts are correlated with the energy, so the enthalpies are
// finite and distinct per component)
void addFluctuatingBlock(PropertyEnthalpy& property, std::size_t block, std::size_t nA, std::size_t nB)
{
  addSample(property, block, nA, nB, -10.0 * static_cast<double>(nA) - 25.0 * static_cast<double>(nB));
  addSample(property, block, nA + 1, nB, -10.0 * static_cast<double>(nA + 1) - 25.0 * static_cast<double>(nB) + 1.0);
  addSample(property, block, nA, nB + 1, -10.0 * static_cast<double>(nA) - 25.0 * static_cast<double>(nB + 1) - 2.0);
  addSample(property, block, nA + 1, nB + 1,
            -10.0 * static_cast<double>(nA + 1) - 25.0 * static_cast<double>(nB + 1) + 0.5);
}

}  // namespace

TEST(ENTHALPY_OF_ADSORPTION, a_block_with_a_singular_covariance_has_no_enthalpy_and_is_skipped)
{
  // one component whose count never changed in the block
  EnthalpyOfAdsorptionTerms constant =
      EnthalpyOfAdsorptionTerms({0}, {3}, -1.0, temperature) + EnthalpyOfAdsorptionTerms({0}, {3}, -2.0, temperature);
  const EnthalpyOfAdsorptionData singular = (constant / 2.0).compositeProperty();
  EXPECT_FALSE(singular.defined);  // not -RT
  EXPECT_FALSE(isFinite(singular));

  // five blocks, the third one without fluctuations: the mean and the error come from the other four
  PropertyEnthalpy property(5, 2);
  for (std::size_t block = 0; block < 5; ++block)
  {
    if (block == 2)
    {
      addSample(property, block, 2, 3, -100.0);
      addSample(property, block, 2, 3, -101.0);
      continue;
    }
    addFluctuatingBlock(property, block, 2, 3);
  }
  const auto [mean, error] = property.average();
  EXPECT_FALSE(property.averaged(2).defined);
  EXPECT_TRUE(mean.defined);
  EXPECT_TRUE(error.defined);
  for (std::size_t k = 0; k < 2; ++k)
  {
    double sum = 0.0;
    for (const std::size_t block : {0uz, 1uz, 3uz, 4uz}) sum += property.averaged(block).values[k];
    EXPECT_NEAR(mean.values[k], sum / 4.0, 1e-9) << "component " << k;
  }
}

TEST(ENTHALPY_OF_ADSORPTION, the_total_is_weighted_by_the_adsorbed_fractions_and_converted_once)
{
  const ForceField forceField({{"X", false, 16.0, 0.0, 0.0, 1, false}}, {{150.0, 3.7}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  std::vector<Component> components{makeComponent(forceField, "a"), makeComponent(forceField, "b")};
  std::vector<std::size_t> swappable{0, 1};

  // component a adsorbs about one molecule, component b about three: weights 1/4 and 3/4, not the bulk halves
  PropertyEnthalpy property(5, 2);
  for (std::size_t block = 0; block < 5; ++block) addFluctuatingBlock(property, block, 1, 3);

  const nlohmann::json status = property.jsonAveragesStatistics(swappable, components);
  const std::vector<double> blockA = status["a"]["block"].get<std::vector<double>>();
  const std::vector<double> blockB = status["b"]["block"].get<std::vector<double>>();
  const std::vector<double> blockTotal = status["total"]["block"].get<std::vector<double>>();
  ASSERT_EQ(blockTotal.size(), 5uz);
  for (std::size_t block = 0; block < 5; ++block)
  {
    const std::vector<double> fractions = property.blockAdsorbedMoleFractions(block);
    EXPECT_NEAR(fractions[0], 1.5 / 5.0, 1e-12);
    EXPECT_NEAR(fractions[1], 3.5 / 5.0, 1e-12);
    // the per-component block values are in Kelvin; the total is their weighted sum, in Kelvin once
    EXPECT_NEAR(blockTotal[block], fractions[0] * blockA[block] + fractions[1] * blockB[block], 1e-9)
        << "block " << block;
    EXPECT_NE(blockTotal[block], 0.5 * blockA[block] + 0.5 * blockB[block]);
  }
  const double meanOfBlocks = std::accumulate(blockTotal.begin(), blockTotal.end(), 0.0) / 5.0;
  EXPECT_NEAR(status["total"]["mean"]["[K]"].get<double>(), meanOfBlocks, 1e-9);
  EXPECT_NEAR(status["total"]["mean"]["[kJ/mol]"].get<double>(),
              meanOfBlocks * Units::EnergyToKJPerMol / Units::EnergyToKelvin, 1e-9);

  // the text report agrees with the JSON
  const std::string text = property.writeAveragesStatistics(swappable, components);
  EXPECT_NE(text.find(std::format("{: .6e}", meanOfBlocks)), std::string::npos);
}
