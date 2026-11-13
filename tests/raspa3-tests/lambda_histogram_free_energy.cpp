#include <gtest/gtest.h>

import std;

import units;
import property_lambda_probability_histogram;

// The Landau free-energy profile of a lambda histogram sampled under the Wang-Landau bias: the bias has to
// be put back, or a flattened histogram reports a flat profile. The profile is reported relative to its
// minimum in whichever unit system is active.

namespace
{

constexpr std::size_t bins = 5;
constexpr double beta = 1.0 / 300.0;  // in internal units (1/K)

// a histogram flattened by Wang-Landau: equal counts in every bin, the true free energy sits in the bias
PropertyLambdaProbabilityHistogram makeFlattenedHistogram(const std::vector<double>& freeEnergyInKelvin)
{
  PropertyLambdaProbabilityHistogram histogram(1, bins);
  for (std::size_t i = 0; i < bins; ++i)
  {
    histogram.biasFactor[i] = beta * freeEnergyInKelvin[i];
    histogram.bookKeepingLambda[0][i] = 1000.0;
  }
  return histogram;
}

class ScopedUnitSystem
{
 public:
  explicit ScopedUnitSystem(Units::System system) : previous_(Units::unitSystem) { Units::unitSystem = system; }
  ~ScopedUnitSystem() { Units::unitSystem = previous_; }

 private:
  Units::System previous_;
};

// the free-energy values of the report, one per bin
std::vector<double> reportedProfile(const std::string& report)
{
  std::vector<double> values;
  std::istringstream lines(report);
  std::string line;
  while (std::getline(lines, line))
  {
    const std::size_t at = line.find("Free energy: ");
    if (at == std::string::npos) continue;
    values.push_back(std::stod(line.substr(at + std::string("Free energy: ").size())));
  }
  return values;
}

}  // namespace

TEST(LAMBDA_HISTOGRAM_FREE_ENERGY, the_bias_is_put_back_into_the_landau_profile)
{
  const std::vector<double> freeEnergy{0.0, 400.0, 900.0, 300.0, -200.0};  // [K]
  const PropertyLambdaProbabilityHistogram histogram = makeFlattenedHistogram(freeEnergy);

  const std::vector<double> profile = histogram.averagedLandauFreeEnergyHistogram(beta);
  ASSERT_EQ(profile.size(), bins);
  for (std::size_t i = 1; i < bins; ++i)
  {
    EXPECT_NEAR(profile[i] - profile[0], freeEnergy[i] - freeEnergy[0], 1e-9) << "bin " << i;
  }

  // the excess chemical potential (with the bias difference the reports add) is the end-to-end profile step
  const double excess = histogram.averagedExcessChemicalPotential(beta) +
                        (histogram.biasFactor[bins - 1] - histogram.biasFactor[0]) / beta;
  EXPECT_NEAR(excess, profile[bins - 1] - profile[0], 1e-9);
}

TEST(LAMBDA_HISTOGRAM_FREE_ENERGY, the_report_is_relative_to_the_minimum_in_both_unit_systems)
{
  const std::vector<double> freeEnergy{0.0, 400.0, 900.0, 300.0, -200.0};  // [K]
  const PropertyLambdaProbabilityHistogram histogram = makeFlattenedHistogram(freeEnergy);

  {
    ScopedUnitSystem units(Units::System::RASPA);
    const std::vector<double> profile =
        reportedProfile(histogram.writeAveragesStatistics(beta, std::nullopt, std::nullopt));
    ASSERT_EQ(profile.size(), bins);
    EXPECT_NEAR(*std::ranges::min_element(profile), 0.0, 1e-6);
    EXPECT_NEAR(profile[2], Units::EnergyToKelvin * (900.0 + 200.0), 1e-3);
  }
  {
    ScopedUnitSystem units(Units::System::ReducedUnits);
    const std::vector<double> profile =
        reportedProfile(histogram.writeAveragesStatistics(beta, std::nullopt, std::nullopt));
    ASSERT_EQ(profile.size(), bins);
    EXPECT_NEAR(*std::ranges::min_element(profile), 0.0, 1e-6);
    EXPECT_NEAR(profile[2], beta * (900.0 + 200.0), 1e-6);
  }
}

// Blocks with different numbers of samples describe the same probability distribution: each block is normalised by
// its own count, so two flat blocks of different size agree bin by bin and leave no error on the profile.
TEST(LAMBDA_HISTOGRAM_FREE_ENERGY, blocks_of_different_size_share_the_profile)
{
  PropertyLambdaProbabilityHistogram histogram(3, bins);
  for (std::size_t i = 0; i < bins; ++i)
  {
    histogram.bookKeepingLambda[0][i] = 1000.0;
    histogram.bookKeepingLambda[1][i] = 10.0;
    histogram.bookKeepingLambda[2][i] = 100.0;
  }
  const auto [profile, error] = histogram.averageLandauFreeEnergyHistogram(beta);
  ASSERT_EQ(error.size(), bins);
  for (std::size_t i = 0; i < bins; ++i)
  {
    EXPECT_NEAR(error[i], 0.0, 1e-9 * std::max(1.0, std::abs(profile[i]))) << "bin " << i;
  }
  const auto [probability, probabilityError] = histogram.result();
  for (std::size_t i = 0; i < bins; ++i)
  {
    EXPECT_NEAR(probability[i], 1.0, 1e-12) << "bin " << i;
    EXPECT_NEAR(probabilityError[i], 0.0, 1e-12) << "bin " << i;
  }
}

// A histogram without samples, or with an unsampled bin, reports finite numbers (no 0/0, no log of zero).
TEST(LAMBDA_HISTOGRAM_FREE_ENERGY, empty_histograms_and_unsampled_bins_stay_finite)
{
  PropertyLambdaProbabilityHistogram empty(2, bins);
  const auto [probability, probabilityError] = empty.result();
  for (std::size_t i = 0; i < bins; ++i)
  {
    EXPECT_EQ(std::to_string(probability[i]).find("nan"), std::string::npos) << "bin " << i;
    EXPECT_EQ(std::to_string(probabilityError[i]).find("nan"), std::string::npos) << "bin " << i;
  }

  PropertyLambdaProbabilityHistogram partial(2, bins);
  for (std::size_t block = 0; block < 2; ++block)
  {
    for (std::size_t i = 0; i + 1 < bins; ++i) partial.bookKeepingLambda[block][i] = 100.0;
  }
  const auto [profile, profileError] = partial.averageLandauFreeEnergyHistogram(beta);
  for (std::size_t i = 0; i < bins; ++i)
  {
    EXPECT_EQ(std::to_string(profile[i]).find("nan"), std::string::npos) << "bin " << i;
    EXPECT_EQ(std::to_string(profile[i]).find("inf"), std::string::npos) << "bin " << i;
    EXPECT_EQ(std::to_string(profileError[i]).find("nan"), std::string::npos) << "bin " << i;
  }
  const std::string report = partial.writeAveragesStatistics(beta, std::nullopt, std::nullopt);
  EXPECT_EQ(report.find("nan"), std::string::npos) << report;
  EXPECT_EQ(report.find("inf"), std::string::npos) << report;
}

// A bias file without values is an input error, not an out-of-range read.
TEST(LAMBDA_HISTOGRAM_FREE_ENERGY, an_empty_bias_file_is_rejected)
{
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / std::format("raspa3_empty_bias_{}.json", ::getpid());
  {
    std::ofstream stream(path);
    stream << R"({"bias": []})";
  }
  PropertyLambdaProbabilityHistogram histogram(1, bins);
  EXPECT_THROW(histogram.readBiasingFile(path), std::runtime_error);
  std::filesystem::remove(path);
}
