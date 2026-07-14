module;

export module property_elastic_constants_fluctuation;

import std;
import archive;
import json;
import units;
import elastic_tensor;
export import property_block_average;

export struct ElasticFluctuationData
{
  std::array<double, 6> configurationalStress{};
  std::array<double, 6> kineticStress{};
  std::array<double, 36> born{};
  std::array<double, 36> kinetic{};
  std::array<double, 36> fluctuation{};
  std::array<double, 36> stiffness{};
  std::array<double, 36> pressureCorrection{};
  std::array<double, 36> tangentStiffness{};
};

export ElasticFluctuationData operator-(const ElasticFluctuationData& a, const ElasticFluctuationData& b);
export ElasticFluctuationData operator*(const ElasticFluctuationData& a, const ElasticFluctuationData& b);
export ElasticFluctuationData operator*(double a, const ElasticFluctuationData& b);
export ElasticFluctuationData sqrt(const ElasticFluctuationData& a);
export ElasticFluctuationData& operator+=(ElasticFluctuationData& a, const ElasticFluctuationData& b);

export struct ElasticFluctuationTerms
{
  std::uint64_t versionNumber{1};
  std::array<double, 6> configurationalStress{};
  std::array<double, 6> kineticStress{};
  std::array<double, 36> stressProducts{};
  std::array<double, 36> born{};
  double beta{};
  double volume{};
  double kineticEntities{};

  ElasticFluctuationTerms() = default;
  ElasticFluctuationTerms(const std::array<double, 6>& configurationalStress,
                          const std::array<double, 6>& kineticStress, const std::array<double, 36>& born, double beta,
                          double volume, double kineticEntities);

  ElasticFluctuationTerms& operator+=(const ElasticFluctuationTerms& other);
  ElasticFluctuationData compositeProperty() const;
  ElasticFluctuationData zeroCompositeProperty() const { return {}; }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ElasticFluctuationTerms& terms);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ElasticFluctuationTerms& terms);
};

export ElasticFluctuationTerms operator*(double a, const ElasticFluctuationTerms& b);
export ElasticFluctuationTerms operator/(const ElasticFluctuationTerms& a, double b);

export struct PropertyElasticConstantsFluctuation : BlockAverage<ElasticFluctuationTerms>
{
  PropertyElasticConstantsFluctuation() = default;
  explicit PropertyElasticConstantsFluctuation(std::size_t numberOfBlocks)
      : BlockAverage<ElasticFluctuationTerms>(numberOfBlocks)
  {
  }

  std::string writeAveragesStatistics(double eigenvalueTolerance = 1.0e-8) const;
  nlohmann::json jsonAveragesStatistics(double eigenvalueTolerance = 1.0e-8) const;
};
