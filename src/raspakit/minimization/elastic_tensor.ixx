module;

export module elastic_tensor;

import std;

export constexpr std::size_t elasticVoigtSize = 6;
export constexpr std::array<std::size_t, elasticVoigtSize> elasticRegularToVoigt = {0, 3, 5, 4, 2, 1};

export struct ElasticDerivedProperties
{
  std::array<double, 36> compliance{};
  std::array<double, 6> stabilityEigenvalues{};
  std::array<double, 3> youngModuli{};
  std::array<double, 9> poissonRatios{};
  double bulkModulusVoigt{};
  double shearModulusVoigt{};
  double bulkModulusReuss{};
  double shearModulusReuss{};
  double bulkModulusHill{};
  double shearModulusHill{};
  bool complianceAvailable{};
};

export double& elasticAt(std::array<double, 36>& matrix, std::size_t row, std::size_t column);
export double elasticAt(const std::array<double, 36>& matrix, std::size_t row, std::size_t column);
export void symmetrizeElasticTensor(std::span<double> matrix, std::size_t size = elasticVoigtSize);
export ElasticDerivedProperties deriveElasticProperties(const std::array<double, 36>& stiffness,
                                                        double relativeTolerance = 1.0e-8);
export std::string elasticMatrixString(const std::array<double, 36>& matrix, double conversion);
