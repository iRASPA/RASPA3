module;

module elastic_tensor;

import std;
import symmetric_eigensolver;

double& elasticAt(std::array<double, 36>& matrix, std::size_t row, std::size_t column)
{
  return matrix[row * elasticVoigtSize + column];
}

double elasticAt(const std::array<double, 36>& matrix, std::size_t row, std::size_t column)
{
  return matrix[row * elasticVoigtSize + column];
}

void symmetrizeElasticTensor(std::span<double> matrix, std::size_t size)
{
  for (std::size_t row = 0; row < size; ++row)
  {
    for (std::size_t column = row + 1; column < size; ++column)
    {
      const double value = 0.5 * (matrix[row * size + column] + matrix[column * size + row]);
      matrix[row * size + column] = value;
      matrix[column * size + row] = value;
    }
  }
}

ElasticDerivedProperties deriveElasticProperties(const std::array<double, 36>& stiffness, double relativeTolerance)
{
  if (!(relativeTolerance > 0.0) || !std::isfinite(relativeTolerance))
    throw std::invalid_argument("Elastic eigenvalue tolerance must be finite and positive");

  ElasticDerivedProperties result{};
  SymmetricEigenSystem eigensystem = diagonalizeSymmetric(stiffness, elasticVoigtSize);
  std::ranges::copy(eigensystem.eigenvalues, result.stabilityEigenvalues.begin());
  const double spectralScale = std::ranges::fold_left(eigensystem.eigenvalues, 0.0, [](double value, double eigenvalue)
                                                      { return std::max(value, std::abs(eigenvalue)); });
  const double threshold = relativeTolerance * spectralScale;
  result.complianceAvailable = std::ranges::none_of(
      eigensystem.eigenvalues, [threshold](double eigenvalue) { return std::abs(eigenvalue) <= threshold; });
  if (result.complianceAvailable)
  {
    for (std::size_t mode = 0; mode < elasticVoigtSize; ++mode)
      for (std::size_t row = 0; row < elasticVoigtSize; ++row)
        for (std::size_t column = 0; column < elasticVoigtSize; ++column)
          elasticAt(result.compliance, row, column) += eigensystem.eigenvector(row, mode) *
                                                       eigensystem.eigenvector(column, mode) /
                                                       eigensystem.eigenvalues[mode];
  }

  result.bulkModulusVoigt =
      (elasticAt(stiffness, 0, 0) + elasticAt(stiffness, 1, 1) + elasticAt(stiffness, 2, 2) +
       2.0 * (elasticAt(stiffness, 0, 1) + elasticAt(stiffness, 0, 2) + elasticAt(stiffness, 1, 2))) /
      9.0;
  result.shearModulusVoigt =
      (elasticAt(stiffness, 0, 0) + elasticAt(stiffness, 1, 1) + elasticAt(stiffness, 2, 2) -
       elasticAt(stiffness, 0, 1) - elasticAt(stiffness, 0, 2) - elasticAt(stiffness, 1, 2) +
       3.0 * (elasticAt(stiffness, 3, 3) + elasticAt(stiffness, 4, 4) + elasticAt(stiffness, 5, 5))) /
      15.0;

  if (result.complianceAvailable)
  {
    const double bulkDenominator = elasticAt(result.compliance, 0, 0) + elasticAt(result.compliance, 1, 1) +
                                   elasticAt(result.compliance, 2, 2) +
                                   2.0 * (elasticAt(result.compliance, 0, 1) + elasticAt(result.compliance, 0, 2) +
                                          elasticAt(result.compliance, 1, 2));
    const double shearDenominator = 4.0 * (elasticAt(result.compliance, 0, 0) + elasticAt(result.compliance, 1, 1) +
                                           elasticAt(result.compliance, 2, 2) - elasticAt(result.compliance, 0, 1) -
                                           elasticAt(result.compliance, 0, 2) - elasticAt(result.compliance, 1, 2)) +
                                    3.0 * (elasticAt(result.compliance, 3, 3) + elasticAt(result.compliance, 4, 4) +
                                           elasticAt(result.compliance, 5, 5));
    result.bulkModulusReuss = 1.0 / bulkDenominator;
    result.shearModulusReuss = 15.0 / shearDenominator;
    result.bulkModulusHill = 0.5 * (result.bulkModulusVoigt + result.bulkModulusReuss);
    result.shearModulusHill = 0.5 * (result.shearModulusVoigt + result.shearModulusReuss);
    for (std::size_t loading = 0; loading < 3; ++loading)
    {
      result.youngModuli[loading] = 1.0 / elasticAt(result.compliance, loading, loading);
      for (std::size_t transverse = 0; transverse < 3; ++transverse)
        result.poissonRatios[loading * 3 + transverse] =
            loading == transverse
                ? 0.0
                : -elasticAt(result.compliance, transverse, loading) / elasticAt(result.compliance, loading, loading);
    }
  }
  return result;
}

std::string elasticMatrixString(const std::array<double, 36>& matrix, double conversion)
{
  std::string output;
  for (std::size_t row = 0; row < elasticVoigtSize; ++row)
  {
    for (std::size_t column = 0; column < elasticVoigtSize; ++column)
      output += std::format("{: 15.7e}{}", conversion * elasticAt(matrix, row, column),
                            column + 1 == elasticVoigtSize ? "\n" : " ");
  }
  return output;
}
