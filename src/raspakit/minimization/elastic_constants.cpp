module;

module elastic_constants;

import std;

import units;
import double3x3;
import generalized_hessian;
import minimization_cell_layout;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import symmetric_eigensolver;

namespace
{
constexpr std::size_t voigtSize = 6;
constexpr std::array<std::size_t, voigtSize> regularToVoigt = {0, 3, 5, 4, 2, 1};

double& at(std::array<double, 36>& matrix, std::size_t row, std::size_t column)
{
  return matrix[row * voigtSize + column];
}

double at(const std::array<double, 36>& matrix, std::size_t row, std::size_t column)
{
  return matrix[row * voigtSize + column];
}

void symmetrize(std::span<double> matrix, std::size_t size)
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

std::array<double, 6> regularCoordinates(const double3x3& matrix)
{
  return {matrix.ax, matrix.ay + matrix.bx, matrix.az + matrix.cx, matrix.by, matrix.bz + matrix.cy, matrix.cz};
}

std::vector<double> pseudoInverse(const std::vector<double>& matrix, std::size_t size, double relativeTolerance,
                                  std::size_t& discardedModes)
{
  if (size == 0) return {};
  SymmetricEigenSystem eigensystem = diagonalizeSymmetric(matrix, size);
  const double spectralScale = std::ranges::fold_left(eigensystem.eigenvalues, 0.0, [](double value, double eigenvalue)
                                                      { return std::max(value, std::abs(eigenvalue)); });
  const double threshold = relativeTolerance * spectralScale;
  if (eigensystem.eigenvalues.front() < -threshold)
  {
    throw std::runtime_error(
        std::format("Cannot compute elastic constants: internal Hessian has a negative mode ({:.8e})",
                    eigensystem.eigenvalues.front()));
  }

  std::vector<double> inverse(size * size, 0.0);
  discardedModes = 0;
  for (std::size_t mode = 0; mode < size; ++mode)
  {
    const double eigenvalue = eigensystem.eigenvalues[mode];
    if (eigenvalue <= threshold)
    {
      ++discardedModes;
      continue;
    }
    for (std::size_t row = 0; row < size; ++row)
    {
      for (std::size_t column = 0; column < size; ++column)
      {
        inverse[row * size + column] +=
            eigensystem.eigenvector(row, mode) * eigensystem.eigenvector(column, mode) / eigenvalue;
      }
    }
  }
  return inverse;
}

bool invertStiffness(const std::array<double, 36>& stiffness, double relativeTolerance,
                     std::array<double, 36>& compliance, std::array<double, 6>& eigenvalues)
{
  SymmetricEigenSystem eigensystem = diagonalizeSymmetric(stiffness, voigtSize);
  std::ranges::copy(eigensystem.eigenvalues, eigenvalues.begin());
  const double spectralScale = std::ranges::fold_left(eigensystem.eigenvalues, 0.0, [](double value, double eigenvalue)
                                                      { return std::max(value, std::abs(eigenvalue)); });
  const double threshold = relativeTolerance * spectralScale;
  if (std::ranges::any_of(eigensystem.eigenvalues,
                          [threshold](double eigenvalue) { return std::abs(eigenvalue) <= threshold; }))
  {
    return false;
  }

  compliance.fill(0.0);
  for (std::size_t mode = 0; mode < voigtSize; ++mode)
  {
    for (std::size_t row = 0; row < voigtSize; ++row)
    {
      for (std::size_t column = 0; column < voigtSize; ++column)
      {
        at(compliance, row, column) +=
            eigensystem.eigenvector(row, mode) * eigensystem.eigenvector(column, mode) / eigensystem.eigenvalues[mode];
      }
    }
  }
  return true;
}

std::string matrixString(const std::array<double, 36>& matrix, double conversion)
{
  std::string output;
  for (std::size_t row = 0; row < voigtSize; ++row)
  {
    for (std::size_t column = 0; column < voigtSize; ++column)
    {
      output += std::format("{: 15.7e}{}", conversion * at(matrix, row, column), column + 1 == voigtSize ? "\n" : " ");
    }
  }
  return output;
}
}  // namespace

ElasticConstantsResult computeElasticConstants(const System& system, double relativeEigenvalueTolerance)
{
  if (!(relativeEigenvalueTolerance > 0.0) || !std::isfinite(relativeEigenvalueTolerance))
  {
    throw std::invalid_argument("Elastic-constant eigenvalue tolerance must be finite and positive");
  }
  if (system.hasExternalField)
  {
    throw std::runtime_error("Elastic constants do not support an external field");
  }

  const double externalPressure = system.pressure;
  System derivativeSystem = system;
  derivativeSystem.pressure = 0.0;
  derivativeSystem.cellMinimizationType = CellMinimizationType::Regular;

  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(CellMinimizationType::Regular, derivativeSystem.monoclinicAngleType);
  const std::size_t numberOfFlexibleFrameworkAtoms = derivativeSystem.framework && !derivativeSystem.framework->rigid
                                                         ? derivativeSystem.spanOfFrameworkAtoms().size()
                                                         : 0;
  const MinimizationDofLayout layout = buildMinimizationDofLayout(
      derivativeSystem.moleculeData, derivativeSystem.components, numberOfFlexibleFrameworkAtoms, cellLayout.size());
  const std::size_t numberOfPositionDofs = layout.numberOfPositionDofs();
  const std::size_t numberOfDofs = layout.numDofs();

  GeneralizedHessian hessian(numberOfDofs, 0);
  std::vector<double> gradient(numberOfDofs, 0.0);
  DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = true};
  DerivativeResults derivatives{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(derivativeSystem, layout, capabilities, derivatives);

  std::vector<double> internalHessian(numberOfPositionDofs * numberOfPositionDofs);
  std::vector<double> positionStrain(numberOfPositionDofs * voigtSize);
  for (std::size_t row = 0; row < numberOfPositionDofs; ++row)
  {
    for (std::size_t column = 0; column < numberOfPositionDofs; ++column)
    {
      internalHessian[row * numberOfPositionDofs + column] = hessian.positionPosition()[row * numberOfDofs + column];
    }
    for (std::size_t strain = 0; strain < voigtSize; ++strain)
    {
      positionStrain[row * voigtSize + strain] =
          hessian.positionPosition()[row * numberOfDofs + *layout.cellDof(regularToVoigt[strain])];
    }
  }
  symmetrize(internalHessian, numberOfPositionDofs);

  ElasticConstantsResult result{};
  std::vector<double> inverseInternal =
      pseudoInverse(internalHessian, numberOfPositionDofs, relativeEigenvalueTolerance, result.discardedInternalModes);

  const double inverseVolume = 1.0 / derivativeSystem.simulationBox.volume;
  for (std::size_t row = 0; row < voigtSize; ++row)
  {
    const std::size_t regularRow = regularToVoigt[row];
    for (std::size_t column = 0; column < voigtSize; ++column)
    {
      const std::size_t regularColumn = regularToVoigt[column];
      const double logStrainHessian =
          hessian.positionPosition()[*layout.cellDof(regularRow) * numberOfDofs + *layout.cellDof(regularColumn)];
      const std::array<double, 6> secondCoordinates =
          regularCoordinates(cellStrainSecondDerivative(cellLayout, regularRow, regularColumn));
      double geometricTerm = 0.0;
      for (std::size_t coordinate = 0; coordinate < voigtSize; ++coordinate)
      {
        geometricTerm += secondCoordinates[coordinate] * gradient[*layout.cellDof(coordinate)];
      }
      // Cell minimization uses F=exp(L). Remove dU/dL : d2F/dL2 to obtain
      // the infinitesimal-strain (F=I+epsilon) Born tensor used by RASPA2.
      at(result.born, row, column) = (logStrainHessian - geometricTerm) * inverseVolume;

      double relaxation = 0.0;
      for (std::size_t i = 0; i < numberOfPositionDofs; ++i)
      {
        for (std::size_t j = 0; j < numberOfPositionDofs; ++j)
        {
          relaxation += positionStrain[i * voigtSize + row] * inverseInternal[i * numberOfPositionDofs + j] *
                        positionStrain[j * voigtSize + column];
        }
      }
      at(result.relaxation, row, column) = relaxation * inverseVolume;
    }
  }
  symmetrize(result.born, voigtSize);
  symmetrize(result.relaxation, voigtSize);

  for (std::size_t direction = 0; direction < 3; ++direction)
  {
    at(result.pressureCorrection, direction, direction) = -externalPressure;
  }
  // Born-Klein hydrostatic correction in engineering-shear Voigt notation.
  for (std::size_t shear = 3; shear < voigtSize; ++shear)
  {
    at(result.pressureCorrection, shear, shear) = -0.5 * externalPressure;
  }
  for (std::size_t index = 0; index < result.stiffness.size(); ++index)
  {
    result.stiffness[index] = result.born[index] - result.relaxation[index] + result.pressureCorrection[index];
  }
  symmetrize(result.stiffness, voigtSize);

  result.complianceAvailable =
      invertStiffness(result.stiffness, relativeEigenvalueTolerance, result.compliance, result.stabilityEigenvalues);
  result.bulkModulusVoigt =
      (at(result.stiffness, 0, 0) + at(result.stiffness, 1, 1) + at(result.stiffness, 2, 2) +
       2.0 * (at(result.stiffness, 0, 1) + at(result.stiffness, 0, 2) + at(result.stiffness, 1, 2))) /
      9.0;
  result.shearModulusVoigt =
      (at(result.stiffness, 0, 0) + at(result.stiffness, 1, 1) + at(result.stiffness, 2, 2) -
       at(result.stiffness, 0, 1) - at(result.stiffness, 0, 2) - at(result.stiffness, 1, 2) +
       3.0 * (at(result.stiffness, 3, 3) + at(result.stiffness, 4, 4) + at(result.stiffness, 5, 5))) /
      15.0;

  if (result.complianceAvailable)
  {
    const double bulkDenominator =
        at(result.compliance, 0, 0) + at(result.compliance, 1, 1) + at(result.compliance, 2, 2) +
        2.0 * (at(result.compliance, 0, 1) + at(result.compliance, 0, 2) + at(result.compliance, 1, 2));
    const double shearDenominator =
        4.0 * (at(result.compliance, 0, 0) + at(result.compliance, 1, 1) + at(result.compliance, 2, 2) -
               at(result.compliance, 0, 1) - at(result.compliance, 0, 2) - at(result.compliance, 1, 2)) +
        3.0 * (at(result.compliance, 3, 3) + at(result.compliance, 4, 4) + at(result.compliance, 5, 5));
    result.bulkModulusReuss = 1.0 / bulkDenominator;
    result.shearModulusReuss = 15.0 / shearDenominator;
    result.bulkModulusHill = 0.5 * (result.bulkModulusVoigt + result.bulkModulusReuss);
    result.shearModulusHill = 0.5 * (result.shearModulusVoigt + result.shearModulusReuss);
    for (std::size_t loading = 0; loading < 3; ++loading)
    {
      result.youngModuli[loading] = 1.0 / at(result.compliance, loading, loading);
      for (std::size_t transverse = 0; transverse < 3; ++transverse)
      {
        result.poissonRatios[loading * 3 + transverse] =
            loading == transverse
                ? 0.0
                : -at(result.compliance, transverse, loading) / at(result.compliance, loading, loading);
      }
    }
  }
  return result;
}

std::string writeElasticConstants(const ElasticConstantsResult& result)
{
  const bool reduced = Units::unitSystem == Units::System::ReducedUnits;
  const double pressureConversion = reduced ? 1.0 : 1.0e-9 * Units::PressureConversionFactor;
  const double complianceConversion = reduced ? 1.0 : 1.0 / pressureConversion;
  const std::string_view pressureUnit = reduced ? "reduced pressure" : "GPa";
  const std::string_view complianceUnit = reduced ? "inverse reduced pressure" : "GPa^-1";

  std::string output = std::format(
      "\nStatic elastic constants at 0 K\n"
      "Voigt order: xx yy zz yz xz xy; shear strains use engineering convention.\n"
      "Discarded internal zero modes: {}\n\n"
      "Born term [{}]\n{}\n"
      "Internal-relaxation term [{}]\n{}\n"
      "Hydrostatic pressure correction [{}]\n{}\n"
      "Relaxed elastic constants [{}]\n{}\n",
      result.discardedInternalModes, pressureUnit, matrixString(result.born, pressureConversion), pressureUnit,
      matrixString(result.relaxation, pressureConversion), pressureUnit,
      matrixString(result.pressureCorrection, pressureConversion), pressureUnit,
      matrixString(result.stiffness, pressureConversion));

  if (result.complianceAvailable)
  {
    output += std::format(
        "Elastic compliance [{}]\n{}\n"
        "Young moduli: {: .7e} {: .7e} {: .7e} [{}]\n"
        "Poisson ratios (rows: loading x/y/z; columns: transverse x/y/z):\n"
        "{: .7e} {: .7e} {: .7e}\n{: .7e} {: .7e} {: .7e}\n{: .7e} {: .7e} {: .7e}\n"
        "Bulk modulus (Voigt/Reuss/Hill): {: .7e} {: .7e} {: .7e} [{}]\n"
        "Shear modulus (Voigt/Reuss/Hill): {: .7e} {: .7e} {: .7e} [{}]\n",
        complianceUnit, matrixString(result.compliance, complianceConversion),
        pressureConversion * result.youngModuli[0], pressureConversion * result.youngModuli[1],
        pressureConversion * result.youngModuli[2], pressureUnit, result.poissonRatios[0], result.poissonRatios[1],
        result.poissonRatios[2], result.poissonRatios[3], result.poissonRatios[4], result.poissonRatios[5],
        result.poissonRatios[6], result.poissonRatios[7], result.poissonRatios[8],
        pressureConversion * result.bulkModulusVoigt, pressureConversion * result.bulkModulusReuss,
        pressureConversion * result.bulkModulusHill, pressureUnit, pressureConversion * result.shearModulusVoigt,
        pressureConversion * result.shearModulusReuss, pressureConversion * result.shearModulusHill, pressureUnit);
  }
  else
  {
    output += "Elastic compliance and derived moduli are unavailable because the stiffness matrix is singular.\n";
  }
  output += std::format("Born stability eigenvalues [{}]:", pressureUnit);
  for (double eigenvalue : result.stabilityEigenvalues)
  {
    output += std::format(" {: .7e}", pressureConversion * eigenvalue);
  }
  output += "\n";
  return output;
}
