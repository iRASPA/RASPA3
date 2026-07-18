module;

module elastic_constants;

import std;

import units;
import double3x3;
import double3;
import generalized_hessian;
import minimization_cell_layout;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import symmetric_eigensolver;
import elastic_tensor;
import atom;
import atom_dynamics;
import molecule;
import component;

namespace
{
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
}  // namespace

std::array<double, 36> computeAffineBornTensor(const System& system)
{
  if (system.hasExternalField) throw std::runtime_error("Elastic constants do not support an external field");

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
  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = true};
  DerivativeResults derivatives{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(derivativeSystem, layout, capabilities, derivatives);

  std::array<double, 36> born{};
  const double inverseVolume = 1.0 / derivativeSystem.simulationBox.volume;
  for (std::size_t row = 0; row < elasticVoigtSize; ++row)
  {
    const std::size_t regularRow = elasticRegularToVoigt[row];
    for (std::size_t column = 0; column < elasticVoigtSize; ++column)
    {
      const std::size_t regularColumn = elasticRegularToVoigt[column];
      const double logStrainHessian =
          hessian.positionPosition()[*layout.cellDof(regularRow) * layout.numDofs() + *layout.cellDof(regularColumn)];
      const std::array<double, 6> secondCoordinates =
          regularCoordinates(cellStrainSecondDerivative(cellLayout, regularRow, regularColumn));
      double geometricTerm = 0.0;
      for (std::size_t coordinate = 0; coordinate < elasticVoigtSize; ++coordinate)
        geometricTerm += secondCoordinates[coordinate] * gradient[*layout.cellDof(coordinate)];
      elasticAt(born, row, column) = (logStrainHessian - geometricTerm) * inverseVolume;
    }
  }
  symmetrizeElasticTensor(born);
  return born;
}

double3x3 computeMolecularKineticVirial(const System& system)
{
  double3x3 stress{};
  const auto add = [&stress](double mass, const double3& velocity)
  {
    stress.ax += mass * velocity.x * velocity.x;
    stress.ay += mass * velocity.x * velocity.y;
    stress.az += mass * velocity.x * velocity.z;
    stress.bx += mass * velocity.y * velocity.x;
    stress.by += mass * velocity.y * velocity.y;
    stress.bz += mass * velocity.y * velocity.z;
    stress.cx += mass * velocity.z * velocity.x;
    stress.cy += mass * velocity.z * velocity.y;
    stress.cz += mass * velocity.z * velocity.z;
  };
  const std::span<const Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  const std::span<const AtomDynamics> moleculeDynamics = system.spanOfMoleculeDynamics();
  const std::span<const GroupState> groupData = system.spanOfGroupData();
  std::size_t atomIndex = 0;
  std::size_t groupIndex = 0;
  for (const Molecule& molecule : system.moleculeData)
  {
    const Component& component = system.components[molecule.componentId];
    if (component.rigid)
    {
      add(molecule.mass, molecule.velocity);
    }
    else if (component.isSemiFlexible() && !groupData.empty())
    {
      // Semi-flexible molecule: rigid groups contribute their center-of-mass momentum flux,
      // flexible atoms contribute individually (matching the barostat coupling).
      std::size_t rigidRank = 0;
      for (const MoleculeGroup& group : component.groups)
      {
        if (group.rigid)
        {
          add(group.mass, groupData[groupIndex + rigidRank].velocity);
          ++rigidRank;
        }
      }
      for (std::size_t i = 0; i < molecule.numberOfAtoms; ++i)
      {
        if (!component.rigidGroupContaining(i).has_value())
          add(system.forceField.pseudoAtoms[moleculeAtoms[atomIndex + i].type].mass,
              moleculeDynamics[atomIndex + i].velocity);
      }
      groupIndex += component.numberOfRigidGroups();
    }
    else
    {
      for (std::size_t i = 0; i < molecule.numberOfAtoms; ++i)
        add(system.forceField.pseudoAtoms[moleculeAtoms[atomIndex + i].type].mass,
            moleculeDynamics[atomIndex + i].velocity);
    }
    atomIndex += molecule.numberOfAtoms;
  }
  if (system.framework && !system.framework->rigid)
  {
    const std::span<const Atom> atoms = system.spanOfFrameworkAtoms();
    const std::span<const AtomDynamics> dynamics = system.spanOfFrameworkDynamics();
    for (std::size_t i = 0; i < atoms.size(); ++i)
      add(system.forceField.pseudoAtoms[atoms[i].type].mass, dynamics[i].velocity);
  }
  return stress;
}

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
  std::vector<double> positionStrain(numberOfPositionDofs * elasticVoigtSize);
  for (std::size_t row = 0; row < numberOfPositionDofs; ++row)
  {
    for (std::size_t column = 0; column < numberOfPositionDofs; ++column)
    {
      internalHessian[row * numberOfPositionDofs + column] = hessian.positionPosition()[row * numberOfDofs + column];
    }
    for (std::size_t strain = 0; strain < elasticVoigtSize; ++strain)
    {
      positionStrain[row * elasticVoigtSize + strain] =
          hessian.positionPosition()[row * numberOfDofs + *layout.cellDof(elasticRegularToVoigt[strain])];
    }
  }
  symmetrizeElasticTensor(internalHessian, numberOfPositionDofs);

  ElasticConstantsResult result{};
  std::vector<double> inverseInternal =
      pseudoInverse(internalHessian, numberOfPositionDofs, relativeEigenvalueTolerance, result.discardedInternalModes);

  const double inverseVolume = 1.0 / derivativeSystem.simulationBox.volume;
  for (std::size_t row = 0; row < elasticVoigtSize; ++row)
  {
    const std::size_t regularRow = elasticRegularToVoigt[row];
    for (std::size_t column = 0; column < elasticVoigtSize; ++column)
    {
      const std::size_t regularColumn = elasticRegularToVoigt[column];
      const double logStrainHessian =
          hessian.positionPosition()[*layout.cellDof(regularRow) * numberOfDofs + *layout.cellDof(regularColumn)];
      const std::array<double, 6> secondCoordinates =
          regularCoordinates(cellStrainSecondDerivative(cellLayout, regularRow, regularColumn));
      double geometricTerm = 0.0;
      for (std::size_t coordinate = 0; coordinate < elasticVoigtSize; ++coordinate)
      {
        geometricTerm += secondCoordinates[coordinate] * gradient[*layout.cellDof(coordinate)];
      }
      // Cell minimization uses F=exp(L). Remove dU/dL : d2F/dL2 to obtain
      // the infinitesimal-strain (F=I+epsilon) Born tensor used by RASPA2.
      elasticAt(result.born, row, column) = (logStrainHessian - geometricTerm) * inverseVolume;

      double relaxation = 0.0;
      for (std::size_t i = 0; i < numberOfPositionDofs; ++i)
      {
        for (std::size_t j = 0; j < numberOfPositionDofs; ++j)
        {
          relaxation += positionStrain[i * elasticVoigtSize + row] * inverseInternal[i * numberOfPositionDofs + j] *
                        positionStrain[j * elasticVoigtSize + column];
        }
      }
      elasticAt(result.relaxation, row, column) = relaxation * inverseVolume;
    }
  }
  symmetrizeElasticTensor(result.born);
  symmetrizeElasticTensor(result.relaxation);

  for (std::size_t direction = 0; direction < 3; ++direction)
  {
    elasticAt(result.pressureCorrection, direction, direction) = -externalPressure;
  }
  // Born-Klein hydrostatic correction in engineering-shear Voigt notation.
  for (std::size_t shear = 3; shear < elasticVoigtSize; ++shear)
  {
    elasticAt(result.pressureCorrection, shear, shear) = -0.5 * externalPressure;
  }
  for (std::size_t index = 0; index < result.stiffness.size(); ++index)
  {
    result.stiffness[index] = result.born[index] - result.relaxation[index] + result.pressureCorrection[index];
  }
  symmetrizeElasticTensor(result.stiffness);

  const ElasticDerivedProperties derived = deriveElasticProperties(result.stiffness, relativeEigenvalueTolerance);
  result.compliance = derived.compliance;
  result.stabilityEigenvalues = derived.stabilityEigenvalues;
  result.youngModuli = derived.youngModuli;
  result.poissonRatios = derived.poissonRatios;
  result.bulkModulusVoigt = derived.bulkModulusVoigt;
  result.shearModulusVoigt = derived.shearModulusVoigt;
  result.bulkModulusReuss = derived.bulkModulusReuss;
  result.shearModulusReuss = derived.shearModulusReuss;
  result.bulkModulusHill = derived.bulkModulusHill;
  result.shearModulusHill = derived.shearModulusHill;
  result.complianceAvailable = derived.complianceAvailable;
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
      result.discardedInternalModes, pressureUnit, elasticMatrixString(result.born, pressureConversion), pressureUnit,
      elasticMatrixString(result.relaxation, pressureConversion), pressureUnit,
      elasticMatrixString(result.pressureCorrection, pressureConversion), pressureUnit,
      elasticMatrixString(result.stiffness, pressureConversion));

  if (result.complianceAvailable)
  {
    output += std::format(
        "Elastic compliance [{}]\n{}\n"
        "Young moduli: {: .7e} {: .7e} {: .7e} [{}]\n"
        "Poisson ratios (rows: loading x/y/z; columns: transverse x/y/z):\n"
        "{: .7e} {: .7e} {: .7e}\n{: .7e} {: .7e} {: .7e}\n{: .7e} {: .7e} {: .7e}\n"
        "Bulk modulus (Voigt/Reuss/Hill): {: .7e} {: .7e} {: .7e} [{}]\n"
        "Shear modulus (Voigt/Reuss/Hill): {: .7e} {: .7e} {: .7e} [{}]\n",
        complianceUnit, elasticMatrixString(result.compliance, complianceConversion),
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
