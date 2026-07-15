module;

module normal_modes;

import std;

import units;
import double3;
import double3x3;
import atom;
import molecule;
import component;
import forcefield;
import simulationbox;
import skelement;
import generalized_hessian;
import minimization_cell_layout;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import minimization_generalized_coordinates;
import symmetric_eigensolver;

namespace
{
/**
 * Inverse-square-root mass metric in the generalized DOF basis. Cartesian and
 * rigid center-of-mass DOFs carry scalar weights 1/sqrt(m); rigid orientation
 * DOFs carry a 3x3 block, the pseudo-inverse square root of the space-frame
 * inertia tensor.
 */
struct MassMetric
{
  // Zero entries mark DOFs that belong to a 3x3 orientation block.
  std::vector<double> scalarInverseSqrt;
  // Base DOF index plus row-major inverse-square-root inertia block.
  std::vector<std::pair<std::size_t, std::array<double, 9>>> blocks;
  std::size_t discardedRotationalDofs{};
};

double checkedInverseSqrtMass(double mass, std::string_view what)
{
  if (!(mass > 0.0) || !std::isfinite(mass))
  {
    throw std::runtime_error(
        std::format("Normal-mode analysis requires positive masses, but {} has mass {}", what, mass));
  }
  return 1.0 / std::sqrt(mass);
}

std::array<double, 9> inverseSqrtInertia(const double3x3& inertia, std::size_t& discardedRotationalDofs)
{
  const std::array<double, 9> matrix = {inertia.ax, inertia.bx, inertia.cx, inertia.ay, inertia.by,
                                        inertia.cy, inertia.az, inertia.bz, inertia.cz};
  const SymmetricEigenSystem eigensystem = diagonalizeSymmetric(matrix, 3);

  // Same near-zero inertia criterion as the component rotational-DOF count.
  const double trace = inertia.ax + inertia.by + inertia.cz;
  const double threshold = std::max(1.0e-2, trace) * 1.0e-5;

  std::array<double, 9> inverseSqrt{};
  for (std::size_t mode = 0; mode < 3; ++mode)
  {
    const double eigenvalue = eigensystem.eigenvalues[mode];
    if (eigenvalue <= threshold)
    {
      ++discardedRotationalDofs;
      continue;
    }
    const double weight = 1.0 / std::sqrt(eigenvalue);
    for (std::size_t row = 0; row < 3; ++row)
    {
      for (std::size_t column = 0; column < 3; ++column)
      {
        inverseSqrt[row * 3 + column] += weight * eigensystem.eigenvector(row, mode) * eigensystem.eigenvector(column, mode);
      }
    }
  }
  return inverseSqrt;
}

MassMetric buildMassMetric(const System& system, const MinimizationDofLayout& layout)
{
  MassMetric metric{};
  metric.scalarInverseSqrt.assign(layout.numDofs(), 0.0);

  const std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  for (std::size_t atom = 0; atom < layout.numberOfFrameworkAtoms(); ++atom)
  {
    const std::size_t type = static_cast<std::size_t>(frameworkAtoms[atom].type);
    const double weight = checkedInverseSqrtMass(system.forceField.pseudoAtoms[type].mass,
                                                 std::format("framework atom {} ({})", atom,
                                                             system.forceField.pseudoAtoms[type].name));
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      metric.scalarInverseSqrt[*layout.frameworkAtomDof(atom, static_cast<MinimizationDofAxis>(axis))] = weight;
    }
  }

  const std::span<const Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = system.moleculeData[moleculeIndex];
    const Component& component = system.components[molecule.componentId];

    if (component.rigid)
    {
      const double weight =
          checkedInverseSqrtMass(molecule.mass, std::format("rigid molecule {} ({})", moleculeIndex, component.name));
      metric.scalarInverseSqrt[*layout.rigidMoleculeDof(moleculeIndex, RigidDof::ComX)] = weight;
      metric.scalarInverseSqrt[*layout.rigidMoleculeDof(moleculeIndex, RigidDof::ComY)] = weight;
      metric.scalarInverseSqrt[*layout.rigidMoleculeDof(moleculeIndex, RigidDof::ComZ)] = weight;

      // The orientation tangent DOFs are lab-frame rotation-vector components, so the
      // kinetic metric is the space-frame inertia tensor built from the current offsets.
      double3x3 inertia{};
      for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
      {
        const Atom& atom = moleculeAtoms[molecule.atomIndex + localAtom];
        const double mass = system.forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].mass;
        const double3 s = atom.position - molecule.centerOfMassPosition;
        inertia.ax += mass * (s.y * s.y + s.z * s.z);
        inertia.by += mass * (s.x * s.x + s.z * s.z);
        inertia.cz += mass * (s.x * s.x + s.y * s.y);
        inertia.bx -= mass * s.x * s.y;
        inertia.ay -= mass * s.x * s.y;
        inertia.cx -= mass * s.x * s.z;
        inertia.az -= mass * s.x * s.z;
        inertia.cy -= mass * s.y * s.z;
        inertia.bz -= mass * s.y * s.z;
      }
      metric.blocks.emplace_back(*layout.rigidMoleculeDof(moleculeIndex, RigidDof::OriX),
                                 inverseSqrtInertia(inertia, metric.discardedRotationalDofs));
    }
    else
    {
      for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
      {
        const std::size_t type = static_cast<std::size_t>(moleculeAtoms[molecule.atomIndex + localAtom].type);
        const double weight = checkedInverseSqrtMass(
            system.forceField.pseudoAtoms[type].mass,
            std::format("atom {} of molecule {} ({})", localAtom, moleculeIndex, component.name));
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
          metric.scalarInverseSqrt[*layout.flexibleAtomDof(moleculeIndex, localAtom,
                                                           static_cast<MinimizationDofAxis>(axis))] = weight;
        }
      }
    }
  }
  return metric;
}

// D = W H W with W the symmetric block-diagonal inverse-square-root mass metric.
std::vector<double> applyMassMetric(std::span<const double> hessian, std::size_t size, const MassMetric& metric)
{
  std::vector<double> weighted(hessian.begin(), hessian.end());

  // Row transform: block rows and scalar rows are disjoint.
  for (std::size_t row = 0; row < size; ++row)
  {
    const double weight = metric.scalarInverseSqrt[row];
    if (weight == 0.0) continue;
    for (std::size_t column = 0; column < size; ++column)
    {
      weighted[row * size + column] *= weight;
    }
  }
  for (const auto& [base, block] : metric.blocks)
  {
    for (std::size_t column = 0; column < size; ++column)
    {
      const std::array<double, 3> rows = {weighted[(base + 0) * size + column], weighted[(base + 1) * size + column],
                                          weighted[(base + 2) * size + column]};
      for (std::size_t k = 0; k < 3; ++k)
      {
        weighted[(base + k) * size + column] =
            block[k * 3 + 0] * rows[0] + block[k * 3 + 1] * rows[1] + block[k * 3 + 2] * rows[2];
      }
    }
  }

  // Column transform.
  for (std::size_t column = 0; column < size; ++column)
  {
    const double weight = metric.scalarInverseSqrt[column];
    if (weight == 0.0) continue;
    for (std::size_t row = 0; row < size; ++row)
    {
      weighted[row * size + column] *= weight;
    }
  }
  for (const auto& [base, block] : metric.blocks)
  {
    for (std::size_t row = 0; row < size; ++row)
    {
      const std::array<double, 3> columns = {weighted[row * size + base + 0], weighted[row * size + base + 1],
                                             weighted[row * size + base + 2]};
      for (std::size_t k = 0; k < 3; ++k)
      {
        weighted[row * size + base + k] =
            columns[0] * block[0 * 3 + k] + columns[1] * block[1 * 3 + k] + columns[2] * block[2 * 3 + k];
      }
    }
  }
  return weighted;
}

// Symmetric block-diagonal metric applied to a single vector: q = W v.
std::vector<double> metricTimesVector(const MassMetric& metric, std::span<const double> vector, std::size_t size)
{
  std::vector<double> result(vector.begin(), vector.end());
  for (std::size_t i = 0; i < size; ++i)
  {
    const double weight = metric.scalarInverseSqrt[i];
    if (weight != 0.0) result[i] = weight * vector[i];
  }
  for (const auto& [base, block] : metric.blocks)
  {
    for (std::size_t k = 0; k < 3; ++k)
    {
      result[base + k] =
          block[k * 3 + 0] * vector[base + 0] + block[k * 3 + 1] * vector[base + 1] + block[k * 3 + 2] * vector[base + 2];
    }
  }
  return result;
}

void writePdbFrame(std::ostream& stream, const System& system, std::size_t modelNumber)
{
  const SimulationBox& box = system.simulationBox;
  std::print(stream, "MODEL {:>4}\n", modelNumber);
  std::print(stream, "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}\n", box.lengthA, box.lengthB, box.lengthC,
             box.angleAlpha * 180.0 / std::numbers::pi, box.angleBeta * 180.0 / std::numbers::pi,
             box.angleGamma * 180.0 / std::numbers::pi);

  std::size_t serialNumber{1};
  const auto writeAtom = [&](const Atom& atom)
  {
    const std::size_t type = static_cast<std::size_t>(atom.type);
    const std::size_t atomicNumber = system.forceField.pseudoAtoms[type].atomicNumber;
    const std::string name = std::format("{:<4}", system.forceField.pseudoAtoms[type].name);
    const std::string element = PredefinedElements::predefinedElements[atomicNumber]._chemicalSymbol;
    std::print(stream, "ATOM  {:>5} {:4}{:1}{:>3} {:1}{:>4}{:1}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4}{:>2}\n",
               serialNumber, name.substr(0, 4), ' ', " ", ' ', 0, ' ', atom.position.x, atom.position.y,
               atom.position.z, 1.0, 0.0, ' ', element);
    ++serialNumber;
  };

  for (const Atom& atom : system.spanOfFrameworkAtoms()) writeAtom(atom);
  for (const Atom& atom : system.spanOfMoleculeAtoms()) writeAtom(atom);
  stream << "ENDMDL\n";
}
}  // namespace

NormalModesResult computeNormalModes(const System& system, double relativeZeroTolerance)
{
  if (!(relativeZeroTolerance > 0.0) || !std::isfinite(relativeZeroTolerance))
  {
    throw std::invalid_argument("Normal-mode zero tolerance must be finite and positive");
  }

  System derivativeSystem = system;
  // Normal modes are evaluated at fixed cell; the layout below carries no cell DOFs.
  derivativeSystem.cellMinimizationType = CellMinimizationType::Fixed;
  const std::size_t numberOfFlexibleFrameworkAtoms = derivativeSystem.framework && !derivativeSystem.framework->rigid
                                                         ? derivativeSystem.spanOfFrameworkAtoms().size()
                                                         : 0;
  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(derivativeSystem.moleculeData, derivativeSystem.components,
                                 numberOfFlexibleFrameworkAtoms, 0);
  const std::size_t numberOfDofs = layout.numDofs();

  NormalModesResult result{};
  result.numberOfModes = numberOfDofs;
  if (numberOfDofs == 0) return result;

  GeneralizedHessian hessian(numberOfDofs, 0);
  std::vector<double> gradient(numberOfDofs, 0.0);
  DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = true};
  DerivativeResults derivatives{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(derivativeSystem, layout, capabilities, derivatives);

  std::vector<double> symmetrized(hessian.positionPosition().begin(), hessian.positionPosition().end());
  for (std::size_t row = 0; row < numberOfDofs; ++row)
  {
    for (std::size_t column = row + 1; column < numberOfDofs; ++column)
    {
      const double average = 0.5 * (symmetrized[row * numberOfDofs + column] + symmetrized[column * numberOfDofs + row]);
      symmetrized[row * numberOfDofs + column] = average;
      symmetrized[column * numberOfDofs + row] = average;
    }
  }

  const MassMetric metric = buildMassMetric(derivativeSystem, layout);
  const std::vector<double> dynamicalMatrix = applyMassMetric(symmetrized, numberOfDofs, metric);

  SymmetricEigenSystem eigensystem = diagonalizeSymmetric(dynamicalMatrix, numberOfDofs);
  result.eigenvalues = std::move(eigensystem.eigenvalues);
  result.eigenvectors = std::move(eigensystem.eigenvectors);
  result.discardedRotationalDofs = metric.discardedRotationalDofs;

  const double spectralScale = std::ranges::fold_left(result.eigenvalues, 0.0, [](double value, double eigenvalue)
                                                      { return std::max(value, std::abs(eigenvalue)); });
  const double threshold = relativeZeroTolerance * spectralScale;
  for (const double eigenvalue : result.eigenvalues)
  {
    if (eigenvalue < -threshold)
      ++result.negativeModes;
    else if (eigenvalue <= threshold)
      ++result.zeroModes;
  }
  return result;
}

std::vector<double> normalModeFrequencies(const NormalModesResult& result)
{
  const bool reduced = Units::unitSystem == Units::System::ReducedUnits;
  // Internal angular frequency (per internal time unit) to cm^-1.
  const double conversion =
      reduced ? 1.0 : 1.0 / (2.0 * std::numbers::pi * 100.0 * Units::SpeedOfLight * Units::TimeConversionFactor);

  std::vector<double> frequencies;
  frequencies.reserve(result.eigenvalues.size());
  for (const double eigenvalue : result.eigenvalues)
  {
    const double magnitude = std::sqrt(std::abs(eigenvalue)) * conversion;
    frequencies.push_back(eigenvalue < 0.0 ? -magnitude : magnitude);
  }
  return frequencies;
}

std::string writeNormalModes(const NormalModesResult& result)
{
  const bool reduced = Units::unitSystem == Units::System::ReducedUnits;

  std::string output = std::format(
      "\nGamma-point normal-mode analysis\n"
      "Number of modes: {} (negative: {}, zero: {})\n"
      "Discarded rigid-rotation DOFs (zero inertia): {}\n"
      "Negative eigenvalues denote imaginary modes and are printed as negative frequencies.\n\n",
      result.numberOfModes, result.negativeModes, result.zeroModes, result.discardedRotationalDofs);

  if (reduced)
  {
    output += "mode       omega^2 [reduced]      frequency [reduced]\n";
    const std::vector<double> frequencies = normalModeFrequencies(result);
    for (std::size_t mode = 0; mode < result.eigenvalues.size(); ++mode)
    {
      output += std::format("{:5d}  {: 20.10e}  {: 20.10e}\n", mode, result.eigenvalues[mode], frequencies[mode]);
    }
    return output;
  }

  // Internal angular frequency (rad/ps) to SI (rad/s).
  const double toRadiansPerSecond = 1.0 / Units::TimeConversionFactor;
  const double toTHz = toRadiansPerSecond / (2.0 * std::numbers::pi * 1.0e12);
  const double toWavenumber = toRadiansPerSecond / (2.0 * std::numbers::pi * 100.0 * Units::SpeedOfLight);
  const double toMilliEv =
      1.0e3 * (Units::PlanckConstant / (2.0 * std::numbers::pi)) * toRadiansPerSecond / Units::ElectronicChargeUnit;

  output += "mode       omega^2 [ps^-2]       frequency [THz]    wavenumber [cm^-1]         energy [meV]\n";
  for (std::size_t mode = 0; mode < result.eigenvalues.size(); ++mode)
  {
    const double eigenvalue = result.eigenvalues[mode];
    const double sign = eigenvalue < 0.0 ? -1.0 : 1.0;
    const double omega = std::sqrt(std::abs(eigenvalue));
    output += std::format("{:5d}  {: 20.10e}  {: 20.10e}  {: 20.10e}  {: 20.10e}\n", mode, eigenvalue,
                          sign * omega * toTHz, sign * omega * toWavenumber, sign * omega * toMilliEv);
  }
  return output;
}

void writeNormalModeMovies(const System& system, const NormalModesResult& result, std::size_t systemIndex,
                           std::size_t numberOfPeriods, std::size_t pointsPerPeriod, double amplitude,
                           const std::filesystem::path& directory)
{
  if (result.numberOfModes == 0) return;
  if (numberOfPeriods == 0 || pointsPerPeriod == 0)
  {
    throw std::invalid_argument("Normal-mode movies require positive numberOfPeriods and pointsPerPeriod");
  }
  if (!(amplitude > 0.0) || !std::isfinite(amplitude))
  {
    throw std::invalid_argument("Normal-mode movie amplitude must be finite and positive");
  }

  System referenceSystem = system;
  referenceSystem.cellMinimizationType = CellMinimizationType::Fixed;
  const std::size_t numberOfFlexibleFrameworkAtoms = referenceSystem.framework && !referenceSystem.framework->rigid
                                                         ? referenceSystem.spanOfFrameworkAtoms().size()
                                                         : 0;
  const MinimizationDofLayout layout = buildMinimizationDofLayout(
      referenceSystem.moleculeData, referenceSystem.components, numberOfFlexibleFrameworkAtoms, 0);
  const std::size_t numberOfDofs = layout.numDofs();
  if (numberOfDofs == 0) return;
  if (result.eigenvectors.size() != numberOfDofs * numberOfDofs)
  {
    throw std::invalid_argument("Normal-mode movies: eigenvector layout does not match the system");
  }

  const MassMetric metric = buildMassMetric(referenceSystem, layout);
  const std::vector<double> frequencies = normalModeFrequencies(result);
  const std::span<const Atom> moleculeAtoms = referenceSystem.spanOfMoleculeAtoms();

  std::filesystem::create_directories(directory);

  const std::size_t totalFrames = numberOfPeriods * pointsPerPeriod;
  for (std::size_t mode = 0; mode < numberOfDofs; ++mode)
  {
    std::vector<double> eigenvector(numberOfDofs);
    for (std::size_t dof = 0; dof < numberOfDofs; ++dof)
    {
      eigenvector[dof] = result.eigenvectors[dof + numberOfDofs * mode];
    }
    // Mass-weighted eigenvector xi maps to a generalized displacement q = W xi.
    const std::vector<double> generalized = metricTimesVector(metric, eigenvector, numberOfDofs);

    // Largest atomic Cartesian displacement produced by the mode pattern.
    double maxDisplacementSquared = 0.0;
    for (std::size_t atom = 0; atom < layout.numberOfFrameworkAtoms(); ++atom)
    {
      double3 displacement{};
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        if (const auto dof = layout.frameworkAtomDof(atom, static_cast<MinimizationDofAxis>(axis)))
        {
          (&displacement.x)[axis] = generalized[*dof];
        }
      }
      maxDisplacementSquared = std::max(maxDisplacementSquared, double3::dot(displacement, displacement));
    }
    for (std::size_t moleculeIndex = 0; moleculeIndex < referenceSystem.moleculeData.size(); ++moleculeIndex)
    {
      const Molecule& molecule = referenceSystem.moleculeData[moleculeIndex];
      if (layout.molecules()[moleculeIndex].rigid)
      {
        const auto comBase = layout.rigidMoleculeDof(moleculeIndex, RigidDof::ComX);
        const auto orientationBase = layout.rigidMoleculeDof(moleculeIndex, RigidDof::OriX);
        const double3 comDisplacement =
            comBase ? double3(generalized[*comBase + 0], generalized[*comBase + 1], generalized[*comBase + 2])
                    : double3();
        const double3 omega = orientationBase ? double3(generalized[*orientationBase + 0],
                                                        generalized[*orientationBase + 1],
                                                        generalized[*orientationBase + 2])
                                              : double3();
        for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
        {
          const double3 offset = moleculeAtoms[molecule.atomIndex + localAtom].position - molecule.centerOfMassPosition;
          const double3 displacement = comDisplacement + double3::cross(omega, offset);
          maxDisplacementSquared = std::max(maxDisplacementSquared, double3::dot(displacement, displacement));
        }
      }
      else
      {
        for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
        {
          double3 displacement{};
          for (std::size_t axis = 0; axis < 3; ++axis)
          {
            if (const auto dof = layout.flexibleAtomDof(moleculeIndex, localAtom, static_cast<MinimizationDofAxis>(axis)))
            {
              (&displacement.x)[axis] = generalized[*dof];
            }
          }
          maxDisplacementSquared = std::max(maxDisplacementSquared, double3::dot(displacement, displacement));
        }
      }
    }

    const double maxDisplacement = std::sqrt(maxDisplacementSquared);
    const double scale = maxDisplacement > 1.0e-12 ? amplitude / maxDisplacement : 0.0;

    std::ofstream stream(directory / std::format("mode_{:04d}.s{}.pdb", mode, systemIndex));
    std::print(stream, "REMARK normal mode {} omega^2={:.10e} [{}] frequency={:.6f} [{}]\n", mode,
               result.eigenvalues[mode], Units::unitSystem == Units::System::ReducedUnits ? "reduced" : "ps^-2",
               frequencies[mode], Units::unitSystem == Units::System::ReducedUnits ? "reduced" : "cm^-1");

    std::vector<double> displacement(numberOfDofs);
    for (std::size_t frame = 0; frame < totalFrames; ++frame)
    {
      const double theta = 2.0 * std::numbers::pi * static_cast<double>(frame) / static_cast<double>(pointsPerPeriod);
      const double factor = scale * std::sin(theta);
      for (std::size_t dof = 0; dof < numberOfDofs; ++dof)
      {
        displacement[dof] = factor * generalized[dof];
      }

      System frameSystem = referenceSystem;
      applyGeneralizedDisplacement(frameSystem, layout, displacement);
      writePdbFrame(stream, frameSystem, frame + 1);
    }
  }
}
