module;

module minimization_mass_metric;

import std;

import double3;
import double3x3;
import atom;
import molecule;
import component;
import forcefield;
import minimization_dof_layout;
import symmetric_eigensolver;

namespace
{
double checkedInverseSqrtMass(double mass, std::string_view what)
{
  if (!(mass > 0.0) || !std::isfinite(mass))
  {
    throw std::runtime_error(
        std::format("Mass-weighting requires positive masses, but {} has mass {}", what, mass));
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
        inverseSqrt[row * 3 + column] +=
            weight * eigensystem.eigenvector(row, mode) * eigensystem.eigenvector(column, mode);
      }
    }
  }
  return inverseSqrt;
}
}  // namespace

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
