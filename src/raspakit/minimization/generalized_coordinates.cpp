module;

module minimization_generalized_coordinates;

import std;

import double3;
import double3x3;
import simd_quatd;
import atom;
import molecule;
import component;
import simulationbox;
import minimization_cell_layout;

namespace
{
simd_quatd quaternionFromRotationVector(const double3& omega)
{
  const double angle = std::sqrt(double3::dot(omega, omega));
  if (angle < 1.0e-30)
  {
    return simd_quatd(0.0, 0.0, 0.0, 1.0);
  }
  return simd_quatd::fromAxisAngle(angle, omega / angle);
}

void regenerateRigidAtoms(System& system, std::size_t moleculeIndex)
{
  Molecule& molecule = system.moleculeData[moleculeIndex];
  const Component& component = system.components[molecule.componentId];
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
  {
    atoms[molecule.atomIndex + localAtom].position =
        molecule.centerOfMassPosition + rotation * component.atoms[localAtom].position;
  }
}

void wrapMoleculeToPrimaryCell(System& system, std::size_t moleculeIndex)
{
  Molecule& molecule = system.moleculeData[moleculeIndex];
  const Component& component = system.components[molecule.componentId];
  std::span<Atom> atoms = system.spanOfMoleculeAtoms().subspan(molecule.atomIndex, molecule.numberOfAtoms);
  const double3 centerOfMass = component.rigid ? molecule.centerOfMassPosition : component.computeCenterOfMass(atoms);
  const double3 centerOfMassPbc = system.simulationBox.mapToBox(centerOfMass);
  const double3 shift = centerOfMassPbc - centerOfMass;
  if (double3::dot(shift, shift) < 1.0e-30)
  {
    return;
  }

  if (component.rigid)
  {
    molecule.centerOfMassPosition = centerOfMassPbc;
  }
  for (Atom& atom : atoms)
  {
    atom.position += shift;
  }
}

double3x3 matrixExponential(double3x3 matrix)
{
  const double norm = std::max({std::abs(matrix.m11) + std::abs(matrix.m12) + std::abs(matrix.m13),
                                std::abs(matrix.m21) + std::abs(matrix.m22) + std::abs(matrix.m23),
                                std::abs(matrix.m31) + std::abs(matrix.m32) + std::abs(matrix.m33)});
  const int scaling = norm > 0.5 ? std::max(0, static_cast<int>(std::ceil(std::log2(norm / 0.5)))) : 0;
  matrix = matrix * std::ldexp(1.0, -scaling);

  double3x3 result = double3x3::identity();
  double3x3 term = double3x3::identity();
  for (std::size_t order = 1; order <= 24; ++order)
  {
    term = term * matrix;
    term = term * (1.0 / static_cast<double>(order));
    result += term;
  }
  for (int i = 0; i < scaling; ++i)
  {
    result = result * result;
  }
  return result;
}

void applyCellDisplacement(System& system, const MinimizationDofLayout& layout, std::span<const double> displacement)
{
  if (layout.numberOfCellDofs() == 0) return;

  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
  if (cellLayout.size() != layout.numberOfCellDofs())
  {
    throw std::invalid_argument("applyGeneralizedDisplacement: cell layout does not match DOF layout");
  }

  double3x3 logarithmicStrain{};
  for (std::size_t coordinate = 0; coordinate < cellLayout.size(); ++coordinate)
  {
    logarithmicStrain += displacement[*layout.cellDof(coordinate)] * cellLayout.bases[coordinate];
  }
  const double3x3 deformation = matrixExponential(logarithmicStrain);

  for (Atom& atom : system.spanOfFrameworkAtoms())
  {
    atom.position = deformation * atom.position;
  }
  std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  for (Molecule& molecule : system.moleculeData)
  {
    const Component& component = system.components[molecule.componentId];
    if (component.rigid)
    {
      molecule.centerOfMassPosition = deformation * molecule.centerOfMassPosition;
    }
    else
    {
      for (Atom& atom : moleculeAtoms.subspan(molecule.atomIndex, molecule.numberOfAtoms))
      {
        atom.position = deformation * atom.position;
      }
    }
  }

  system.simulationBox = SimulationBox(deformation * system.simulationBox.cell);
}
}  // namespace

void applyGeneralizedDisplacement(System& system, const MinimizationDofLayout& layout,
                                  std::span<const double> displacement)
{
  if (displacement.size() != layout.numDofs())
  {
    throw std::invalid_argument("applyGeneralizedDisplacement: displacement extent does not match DOF layout");
  }

  std::span<Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  for (std::size_t atom = 0; atom < layout.numberOfFrameworkAtoms(); ++atom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      if (const auto dof = layout.frameworkAtomDof(atom, static_cast<MinimizationDofAxis>(axis)))
      {
        (&frameworkAtoms[atom].position.x)[axis] += displacement[*dof];
      }
    }
  }

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    Molecule& molecule = system.moleculeData[moleculeIndex];
    if (!layout.molecules()[moleculeIndex].rigid)
    {
      for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
      {
        double3 delta{};
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
          const auto dof = layout.flexibleAtomDof(moleculeIndex, localAtom, static_cast<MinimizationDofAxis>(axis));
          if (dof)
          {
            (&delta.x)[axis] = displacement[*dof];
          }
        }
        atoms[molecule.atomIndex + localAtom].position += delta;
      }
      wrapMoleculeToPrimaryCell(system, moleculeIndex);
      continue;
    }

    const auto comBase = layout.rigidMoleculeDof(moleculeIndex, RigidDof::ComX);
    const auto orientationBase = layout.rigidMoleculeDof(moleculeIndex, RigidDof::OriX);
    if (comBase)
    {
      molecule.centerOfMassPosition +=
          double3(displacement[*comBase + 0], displacement[*comBase + 1], displacement[*comBase + 2]);
    }
    if (orientationBase)
    {
      const double3 omega(displacement[*orientationBase + 0], displacement[*orientationBase + 1],
                          displacement[*orientationBase + 2]);
      molecule.orientation = (quaternionFromRotationVector(omega) * molecule.orientation).normalized();
    }
    regenerateRigidAtoms(system, moleculeIndex);
    wrapMoleculeToPrimaryCell(system, moleculeIndex);
  }

  if (layout.numberOfCellDofs() != 0)
  {
    applyCellDisplacement(system, layout, displacement);
    for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
    {
      if (layout.molecules()[moleculeIndex].rigid)
      {
        regenerateRigidAtoms(system, moleculeIndex);
      }
      wrapMoleculeToPrimaryCell(system, moleculeIndex);
    }
    // Keep cutoff, Ewald alpha, and reciprocal integer bounds fixed during a minimization.
    // Re-selecting numerical parameters would make the objective discontinuous as the cell changes.
    system.precomputeTotalRigidEnergy();
  }
}
