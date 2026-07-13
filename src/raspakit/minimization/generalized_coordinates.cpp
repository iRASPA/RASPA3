module;

module minimization_generalized_coordinates;

import std;

import double3;
import double3x3;
import simd_quatd;
import atom;
import molecule;
import component;

namespace
{
simd_quatd quaternionFromRotationVector(const double3 &omega)
{
  const double angle = std::sqrt(double3::dot(omega, omega));
  if (angle < 1.0e-30)
  {
    return simd_quatd(0.0, 0.0, 0.0, 1.0);
  }
  return simd_quatd::fromAxisAngle(angle, omega / angle);
}

void regenerateRigidAtoms(System &system, std::size_t moleculeIndex)
{
  Molecule &molecule = system.moleculeData[moleculeIndex];
  const Component &component = system.components[molecule.componentId];
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
  {
    atoms[molecule.atomIndex + localAtom].position =
        molecule.centerOfMassPosition + rotation * component.atoms[localAtom].position;
  }
}
}  // namespace

void applyGeneralizedDisplacement(System &system, const MinimizationDofLayout &layout,
                                  std::span<const double> displacement)
{
  if (displacement.size() != layout.numDofs())
  {
    throw std::invalid_argument("applyGeneralizedDisplacement: displacement extent does not match DOF layout");
  }

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    Molecule &molecule = system.moleculeData[moleculeIndex];
    if (!layout.molecules()[moleculeIndex].rigid)
    {
      for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
      {
        double3 delta{};
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
          const auto dof =
              layout.flexibleAtomDof(moleculeIndex, localAtom, static_cast<MinimizationDofAxis>(axis));
          if (dof)
          {
            (&delta.x)[axis] = displacement[*dof];
          }
        }
        atoms[molecule.atomIndex + localAtom].position += delta;
      }
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
  }
}
