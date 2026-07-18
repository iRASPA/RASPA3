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
  const bool authoritativeGroups = !system.groupData.empty();
  std::span<GroupState> groupData = system.spanOfGroupData();
  std::size_t globalRigidGroupIndex{};
  for (Molecule& molecule : system.moleculeData)
  {
    const Component& component = system.components[molecule.componentId];
    if (component.rigid)
    {
      molecule.centerOfMassPosition = deformation * molecule.centerOfMassPosition;
    }
    else if (component.isSemiFlexible())
    {
      // The cell drives each rigid group's center of mass (the internal group geometry does not
      // scale) and the flexible atoms individually.
      std::span<Atom> span = moleculeAtoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
      std::size_t rigidRank{};
      for (std::size_t groupIndex = 0; groupIndex < component.groups.size(); ++groupIndex)
      {
        if (!component.groups[groupIndex].rigid)
        {
          continue;
        }
        const std::size_t rigidSlot = rigidRank++;
        GroupState localState;
        if (!authoritativeGroups)
        {
          localState = component.deriveGroupState(groupIndex, span);
        }
        GroupState& state = authoritativeGroups ? groupData[globalRigidGroupIndex + rigidSlot] : localState;
        state.centerOfMassPosition = deformation * state.centerOfMassPosition;
        component.regenerateGroupAtoms(state, groupIndex, span);
      }
      for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
      {
        if (!component.rigidGroupContaining(localAtom).has_value())
        {
          span[localAtom].position = deformation * span[localAtom].position;
        }
      }
      globalRigidGroupIndex += component.numberOfRigidGroups();
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

  // Authoritative rigid-group state: when 'System::groupData' is populated (minimization keeps it so
  // via 'System::initializeGroupData'), each rigid group's orientation quaternion is carried across
  // steps instead of being re-fitted from the atom positions with an orthogonal Procrustes solve on
  // every displacement. This is both cheaper and free of the small chart noise the repeated fit
  // introduces. Callers that do not maintain 'groupData' (empty span) still get the exact same result
  // through the derive-from-atoms fallback below.
  const bool authoritativeGroups = !system.groupData.empty();
  std::span<GroupState> groupData = system.spanOfGroupData();
  std::size_t globalRigidGroupIndex{};

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    Molecule& molecule = system.moleculeData[moleculeIndex];
    const Component& component = system.components[molecule.componentId];
    if (component.isSemiFlexible())
    {
      std::span<Atom> moleculeAtoms = atoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
      // Flexible atoms move along their Cartesian degrees of freedom.
      for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
      {
        double3 delta{};
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
          if (const auto dof = layout.flexibleAtomDof(moleculeIndex, localAtom, static_cast<MinimizationDofAxis>(axis)))
          {
            (&delta.x)[axis] = displacement[*dof];
          }
        }
        moleculeAtoms[localAtom].position += delta;
      }
      // Each rigid group advances its center of mass and orientation, then rebuilds its atoms.
      std::size_t rigidRank{};
      for (std::size_t groupIndex = 0; groupIndex < component.groups.size(); ++groupIndex)
      {
        if (!component.groups[groupIndex].rigid)
        {
          continue;
        }
        const std::size_t rigidSlot = rigidRank++;
        const auto comBase = layout.atomRigidComDof(moleculeIndex, component.groups[groupIndex].atoms.front());
        if (!comBase)
        {
          continue;
        }
        // Advance the authoritative group state in place, or derive it from the atoms as a fallback.
        GroupState localState;
        if (!authoritativeGroups)
        {
          localState = component.deriveGroupState(groupIndex, moleculeAtoms);
        }
        GroupState& state = authoritativeGroups ? groupData[globalRigidGroupIndex + rigidSlot] : localState;
        state.centerOfMassPosition +=
            double3(displacement[*comBase + 0], displacement[*comBase + 1], displacement[*comBase + 2]);
        const double3 omega(displacement[*comBase + 3], displacement[*comBase + 4], displacement[*comBase + 5]);
        state.orientation = (quaternionFromRotationVector(omega) * state.orientation).normalized();
        component.regenerateGroupAtoms(state, groupIndex, moleculeAtoms);
      }
      // Wrap the molecule into the primary cell. In authoritative mode the same rigid shift must be
      // folded into the stored group centers of mass so the next regeneration does not undo it.
      const double3 com = component.computeCenterOfMass(moleculeAtoms);
      const double3 comPbc = system.simulationBox.mapToBox(com);
      const double3 shift = comPbc - com;
      if (double3::dot(shift, shift) >= 1.0e-30)
      {
        for (Atom& atom : moleculeAtoms) atom.position += shift;
        if (authoritativeGroups)
        {
          for (std::size_t r = 0; r < component.numberOfRigidGroups(); ++r)
          {
            groupData[globalRigidGroupIndex + r].centerOfMassPosition += shift;
          }
        }
      }
      globalRigidGroupIndex += component.numberOfRigidGroups();
      continue;
    }
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
    std::size_t wrapRigidGroupIndex{};
    for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
    {
      if (layout.molecules()[moleculeIndex].rigid)
      {
        regenerateRigidAtoms(system, moleculeIndex);
      }
      const Component& component = system.components[system.moleculeData[moleculeIndex].componentId];
      if (component.isSemiFlexible())
      {
        // Wrap group-aware: the rigid shift must also be folded into the stored group centers of
        // mass so the next regeneration does not undo it.
        Molecule& molecule = system.moleculeData[moleculeIndex];
        std::span<Atom> moleculeAtoms =
            system.spanOfMoleculeAtoms().subspan(molecule.atomIndex, molecule.numberOfAtoms);
        const double3 com = component.computeCenterOfMass(moleculeAtoms);
        const double3 shift = system.simulationBox.mapToBox(com) - com;
        if (double3::dot(shift, shift) >= 1.0e-30)
        {
          for (Atom& atom : moleculeAtoms) atom.position += shift;
          if (authoritativeGroups)
          {
            for (std::size_t r = 0; r < component.numberOfRigidGroups(); ++r)
            {
              groupData[wrapRigidGroupIndex + r].centerOfMassPosition += shift;
            }
          }
        }
        wrapRigidGroupIndex += component.numberOfRigidGroups();
        continue;
      }
      wrapMoleculeToPrimaryCell(system, moleculeIndex);
    }
    // Keep cutoff, Ewald alpha, and reciprocal integer bounds fixed during a minimization.
    // Re-selecting numerical parameters would make the objective discontinuous as the cell changes.
    system.precomputeTotalRigidEnergy();
  }
}
