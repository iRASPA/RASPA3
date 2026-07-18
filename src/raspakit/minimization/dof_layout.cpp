module;

module minimization_dof_layout;

import std;

MinimizationDofLayout buildMinimizationDofLayout(std::span<const Molecule> moleculeData,
                                                 std::span<const Component> components,
                                                 std::size_t numberOfFlexibleFrameworkAtoms,
                                                 std::size_t numberOfCellDofs)
{
  MinimizationDofLayout layout;
  layout._numberOfFrameworkAtoms = numberOfFlexibleFrameworkAtoms;
  layout._numberOfCellDofs = numberOfCellDofs;
  layout._molecules.resize(moleculeData.size());

  for (const Molecule& molecule : moleculeData)
  {
    layout._maxAtomsPerMolecule = std::max(layout._maxAtomsPerMolecule, molecule.numberOfAtoms);
  }

  layout._flexibleAtomDof.assign(moleculeData.size() * layout._maxAtomsPerMolecule * 3, -1);
  layout._rigidMoleculeDof.assign(moleculeData.size() * 6, -1);
  layout._atomRigidBodyDof.assign(moleculeData.size() * layout._maxAtomsPerMolecule, -1);
  layout._frameworkAtomDof.assign(numberOfFlexibleFrameworkAtoms * 3, -1);

  std::size_t nextDof{};
  for (std::size_t atom = 0; atom < numberOfFlexibleFrameworkAtoms; ++atom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      layout._frameworkAtomDof[atom * 3 + axis] = static_cast<std::int32_t>(nextDof++);
    }
  }
  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    const Component& component = components[molecule.componentId];
    MolDofInfo& info = layout._molecules[moleculeIndex];
    info.rigid = component.rigid;
    info.nAtoms = molecule.numberOfAtoms;
    info.base = nextDof;

    if (component.rigid)
    {
      const std::size_t base = nextDof;
      for (std::size_t dof = 0; dof < 6; ++dof)
      {
        layout._rigidMoleculeDof[moleculeIndex * 6 + dof] = static_cast<std::int32_t>(nextDof++);
      }
      layout._rigidBodies.push_back(
          RigidBodyDofInfo{.moleculeIndex = moleculeIndex, .groupIndex = 0, .wholeMolecule = true, .base = base});
      for (std::size_t atom = 0; atom < molecule.numberOfAtoms; ++atom)
      {
        layout._atomRigidBodyDof[moleculeIndex * layout._maxAtomsPerMolecule + atom] = static_cast<std::int32_t>(base);
      }
    }
    else if (component.isSemiFlexible())
    {
      // Flexible atoms carry Cartesian DOFs; each rigid group carries a single six-DOF block.
      for (std::size_t atom = 0; atom < molecule.numberOfAtoms; ++atom)
      {
        if (component.rigidGroupContaining(atom).has_value())
        {
          continue;
        }
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
          const std::size_t offset = (moleculeIndex * layout._maxAtomsPerMolecule + atom) * 3 + axis;
          layout._flexibleAtomDof[offset] = static_cast<std::int32_t>(nextDof++);
        }
      }
      for (std::size_t groupIndex = 0; groupIndex < component.groups.size(); ++groupIndex)
      {
        const MoleculeGroup& group = component.groups[groupIndex];
        if (!group.rigid)
        {
          continue;
        }
        const std::size_t base = nextDof;
        nextDof += 6;
        layout._rigidBodies.push_back(RigidBodyDofInfo{
            .moleculeIndex = moleculeIndex, .groupIndex = groupIndex, .wholeMolecule = false, .base = base});
        for (const std::size_t atom : group.atoms)
        {
          layout._atomRigidBodyDof[moleculeIndex * layout._maxAtomsPerMolecule + atom] =
              static_cast<std::int32_t>(base);
        }
      }
    }
    else
    {
      for (std::size_t atom = 0; atom < molecule.numberOfAtoms; ++atom)
      {
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
          const std::size_t offset = (moleculeIndex * layout._maxAtomsPerMolecule + atom) * 3 + axis;
          layout._flexibleAtomDof[offset] = static_cast<std::int32_t>(nextDof++);
        }
      }
    }
  }

  layout._numberOfPositionDofs = nextDof;
  layout._numDofs = nextDof + numberOfCellDofs;
  return layout;
}

std::optional<std::size_t> MinimizationDofLayout::cellDof(std::size_t cellCoordinate) const noexcept
{
  if (cellCoordinate >= _numberOfCellDofs) return std::nullopt;
  return _numberOfPositionDofs + cellCoordinate;
}

std::optional<std::size_t> MinimizationDofLayout::frameworkAtomDof(std::size_t atom,
                                                                   MinimizationDofAxis axis) const noexcept
{
  if (atom >= _numberOfFrameworkAtoms) return std::nullopt;
  const std::int32_t index = _frameworkAtomDof[atom * 3 + static_cast<std::size_t>(axis)];
  if (index < 0) return std::nullopt;
  return static_cast<std::size_t>(index);
}

std::optional<std::size_t> MinimizationDofLayout::flexibleAtomDof(std::size_t moleculeIndex, std::size_t localAtom,
                                                                  MinimizationDofAxis axis) const noexcept
{
  if (moleculeIndex >= _molecules.size() || _molecules[moleculeIndex].rigid ||
      localAtom >= _molecules[moleculeIndex].nAtoms)
  {
    return std::nullopt;
  }

  const std::size_t offset = (moleculeIndex * _maxAtomsPerMolecule + localAtom) * 3 + static_cast<std::size_t>(axis);
  const std::int32_t index = _flexibleAtomDof[offset];
  if (index < 0)
  {
    return std::nullopt;
  }
  return static_cast<std::size_t>(index);
}

std::optional<std::size_t> MinimizationDofLayout::rigidMoleculeDof(std::size_t moleculeIndex,
                                                                   RigidDof dof) const noexcept
{
  if (moleculeIndex >= _molecules.size() || !_molecules[moleculeIndex].rigid)
  {
    return std::nullopt;
  }

  const std::int32_t index = _rigidMoleculeDof[moleculeIndex * 6 + static_cast<std::size_t>(dof)];
  if (index < 0)
  {
    return std::nullopt;
  }
  return static_cast<std::size_t>(index);
}

std::optional<std::size_t> MinimizationDofLayout::atomRigidComDof(std::size_t moleculeIndex,
                                                                  std::size_t localAtom) const noexcept
{
  if (moleculeIndex >= _molecules.size() || localAtom >= _molecules[moleculeIndex].nAtoms)
  {
    return std::nullopt;
  }
  const std::int32_t index = _atomRigidBodyDof[moleculeIndex * _maxAtomsPerMolecule + localAtom];
  if (index < 0)
  {
    return std::nullopt;
  }
  return static_cast<std::size_t>(index);
}

std::optional<std::size_t> MinimizationDofLayout::atomRigidOrientationDof(std::size_t moleculeIndex,
                                                                          std::size_t localAtom) const noexcept
{
  const std::optional<std::size_t> base = atomRigidComDof(moleculeIndex, localAtom);
  if (!base)
  {
    return std::nullopt;
  }
  return *base + 3;
}

std::size_t MinimizationDofLayout::flexibleAtomDofBase(std::size_t moleculeIndex, std::size_t localAtom) const
{
  return *flexibleAtomDof(moleculeIndex, localAtom, MinimizationDofAxis::X);
}

std::size_t MinimizationDofLayout::rigidMoleculeDofBase(std::size_t moleculeIndex) const
{
  return _molecules.at(moleculeIndex).base;
}
