module;

module system;

import std;

import randomnumbers;
import atom;
import atom_dynamics;
import molecule;
import component;
import double3;
import simd_quatd;

// System layout: spans, iterators, and molecule indexing over contiguous storage.

std::size_t System::randomMoleculeOfComponent(RandomNumber& random, std::size_t selectedComponent)
{
  return std::size_t(random.uniform() * static_cast<double>(numberOfMoleculesPerComponent[selectedComponent]));
}

std::size_t System::randomIntegerMoleculeOfComponent(RandomNumber& random, std::size_t selectedComponent)
{
  return numberOfFractionalMoleculesPerComponent[selectedComponent] +
         std::size_t(random.uniform() * static_cast<double>(numberOfIntegerMoleculesPerComponent[selectedComponent]));
}

std::vector<Atom>::iterator System::iteratorForMolecule(std::size_t selectedComponent, std::size_t selectedMolecule)
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    std::size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  std::size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule + numberOfFrameworkAtoms;
  return atomData.begin() + static_cast<std::vector<Atom>::difference_type>(index);
}

std::vector<double3>::iterator System::iteratorForElectricField(std::size_t selectedComponent,
                                                                std::size_t selectedMolecule)
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    std::size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  std::size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule + numberOfFrameworkAtoms;
  return electricField.begin() + static_cast<std::vector<double3>::difference_type>(index);
}

std::vector<Molecule>::iterator System::indexForMolecule(std::size_t selectedComponent, std::size_t selectedMolecule)
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    index += numberOfMoleculesPerComponent[i];
  }
  index += selectedMolecule;
  return moleculeData.begin() + static_cast<std::vector<Atom>::difference_type>(index);
}

std::size_t System::moleculeIndexOfComponent(std::size_t selectedComponent, std::size_t selectedMolecule)
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    index += numberOfMoleculesPerComponent[i];
  }
  index += selectedMolecule;
  return index;
}

std::span<const Atom> System::spanOfFrameworkAtoms() const
{
  return std::span(atomData.begin(), numberOfFrameworkAtoms);
}

std::span<Atom> System::spanOfFrameworkAtoms() { return std::span(atomData.begin(), numberOfFrameworkAtoms); }

std::span<const Atom> System::spanOfRigidFrameworkAtoms() const
{
  return std::span(atomData.begin(), numberOfRigidFrameworkAtoms);
}

std::span<const Atom> System::spanOfFlexibleAtoms() const
{
  return std::span(atomData.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfFrameworkAtoms),
                   atomData.end());
}

std::span<const Atom> System::spanOfMoleculeAtoms() const
{
  return std::span(atomData.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfFrameworkAtoms),
                   atomData.end());
}

std::span<Atom> System::spanOfMoleculeAtoms()
{
  return std::span(atomData.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfFrameworkAtoms),
                   atomData.end());
}

std::span<const AtomDynamics> System::spanOfFrameworkDynamics() const
{
  return std::span(atomDynamics.begin(), numberOfFrameworkAtoms);
}

std::span<AtomDynamics> System::spanOfFrameworkDynamics()
{
  return std::span(atomDynamics.begin(), numberOfFrameworkAtoms);
}

std::span<const AtomDynamics> System::spanOfMoleculeDynamics() const
{
  return std::span(
      atomDynamics.begin() + static_cast<std::vector<AtomDynamics>::difference_type>(numberOfFrameworkAtoms),
      atomDynamics.end());
}

std::span<AtomDynamics> System::spanOfMoleculeDynamics()
{
  return std::span(
      atomDynamics.begin() + static_cast<std::vector<AtomDynamics>::difference_type>(numberOfFrameworkAtoms),
      atomDynamics.end());
}

void System::initializeGroupData()
{
  std::size_t totalRigidGroups{};
  for (const Molecule& molecule : moleculeData)
  {
    totalRigidGroups += components[molecule.componentId].numberOfRigidGroups();
  }
  // When the layout is unchanged (e.g. entering production after equilibration, where the group states
  // are kept aligned across particle exchanges) the velocities and orientation momenta are preserved;
  // only the kinematic state is re-derived from the (authoritative) atom positions. A fresh or
  // resized groupData starts with zero dynamic state.
  const bool preserveDynamics = groupData.size() == totalRigidGroups;
  if (!preserveDynamics)
  {
    groupData.assign(totalRigidGroups, GroupState{});
  }
  if (totalRigidGroups == 0) return;

  std::span<const Atom> atoms = spanOfMoleculeAtoms();
  std::size_t atomIndex{};
  std::size_t groupIndex{};
  for (const Molecule& molecule : moleculeData)
  {
    const Component& component = components[molecule.componentId];
    if (component.isSemiFlexible())
    {
      std::span<const Atom> span = atoms.subspan(atomIndex, molecule.numberOfAtoms);
      std::size_t rigidRank{};
      for (std::size_t g = 0; g != component.groups.size(); ++g)
      {
        if (component.groups[g].rigid)
        {
          GroupState& state = groupData[groupIndex + rigidRank];
          GroupState derived = component.deriveGroupState(g, span);
          if (preserveDynamics)
          {
            // q and -q describe the same rotation; keep the sign continuous with the stored
            // orientation because the preserved orientation momentum is paired with it
            const double overlap = derived.orientation.r * state.orientation.r +
                                   derived.orientation.ix * state.orientation.ix +
                                   derived.orientation.iy * state.orientation.iy +
                                   derived.orientation.iz * state.orientation.iz;
            if (overlap < 0.0)
            {
              derived.orientation = simd_quatd(-derived.orientation.ix, -derived.orientation.iy,
                                               -derived.orientation.iz, -derived.orientation.r);
            }
          }
          state.centerOfMassPosition = derived.centerOfMassPosition;
          state.orientation = derived.orientation;
          state.gradient = derived.gradient;
          state.orientationGradient = derived.orientationGradient;
          if (!preserveDynamics)
          {
            state.velocity = derived.velocity;
            state.orientationMomentum = derived.orientationMomentum;
          }
          ++rigidRank;
        }
      }
      groupIndex += component.numberOfRigidGroups();
    }
    atomIndex += molecule.numberOfAtoms;
  }
}

std::span<GroupState> System::spanOfGroupData() { return std::span(groupData); }

std::span<const GroupState> System::spanOfGroupData() const { return std::span(groupData); }

std::span<double> System::spanOfMoleculeElectrostaticPotential()
{
  return std::span(
      electricPotential.begin() + static_cast<std::vector<double3>::difference_type>(numberOfFrameworkAtoms),
      electricPotential.end());
}

std::span<double3> System::spanOfMoleculeElectricField()
{
  return std::span(electricField.begin() + static_cast<std::vector<double3>::difference_type>(numberOfFrameworkAtoms),
                   electricField.end());
}

std::span<double3> System::spanOfMoleculeElectricFieldNew()
{
  return std::span(
      electricFieldNew.begin() + static_cast<std::vector<double3>::difference_type>(numberOfFrameworkAtoms),
      electricFieldNew.end());
}

std::span<Atom> System::spanOfMolecule(std::size_t selectedComponent, std::size_t selectedMolecule)
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    std::size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  std::size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule;
  return std::span(&atomData[index + numberOfFrameworkAtoms], size);
}

const std::span<const Atom> System::spanOfMolecule(std::size_t selectedComponent, std::size_t selectedMolecule) const
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    std::size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  std::size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule;
  return std::span(&atomData[index + numberOfFrameworkAtoms], size);
}

const std::span<const Atom> System::spanOfIntegerAtomsOfComponent(std::size_t selectedComponent) const
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    std::size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  std::size_t size = components[selectedComponent].atoms.size();
  index += size * numberOfFractionalMoleculesPerComponent[selectedComponent];
  std::size_t number_of_atoms = size * (numberOfMoleculesPerComponent[selectedComponent] -
                                        numberOfFractionalMoleculesPerComponent[selectedComponent]);
  return std::span(&atomData[index + numberOfFrameworkAtoms], number_of_atoms);
}

std::span<double3> System::spanElectricFieldNew(std::size_t selectedComponent, std::size_t selectedMolecule)
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    std::size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  std::size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule;
  return std::span(&electricFieldNew[index + numberOfFrameworkAtoms], size);
}

const std::span<const double3> System::spanElectricFieldNew(std::size_t selectedComponent,
                                                            std::size_t selectedMolecule) const
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    std::size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  std::size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule;
  return std::span(&electricFieldNew[index + numberOfFrameworkAtoms], size);
}

std::span<double3> System::spanElectricFieldOld(std::size_t selectedComponent, std::size_t selectedMolecule)
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    std::size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  std::size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule;
  return std::span(&electricField[index + numberOfFrameworkAtoms], size);
}

const std::span<const double3> System::spanElectricFieldOld(std::size_t selectedComponent,
                                                            std::size_t selectedMolecule) const
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    std::size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  std::size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule;
  return std::span(&electricField[index + numberOfFrameworkAtoms], size);
}

std::size_t System::globalIndexOfComponentAndMolecule(std::size_t selectedComponent, std::size_t selectedMolecule)
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    index += numberOfMoleculesPerComponent[i];
  }
  index += selectedMolecule;
  return index;
}


std::size_t System::indexOfFirstMolecule(std::size_t selectedComponent)
{
  std::size_t index{0};
  for (std::size_t i = 0; i < selectedComponent; ++i)
  {
    std::size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  return index + numberOfFrameworkAtoms;
}
