module;

module system;

import <cstddef>;

import double3;

import atom;
import component;

// system_reservoir.cpp


// TODO: add equilibration for flexible molecules

std::vector<Atom> System::equilibratedMoleculeRandomInBox(RandomNumber &random, size_t selectedComponent, std::span<Atom> molecule, double scaling, size_t moleculeId) const
{
  size_t startingBead = components[selectedComponent].startingBead;
  double3 center = molecule[startingBead].position;
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());

  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  double3 position = simulationBox.randomPosition(random);

  for (size_t i = 0; i != copied_atoms.size(); ++i)
  {
    copied_atoms[i].setScaling(scaling);
    copied_atoms[i].position = position + randomRotationMatrix * (molecule[i].position - center);
    copied_atoms[i].moleculeId = static_cast<uint32_t>(moleculeId);
  }
  return copied_atoms;
}
