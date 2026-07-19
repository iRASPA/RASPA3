module;

module minimization_rigid_kinematics;

import std;
import component;
import framework;

namespace
{
// Space-frame offset of the atom from the molecule's center of mass.
double3 atomOffsetFromCenterOfMass(const Molecule &molecule, const Component &component, std::size_t localAtom)
{
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  return rotation * component.atoms[localAtom].position;
}

// A quaternion-tangent perturbation about lab axis e_k rotates the space-frame offset s,
// so dp/domega_k = e_k x s. The second derivative of the exponential map at omega = 0 is
// the symmetrized product (1/2)(E_alpha E_beta + E_beta E_alpha) s with E_k = [e_k]x,
// giving (1/2)(e_alpha (e_beta . s) + e_beta (e_alpha . s)) - delta(alpha,beta) s.
void fillDerivativesFromOffset(Minimization::RigidAtomDerivatives &derivatives, const double3 &s)
{
  derivatives.dVecX = double3(0.0, -s.z, s.y);
  derivatives.dVecY = double3(s.z, 0.0, -s.x);
  derivatives.dVecZ = double3(-s.y, s.x, 0.0);

  derivatives.ddVecAX = double3(0.0, -s.y, -s.z);
  derivatives.ddVecBY = double3(-s.x, 0.0, -s.z);
  derivatives.ddVecCZ = double3(-s.x, -s.y, 0.0);
  derivatives.ddVecAY = double3(0.5 * s.y, 0.5 * s.x, 0.0);
  derivatives.ddVecAZ = double3(0.5 * s.z, 0.0, 0.5 * s.x);
  derivatives.ddVecBZ = double3(0.0, 0.5 * s.z, 0.5 * s.y);
}
}  // namespace

bool Minimization::RigidDerivativeCache::moleculeIsRigid(std::size_t moleculeIndex) const noexcept
{
  return moleculeIndex < _molecules.size() && !_molecules[moleculeIndex].empty();
}

const Minimization::RigidAtomDerivatives &Minimization::RigidDerivativeCache::atom(std::size_t moleculeIndex,
                                                                                   std::size_t localAtomIndex) const noexcept
{
  return _molecules[moleculeIndex][localAtomIndex];
}

const double3 &Minimization::RigidDerivativeCache::bodyCenterOfMass(std::size_t moleculeIndex,
                                                                    std::size_t localAtomIndex) const noexcept
{
  return _bodyCentersOfMass[moleculeIndex][localAtomIndex];
}

bool Minimization::RigidDerivativeCache::frameworkAtomIsRigid(std::size_t frameworkAtom) const noexcept
{
  return frameworkAtom < _frameworkAtomRigid.size() && _frameworkAtomRigid[frameworkAtom] != 0;
}

const Minimization::RigidAtomDerivatives &Minimization::RigidDerivativeCache::frameworkAtom(
    std::size_t frameworkAtom) const noexcept
{
  return _framework[frameworkAtom];
}

const double3 &Minimization::RigidDerivativeCache::frameworkBodyCenterOfMass(std::size_t frameworkAtom) const noexcept
{
  return _frameworkBodyCentersOfMass[frameworkAtom];
}

Minimization::RigidDerivativeCache Minimization::RigidDerivativeCache::build(std::span<const Molecule> moleculeData,
                                                                             std::span<const Component> components,
                                                                             std::span<const Atom> moleculeAtoms,
                                                                             const Framework *framework,
                                                                             std::span<const Atom> frameworkAtoms)
{
  RigidDerivativeCache cache{};
  cache._molecules.resize(moleculeData.size());
  cache._bodyCentersOfMass.resize(moleculeData.size());

  // Framework rigid groups (mixed Fixed/Rigid/Flexible host): the space-frame offsets follow directly
  // from the current laboratory positions relative to each group's mass-weighted center of mass.
  if (framework && framework->isMixed() && !frameworkAtoms.empty())
  {
    cache._framework.resize(frameworkAtoms.size());
    cache._frameworkBodyCentersOfMass.assign(frameworkAtoms.size(), double3{});
    cache._frameworkAtomRigid.assign(frameworkAtoms.size(), 0);
    for (const FrameworkGroup &group : framework->groups)
    {
      if (!group.isRigidBody()) continue;
      double3 centerOfMass{};
      double totalMass{};
      for (std::size_t k = 0; k != group.atoms.size(); ++k)
      {
        centerOfMass += group.atomMasses[k] * frameworkAtoms[group.atoms[k]].position;
        totalMass += group.atomMasses[k];
      }
      centerOfMass = centerOfMass / totalMass;
      for (const std::size_t atom : group.atoms)
      {
        cache._frameworkBodyCentersOfMass[atom] = centerOfMass;
        cache._frameworkAtomRigid[atom] = 1;
        fillDerivativesFromOffset(cache._framework[atom], frameworkAtoms[atom].position - centerOfMass);
      }
    }
  }

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    const Component &component = components[molecule.componentId];

    if (component.rigid)
    {
      cache._molecules[moleculeIndex].resize(molecule.numberOfAtoms);
      cache._bodyCentersOfMass[moleculeIndex].assign(molecule.numberOfAtoms, molecule.centerOfMassPosition);
      for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
      {
        const double3 s = atomOffsetFromCenterOfMass(molecule, component, localAtom);
        fillDerivativesFromOffset(cache._molecules[moleculeIndex][localAtom], s);
      }
    }
    else if (component.isSemiFlexible())
    {
      // Rigid groups of a semi-flexible molecule: the space-frame offsets follow directly from
      // the current laboratory positions relative to the group's mass-weighted center of mass.
      cache._molecules[moleculeIndex].resize(molecule.numberOfAtoms);
      cache._bodyCentersOfMass[moleculeIndex].assign(molecule.numberOfAtoms, double3{});
      const std::span<const Atom> atoms = moleculeAtoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
      for (const Fragment &group : component.fragmentGraph.fragments)
      {
        if (!group.isRigidBody())
        {
          continue;
        }
        double3 centerOfMass{};
        for (const std::size_t localAtom : group.atoms)
        {
          centerOfMass += component.definedAtoms[localAtom].second * atoms[localAtom].position;
        }
        centerOfMass = centerOfMass / group.mass;
        for (const std::size_t localAtom : group.atoms)
        {
          cache._bodyCentersOfMass[moleculeIndex][localAtom] = centerOfMass;
          fillDerivativesFromOffset(cache._molecules[moleculeIndex][localAtom],
                                    atoms[localAtom].position - centerOfMass);
        }
      }
    }
  }

  return cache;
}
