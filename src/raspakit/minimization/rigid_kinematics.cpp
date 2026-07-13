module;

module minimization_rigid_kinematics;

import std;
import component;

namespace
{
// Space-frame offset of the atom from the molecule's center of mass.
double3 atomOffsetFromCenterOfMass(const Molecule &molecule, const Component &component, std::size_t localAtom)
{
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  return rotation * component.atoms[localAtom].position;
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

Minimization::RigidDerivativeCache Minimization::RigidDerivativeCache::build(std::span<const Molecule> moleculeData,
                                                                             std::span<const Component> components,
                                                                             std::span<const Atom> atoms)
{
  RigidDerivativeCache cache{};
  cache._molecules.resize(moleculeData.size());

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    const Component &component = components[molecule.componentId];
    if (!component.rigid)
    {
      continue;
    }

    cache._molecules[moleculeIndex].resize(molecule.numberOfAtoms);
    for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
    {
      // A quaternion-tangent perturbation about lab axis e_k rotates the space-frame offset s,
      // so dp/domega_k = e_k x s. The second derivative of the exponential map at omega = 0 is
      // the symmetrized product (1/2)(E_alpha E_beta + E_beta E_alpha) s with E_k = [e_k]x,
      // giving (1/2)(e_alpha (e_beta . s) + e_beta (e_alpha . s)) - delta(alpha,beta) s.
      const double3 s = atomOffsetFromCenterOfMass(molecule, component, localAtom);

      RigidAtomDerivatives &derivatives = cache._molecules[moleculeIndex][localAtom];
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
  }

  (void)atoms;
  return cache;
}
