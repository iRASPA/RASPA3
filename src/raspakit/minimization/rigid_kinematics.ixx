module;

export module minimization_rigid_kinematics;

import std;

import double3;
import double3x3;
import simd_quatd;
import atom;
import molecule;
import component;

export namespace Minimization
{
struct RigidAtomDerivatives
{
  double3 dVecX{};
  double3 dVecY{};
  double3 dVecZ{};
  double3 ddVecAX{};
  double3 ddVecBY{};
  double3 ddVecCZ{};
  double3 ddVecAY{};
  double3 ddVecAZ{};
  double3 ddVecBZ{};
};

/**
 * Precomputed orientation Jacobians for rigid molecules (RASPA2 DVec / DDVec tables).
 */
class RigidDerivativeCache
{
 public:
  bool moleculeIsRigid(std::size_t moleculeIndex) const noexcept;

  const RigidAtomDerivatives &atom(std::size_t moleculeIndex, std::size_t localAtomIndex) const noexcept;

  static RigidDerivativeCache build(std::span<const Molecule> moleculeData, std::span<const Component> components,
                                    std::span<const Atom> atoms);

 private:
  std::vector<std::vector<RigidAtomDerivatives>> _molecules;
};

}  // namespace Minimization
