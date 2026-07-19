module;

export module minimization_rigid_kinematics;

import std;

import double3;
import double3x3;
import simd_quatd;
import atom;
import molecule;
import component;
import framework;

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
 * Precomputed orientation Jacobians for rigid bodies (RASPA2 DVec / DDVec tables). A rigid body is
 * either a whole rigid molecule or a rigid group inside a semi-flexible molecule; the tables are
 * indexed per atom and are populated only for atoms driven by a rigid body.
 */
class RigidDerivativeCache
{
 public:
  bool moleculeIsRigid(std::size_t moleculeIndex) const noexcept;

  const RigidAtomDerivatives &atom(std::size_t moleculeIndex, std::size_t localAtomIndex) const noexcept;

  /** Laboratory-frame center of mass of the rigid body driving this atom (molecule or group). */
  const double3 &bodyCenterOfMass(std::size_t moleculeIndex, std::size_t localAtomIndex) const noexcept;

  /** True when the (global) framework atom is driven by a mixed-framework rigid group. */
  bool frameworkAtomIsRigid(std::size_t frameworkAtom) const noexcept;

  /** Orientation Jacobian of the framework rigid group driving this (global) framework atom. */
  const RigidAtomDerivatives &frameworkAtom(std::size_t frameworkAtom) const noexcept;

  /** Laboratory-frame center of mass of the framework rigid group driving this atom. */
  const double3 &frameworkBodyCenterOfMass(std::size_t frameworkAtom) const noexcept;

  static RigidDerivativeCache build(std::span<const Molecule> moleculeData, std::span<const Component> components,
                                    std::span<const Atom> moleculeAtoms, const Framework *framework = nullptr,
                                    std::span<const Atom> frameworkAtoms = {});

 private:
  std::vector<std::vector<RigidAtomDerivatives>> _molecules;
  std::vector<std::vector<double3>> _bodyCentersOfMass;
  std::vector<RigidAtomDerivatives> _framework;        ///< Per (global) framework atom; only rigid-group atoms set.
  std::vector<double3> _frameworkBodyCentersOfMass;    ///< Group center of mass per rigid-group framework atom.
  std::vector<std::uint8_t> _frameworkAtomRigid;       ///< 1 when the framework atom is driven by a rigid group.
};

}  // namespace Minimization
