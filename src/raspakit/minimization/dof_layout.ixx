module;

export module minimization_dof_layout;

import std;

import molecule;
import component;
import framework;

export enum class MinimizationDofAxis : std::size_t {
  X = 0,
  Y = 1,
  Z = 2,
};

export enum class RigidDof : std::size_t {
  ComX = 0,
  ComY = 1,
  ComZ = 2,
  OriX = 3,
  OriY = 4,
  OriZ = 5,
};

export struct MolDofInfo
{
  bool rigid{false};
  std::size_t base{};
  std::size_t nAtoms{};
};

/**
 * A rigid body carrying six degrees of freedom (center of mass + orientation tangent). It is
 * either a whole rigid molecule or a rigid group inside a semi-flexible molecule.
 */
export struct RigidBodyDofInfo
{
  std::size_t moleculeIndex{};
  std::size_t groupIndex{};   ///< Rigid-group index within the component (unused for whole molecules).
  bool wholeMolecule{false};  ///< True for a fully rigid molecule, false for a rigid group.
  bool isFramework{false};    ///< True for a rigid group of a mixed framework.
  std::size_t base{};         ///< Center-of-mass DOF base; orientation base is 'base + 3'.
};

/**
 * Unordered minimization degree-of-freedom layout.
 *
 * Molecules keep simulation order. Each flexible atom contributes three Cartesian DOFs;
 * each rigid body (a rigid molecule, or a rigid group inside a semi-flexible molecule)
 * contributes six (center of mass and orientation tangent components).
 */
export class MinimizationDofLayout
{
 public:
  std::size_t numDofs() const noexcept { return _numDofs; }
  std::size_t numberOfPositionDofs() const noexcept { return _numberOfPositionDofs; }
  std::size_t numberOfCellDofs() const noexcept { return _numberOfCellDofs; }
  std::optional<std::size_t> cellDof(std::size_t cellCoordinate) const noexcept;

  std::span<const MolDofInfo> molecules() const noexcept { return _molecules; }

  std::size_t numberOfFrameworkAtoms() const noexcept { return _numberOfFrameworkAtoms; }
  std::optional<std::size_t> frameworkAtomDof(std::size_t atom, MinimizationDofAxis axis) const noexcept;

  /** COM DOF base of the framework rigid group driving this atom, if any. */
  std::optional<std::size_t> frameworkAtomRigidComDof(std::size_t atom) const noexcept;

  std::optional<std::size_t> flexibleAtomDof(std::size_t moleculeIndex, std::size_t localAtom,
                                             MinimizationDofAxis axis) const noexcept;

  std::optional<std::size_t> rigidMoleculeDof(std::size_t moleculeIndex, RigidDof dof) const noexcept;

  /**
   * Center-of-mass DOF base of the rigid body driving this atom (whole rigid molecule or rigid
   * group of a semi-flexible molecule); nullopt when the atom carries Cartesian DOFs instead.
   */
  std::optional<std::size_t> atomRigidComDof(std::size_t moleculeIndex, std::size_t localAtom) const noexcept;

  /** Orientation DOF base of the rigid body driving this atom (com base + 3), if any. */
  std::optional<std::size_t> atomRigidOrientationDof(std::size_t moleculeIndex, std::size_t localAtom) const noexcept;

  /** All rigid bodies (rigid molecules and rigid groups) in DOF order. */
  std::span<const RigidBodyDofInfo> rigidBodies() const noexcept { return _rigidBodies; }

  std::size_t flexibleAtomDofBase(std::size_t moleculeIndex, std::size_t localAtom) const;

  std::size_t rigidMoleculeDofBase(std::size_t moleculeIndex) const;

 private:
  friend MinimizationDofLayout buildMinimizationDofLayout(std::span<const Molecule> moleculeData,
                                                          std::span<const Component> components,
                                                          std::size_t numberOfFrameworkAtoms,
                                                          std::size_t numberOfCellDofs,
                                                          const Framework* framework);

  std::size_t _numDofs{};
  std::size_t _numberOfPositionDofs{};
  std::size_t _numberOfCellDofs{};
  std::size_t _numberOfFrameworkAtoms{};
  std::size_t _maxAtomsPerMolecule{};
  std::vector<std::int32_t> _frameworkAtomDof;
  std::vector<std::int32_t> _frameworkAtomRigidBodyDof;
  std::vector<MolDofInfo> _molecules;
  std::vector<std::int32_t> _flexibleAtomDof;
  std::vector<std::int32_t> _rigidMoleculeDof;
  std::vector<std::int32_t> _atomRigidBodyDof;
  std::vector<RigidBodyDofInfo> _rigidBodies;
};

export MinimizationDofLayout buildMinimizationDofLayout(std::span<const Molecule> moleculeData,
                                                        std::span<const Component> components,
                                                        std::size_t numberOfFrameworkAtoms = 0,
                                                        std::size_t numberOfCellDofs = 0,
                                                        const Framework* framework = nullptr);
