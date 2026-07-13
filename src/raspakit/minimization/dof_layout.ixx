module;

export module minimization_dof_layout;

import std;

import molecule;
import component;

export enum class MinimizationDofAxis : std::size_t
{
  X = 0,
  Y = 1,
  Z = 2,
};

export enum class RigidDof : std::size_t
{
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
 * Unordered minimization degree-of-freedom layout.
 *
 * Molecules keep simulation order. Each flexible atom contributes three Cartesian DOFs;
 * each rigid molecule contributes six (center of mass and orientation tangent components).
 */
export class MinimizationDofLayout
{
 public:
  std::size_t numDofs() const noexcept { return _numDofs; }

  std::span<const MolDofInfo> molecules() const noexcept { return _molecules; }

  std::optional<std::size_t> flexibleAtomDof(std::size_t moleculeIndex, std::size_t localAtom,
                                             MinimizationDofAxis axis) const noexcept;

  std::optional<std::size_t> rigidMoleculeDof(std::size_t moleculeIndex, RigidDof dof) const noexcept;

  std::size_t flexibleAtomDofBase(std::size_t moleculeIndex, std::size_t localAtom) const;

  std::size_t rigidMoleculeDofBase(std::size_t moleculeIndex) const;

 private:
  friend MinimizationDofLayout buildMinimizationDofLayout(std::span<const Molecule> moleculeData,
                                                          std::span<const Component> components);

  std::size_t _numDofs{};
  std::size_t _maxAtomsPerMolecule{};
  std::vector<MolDofInfo> _molecules;
  std::vector<std::int32_t> _flexibleAtomDof;
  std::vector<std::int32_t> _rigidMoleculeDof;
};

export MinimizationDofLayout buildMinimizationDofLayout(std::span<const Molecule> moleculeData,
                                                        std::span<const Component> components);
