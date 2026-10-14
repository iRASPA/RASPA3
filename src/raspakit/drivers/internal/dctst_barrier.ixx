module;

export module dctst_barrier;

import std;

import int3;
import double3;
import double3x3;
import simulationbox;
import randomnumbers;

/**
 * \brief The dividing surface of a dcTST calculation and the cage lattice it generates.
 *
 * The barrier is the plane through 'point' with unit normal 'normal'. Its periodic images at
 * q = k lambda (k integer) along the normal are the dividing surfaces between consecutive cages;
 * the reaction coordinate of a bead at r is q = (r - point) . normal and its cage index is
 * floor(q / lambda). The barrier plane itself separates cage A (index -1, q < 0) from cage B
 * (index 0, q > 0).
 *
 * All projections apply the minimum-image convention to the vector from 'point' to the bead, so
 * wrapped and unwrapped positions give the same result as long as the box is at least two cages
 * wide along the normal. Net cage changes along a trajectory are accumulated from the
 * minimum-image displacements between consecutive integration steps ('displacement').
 */
export struct DcTSTBarrier
{
  DcTSTBarrier() = default;
  DcTSTBarrier(const SimulationBox& simulationBox, double3 fractionalPosition, double3 cartesianNormal,
               double latticeDistance);

  double3 point{};            ///< A point on the barrier plane (Cartesian).
  double3 normal{1, 0, 0};    ///< Unit normal of the barrier plane (Cartesian).
  double3 inPlaneU{0, 1, 0};  ///< First unit vector spanning the plane.
  double3 inPlaneV{0, 0, 1};  ///< Second unit vector spanning the plane.
  double lambda{1.0};         ///< Distance between equivalent barrier planes (the cage-to-cage distance).

  /// Reaction coordinate of a bead at 'r': the minimum-image distance to the barrier plane along
  /// the normal.
  [[nodiscard]] double projection(const SimulationBox& simulationBox, double3 r) const;

  /// The component of the minimum-image displacement from 'rOld' to 'rNew' along the normal.
  [[nodiscard]] double displacement(const SimulationBox& simulationBox, double3 rOld, double3 rNew) const;

  /// Cage index of reaction coordinate 'q': floor(q / lambda); cage B is 0, cage A is -1.
  [[nodiscard]] std::int64_t cageIndex(double q) const;

  /// Distance (in cages) of cage 'cage' from the pair {A, B}: 0 for A and B themselves, 1 for their
  /// neighbors, and so on.
  [[nodiscard]] static std::size_t shellOfCage(std::int64_t cage);

  /// 'q' reduced into the period [-lambda/2, lambda/2) around the barrier.
  [[nodiscard]] double reducedCoordinate(double q) const;

  /// Bin index of 'q' on a histogram of 'numberOfBins' bins over [-lambda/2, lambda/2) with bin 0
  /// centered on the barrier.
  [[nodiscard]] std::size_t bin(double q, std::size_t numberOfBins) const;

  /// Reaction coordinate of the center of bin 'index' of such a histogram.
  [[nodiscard]] double binCenter(std::size_t index, std::size_t numberOfBins) const;

  /// A uniformly distributed random point on the plane q = 'offset' (within the periodic box).
  [[nodiscard]] double3 randomPointOnPlane(RandomNumber& random, const SimulationBox& simulationBox,
                                           double offset = 0.0) const;

  /// Shifts 'r' along the normal so that it lies on the plane q = 'offset'.
  [[nodiscard]] double3 projectOntoPlane(const SimulationBox& simulationBox, double3 r, double offset = 0.0) const;
};

/**
 * \brief The three-dimensional cage lattice used by the occupancy-resolved mode.
 *
 * One cage per unit cell of the framework: cage (i, j, k) is the unit-cell parallelepiped centered on
 * 'origin' + cell (i, j, k), with (i, j, k) taken modulo the number of unit cells of the supercell. A
 * bead belongs to the cage whose center is nearest in fractional coordinates, so that along the barrier
 * normal the cage boundaries coincide with the periodic images of the barrier plane.
 */
export struct DcTSTCageLattice
{
  DcTSTCageLattice() = default;
  DcTSTCageLattice(const SimulationBox& unitCell, int3 numberOfUnitCells, double3 cageCenter);

  double3x3 cell{};         ///< Unit cell (columns are the lattice vectors).
  double3x3 inverseCell{};  ///< Its inverse.
  int3 numberOfUnitCells{1, 1, 1};
  double3 origin{};  ///< Cartesian center of cage (0, 0, 0).

  [[nodiscard]] std::size_t numberOfCages() const;

  /// Linear index of the cage containing 'r' (periodic in the supercell).
  [[nodiscard]] std::size_t cageOf(double3 r) const;

  /// Cartesian center of cage 'cage'.
  [[nodiscard]] double3 center(std::size_t cage) const;

  /// The cage whose center is 'shift' away from the center of 'cage'.
  [[nodiscard]] std::size_t neighbor(std::size_t cage, double3 shift) const;

  /// A uniformly distributed random point inside cage 'cage'.
  [[nodiscard]] double3 randomPointInCage(RandomNumber& random, std::size_t cage) const;

  [[nodiscard]] int3 indices(std::size_t cage) const;
  [[nodiscard]] std::size_t index(int3 n) const;

  /// The neighbor of 'cage' one cell along lattice axis 'axis' (0, 1, 2) in direction 'step' (+1 or -1).
  [[nodiscard]] std::size_t neighborAlongAxis(std::size_t cage, int axis, int step) const;

  /// Fractional coordinates (unit-cell basis) of 'r' relative to the center of 'cage', with the minimum-image
  /// convention of the supercell 'simulationBox'; a bead inside its cage has all components in [-1/2, 1/2).
  [[nodiscard]] double3 fractionalOffset(const SimulationBox& simulationBox, double3 r, std::size_t cage) const;

  /// The lattice axis most aligned with 'direction', the sign of the alignment (+1 if the lattice vector
  /// points along 'direction') and |cos| of the angle between them.
  struct AxisAlignment
  {
    int axis{0};
    int sign{1};
    double cosine{1.0};
  };
  [[nodiscard]] AxisAlignment axisAlongDirection(double3 direction) const;

  /// Length of lattice vector 'axis'.
  [[nodiscard]] double axisLength(int axis) const;
};
