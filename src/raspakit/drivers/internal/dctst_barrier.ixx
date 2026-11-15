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
 * The normal must be a crystallographic direction of the periodic box: normal ~ h1 b1 + h2 b2 + h3 b3
 * with integer (h1, h2, h3) and b_i the reciprocal vectors of the box, so that the barrier and its
 * periodic images are lattice planes s . h = const (s the fractional coordinates). Along the normal the
 * box then repeats every 'period' = 1 / |h1 b1 + h2 b2 + h3 b3| (coprime h), which must be a whole
 * number of lattice distances. The reaction coordinate is the scalar distance from the plane wrapped
 * into [-period/2, period/2): it is continuous across the faces of the box for every such normal
 * (a component-wise minimum image followed by a dot product jumps by O(L) for an oblique normal), and
 * the projection onto the plane maps a uniform point of the box onto a uniform point of the periodic
 * plane (every point of the plane collects the same segment of length 'period'). Net cage changes
 * along a trajectory are accumulated from the minimum-image displacements between consecutive
 * integration steps ('displacement').
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
  int3 millerIndices{1, 0, 0};  ///< Coprime (h1, h2, h3) of the barrier plane in the basis of the box.
  double period{1.0};           ///< Repeat distance of the periodic box along the normal (a multiple of lambda).

  /// Reaction coordinate of a bead at 'r': the signed distance to the barrier plane along the normal,
  /// wrapped into [-period/2, period/2).
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

  /// Whether 'x' lies in the periodic interval of radius 'halfWidth' about 'center' on the reduced
  /// coordinate (the reinsertion support of an unconditioned window).
  [[nodiscard]] bool inPeriodicWindow(double x, double center, double halfWidth) const;

  /// A uniformly random point of that periodic interval. When 'halfWidth' is at least half a period
  /// the sample is uniform on the whole period.
  [[nodiscard]] double samplePeriodicWindow(RandomNumber& random, double center, double halfWidth) const;

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
 * The unit cell of the framework is divided into 'cagesPerUnitCell' (c1, c2, c3) cages along its
 * axes (one cage per cell by default; e.g. [2, 2, 2] for the eight alpha cages of the conventional
 * LTA cell): cage (i, j, k) is the parallelepiped of the cage cell (the unit-cell vectors divided by
 * c) centered on 'origin' + cell (i, j, k), with (i, j, k) taken modulo the number of cages of the
 * supercell. A bead belongs to the cage whose center is nearest in fractional coordinates. The cage
 * boundaries along axis i are the lattice planes s_i = const, perpendicular to the reciprocal vector
 * b_i and 'planeSpacing(i)' apart: the cages are consistent with a barrier only when its normal is
 * b_i / |b_i| for some axis and its lattice distance is that spacing ('inconsistencyWith').
 */
export struct DcTSTCageLattice
{
  DcTSTCageLattice() = default;
  DcTSTCageLattice(const SimulationBox& unitCell, int3 numberOfUnitCells, double3 cageCenter,
                   int3 cagesPerUnitCell = int3(1, 1, 1));

  double3x3 cell{};         ///< Cage cell (columns are the unit-cell vectors divided by 'cagesPerUnitCell').
  double3x3 inverseCell{};  ///< Its inverse.
  int3 numberOfUnitCells{1, 1, 1};  ///< Number of cages along each axis of the supercell.
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

  /// The cage axis whose boundary planes are most nearly perpendicular to 'direction' (the reciprocal vector
  /// b_axis most aligned with it), the sign of the alignment (+1 if b_axis points along 'direction', so that a
  /// step +sign along the axis moves along 'direction') and |cos| of the angle between them.
  struct AxisAlignment
  {
    int axis{0};
    int sign{1};
    double cosine{1.0};
  };
  [[nodiscard]] AxisAlignment axisAlongDirection(double3 direction) const;

  /// Distance between consecutive cage boundary planes s_axis = const (1 / |b_axis|; the length of the cage
  /// vector for an orthogonal cell).
  [[nodiscard]] double planeSpacing(int axis) const;

  /// Why the cages do not partition space along the barrier ('normal', 'lambda'), or nothing when they do: the
  /// normal is perpendicular to the cage planes of one axis, the lattice distance is their spacing, and the
  /// barrier plane passes midway between cage centers.
  [[nodiscard]] std::optional<std::string> inconsistencyWith(const SimulationBox& simulationBox,
                                                             const DcTSTBarrier& barrier) const;
};
