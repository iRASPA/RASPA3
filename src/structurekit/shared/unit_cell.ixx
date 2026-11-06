module;

export module unit_cell;

import std;

import double3;
import int3;
import double3x3;

#if defined(__GNUC__)
#define ALWAYS_INLINE __attribute__((__always_inline__))
#elif defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#endif

// The periodic cell a structure is measured in, and the only thing here that knows about periodicity.
//
// This is the cell as an analysis needs it and not as a simulation does: a shape to wrap distances with and
// to convert between fractional and Cartesian coordinates. Nothing here can be scaled by a barostat,
// restored from an archive, averaged over a run or printed as part of a system, because none of that is what
// a pore diameter or a surface area is measured against. Keeping it separate from the engine's box is what
// lets the whole library be built and tested without the engine, in the way the symmetry code already is.
//
// The lengths and angles are carried alongside the matrix rather than derived on demand, so that a cell built
// from them reports the ones it was given and not those a matrix returns a rounding out.
export struct UnitCell
{
  // Whether the cell is a box lying along the axes: three right angles, and the matrix the three lengths down
  // its diagonal. It only decides which form of the minimum image convention is used, the cheap one reading
  // nothing but that diagonal. A box that is turned, or has its axes in another order, or has one pointing
  // backwards, has the right angles and not the diagonal, and is wrapped the general way.
  enum class Type : int
  {
    Rectangular = 0,
    Triclinic = 1
  };

  UnitCell() = default;

  // A box of the given edge lengths, which is what a test wants when the shape of the cell is not the thing
  // being tested.
  explicit UnitCell(double a, double b, double c);

  // From lengths and angles, in radians, as a crystallographic description gives them. The type is read off
  // the angles.
  explicit UnitCell(double a, double b, double c, double alpha, double beta, double gamma);
  explicit UnitCell(double a, double b, double c, double alpha, double beta, double gamma, Type type);

  // From a cell matrix, with the type read off it: a cell whose three angles are right to within 1e-5 radians
  // and whose vectors lie along x, y and z to the same, pointing forwards, is treated as rectangular.
  explicit UnitCell(double3x3 cell);

  // From a cell matrix whose type is already known, for a caller that has been told rather than having to
  // guess. Told rectangular of a cell that is not a box along the axes, it is taken as triclinic, the quick
  // form being wrong for it.
  explicit UnitCell(double3x3 cell, Type type);

  bool operator==(UnitCell const&) const = default;

  double3x3 cell{};
  double3x3 inverseCell{};
  double volume{0.0};  // Å³, never negative, whichever hand the three vectors are

  double lengthA{0.0};
  double lengthB{0.0};
  double lengthC{0.0};
  double angleAlpha{0.0};
  double angleBeta{0.0};
  double angleGamma{0.0};

  Type type{Type::Triclinic};

  // The same lattice by three vectors that, with minus their sum for a fourth, meet one another at no acute
  // angle: an obtuse superbase, which every lattice in three dimensions has (Selling). It is what the shortest
  // image is found in. In fractions of the cell as given, the image nearest in fractions is the one inside the
  // cell-shaped box around the origin, and that box is not the set of points nearer the origin than any other
  // lattice point unless the cell is a rectangular one: with two axes at sixty degrees the corner of it is half
  // as far again from the origin as from the lattice point next door. In this base the lattice points that
  // bound that set are known outright, the three vectors, their sums in pairs and the sum of all three, and
  // their opposites (Conway and Sloane, Proc. R. Soc. A 436, 55 (1992)), so a vector that none of those seven
  // shortens is the shortest there is, whatever cell of the lattice was handed in.
  double3x3 obtuseCell{};
  double3x3 obtuseInverse{};
  double3x3 obtuseSteps{};                       // whole numbers: obtuseCell = cell * obtuseSteps
  std::array<double, 7> obtuseHalfSquares{};     // half the squared length of each of the seven

  // The shortest of the images of `dr`, which is the separation to use for a pair when the structure repeats
  // for ever, and it is the shortest in Ångström, in a cell of any shape.
  ALWAYS_INLINE inline double3 applyPeriodicBoundaryConditions(const double3& dr) const
  {
    switch (type)
    {
      case Type::Rectangular:
      {
        double3 s;
        s.x = dr.x - static_cast<double>(static_cast<std::make_signed_t<std::size_t>>(dr.x * inverseCell.ax +
                                                                                      ((dr.x >= 0.0) ? 0.5 : -0.5))) *
                         cell.ax;
        s.y = dr.y - static_cast<double>(static_cast<std::make_signed_t<std::size_t>>(dr.y * inverseCell.by +
                                                                                      ((dr.y >= 0.0) ? 0.5 : -0.5))) *
                         cell.by;
        s.z = dr.z - static_cast<double>(static_cast<std::make_signed_t<std::size_t>>(dr.z * inverseCell.cz +
                                                                                      ((dr.z >= 0.0) ? 0.5 : -0.5))) *
                         cell.cz;
        return s;
      }
      default:
      {
        // compute value in variable first to avoid microsoft compiler bug
        double3 r = shortestImage(dr, nullptr);
        return r;
      }
    }
  }

  // The same, for a cell of either type, with which image it is: the whole cell vectors taken off `dr` to
  // leave it, so that `dr - cell * translation` is what comes back.
  double3 shortestImage(const double3& dr, int3* translation) const;

  // Those whole cell vectors alone, for a caller that has to know where the other end went.
  int3 translationToShortestImage(const double3& dr) const;

  // The distance between the two faces of the cell perpendicular to each axis. A fractional displacement of
  // t along an axis is at least |t| times this in Cartesian terms, which is how far a search over images has
  // to go and how much room a grid spacing stands for.
  double3 perpendicularWidths() const;

  // The longest body diagonal of the box spanned by the given fraction of each cell vector: (1, 1, 1) is the
  // cell itself, and (1/nx, 1/ny, 1/nz) is one voxel of a grid over it. Half of it is how far a point of the
  // box can be from the nearest of its corners, whatever the shape of the box. The sum of the three vectors
  // is one of the four diagonals and not always the longest: wherever two of them meet at an obtuse angle
  // their difference is longer than their sum, so it is the longest of the four that has to be taken.
  double longestDiagonal(double3 fractions = double3(1.0, 1.0, 1.0)) const;

  // How many unit cells are needed along each axis before a cutoff of `cutOff` fits inside half the cell,
  // which is the condition for the minimum image convention to be the whole interaction.
  int3 smallestNumberOfUnitCellsForMinimumImagesConvention(double cutOff) const;

  // The same cell repeated `scale` times along each axis: its own three vectors, each that many times over,
  // so it lies as the cell lies and is the same hand, and a point of the cell is where it was.
  UnitCell scaled(int3 scale) const;

 private:
  // What every constructor ends with: the type held to what the matrix can bear, and the superbase.
  void settle();
};
