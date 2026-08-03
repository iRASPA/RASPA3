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
// The lengths and angles are carried alongside the matrix rather than derived on demand because `scaled`
// rebuilds the cell from them, in the canonical orientation with a along x and b in the xy-plane.
export struct UnitCell
{
  // Whether the cell is a box. It only decides which form of the minimum image convention is used, the
  // cheap one being available when the axes are orthogonal.
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

  // From a cell matrix, with the type read off the angles: a cell whose three angles are right to within
  // 1e-5 radians is treated as rectangular.
  explicit UnitCell(double3x3 cell);

  // From a cell matrix whose type is already known, for a caller that has been told rather than having to
  // guess.
  explicit UnitCell(double3x3 cell, Type type);

  bool operator==(UnitCell const&) const = default;

  double3x3 cell{};
  double3x3 inverseCell{};
  double volume{0.0};

  double lengthA{0.0};
  double lengthB{0.0};
  double lengthC{0.0};
  double angleAlpha{0.0};
  double angleBeta{0.0};
  double angleGamma{0.0};

  Type type{Type::Triclinic};

  // The shortest of the images of `dr`, which is the separation to use for a pair when the structure repeats
  // for ever. Hot enough to be worth writing out for both cell shapes.
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
        double3 s = inverseCell * dr;
        s.x -= static_cast<double>(static_cast<std::make_signed_t<std::size_t>>(s.x + ((s.x >= 0.0) ? 0.5 : -0.5)));
        s.y -= static_cast<double>(static_cast<std::make_signed_t<std::size_t>>(s.y + ((s.y >= 0.0) ? 0.5 : -0.5)));
        s.z -= static_cast<double>(static_cast<std::make_signed_t<std::size_t>>(s.z + ((s.z >= 0.0) ? 0.5 : -0.5)));

        // compute value in variable first to avoid microsoft compiler bug
        double3 r = cell * s;
        return r;
      }
    }
  }

  // The distance between the two faces of the cell perpendicular to each axis. A fractional displacement of
  // t along an axis is at least |t| times this in Cartesian terms, which is how far a search over images has
  // to go and how much room a grid spacing stands for.
  double3 perpendicularWidths() const;

  // How many unit cells are needed along each axis before a cutoff of `cutOff` fits inside half the cell,
  // which is the condition for the minimum image convention to be the whole interaction.
  int3 smallestNumberOfUnitCellsForMinimumImagesConvention(double cutOff) const;

  // The same cell repeated `scale` times along each axis.
  UnitCell scaled(int3 scale) const;
};
