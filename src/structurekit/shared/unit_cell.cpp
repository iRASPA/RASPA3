module;

module unit_cell;

import std;

import double3;
import int3;
import double3x3;

namespace
{

constexpr double rightAngle = 0.5 * std::numbers::pi;
constexpr double angleTolerance = 1.0e-5;  // radians, and the same as a share of a length

bool hasRightAngles(const UnitCell &unitCell)
{
  return std::fabs(unitCell.angleAlpha - rightAngle) < angleTolerance &&
         std::fabs(unitCell.angleBeta - rightAngle) < angleTolerance &&
         std::fabs(unitCell.angleGamma - rightAngle) < angleTolerance;
}

// Whether the cell is the matrix with three positive lengths down its diagonal, to the tolerance the angles
// are held to, which is what the quick form of the minimum image convention takes it for.
bool liesAlongTheAxes(const double3x3 &cell)
{
  if (!(cell.ax > 0.0 && cell.by > 0.0 && cell.cz > 0.0)) return false;
  return std::fabs(cell.ay) <= angleTolerance * cell.ax && std::fabs(cell.az) <= angleTolerance * cell.ax &&
         std::fabs(cell.bx) <= angleTolerance * cell.by && std::fabs(cell.bz) <= angleTolerance * cell.by &&
         std::fabs(cell.cx) <= angleTolerance * cell.cz && std::fabs(cell.cy) <= angleTolerance * cell.cz;
}

// The angle between two cell vectors. The cosine of two vectors along one line comes out a rounding error past
// one as often as not, and that of a vector of no length is nothing over nothing; neither is left to the arc
// cosine, which has no answer for them. A vector of no length makes no angle and is given none.
double angleBetween(const double3 &u, const double3 &v)
{
  const double lengths = u.length() * v.length();
  if (!(lengths > 0.0)) return 0.0;
  return std::acos(std::clamp(double3::dot(u, v) / lengths, -1.0, 1.0));
}

// The seven lattice vectors of an obtuse superbase that can shorten a vector, as whole numbers of its three.
constexpr std::array<std::array<int, 3>, 7> shorteners{
    {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}};

}  // namespace

UnitCell::UnitCell(double a, double b, double c)
    : lengthA(a),
      lengthB(b),
      lengthC(c),
      angleAlpha(0.5 * std::numbers::pi),
      angleBeta(0.5 * std::numbers::pi),
      angleGamma(0.5 * std::numbers::pi),
      type(Type::Rectangular)
{
  this->cell = double3x3(double3(a, 0.0, 0.0), double3(0.0, b, 0.0), double3(0.0, 0.0, c));
  if (a != 0.0 && b != 0.0 && c != 0.0)
  {
    this->inverseCell = this->cell.inverse();
    this->volume = std::fabs(this->cell.determinant());
  }
  settle();
}

UnitCell::UnitCell(double a, double b, double c, double alpha, double beta, double gamma)
    : UnitCell(a, b, c, alpha, beta, gamma, Type::Triclinic)
{
  if (hasRightAngles(*this) && liesAlongTheAxes(this->cell)) type = Type::Rectangular;
}

UnitCell::UnitCell(double a, double b, double c, double alpha, double beta, double gamma, Type type)
    : lengthA(a), lengthB(b), lengthC(c), angleAlpha(alpha), angleBeta(beta), angleGamma(gamma), type(type)
{
  // exact zeros for right angles (see 'double3x3::fromCellParameters')
  this->cell = double3x3::fromCellParameters(a, b, c, alpha, beta, gamma);
  if (a != 0.0 && b != 0.0 && c != 0.0)
  {
    this->inverseCell = this->cell.inverse();
    this->volume = std::fabs(this->cell.determinant());
  }
  settle();
}

UnitCell::UnitCell(double3x3 cell) : UnitCell(cell, Type::Triclinic)
{
  if (hasRightAngles(*this) && liesAlongTheAxes(this->cell)) type = Type::Rectangular;
}

UnitCell::UnitCell(double3x3 cell, Type type) : cell(cell), type(type)
{
  // The size of the determinant, which is the volume. Its sign is which hand the three vectors are, and that
  // is the order they were named in and nothing about the cell: a volume that came out below zero for half the
  // orders made every width, spacing and density taken from it negative as well.
  this->volume = std::fabs(cell.determinant());
  // a cell of no volume has no inverse, and dividing by its determinant would fill one with what is not a number
  if (std::fabs(this->volume) > 0.0) this->inverseCell = cell.inverse();

  double3 column1 = this->cell[0];
  double3 column2 = this->cell[1];
  double3 column3 = this->cell[2];
  this->lengthA = column1.length();
  this->lengthB = column2.length();
  this->lengthC = column3.length();

  this->angleAlpha = angleBetween(column2, column3);
  this->angleBeta = angleBetween(column1, column3);
  this->angleGamma = angleBetween(column1, column2);
  settle();
}

void UnitCell::settle()
{
  // The quick form of the convention reads the diagonal of the matrix and nothing else, so it is for a matrix
  // that is its diagonal, whatever the caller said of it.
  if (this->type == Type::Rectangular && !liesAlongTheAxes(this->cell)) this->type = Type::Triclinic;

  this->obtuseCell = this->cell;
  this->obtuseInverse = this->inverseCell;
  this->obtuseSteps = double3x3(double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0));
  this->obtuseHalfSquares = {};
  if (!(std::fabs(this->volume) > 0.0)) return;

  // Selling's reduction. The three cell vectors and minus their sum are a superbase, four vectors adding to
  // nothing of which any three span the lattice. While two of them meet at an acute angle, one of the two is
  // turned round and added to the other two that are left, which is again a superbase of the same lattice and
  // shorter in the sum of its squares, so it ends, and it ends with no acute angle among them. What is carried
  // alongside is each vector in whole numbers of the cell as given.
  std::array<double3, 4> base{this->cell[0], this->cell[1], this->cell[2],
                              -1.0 * (this->cell[0] + this->cell[1] + this->cell[2])};
  std::array<double3, 4> steps{double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0),
                               double3(-1.0, -1.0, -1.0)};
  for (int pass = 0; pass < 10000; ++pass)
  {
    std::size_t first = 4;
    std::size_t second = 4;
    double sharpest = 0.0;
    for (std::size_t i = 0; i < 4; ++i)
    {
      for (std::size_t j = i + 1; j < 4; ++j)
      {
        // a right angle is not acute, and one that is right to rounding is not to be chased for ever
        const double along = double3::dot(base[i], base[j]);
        const double cosine = along / (base[i].length() * base[j].length());
        if (cosine > 1.0e-12 && along > sharpest)
        {
          sharpest = along;
          first = i;
          second = j;
        }
      }
    }
    if (first == 4) break;

    for (std::size_t k = 0; k < 4; ++k)
    {
      if (k == first || k == second) continue;
      base[k] += base[first];
      steps[k] += steps[first];
    }
    base[first] = -1.0 * base[first];
    steps[first] = -1.0 * steps[first];
  }

  this->obtuseCell = double3x3(base[0], base[1], base[2]);
  this->obtuseInverse = this->obtuseCell.inverse();
  this->obtuseSteps = double3x3(steps[0], steps[1], steps[2]);
  for (std::size_t k = 0; k < shorteners.size(); ++k)
  {
    const double3 vector = this->obtuseCell * double3(shorteners[k][0], shorteners[k][1], shorteners[k][2]);
    this->obtuseHalfSquares[k] = 0.5 * double3::dot(vector, vector);
  }
}

double3 UnitCell::shortestImage(const double3 &dr, int3 *translation) const
{
  if (translation != nullptr) *translation = int3(0, 0, 0);

  if (this->type == Type::Rectangular)
  {
    const double3 whole(std::round(dr.x * this->inverseCell.ax), std::round(dr.y * this->inverseCell.by),
                        std::round(dr.z * this->inverseCell.cz));
    if (translation != nullptr)
    {
      *translation = int3(static_cast<std::int32_t>(whole.x), static_cast<std::int32_t>(whole.y),
                          static_cast<std::int32_t>(whole.z));
    }
    return double3(dr.x - whole.x * this->cell.ax, dr.y - whole.y * this->cell.by, dr.z - whole.z * this->cell.cz);
  }

  // A cell of no volume repeats nothing, and a displacement in it has no other image.
  if (!(std::fabs(this->volume) > 0.0)) return dr;

  // The image nearest in fractions of the obtuse superbase to begin with, which is the shortest outright unless
  // it lies out towards a corner of that cell, and then whichever of the seven vectors shortens it most is
  // taken off, until none of them does. Every step shortens it, so this ends, and where it ends no lattice
  // point is nearer than the origin.
  const double3 fractional = this->obtuseInverse * dr;
  double3 whole(std::round(fractional.x), std::round(fractional.y), std::round(fractional.z));
  double3 image = dr - this->obtuseCell * whole;

  const double3 first = this->obtuseCell[0];
  const double3 second = this->obtuseCell[1];
  const double3 third = this->obtuseCell[2];
  for (int step = 0; step < 64; ++step)
  {
    const double alongFirst = double3::dot(image, first);
    const double alongSecond = double3::dot(image, second);
    const double alongThird = double3::dot(image, third);
    const std::array<double, 7> along{alongFirst,
                                      alongSecond,
                                      alongThird,
                                      alongFirst + alongSecond,
                                      alongFirst + alongThird,
                                      alongSecond + alongThird,
                                      alongFirst + alongSecond + alongThird};

    // Taking a vector L off shortens the square of the image by 2 image.L - L.L, and adding it by the same with
    // the sign of image.L turned, so half the gain either way is |image.L| less half of L.L. One that gains
    // nothing to rounding is left: two images of one length are both the shortest.
    std::size_t best = shorteners.size();
    double mostGained = 0.0;
    for (std::size_t k = 0; k < shorteners.size(); ++k)
    {
      const double gained = std::fabs(along[k]) - this->obtuseHalfSquares[k];
      if (gained > mostGained && gained > 1.0e-12 * this->obtuseHalfSquares[k])
      {
        mostGained = gained;
        best = k;
      }
    }
    if (best == shorteners.size()) break;

    const double sign = along[best] > 0.0 ? 1.0 : -1.0;
    whole += sign * double3(shorteners[best][0], shorteners[best][1], shorteners[best][2]);
    image = dr - this->obtuseCell * whole;
  }

  if (translation != nullptr)
  {
    const double3 given = this->obtuseSteps * whole;
    *translation =
        int3(static_cast<std::int32_t>(std::lround(given.x)), static_cast<std::int32_t>(std::lround(given.y)),
             static_cast<std::int32_t>(std::lround(given.z)));
  }
  return image;
}

int3 UnitCell::translationToShortestImage(const double3 &dr) const
{
  int3 translation(0, 0, 0);
  shortestImage(dr, &translation);
  return translation;
}

double3 UnitCell::perpendicularWidths() const
{
  double3 c1 = double3::cross(cell[0], cell[1]);
  double3 c2 = double3::cross(cell[1], cell[2]);
  double3 c3 = double3::cross(cell[0], cell[2]);

  double v = std::fabs(cell[0].x * c2.x + cell[0].y * c2.y + cell[0].z * c2.z);

  return double3(v / std::sqrt(double3::dot(c2, c2)), v / std::sqrt(double3::dot(c3, c3)),
                 v / std::sqrt(double3::dot(c1, c1)));
}

double UnitCell::longestDiagonal(double3 fractions) const
{
  double longest = 0.0;
  for (double sy : {-1.0, 1.0})
  {
    for (double sz : {-1.0, 1.0})
    {
      longest = std::max(longest, (cell * double3(fractions.x, sy * fractions.y, sz * fractions.z)).length());
    }
  }
  return longest;
}

int3 UnitCell::smallestNumberOfUnitCellsForMinimumImagesConvention(double cutOff) const
{
  double3 widths = perpendicularWidths();
  return int3(static_cast<std::int32_t>(std::ceil(2.0 * cutOff / widths.x)),
              static_cast<std::int32_t>(std::ceil(2.0 * cutOff / widths.y)),
              static_cast<std::int32_t>(std::ceil(2.0 * cutOff / widths.z)));
}

UnitCell UnitCell::scaled(int3 scale) const
{
  // The cell's own vectors, each so many times over, and built by the constructor so that what the constructors
  // work out of a matrix is worked out of this one.
  //
  // It was rebuilt from the lengths and the angles before, which say what shape a cell is and nothing of how
  // it lies: the repeated cell came out with its first vector along x and its second in the xy-plane, and
  // right-handed, whatever the cell itself was. For a cell that lies that way already it is the same cell.
  // For one given by a matrix turned any other way, or left-handed, it was the right shape somewhere else, and
  // an atom placed in it by its position in space was placed in the wrong cell.
  const double3 counts(static_cast<double>(scale.x), static_cast<double>(scale.y), static_cast<double>(scale.z));
  UnitCell repeated(double3x3(counts.x * cell[0], counts.y * cell[1], counts.z * cell[2]), type);

  // The lengths and the angles are known without reading them back off the matrix, which would return them
  // a rounding out: the lengths are the cell's so many times over, and repeating a cell turns nothing.
  repeated.lengthA = counts.x * lengthA;
  repeated.lengthB = counts.y * lengthB;
  repeated.lengthC = counts.z * lengthC;
  repeated.angleAlpha = angleAlpha;
  repeated.angleBeta = angleBeta;
  repeated.angleGamma = angleGamma;
  return repeated;
}
