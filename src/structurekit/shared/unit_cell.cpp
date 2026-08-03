module;

module unit_cell;

import std;

import double3;
import int3;
import double3x3;

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
    this->volume = this->cell.determinant();
  }
}

UnitCell::UnitCell(double a, double b, double c, double alpha, double beta, double gamma)
    : UnitCell(a, b, c, alpha, beta, gamma, Type::Triclinic)
{
  if (((std::fabs(angleAlpha - 0.5 * std::numbers::pi)) < 1e-5) &&
      ((std::fabs(angleBeta - 0.5 * std::numbers::pi)) < 1e-5) &&
      ((std::fabs(angleGamma - 0.5 * std::numbers::pi)) < 1e-5))
  {
    type = Type::Rectangular;
  }
}

UnitCell::UnitCell(double a, double b, double c, double alpha, double beta, double gamma, Type type)
    : lengthA(a), lengthB(b), lengthC(c), angleAlpha(alpha), angleBeta(beta), angleGamma(gamma), type(type)
{
  double temp = (std::cos(alpha) - std::cos(gamma) * std::cos(beta)) / std::sin(gamma);

  double3 v1 = double3(a, 0.0, 0.0);
  double3 v2 = double3(b * std::cos(gamma), b * std::sin(gamma), 0.0);
  double3 v3 =
      double3(c * std::cos(beta), c * temp, c * std::sqrt(1.0 - std::cos(beta) * std::cos(beta) - temp * temp));
  this->cell = double3x3(v1, v2, v3);
  if (a != 0.0 && b != 0.0 && c != 0.0)
  {
    this->inverseCell = this->cell.inverse();
    this->volume = this->cell.determinant();
  }
}

UnitCell::UnitCell(double3x3 cell) : UnitCell(cell, Type::Triclinic)
{
  if (((std::fabs(angleAlpha - 0.5 * std::numbers::pi)) < 1e-5) &&
      ((std::fabs(angleBeta - 0.5 * std::numbers::pi)) < 1e-5) &&
      ((std::fabs(angleGamma - 0.5 * std::numbers::pi)) < 1e-5))
  {
    type = Type::Rectangular;
  }
}

UnitCell::UnitCell(double3x3 cell, Type type) : cell(cell), type(type)
{
  this->inverseCell = cell.inverse();
  this->volume = cell.determinant();

  double3 column1 = this->cell[0];
  double3 column2 = this->cell[1];
  double3 column3 = this->cell[2];
  this->lengthA = column1.length();
  this->lengthB = column2.length();
  this->lengthC = column3.length();

  this->angleAlpha = std::acos(double3::dot(column2, column3) / (this->lengthB * this->lengthC));
  this->angleBeta = std::acos(double3::dot(column1, column3) / (this->lengthA * this->lengthC));
  this->angleGamma = std::acos(double3::dot(column1, column2) / (this->lengthA * this->lengthB));
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

int3 UnitCell::smallestNumberOfUnitCellsForMinimumImagesConvention(double cutOff) const
{
  double3 widths = perpendicularWidths();
  return int3(static_cast<std::int32_t>(std::ceil(2.0 * cutOff / widths.x)),
              static_cast<std::int32_t>(std::ceil(2.0 * cutOff / widths.y)),
              static_cast<std::int32_t>(std::ceil(2.0 * cutOff / widths.z)));
}

UnitCell UnitCell::scaled(int3 scale) const
{
  UnitCell v;

  v.lengthA = static_cast<double>(scale.x) * lengthA;
  v.lengthB = static_cast<double>(scale.y) * lengthB;
  v.lengthC = static_cast<double>(scale.z) * lengthC;
  v.angleAlpha = angleAlpha;
  v.angleBeta = angleBeta;
  v.angleGamma = angleGamma;

  double temp = (std::cos(v.angleAlpha) - std::cos(v.angleGamma) * std::cos(v.angleBeta)) / std::sin(v.angleGamma);
  double3 v1 = double3(v.lengthA, 0.0, 0.0);
  double3 v2 = double3(v.lengthB * std::cos(v.angleGamma), v.lengthB * std::sin(v.angleGamma), 0.0);
  double3 v3 = double3(v.lengthC * std::cos(v.angleBeta), v.lengthC * temp,
                       v.lengthC * std::sqrt(1.0 - std::cos(v.angleBeta) * std::cos(v.angleBeta) - temp * temp));
  v.cell = double3x3(v1, v2, v3);
  v.inverseCell = v.cell.inverse();
  v.volume = v.cell.determinant();
  v.type = type;

  return v;
}
