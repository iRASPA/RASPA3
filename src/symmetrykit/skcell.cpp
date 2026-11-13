module;

module skcell;

import std;

import int3;
import bool3;
import double3;
import float4;
import double3x3;
import skboundingbox;

namespace
{
// exact 0 and 1 for right angles (see 'double3x3::cosCellAngle')
// The transformation that takes the canonical frame of double3x3::fromCellParameters (a along x, b in the xy-plane)
// to the frame the cell is actually in: the identity for a cell in the canonical frame (kept exactly, so that its
// zeros stay exact), otherwise the cell times the inverse of its canonical form. Axes rebuilt from new lattice
// parameters are carried over with it, so that a rotated cell keeps its orientation.
double3x3 orientation(const SKCell& cell)
{
  const double3x3 unitCell = cell.unitCell();
  // only a triangular cell with a positive diagonal is the canonical form: a triangular one with a negative diagonal
  // entry is rotated (by 180 degrees about an axis) or left-handed, and keeps that orientation through the product
  if (unitCell.ay == 0.0 && unitCell.az == 0.0 && unitCell.bz == 0.0 && unitCell.ax > 0.0 && unitCell.by > 0.0 &&
      unitCell.cz > 0.0)
  {
    return double3x3::identity();
  }
  return unitCell *
         double3x3::fromCellParameters(cell.a(), cell.b(), cell.c(), cell.alpha(), cell.beta(), cell.gamma()).inverse();
}
}  // namespace

SKCell::SKCell()
{
  double a = 20.0;
  double b = 20.0;
  double c = 20.0;
  double alpha = 90.0 * std::numbers::pi / 180.0;
  double beta = 90.0 * std::numbers::pi / 180.0;
  double gamma = 90.0 * std::numbers::pi / 180.0;

  _unitCell = double3x3::fromCellParameters(a, b, c, alpha, beta, gamma);

  _inverseUnitCell = _unitCell.inverse();
  _fullCell = _unitCell;

  int dx = _maximumReplica.x - _minimumReplica.x + 1;
  int dy = _maximumReplica.y - _minimumReplica.y + 1;
  int dz = _maximumReplica.z - _minimumReplica.z + 1;

  _fullCell.ax *= dx;
  _fullCell.bx *= dy;
  _fullCell.cx *= dz;
  _fullCell.ay *= dx;
  _fullCell.by *= dy;
  _fullCell.cy *= dz;
  _fullCell.az *= dx;
  _fullCell.bz *= dy;
  _fullCell.cz *= dz;

  _inverseFullCell = _fullCell.inverse();

  _boundingBox = enclosingBoundingBox();
}

SKCell::SKCell(double3 v1, double3 v2, double3 v3)
{
  _unitCell = double3x3(v1, v2, v3);
  _inverseUnitCell = _unitCell.inverse();

  _fullCell = _unitCell;
  _inverseFullCell = _fullCell.inverse();

  _boundingBox = enclosingBoundingBox();
}

SKCell::SKCell(double3x3 m)
{
  _unitCell = m;
  _inverseUnitCell = _unitCell.inverse();

  _fullCell = m;
  _inverseFullCell = _fullCell.inverse();

  _boundingBox = enclosingBoundingBox();
}

SKCell::SKCell(double a, double b, double c, double alpha, double beta, double gamma)
{
  // exact zeros for right angles (see 'double3x3::fromCellParameters')
  _unitCell = double3x3::fromCellParameters(a, b, c, alpha, beta, gamma);

  _inverseUnitCell = _unitCell.inverse();
  _fullCell = _unitCell;

  int dx = _maximumReplica.x - _minimumReplica.x + 1;
  int dy = _maximumReplica.y - _minimumReplica.y + 1;
  int dz = _maximumReplica.z - _minimumReplica.z + 1;

  _fullCell.ax *= dx;
  _fullCell.bx *= dy;
  _fullCell.cx *= dz;
  _fullCell.ay *= dx;
  _fullCell.by *= dy;
  _fullCell.cy *= dz;
  _fullCell.az *= dx;
  _fullCell.bz *= dy;
  _fullCell.cz *= dz;

  _inverseFullCell = _fullCell.inverse();

  _boundingBox = enclosingBoundingBox();
}

SKCell::SKCell(SKBoundingBox boundingBox)
{
  _minimumReplica = int3(0, 0, 0);
  _maximumReplica = int3(0, 0, 0);

  double3 v1 = double3(std::fabs(boundingBox.maximum().x - boundingBox.minimum().x), 0.0, 0.0);
  double3 v2 = double3(0.0, std::fabs(boundingBox.maximum().y - boundingBox.minimum().y), 0.0);
  double3 v3 = double3(0.0, 0.0, std::fabs(boundingBox.maximum().z - boundingBox.minimum().z));

  _unitCell = double3x3(v1, v2, v3);
  _inverseUnitCell = _unitCell.inverse();
  _fullCell = _unitCell;
  _inverseFullCell = _fullCell.inverse();
  _boundingBox = enclosingBoundingBox();
}

SKCell::SKCell(SKCell superCell, int3 minimumReplica, int3 maximumReplica)
{
  double3 v1 = double(maximumReplica.x - minimumReplica.x + 1) * superCell._unitCell[0];
  double3 v2 = double(maximumReplica.y - minimumReplica.y + 1) * superCell._unitCell[1];
  double3 v3 = double(maximumReplica.z - minimumReplica.z + 1) * superCell._unitCell[2];
  _unitCell = double3x3(v1, v2, v3);
  _inverseUnitCell = _unitCell.inverse();
  _fullCell = _unitCell;
  _inverseFullCell = _fullCell.inverse();

  _minimumReplica = int3(0, 0, 0);
  _maximumReplica = int3(0, 0, 0);
  _boundingBox = enclosingBoundingBox();
}

std::shared_ptr<SKCell> SKCell::superCell() const
{
  double3 v1 = double(_maximumReplica.x - _minimumReplica.x + 1) * _unitCell[0];
  double3 v2 = double(_maximumReplica.y - _minimumReplica.y + 1) * _unitCell[1];
  double3 v3 = double(_maximumReplica.z - _minimumReplica.z + 1) * _unitCell[2];
  return std::make_shared<SKCell>(v1, v2, v3);
}

void SKCell::setUnitCell(const double3x3& unitCell)
{
  _unitCell = unitCell;

  _inverseUnitCell = _unitCell.inverse();
  updateFullCell();
}

int3 SKCell::minimumReplicas() const { return _minimumReplica; }

void SKCell::setMinimumReplicas(const int3& minimumReplicas)
{
  _minimumReplica = minimumReplicas;
  updateFullCell();
}

int3 SKCell::maximumReplicas() const { return _maximumReplica; }

void SKCell::setMaximumReplicas(const int3& maximumReplicas)
{
  _maximumReplica = maximumReplicas;
  updateFullCell();
}

SKBoundingBox SKCell::enclosingBoundingBox()
{
  double3 c0 = _unitCell * double3(double(_minimumReplica.x), double(_minimumReplica.y), double(_minimumReplica.z));
  double3 c1 = _unitCell * double3(double(_maximumReplica.x + 1), double(_minimumReplica.y), double(_minimumReplica.z));
  double3 c2 =
      _unitCell * double3(double(_maximumReplica.x + 1), double(_maximumReplica.y + 1), double(_minimumReplica.z));
  double3 c3 = _unitCell * double3(double(_minimumReplica.x), double(_maximumReplica.y + 1), double(_minimumReplica.z));
  double3 c4 = _unitCell * double3(double(_minimumReplica.x), double(_minimumReplica.y), double(_maximumReplica.z + 1));
  double3 c5 =
      _unitCell * double3(double(_maximumReplica.x + 1), double(_minimumReplica.y), double(_maximumReplica.z + 1));
  double3 c6 =
      _unitCell * double3(double(_maximumReplica.x + 1), double(_maximumReplica.y + 1), double(_maximumReplica.z + 1));
  double3 c7 =
      _unitCell * double3(double(_minimumReplica.x), double(_maximumReplica.y + 1), double(_maximumReplica.z + 1));

  double listx[] = {c0.x, c1.x, c2.x, c3.x, c4.x, c5.x, c6.x, c7.x};
  double listy[] = {c0.y, c1.y, c2.y, c3.y, c4.y, c5.y, c6.y, c7.y};
  double listz[] = {c0.z, c1.z, c2.z, c3.z, c4.z, c5.z, c6.z, c7.z};
  double3 minimum = double3(*std::min_element(listx, listx + 8), *std::min_element(listy, listy + 8),
                            *std::min_element(listz, listz + 8));
  double3 maximum = double3(*std::max_element(listx, listx + 8), *std::max_element(listy, listy + 8),
                            *std::max_element(listz, listz + 8));

  return SKBoundingBox(minimum, maximum);
}

std::array<double3, 8> SKCell::corners()
{
  std::array<double3, 8> temp = {{double3(0.0, 0.0, 0.0), double3(1.0, 0.0, 0.0), double3(1.0, 1.0, 0.0),
                                  double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0), double3(1.0, 0.0, 1.0),
                                  double3(1.0, 1.0, 1.0), double3(0.0, 1.0, 1.0)}};
  return temp;
}

int3 SKCell::replicaFromIndex(int index)
{
  int dx = _maximumReplica.x - _minimumReplica.x + 1;
  int dy = _maximumReplica.y - _minimumReplica.y + 1;
  int dz = _maximumReplica.z - _minimumReplica.z + 1;
  int n = dx * dy * dz;

  // the order of renderTranslationVectors: x fastest, z slowest (it used to be decoded the other way round)
  int k1 = (index % n) % dx + _minimumReplica.x;
  int k2 = ((index % n) / dx) % dy + _minimumReplica.y;
  int k3 = ((index % n) / (dx * dy)) + _minimumReplica.z;

  return int3(k1, k2, k3);
}

double3x3 SKCell::box() const { return _fullCell; }
void SKCell::setBox(const double3x3& fullCell)
{
  _fullCell = fullCell;
  _inverseFullCell = _fullCell.inverse();

  int dx = _maximumReplica.x - _minimumReplica.x + 1;
  int dy = _maximumReplica.y - _minimumReplica.y + 1;
  int dz = _maximumReplica.z - _minimumReplica.z + 1;

  // the unit cell is the box divided over the replicas
  _unitCell = fullCell;
  _unitCell.ax /= double(dx);
  _unitCell.bx /= double(dy);
  _unitCell.cx /= double(dz);
  _unitCell.ay /= double(dx);
  _unitCell.by /= double(dy);
  _unitCell.cy /= double(dz);
  _unitCell.az /= double(dx);
  _unitCell.bz /= double(dy);
  _unitCell.cz /= double(dz);

  _inverseUnitCell = _unitCell.inverse();
  _boundingBox = enclosingBoundingBox();
}

int SKCell::totalNumberOfReplicas() const
{
  int dx = _maximumReplica.x - _minimumReplica.x + 1;
  int dy = _maximumReplica.y - _minimumReplica.y + 1;
  int dz = _maximumReplica.z - _minimumReplica.z + 1;

  return dx * dy * dz;
}

int SKCell::minimumReplicaX() const { return _minimumReplica.x; }

void SKCell::setMinimumReplicaX(const int newValue)
{
  _minimumReplica.x = newValue;
  updateFullCell();
}

int SKCell::maximumReplicaX() const { return _maximumReplica.x; }

void SKCell::setMaximumReplicaX(const int newValue)
{
  _maximumReplica.x = newValue;
  updateFullCell();
}

int SKCell::minimumReplicaY() const { return _minimumReplica.y; }

void SKCell::setMinimumReplicaY(const int newValue)
{
  _minimumReplica.y = newValue;
  updateFullCell();
}

int SKCell::maximumReplicaY() const { return _maximumReplica.y; }

void SKCell::setMaximumReplicaY(const int newValue)
{
  _maximumReplica.y = newValue;
  updateFullCell();
}

int SKCell::minimumReplicaZ() const { return _minimumReplica.z; }

void SKCell::setMinimumReplicaZ(const int newValue)
{
  _minimumReplica.z = newValue;
  updateFullCell();
}

int SKCell::maximumReplicaZ() const { return _maximumReplica.z; }

void SKCell::setMaximumReplicaZ(const int newValue)
{
  _maximumReplica.z = newValue;
  updateFullCell();
}

int3 SKCell::numberOfReplicas()
{
  int dx = _maximumReplica.x - _minimumReplica.x + 1;
  int dy = _maximumReplica.y - _minimumReplica.y + 1;
  int dz = _maximumReplica.z - _minimumReplica.z + 1;
  return int3(dx, dy, dz);
}

double SKCell::a() const { return _unitCell[0].length(); }

// The full cell is the unit cell scaled by the replica counts, and the bounding box encloses it; every change of
// the unit cell or the replicas has to carry both along (the a/b/c and alpha/beta/gamma setters failed to do so for
// the full cell, and no setter did for the bounding box).
void SKCell::updateFullCell()
{
  int dx = _maximumReplica.x - _minimumReplica.x + 1;
  int dy = _maximumReplica.y - _minimumReplica.y + 1;
  int dz = _maximumReplica.z - _minimumReplica.z + 1;
  _fullCell = double3x3(double(dx) * _unitCell[0], double(dy) * _unitCell[1], double(dz) * _unitCell[2]);
  _inverseFullCell = _fullCell.inverse();
  _boundingBox = enclosingBoundingBox();
}

void SKCell::a(const double& newValue)
{
  double multiplier = newValue / a();
  _unitCell = double3x3(multiplier * _unitCell[0], _unitCell[1], _unitCell[2]);
  _inverseUnitCell = _unitCell.inverse();
  updateFullCell();
}

void SKCell::setLengthA(double value)
{
  double multiplier = value / a();
  _unitCell = double3x3(multiplier * _unitCell[0], _unitCell[1], _unitCell[2]);
  _inverseUnitCell = _unitCell.inverse();

  updateFullCell();
}

double SKCell::b() const { return _unitCell[1].length(); }
void SKCell::b(const double& newValue)
{
  double multiplier = newValue / b();
  _unitCell = double3x3(_unitCell[0], multiplier * _unitCell[1], _unitCell[2]);
  _inverseUnitCell = _unitCell.inverse();
  updateFullCell();
}

void SKCell::setLengthB(double value)
{
  double multiplier = value / b();
  _unitCell = double3x3(_unitCell[0], multiplier * _unitCell[1], _unitCell[2]);
  _inverseUnitCell = _unitCell.inverse();

  updateFullCell();
}

double SKCell::c() const { return _unitCell[2].length(); }

void SKCell::setLengthC(double value)
{
  double multiplier = value / c();
  _unitCell = double3x3(_unitCell[0], _unitCell[1], multiplier * _unitCell[2]);
  _inverseUnitCell = _unitCell.inverse();

  updateFullCell();
}

void SKCell::c(const double& newValue)
{
  double multiplier = newValue / c();
  _unitCell = double3x3(_unitCell[0], _unitCell[1], multiplier * _unitCell[2]);
  _inverseUnitCell = _unitCell.inverse();
  updateFullCell();
}
double SKCell::alpha() const
{
  double3 column2 = _unitCell[1];
  double3 column3 = _unitCell[2];
  double length2 = column2.length();
  double length3 = column3.length();

  // clamped: round-off puts the cosine of (nearly) collinear axes a little beyond +-1, where acos gives NaN
  return std::acos(std::clamp(double3::dot(column2, column3) / (length2 * length3), -1.0, 1.0));
}
void SKCell::alpha(const double& newValue) { setAlphaAngle(newValue); }

void SKCell::setAlphaAngle(double value)
{
  // all three axes from the lattice parameters, carried into the frame the cell is in: the setter used to put
  // the new c (and b) axis in the canonical frame while leaving the a axis where it was, which bent the angles
  // of a rotated cell
  _unitCell = orientation(*this) * double3x3::fromCellParameters(a(), b(), c(), value, beta(), gamma());
  _inverseUnitCell = _unitCell.inverse();
  updateFullCell();
}

double SKCell::beta() const
{
  double3 column1 = _unitCell[0];
  double3 column3 = _unitCell[2];
  double length1 = column1.length();
  double length3 = column3.length();

  // clamped: round-off puts the cosine of (nearly) collinear axes a little beyond +-1, where acos gives NaN
  return std::acos(std::clamp(double3::dot(column1, column3) / (length1 * length3), -1.0, 1.0));
}

void SKCell::setBetaAngle(double value)
{
  // see setAlphaAngle
  _unitCell = orientation(*this) * double3x3::fromCellParameters(a(), b(), c(), alpha(), value, gamma());
  _inverseUnitCell = _unitCell.inverse();
  updateFullCell();
}

void SKCell::beta(const double& newValue) { setBetaAngle(newValue); }
double SKCell::gamma() const
{
  double3 column1 = _unitCell[0];
  double3 column2 = _unitCell[1];
  double length1 = column1.length();
  double length2 = column2.length();

  // clamped: round-off puts the cosine of (nearly) collinear axes a little beyond +-1, where acos gives NaN
  return std::acos(std::clamp(double3::dot(column1, column2) / (length1 * length2), -1.0, 1.0));
}

void SKCell::setGammaAngle(double value)
{
  // see setAlphaAngle
  _unitCell = orientation(*this) * double3x3::fromCellParameters(a(), b(), c(), alpha(), beta(), value);
  _inverseUnitCell = _unitCell.inverse();
  updateFullCell();
}

void SKCell::gamma(const double& newValue) { setGammaAngle(newValue); }

bool SKCell::orthorhombic() const
{
  return (std::fabs(alpha() - std::numbers::pi / 2.0) < 0.001) &&
         (std::fabs(beta() - std::numbers::pi / 2.0) < 0.001) && (std::fabs(gamma() - std::numbers::pi / 2.0) < 0.001);
}

double SKCell::volume() const
{
  double3 column1 = _unitCell[0];
  double3 column2 = _unitCell[1];
  double3 column3 = _unitCell[2];

  double3 v2 = double3::cross(column2, column3);

  // the triple product is signed: a left-handed cell has a negative determinant but the same volume
  return std::fabs(double3::dot(column1, v2));
}

double3 SKCell::perpendicularWidths() const
{
  double3 column1 = _unitCell[0];
  double3 column2 = _unitCell[1];
  double3 column3 = _unitCell[2];

  double3 v1 = double3::cross(column1, column2);
  double3 v2 = double3::cross(column2, column3);
  double3 v3 = double3::cross(column3, column1);

  // the widths are distances between opposite faces, so the signed triple product of a left-handed cell is not
  // taken as it is (the widths came out negative)
  double volume = std::fabs(double3::dot(column1, v2));

  return double3(volume / v2.length(), volume / v3.length(), volume / v1.length());
}

double3 SKCell::applyFullCellBoundaryCondition(double3 dr)
{
  // convert from xyz to abc
  double3 s = _inverseFullCell * dr;

  // apply boundary condition
  s.x -= std::rint(s.x);
  s.y -= std::rint(s.y);
  s.z -= std::rint(s.z);

  // convert from abc to xyz and return value
  return _fullCell * s;
}

double3 SKCell::applyUnitCellBoundaryCondition(double3 dr)
{
  // convert from xyz to abc
  double3 s = _inverseUnitCell * dr;

  // apply boundary condition
  s.x -= std::rint(s.x);
  s.y -= std::rint(s.y);
  s.z -= std::rint(s.z);

  // convert from abc to xyz and return value
  return _unitCell * s;
}

int3 SKCell::numberOfReplicas(double cutoff)
{
  double3 column1 = _unitCell[0];
  double3 column2 = _unitCell[1];
  double3 column3 = _unitCell[2];

  double3 v1 = double3::cross(column1, column2);
  double3 v2 = double3::cross(column2, column3);
  double3 v3 = double3::cross(column3, column1);

  // as in perpendicularWidths: a left-handed cell gave negative widths and so zero or negative replica counts
  double volume = std::fabs(double3::dot(column1, v2));

  double3 perpendicularWith = double3(volume / v2.length(), volume / v3.length(), volume / v1.length());

  int3 replicas = int3(int(std::ceil(2.0 * cutoff / (perpendicularWith.x + 0.000001))),
                       int(std::ceil(2.0 * cutoff / (perpendicularWith.y + 0.000001))),
                       int(std::ceil(2.0 * cutoff / (perpendicularWith.z + 0.000001))));

  return replicas;
}

std::vector<float4> SKCell::renderTranslationVectors()
{
  std::vector<float4> vectors{};

  for (int k3 = _minimumReplica.z; k3 <= _maximumReplica.z; k3++)
  {
    for (int k2 = _minimumReplica.y; k2 <= _maximumReplica.y; k2++)
    {
      for (int k1 = _minimumReplica.x; k1 <= _maximumReplica.x; k1++)
      {
        vectors.push_back(float4{float(k1), float(k2), float(k3), 0.0f});
      }
    }
  }

  return vectors;
}
