module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <numbers>
#include <vector>
#endif

module skcell;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import bool3;
import double3;
import float4;
import double3x3;
import skboundingbox;

SKCell::SKCell()
{
  double a = 20.0;
  double b = 20.0;
  double c = 20.0;
  double alpha = 90.0 * std::numbers::pi / 180.0;
  double beta = 90.0 * std::numbers::pi / 180.0;
  double gamma = 90.0 * std::numbers::pi / 180.0;

  double temp = (std::cos(alpha) - std::cos(gamma) * std::cos(beta)) / std::sin(gamma);

  double3 v1 = double3(a, 0.0, 0.0);
  double3 v2 = double3(b * std::cos(gamma), b * std::sin(gamma), 0.0);
  double3 v3 =
      double3(c * std::cos(beta), c * temp, c * std::sqrt(1.0 - std::cos(beta) * std::cos(beta) - temp * temp));
  _unitCell = double3x3(v1, v2, v3);

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
  double temp = (std::cos(alpha) - std::cos(gamma) * std::cos(beta)) / std::sin(gamma);

  double3 v1 = double3(a, 0.0, 0.0);
  double3 v2 = double3(b * std::cos(gamma), b * std::sin(gamma), 0.0);
  double3 v3 =
      double3(c * std::cos(beta), c * temp, c * std::sqrt(1.0 - std::cos(beta) * std::cos(beta) - temp * temp));
  _unitCell = double3x3(v1, v2, v3);

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
}

int3 SKCell::minimumReplicas() const { return _minimumReplica; }

void SKCell::setMinimumReplicas(const int3& minimumReplicas)
{
  _minimumReplica = minimumReplicas;

  int dx = _maximumReplica.x - _minimumReplica.x + 1;

  _fullCell.ax = _unitCell.ax * double(dx);
  _fullCell.ay = _unitCell.ay * double(dx);
  _fullCell.az = _unitCell.az * double(dx);

  int dy = _maximumReplica.y - _minimumReplica.y + 1;

  _fullCell.bx = _unitCell.bx * double(dy);
  _fullCell.by = _unitCell.by * double(dy);
  _fullCell.bz = _unitCell.bz * double(dy);

  int dz = _maximumReplica.z - _minimumReplica.z + 1;

  _fullCell.cx = _unitCell.cx * double(dz);
  _fullCell.cy = _unitCell.cy * double(dz);
  _fullCell.cz = _unitCell.cz * double(dz);

  _inverseFullCell = _fullCell.inverse();
}

int3 SKCell::maximumReplicas() const { return _maximumReplica; }

void SKCell::setMaximumReplicas(const int3& maximumReplicas)
{
  _maximumReplica = maximumReplicas;

  int dx = _maximumReplica.x - _minimumReplica.x + 1;

  _fullCell.ax = _unitCell.ax * double(dx);
  _fullCell.ay = _unitCell.ay * double(dx);
  _fullCell.az = _unitCell.az * double(dx);

  int dy = _maximumReplica.y - _minimumReplica.y + 1;

  _fullCell.bx = _unitCell.bx * double(dy);
  _fullCell.by = _unitCell.by * double(dy);
  _fullCell.bz = _unitCell.bz * double(dy);

  int dz = _maximumReplica.z - _minimumReplica.z + 1;

  _fullCell.cx = _unitCell.cx * double(dz);
  _fullCell.cy = _unitCell.cy * double(dz);
  _fullCell.cz = _unitCell.cz * double(dz);

  _inverseFullCell = _fullCell.inverse();
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

  int k1 = ((index % n) / (dz * dy)) + _minimumReplica.x;
  int k2 = ((index % n) / dz) % dy + _minimumReplica.y;
  int k3 = (index % n) % dz + _minimumReplica.z;

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

  int dx = _maximumReplica.x - _minimumReplica.x + 1;

  _fullCell.ax = _unitCell.ax * double(dx);
  _fullCell.ay = _unitCell.ay * double(dx);
  _fullCell.az = _unitCell.az * double(dx);

  _inverseFullCell = _fullCell.inverse();
}

int SKCell::maximumReplicaX() const { return _maximumReplica.x; }

void SKCell::setMaximumReplicaX(const int newValue)
{
  _maximumReplica.x = newValue;

  int dx = _maximumReplica.x - _minimumReplica.x + 1;

  _fullCell.ax = _unitCell.ax * double(dx);
  _fullCell.ay = _unitCell.ay * double(dx);
  _fullCell.az = _unitCell.az * double(dx);

  _inverseFullCell = _fullCell.inverse();
}

int SKCell::minimumReplicaY() const { return _minimumReplica.y; }

void SKCell::setMinimumReplicaY(const int newValue)
{
  _minimumReplica.y = newValue;

  int dy = _maximumReplica.y - _minimumReplica.y + 1;

  _fullCell.bx = _unitCell.bx * double(dy);
  _fullCell.by = _unitCell.by * double(dy);
  _fullCell.bz = _unitCell.bz * double(dy);

  _inverseFullCell = _fullCell.inverse();
}

int SKCell::maximumReplicaY() const { return _maximumReplica.y; }

void SKCell::setMaximumReplicaY(const int newValue)
{
  _maximumReplica.y = newValue;

  int dy = _maximumReplica.y - _minimumReplica.y + 1;

  _fullCell.bx = _unitCell.bx * double(dy);
  _fullCell.by = _unitCell.by * double(dy);
  _fullCell.bz = _unitCell.bz * double(dy);

  _inverseFullCell = _fullCell.inverse();
}

int SKCell::minimumReplicaZ() const { return _minimumReplica.z; }

void SKCell::setMinimumReplicaZ(const int newValue)
{
  _minimumReplica.z = newValue;

  int dz = _maximumReplica.z - _minimumReplica.z + 1;

  _fullCell.cx = _unitCell.cx * double(dz);
  _fullCell.cy = _unitCell.cy * double(dz);
  _fullCell.cz = _unitCell.cz * double(dz);

  _inverseFullCell = _fullCell.inverse();
}

int SKCell::maximumReplicaZ() const { return _maximumReplica.z; }

void SKCell::setMaximumReplicaZ(const int newValue)
{
  _maximumReplica.z = newValue;

  int dz = _maximumReplica.z - _minimumReplica.z + 1;

  _fullCell.cx = _unitCell.cx * double(dz);
  _fullCell.cy = _unitCell.cy * double(dz);
  _fullCell.cz = _unitCell.cz * double(dz);

  _inverseFullCell = _fullCell.inverse();
}

int3 SKCell::numberOfReplicas()
{
  int dx = _maximumReplica.x - _minimumReplica.x + 1;
  int dy = _maximumReplica.y - _minimumReplica.y + 1;
  int dz = _maximumReplica.z - _minimumReplica.z + 1;
  return int3(dx, dy, dz);
}

double SKCell::a() const { return _unitCell[0].length(); }

void SKCell::a(const double& newValue)
{
  double multiplier = newValue / a();
  _unitCell = double3x3(multiplier * _unitCell[0], _unitCell[1], _unitCell[2]);
}

void SKCell::setLengthA(double value)
{
  double multiplier = value / a();
  _unitCell = double3x3(multiplier * _unitCell[0], _unitCell[1], _unitCell[2]);
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
}

double SKCell::b() const { return _unitCell[1].length(); }
void SKCell::b(const double& newValue)
{
  double multiplier = newValue / b();
  _unitCell = double3x3(_unitCell[0], multiplier * _unitCell[1], _unitCell[2]);
}

void SKCell::setLengthB(double value)
{
  double multiplier = value / b();
  _unitCell = double3x3(_unitCell[0], multiplier * _unitCell[1], _unitCell[2]);
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
}

double SKCell::c() const { return _unitCell[2].length(); }

void SKCell::setLengthC(double value)
{
  double multiplier = value / c();
  _unitCell = double3x3(_unitCell[0], _unitCell[1], multiplier * _unitCell[2]);
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
}

void SKCell::c(const double& newValue)
{
  double multiplier = newValue / c();
  _unitCell = double3x3(_unitCell[0], _unitCell[1], multiplier * _unitCell[2]);
}
double SKCell::alpha() const
{
  double3 column2 = _unitCell[1];
  double3 column3 = _unitCell[2];
  double length2 = column2.length();
  double length3 = column3.length();

  return std::acos(double3::dot(column2, column3) / (length2 * length3));
}
void SKCell::alpha(const double& newValue)
{
  double betaAngle = beta();
  double gammaAngle = gamma();
  double temp = (std::cos(newValue) - std::cos(gammaAngle) * std::cos(betaAngle)) / std::sin(gammaAngle);
  double lengthC = c();

  double3 v3 = lengthC * double3(std::cos(betaAngle), temp,
                                 std::sqrt(1.0 - std::cos(betaAngle) * std::cos(betaAngle) - temp * temp));
  _unitCell = double3x3(_unitCell[0], _unitCell[1], v3);
}

void SKCell::setAlphaAngle(double value)
{
  double betaAngle = beta();
  double gammaAngle = gamma();
  double temp = (std::cos(value) - std::cos(gammaAngle) * std::cos(betaAngle)) / std::sin(gammaAngle);
  double lengthC = c();

  double3 v3 = lengthC * double3(std::cos(betaAngle), temp,
                                 std::sqrt(1.0 - std::cos(betaAngle) * std::cos(betaAngle) - temp * temp));
  _unitCell = double3x3(_unitCell[0], _unitCell[1], v3);

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
}

double SKCell::beta() const
{
  double3 column1 = _unitCell[0];
  double3 column3 = _unitCell[2];
  double length1 = column1.length();
  double length3 = column3.length();

  return std::acos(double3::dot(column1, column3) / (length1 * length3));
}

void SKCell::setBetaAngle(double value)
{
  double alphaAngle = alpha();
  double gammaAngle = gamma();
  double temp = (std::cos(alphaAngle) - std::cos(gammaAngle) * std::cos(value)) / std::sin(gammaAngle);
  double lengthC = c();

  double3 v3 =
      lengthC * double3(std::cos(value), temp, std::sqrt(1.0 - std::cos(value) * std::cos(value) - temp * temp));
  _unitCell = double3x3(_unitCell[0], _unitCell[1], v3);

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
}

void SKCell::beta(const double& newValue)
{
  double alphaAngle = alpha();
  double gammaAngle = gamma();
  double temp = (std::cos(alphaAngle) - std::cos(gammaAngle) * std::cos(newValue)) / std::sin(gammaAngle);
  double lengthC = c();

  double3 v3 = lengthC * double3(std::cos(newValue), temp,
                                 std::sqrt(1.0 - std::cos(newValue) * std::cos(newValue) - temp * temp));
  _unitCell = double3x3(_unitCell[0], _unitCell[1], v3);
}
double SKCell::gamma() const
{
  double3 column1 = _unitCell[0];
  double3 column2 = _unitCell[1];
  double length1 = column1.length();
  double length2 = column2.length();

  return std::acos(double3::dot(column1, column2) / (length1 * length2));
}

void SKCell::setGammaAngle(double value)
{
  double alphaAngle = alpha();
  double betaAngle = beta();
  double lengthB = b();
  double lengthC = c();
  double temp = (std::cos(alphaAngle) - std::cos(value) * std::cos(betaAngle)) / std::sin(value);

  double3 v2 = lengthB * double3(std::cos(value), std::sin(value), 0.0);
  double3 v3 = lengthC * double3(std::cos(betaAngle), temp,
                                 std::sqrt(1.0 - std::cos(betaAngle) * std::cos(betaAngle) - temp * temp));
  _unitCell = double3x3(_unitCell[0], v2, v3);

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
}

void SKCell::gamma(const double& newValue)
{
  double alphaAngle = alpha();
  double betaAngle = beta();
  double lengthB = b();
  double lengthC = c();
  double temp = (std::cos(alphaAngle) - std::cos(newValue) * std::cos(betaAngle)) / std::sin(newValue);

  double3 v2 = lengthB * double3(std::cos(newValue), std::sin(newValue), 0.0);
  double3 v3 = lengthC * double3(std::cos(betaAngle), temp,
                                 std::sqrt(1.0 - std::cos(betaAngle) * std::cos(betaAngle) - temp * temp));
  _unitCell = double3x3(_unitCell[0], v2, v3);
}

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

  return double3::dot(column1, v2);
}

double3 SKCell::perpendicularWidths() const
{
  double3 column1 = _unitCell[0];
  double3 column2 = _unitCell[1];
  double3 column3 = _unitCell[2];

  double3 v1 = double3::cross(column1, column2);
  double3 v2 = double3::cross(column2, column3);
  double3 v3 = double3::cross(column3, column1);

  double volume = double3::dot(column1, v2);

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

  double volume = double3::dot(column1, v2);

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
