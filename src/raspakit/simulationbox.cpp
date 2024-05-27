module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <ostream>
#include <source_location>
#include <sstream>
#include <string>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#endif

module simulationbox;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <numbers>;
import <string>;
import <iostream>;
import <ostream>;
import <sstream>;
import <fstream>;
import <exception>;
import <source_location>;
import <complex>;
import <array>;
import <map>;
import <algorithm>;
#if defined(__has_include) && __has_include(<print>)
import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
import print;
#endif

import randomnumbers;
import double3x3;
import double3;
import archive;
import units;
import stringutils;
import hdf5;

SimulationBox::SimulationBox(double a, double b, double c, Type type) : type(type)
{
  lengthA = a;
  lengthB = b;
  lengthC = c;
  angleAlpha = std::numbers::pi * 90.0 / 180.0;
  angleBeta = std::numbers::pi * 90.0 / 180.0;
  angleGamma = std::numbers::pi * 90.0 / 180.0;

  double3 v1 = double3(a, 0.0, 0.0);
  double3 v2 = double3(0.0, b, 0.0);
  double3 v3 = double3(0.0, 0.0, c);
  cell = double3x3(v1, v2, v3);
  if (a != 0.0 && b != 0.0 && c != 0.0)
  {
    inverseCell = cell.inverse();
    volume = cell.determinant();
  }
  else
  {
    volume = 0.0;
  }
}

SimulationBox::SimulationBox(double a, double b, double c, double alpha, double beta, double gamma, Type type)
    : lengthA(a), lengthB(b), lengthC(c), angleAlpha(alpha), angleBeta(beta), angleGamma(gamma), type(type)
{
  double temp = (cos(alpha) - cos(gamma) * cos(beta)) / sin(gamma);

  double3 v1 = double3(a, 0.0, 0.0);
  double3 v2 = double3(b * cos(gamma), b * sin(gamma), 0.0);
  double3 v3 = double3(c * cos(beta), c * temp, c * sqrt(1.0 - cos(beta) * cos(beta) - temp * temp));
  cell = double3x3(v1, v2, v3);
  if (a != 0.0 && b != 0.0 && c != 0.0)
  {
    inverseCell = cell.inverse();
    volume = cell.determinant();
  }
  else
  {
    volume = 0.0;
  }
}

SimulationBox::SimulationBox(double3x3 m, Type type) : type(type)
{
  this->cell = m;
  this->inverseCell = m.inverse();
  this->volume = m.determinant();

  double3 column1 = cell[0];
  double3 column2 = cell[1];
  double3 column3 = cell[2];
  this->lengthA = column1.length();
  this->lengthB = column2.length();
  this->lengthC = column3.length();

  angleAlpha = std::acos(double3::dot(column2, column3) / (this->lengthB * this->lengthC));
  angleBeta = std::acos(double3::dot(column1, column3) / (this->lengthA * this->lengthC));
  angleGamma = std::acos(double3::dot(column1, column2) / (this->lengthA * this->lengthB));
}

double3 SimulationBox::randomPosition(RandomNumber &random) const
{
  return cell * double3(random.uniform(), random.uniform(), random.uniform());
}

double3 SimulationBox::perpendicularWidths() const
{
  double3 c1 = double3::cross(cell[0], cell[1]);
  double3 c2 = double3::cross(cell[1], cell[2]);
  double3 c3 = double3::cross(cell[0], cell[2]);

  double v = std::fabs(cell[0].x * c2.x + cell[0].y * c2.y + cell[0].z * c2.z);

  // calculate cell perpendicular widths
  return double3(v / std::sqrt(double3::dot(c2, c2)), v / std::sqrt(double3::dot(c3, c3)),
                 v / std::sqrt(double3::dot(c1, c1)));
}

void SimulationBox::setBoxLengths(double3 lengths)
{
  lengthA = lengths.x;
  lengthB = lengths.y;
  lengthC = lengths.z;
  double temp = (cos(angleAlpha) - cos(angleGamma) * cos(angleBeta)) / sin(angleGamma);
  double3 v1 = double3(lengths.x, 0.0, 0.0);
  double3 v2 = double3(lengths.y * cos(angleGamma), lengths.y * sin(angleGamma), 0.0);
  double3 v3 = double3(lengths.z * cos(angleBeta), lengths.z * temp,
                       lengths.z * sqrt(1.0 - cos(angleBeta) * cos(angleBeta) - temp * temp));
  cell = double3x3(v1, v2, v3);
  inverseCell = cell.inverse();
  volume = cell.determinant();
}

void SimulationBox::setBoxAngles(double3 angles)
{
  angleAlpha = angles.x;
  angleBeta = angles.y;
  angleGamma = angles.z;
  double3 lengths = this->lengths();
  lengthA = lengths.x;
  lengthB = lengths.y;
  lengthC = lengths.z;
  double temp = (cos(angles.x) - cos(angles.z) * cos(angles.y)) / sin(angles.z);
  double3 v1 = double3(lengths.x, 0.0, 0.0);
  double3 v2 = double3(lengths.y * cos(angles.z), lengths.y * sin(angles.z), 0.0);
  double3 v3 = double3(lengths.z * cos(angles.y), lengths.z * temp,
                       lengths.z * sqrt(1.0 - cos(angles.y) * cos(angles.y) - temp * temp));
  cell = double3x3(v1, v2, v3);
  inverseCell = cell.inverse();
  volume = cell.determinant();
}

double3 SimulationBox::lengths() { return double3(cell[0].length(), cell[1].length(), cell[2].length()); }

double3 SimulationBox::angles()
{
  double3 column1 = cell[0];
  double3 column2 = cell[1];
  double3 column3 = cell[2];
  double length1 = column1.length();
  double length2 = column2.length();
  double length3 = column3.length();

  return double3(acos(double3::dot(column2, column3) / (length2 * length3)),
                 acos(double3::dot(column1, column3) / (length1 * length3)),
                 acos(double3::dot(column1, column2) / (length1 * length2)));
}

std::string SimulationBox::printStatus() const
{
  std::ostringstream stream;

  std::print(stream, "Box:     {:9.5f} {:9.5f} {:9.5f}\n", cell.ax, cell.bx, cell.cx);
  std::print(stream, "         {:9.5f} {:9.5f} {:9.5f}\n", cell.ay, cell.by, cell.cy);
  std::print(stream, "         {:9.5f} {:9.5f} {:9.5f}\n", cell.az, cell.bz, cell.cz);
  std::print(stream, "Lengths: {:9.5f} {:9.5f} {:9.5f}\n", lengthA, lengthB, lengthC);
  double conv = 180.0 / std::numbers::pi;
  std::print(stream, "Angles:  {:9.5f} {:9.5f} {:9.5f}\n", conv * angleAlpha, conv * angleBeta, conv * angleGamma);

  return stream.str();
}

void SimulationBox::logStatus(HDF5Writer &hdf5) const
{
  hdf5.createDataset<double>("initial_conditions", "box", std::vector<size_t>{3, 3}, {{"dimensions", "(xyz, xyz)"}});
  hdf5.writeVector(
      "initial_conditions", "box",
      std::vector<double>{cell.ax, cell.bx, cell.cx, cell.ay, cell.by, cell.cy, cell.az, cell.bz, cell.cz});
  hdf5.createDataset<double>("initial_conditions", "boxLengths", std::vector<size_t>{3}, {{"dimensions", "(xyz, )"}});
  hdf5.writeVector("initial_conditions", "boxLengths", std::vector<double>{lengthA, lengthB, lengthC});
  hdf5.createDataset<double>("initial_conditions", "boxAngles", std::vector<size_t>{3}, {{"dimensions", "(xyz, )"}});
}

std::string SimulationBox::printStatus(const SimulationBox &average, [[maybe_unused]] const SimulationBox &error) const
{
  std::ostringstream stream;
  std::print(stream, "Box:     {:9.5f} {:9.5f} {:9.5f}  Average: {:9.5f} {:9.5f} {:9.5f}\n", cell.ax, cell.bx, cell.cx,
             average.cell.ax, average.cell.bx, average.cell.cx);
  std::print(stream, "         {:9.5f} {:9.5f} {:9.5f}           {:9.5f} {:9.5f} {:9.5f}\n", cell.ay, cell.by, cell.cy,
             average.cell.ay, average.cell.by, average.cell.cy);
  std::print(stream, "         {:9.5f} {:9.5f} {:9.5f}           {:9.5f} {:9.5f} {:9.5f}\n", cell.az, cell.bz, cell.cz,
             average.cell.az, average.cell.bz, average.cell.cz);
  std::print(stream, "Lengths: {:9.5f} {:9.5f} {:9.5f}  Average: {:9.5f} {:9.5f} {:9.5f}\n", lengthA, lengthB, lengthC,
             average.lengthA, average.lengthB, average.lengthC);

  double conv = 180.0 / std::numbers::pi;
  std::print(stream, "Angles:  {:9.5f} {:9.5f} {:9.5f}  Average: {:9.5f} {:9.5f} {:9.5f}\n", conv * angleAlpha,
             conv * angleBeta, conv * angleGamma, conv * average.angleAlpha, conv * average.angleBeta,
             conv * average.angleGamma);

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const SimulationBox &box)
{
  archive << box.versionNumber;

  archive << box.lengthA;
  archive << box.lengthB;
  archive << box.lengthC;
  archive << box.angleAlpha;
  archive << box.angleBeta;
  archive << box.angleGamma;
  archive << box.cell;
  archive << box.inverseCell;
  archive << box.volume;
  archive << box.type;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, SimulationBox &box)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > box.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'SimulationBox' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> box.lengthA;
  archive >> box.lengthB;
  archive >> box.lengthC;
  archive >> box.angleAlpha;
  archive >> box.angleBeta;
  archive >> box.angleGamma;
  archive >> box.cell;
  archive >> box.inverseCell;
  archive >> box.volume;
  archive >> box.type;

  return archive;
}
