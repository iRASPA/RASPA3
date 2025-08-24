module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iostream>
#include <istream>
#include <map>
#include <numbers>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#endif

export module simulationbox;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;

import double3;
import int3;
import double3x3;
import randomnumbers;
import units;
import json;

#if defined(__GNUC__)
#define ALWAYS_INLINE __attribute__((__always_inline__))
#elif defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#endif

/**
 * \brief Represents a simulation box used in simulations.
 *
 * The SimulationBox struct defines the dimensions and angles of the simulation box.
 * It supports different types, such as rectangular and triclinic, and contains
 * information about the cell matrix and its inverse, as well as the volume of the box.
 */
export struct SimulationBox
{
  /**
   * \brief Enumeration for the type of the simulation box.
   */
  enum class Type : int
  {
    Rectangular = 0,
    Triclinic = 1
  };

  /**
   * \brief Default constructor for the SimulationBox class.
   *
   * Initializes a SimulationBox with all lengths and angles set to zero,
   * and sets the cell matrices to zero matrices.
   */
  SimulationBox()
      : lengthA(0.0),
        lengthB(0.0),
        lengthC(0.0),
        angleAlpha(0.0),
        angleBeta(0.0),
        angleGamma(0.0),
        cell(double3x3(double3(0.0, 0.0, 0.0), double3(0.0, 0.0, 0.0), double3(0.0, 0.0, 0.0))),
        inverseCell(double3x3(double3(0.0, 0.0, 0.0), double3(0.0, 0.0, 0.0), double3(0.0, 0.0, 0.0))),
        volume(0.0) {};

  bool operator==(SimulationBox const&) const = default;

  /**
   * \brief Constructs a rectangular SimulationBox with given lengths.
   *
   * Initializes a rectangular SimulationBox with the specified lengths for the sides.
   * Sets the angles to 90 degrees and computes the cell matrices and volume.
   *
   * \param a Length of side A.
   * \param b Length of side B.
   * \param c Length of side C.
   */
  explicit SimulationBox(double a, double b, double c);

  /**
   * \brief Constructs a SimulationBox with given lengths and angles.
   *
   * Initializes a SimulationBox with the specified lengths for the sides and angles.
   * Computes the cell matrices and volume based on the provided dimensions.
   *
   * \param a Length of side A.
   * \param b Length of side B.
   * \param c Length of side C.
   * \param alpha Angle between sides B and C (in radians).
   * \param beta Angle between sides A and C (in radians).
   * \param gamma Angle between sides A and B (in radians).
   */
  explicit SimulationBox(double a, double b, double c, double alpha, double beta, double gamma);
  /**
   * \brief Constructs a SimulationBox with given lengths, angles, and type.
   *
   * Initializes a SimulationBox with the specified lengths, angles, and box type.
   * Computes the cell matrices and volume based on the provided dimensions.
   *
   * \param a Length of side A.
   * \param b Length of side B.
   * \param c Length of side C.
   * \param alpha Angle between sides B and C (in radians).
   * \param beta Angle between sides A and C (in radians).
   * \param gamma Angle between sides A and B (in radians).
   * \param type The type of the simulation box (Rectangular or Triclinic).
   */
  explicit SimulationBox(double a, double b, double c, double alpha, double beta, double gamma, Type type);

  /**
   * \brief Constructs a SimulationBox with a given cell matrix.
   *
   * Initializes a SimulationBox with the specified cell matrix.
   * Computes the lengths and angles from the cell matrix, as well as the inverse cell matrix and volume.
   *
   * \param m The cell matrix defining the box dimensions and angles.
   */
  explicit SimulationBox(double3x3 m);
  /**
   * \brief Constructs a SimulationBox with a given cell matrix and type.
   *
   * Initializes a SimulationBox with the specified cell matrix and box type.
   * Computes the lengths and angles from the cell matrix, as well as the inverse cell matrix and volume.
   *
   * \param m The cell matrix defining the box dimensions and angles.
   * \param type The type of the simulation box (Rectangular or Triclinic).
   */
  explicit SimulationBox(double3x3 m, Type type);

  /**
   * \brief Applies periodic boundary conditions to a displacement vector.
   *
   * Adjusts the displacement vector \p dr to account for periodic boundary conditions,
   * ensuring that it lies within the primary simulation cell.
   * For rectangular boxes, the method is optimized for performance.
   *
   * \param dr The displacement vector to be adjusted.
   * \return The adjusted displacement vector within the primary simulation cell.
   */
  ALWAYS_INLINE inline double3 applyPeriodicBoundaryConditions(const double3& dr) const
  {
    switch (type)
    {
      case SimulationBox::Type::Rectangular:
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
        /*
        s.x = dr.x - cell.ax * std::rint(dr.x * inverseCell.ax);
        s.y = dr.y - cell.by * std::rint(dr.y * inverseCell.by);
        s.z = dr.z - cell.cz * std::rint(dr.z * inverseCell.cz);
        */
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

  /**
   * \brief Sets the lengths of the simulation box.
   *
   * Updates the lengths of the box and recalculates the cell matrix,
   * inverse cell matrix, and volume based on the new lengths.
   *
   * \param lengths A double3 containing the new lengths (lengthA, lengthB, lengthC).
   */
  void setBoxLengths(double3 lengths);
  /**
   * \brief Sets the angles of the simulation box.
   *
   * Updates the angles of the box and recalculates the cell matrix,
   * inverse cell matrix, and volume based on the new angles.
   * Angles are provided in radians.
   *
   * \param angles A double3 containing the new angles (angleAlpha, angleBeta, angleGamma).
   */
  void setBoxAngles(double3 angles);
  /**
   * \brief Retrieves the lengths of the simulation box.
   *
   * Returns the lengths of the sides of the simulation box.
   *
   * \return A double3 containing the lengths (lengthA, lengthB, lengthC).
   */
  double3 lengths() const;
  /**
   * \brief Retrieves the angles of the simulation box.
   *
   * Returns the angles between the sides of the simulation box.
   *
   * \return A double3 containing the angles (angleAlpha, angleBeta, angleGamma) in radians.
   */
  double3 angles() const;
  /**
   * \brief Generates a random position within the simulation box.
   *
   * Generates a random position within the simulation box using the provided random number generator.
   *
   * \param random A random number generator.
   * \return A double3 representing a random position within the box.
   */
  double3 randomPosition(RandomNumber& random) const;
  /**
   * \brief Calculates the perpendicular widths of the simulation box.
   *
   * Computes the perpendicular widths (shortest distances) along each axis of the box.
   *
   * \return A double3 containing the perpendicular widths along x, y, and z axes.
   */
  double3 perpendicularWidths() const;

  int3 smallestNumberOfUnitCellsForMinimumImagesConvention(double cutOff) const;

  /**
   * \brief Returns a string representation of the simulation box.
   *
   * Generates a string that includes the cell matrix, lengths, angles, and boundary conditions type.
   *
   * \return A string representing the simulation box.
   */
  std::string printStatus() const;
  /**
   * \brief Returns a string representation of the simulation box with averages and errors.
   *
   * Generates a string that includes the cell matrix, lengths, angles, and boundary conditions type,
   * along with average values and errors for comparison.
   *
   * \param average The average simulation box for comparison.
   * \param error The error simulation box for comparison.
   * \return A string representing the simulation box with averages and errors.
   */
  std::string printStatus(const SimulationBox& average, const SimulationBox& error) const;

  /**
   * \brief Returns a JSON representation of the simulation box.
   *
   * Generates a JSON object containing the cell matrix, lengths, and angles of the box.
   *
   * \return A JSON object representing the simulation box.
   */
  nlohmann::json jsonStatus() const;

  std::uint64_t versionNumber{1};
  double lengthA;
  double lengthB;
  double lengthC;
  double angleAlpha;
  double angleBeta;
  double angleGamma;
  double3x3 cell;
  double3x3 inverseCell;
  double volume;
  Type type;

  inline SimulationBox& operator+=(const SimulationBox& b)
  {
    lengthA += b.lengthA;
    lengthB += b.lengthB;
    lengthC += b.lengthC;
    angleAlpha += b.angleAlpha;
    angleBeta += b.angleBeta;
    angleGamma += b.angleGamma;
    cell += b.cell;
    volume += volume;

    return *this;
  }

  /**
   * \brief Scales the simulation box by a scalar factor.
   *
   * Multiplies the lengths of the box by the given scalar and recalculates
   * the cell matrix, inverse cell matrix, and volume.
   *
   * \param scale The scalar factor to scale the box by.
   */
  void scale(double scale)
  {
    lengthA *= scale;
    lengthB *= scale;
    lengthC *= scale;

    double temp = (std::cos(angleAlpha) - std::cos(angleGamma) * std::cos(angleBeta)) / std::sin(angleGamma);
    double3 v1 = double3(lengthA, 0.0, 0.0);
    double3 v2 = double3(lengthB * std::cos(angleGamma), lengthB * std::sin(angleGamma), 0.0);
    double3 v3 = double3(lengthC * std::cos(angleBeta), lengthC * temp,
                         lengthC * std::sqrt(1.0 - std::cos(angleBeta) * std::cos(angleBeta) - temp * temp));
    cell = double3x3(v1, v2, v3);
    inverseCell = cell.inverse();
    volume = cell.determinant();
  }

  /**
   * \brief Scales the simulation box by integer factors along each axis.
   *
   * Multiplies the lengths of the box by the given integer factors along x, y, and z axes,
   * and recalculates the cell matrix, inverse cell matrix, and volume.
   *
   * \param scale An int3 containing the scaling factors along each axis.
   */
  void scale(int3 scale)
  {
    lengthA *= static_cast<double>(scale.x);
    lengthB *= static_cast<double>(scale.y);
    lengthC *= static_cast<double>(scale.z);

    double temp = (std::cos(angleAlpha) - std::cos(angleGamma) * std::cos(angleBeta)) / std::sin(angleGamma);
    double3 v1 = double3(lengthA, 0.0, 0.0);
    double3 v2 = double3(lengthB * std::cos(angleGamma), lengthB * std::sin(angleGamma), 0.0);
    double3 v3 = double3(lengthC * std::cos(angleBeta), lengthC * temp,
                         lengthC * std::sqrt(1.0 - std::cos(angleBeta) * std::cos(angleBeta) - temp * temp));
    cell = double3x3(v1, v2, v3);
    inverseCell = cell.inverse();
    volume = cell.determinant();
  }

  /**
   * \brief Returns a new SimulationBox scaled by a scalar factor.
   *
   * Creates a new SimulationBox with lengths multiplied by the given scalar,
   * and recalculates the cell matrix, inverse cell matrix, and volume.
   *
   * \param scale The scalar factor to scale the box by.
   * \return A new scaled SimulationBox.
   */
  SimulationBox scaled(double scale)
  {
    SimulationBox v;

    v.lengthA = scale * lengthA;
    v.lengthB = scale * lengthB;
    v.lengthC = scale * lengthC;
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

  /**
   * \brief Returns a new SimulationBox scaled by integer factors along each axis.
   *
   * Creates a new SimulationBox with lengths multiplied by the given integer factors,
   * and recalculates the cell matrix, inverse cell matrix, and volume.
   *
   * \param scale An int3 containing the scaling factors along each axis.
   * \return A new scaled SimulationBox.
   */
  SimulationBox scaled(int3 scale) const
  {
    SimulationBox v;

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

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const SimulationBox& box);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, SimulationBox& box);
  friend void to_json(nlohmann::json&, const SimulationBox&);
  friend void from_json(const nlohmann::json&, SimulationBox&);
};

export inline SimulationBox operator+(const SimulationBox& a, const SimulationBox& b)
{
  SimulationBox m;

  m.lengthA = a.lengthA + b.lengthA;
  m.lengthB = a.lengthB + b.lengthB;
  m.lengthC = a.lengthC + b.lengthC;
  m.angleAlpha = a.angleAlpha + b.angleAlpha;
  m.angleBeta = a.angleBeta + b.angleBeta;
  m.angleGamma = a.angleGamma + b.angleGamma;
  m.cell = a.cell + b.cell;
  m.volume = a.volume + b.volume;

  return m;
}

export inline SimulationBox operator-(const SimulationBox& a, const SimulationBox& b)
{
  SimulationBox m;

  m.lengthA = a.lengthA - b.lengthA;
  m.lengthB = a.lengthB - b.lengthB;
  m.lengthC = a.lengthC - b.lengthC;
  m.angleAlpha = a.angleAlpha - b.angleAlpha;
  m.angleBeta = a.angleBeta - b.angleBeta;
  m.angleGamma = a.angleGamma - b.angleGamma;
  m.cell = a.cell - b.cell;
  m.volume = a.volume - b.volume;

  return m;
}

export inline SimulationBox operator*(const SimulationBox& a, const SimulationBox& b)
{
  SimulationBox m;

  m.lengthA = a.lengthA * b.lengthA;
  m.lengthB = a.lengthB * b.lengthB;
  m.lengthC = a.lengthC * b.lengthC;
  m.angleAlpha = a.angleAlpha * b.angleAlpha;
  m.angleBeta = a.angleBeta * b.angleBeta;
  m.angleGamma = a.angleGamma * b.angleGamma;
  m.cell = a.cell * b.cell;
  m.volume = a.volume * b.volume;
  return m;
}

export inline SimulationBox operator*(const double& a, const SimulationBox& b)
{
  SimulationBox m;

  m.lengthA = a * b.lengthA;
  m.lengthB = a * b.lengthB;
  m.lengthC = a * b.lengthC;
  m.angleAlpha = a * b.angleAlpha;
  m.angleBeta = a * b.angleBeta;
  m.angleGamma = a * b.angleGamma;
  m.cell = a * b.cell;
  m.volume = a * b.volume;
  return m;
}

export inline SimulationBox operator/(const SimulationBox& a, const double& b)
{
  SimulationBox m;

  double temp = 1.0 / b;
  m.lengthA = a.lengthA * temp;
  m.lengthB = a.lengthB * temp;
  m.lengthC = a.lengthC * temp;
  m.angleAlpha = a.angleAlpha * temp;
  m.angleBeta = a.angleBeta * temp;
  m.angleGamma = a.angleGamma * temp;
  m.cell = a.cell * temp;
  m.volume = a.volume * temp;
  return m;
}

/**
 * \brief Computes the square root of each component of a SimulationBox.
 *
 * Returns a new SimulationBox where each component is the square root of the corresponding component of \p a.
 *
 * \param a The SimulationBox.
 * \return A new SimulationBox with square-rooted components.
 */
export inline SimulationBox sqrt(const SimulationBox& a)
{
  SimulationBox m;

  m.lengthA = std::sqrt(a.lengthA);
  m.lengthB = std::sqrt(a.lengthB);
  m.lengthC = std::sqrt(a.lengthC);
  m.angleAlpha = std::sqrt(a.angleAlpha);
  m.angleBeta = std::sqrt(a.angleBeta);
  m.angleGamma = std::sqrt(a.angleGamma);
  m.cell = sqrt(a.cell);
  m.volume = std::sqrt(a.volume);

  return m;
}

/**
 * \brief Computes the component-wise maximum of two SimulationBoxes.
 *
 * Returns a new SimulationBox where each component is the maximum of the corresponding components of \p a and \p b.
 *
 * \param a The first SimulationBox.
 * \param b The second SimulationBox.
 * \return A new SimulationBox with maximum components.
 */
export inline SimulationBox max(const SimulationBox& a, const SimulationBox& b)
{
  SimulationBox c = SimulationBox(std::max(a.lengthA, b.lengthA), std::max(a.lengthB, b.lengthB),
                                  std::max(a.lengthC, b.lengthC), std::max(a.angleAlpha, b.angleAlpha),
                                  std::max(a.angleBeta, b.angleBeta), std::max(a.angleGamma, b.angleGamma));

  if (a.type == SimulationBox::Type::Triclinic || b.type == SimulationBox::Type::Triclinic)
  {
    c.type = SimulationBox::Type::Triclinic;
  }

  return c;
}

export inline SimulationBox operator*(const double3x3& a, const SimulationBox& b) { return SimulationBox(a * b.cell); }

export void to_json(nlohmann::json& j, const SimulationBox& b)
{
  j = nlohmann::json{{"length-a", b.lengthA},
                     {"length-b", b.lengthB},
                     {"length-c", b.lengthC},
                     {"angle-alpha", b.angleAlpha * Units::RadiansToDegrees},
                     {"angle-beta", b.angleBeta * Units::RadiansToDegrees},
                     {"angle-gamma", b.angleGamma * Units::RadiansToDegrees}};
}

export void from_json(const nlohmann::json& j, SimulationBox& b)
{
  j.at("length-a").get_to(b.lengthA);
  j.at("length-b").get_to(b.lengthB);
  j.at("length-c").get_to(b.lengthC);
  j.at("angle-alpha").get_to(b.angleAlpha);
  j.at("angle-beta").get_to(b.angleBeta);
  j.at("angle-gamma").get_to(b.angleGamma);
  b.angleAlpha *= Units::DegreesToRadians;
  b.angleBeta *= Units::DegreesToRadians;
  b.angleGamma *= Units::DegreesToRadians;

  double temp = (std::cos(b.angleAlpha) - std::cos(b.angleGamma) * std::cos(b.angleBeta)) / std::sin(b.angleGamma);

  double3 v1 = double3(b.lengthA, 0.0, 0.0);
  double3 v2 = double3(b.lengthB * std::cos(b.angleGamma), b.lengthB * std::sin(b.angleGamma), 0.0);
  double3 v3 = double3(b.lengthC * std::cos(b.angleBeta), b.lengthC * temp,
                       b.lengthC * std::sqrt(1.0 - std::cos(b.angleBeta) * std::cos(b.angleBeta) - temp * temp));
  b.cell = double3x3(v1, v2, v3);
  if (b.lengthA != 0.0 && b.lengthB != 0.0 && b.lengthC != 0.0)
  {
    b.inverseCell = b.cell.inverse();
    b.volume = b.cell.determinant();
  }
  else
  {
    b.volume = 0.0;
  }
}
