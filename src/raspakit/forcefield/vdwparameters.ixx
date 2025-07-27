module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <optional>
#include <ostream>
#include <string>
#include <vector>
#endif

export module vdwparameters;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import double4;
import units;

/**
 * \brief Represents the van der Waals parameters for particle interactions.
 *
 * The VDWParameters struct encapsulates parameters for different types of van der Waals potentials,
 * including Lennard-Jones, Buckingham, Morse, Feynmann-Hibbs, MM3, and Born-Huggins-Meyer potentials.
 * It includes the potential parameters, a shift value calculated at the cutoff distance,
 * tail correction energy, and the type of potential.
 */
export struct VDWParameters
{
  /**
   * \brief Enumeration of van der Waals potential types.
   */
  enum class Type : std::size_t
  {
    None = 0,
    LennardJones = 1,
    BuckingHam = 2,
    Morse = 3,
    FeynmannHibbs = 4,
    MM3 = 5,
    BornHugginsMeyer = 6,
    RepulsiveHarmonic = 100
  };

  double4 parameters;             ///< The potential parameters. For LJ: epsilon, sigma; for Buckingham: 3 parameters.
  double shift;                   ///< The potential energy shift calculated at the cutoff distance.
  double tailCorrectionEnergy;    ///< The tail correction energy for the potential.
  double tailCorrectionPressure;  ///< The tail correction energy for the potential.
  Type type{0};                   ///< The type of van der Waals potential.

  /**
   * \brief Default constructor for VDWParameters.
   *
   * Initializes the parameters with zero values.
   */
  VDWParameters()
      : parameters(0.0, 0.0, 0.0, 0.0),
        shift(0.0),
        tailCorrectionEnergy(0.0),
        tailCorrectionPressure(0.0),
        type(Type::LennardJones)
  {
  }

  VDWParameters(double4 parameters, double shift, double tailCorrectionEnergy, double tailCorrectionPressure, Type type)
      : parameters(parameters),
        shift(shift),
        tailCorrectionEnergy(tailCorrectionEnergy),
        tailCorrectionPressure(tailCorrectionPressure),
        type(type)
  {
  }

  /**
   * \brief Constructs a Lennard-Jones VDWParameter structure.
   *
   * Initializes the VDWParameters with the given epsilon and sigma values for a Lennard-Jones potential.
   *
   * \param epsilon The strength parameter of the potential in units of Kelvin.
   * \param sigma The size parameter of the potential in units of Angstrom.
   */
  VDWParameters(double epsilon, double sigma)
      : parameters(double4(epsilon * Units::KelvinToEnergy, sigma, 0.0, 0.0)),
        shift(0.0),
        tailCorrectionEnergy(0.0),
        tailCorrectionPressure(0.0),
        type(Type::LennardJones)
  {
  }

  VDWParameters(double epsilon, double sigma, Type type)
      : parameters(double4(epsilon * Units::KelvinToEnergy, sigma, 0.0, 0.0)),
        shift(0.0),
        tailCorrectionEnergy(0.0),
        tailCorrectionPressure(0.0),
        type(type)
  {
  }

  /**
   * \brief Returns interaction type from given string
   *
   * \param interactionType string to convert to interaction type
   */
  static Type stringToEnum(std::string interactionType);

  /**
   * \brief Computes the potential energy shift at the cutoff distance.
   *
   * Calculates the shift in potential energy at the specified cutoff distance to ensure continuity
   * of the potential.
   *
   * \param cutOff The cutoff distance at which to compute the shift.
   */
  void computeShiftAtCutOff(double cutOff)
  {
    shift = 0.0;
    double scaling = 1.0;
    double arg1 = parameters.x;
    double arg2 = parameters.y * parameters.y;
    double rr = cutOff * cutOff;
    double temp = (rr / arg2);
    double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
    shift = scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)));
  }

  double sizeParameter() const
  {
    switch (type)
    {
      case Type::LennardJones:
        return parameters.y;
      default:
        return 0.0;
    }
  };

  bool operator==(const VDWParameters &other) const;
  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VDWParameters &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VDWParameters &p);
};
