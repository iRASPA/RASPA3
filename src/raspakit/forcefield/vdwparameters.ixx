module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <ostream>
#include <fstream>
#include <optional>
#endif

export module vdwparameters;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <string>;
import <algorithm>;
import <iostream>;
import <ostream>;
import <fstream>;
import <optional>;
#endif

import archive;
import double4;

export struct VDWParameters
{
  enum class Type : int
  {
    LennardJones = 0,
    BuckingHam = 1,
    Morse = 2,
    FeynmannHibbs = 3,
    MM3 = 4,
    BornHugginsMeyer = 5
  };

  double4 parameters;      // for LJ: epsilon, sigma, for Buckingham: 3 parameters
  double shift;
  double tailCorrectionEnergy;
  Type type{ 0 };

  VDWParameters(): parameters(double4(158.5/ 1.2027242847,3.72,0.0,0.0)), shift(-0.56217796) {}

  VDWParameters(double epsilon, double sigma) : parameters(double4(epsilon, sigma, 0.0, 0.0)), shift(0.0), 
                                                tailCorrectionEnergy(0.0), type(Type::LennardJones)
  {
  }

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

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VDWParameters &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VDWParameters &p);
};

