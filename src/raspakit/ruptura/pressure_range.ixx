module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <vector>
#endif

export module pressure_range;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;

export struct PressureRange
{
  enum class Scale : std::size_t
  {
    Log = 0,
    Linear = 1
  };

  bool operator==(PressureRange const &) const = default;

  double pressureStart{1e-2};
  double pressureEnd{1e7};
  std::size_t numberOfPoints{100};
  Scale scale{Scale::Log};

  inline std::vector<double> pressures() const
  {
    std::vector<double> p(numberOfPoints);
    switch (scale)
    {
      case Scale::Log:
      default:
        for (std::size_t i = 0; i < numberOfPoints; ++i)
        {
          p[i] = std::pow(
              10, std::log10(pressureStart) + ((std::log10(pressureEnd) - std::log10(pressureStart)) *
                                               (static_cast<double>(i) / static_cast<double>(numberOfPoints - 1))));
        }
        break;
      case Scale::Linear:
        for (std::size_t i = 0; i < numberOfPoints; ++i)
        {
          p[i] = pressureStart +
                 (pressureEnd - pressureStart) * (static_cast<double>(i) / static_cast<double>(numberOfPoints - 1));
        }
        break;
    }
    return p;
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PressureRange &r);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PressureRange &r);
};
