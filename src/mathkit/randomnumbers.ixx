module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <tuple>
#include <utility>
#endif

export module randomnumbers;

#ifndef USE_LEGACY_HEADERS
import <tuple>;
import <random>;
import <cmath>;
import <utility>;
import <fstream>;
import <iostream>;
import <optional>;
#endif

import archive;
import double3;
import double3x3;
import simd_quatd;

export struct RandomNumber
{
  RandomNumber(std::optional<size_t> s)
  {
    std::random_device rd;
    seed = s.has_value() ? s.value() : rd();
    mt = std::mt19937_64(seed);
    uniformDistribution = std::uniform_real_distribution<double>(0.0, 1.0);
    normalDistribution = std::normal_distribution<double>();
  }

  bool operator==(RandomNumber const &rhs) const
  {
    return (mt == rhs.mt) && (seed == rhs.seed) && (count == rhs.count);
  }

  std::mt19937_64 mt;
  size_t seed{1400};
  size_t count{0};
  std::uniform_real_distribution<double> uniformDistribution;
  std::normal_distribution<double> normalDistribution;

  double uniform()
  {
    ++count;
    return uniformDistribution(mt);
  }

  double Gaussian()
  {
    ++count;
    return normalDistribution(mt);
  }

  size_t integer(size_t i, size_t j) { return i + static_cast<size_t>(static_cast<double>(j + 1 - i) * uniform()); }

  std::pair<size_t, size_t> randomPairAdjacentIntegers(size_t size)
  {
    if (size <= 1) return std::make_pair(0, 0);
    size_t first = static_cast<size_t>(static_cast<double>(size - 1) * uniform());
    size_t second = first + 1;
    if (uniform() < 0.5) std::swap(first, second);
    return std::make_pair(first, second);
  }

  inline double3 UnitSphere()
  {
    double ran1, ran2, ranh, ransq;

    do
    {
      ran1 = 2.0 * uniform();
      ran2 = 2.0 * uniform();
      ransq = ran1 * ran1 + ran2 * ran2;
    } while (ransq >= 1.0);

    ranh = 2.0 * std::sqrt(1.0 - ransq);

    return double3(ran1 * ranh, ran2 * ranh, 1.0 - 2.0 * ransq);
  }

  double3 randomVectorOnUnitSphere();
  double3x3 randomRotationMatrix();
  double3x3 randomRotationAroundX(double angle);
  double3x3 randomRotationAroundY(double angle);
  double3x3 randomRotationAroundZ(double angle);
  simd_quatd randomSimdQuatd();
  simd_quatd smallRandomQuaternion(double angleRange);

  size_t categoricalDistribution(const std::vector<double> &probabilities)
  {
    std::discrete_distribution<size_t> d(probabilities.begin(), probabilities.end());
    return d(mt);
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const RandomNumber &r);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, RandomNumber &r);
};
