module;

export module pressures;

import std;

import archive;
import double3x3;

export struct Pressures
{
  Pressures():
    totalPressureTensor(double3x3()),
    excessPressureTensor(double3x3()),
    idealGasPressureTensor(double3x3()),
    totalPressure(0.0),
    excessPressure(0.0),
    idealGasPressure(0.0)
  {
  };

  Pressures(double3x3 totalPressureTensor, double3x3 excessPressureTensor, double3x3 idealGasPressureTensor,
            double totalPressure, double excessPressure, double idealGasPressure):
            totalPressureTensor(totalPressureTensor),
            excessPressureTensor(excessPressureTensor),
            idealGasPressureTensor(idealGasPressureTensor),
            totalPressure(totalPressure),
            excessPressure(excessPressure),
            idealGasPressure(idealGasPressure)
  {
  }

  bool operator==(Pressures const&) const = default;

  inline Pressures& operator+=(const Pressures& b)
  {
    totalPressureTensor += b.totalPressureTensor;
    excessPressureTensor += b.excessPressureTensor;
    idealGasPressureTensor += b.idealGasPressureTensor;

    totalPressure += b.totalPressure;
    excessPressure += b.excessPressure;
    idealGasPressure += b.idealGasPressure;

    return *this;
  }

  std::uint64_t versionNumber{1};

  double3x3 totalPressureTensor{};
  double3x3 excessPressureTensor{};
  double3x3 idealGasPressureTensor{};
  double totalPressure{};
  double excessPressure{};
  double idealGasPressure{};

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Pressures& l);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Pressures& l);
};

export inline Pressures operator+(const Pressures& a, const Pressures& b)
{
  Pressures m{}; 

  m.totalPressureTensor = a.totalPressureTensor + b.totalPressureTensor;
  m.excessPressureTensor = a.excessPressureTensor + b.excessPressureTensor;
  m.idealGasPressureTensor = a.idealGasPressureTensor + b.idealGasPressureTensor;

  m.totalPressure = a.totalPressure + b.totalPressure;
  m.excessPressure = a.excessPressure + b.excessPressure;
  m.idealGasPressure = a.idealGasPressure + b.idealGasPressure;

  return m;
}

export inline Pressures operator-(const Pressures& a, const Pressures& b)
{
  Pressures m{};

  m.totalPressureTensor = a.totalPressureTensor - b.totalPressureTensor;
  m.excessPressureTensor = a.excessPressureTensor - b.excessPressureTensor;
  m.idealGasPressureTensor = a.idealGasPressureTensor - b.idealGasPressureTensor;

  m.totalPressure = a.totalPressure - b.totalPressure;
  m.excessPressure = a.excessPressure - b.excessPressure;
  m.idealGasPressure = a.idealGasPressure - b.idealGasPressure;

  return m;
}

export inline Pressures operator*(const Pressures& a, const Pressures& b)
{
  Pressures m{};

  m.totalPressureTensor = double3x3::multiplyElementWise(a.totalPressureTensor, b.totalPressureTensor);
  m.excessPressureTensor = double3x3::multiplyElementWise(a.excessPressureTensor, b.excessPressureTensor);
  m.idealGasPressureTensor = double3x3::multiplyElementWise(a.idealGasPressureTensor, b.idealGasPressureTensor);

  m.totalPressure = a.totalPressure * b.totalPressure;
  m.excessPressure = a.excessPressure * b.excessPressure;
  m.idealGasPressure = a.idealGasPressure * b.idealGasPressure;

  return m;
}

export inline Pressures operator*(const double& a, const Pressures& b)
{
  Pressures m{};

  m.totalPressureTensor = a * b.totalPressureTensor;
  m.excessPressureTensor = a * b.excessPressureTensor;
  m.idealGasPressureTensor = a * b.idealGasPressureTensor;

  m.totalPressure = a * b.totalPressure;
  m.excessPressure = a * b.excessPressure;
  m.idealGasPressure = a * b.idealGasPressure;

  return m;
}

export inline Pressures operator/(const Pressures& a, const double& b)
{
  Pressures m{};

  double inv_b = 1.0 / b;

  m.totalPressureTensor = inv_b * a.totalPressureTensor;
  m.excessPressureTensor = inv_b * a.excessPressureTensor;
  m.idealGasPressureTensor = inv_b * a.idealGasPressureTensor;

  m.totalPressure = inv_b * a.totalPressure;
  m.excessPressure = inv_b * a.excessPressure;
  m.idealGasPressure = inv_b * a.idealGasPressure;

  return m;
}

export inline Pressures sqrt(const Pressures& a)
{
  Pressures m{};

  m.totalPressureTensor = sqrt(a.totalPressureTensor);
  m.excessPressureTensor = sqrt(a.excessPressureTensor);
  m.idealGasPressureTensor = sqrt(a.idealGasPressureTensor);

  m.totalPressure = std::sqrt(a.totalPressure);
  m.excessPressure = std::sqrt(a.excessPressure);
  m.idealGasPressure = std::sqrt(a.idealGasPressure);

  return m;
}
