module;

export module property_pressure;

import std;

import archive;
import double3x3;
import averages;
import units;
import json;
export import property_block_average;

export struct PressureData
{
  PressureData():
    totalPressureTensor(double3x3()),
    excessPressureTensor(double3x3()),
    idealGasPressureTensor(double3x3()),
    totalPressure(0.0),
    excessPressure(0.0),
    idealGasPressure(0.0)
  {
  };

  PressureData(double3x3 totalPressureTensor, double3x3 excessPressureTensor, double3x3 idealGasPressureTensor,
               double totalPressure, double excessPressure, double idealGasPressure):
               totalPressureTensor(totalPressureTensor),
               excessPressureTensor(excessPressureTensor),
               idealGasPressureTensor(idealGasPressureTensor),
               totalPressure(totalPressure),
               excessPressure(excessPressure),
               idealGasPressure(idealGasPressure)
  {
  }

  inline PressureData& operator+=(const PressureData& b)
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

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const PressureData& l);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, PressureData& l);
};

export inline PressureData operator+(const PressureData &a, const PressureData &b)
{
  PressureData m{}; 

  m.totalPressureTensor = a.totalPressureTensor + b.totalPressureTensor;
  m.excessPressureTensor = a.excessPressureTensor + b.excessPressureTensor;
  m.idealGasPressureTensor = a.idealGasPressureTensor + b.idealGasPressureTensor;

  m.totalPressure = a.totalPressure + b.totalPressure;
  m.excessPressure = a.excessPressure + b.excessPressure;
  m.idealGasPressure = a.idealGasPressure + b.idealGasPressure;

  return m;
}

export inline PressureData operator-(const PressureData &a, const PressureData &b)
{
  PressureData m{};

  m.totalPressureTensor = a.totalPressureTensor - b.totalPressureTensor;
  m.excessPressureTensor = a.excessPressureTensor - b.excessPressureTensor;
  m.idealGasPressureTensor = a.idealGasPressureTensor - b.idealGasPressureTensor;

  m.totalPressure = a.totalPressure - b.totalPressure;
  m.excessPressure = a.excessPressure - b.excessPressure;
  m.idealGasPressure = a.idealGasPressure - b.idealGasPressure;

  return m;
}

export inline PressureData operator*(const PressureData &a, const PressureData &b)
{
  PressureData m{};

  m.totalPressureTensor = double3x3::multiplyElementWise(a.totalPressureTensor, b.totalPressureTensor);
  m.excessPressureTensor = double3x3::multiplyElementWise(a.excessPressureTensor, b.excessPressureTensor);
  m.idealGasPressureTensor = double3x3::multiplyElementWise(a.idealGasPressureTensor, b.idealGasPressureTensor);

  m.totalPressure = a.totalPressure * b.totalPressure;
  m.excessPressure = a.excessPressure * b.excessPressure;
  m.idealGasPressure = a.idealGasPressure * b.idealGasPressure;

  return m;
}

export inline PressureData operator*(const double& a, const PressureData &b)
{
  PressureData m{};

  m.totalPressureTensor = a * b.totalPressureTensor;
  m.excessPressureTensor = a * b.excessPressureTensor;
  m.idealGasPressureTensor = a * b.idealGasPressureTensor;

  m.totalPressure = a * b.totalPressure;
  m.excessPressure = a * b.excessPressure;
  m.idealGasPressure = a * b.idealGasPressure;

  return m;
}

export inline PressureData operator/(const PressureData &a, const double& b)
{
  PressureData m{};

  double inv_b = 1.0 / b;

  m.totalPressureTensor = inv_b * a.totalPressureTensor;
  m.excessPressureTensor = inv_b * a.excessPressureTensor;
  m.idealGasPressureTensor = inv_b * a.idealGasPressureTensor;

  m.totalPressure = inv_b * a.totalPressure;
  m.excessPressure = inv_b * a.excessPressure;
  m.idealGasPressure = inv_b * a.idealGasPressure;

  return m;
}

export inline PressureData sqrt(const PressureData &a)
{
  PressureData m{};

  m.totalPressureTensor = sqrt(a.totalPressureTensor);
  m.excessPressureTensor = sqrt(a.excessPressureTensor);
  m.idealGasPressureTensor = sqrt(a.idealGasPressureTensor);

  m.totalPressure = std::sqrt(a.totalPressure);
  m.excessPressure = std::sqrt(a.excessPressure);
  m.idealGasPressure = std::sqrt(a.idealGasPressure);

  return m;
}

export struct PropertyPressure : BlockAverage<PressureData>
{
  PropertyPressure() = default;

  PropertyPressure(std::size_t numberOfBlocks) : BlockAverage<PressureData>(numberOfBlocks) {}

  bool operator==(PropertyPressure const &) const = default;

  using BlockAverage<PressureData>::addSample;

  /// Builds the full PressureData sample (tensors and scalar traces) from the sampled
  /// ideal-gas pressure and excess pressure tensor.
  inline void addSample(std::size_t blockIndex, double idealGasPressureValue, double3x3 excessPressureTensor,
                        double weight)
  {
    double3x3 ideal_gas_pressure_tensor = double3x3(idealGasPressureValue, idealGasPressureValue, idealGasPressureValue);
    double excess_pressure = excessPressureTensor.trace() / 3.0;

    PressureData sample;
    sample.totalPressureTensor = excessPressureTensor + ideal_gas_pressure_tensor;
    sample.excessPressureTensor = excessPressureTensor;
    sample.idealGasPressureTensor = ideal_gas_pressure_tensor;
    sample.totalPressure = idealGasPressureValue + excess_pressure;
    sample.excessPressure = excess_pressure;
    sample.idealGasPressure = idealGasPressureValue;

    addSample(blockIndex, sample, weight);
  }

  std::string writeAveragesStatistics() const;
  nlohmann::json jsonAveragesStatistics() const;
};
