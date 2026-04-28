module;

export module velocity_autocorrelation_function_data;

import std;

import archive;

export struct VelocityAutoCorrelationFunctionData
{
  VelocityAutoCorrelationFunctionData():
    time(0.0),
    xyz(0.0),
    x(0.0),
    y(0.0),
    z(0.0),
    numberOfSamples(0.0)
  {
  };

  VelocityAutoCorrelationFunctionData(double time, double xyz, double x, double y, double z, double numberOfSamples):
    time(time),
    xyz(xyz),
    x(x),
    y(y),
    z(z),
    numberOfSamples(numberOfSamples)
  {
  };

  std::uint64_t versionNumber{1};

  double time{};
  double xyz{};
  double x{};
  double y{};
  double z{};
  double numberOfSamples{};

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const VelocityAutoCorrelationFunctionData& l);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, VelocityAutoCorrelationFunctionData& l);
};

