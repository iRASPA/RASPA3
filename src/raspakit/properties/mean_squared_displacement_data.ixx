module;

export module mean_squared_displacement_data;

import std;

import archive;

export struct MeanSquaredDisplacementData
{
  MeanSquaredDisplacementData():
    time(0.0),
    xyz(0.0),
    x(0.0),
    y(0.0),
    z(0.0),
    numberOfSamples(0.0)
  {
  };

  MeanSquaredDisplacementData(double time, double xyz, double x, double y, double z, double numberOfSamples):
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

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MeanSquaredDisplacementData& l);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MeanSquaredDisplacementData& l);
};

