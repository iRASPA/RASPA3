module;

export module loading_data;

import std;

import archive;
import int3;
import simulationbox;
import component;

export struct LoadingData
{
  LoadingData() {};

  LoadingData(std::size_t size)
      : size(size),
        numberOfMolecules(std::vector<double>(size)),
        numberDensities(std::vector<double>(size)),
        inverseNumberDensities(std::vector<double>(size))
  {
  }

  bool operator==(LoadingData const&) const = default;

  void resize(std::size_t numberOfComponents)
  {
    numberOfMolecules.resize(numberOfComponents);
    numberDensities.resize(numberOfComponents);
    inverseNumberDensities.resize(numberOfComponents);
  }

  LoadingData(std::size_t size, const std::vector<std::size_t>& numberOfIntegerMolecules, const SimulationBox& box)
      : size(size), numberOfMolecules(size), numberDensities(size), inverseNumberDensities(size)
  {
    double inverseVolume = 1.0 / box.volume;

    totalNumberOfMolecules = 0.0;
    totalDensity = 0.0;
    for (std::size_t i = 0; i < this->numberOfMolecules.size(); ++i)
    {
      this->numberOfMolecules[i] = static_cast<double>(numberOfIntegerMolecules[i]);
      totalNumberOfMolecules += static_cast<double>(numberOfIntegerMolecules[i]);
      ;
      this->numberDensities[i] = static_cast<double>(numberOfIntegerMolecules[i]) * inverseVolume;
      totalDensity += static_cast<double>(numberOfIntegerMolecules[i]) * inverseVolume;
      this->inverseNumberDensities[i] = box.volume / static_cast<double>(numberOfIntegerMolecules[i]);
    }
  }

  double& operator()(std::size_t compA) { return numberOfMolecules[compA]; }

  std::string printStatus(std::size_t componentId, const std::string &componentName, double componentTotalMass,
                          double componentAmountOfExcessMolecules, std::optional<double> frameworkMass,
                          std::optional<int3> numberOfUnitCells) const;
  std::string printStatus(std::size_t componentId, const std::string &componentName, double componentTotalMass,
                          double componentAmountOfExcessMolecules, const LoadingData& average, const LoadingData& error,
                          std::optional<double> frameworkMass, std::optional<int3> numberOfUnitCells) const;

  inline LoadingData& operator+=(const LoadingData& b)
  {
    totalNumberOfMolecules += b.totalNumberOfMolecules;
    totalDensity += b.totalDensity;
    for (std::size_t i = 0; i < this->numberOfMolecules.size(); ++i)
    {
      numberOfMolecules[i] += b.numberOfMolecules[i];
      numberDensities[i] += b.numberDensities[i];
      inverseNumberDensities[i] += b.inverseNumberDensities[i];
    }
    return *this;
  }

  std::uint64_t versionNumber{1};

  std::size_t size;
  double totalNumberOfMolecules;
  double totalDensity;
  std::vector<double> numberOfMolecules;
  std::vector<double> numberDensities;
  std::vector<double> inverseNumberDensities;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const LoadingData& l);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, LoadingData& l);
};

export inline LoadingData operator+(const LoadingData& a, const LoadingData& b)
{
  LoadingData m(a.numberOfMolecules.size());

  m.totalNumberOfMolecules = a.totalNumberOfMolecules + b.totalNumberOfMolecules;
  m.totalDensity = a.totalDensity + b.totalDensity;
  for (std::size_t i = 0; i != a.numberOfMolecules.size(); ++i)
  {
    m.numberOfMolecules[i] = a.numberOfMolecules[i] + b.numberOfMolecules[i];
    m.numberDensities[i] = a.numberDensities[i] + b.numberDensities[i];
    m.inverseNumberDensities[i] = a.inverseNumberDensities[i] + b.inverseNumberDensities[i];
  }
  return m;
}

export inline LoadingData operator-(const LoadingData& a, const LoadingData& b)
{
  LoadingData m(a.numberOfMolecules.size());

  m.totalNumberOfMolecules = a.totalNumberOfMolecules - b.totalNumberOfMolecules;
  m.totalDensity = a.totalDensity - b.totalDensity;
  for (std::size_t i = 0; i != a.numberOfMolecules.size(); ++i)
  {
    m.numberOfMolecules[i] = a.numberOfMolecules[i] - b.numberOfMolecules[i];
    m.numberDensities[i] = a.numberDensities[i] - b.numberDensities[i];
    m.inverseNumberDensities[i] = a.inverseNumberDensities[i] - b.inverseNumberDensities[i];
  }
  return m;
}

export inline LoadingData operator*(const LoadingData& a, const LoadingData& b)
{
  LoadingData m(a.numberOfMolecules.size());

  m.totalNumberOfMolecules = a.totalNumberOfMolecules * b.totalNumberOfMolecules;
  m.totalDensity = a.totalDensity * b.totalDensity;
  for (std::size_t i = 0; i != a.numberOfMolecules.size(); ++i)
  {
    m.numberOfMolecules[i] = a.numberOfMolecules[i] * b.numberOfMolecules[i];
    m.numberDensities[i] = a.numberDensities[i] * b.numberDensities[i];
    m.inverseNumberDensities[i] = a.inverseNumberDensities[i] * b.inverseNumberDensities[i];
  }
  return m;
}

export inline LoadingData operator*(const double& a, const LoadingData& b)
{
  LoadingData m(b.numberOfMolecules.size());

  m.totalNumberOfMolecules = a * b.totalNumberOfMolecules;
  m.totalDensity = a * b.totalDensity;
  for (std::size_t i = 0; i != b.numberOfMolecules.size(); ++i)
  {
    m.numberOfMolecules[i] = a * b.numberOfMolecules[i];
    m.numberDensities[i] = a * b.numberDensities[i];
    m.inverseNumberDensities[i] = a * b.inverseNumberDensities[i];
  }
  return m;
}

export inline LoadingData operator/(const LoadingData& a, const double& b)
{
  LoadingData m(a.numberOfMolecules.size());

  double temp = 1.0 / b;
  m.totalNumberOfMolecules = a.totalNumberOfMolecules * temp;
  m.totalDensity = a.totalDensity * temp;
  for (std::size_t i = 0; i != a.numberOfMolecules.size(); ++i)
  {
    m.numberOfMolecules[i] = a.numberOfMolecules[i] * temp;
    m.numberDensities[i] = a.numberDensities[i] * temp;
    m.inverseNumberDensities[i] = a.inverseNumberDensities[i] * temp;
  }
  return m;
}

export inline LoadingData sqrt(const LoadingData& a)
{
  LoadingData m(a.numberOfMolecules.size());

  m.totalNumberOfMolecules = std::sqrt(a.totalNumberOfMolecules);
  m.totalDensity = std::sqrt(a.totalDensity);
  for (std::size_t i = 0; i != a.numberOfMolecules.size(); ++i)
  {
    m.numberOfMolecules[i] = std::sqrt(a.numberOfMolecules[i]);
    m.numberDensities[i] = std::sqrt(a.numberDensities[i]);
    m.inverseNumberDensities[i] = std::sqrt(a.inverseNumberDensities[i]);
  }
  return m;
}
