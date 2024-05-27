module;

#ifdef USE_LEGACY_HEADERS
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <map>
#include <optional>
#include <ostream>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#endif

export module framework;

#ifndef USE_LEGACY_HEADERS
import <format>;
import <tuple>;
import <vector>;
import <string>;
import <chrono>;
import <cstdint>;
import <fstream>;
import <sstream>;
import <ostream>;
import <vector>;
import <array>;
import <map>;
import <optional>;
import <span>;
#if defined(__has_include) && __has_include(<print>)
import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
import print;
#endif

import stringutils;
import archive;
import randomnumbers;
import int3;
import double3;
import averages;
import atom;
import forcefield;
import simulationbox;
import property_widom;
import isotherm;
import multi_site_isotherm;
import bond_potential;
import hdf5;

export struct Framework
{
  Framework();
  Framework(size_t currentComponent, const ForceField &forceField, const std::string &componentName,
            std::optional<const std::string> fileName, int3 numberOfUnitCells) noexcept(false);
  Framework(size_t componentId, const ForceField &forceField, std::string componentName, SimulationBox simulationBox,
            size_t spaceGroupHallNumber, std::vector<Atom> definedAtoms, int3 numberOfUnitCells) noexcept(false);

  uint64_t versionNumber{1};

  SimulationBox simulationBox;
  size_t spaceGroupHallNumber{1};
  int3 numberOfUnitCells{1, 1, 1};

  size_t frameworkId{0};
  std::string name{};
  std::optional<std::string> filenameData{};
  std::string filename{};

  bool rigid{true};

  double mass{0.0};
  double unitCellMass{0.0};
  double netCharge{0.0};
  std::vector<Atom> definedAtoms{};
  std::vector<Atom> atoms{};
  std::vector<Atom> unitCellAtoms;

  std::vector<size_t> chiralCenters{};
  std::vector<BondPotential> bonds{};
  std::vector<std::pair<size_t, size_t>> bondDipoles{};
  std::vector<std::tuple<size_t, size_t, size_t>> bends{};
  std::vector<std::pair<size_t, size_t>> UreyBradley{};
  std::vector<std::tuple<size_t, size_t, size_t, size_t>> inversionBends{};
  std::vector<std::tuple<size_t, size_t, size_t, size_t>> Torsion{};
  std::vector<std::tuple<size_t, size_t, size_t, size_t>> ImproperTorsions{};
  std::vector<std::tuple<size_t, size_t, size_t>> bondBonds{};
  std::vector<std::tuple<size_t, size_t, size_t>> stretchBends{};
  std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendBends{};
  std::vector<std::tuple<size_t, size_t, size_t, size_t>> stretchTorsions{};
  std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendTorsions{};
  std::vector<std::pair<size_t, size_t>> intraVDW{};
  std::vector<std::pair<size_t, size_t>> intraCoulomb{};
  std::vector<std::pair<size_t, size_t>> excludedIntraCoulomb{};

  void readFramework(const ForceField &forceField, const std::string &fileName);

  void expandDefinedAtomsToUnitCell();
  void makeSuperCell();

  std::string printStatus(const ForceField &forceField) const;
  std::string printBreakthroughStatus() const;
  void logStatus(HDF5Writer &hdf5, const ForceField &forceField) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Framework &c);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Framework &c);

  std::string repr() const;
};

template <typename T>
std::vector<T> parseListOfParameters(const std::string &arguments, size_t lineNumber)
{
  std::vector<T> list{};

  std::string str;
  std::istringstream ss(arguments);

  while (ss >> str)
  {
    if (trim(str).rfind("//", 0) == 0)
    {
      if (list.empty())
      {
        throw std::runtime_error(std::format("No values could be read at line: {}\n", lineNumber));
      }
      return list;
    }
    T value;
    std::istringstream s(str);
    if (s >> value)
    {
      list.push_back(value);
    }
    else
    {
      if (list.empty())
      {
        throw std::runtime_error(std::format("No values could be read at line: {}\n", lineNumber));
      }
      return list;
    }
  };

  return list;
}
