export module framework;

import <format>;
import <tuple>;
import <vector>;
import <string>;
import <chrono>;
import <cstdint>;
import <fstream>;
import <ostream>;
import <vector>;
import <array>;
import <map>;
import <optional>;
import <span>;

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

export struct Framework
{
  Framework();
  Framework(size_t currentComponent, const ForceField& forceField, const std::string &componentName, 
            std::optional<const std::string> fileName) noexcept(false);
  Framework(size_t componentId, std::string componentName, double mass, SimulationBox simulationBox, 
            size_t spaceGroupHallNumber, std::vector<Atom> definedAtoms, int3 numberOfUnitCells) noexcept(false);

  uint64_t versionNumber{ 1 };

  SimulationBox simulationBox;
  size_t spaceGroupHallNumber{ 1 };
  int3 numberOfUnitCells{ 1, 1, 1};

  size_t frameworkId{ 0 };
  std::string name{};
  std::optional<std::string> filenameData{};
  std::string filename{};

  bool rigid { true };

  double mass{ 0.0 };
  double netCharge{ 0.0 };
  std::vector<Atom> definedAtoms{};
  std::vector<Atom> atoms{};
  std::vector<Atom> unitCellAtoms;

  std::vector<size_t> chiralCenters{};
  std::vector<BondPotential> bonds{};
  std::vector<std::pair<size_t, size_t>> bondDipoles{};
  std::vector<std::tuple<size_t, size_t, size_t>> bends{};
  std::vector<std::pair<size_t, size_t>>  UreyBradley{};
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

  std::vector<Atom> frameworkAtoms() const;
  void readFramework(const ForceField& forceField, const std::string& fileName);

  std::string printStatus(const ForceField& forceField) const;
  std::string printBreakthroughStatus() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Framework &c);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Framework &c);

  std::string repr() const;
};
