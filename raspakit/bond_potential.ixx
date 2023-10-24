export module bond_potential;

import <string>;
import <map>;
import <vector>;
import <array>;
import <print>;

import stringutils;

export const size_t maximumNumberOfBondParameters{ 4 };

export enum class BondType : size_t
{
    Undefined = 0,
    Fixed = 1,
    Rigid = 2,
    Harmonic = 3
};

export struct BondPotential
{
  BondType bondType;
  std::pair<size_t, size_t> bondIds;
  std::array<double, maximumNumberOfBondParameters> parameters;

  BondPotential(): bondType(BondType::Undefined), bondIds(0,0) {}
  BondPotential(const BondType bondType, std::pair<size_t, size_t> bondIds): bondType(bondType), bondIds(bondIds) {}

  bool operator==(BondPotential const&) const = default;

  std::string print() const
  {
      switch (bondType)
      {
      case BondType::Fixed:
          return std::format("FIXED_BOND ({}-{})\n", bondIds.first, bondIds.second);
      case BondType::Rigid:
          return std::format("RIGID_BOND ({}-{})\n", bondIds.first, bondIds.second);
      case BondType::Harmonic:
          return std::format("HARMONIC_BOND ({}-{}): p_0/k_B={} [K/A^2], p_1={} [A]\n", bondIds.first, bondIds.second, parameters[0], parameters[1]);
      default:
          return "Unknown potential";
      }
  }

  static inline std::vector<size_t> numberOfBondParameters{ 0, 0, 2 };


  struct comp {
      bool operator() (const std::string& lhs, const std::string& rhs) const {
#if defined(_WIN32)
          return _stricmp(lhs.c_str(), rhs.c_str()) < 0;
#else
          return strcasecmp(lhs.c_str(), rhs.c_str()) < 0;
#endif

      }
  };

  static inline std::map<std::string, BondType, comp> bondDefinitionForString{
    {"FIXED_BOND", BondType::Fixed} ,
    {"RIGID_BOND", BondType::Rigid} ,
    {"HARMONIC_BOND", BondType::Harmonic} };

};


