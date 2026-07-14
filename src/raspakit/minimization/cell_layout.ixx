module;

export module minimization_cell_layout;

import std;

import double3x3;

export enum class CellMinimizationType : std::uint8_t {
  Fixed,
  Isotropic,
  Anisotropic,
  Monoclinic,
  Regular,
  RegularUpperTriangle,
  MonoclinicUpperTriangle
};

export enum class MonoclinicAngleType : std::uint8_t { Alpha, Beta, Gamma };

export struct CellMinimizationLayout
{
  CellMinimizationType requestedType{CellMinimizationType::Fixed};
  MonoclinicAngleType monoclinicAngle{MonoclinicAngleType::Beta};
  std::vector<double3x3> bases{};

  std::size_t size() const noexcept { return bases.size(); }
  bool empty() const noexcept { return bases.empty(); }
  bool isCompatibilityAlias() const noexcept;
  std::string name() const;
};

export std::optional<CellMinimizationType> cellMinimizationTypeFromString(std::string_view value);
export std::optional<MonoclinicAngleType> monoclinicAngleTypeFromString(std::string_view value);
export std::string cellMinimizationTypeName(CellMinimizationType type);
export std::string monoclinicAngleTypeName(MonoclinicAngleType type);
export CellMinimizationLayout makeCellMinimizationLayout(CellMinimizationType type, MonoclinicAngleType angle);
export double3x3 cellStrainSecondDerivative(const CellMinimizationLayout& layout, std::size_t a, std::size_t b);
