module;

module minimization_cell_layout;

import std;

import double3x3;

namespace
{
bool equalCaseInsensitive(std::string_view lhs, std::string_view rhs)
{
  return lhs.size() == rhs.size() &&
         std::ranges::equal(
             lhs, rhs, [](char a, char b)
             { return std::tolower(static_cast<unsigned char>(a)) == std::tolower(static_cast<unsigned char>(b)); });
}

double3x3 diagonalBasis(std::size_t axis)
{
  double3x3 basis{};
  basis.mm[axis][axis] = 1.0;
  return basis;
}

double3x3 shearBasis(std::size_t a, std::size_t b)
{
  double3x3 basis{};
  basis.mm[a][b] = 0.5;
  basis.mm[b][a] = 0.5;
  return basis;
}
}  // namespace

bool CellMinimizationLayout::isCompatibilityAlias() const noexcept
{
  return requestedType == CellMinimizationType::RegularUpperTriangle ||
         requestedType == CellMinimizationType::MonoclinicUpperTriangle;
}

std::string CellMinimizationLayout::name() const { return cellMinimizationTypeName(requestedType); }

std::optional<CellMinimizationType> cellMinimizationTypeFromString(std::string_view value)
{
  if (equalCaseInsensitive(value, "REGULAR_UPPER_TRIANGLE")) return CellMinimizationType::RegularUpperTriangle;
  if (equalCaseInsensitive(value, "MONOCLINIC_UPPER_TRIANGLE")) return CellMinimizationType::MonoclinicUpperTriangle;
  for (CellMinimizationType type :
       {CellMinimizationType::Fixed, CellMinimizationType::Isotropic, CellMinimizationType::Anisotropic,
        CellMinimizationType::Monoclinic, CellMinimizationType::Regular, CellMinimizationType::RegularUpperTriangle,
        CellMinimizationType::MonoclinicUpperTriangle})
  {
    if (equalCaseInsensitive(value, cellMinimizationTypeName(type))) return type;
  }
  return std::nullopt;
}

std::optional<MonoclinicAngleType> monoclinicAngleTypeFromString(std::string_view value)
{
  for (MonoclinicAngleType type : {MonoclinicAngleType::Alpha, MonoclinicAngleType::Beta, MonoclinicAngleType::Gamma})
  {
    if (equalCaseInsensitive(value, monoclinicAngleTypeName(type))) return type;
  }
  return std::nullopt;
}

std::string cellMinimizationTypeName(CellMinimizationType type)
{
  switch (type)
  {
    case CellMinimizationType::Fixed:
      return "Fixed";
    case CellMinimizationType::Isotropic:
      return "Isotropic";
    case CellMinimizationType::Anisotropic:
      return "Anisotropic";
    case CellMinimizationType::Monoclinic:
      return "Monoclinic";
    case CellMinimizationType::Regular:
      return "Regular";
    case CellMinimizationType::RegularUpperTriangle:
      return "RegularUpperTriangle";
    case CellMinimizationType::MonoclinicUpperTriangle:
      return "MonoclinicUpperTriangle";
  }
  std::unreachable();
}

std::string monoclinicAngleTypeName(MonoclinicAngleType type)
{
  switch (type)
  {
    case MonoclinicAngleType::Alpha:
      return "Alpha";
    case MonoclinicAngleType::Beta:
      return "Beta";
    case MonoclinicAngleType::Gamma:
      return "Gamma";
  }
  std::unreachable();
}

CellMinimizationLayout makeCellMinimizationLayout(CellMinimizationType type, MonoclinicAngleType angle)
{
  CellMinimizationLayout layout{.requestedType = type, .monoclinicAngle = angle};
  switch (type)
  {
    case CellMinimizationType::Fixed:
      break;
    case CellMinimizationType::Isotropic:
      layout.bases.emplace_back(1.0, 1.0, 1.0);
      break;
    case CellMinimizationType::Anisotropic:
      layout.bases = {diagonalBasis(0), diagonalBasis(1), diagonalBasis(2)};
      break;
    case CellMinimizationType::Monoclinic:
    case CellMinimizationType::MonoclinicUpperTriangle:
      layout.bases = {diagonalBasis(0), diagonalBasis(1), diagonalBasis(2)};
      switch (angle)
      {
        case MonoclinicAngleType::Alpha:
          layout.bases.push_back(shearBasis(1, 2));
          break;
        case MonoclinicAngleType::Beta:
          layout.bases.push_back(shearBasis(0, 2));
          break;
        case MonoclinicAngleType::Gamma:
          layout.bases.push_back(shearBasis(0, 1));
          break;
      }
      break;
    case CellMinimizationType::Regular:
    case CellMinimizationType::RegularUpperTriangle:
      layout.bases = {diagonalBasis(0), shearBasis(0, 1), shearBasis(0, 2),
                      diagonalBasis(1), shearBasis(1, 2), diagonalBasis(2)};
      break;
  }
  return layout;
}

double3x3 cellStrainSecondDerivative(const CellMinimizationLayout& layout, std::size_t a, std::size_t b)
{
  return 0.5 * (layout.bases.at(a) * layout.bases.at(b) + layout.bases.at(b) * layout.bases.at(a));
}
