module;

module blocking_pockets;

import std;

import stringutils;
import double3;
import double4;
import atom;
import forcefield;
import framework;
import structure_input;

import unit_cell;
import voronoi_blocking_spheres;

namespace BlockingPockets
{

Specification parse(const nlohmann::basic_json<nlohmann::raspa_map> &item)
{
  if (item.is_string())
  {
    const std::string value = item.get<std::string>();
    if (caseInSensStringCompare(value, "auto"))
    {
      return Specification{.pockets = {}, .automatic = true};
    }
    return Specification{.pockets = readBlockingPocketFile(value), .automatic = false};
  }

  Specification specification{};
  for (const auto &[_, sphere] : item.items())
  {
    if (!sphere.is_array() || sphere.size() != 4)
    {
      throw std::runtime_error(
          std::format("[Blocking pockets]: item {} must be an array with four elements, the fractional positions "
                      "s_x, s_y, s_z and a radius in Angstrom; 'BlockingPockets' as a whole must be a list of "
                      "those, or the string 'auto', or the name of a '.block' file\n",
                      sphere.dump()));
    }

    const std::vector<double> data = sphere.get<std::vector<double>>();
    specification.pockets.push_back(double4(data[0], data[1], data[2], data[3]));
  }

  return specification;
}

std::vector<double4> readBlockingPocketFile(const std::string &fileName)
{
  std::istringstream stream{readFileContent(fileName, ".block")};

  std::size_t numberOfSpheres{};
  if (!(stream >> numberOfSpheres))
  {
    throw std::runtime_error(std::format(
        "[Blocking pockets]: file '{}' does not start with the number of spheres in it\n", fileName));
  }

  std::vector<double4> pockets;
  pockets.reserve(numberOfSpheres);
  for (std::size_t i = 0; i != numberOfSpheres; ++i)
  {
    double3 center{};
    double radius{};
    if (!(stream >> center.x >> center.y >> center.z >> radius))
    {
      throw std::runtime_error(
          std::format("[Blocking pockets]: file '{}' says it holds {} spheres but only {} could be read; each "
                      "sphere is a line of the fractional positions s_x, s_y, s_z and a radius in Angstrom\n",
                      fileName, numberOfSpheres, i));
    }
    pockets.push_back(double4(center.x, center.y, center.z, radius));
  }

  return pockets;
}

std::vector<double4> compute(const Framework &framework, const ForceField &forceField)
{
  const UnitCell unitCell = StructureInput::makeUnitCell(framework.simulationBox);

  std::vector<double3> fractionalPositions = framework.fractionalAtomPositionsUnitCell();
  std::vector<double> radii;
  radii.reserve(framework.unitCellAtoms.size());
  for (const Atom &atom : framework.unitCellAtoms)
  {
    std::size_t type = static_cast<std::size_t>(atom.type);
    radii.push_back(0.5 * forceField(type, type).sizeParameter());
  }

  // A framework of points has no pocket in it, and reporting none would be read as a framework with none.
  // The force fields this happens with are not broken: giving a framework atom no size of its own and naming
  // every framework-guest pair outright is how several of the zeolite models are written, and they simulate
  // correctly. It is only the geometry that cannot be taken from them, and only sizes it is missing, so what
  // is wanted is the sentence that says which.
  if (std::ranges::all_of(radii, [](double radius) { return radius <= 0.0; }))
  {
    throw std::runtime_error(std::format(
        "[Blocking pockets]: the force field gives no van der Waals size to any atom of framework '{}', so it "
        "describes a framework of points, in which no cavity is closed to anything. Force fields that leave the "
        "framework atoms without self-interactions and name every framework-guest pair outright are of this kind. "
        "Give the framework atoms self-interactions, or list the blocking pockets instead of asking for 'auto'\n",
        framework.name));
  }

  // Nitrogen, at the size the force field gives it where it has one and at the standard size where it does
  // not, which is the case for every force field written for a simulation rather than for an analysis.
  std::optional<std::size_t> probeType = forceField.findPseudoAtom(std::string(nitrogenProbeName));
  double probeRadius = 0.5 * (probeType.has_value()
                                  ? forceField(probeType.value(), probeType.value()).sizeParameter()
                                  : nitrogenProbeSizeParameter);

  VoronoiBlockingSpheres blocking{};
  blocking.compute(unitCell, fractionalPositions, radii, probeRadius);

  std::vector<double4> pockets;
  pockets.reserve(blocking.spheres.size());
  for (const BlockingSphere &sphere : blocking.spheres)
  {
    pockets.push_back(
        double4(sphere.centerFractional.x, sphere.centerFractional.y, sphere.centerFractional.z, sphere.radius));
  }

  return pockets;
}

}  // namespace BlockingPockets
