module;

module apollonius_pore_size_distribution;

import std;

import double3;
import crystal;
import pair_interactions;
import apollonius_accessibility;
import exact_pore_size_distribution;
import voronoi_pore_size_distribution;

void ApolloniusPoreSizeDistribution::run(const PairInteractions& interactions, const Crystal& framework,
                                         std::string probePseudoAtom, std::optional<double> maximumDiameter,
                                         std::optional<std::size_t> numberOfBins, std::size_t subdivisions)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("ApolloniusPoreSizeDistribution: Unknown probe-atom type\n");
  }
  const double probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  fractionalPositions.reserve(framework.atoms.size());
  radii.reserve(framework.atoms.size());
  for (const CrystalAtom& atom : framework.atoms)
  {
    fractionalPositions.push_back(framework.unitCell.inverseCell * atom.position);
    std::size_t type = atom.type;
    radii.push_back(0.5 * interactions(type, type).sizeParameter);
  }

  auto build = [&](double probeRadius)
  {
    return ApolloniusAccessibility::create(framework.unitCell, fractionalPositions, radii, probeRadius)
        .accessibility;
  };

  curve = exactPoreSizeDistribution(build, framework.unitCell.volume, maximumDiameter.value_or(20.0),
                                    numberOfBins.value_or(100), subdivisions, probeRadius);

  writePoreSizeDistribution(framework, "apollonius", probePseudoAtom, curve);
}
