module;

module apollonius_pore_size_distribution;

import std;

import double3;
import atom;
import framework;
import forcefield;
import apollonius_accessibility;
import exact_pore_size_distribution;
import voronoi_pore_size_distribution;

void ApolloniusPoreSizeDistribution::run(const ForceField& forceField, const Framework& framework,
                                         std::string probePseudoAtom, std::optional<double> maximumDiameter,
                                         std::optional<std::size_t> numberOfBins, std::size_t subdivisions)
{
  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("ApolloniusPoreSizeDistribution: Unknown probe-atom type\n");
  }
  const double probeRadius = 0.5 * forceField[probeType.value()].sizeParameter();

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  fractionalPositions.reserve(framework.unitCellAtoms.size());
  radii.reserve(framework.unitCellAtoms.size());
  for (const Atom& atom : framework.unitCellAtoms)
  {
    fractionalPositions.push_back(framework.simulationBox.inverseCell * atom.position);
    std::size_t type = static_cast<std::size_t>(atom.type);
    radii.push_back(0.5 * forceField(type, type).sizeParameter());
  }

  auto build = [&](double probeRadius)
  {
    return ApolloniusAccessibility::create(framework.simulationBox, fractionalPositions, radii, probeRadius)
        .accessibility;
  };

  curve = exactPoreSizeDistribution(build, framework.simulationBox.volume, maximumDiameter.value_or(20.0),
                                    numberOfBins.value_or(100), subdivisions, probeRadius);

  writePoreSizeDistribution(framework, "apollonius", probePseudoAtom, curve);
}
