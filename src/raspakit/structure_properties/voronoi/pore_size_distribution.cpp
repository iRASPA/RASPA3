module;

module voronoi_pore_size_distribution;

import std;

import double3;
import units;
import atom;
import framework;
import forcefield;
import skspacegroupdatabase;
import voronoi_accessibility;
import exact_pore_size_distribution;

// The bare radii of the framework's atoms, and where they are in fractional coordinates.
std::pair<std::vector<double3>, std::vector<double>> frameworkSpheres(const ForceField& forceField,
                                                                     const Framework& framework)
{
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
  return {fractionalPositions, radii};
}


void writePoreSizeDistribution(const Framework& framework, const std::string& diagramName,
                               const PoreSizeDistributionCurve& curve)
{
  std::ofstream report;
  report.open(std::format("{}.{}.psd.txt", framework.name, diagramName));

  std::print(report, "# Pore-size distribution ({}, exact)\n", diagramName);
  std::print(report, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(report, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(report, "# Number of framework atoms: {}\n", framework.unitCellAtoms.size());
  std::print(report, "# Framework volume: {} [Å³]\n", curve.cellVolume);
  std::print(report, "# Framework mass: {} [g/mol]\n", framework.unitCellMass);
  std::print(report, "# Void volume: {} [Å³], void fraction {}\n", curve.voidVolume,
             (curve.cellVolume > 0.0) ? curve.voidVolume / curve.cellVolume : 0.0);
  std::print(report, "# Diameters evaluated: {}, up to {} [Å]\n", curve.points.size(),
             curve.points.empty() ? 0.0 : curve.points.back().diameter);
  std::print(report, "# Seconds: {}, over {} evaluations of the surface\n", curve.seconds,
             curve.numberOfEvaluations);
  std::print(report, "#\n");
  std::print(report, "# The pore size at a point of the void is the diameter of the largest sphere that holds\n");
  std::print(report, "# the point and fits in the void, and the cumulative volume is the volume where that\n");
  std::print(report, "# diameter is at least d. That volume is the cell less the volume enclosed by the\n");
  std::print(report, "# solvent excluded surface of the framework at probe radius d/2, and the distribution\n");
  std::print(report, "# is its derivative, which is an integral over the reentrant part of that surface\n");
  std::print(report, "# alone. Both are closed forms over the same patches, arcs and vertices: no point of\n");
  std::print(report, "# the void is tested and no trial sphere is drawn, so no row of this file is an\n");
  std::print(report, "# estimate of another row's value.\n");
  std::print(report, "#\n");
  std::print(report, "# The distribution is not a function of d alone. A pore holds the whole of its volume at\n");
  std::print(report, "# the diameter of the largest sphere that fits in it -- a cavity that is a ball of radius\n");
  std::print(report, "# a is filled by the positions of any probe up to that radius, so every point of it has\n");
  std::print(report, "# the one pore size 2a -- and only the corrugation of the walls puts volume at smaller\n");
  std::print(report, "# sizes. So there is a spike at the largest sphere of every family of pores the void\n");
  std::print(report, "# falls into, and the rows below are the continuous part between them.\n");
  std::print(report, "#\n");
  std::print(report, "# The spikes are found by bisection, wherever the volume lost across a step exceeds what\n");
  std::print(report, "# the derivative accounts for. A corner of the continuous part loses that excess under\n");
  std::print(report, "# refinement and a cliff keeps it, so what is listed here survived the narrowing of its\n");
  std::print(report, "# own interval by a factor of a thousand or more:\n");
  std::print(report, "#\n");
  std::print(report, "#      d [Å]        weight     cornered within [Å]\n");
  for (const PoreSizeSpike& spike : curve.spikes)
  {
    std::print(report, "#   {:9.5f}    {:10.6f}    {:.2e}\n", spike.diameter, spike.weight, spike.bracket);
  }
  std::print(report, "#\n");
  std::print(report, "# The continuous part integrates to {:.6f} and the spikes weigh {:.6f}, which come to\n",
             curve.integral, curve.singularWeight);
  std::print(report, "# {:.6f} against one. The void has some pore size everywhere, so the two have to make up\n",
             curve.integral + curve.singularWeight);
  std::print(report, "# the whole of it, and they are computed independently: how near that comes is the check\n");
  std::print(report, "# of the volume against its own derivative over the whole range.\n");
  std::print(report, "#\n");
  std::print(report, "# The largest sphere that fits in the void is {:.5f} [Å] across, which is where the last\n",
             curve.largestDiameter);
  std::print(report, "# spike stands. Beyond it {:.2e} of the void is left, which ought to be nothing: where it\n",
             curve.truncatedWeight);
  std::print(report, "# is not, the range asked for stopped short of the largest pore.\n");
  std::print(report, "#\n");
  std::print(report, "# column 1: diameter d [Å]\n");
  std::print(report, "# column 2: P(d) [1/Å], normalised to integrate to one\n");
  std::print(report, "# column 3: cumulative pore volume, as a fraction of the void volume\n");
  std::print(report, "# column 4: cumulative pore volume [Å³]\n");
  std::print(report, "# column 5: cumulative pore volume [cm³/g]\n");
  std::print(report, "# columns 6-8: the share of column 2 on surfaces the probe can reach, on surfaces\n");
  std::print(report, "#              sealed off from it, and on surfaces the network could not place\n");
  std::print(report, "# columns 9-13: arcs, of them cusped, vertices, of them clipped, of them degenerate\n");

  // Å³ per unit cell to cm³ per gram.
  const double gramsPerCell = framework.unitCellMass / Units::AvogadroConstant;
  const double toVolumePerMass = (gramsPerCell > 0.0) ? 1.0e-24 / gramsPerCell : 0.0;

  for (const PoreSizeDistributionPoint& point : curve.points)
  {
    // The diameter to more figures than a plot needs, because the columns are meant to be differenced
    // against one another and a spacing read back from a rounded abscissa is not the spacing used.
    std::print(report, "{:11.6f} {:14.8f} {:12.8f} {:14.5f} {:12.6f} {:14.8f} {:14.8f} {:14.8f} {:8} {:8} {:8} {:8} {:8}\n",
               point.diameter, point.distribution, point.cumulative, point.poreVolume,
               point.poreVolume * toVolumePerMass, point.accessible, point.inaccessible, point.undecided,
               point.numberOfArcs, point.cuspedArcs, point.numberOfVertices, point.clippedVertices,
               point.degenerateVertices);
  }
  report.close();
}


void VoronoiPoreSizeDistribution::run(const ForceField& forceField, const Framework& framework,
                                      std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins,
                                      std::size_t subdivisions)
{
  auto [fractionalPositions, radii] = frameworkSpheres(forceField, framework);

  auto build = [&](double probeRadius)
  { return VoronoiAccessibility::create(framework.simulationBox, fractionalPositions, radii, probeRadius); };

  curve = exactPoreSizeDistribution(build, framework.simulationBox.volume, maximumDiameter.value_or(20.0),
                                    numberOfBins.value_or(100), subdivisions);

  writePoreSizeDistribution(framework, "voronoi", curve);
}
