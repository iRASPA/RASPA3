module;

module energy_shared_pore_size_distribution;

import std;

import uint3;
import double3;
import unit_cell;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import grid_pore_size;
import grid_percolation;
import grid_connected_components;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;


EnergyPoreSizeDistribution::EnergyPoreSizeDistribution() {}


EnergyPoreSizeDistribution::~EnergyPoreSizeDistribution() {}


namespace
{
// What one landscape yields before it is binned. Both reductions of the landscape go through here, so the
// two curves cannot be measured differently; they are binned together afterwards, over a range wide enough
// for both, so that they can be read off against each other row by row.
struct SweptField
{
  std::vector<float> distance;
  std::vector<float> poreRadius;
  std::vector<std::uint8_t> accessible;
  int dimensionality{0};

  double largestDiameter() const
  {
    double largest = 0.0;
    for (std::size_t v = 0; v < this->distance.size(); ++v)
    {
      if (this->distance[v] < 0.0f) continue;
      largest = std::max(largest, 2.0 * static_cast<double>(this->poreRadius[v]));
    }
    return largest;
  }
};

// The escape barrier out of every piece of a landscape, which does not depend on the level a region is cut
// at and so is computed once per field and handed to each cut of it.
std::vector<float> escapeOf(uint3 gridSize, std::span<const float> openness)
{
  GridPercolation swept = sweepPercolation(gridSize, openness, std::numeric_limits<float>::lowest(), true);
  return std::move(swept.escapeOpenness);
}

SweptField sweepField(const EnergyBackend &backend, const UnitCell &unitCell, uint3 gridSize,
                      std::span<const float> energy, double isoValue, double slack,
                      std::span<const float> escapeOpenness, double kT, double thresholdInKT)
{
  SweptField result;
  const std::size_t numberOfVoxels = energy.size();
  if (numberOfVoxels == 0) return result;

  // Everything downstream reads a field as an openness, larger meaning more open, so the landscape is turned
  // over once here and the level with it. A molecule is more welcome where the energy is lower.
  std::vector<float> openness(numberOfVoxels);
  for (std::size_t v = 0; v < numberOfVoxels; ++v) openness[v] = -energy[v];
  const double level = -isoValue;

  result.distance = distanceToIsosurface(gridSize, unitCell, openness, level);
  result.poreRadius = backend.poreRadiusField(gridSize, unitCell, result.distance, slack);

  // A cavity the molecule cannot get into contributes nothing to an uptake however wide it is, so the void
  // divides by whether a point can be reached. On a landscape that is not the same question as whether a
  // path exists at the level: a cage closed at the level but open a few kT above it is filled in any
  // experiment, and DDR is exactly that for methane. So a piece counts as reached if it runs through the
  // crystal at the level, or if the climb out of it is one the molecule will make.
  GridComponents components = GridComponents::compute(gridSize, openness, level);

  const float noWayOut = std::numeric_limits<float>::lowest();
  std::vector<float> escapeOfPore(components.pores.size(), noWayOut);
  if (!escapeOpenness.empty())
  {
    for (std::size_t v = 0; v < numberOfVoxels; ++v)
    {
      const std::int32_t pore = components.voxelPore[v];
      if (pore < 0) continue;
      escapeOfPore[static_cast<std::size_t>(pore)] =
          std::max(escapeOfPore[static_cast<std::size_t>(pore)], escapeOpenness[v]);
    }
  }

  std::vector<std::uint8_t> poreIsReached(components.pores.size(), 0);
  for (std::size_t poreIndex = 0; poreIndex < components.pores.size(); ++poreIndex)
  {
    const GridPore &pore = components.pores[poreIndex];
    if (pore.isChannel)
    {
      poreIsReached[poreIndex] = 1;
      continue;
    }

    const float escape = escapeOfPore[poreIndex];
    if (escape == noWayOut || kT <= 0.0) continue;

    const double barrier = (-static_cast<double>(escape) - (-pore.largestOpenness)) / kT;
    if (barrier <= thresholdInKT) poreIsReached[poreIndex] = 1;
  }

  result.accessible.assign(numberOfVoxels, 0);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    const std::int32_t pore = components.voxelPore[v];
    if (pore < 0) continue;
    if (poreIsReached[static_cast<std::size_t>(pore)]) result.accessible[v] = 1;
  }
  result.dimensionality = components.dimensionality;

  return result;
}

// How much of its time the molecule spends at each point of the void, up to a constant. Only points inside
// the void are given a weight, since the ones outside are not part of the curve at all; the constant is
// fixed by handing the deepest point of the cell a weight of one, which cancels out of every ratio the curve
// is made of and keeps exp() away from its limits at low temperature.
struct Occupancy
{
  std::vector<double> weight;

  // Over the whole cell, and over the part of it the sizes were measured on. The second is smaller because
  // the reference contour has to close somewhere, and the molecule does spend a little of its time on the
  // repulsive side of it. Their ratio says how much of the answer the curve is entitled to speak for.
  double total{0.0};
  double inside{0.0};
};

Occupancy boltzmannWeights(std::span<const float> energy, std::span<const float> distance, double temperature)
{
  Occupancy occupancy;
  const std::size_t numberOfVoxels = energy.size();
  occupancy.weight.assign(numberOfVoxels, 0.0);
  if (temperature <= 0.0) return occupancy;

  double deepest = std::numeric_limits<double>::max();
  for (std::size_t v = 0; v < numberOfVoxels; ++v) deepest = std::min(deepest, static_cast<double>(energy[v]));
  if (deepest == std::numeric_limits<double>::max()) return occupancy;

  const double beta = 1.0 / (Units::KB * temperature);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    const double w = std::exp(-beta * (static_cast<double>(energy[v]) - deepest));
    occupancy.total += w;
    if (distance[v] < 0.0f) continue;
    occupancy.inside += w;
    occupancy.weight[v] = w;
  }
  return occupancy;
}
}  // namespace


void EnergyPoreSizeDistribution::run(const EnergyBackend &backend, const PairInteractions &interactions,
                                     const Crystal &framework, const LinearProbe &probe, double isoValue,
                                     uint3 gridSize, std::size_t numberOfOrientations, double temperature,
                                     double thresholdInKT, bool useElectrostatics, double relativePrecision,
                                     std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins)
{
  MolecularField field = buildMolecularField(backend, interactions, framework, probe, gridSize, numberOfOrientations,
                                             temperature, useElectrostatics, relativePrecision);
  this->grid = field.grid;
  this->potential = field.potential;
  this->isoValue = isoValue;
  this->temperature = temperature;
  this->thresholdInKT = thresholdInKT;

  if (this->grid.numberOfVoxels() == 0) return;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  this->slack = coveringSlack(framework.unitCell, gridSize);

  const double kT = Units::KB * temperature;

  // What it costs to leave is a property of the landscape and not of the level, so one pass over each of the
  // two reductions serves every cut of them below.
  std::vector<float> openness(this->grid.numberOfVoxels());
  for (std::size_t v = 0; v < openness.size(); ++v) openness[v] = -this->grid.freeEnergy[v];
  const std::vector<float> freeEscape = escapeOf(gridSize, openness);

  for (std::size_t v = 0; v < openness.size(); ++v) openness[v] = -this->grid.minimumEnergy[v];
  const std::vector<float> leastEscape = escapeOf(gridSize, openness);

  SweptField free = sweepField(backend, framework.unitCell, gridSize, this->grid.freeEnergy, isoValue,
                               this->slack, freeEscape, kT, this->thresholdInKT);
  SweptField least = sweepField(backend, framework.unitCell, gridSize, this->grid.minimumEnergy, isoValue,
                                this->slack, leastEscape, kT, this->thresholdInKT);

  // The weighted curve is measured against the reference contour rather than against the chosen level, and
  // this is the whole reason it is worth having. Weighting the counting alone does not free the curve of the
  // level: the size at a point is the largest sphere that fits in the void, so moving the level moves every
  // size, and the axis slides out from under the weights. Pinning the contour holds the axis still, and the
  // curve is then a statement about the framework at a temperature and not about a number someone chose.
  SweptField reference = (isoValue == 0.0) ? SweptField{}
                                           : sweepField(backend, framework.unitCell, gridSize,
                                                        this->grid.freeEnergy, 0.0, this->slack, freeEscape, kT,
                                                        this->thresholdInKT);
  const SweptField &onReference = (isoValue == 0.0) ? free : reference;

  // The minimum-energy void contains the free-energy one, so its pores are the wider, and a level below zero
  // puts the reference contour outside them both. All the curves are binned over whichever reaches furthest,
  // since curves on different bins cannot be read off against each other.
  const double widest =
      std::max({free.largestDiameter(), least.largestDiameter(), onReference.largestDiameter(), 1.0});
  PoreSizeCurve freeCurve = binPoreSizes(free.poreRadius, free.distance, free.accessible, {},
                                         maximumDiameter.value_or(widest), numberOfBins);
  PoreSizeCurve leastCurve = binPoreSizes(least.poreRadius, least.distance, least.accessible, {},
                                          maximumDiameter.value_or(widest), numberOfBins);

  Occupancy occupancy = boltzmannWeights(this->grid.freeEnergy, onReference.distance, temperature);
  PoreSizeCurve occupancyCurve =
      binPoreSizes(onReference.poreRadius, onReference.distance, onReference.accessible, occupancy.weight,
                   maximumDiameter.value_or(widest), numberOfBins);

  // The volume the molecules are being compared against has to be the volume they are weighted over, or the
  // two means are answers to questions about different regions.
  PoreSizeCurve referenceCurve =
      (isoValue == 0.0) ? freeCurve
                        : binPoreSizes(onReference.poreRadius, onReference.distance, onReference.accessible, {},
                                       maximumDiameter.value_or(widest), numberOfBins);

  const double total = static_cast<double>(this->grid.numberOfVoxels());
  this->points = freeCurve.points;
  this->largestDiameter = freeCurve.largestDiameter;
  this->voidFraction = static_cast<double>(freeCurve.numberOfVoidVoxels) / total;
  this->accessibleVoidFraction = static_cast<double>(freeCurve.numberOfAccessibleVoxels) / total;
  this->dimensionality = free.dimensionality;

  this->occupancyPoints = occupancyCurve.points;
  this->volumetricMeanDiameter = referenceCurve.meanDiameter;
  this->occupancyMeanDiameter = occupancyCurve.meanDiameter;
  this->reachableOccupancyFraction =
      (occupancyCurve.totalWeight > 0.0) ? occupancyCurve.totalAccessibleWeight / occupancyCurve.totalWeight : 0.0;
  this->occupancyInsideReference = (occupancy.total > 0.0) ? occupancy.inside / occupancy.total : 0.0;

  this->minimumEnergyPoints = leastCurve.points;
  this->minimumEnergyLargestDiameter = leastCurve.largestDiameter;
  this->minimumEnergyVoidFraction = static_cast<double>(leastCurve.numberOfVoidVoxels) / total;

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  this->sweepSeconds = elapsed.count();
  this->seconds = this->sweepSeconds;

  double3 spacing = this->grid.spacing();
  const double volume = framework.unitCell.volume;

  std::ofstream myfile;
  myfile.open(framework.name + "." + probe.name + ".energy.psd." + this->grid.backend + ".txt");
  std::print(myfile, "# Energy-based pore-size distribution (Gelb-Gubbins on an iso-surface of the landscape)\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", volume);
  std::print(myfile, "# Crystal mass: {} [g/mol]\n", framework.mass);
  std::print(myfile, "# Molecule: {}, {} sites, {:.4f} [Å] end to end\n", probe.name, probe.sites.size(),
             probe.length());
  std::print(myfile, "# Orientations sampled: {} over the {}\n", this->grid.numberOfOrientations,
             this->grid.overHemisphere ? "hemisphere, the molecule being the same end for end" : "whole sphere");
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å], far face left out\n",
             gridSize.x, gridSize.y, gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Temperature: {} [K]\n", temperature);
  std::print(myfile, "# Cutoff: {} [Å]\n", this->grid.cutOff);

  if (this->grid.chargesIncluded)
  {
    std::print(myfile, "# Electrostatics: Ewald, split at alpha = {:.5f} [1/Å], {} wave vectors\n",
               this->potential.alpha, this->potential.numberOfWaveVectors);
  }
  if (this->grid.chargesIgnored)
  {
    std::print(myfile, "#\n");
    std::print(myfile, "# WARNING: this molecule carries partial charges and they have not been acted on.\n");
  }

  std::print(myfile, "# Iso-value: {} [internal], {:.4f} [K]\n", this->isoValue,
             this->isoValue * Units::EnergyToKelvin);
  std::print(myfile, "# A sphere covers a point when it comes within {:.5f} [Å] of it, half a voxel diagonal "
                     "twice\n", this->slack);
  std::print(myfile, "# Diameters evaluated: {}, up to {} [Å]\n", this->points.size(),
             this->points.empty() ? 0.0 : this->points.back().diameter + 0.5 * freeCurve.binWidth);
  std::print(myfile, "# Timing ({}): {} [s] for the landscape\n", this->grid.backend, this->grid.seconds);
  std::print(myfile, "# Timing ({}): {} [s] for the two sweeps\n", backend.name, this->sweepSeconds);
  std::print(myfile, "#\n");

  std::print(myfile, "{:44} {:14.8f}\n", "Void fraction at this level (free energy):", this->voidFraction);
  std::print(myfile, "{:44} {:14.8f}\n", "Of the cell, the molecule can reach:", this->accessibleVoidFraction);
  std::print(myfile, "# Reached means the piece runs through the crystal at this level, or the climb out of\n");
  std::print(myfile, "# it is under {} kT and the molecule makes it. The second half matters: a cage shut at\n",
             this->thresholdInKT);
  std::print(myfile, "# the level but open a little above it holds gas in any experiment, and counting it out\n");
  std::print(myfile, "# because no path exists at one chosen level would throw that uptake away.\n");
  std::print(myfile, "{:44} {:14.5f} [Å]\n", "Largest sphere the void holds:", this->largestDiameter);
  std::print(myfile, "{:44} {:14}\n", "Directions the void runs in:", this->dimensionality);
  std::print(myfile, "{:44} {:14.8f}\n",
             "Void fraction at this level (minimum energy):", this->minimumEnergyVoidFraction);
  std::print(myfile, "{:44} {:14.5f} [Å]\n",
             "Largest sphere that void holds:", this->minimumEnergyLargestDiameter);
  std::print(myfile, "\n");

  std::print(myfile, "{:44} {:14.5f} [Å]\n", "Mean pore, per unit of volume:", this->volumetricMeanDiameter);
  std::print(myfile, "{:44} {:14.5f} [Å]\n", "Mean pore, per molecule:", this->occupancyMeanDiameter);
  std::print(myfile, "{:44} {:14.8f}\n",
             "Of the molecules, in reachable void:", this->reachableOccupancyFraction);
  std::print(myfile, "{:44} {:14.8f}\n",
             "Of the molecules, inside the contour:", this->occupancyInsideReference);
  std::print(myfile, "\n");

  std::print(myfile, "# The size at a point is the largest sphere that fits in the void and covers the point, so\n");
  std::print(myfile, "# a narrow neck between two wide cavities is given the width of a cavity rather than its\n");
  std::print(myfile, "# own: that is what the definition says, and it is why the curve has no weight at all\n");
  std::print(myfile, "# below the narrowest sphere the void can hold anywhere.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# The void here is bounded by a level set of the energy and not by an atom's surface, so\n");
  std::print(myfile, "# the room at each point has to be measured off the field rather than known exactly. It is\n");
  std::print(myfile, "# found by walking out to where the field crosses the level and interpolating the crossing\n");
  std::print(myfile, "# between the two grid points that straddle it, the same construction marching cubes puts\n");
  std::print(myfile, "# its vertices by. That leaves an error of up to about half a voxel in each distance,\n");
  std::print(myfile, "# which the geometric route does not have, and refining the grid is the only thing that\n");
  std::print(myfile, "# removes it.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# Columns 8 and 9 count each point by exp(-A/kT) at {} K instead of by the volume it takes\n",
             this->temperature);
  std::print(myfile, "# up. Columns 2 to 7 say how the room is divided among pore widths; columns 8 and 9 say how\n");
  std::print(myfile, "# the molecules are, and it is the latter that a measured isotherm is a statement about.\n");
  std::print(myfile, "# Where the two differ most is at the narrow end, which holds little volume and, being\n");
  std::print(myfile, "# where the walls close in from both sides, much of the binding.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# Columns 8 and 9 carry no iso-value. They are measured against the fixed A = 0 contour,\n");
  std::print(myfile, "# where the molecule is no better off than in the gas, and not against the level above.\n");
  std::print(myfile, "# Weighting alone would not have been enough: the size at a point is the width of the\n");
  std::print(myfile, "# largest sphere that fits, so raising the level widens every sphere and slides the axis,\n");
  std::print(myfile, "# and weights cannot hold still an axis that has moved. Pinning the contour is what makes\n");
  std::print(myfile, "# these two columns depend on the temperature and on nothing else that was chosen.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# The molecule does spend a little of its time outside that contour, pressed into the\n");
  std::print(myfile, "# wall, where there is no size to give it. That share is reported above, and columns 8\n");
  std::print(myfile, "# and 9 are normalised over the rest.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# column 1: diameter d [Å]\n");
  std::print(myfile, "# column 2: share of the void in pores at least d across, free energy\n");
  std::print(myfile, "# column 3: the distribution, the derivative of column 2 [1/Å]\n");
  std::print(myfile, "# column 4: share of the reachable void in pores at least d across\n");
  std::print(myfile, "# column 5: the derivative of column 4 [1/Å]\n");
  std::print(myfile, "# column 6: share of the void in pores at least d across, minimum energy\n");
  std::print(myfile, "# column 7: the derivative of column 6 [1/Å]\n");
  std::print(myfile, "# column 8: share of the molecules in pores at least d across\n");
  std::print(myfile, "# column 9: the derivative of column 8 [1/Å]\n");
  std::print(myfile, "#      d [Å]     cumulative   distribution     accessible   distribution      least cum. "
                     "  least dist.      occupancy    occ. dist.\n");
  for (std::size_t bin = 0; bin < this->points.size(); ++bin)
  {
    const PoreSizeCurvePoint &point = this->points[bin];
    const double leastCumulative =
        (bin < this->minimumEnergyPoints.size()) ? this->minimumEnergyPoints[bin].cumulativeVoidFraction : 0.0;
    const double leastDistribution =
        (bin < this->minimumEnergyPoints.size()) ? this->minimumEnergyPoints[bin].distribution : 0.0;
    const double occupancyCumulative =
        (bin < this->occupancyPoints.size()) ? this->occupancyPoints[bin].cumulativeVoidFraction : 0.0;
    const double occupancyDistribution =
        (bin < this->occupancyPoints.size()) ? this->occupancyPoints[bin].distribution : 0.0;
    std::print(myfile, "{:11.6f} {:14.8f} {:14.8f} {:14.8f} {:14.8f} {:14.8f} {:14.8f} {:14.8f} {:14.8f}\n",
               point.diameter, point.cumulativeVoidFraction, point.distribution, point.cumulativeAccessibleFraction,
               point.accessibleDistribution, leastCumulative, leastDistribution, occupancyCumulative,
               occupancyDistribution);
  }

  myfile.close();
}
