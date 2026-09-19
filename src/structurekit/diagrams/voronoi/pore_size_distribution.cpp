module;

module voronoi_pore_size_distribution;

import std;

import double3;
import units;
import crystal;
import pair_interactions;
import skspacegroupdatabase;
import pore_accessibility;
import exact_pore_size_distribution;
import exact_void_split;
import voronoi_blocking_spheres;

// The bare radii of the framework's atoms, and where they are in fractional coordinates.
std::pair<std::vector<double3>, std::vector<double>> frameworkSpheres(const PairInteractions& interactions,
                                                                     const Crystal& framework)
{
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
  return {fractionalPositions, radii};
}


// A quantity that rounds away to nothing, printed as nothing rather than as a signed nothing.
//
// Past the largest pore of a framework there is no volume left, and what the arithmetic hands back there is the
// difference of two numbers some twelve orders of magnitude larger: a few parts in 1e11, of either sign,
// according to which way the last bits of the two fell. Every column below is printed to eight decimals at the
// finest, so such a number is already nothing in the table; all that reaches the reader is the minus sign in
// front of it, which is round-off wearing the look of a measurement. Nothing else is touched, the threshold
// being smaller than the least the table can show.
double asPrinted(double value)
{
  constexpr double finestColumn = 1.0e-8;
  return (std::abs(value) < 0.5 * finestColumn) ? 0.0 : value;
}


void writePoreSizeDistribution(const Crystal& framework, const std::string& diagramName,
                               const std::string& probePseudoAtom, const PoreSizeDistributionCurve& curve)
{
  std::ofstream report;
  report.open(std::format("{}.{}.psd.txt", framework.name, diagramName));

  std::print(report, "# Pore-size distribution ({}, exact)\n", diagramName);
  std::print(report, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(report, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(report, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(report, "# Crystal volume: {} [Å³]\n", curve.cellVolume);
  std::print(report, "# Crystal mass: {} [g/mol]\n", framework.mass);
  std::print(report, "# Void volume: {} [Å³], void fraction {}\n", curve.voidVolume,
             (curve.cellVolume > 0.0) ? curve.voidVolume / curve.cellVolume : 0.0);
  std::print(report, "# Probe (blocking / accessibility): {}, radius {} [Å], diameter {} [Å]\n", probePseudoAtom,
             curve.probeRadius, 2.0 * curve.probeRadius);
  std::print(report, "# Floor radius: {} [Å] ({})\n", curve.floorRadius,
             (curve.floorRadius <= 0.0) ? "hybrid: He pockets filled by blocking spheres, bare PSD of the rest"
                                       : "probe-occupiable: flat below this diameter");
  std::print(report, "# Volume of the primary curve: {} [Å³], {} of the void\n", curve.probeAccessibleVolume,
             (curve.voidVolume > 0.0) ? curve.probeAccessibleVolume / curve.voidVolume : 0.0);
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
  std::print(report, "# Pockets sealed to the blocking probe are excluded from columns 2 to 5. By default the\n");
  std::print(report, "# floor is zero (hybrid): those pockets are filled by their helium blocking spheres and the\n");
  std::print(report, "# bare pore-size distribution of the remaining void is reported, so wall corrugation below\n");
  std::print(report, "# the blocking probe's diameter is still shown. With a positive floor the primary curve is\n");
  std::print(report, "# instead the probe-occupiable distribution, flat below twice that radius.\n");
  std::print(report, "#\n");
  if (curve.floorRadius > 0.0)
  {
    std::print(report, "# The primary curve is flat below the floor diameter {:.5f} [Å] and normalised by the\n",
               2.0 * curve.floorRadius);
    std::print(report, "# volume left after blocking. Columns 6 to 8 divide each row by what a probe of that\n");
    std::print(report, "# row's own diameter can reach, which moves along the curve where the blocking probe\n");
    std::print(report, "# stands still. Columns 14 to 17 hold the bare whole-void curve from the same sweep.\n");
  }
  else
  {
    std::print(report, "# The primary curve has no diameter floor: helium-sealed pockets are filled by blocking\n");
    std::print(report, "# spheres and the remaining open network is given its bare Gelb--Gubbins distribution,\n");
    std::print(report, "# normalised by that remaining void. The open framework without blocking is not\n");
    std::print(report, "# reported: sealed cages are not an adsorption observable. Columns 14 to 17 mirror\n");
    std::print(report, "# columns 2 to 5.\n");
  }
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
  std::print(report, "# own interval by a factor of a thousand or more. The list below is for the blocked curve:\n");
  std::print(report, "#\n");
  std::print(report, "#      d [Å]        weight     cornered within [Å]\n");
  if (curve.probeAccessibleVolume > 0.0)
  {
    for (const PoreSizeSpike& spike : curve.probeAccessibleSpikes)
    {
      std::print(report, "#   {:9.5f}    {:10.6f}    {:.2e}\n", spike.diameter, spike.weight, spike.bracket);
    }
    std::print(report, "#\n");
    std::print(report, "# The continuous part integrates to {:.6f} and the spikes weigh {:.6f}, which come to\n",
               curve.probeAccessibleIntegral, curve.probeAccessibleSingularWeight);
    std::print(report, "# {:.6f} against one over the volume the probe can reach.\n",
               curve.probeAccessibleIntegral + curve.probeAccessibleSingularWeight);
    std::print(report, "# The largest sphere that fits anywhere the probe can reach is {:.5f} [Å] across, with\n",
               curve.probeAccessibleLargestDiameter);
    std::print(report, "# {:.2e} of that volume left beyond the end of the range.\n",
               curve.probeAccessibleTruncatedWeight);
  }
  else
  {
    std::print(report, "#   (none: every pore is sealed to this probe, so columns 2 to 5 are zero throughout)\n");
    std::print(report, "#\n");
    std::print(report, "# The probe reaches nothing of this framework at that size. That is a statement about the\n");
    std::print(report, "# framework and not a failure of the analysis.\n");
  }
  std::print(report, "#\n");
  if (curve.floorRadius > 0.0)
  {
    std::print(report, "# The same three things for the whole void, with pockets not blocked:\n");
    std::print(report, "#\n");
    std::print(report, "#   whole void   d [Å]        weight     cornered within [Å]\n");
    for (const PoreSizeSpike& spike : curve.spikes)
    {
      std::print(report, "#   {:9.5f}    {:10.6f}    {:.2e}\n", spike.diameter, spike.weight, spike.bracket);
    }
    std::print(report, "#\n");
    std::print(report, "# Its continuous part comes to {:.6f} and its spikes to {:.6f}, together {:.6f} against one.\n",
               curve.integral, curve.singularWeight, curve.integral + curve.singularWeight);
    std::print(report, "# The largest sphere that fits in the void is {:.5f} [Å] across, with {:.2e} left beyond\n",
               curve.largestDiameter, curve.truncatedWeight);
    std::print(report, "# the end of the range.\n");
  }
  else
  {
    std::print(report, "# The open whole-void distribution is not computed for the hybrid (no floor).\n");
  }
  std::print(report, "#\n");
  std::print(report, "# column 1: diameter d [Å]\n");
  std::print(report, "# column 2: P(d) [1/Å] over the open network (pockets blocked), normalised by its volume\n");
  std::print(report, "# column 3: that volume at this diameter, as a fraction of the whole of it\n");
  std::print(report, "# column 4: that volume at this diameter [Å³]\n");
  std::print(report, "# column 5: the same [cm³/g]\n");
  std::print(report, "# columns 6-8: the share of P(d) on surfaces a probe of that row's own diameter can\n");
  std::print(report, "#              reach, on surfaces sealed off from it, and on surfaces the network could\n");
  std::print(report, "#              not place\n");
  std::print(report, "# columns 9-13: arcs, of them cusped, vertices, of them clipped, of them degenerate\n");
  if (curve.floorRadius > 0.0)
  {
    std::print(report, "# column 14: P(d) [1/Å] over the whole void, normalised to integrate to one\n");
    std::print(report, "# column 15: cumulative pore volume, as a fraction of the void volume\n");
    std::print(report, "# column 16: cumulative pore volume [Å³]\n");
    std::print(report, "# column 17: the same [cm³/g]\n");
  }
  else
  {
    std::print(report, "# columns 14-17: mirror of columns 2-5 (no separate whole-void curve)\n");
  }

  // Å³ per unit cell to cm³ per gram.
  const double gramsPerCell = framework.mass / Units::AvogadroConstant;
  const double toVolumePerMass = (gramsPerCell > 0.0) ? 1.0e-24 / gramsPerCell : 0.0;

  for (const PoreSizeDistributionPoint& point : curve.points)
  {
    // The diameter to more figures than a plot needs, because the columns are meant to be differenced
    // against one another and a spacing read back from a rounded abscissa is not the spacing used.
    std::print(report,
               "{:11.6f} {:14.8f} {:12.8f} {:14.5f} {:12.6f} {:14.8f} {:14.8f} {:14.8f} {:8} {:8} {:8} {:8} {:8}"
               " {:14.8f} {:12.8f} {:14.5f} {:12.6f}\n",
               point.diameter, asPrinted(point.probeAccessibleDistribution),
               asPrinted(point.probeAccessibleCumulative), asPrinted(point.probeAccessiblePoreVolume),
               asPrinted(point.probeAccessiblePoreVolume * toVolumePerMass), asPrinted(point.accessible),
               asPrinted(point.inaccessible), asPrinted(point.undecided), point.numberOfArcs, point.cuspedArcs,
               point.numberOfVertices, point.clippedVertices, point.degenerateVertices,
               asPrinted(point.distribution), asPrinted(point.cumulative), asPrinted(point.poreVolume),
               asPrinted(point.poreVolume * toVolumePerMass));
  }
  report.close();
}



void writePoreSizePeaks(const Crystal& framework, const std::string& diagramName,
                        const std::string& probePseudoAtom, const PoreSizeDistributionCurve& curve)
{
  std::ofstream report;
  report.open(std::format("{}.{}.psd-peaks.txt", framework.name, diagramName));

  std::print(report, "# Pore-size peaks ({}, exact): room diameters and void fraction at each\n", diagramName);
  std::print(report, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(report, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(report, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(report, "# Crystal volume: {} [Å³]\n", curve.cellVolume);
  std::print(report, "# Crystal mass: {} [g/mol]\n", framework.mass);
  std::print(report, "# Void volume: {} [Å³], void fraction {}\n", curve.voidVolume,
             (curve.cellVolume > 0.0) ? curve.voidVolume / curve.cellVolume : 0.0);
  std::print(report, "# Probe (blocking / accessibility): {}, radius {} [Å], diameter {} [Å]\n", probePseudoAtom,
             curve.probeRadius, 2.0 * curve.probeRadius);
  std::print(report, "# Floor radius: {} [Å] ({})\n", curve.floorRadius,
             (curve.floorRadius <= 0.0) ? "hybrid: He pockets filled by blocking spheres, bare peaks of the rest"
                                       : "probe-occupiable: peaks only above this diameter");
  std::print(report, "# Volume of the primary curve: {} [Å³], {} of the void\n", curve.probeAccessibleVolume,
             (curve.voidVolume > 0.0) ? curve.probeAccessibleVolume / curve.voidVolume : 0.0);
  std::print(report, "# Seconds: {}, over {} evaluations of the surface\n", curve.seconds, curve.numberOfEvaluations);
  std::print(report, "#\n");
  std::print(report, "# Only the discrete room sizes of the Gelb--Gubbins measure: each spike is the diameter of\n");
  std::print(report, "# a family of maximal inscribed spheres and the fraction of the (reachable) void that sits\n");
  std::print(report, "# at exactly that size. Wall corrugation between rooms is omitted; use --pore-size-distribution\n");
  std::print(report, "# for the continuous curve as well. Diameters come from the pore-network maxima; each distinct\n");
  std::print(report, "# size is sampled once just below and once just above the cliff (no dense diameter sweep).\n");
  std::print(report, "#\n");
  std::print(report, "# Rows sorted by weight descending.\n");
  std::print(report, "# column 1: diameter d [Å]\n");
  std::print(report, "# column 2: weight (fraction of the primary void at this pore size)\n");
  std::print(report, "# column 3: volume at this pore size [Å³]\n");
  std::print(report, "# column 4: the same [cm³/g]\n");
  std::print(report, "# column 5: bracket width the diameter was cornered within [Å]\n");

  const double gramsPerCell = framework.mass / Units::AvogadroConstant;
  const double toVolumePerMass = (gramsPerCell > 0.0) ? 1.0e-24 / gramsPerCell : 0.0;
  const double normVolume = curve.probeAccessibleVolume;

  auto byWeightDescending = [](const PoreSizeSpike& a, const PoreSizeSpike& b) { return a.weight > b.weight; };

  if (normVolume <= 0.0)
  {
    std::print(report, "# (none: every pore is sealed to this probe)\n");
  }
  else
  {
    std::vector<PoreSizeSpike> peaks = curve.probeAccessibleSpikes;
    std::ranges::sort(peaks, byWeightDescending);
    for (const PoreSizeSpike& spike : peaks)
    {
      const double volume = spike.weight * normVolume;
      std::print(report, "{:11.6f} {:14.8f} {:14.5f} {:12.6f} {:.2e}\n", spike.diameter, spike.weight, volume,
                 volume * toVolumePerMass, spike.bracket);
    }
    std::print(report, "#\n");
    std::print(report, "# {} peaks, singular weight {:.6f} of the primary void. Largest room {:.5f} [Å].\n",
               curve.probeAccessibleSpikes.size(), curve.probeAccessibleSingularWeight,
               curve.probeAccessibleLargestDiameter);
  }

  if (curve.floorRadius > 0.0)
  {
    std::print(report, "#\n");
    std::print(report, "# Whole void (pockets not blocked):\n");
    std::print(report, "# Rows sorted by weight descending.\n");
    std::print(report, "# column 1: diameter d [Å]\n");
    std::print(report, "# column 2: weight (fraction of the void)\n");
    std::print(report, "# column 3: volume [Å³]\n");
    std::print(report, "# column 4: [cm³/g]\n");
    std::print(report, "# column 5: bracket [Å]\n");
    std::vector<PoreSizeSpike> peaks = curve.spikes;
    std::ranges::sort(peaks, byWeightDescending);
    for (const PoreSizeSpike& spike : peaks)
    {
      const double volume = spike.weight * curve.voidVolume;
      std::print(report, "{:11.6f} {:14.8f} {:14.5f} {:12.6f} {:.2e}\n", spike.diameter, spike.weight, volume,
                 volume * toVolumePerMass, spike.bracket);
    }
    std::print(report, "#\n");
    std::print(report, "# {} whole-void peaks, singular weight {:.6f}. Largest {:.5f} [Å].\n", curve.spikes.size(),
               curve.singularWeight, curve.largestDiameter);
  }

  report.close();
}


PoreSizeDistributionCurve hybridFromBlockedCurve(PoreSizeDistributionCurve blocked, double accessibilityRadius)
{
  PoreSizeDistributionCurve curve = std::move(blocked);
  curve.probeRadius = accessibilityRadius;
  curve.floorRadius = 0.0;
  // Bare PSD of the blocked network is the only curve: mirror it into the primary columns.
  // voidVolume is the remaining open-network volume after blocking spheres (zero-probe), which
  // also normalises the distribution; the He void-split volume is a different quantity and is not
  // substituted here.
  curve.probeAccessibleVolume = curve.voidVolume;
  curve.probeAccessibleIntegral = curve.integral;
  curve.probeAccessibleSingularWeight = curve.singularWeight;
  curve.probeAccessibleSpikes = curve.spikes;
  curve.probeAccessibleTruncatedWeight = curve.truncatedWeight;
  curve.probeAccessibleLargestDiameter = curve.largestDiameter;
  for (PoreSizeDistributionPoint& point : curve.points)
  {
    point.probeAccessiblePoreVolume = point.poreVolume;
    point.probeAccessibleCumulative = point.cumulative;
    point.probeAccessibleDistribution = point.distribution;
  }
  return curve;
}


void VoronoiPoreSizeDistribution::run(const PairInteractions& interactions, const Crystal& framework,
                                      std::string probePseudoAtom, std::optional<double> maximumDiameter,
                                      std::optional<std::size_t> numberOfBins, std::size_t subdivisions,
                                      double floorRadius, bool peaksOnly)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("VoronoiPoreSizeDistribution: Unknown probe-atom type\n");
  }
  const double probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

  auto [fractionalPositions, radii] = frameworkSpheres(interactions, framework);

  auto build = [&](double inflation)
  { return PoreAccessibility::create(framework.unitCell, fractionalPositions, radii, inflation); };

  const double cellVolume = framework.unitCell.volume;
  const double maxDiameter = maximumDiameter.value_or(20.0);
  const std::size_t bins = numberOfBins.value_or(100);

  if (floorRadius > 0.0)
  {
    curve = peaksOnly ? exactPoreSizePeaks(build, cellVolume, subdivisions, probeRadius, floorRadius)
                      : exactPoreSizeDistribution(build, cellVolume, maxDiameter, bins, subdivisions, probeRadius,
                                                  floorRadius);
  }
  else
  {
    // Hybrid: fill helium-sealed pockets with their blocking spheres, then take the bare PSD of what remains.
    PoreAccessibility he = PoreAccessibility::create(framework.unitCell, fractionalPositions, radii, probeRadius);
    ExactVoidSplit split = exactVoidSplitByComponents(he, cellVolume, subdivisions);
    std::vector<BlockingSphere> spheres;
    if (measuredSpheresRefused(split).empty())
    {
      spheres = exactBlockingSpheres(split);
    }
    else
    {
      spheres = computeBlockingSpheres(he, static_cast<std::size_t>(200.0 * cellVolume));
    }

    std::vector<double3> blockedPositions = fractionalPositions;
    std::vector<double> blockedRadii = radii;
    blockedPositions.reserve(fractionalPositions.size() + spheres.size());
    blockedRadii.reserve(radii.size() + spheres.size());
    for (const BlockingSphere& sphere : spheres)
    {
      blockedPositions.push_back(sphere.centerFractional);
      blockedRadii.push_back(sphere.radius);
    }

    auto buildBlocked = [&](double inflation)
    { return PoreAccessibility::create(framework.unitCell, blockedPositions, blockedRadii, inflation); };

    PoreSizeDistributionCurve blocked =
        peaksOnly ? exactPoreSizePeaks(buildBlocked, cellVolume, subdivisions, 0.0, 0.0)
                  : exactPoreSizeDistribution(buildBlocked, cellVolume, maxDiameter, bins, subdivisions, 0.0, 0.0);
    curve = hybridFromBlockedCurve(std::move(blocked), probeRadius);
  }

  if (peaksOnly)
    writePoreSizePeaks(framework, "voronoi", probePseudoAtom, curve);
  else
    writePoreSizeDistribution(framework, "voronoi", probePseudoAtom, curve);
}
