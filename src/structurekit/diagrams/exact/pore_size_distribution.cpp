module;

module exact_pore_size_distribution;

import std;

import double3;
import unit_cell;
import pore_accessibility;
import exact_surface_patches;
import exact_union_volume;
import exact_boundary_components;
import structure_parallel;
import exact_solvent_excluded;
import voronoi_network;

namespace
{

// One diameter, and the two things known there: the fraction with a pore size at least this, and the derivative
// of that fraction. They are computed by routes with nothing in common, which is what makes the comparison below
// worth making. Both are kept twice, once over the whole void and once over the part left after blocking.
struct Series
{
  double cumulative{0.0};
  double distribution{0.0};
};

struct Sample
{
  double diameter{0.0};
  Series whole;
  Series reachable;
};

// A spike has to be worth something to be worth looking for, and below this it would be lost in the trapezium
// error anyway.
constexpr double spikeFloor = 1.0e-5;

// The volume the derivative accounts for between two samples, by the trapezium rule, and the volume actually
// lost there. A smooth stretch has the two agreeing to the truncation error of the rule; a cliff does not.
double unaccounted(const Sample& left, const Sample& right, Series Sample::* of)
{
  const Series& a = left.*of;
  const Series& b = right.*of;
  double lost = a.cumulative - b.cumulative;
  double accounted = 0.5 * (a.distribution + b.distribution) * (right.diameter - left.diameter);
  return lost - accounted;
}

// Everything one series of the sweep collects: the integral of its continuous part, the spikes, and the end of
// the last interval a spike was found in, which is what tells a cliff the bisection straddled from two of them.
struct Collector
{
  double integral{0.0};
  double singularWeight{0.0};
  std::vector<PoreSizeSpike> spikes;
  double lastEnd{-1.0};
};

// The trapezium of one refined interval, and the spike where an excess survived the refinement. The excess is
// not allowed past `allowance`, and the weight the spike was given is returned so that one series can hold
// another to it. Nothing is capped in the whole void, whose allowance is unbounded.
double collect(Collector& into, const Sample& left, const Sample& right, Series Sample::* of, double allowance)
{
  const double gap = std::min(unaccounted(left, right, of), allowance);
  into.integral += 0.5 * ((left.*of).distribution + (right.*of).distribution) * (right.diameter - left.diameter);
  if (gap <= spikeFloor) return 0.0;

  PoreSizeSpike spike;
  spike.diameter = 0.5 * (left.diameter + right.diameter);
  spike.weight = gap;
  spike.bracket = right.diameter - left.diameter;

  // A cliff the bisection has straddled leaves an excess on either side of itself, and those two are one
  // cliff rather than two. Anything within sixteen brackets of the last is taken to be the same one, which is
  // a stated resolution rather than a tolerance: two pore sizes that near each other are the same pore size,
  // and two families that are really distinct are further apart than this by orders of magnitude.
  if (!into.spikes.empty() &&
      (into.lastEnd == left.diameter ||
       spike.diameter - into.spikes.back().diameter < 16.0 * (spike.bracket + into.spikes.back().bracket)))
  {
    PoreSizeSpike& previous = into.spikes.back();
    double total = previous.weight + spike.weight;
    previous.diameter = (previous.diameter * previous.weight + spike.diameter * spike.weight) / total;
    previous.bracket += spike.bracket;
    previous.weight = total;
  }
  else
  {
    into.spikes.push_back(spike);
  }
  into.lastEnd = right.diameter;
  into.singularWeight += gap;
  return gap;
}

// What a diameter is measured against, and the measurement itself.
//
// Two volumes normalise the curve, and neither is known until it has been measured: the whole of the void, at a
// vanishing probe, and the part left after blocking pockets sealed to the accessibility probe. The sampler is
// filled in that order --- `measureVoidVolume`, then `measureProbe`, then any number of calls to `at` --- and
// each of the first two leaves behind what the ones after it need. Past that the samples are independent of one
// another: `at` reads the sampler and writes nothing but the count of evaluations.
//
// `accessibilityRadius` is the probe that decides which pockets are blocked. `floorRadius` is the diameter floor
// of the primary curve: below 2 floorRadius that curve is held flat (the probe-occupiable convention). Zero
// floor is the hybrid: blocked pockets are still excluded, but the remaining void keeps its bare pore sizes.
struct Sampler
{
  Sampler(const std::function<PoreAccessibility(double)>& build, std::size_t subdivisions,
          double accessibilityRadius, double floorRadius, bool openNetwork)
      : build(build),
        subdivisions(subdivisions),
        accessibilityRadius(std::max(0.0, accessibilityRadius)),
        floorRadius(std::max(0.0, floorRadius)),
        openNetwork(openNetwork)
  {
  }

  const std::function<PoreAccessibility(double)>& build;
  std::size_t subdivisions{1};
  double accessibilityRadius{0.0};
  double floorRadius{0.0};

  // Whether the void is taken to be the open network alone: what a point can leave through the crystal, the
  // void shut in at a vanishing probe being left out of every volume here, the one that normalises included.
  // The hybrid asks for it. Its sealed pockets are filled before the sweep, but only the ones the probe's
  // centre can stand in: a cavity narrower than the probe holds no position of it, is no pocket of its
  // network and is given no sphere, and swept bare it came back as void the probe reaches.
  bool openNetwork{false};

  // The network at a vanishing probe, kept when the open network is asked for: it is what says, at every
  // diameter, whether the void a surface faces runs away through the crystal or is shut in. It is at a probe
  // no larger than any row's, which is all a reference has to be.
  PoreAccessibility bare;

  // The network of the accessibility probe. It is built once and kept: every diameter above the floor asks
  // the same network which of its pores each surface faces.
  PoreAccessibility reference;

  // The diameter of the largest sphere that fits in the void, from the network at vanishing probe. Rows of the
  // report past that diameter hold no volume, and evaluating them would rebuild the whole geometry for nothing
  // but a signed zero: they are left at the zeros they start as.
  double largestIncludedDiameter{0.0};

  // The largest sphere that fits where the primary curve's void is, from the network that decides which void
  // that is. It is what the curve's largest pore is reported as when the range ends before it.
  double largestReachableDiameter{0.0};

  // Distinct maximal-sphere diameters from that same vanishing-probe network (peaks-only route).
  std::vector<double> candidateDiameters;
  // The smallest diameter merged into each candidate, which is the largest of its family.
  std::vector<double> candidateLowestDiameters;

  double scale{0.0};                  // one over the void volume, which the whole curve is normalised by
  double probeAccessibleVolume{0.0};  // and the blocked remaining volume, which the primary curve is

  // Atomic because the diameters may be evaluated several at a time, and this is the one thing they all
  // write to. It is a count and not a sum, so how the increments interleave cannot change it.
  std::atomic<std::size_t> evaluations{0};

  double reachableScale() const { return (probeAccessibleVolume > 0.0) ? 1.0 / probeAccessibleVolume : 0.0; }

  // The two diameters a candidate's cliff lies between: from below the smallest member of its family to above
  // the largest, by a half width that is kept clear of the families either side. Nothing where the candidate
  // is past the largest sphere in the void.
  std::optional<std::pair<double, double>> candidateBracket(std::size_t index) const;

  // Volume behind components sealed or unplaced relative to the network `against`, and the reentrant
  // derivative on the reachable side. Shared by the normalisation sample and every diameter of the sweep.
  void splitByReference(const PoreAccessibility& accessibility, const BoundaryComponents& components,
                        const SolventExcludedGeometry& geometry, const PoreAccessibility& against, double& sealed,
                        double& unplaced, double& reachableDerivative) const;

  // The pore volume and its derivative that the sweep is over: all of it, or the open network's share.
  struct Swept
  {
    double poreVolume{0.0};
    double distribution{0.0};
  };
  Swept swept(const PoreAccessibility& accessibility, const BoundaryComponents& components,
              const SolventExcludedGeometry& geometry) const;

  double measureVoidVolume(double cellVolume);
  Sample measureProbe(double voidVolume);
  Sample at(double diameter, PoreSizeDistributionPoint* row);
};

// Symmetry-equivalent network maxima can differ by a few 1e-3 Å; those are one room, not many.
constexpr double peaksDiameterMergeTol = 1.0e-2;

// Tight brackets: the network already supplies the diameter, so the two samples only need to sit on
// opposite sides of the cliff. Continuous loss over a few 1e-3 Å is below the spike floor; a true room
// loses a finite fraction at once and is kept. One family after merging → one bracket, so no double count.
constexpr double candidateHalfWidth = 5.0e-3;

std::optional<std::pair<double, double>> Sampler::candidateBracket(std::size_t index) const
{
  const double diameter = candidateDiameters[index];
  if (diameter <= 0.0 || diameter > largestIncludedDiameter + candidateHalfWidth) return std::nullopt;

  const double lowest = candidateLowestDiameters[index];
  double half = candidateHalfWidth;
  if (index > 0) half = std::min(half, 0.25 * (lowest - candidateDiameters[index - 1]));
  if (index + 1 < candidateDiameters.size())
    half = std::min(half, 0.25 * (candidateLowestDiameters[index + 1] - diameter));
  half = std::max(half, 1.0e-4);
  return std::pair{lowest - half, diameter + half};
}

// Below the probe's diameter the surfaces cannot divide the void by that probe: a pocket whose windows are
// still open at the row's radius is one surface with the channel outside it, kept whole until the inflated
// atoms close the windows and dropped whole after. So a floor short of a positive probe is a curve with a
// cliff at the windows' size and nothing above it to stand on, and it is refused here rather than drawn. The
// drivers reach that range through a blocked structure instead, at a vanishing probe.
void requireFloorAtProbe(const char* route, double probeRadius, double floorRadius)
{
  if (probeRadius > 0.0 && floorRadius < probeRadius)
  {
    throw std::invalid_argument(std::format(
        "{}: the floor radius {} lies below the probe radius {}; below the probe's diameter the surfaces cannot "
        "divide the void by that probe, so the hybrid is taken through a blocked structure at a vanishing probe "
        "(see the pore-size-distribution drivers) and a direct call needs a floor at or above the probe\n",
        route, floorRadius, probeRadius));
  }
}

void Sampler::splitByReference(const PoreAccessibility& accessibility, const BoundaryComponents& components,
                               const SolventExcludedGeometry& geometry, const PoreAccessibility& against,
                               double& sealed, double& unplaced, double& reachableDerivative) const
{
  sealed = 0.0;
  unplaced = 0.0;
  reachableDerivative = 0.0;

  // Not `surfaceSides`, and the difference is the whole point of this curve rather than an oversight.
  // That division asks whether a surface seals void, and answers from the sign of the volume the
  // surface encloses before it consults anything: a surface that closes around void has sealed it, at
  // this radius. Here the question is whether the accessibility probe can reach that void, which is a
  // question about that probe's network and not about this surface, so a closed surface is put to the
  // network like any other. A pocket nothing of this size can leave may still stand in a channel the
  // accessibility probe moves freely along, and it is exactly that case the primary curve is drawn to show.
  //
  // What survives of the geometric argument is the one case the network cannot improve on: a surface
  // running away through the crystal walls a channel at this radius, so it walls one at the accessibility
  // probe's radius too, and there is nothing to ask.
  const std::vector<ComponentVerdict> referenceVerdicts =
      boundaryComponentVerdicts(accessibility, components, &against);
  for (std::size_t component = 0; component < components.numberOfComponents; ++component)
  {
    const bool percolates = components.componentPercolates[component] != 0;
    const ComponentVerdict& verdict = referenceVerdicts[component];
    const int side = percolates ? 1 : (!verdict.decided ? 0 : (verdict.accessible ? 1 : -1));
    if (side > 0)
    {
      reachableDerivative += geometry.componentDistribution[component];
    }
    else if (side < 0)
    {
      sealed += geometry.componentEnclosedVolume[component] + geometry.componentShellVolume[component];
    }
    else
    {
      unplaced += geometry.componentEnclosedVolume[component] + geometry.componentShellVolume[component];
    }
  }
}

Sampler::Swept Sampler::swept(const PoreAccessibility& accessibility, const BoundaryComponents& components,
                              const SolventExcludedGeometry& geometry) const
{
  if (!openNetwork) return {geometry.poreVolume, geometry.distribution};

  // Put to the bare network as the probe's curve puts its surfaces to the probe's: a surface running away
  // through the crystal walls a channel, and any other is asked which pore of that network it faces. What
  // faces a pore that runs nowhere at a vanishing probe is void nothing reaches at any, and is taken out.
  double sealed = 0.0;
  double unplaced = 0.0;
  double openDerivative = 0.0;
  splitByReference(accessibility, components, geometry, bare, sealed, unplaced, openDerivative);
  return {std::clamp(geometry.poreVolume - sealed - unplaced, 0.0, std::max(0.0, geometry.poreVolume)),
          openDerivative};
}

// The void volume, which the whole curve is normalised by. It is the pore volume at zero probe radius, where
// the excluded surface is the surface of the bare atoms and the sweep is the void fraction's own.
double Sampler::measureVoidVolume(double cellVolume)
{
  bare = build(0.0);
  double voidVolume = cellVolume - unionOfBallsVolume(bare, subdivisions);
  ++evaluations;

  // Of the open network, the void less what is shut in at this same vanishing probe. A remainder that is the
  // round-off of that subtraction is no void, as below for the probe's.
  if (openNetwork)
  {
    BoundaryComponents components = boundaryComponents(bare);
    std::vector<ComponentVerdict> verdicts = boundaryComponentVerdicts(bare, components);
    SolventExcludedGeometry geometry = solventExcludedGeometry(bare, 0.0, components, verdicts, subdivisions);
    double sealed = 0.0;
    double unplaced = 0.0;
    double openDerivative = 0.0;
    splitByReference(bare, components, geometry, bare, sealed, unplaced, openDerivative);
    const double open = voidVolume - sealed - unplaced;
    voidVolume = (open > 1.0e-6 * std::max(voidVolume, 1.0)) ? open : 0.0;
  }

  // The largest sphere in the void, which is where the cumulative hits zero and past which a row of the report
  // would be a full analysis spent on nothing. Taken here, where the network is already in hand, rather than
  // found by the first zero of the sweep: that zero is the same number to a bin width, and costs a hundred
  // analyses to reach.
  largestIncludedDiameter = bare.network.largestIncludedSphereDiameter();

  candidateDiameters.clear();
  candidateDiameters.reserve(bare.network.nodes.size());
  // A room is where the clearance peaks, and a node the clearance still rises from is not one: it stands in
  // a window or on a ridge, where the volume is continuous and its derivative is not, and the excess the
  // trapezium leaves across a bracket about such a point was reported as a room of that size. The passage
  // the clearance rises along is as wide at its narrowest as the node it leaves, which is how it is known.
  std::vector<char> risesFrom(bare.network.nodes.size(), 0);
  for (const VoronoiEdge& edge : bare.network.edges)
  {
    for (std::size_t end : {edge.from, edge.to})
    {
      const double peak = bare.network.nodes[end].maximalRadius;
      if (edge.from != edge.to && edge.radius >= peak - 1.0e-9 * std::max(1.0, peak)) risesFrom[end] = 1;
    }
  }
  for (std::size_t index = 0; index < bare.network.nodes.size(); ++index)
  {
    const VoronoiNode& node = bare.network.nodes[index];
    if (node.maximalRadius > 0.0 && risesFrom[index] == 0) candidateDiameters.push_back(2.0 * node.maximalRadius);
  }
  std::sort(candidateDiameters.begin(), candidateDiameters.end());
  // The family keeps its smallest diameter beside its largest. The tolerance is twice the half width of the
  // bracket about a candidate, and a family chains on for as long as each member is within it of the last: kept
  // as its largest member alone, the bracket about that started above the smaller members, past the cliffs of
  // rooms 0.008 Å narrower, whose volume then went into no spike.
  std::vector<double> merged;
  std::vector<double> lowest;
  merged.reserve(candidateDiameters.size());
  lowest.reserve(candidateDiameters.size());
  for (double diameter : candidateDiameters)
  {
    if (merged.empty() || diameter - merged.back() > peaksDiameterMergeTol)
    {
      merged.push_back(diameter);
      lowest.push_back(diameter);
    }
    else
    {
      merged.back() = std::max(merged.back(), diameter);
    }
  }
  candidateDiameters = std::move(merged);
  candidateLowestDiameters = std::move(lowest);

  scale = (voidVolume > 0.0) ? 1.0 / voidVolume : 0.0;
  return voidVolume;
}

// The volume left after blocking pockets sealed to the accessibility probe, which normalises the primary
// curve, and the sample at the floor diameter. The volume is taken at the floor radius: zero for the hybrid
// (bare void of the open network), the accessibility probe's own radius for the probe-occupiable curve.
//
// When the floor is positive, the same evaluation is a sample of both curves, and one that is wanted: the
// primary curve has a corner at the floor diameter, above which volume leaves it and below which none does,
// and a trapezium laid across such a corner would take it for a cliff. So the diameter is made a point of the
// sweep whether or not a row of the table falls on it. At a zero floor the corner is the origin itself.
Sample Sampler::measureProbe(double voidVolume)
{
  reference = build(accessibilityRadius);

  const double radius = floorRadius;
  PoreAccessibility accessibility = (std::abs(radius - accessibilityRadius) < 1.0e-15) ? reference : build(radius);
  BoundaryComponents components = boundaryComponents(accessibility);
  std::vector<ComponentVerdict> verdicts = boundaryComponentVerdicts(accessibility, components);
  SolventExcludedGeometry geometry =
      solventExcludedGeometry(accessibility, radius, components, verdicts, subdivisions);
  ++evaluations;

  const Swept all = swept(accessibility, components, geometry);

  double sealed = 0.0;
  double unplaced = 0.0;
  double reachableDerivative = 0.0;
  double reachableVolume = all.poreVolume;
  if (accessibilityRadius > 0.0)
  {
    splitByReference(accessibility, components, geometry, reference, sealed, unplaced, reachableDerivative);
    reachableVolume = std::clamp(geometry.poreVolume - sealed - unplaced, 0.0, std::max(0.0, geometry.poreVolume));
  }
  else
  {
    // No blocking probe: the primary curve is the whole of what is swept.
    reachableDerivative = all.distribution;
  }

  // The largest sphere where that curve's void is: among the pores the probe's network runs through the
  // crystal, each of whose nodes has the probe's radius less room in that network than in the bare one.
  const PoreAccessibility& decides = (accessibilityRadius > 0.0) ? reference : bare;
  largestReachableDiameter = 0.0;
  if (accessibilityRadius > 0.0 || openNetwork)
  {
    for (std::size_t node = 0; node < decides.network.nodes.size(); ++node)
    {
      if (decides.nodeAccessible[node] == 0) continue;
      largestReachableDiameter = std::max(
          largestReachableDiameter, 2.0 * (decides.network.nodes[node].maximalRadius + decides.probeRadius));
    }
  }
  else
  {
    largestReachableDiameter = largestIncludedDiameter;
  }

  // A volume this small is the round-off left by subtracting the sealed pores from a total they make up the
  // whole of, and not a pore anything could be in: a framework whose void is nothing but sealed cages ends
  // with parts in 1e10 of the cell here, and a primary curve normalised by that would be a distribution
  // made out of nothing. Below the slack the probe reaches nothing and there is no primary curve to draw.
  const double slack = 1.0e-6 * std::max(voidVolume, 1.0);
  probeAccessibleVolume = (reachableVolume > slack) ? reachableVolume : 0.0;

  Sample sample;
  sample.diameter = 2.0 * radius;
  sample.whole.cumulative = all.poreVolume * scale;
  sample.whole.distribution = 0.5 * all.distribution * scale;
  sample.reachable.cumulative = probeAccessibleVolume * reachableScale();

  // At a positive floor the derivative from above is the side of the corner the primary curve has. At a zero
  // floor this sample coincides with the origin, whose derivative is zero on the bare convex surface.
  sample.reachable.distribution = (radius > 0.0) ? 0.5 * reachableDerivative * reachableScale() : 0.0;
  return sample;
}

// What is known at a diameter, and the row of the report that goes with it. The distribution is asked for per
// unit diameter where the closed form gives the derivative in the probe radius, hence the halving.
Sample Sampler::at(double diameter, PoreSizeDistributionPoint* row)
{
  const double radius = 0.5 * diameter;
  PoreAccessibility accessibility = build(radius);
  BoundaryComponents components = boundaryComponents(accessibility);
  std::vector<ComponentVerdict> verdicts = boundaryComponentVerdicts(accessibility, components);
  SolventExcludedGeometry geometry = solventExcludedGeometry(accessibility, radius, components, verdicts, subdivisions);
  ++evaluations;

  // Below the floor the primary region is treated as a union of balls of that floor radius, so every point of
  // it has a pore size of at least the floor diameter and none of it has gone anywhere yet. The hybrid uses a
  // zero floor and asks the accessibility network at every diameter, including the wall-corrugation range.
  const Swept all = swept(accessibility, components, geometry);

  double sealed = 0.0;
  double unplaced = 0.0;
  double reachableDerivative = 0.0;
  const bool divide = radius >= floorRadius;

  // The reachable volume is what the sealed pores leave of the total, which is a difference of volumes of the
  // cell and so can end up a rounding either side of its bounds where one of the two is the whole of it.
  double reachableVolume = probeAccessibleVolume;
  if (divide)
  {
    if (accessibilityRadius > 0.0)
    {
      splitByReference(accessibility, components, geometry, reference, sealed, unplaced, reachableDerivative);
      reachableVolume =
          std::clamp(geometry.poreVolume - sealed - unplaced, 0.0, std::max(0.0, geometry.poreVolume));
    }
    else
    {
      reachableDerivative = all.distribution;
      reachableVolume = all.poreVolume;
    }
  }

  Sample sample;
  sample.diameter = diameter;
  sample.whole.cumulative = all.poreVolume * scale;
  sample.whole.distribution = 0.5 * all.distribution * scale;
  sample.reachable.cumulative = reachableVolume * reachableScale();
  sample.reachable.distribution = 0.5 * reachableDerivative * reachableScale();

  if (row != nullptr)
  {
    row->diameter = diameter;
    row->poreVolume = all.poreVolume;
    row->cumulative = sample.whole.cumulative;
    row->distribution = sample.whole.distribution;
    row->accessible = 0.5 * geometry.accessibleDistribution * scale;
    row->inaccessible = 0.5 * geometry.inaccessibleDistribution * scale;
    row->undecided = 0.5 * geometry.undecidedDistribution * scale;
    row->probeAccessiblePoreVolume = reachableVolume;
    row->probeAccessibleCumulative = sample.reachable.cumulative;
    row->probeAccessibleDistribution = sample.reachable.distribution;
    row->numberOfArcs = geometry.diagnostics.numberOfArcs;
    row->cuspedArcs = geometry.diagnostics.cuspedArcs;
    row->numberOfVertices = geometry.diagnostics.numberOfVertices;
    row->clippedVertices = geometry.diagnostics.clippedVertices;
    row->degenerateVertices = geometry.diagnostics.degenerateVertices;
  }
  return sample;
}

// The rows asked for, at the midpoint of each step. No evaluation reads anything another one writes: each
// builds its own network from the same fixed input and fills its own row, so the loop is over independent
// work and nothing but the order of the rows ties them together.
//
// Rows past the largest sphere in the void are not evaluated. The cumulative is already zero there by the
// definition of that diameter, and rebuilding the geometry for a diameter that holds no volume is the most
// expensive thing the sweep does --- every atom overlapping every other --- for an answer known in advance.
// One step past that diameter is still taken, so the cliff that puts the last of the void over is bracketed
// and the spike search has an interval to work in; everything beyond is left at the zeros the rows start as.
//
// Where that diameter comes from is the network, which is an estimate and not the definition: its nodes are
// walked up to the peaks of the clearance, and a peak the walk falls short of leaves the estimate below the
// largest sphere. A row left at zero below that sphere would be a cliff that is not there, and the search
// would corner it as a spike at the row's own diameter, with the true one beyond it never looked at. So the
// estimate is where the sweep starts to look for the end and not where it stops: while the last row evaluated
// still holds volume, the rows after it are evaluated too, a batch at a time, until one comes back empty or
// the table runs out. On a sound estimate the first pass already ends on an empty row and nothing more is done.
//
// This is where the threads go. It is nearly all of the cost of the curve, one diameter costing a whole
// geometry, and each of them lands in a row of its own with nothing summed across them --- so the rows come
// out the same on any number of threads, to the bit, which is not true of every loop that is threaded here.
std::vector<Sample> evaluateRows(Sampler& sampler, std::vector<PoreSizeDistributionPoint>& points, double step)
{
  std::vector<Sample> rows(points.size());
  for (std::size_t bin = 0; bin < points.size(); ++bin)
  {
    const double diameter = step * (static_cast<double>(bin) + 0.5);
    points[bin].diameter = diameter;
    rows[bin].diameter = diameter;
  }

  const std::size_t workers = std::max<std::size_t>(1, workersAvailable());
  auto evaluate = [&](std::size_t from, std::size_t to)
  {
    forEachIndex(to - from, workers,
                 [&](std::size_t, std::size_t i)
                 { rows[from + i] = sampler.at(rows[from + i].diameter, &points[from + i]); });
  };

  const double cutoff = sampler.largestIncludedDiameter + step;
  std::size_t evaluated = 0;
  while (evaluated < rows.size() && rows[evaluated].diameter <= cutoff) ++evaluated;
  evaluate(0, evaluated);

  while (evaluated > 0 && evaluated < rows.size() && rows[evaluated - 1].whole.cumulative > spikeFloor)
  {
    const std::size_t to = std::min(rows.size(), evaluated + workers);
    evaluate(evaluated, to);
    evaluated = to;
  }
  return rows;
}

// The samples in order of diameter: the origin, the rows, and the floor diameter put in its place among
// them. That one is a point of the sweep and not a row of the table, the rows being the evenly spaced ones the
// caller asked for; where a row falls on it already there is nothing to insert. A zero floor coincides with
// the origin and is not inserted again.
std::vector<Sample> sweepGrid(const Sample& origin, const std::vector<Sample>& rows, const Sample& probeSample,
                              double maximumDiameter)
{
  std::vector<Sample> grid;
  grid.reserve(rows.size() + 2);
  grid.push_back(origin);

  const double probeDiameter = probeSample.diameter;
  bool probePlaced = !(probeDiameter > 0.0 && probeDiameter < maximumDiameter);
  for (const Sample& row : rows)
  {
    if (!probePlaced && probeDiameter <= row.diameter)
    {
      if (probeDiameter < row.diameter) grid.push_back(probeSample);
      probePlaced = true;
    }
    grid.push_back(row);
  }
  // Past the last row, which stands half a step short of the largest diameter: the floor is a point of the sweep
  // there as anywhere, and left out no interval started at or above it and the primary curve came out empty.
  if (!probePlaced) grid.push_back(probeSample);
  return grid;
}

// Whether the trapezium across an interval is to be believed, or the interval has to be narrowed.
//
// `bothSigns` asks it of a shortfall as of an excess. A cliff leaves an excess, and so does a window: where a
// neck pinches off the derivative runs up as one over the root of the distance to it, on the one side, and a
// trapezium with an end just above such a point overstates the integral across it by many times the floor. Left
// alone that put the curve's total over one by a per cent, and it did worse to the interval's own account,
// taking the weight of any room inside the same step off against the overstatement so that the room showed no
// excess and was never looked for. The peaks-only route looks for excesses alone, having no integral to keep.
//
// Both series are asked. A cliff of the primary curve is one of the whole void too, but the two are measured
// against different volumes and the floor is the same number for both: where the probe reaches a small part of
// the void, a room that is most of that part is under the floor of the whole and was not narrowed on its account.
bool wantsNarrowing(const Sample& left, const Sample& right, double floorDiameter, bool bothSigns)
{
  auto beyond = [&](Series Sample::* of)
  {
    const double gap = unaccounted(left, right, of);
    return (bothSigns ? std::abs(gap) : gap) > spikeFloor;
  };
  return beyond(&Sample::whole) || (left.diameter >= floorDiameter && beyond(&Sample::reachable));
}

// The ends of the bracket about every room the network knows of, put among the samples of the grid.
//
// The search below finds a cliff by what the trapezium across a step fails to account for, and a step of the
// table is wide enough for that to fail the other way. Past a room the derivative is not what it was before
// it, the surface the room's going leaves exposed being reentrant surface that was not there, and a trapezium
// laid from the one value to the other across the cliff overstates the continuous part by the rise times the
// distance of the cliff from the middle of the step. In FAU that overstatement came to within three parts in a
// thousand of the weight of the room beside the twelve-ring, 9 Å³ at 6.589 Å, and the step showed no excess at
// all. The network says where the rooms are, so each is given a step of its own, a hundredth of an Angstrom
// wide, in which nothing of that size can hide it.
void bracketRooms(Sampler& sampler, std::vector<Sample>& grid)
{
  if (grid.size() < 2) return;
  const double first = grid.front().diameter;
  const double last = grid.back().diameter;

  std::vector<double> wanted;
  for (std::size_t index = 0; index < sampler.candidateDiameters.size(); ++index)
  {
    const std::optional<std::pair<double, double>> ends = sampler.candidateBracket(index);
    // a room past the end of the range is truncated with everything else there
    if (!ends.has_value() || ends->second >= last) continue;
    if (ends->first > first) wanted.push_back(ends->first);
    wanted.push_back(ends->second);
  }
  // not where the grid has a sample already
  std::erase_if(wanted, [&](double diameter)
                {
                  return std::ranges::any_of(
                      grid, [&](const Sample& sample) { return std::abs(sample.diameter - diameter) < 1.0e-9; });
                });
  if (wanted.empty()) return;

  std::vector<Sample> extra(wanted.size());
  forEachIndex(wanted.size(), workersAvailable(),
               [&](std::size_t, std::size_t i) { extra[i] = sampler.at(wanted[i], nullptr); });

  grid.insert(grid.end(), extra.begin(), extra.end());
  std::ranges::stable_sort(grid, [](const Sample& a, const Sample& b) { return a.diameter < b.diameter; });
}

// One interval of the grid, bisected towards whatever cliff it holds.
std::vector<Sample> narrowInterval(Sampler& sampler, const Sample& left, const Sample& right,
                                   std::size_t refinements, bool bothSigns)
{
  const double floorDiameter = 2.0 * sampler.floorRadius;
  std::vector<Sample> pending = {left, right};
  if (!wantsNarrowing(pending[0], pending[1], floorDiameter, bothSigns)) return pending;

  // Bisect towards it. The half that keeps an excess holds the cliff; a half that loses it was a corner of the
  // continuous part, and its own trapezium is then good enough to add and be done with.
  for (std::size_t level = 0; level < refinements; ++level)
  {
    std::vector<Sample> next;
    next.reserve(pending.size() + 2);
    bool refined = false;
    for (std::size_t k = 0; k + 1 < pending.size(); ++k)
    {
      next.push_back(pending[k]);
      if (wantsNarrowing(pending[k], pending[k + 1], floorDiameter, bothSigns))
      {
        next.push_back(sampler.at(0.5 * (pending[k].diameter + pending[k + 1].diameter), nullptr));
        refined = true;
      }
    }
    next.push_back(pending.back());
    pending = std::move(next);
    if (!refined) break;
  }
  return pending;
}

// Every interval, narrowed. Narrowing one needs nothing from outside that interval, so they are done
// independently and, where there are threads to be had, at the same time. What could not be done that way is
// the collection afterwards, which carries a running account of the spikes from one interval into the next and
// so has to see them in order; it does no geometry and costs nothing. Keeping the two apart is what leaves all
// of the expense here in work that has nothing shared in it, and each interval writing only its own entry is
// what makes the result the same however many threads did it.
std::vector<std::vector<Sample>> narrowIntervals(Sampler& sampler, const std::vector<Sample>& grid,
                                                 std::size_t refinements)
{
  std::vector<std::vector<Sample>> narrowed(grid.size() - 1);
  forEachIndex(narrowed.size(), workersAvailable(), [&](std::size_t, std::size_t i)
               { narrowed[i] = narrowInterval(sampler, grid[i], grid[i + 1], refinements, true); });
  return narrowed;
}

// The trapezia and the spikes of both series, over the narrowed intervals in order of diameter.
void collectSeries(const std::vector<std::vector<Sample>>& narrowed, double floorDiameter, double allowance,
                   Collector& whole, Collector& reachable)
{
  // No bound at all for the whole void, which is the series the other one is held to. The largest finite double
  // rather than an infinity: the build turns infinities off, and a cap no excess can reach is a cap of none.
  const double unbounded = std::numeric_limits<double>::max();

  for (const std::vector<Sample>& pending : narrowed)
  {
    for (std::size_t k = 0; k + 1 < pending.size(); ++k)
    {
      collect(whole, pending[k], pending[k + 1], &Sample::whole, unbounded);

      // (what the whole void lost beyond its derivative's account, whether or not that came to a spike of its
      // own: held to the spike, a room of a region that is a small part of the void was held to nothing)
      const double cliff = std::max(0.0, unaccounted(pending[k], pending[k + 1], &Sample::whole));

      // With a positive floor the primary curve is flat below that diameter and starts there. Since that
      // diameter is a point of the sweep, every interval lies on one side of it, and the ones below hold
      // nothing to collect. With a zero floor every interval is collected.
      //
      // What the primary curve may lose over the interval is held to what the whole void lost there, which
      // is an inequality rather than a safeguard: the remaining region is part of the void, so the volume it
      // loses between two diameters is part of the volume the void loses, and in cubic Angstrom the one is at
      // most the other. The two are normalised by different volumes, whence the ratio.
      if (pending[k].diameter >= floorDiameter)
      {
        collect(reachable, pending[k], pending[k + 1], &Sample::reachable, cliff * allowance);
      }
    }
  }
}

}  // namespace

PoreSizeDistributionCurve exactPoreSizeDistribution(const std::function<PoreAccessibility(double)>& build,
                                                    double cellVolume, double maximumDiameter, std::size_t numberOfBins,
                                                    std::size_t subdivisions, double probeRadius, double floorRadius,
                                                    std::size_t refinements, bool openNetwork)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  PoreSizeDistributionCurve curve;
  curve.cellVolume = cellVolume;
  curve.probeRadius = std::max(0.0, probeRadius);
  curve.floorRadius = std::max(0.0, floorRadius);
  requireFloorAtProbe("exactPoreSizeDistribution", curve.probeRadius, curve.floorRadius);

  const std::size_t bins = std::max<std::size_t>(1, numberOfBins);
  const double step = maximumDiameter / static_cast<double>(bins);

  Sampler sampler{build, subdivisions, curve.probeRadius, curve.floorRadius, openNetwork};

  curve.voidVolume = sampler.measureVoidVolume(cellVolume);
  const Sample probeSample = sampler.measureProbe(curve.voidVolume);
  curve.probeAccessibleVolume = sampler.probeAccessibleVolume;

  // At a vanishing probe the excluded surface is the surface of the bare atoms, which is convex everywhere it
  // is exposed, so nothing of it is reentrant and the derivative is zero. The whole void is open, so the
  // cumulative is one. Neither is evaluated; both are what the definitions say. With a positive floor the
  // primary curve is one there too, its own region being untouched until the floor diameter is passed; with a
  // zero floor it starts losing volume as soon as the wall corrugation does.
  Sample origin;
  origin.diameter = 0.0;
  // (of a void that is there: a cell that is solid throughout has no cumulative to be one, and a one set against
  // the nought of its first row was a spike of the whole weight at no diameter)
  origin.whole.cumulative = (curve.voidVolume > 0.0) ? 1.0 : 0.0;
  origin.whole.distribution = 0.0;
  origin.reachable.cumulative = (curve.probeAccessibleVolume > 0.0) ? 1.0 : 0.0;
  origin.reachable.distribution = 0.0;

  curve.points.assign(bins, PoreSizeDistributionPoint{});
  const std::vector<Sample> rows = evaluateRows(sampler, curve.points, step);
  std::vector<Sample> grid = sweepGrid(origin, rows, probeSample, maximumDiameter);
  bracketRooms(sampler, grid);
  const std::vector<std::vector<Sample>> narrowed = narrowIntervals(sampler, grid, refinements);

  // What a weight of the whole void is worth to the primary curve, the two being normalised by different
  // volumes.
  const double allowance = (curve.probeAccessibleVolume > 0.0) ? curve.voidVolume / curve.probeAccessibleVolume : 0.0;

  Collector whole;
  Collector reachable;
  collectSeries(narrowed, probeSample.diameter, allowance, whole, reachable);

  curve.integral = whole.integral;
  curve.singularWeight = whole.singularWeight;
  curve.spikes = std::move(whole.spikes);
  curve.probeAccessibleIntegral = reachable.integral;
  curve.probeAccessibleSingularWeight = reachable.singularWeight;
  curve.probeAccessibleSpikes = std::move(reachable.spikes);

  // The largest pore is the last spike where the range reached it. Where volume is left beyond the end of the
  // range it did not: the last spike found is then some smaller room, or there is none, and the largest sphere
  // is the network's, which is where the sweep would have found the cliff had it been taken that far.
  curve.truncatedWeight = std::max(0.0, grid.back().whole.cumulative);
  curve.largestDiameter = curve.spikes.empty() ? 0.0 : curve.spikes.back().diameter;
  if (curve.truncatedWeight > spikeFloor)
  {
    curve.largestDiameter = openNetwork ? sampler.largestReachableDiameter : sampler.largestIncludedDiameter;
  }
  curve.probeAccessibleTruncatedWeight = std::max(0.0, grid.back().reachable.cumulative);
  curve.probeAccessibleLargestDiameter =
      curve.probeAccessibleSpikes.empty() ? 0.0 : curve.probeAccessibleSpikes.back().diameter;
  if (curve.probeAccessibleTruncatedWeight > spikeFloor)
  {
    curve.probeAccessibleLargestDiameter = sampler.largestReachableDiameter;
  }

  curve.numberOfEvaluations = sampler.evaluations;

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - begin;
  curve.seconds = elapsed.count();
  return curve;
}

PoreSizeDistributionCurve exactPoreSizePeaks(const std::function<PoreAccessibility(double)>& build,
                                             double cellVolume, std::size_t subdivisions, double probeRadius,
                                             double floorRadius, std::size_t refinements, bool openNetwork)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  PoreSizeDistributionCurve curve;
  curve.cellVolume = cellVolume;
  curve.probeRadius = std::max(0.0, probeRadius);
  curve.floorRadius = std::max(0.0, floorRadius);
  requireFloorAtProbe("exactPoreSizePeaks", curve.probeRadius, curve.floorRadius);

  Sampler sampler{build, subdivisions, curve.probeRadius, curve.floorRadius, openNetwork};

  curve.voidVolume = sampler.measureVoidVolume(cellVolume);
  const Sample probeSample = sampler.measureProbe(curve.voidVolume);
  curve.probeAccessibleVolume = sampler.probeAccessibleVolume;

  Sample origin;
  origin.diameter = 0.0;
  // (of a void that is there: a cell that is solid throughout has no cumulative to be one, and a one set against
  // the nought of its first row was a spike of the whole weight at no diameter)
  origin.whole.cumulative = (curve.voidVolume > 0.0) ? 1.0 : 0.0;
  origin.whole.distribution = 0.0;
  origin.reachable.cumulative = (curve.probeAccessibleVolume > 0.0) ? 1.0 : 0.0;
  origin.reachable.distribution = 0.0;

  const std::vector<double>& candidates = sampler.candidateDiameters;
  std::vector<std::vector<Sample>> brackets(candidates.size());

  forEachIndex(candidates.size(), workersAvailable(),
               [&](std::size_t, std::size_t index)
               {
                 // the bracket spans the family, from below its smallest member to above its largest
                 const std::optional<std::pair<double, double>> ends = sampler.candidateBracket(index);
                 if (!ends.has_value()) return;

                 Sample left = (ends->first <= 0.0) ? origin : sampler.at(ends->first, nullptr);
                 Sample right = sampler.at(ends->second, nullptr);
                 if (refinements == 0)
                   brackets[index] = {left, right};
                 else
                   brackets[index] = narrowInterval(sampler, left, right, refinements, false);
               });

  const double allowance = (curve.probeAccessibleVolume > 0.0) ? curve.voidVolume / curve.probeAccessibleVolume : 0.0;

  // Do not reuse collectSeries here: its "sixteen brackets" merge is meant for refined intervals of ~1e-5 Å.
  // With a ±5e-3 Å peaks bracket that window is ~0.3 Å and chains every nearby MFI room into one spike.
  //
  // A bracket is one room instead, however its refinement fell: the candidates are further apart than the
  // brackets are wide, and a cliff the bisection has straddled leaves an excess on either side of itself. So what
  // a bracket's intervals hold is added up, and placed where the bisection found it.
  auto addInterval = [](PoreSizeSpike& spike, const Sample& left, const Sample& right, double weight)
  {
    const double total = spike.weight + weight;
    spike.diameter = (spike.diameter * spike.weight + 0.5 * (left.diameter + right.diameter) * weight) / total;
    spike.bracket += right.diameter - left.diameter;
    spike.weight = total;
  };

  // Left unrefined the bracket says only which side of the candidate the cliff is on, and the network's own
  // diameter is the better answer.
  auto keepSpike = [&](Collector& into, PoreSizeSpike spike, std::size_t index)
  {
    if (spike.weight <= 0.0) return;
    if (refinements == 0) spike.diameter = candidates[index];
    into.spikes.push_back(spike);
    into.singularWeight += spike.weight;
  };

  Collector whole;
  Collector reachable;
  for (std::size_t index = 0; index < brackets.size(); ++index)
  {
    const std::vector<Sample>& pending = brackets[index];
    if (pending.size() < 2) continue;

    PoreSizeSpike wholeSpike;
    PoreSizeSpike reachableSpike;
    for (std::size_t k = 0; k + 1 < pending.size(); ++k)
    {
      const Sample& left = pending[k];
      const Sample& right = pending[k + 1];
      const double wholeGap = unaccounted(left, right, &Sample::whole);
      if (wholeGap > spikeFloor) addInterval(wholeSpike, left, right, wholeGap);

      // The primary curve turns a corner at the floor diameter: flat below it, losing volume above. An interval
      // wholly below the floor holds nothing of it. One that reaches across the floor holds what is lost from
      // the floor up, and for that the floor is where it starts: the sample there is in hand, and is the one the
      // sweep makes a point of its grid for the same reason. Passing such an interval over, as one below the
      // floor, lost the spike of a room only just wider than the floor, whose bracket starts under it.
      if (right.diameter <= probeSample.diameter) continue;
      const Sample& from = (left.diameter < probeSample.diameter) ? probeSample : left;
      // (held to what the whole void lost there beyond its derivative's account, a spike of its own or not)
      const double reachableGap =
          std::min(unaccounted(from, right, &Sample::reachable), std::max(0.0, wholeGap) * allowance);
      if (reachableGap <= spikeFloor) continue;
      addInterval(reachableSpike, from, right, reachableGap);
    }
    keepSpike(whole, wholeSpike, index);
    keepSpike(reachable, reachableSpike, index);
  }

  curve.integral = 0.0;
  curve.singularWeight = whole.singularWeight;
  curve.spikes = std::move(whole.spikes);
  curve.probeAccessibleIntegral = 0.0;
  curve.probeAccessibleSingularWeight = reachable.singularWeight;
  curve.probeAccessibleSpikes = std::move(reachable.spikes);
  curve.truncatedWeight = 0.0;
  curve.largestDiameter = curve.spikes.empty() ? 0.0 : curve.spikes.back().diameter;
  curve.probeAccessibleTruncatedWeight = 0.0;
  curve.probeAccessibleLargestDiameter =
      curve.probeAccessibleSpikes.empty() ? 0.0 : curve.probeAccessibleSpikes.back().diameter;
  curve.numberOfEvaluations = sampler.evaluations;
  curve.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
  return curve;
}
