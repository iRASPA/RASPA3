module;

module exact_pore_size_distribution;

import std;

import double3;
import simulationbox;
import pore_accessibility;
import exact_surface_patches;
import exact_union_volume;
import exact_boundary_components;
import exact_parallel;
import exact_solvent_excluded;

namespace
{

// One diameter, and the two things known there: the fraction with a pore size at least this, and the derivative
// of that fraction. They are computed by routes with nothing in common, which is what makes the comparison below
// worth making. Both are kept twice, once over the whole void and once over the part one fixed probe can reach.
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
// vanishing probe, and the part of it the fixed probe reaches, at the probe's own diameter. So the sampler is
// filled in that order --- `measureVoidVolume`, then `measureProbe`, then any number of calls to `at` --- and
// each of the first two leaves behind what the ones after it need. Past that the samples are independent of one
// another: `at` reads the sampler and writes nothing but the count of evaluations.
struct Sampler
{
  Sampler(const std::function<PoreAccessibility(double)>& build, std::size_t subdivisions, double probeRadius)
      : build(build), subdivisions(subdivisions), probeRadius(probeRadius)
  {
  }

  const std::function<PoreAccessibility(double)>& build;
  std::size_t subdivisions{1};
  double probeRadius{0.0};

  // The network of the fixed probe. It is built once and kept: every diameter above the probe's own asks the
  // same network which of its pores each surface faces.
  PoreAccessibility reference;

  // The diameter of the largest sphere that fits in the void, from the network at vanishing probe. Rows of the
  // report past that diameter hold no volume, and evaluating them would rebuild the whole geometry for nothing
  // but a signed zero: they are left at the zeros they start as.
  double largestIncludedDiameter{0.0};

  double scale{0.0};                  // one over the void volume, which the whole curve is normalised by
  double probeAccessibleVolume{0.0};  // and the volume the probe reaches, which the accessible curve is

  // Atomic because the diameters may be evaluated several at a time, and this is the one thing they all
  // write to. It is a count and not a sum, so how the increments interleave cannot change it.
  std::atomic<std::size_t> evaluations{0};

  double reachableScale() const { return (probeAccessibleVolume > 0.0) ? 1.0 / probeAccessibleVolume : 0.0; }

  double measureVoidVolume(double cellVolume);
  Sample measureProbe(double voidVolume);
  Sample at(double diameter, PoreSizeDistributionPoint* row);
};

// The void volume, which the whole curve is normalised by. It is the pore volume at zero probe radius, where
// the excluded surface is the surface of the bare atoms and the sweep is the void fraction's own.
double Sampler::measureVoidVolume(double cellVolume)
{
  PoreAccessibility bare = build(0.0);
  const double voidVolume = cellVolume - unionOfBallsVolume(bare, subdivisions);
  ++evaluations;

  // The largest sphere in the void, which is where the cumulative hits zero and past which a row of the report
  // would be a full analysis spent on nothing. Taken here, where the network is already in hand, rather than
  // found by the first zero of the sweep: that zero is the same number to a bin width, and costs a hundred
  // analyses to reach.
  largestIncludedDiameter = bare.network.largestIncludedSphereDiameter();

  scale = (voidVolume > 0.0) ? 1.0 / voidVolume : 0.0;
  return voidVolume;
}

// The volume the fixed probe can reach, which normalises the accessible curve, and the sample at the probe's
// own diameter. The volume is the pore volume there less the pores that probe cannot get into, so it is the
// void fraction's accessible volume, taken as room for the whole probe rather than for its centre.
//
// The same evaluation is a sample of both curves, and one that is wanted: the accessible curve has a corner at
// the probe's diameter, above which volume leaves it and below which none does, and a trapezium laid across
// such a corner would take it for a cliff. So the diameter is made a point of the sweep whether or not a row of
// the table falls on it.
Sample Sampler::measureProbe(double voidVolume)
{
  reference = build(probeRadius);

  BoundaryComponents components = boundaryComponents(reference);
  std::vector<ComponentVerdict> verdicts = boundaryComponentVerdicts(reference, components);
  SolventExcludedGeometry geometry =
      solventExcludedGeometry(reference, probeRadius, components, verdicts, subdivisions);
  ++evaluations;

  // A volume this small is the round-off left by subtracting the sealed pores from a total they make up the
  // whole of, and not a pore anything could be in: a framework whose void is nothing but sealed cages ends
  // with parts in 1e10 of the cell here, and an accessible curve normalised by that would be a distribution
  // made out of nothing. Below the slack the probe reaches nothing and there is no accessible curve to draw.
  const double slack = 1.0e-6 * std::max(voidVolume, 1.0);
  probeAccessibleVolume = (geometry.accessiblePoreVolume > slack) ? geometry.accessiblePoreVolume : 0.0;

  Sample sample;
  sample.diameter = 2.0 * probeRadius;
  sample.whole.cumulative = geometry.poreVolume * scale;
  sample.whole.distribution = 0.5 * geometry.distribution * scale;
  sample.reachable.cumulative = probeAccessibleVolume * reachableScale();

  // The derivative here is the one from above, which is the side of the corner the accessible curve has. At
  // this radius the two divisions are the same one, the pores of the network being the pores of the boundary.
  sample.reachable.distribution = 0.5 * geometry.accessibleDistribution * reachableScale();
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

  // The same boundary, divided by the pores of the fixed probe instead of by the pores of this diameter.
  // Above the probe's radius the two questions are different ones: a surface may close on a pore nothing of
  // this size can leave and yet stand in a channel the probe moves along, and the volume behind it is then
  // reachable although it is sealed at its own size.
  //
  // Below the probe's radius the question is not asked at all: the region that probe can occupy is a union of
  // balls of its radius, so every point of it has a pore size of at least the probe's diameter and none of it
  // has gone anywhere yet. The walk could not be taken there either, being clear of these atoms and not of the
  // larger ones the probe inflates them to.
  double sealed = 0.0;
  double unplaced = 0.0;
  double reachableDerivative = 0.0;
  const bool divide = radius >= probeRadius;
  if (divide)
  {
    const std::vector<ComponentVerdict> referenceVerdicts =
        boundaryComponentVerdicts(accessibility, components, &reference);
    for (std::size_t component = 0; component < components.numberOfComponents; ++component)
    {
      // Not `surfaceSides`, and the difference is the whole point of this curve rather than an oversight.
      // That division asks whether a surface seals void, and answers from the sign of the volume the
      // surface encloses before it consults anything: a surface that closes around void has sealed it, at
      // this radius. Here the question is whether the fixed probe can reach that void, which is a question
      // about the probe's network and not about this surface, so a closed surface is put to the network
      // like any other. A pocket nothing of this size can leave may still stand in a channel the smaller
      // probe moves freely along, and it is exactly that case the accessible curve is drawn to show.
      //
      // What survives of the geometric argument is the one case the network cannot improve on: a surface
      // running away through the crystal walls a channel at this radius, so it walls one at the probe's
      // smaller radius too, and there is nothing to ask.
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

  // The reachable volume is what the sealed pores leave of the total, which is a difference of volumes of the
  // cell and so can end up a rounding either side of its bounds where one of the two is the whole of it.
  const double reachableVolume =
      divide ? std::clamp(geometry.poreVolume - sealed - unplaced, 0.0, std::max(0.0, geometry.poreVolume))
             : probeAccessibleVolume;

  Sample sample;
  sample.diameter = diameter;
  sample.whole.cumulative = geometry.poreVolume * scale;
  sample.whole.distribution = 0.5 * geometry.distribution * scale;
  sample.reachable.cumulative = reachableVolume * reachableScale();
  sample.reachable.distribution = 0.5 * reachableDerivative * reachableScale();

  if (row != nullptr)
  {
    row->diameter = diameter;
    row->poreVolume = geometry.poreVolume;
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
// This is where the threads go. It is nearly all of the cost of the curve, one diameter costing a whole
// geometry, and each of them lands in a row of its own with nothing summed across them --- so the rows come
// out the same on any number of threads, to the bit, which is not true of every loop that is threaded here.
std::vector<Sample> evaluateRows(Sampler& sampler, std::vector<PoreSizeDistributionPoint>& points, double step)
{
  std::vector<Sample> rows(points.size());
  const double cutoff = sampler.largestIncludedDiameter + step;
  forEachIndex(points.size(), workersAvailable(),
               [&](std::size_t, std::size_t bin)
               {
                 const double diameter = step * (static_cast<double>(bin) + 0.5);
                 points[bin].diameter = diameter;
                 rows[bin].diameter = diameter;
                 if (diameter <= cutoff) rows[bin] = sampler.at(diameter, &points[bin]);
               });
  return rows;
}

// The samples in order of diameter: the origin, the rows, and the probe's own diameter put in its place among
// them. That one is a point of the sweep and not a row of the table, the rows being the evenly spaced ones the
// caller asked for; where a row falls on it already there is nothing to insert.
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
  return grid;
}

// One interval of the grid, bisected towards whatever cliff it holds.
std::vector<Sample> narrowInterval(Sampler& sampler, const Sample& left, const Sample& right,
                                   std::size_t refinements)
{
  // The refinement follows the whole void. It does not need following twice: volume that goes over a cliff in
  // a pore the probe can reach goes over the same cliff in the total, so every interval holding a spike of the
  // accessible curve holds one of this curve as well and is narrowed on its account.
  std::vector<Sample> pending = {left, right};
  if (unaccounted(pending[0], pending[1], &Sample::whole) <= spikeFloor) return pending;

  // Bisect towards it. The half that keeps an excess holds the cliff; a half that loses it was a corner of the
  // continuous part, and its own trapezium is then good enough to add and be done with.
  for (std::size_t level = 0; level < refinements; ++level)
  {
    std::vector<Sample> next;
    next.reserve(pending.size() + 2);
    bool refined = false;
    for (std::size_t k = 0; k + 1 < pending.size(); ++k)
    {
      double gap = unaccounted(pending[k], pending[k + 1], &Sample::whole);
      next.push_back(pending[k]);
      if (gap > spikeFloor)
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
               { narrowed[i] = narrowInterval(sampler, grid[i], grid[i + 1], refinements); });
  return narrowed;
}

// The trapezia and the spikes of both series, over the narrowed intervals in order of diameter.
void collectSeries(const std::vector<std::vector<Sample>>& narrowed, double probeDiameter, double allowance,
                   Collector& whole, Collector& reachable)
{
  // No bound at all for the whole void, which is the series the other one is held to. The largest finite double
  // rather than an infinity: the build turns infinities off, and a cap no excess can reach is a cap of none.
  const double unbounded = std::numeric_limits<double>::max();

  for (const std::vector<Sample>& pending : narrowed)
  {
    for (std::size_t k = 0; k + 1 < pending.size(); ++k)
    {
      const double cliff = collect(whole, pending[k], pending[k + 1], &Sample::whole, unbounded);

      // The accessible curve is flat below the probe's own diameter and starts there. Since that diameter is a
      // point of the sweep, every interval lies on one side of it, and the ones below hold nothing to collect.
      //
      // What the accessible curve may lose over the interval is held to what the whole void lost there, which
      // is an inequality rather than a safeguard: the accessible region is part of the void, so the volume it
      // loses between two diameters is part of the volume the void loses, and in cubic Angstrom the one is at
      // most the other. The two are normalised by different volumes, whence the ratio.
      if (pending[k].diameter >= probeDiameter)
      {
        collect(reachable, pending[k], pending[k + 1], &Sample::reachable, cliff * allowance);
      }
    }
  }
}

}  // namespace

PoreSizeDistributionCurve exactPoreSizeDistribution(const std::function<PoreAccessibility(double)>& build,
                                                    double cellVolume, double maximumDiameter, std::size_t numberOfBins,
                                                    std::size_t subdivisions, double probeRadius,
                                                    std::size_t refinements)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  PoreSizeDistributionCurve curve;
  curve.cellVolume = cellVolume;
  curve.probeRadius = std::max(0.0, probeRadius);

  const std::size_t bins = std::max<std::size_t>(1, numberOfBins);
  const double step = maximumDiameter / static_cast<double>(bins);

  Sampler sampler{build, subdivisions, curve.probeRadius};

  curve.voidVolume = sampler.measureVoidVolume(cellVolume);
  const Sample probeSample = sampler.measureProbe(curve.voidVolume);
  curve.probeAccessibleVolume = sampler.probeAccessibleVolume;

  // At a vanishing probe the excluded surface is the surface of the bare atoms, which is convex everywhere it
  // is exposed, so nothing of it is reentrant and the derivative is zero. The whole void is open, so the
  // cumulative is one. Neither is evaluated; both are what the definitions say. The accessible curve is one
  // there too, its own region being untouched until the probe's diameter is passed.
  Sample origin;
  origin.diameter = 0.0;
  origin.whole.cumulative = 1.0;
  origin.whole.distribution = 0.0;
  origin.reachable.cumulative = (curve.probeAccessibleVolume > 0.0) ? 1.0 : 0.0;
  origin.reachable.distribution = 0.0;

  curve.points.assign(bins, PoreSizeDistributionPoint{});
  const std::vector<Sample> rows = evaluateRows(sampler, curve.points, step);
  const std::vector<Sample> grid = sweepGrid(origin, rows, probeSample, maximumDiameter);
  const std::vector<std::vector<Sample>> narrowed = narrowIntervals(sampler, grid, refinements);

  // What a weight of the whole void is worth to the accessible curve, the two being normalised by different
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

  curve.truncatedWeight = std::max(0.0, grid.back().whole.cumulative);
  curve.largestDiameter = curve.spikes.empty() ? 0.0 : curve.spikes.back().diameter;
  curve.probeAccessibleTruncatedWeight = std::max(0.0, grid.back().reachable.cumulative);
  curve.probeAccessibleLargestDiameter =
      curve.probeAccessibleSpikes.empty() ? 0.0 : curve.probeAccessibleSpikes.back().diameter;

  curve.numberOfEvaluations = sampler.evaluations;

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - begin;
  curve.seconds = elapsed.count();
  return curve;
}
