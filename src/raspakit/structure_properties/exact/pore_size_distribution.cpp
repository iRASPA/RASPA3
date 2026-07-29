module;

module exact_pore_size_distribution;

import std;

import double3;
import simulationbox;
import voronoi_accessibility;
import exact_surface_patches;
import exact_union_volume;
import exact_boundary_components;
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
double collect(Collector& into, const Sample& left, const Sample& right, Series Sample::* of, double floor,
               double allowance)
{
  const double gap = std::min(unaccounted(left, right, of), allowance);
  into.integral += 0.5 * ((left.*of).distribution + (right.*of).distribution) * (right.diameter - left.diameter);
  if (gap <= floor) return 0.0;

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

}  // namespace

PoreSizeDistributionCurve exactPoreSizeDistribution(const std::function<VoronoiAccessibility(double)>& build,
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

  std::size_t evaluations = 0;

  // The void volume, which the whole curve is normalised by. It is the pore volume at zero probe radius, where
  // the excluded surface is the surface of the bare atoms and the sweep is the void fraction's own.
  {
    VoronoiAccessibility bare = build(0.0);
    curve.voidVolume = cellVolume - unionOfBallsVolume(bare, subdivisions);
    ++evaluations;
  }
  const double scale = (curve.voidVolume > 0.0) ? 1.0 / curve.voidVolume : 0.0;

  // The network of the fixed probe, which is what the accessible curve is divided by. It is built once and
  // kept: every diameter above the probe's own asks the same network which of its pores each surface faces.
  const VoronoiAccessibility reference = build(curve.probeRadius);

  // What is known at a diameter, and the row of the report that goes with it. The distribution is asked for
  // per unit diameter where the closed form gives the derivative in the probe radius, hence the halving.
  auto evaluate = [&](double diameter, PoreSizeDistributionPoint* row) -> Sample
  {
    const double radius = 0.5 * diameter;
    VoronoiAccessibility accessibility = build(radius);
    BoundaryComponents components = boundaryComponents(accessibility);
    std::vector<ComponentLabel> labels = labelBoundaryComponents(accessibility, components);
    SolventExcludedGeometry geometry = solventExcludedGeometry(accessibility, radius, components, labels, subdivisions);
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
    const bool divide = radius >= curve.probeRadius;
    if (divide)
    {
      const std::vector<ComponentLabel> referenceLabels =
          labelBoundaryComponents(accessibility, components, &reference);
      for (std::size_t label = 0; label < components.numberOfComponents; ++label)
      {
        // A surface running away through the crystal walls a channel at this radius, so it walls one at the
        // probe's smaller radius too, and there is nothing to ask.
        const bool percolates = components.componentPercolates[label] != 0;
        const int side =
            percolates ? 1 : (referenceLabels[label].decided ? (referenceLabels[label].accessible ? 1 : -1) : 0);
        if (side > 0)
        {
          reachableDerivative += geometry.componentDistribution[label];
        }
        else if (side < 0)
        {
          sealed += geometry.componentEnclosedVolume[label] + geometry.componentShellVolume[label];
        }
        else
        {
          unplaced += geometry.componentEnclosedVolume[label] + geometry.componentShellVolume[label];
        }
      }
    }

    // The reachable volume is what the sealed pores leave of the total, which is a difference of volumes of the
    // cell and so can end up a rounding either side of its bounds where one of the two is the whole of it.
    const double reachableVolume =
        divide ? std::clamp(geometry.poreVolume - sealed - unplaced, 0.0, std::max(0.0, geometry.poreVolume))
               : curve.probeAccessibleVolume;
    const double reachableScale = (curve.probeAccessibleVolume > 0.0) ? 1.0 / curve.probeAccessibleVolume : 0.0;

    Sample sample;
    sample.diameter = diameter;
    sample.whole.cumulative = geometry.poreVolume * scale;
    sample.whole.distribution = 0.5 * geometry.distribution * scale;
    sample.reachable.cumulative = reachableVolume * reachableScale;
    sample.reachable.distribution = 0.5 * reachableDerivative * reachableScale;

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
      row->numberOfArcs = geometry.numberOfArcs;
      row->cuspedArcs = geometry.cuspedArcs;
      row->numberOfVertices = geometry.numberOfVertices;
      row->clippedVertices = geometry.clippedVertices;
      row->degenerateVertices = geometry.degenerateVertices;
    }
    return sample;
  };

  // The volume the fixed probe can reach, which normalises the accessible curve. It is the pore volume at the
  // probe's own diameter less the pores that probe cannot get into, so it is measured at that one diameter and
  // is the void fraction's accessible volume, taken as room for the whole probe rather than for its centre.
  //
  // The same evaluation is a sample of both curves, and one that is wanted: the accessible curve has a corner
  // at the probe's diameter, above which volume leaves it and below which none does, and a trapezium laid
  // across such a corner would take it for a cliff. So the diameter is made a point of the sweep whether or
  // not a row of the table falls on it.
  const double probeDiameter = 2.0 * curve.probeRadius;
  Sample probeSample;
  {
    BoundaryComponents components = boundaryComponents(reference);
    std::vector<ComponentLabel> labels = labelBoundaryComponents(reference, components);
    SolventExcludedGeometry geometry =
        solventExcludedGeometry(reference, curve.probeRadius, components, labels, subdivisions);
    ++evaluations;

    // A volume this small is the round-off left by subtracting the sealed pores from a total they make up the
    // whole of, and not a pore anything could be in: a framework whose void is nothing but sealed cages ends
    // with parts in 1e10 of the cell here, and an accessible curve normalised by that would be a distribution
    // made out of nothing. Below the slack the probe reaches nothing and there is no accessible curve to draw.
    const double slack = 1.0e-6 * std::max(curve.voidVolume, 1.0);
    curve.probeAccessibleVolume = (geometry.accessiblePoreVolume > slack) ? geometry.accessiblePoreVolume : 0.0;

    const double reachableScale = (curve.probeAccessibleVolume > 0.0) ? 1.0 / curve.probeAccessibleVolume : 0.0;

    probeSample.diameter = probeDiameter;
    probeSample.whole.cumulative = geometry.poreVolume * scale;
    probeSample.whole.distribution = 0.5 * geometry.distribution * scale;
    probeSample.reachable.cumulative = curve.probeAccessibleVolume * reachableScale;

    // The derivative here is the one from above, which is the side of the corner the accessible curve has. At
    // this radius the two divisions are the same one, the pores of the network being the pores of the boundary.
    probeSample.reachable.distribution = 0.5 * geometry.accessibleDistribution * reachableScale;
  }

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

  // The rows asked for, at the midpoint of each step. No evaluation reads anything another one writes: each
  // builds its own network from the same fixed input and fills its own row, so the loop is over independent
  // work and nothing but the order of the rows ties them together.
  curve.points.assign(bins, PoreSizeDistributionPoint{});
  std::vector<Sample> rows(bins);
  for (std::size_t bin = 0; bin < bins; ++bin)
  {
    rows[bin] = evaluate(step * (static_cast<double>(bin) + 0.5), &curve.points[bin]);
  }

  // The samples in order of diameter, with the probe's own put in its place among them. It is a point of the
  // sweep and not a row of the table, the rows being the evenly spaced ones the caller asked for; where a row
  // falls on it already there is nothing to insert.
  std::vector<Sample> grid;
  grid.reserve(bins + 2);
  grid.push_back(origin);
  bool probePlaced = !(probeDiameter > 0.0 && probeDiameter < maximumDiameter);
  for (std::size_t bin = 0; bin < bins; ++bin)
  {
    if (!probePlaced && probeDiameter <= rows[bin].diameter)
    {
      if (probeDiameter < rows[bin].diameter) grid.push_back(probeSample);
      probePlaced = true;
    }
    grid.push_back(rows[bin]);
  }

  // A spike has to be worth something to be worth looking for, and below this it would be lost in the
  // trapezium error anyway.
  const double floor = 1.0e-5;

  // What a weight of the whole void is worth to the accessible curve, the two being normalised by different
  // volumes, and no bound at all for the whole void itself.
  const double unbounded = std::numeric_limits<double>::infinity();
  const double allowance = (curve.probeAccessibleVolume > 0.0) ? curve.voidVolume / curve.probeAccessibleVolume : 0.0;

  Collector whole;
  Collector reachable;

  // Narrowing an interval towards a cliff needs nothing from outside that interval, so the intervals are
  // narrowed one at a time and independently. What could not be done that way is the collection below, which
  // carries a running account of the spikes from one interval into the next and so has to see them in order;
  // it does no geometry and costs nothing. Keeping the two apart is what leaves all of the expense here in
  // work that has nothing shared in it.
  std::vector<std::vector<Sample>> narrowed(grid.size() - 1);
  for (std::size_t i = 0; i + 1 < grid.size(); ++i)
  {
    // The refinement follows the whole void. It does not need following twice: volume that goes over a cliff
    // in a pore the probe can reach goes over the same cliff in the total, so every interval holding a spike
    // of the accessible curve holds one of this curve as well and is narrowed on its account.
    std::vector<Sample> pending = {grid[i], grid[i + 1]};
    if (unaccounted(pending[0], pending[1], &Sample::whole) > floor)
    {
      // Bisect towards it. The half that keeps an excess holds the cliff; a half that loses it was a corner of
      // the continuous part, and its own trapezium is then good enough to add and be done with.
      for (std::size_t level = 0; level < refinements; ++level)
      {
        std::vector<Sample> next;
        next.reserve(pending.size() + 2);
        bool refined = false;
        for (std::size_t k = 0; k + 1 < pending.size(); ++k)
        {
          double gap = unaccounted(pending[k], pending[k + 1], &Sample::whole);
          next.push_back(pending[k]);
          if (gap > floor)
          {
            next.push_back(evaluate(0.5 * (pending[k].diameter + pending[k + 1].diameter), nullptr));
            refined = true;
          }
        }
        next.push_back(pending.back());
        pending = std::move(next);
        if (!refined) break;
      }
    }
    narrowed[i] = std::move(pending);
  }

  for (const std::vector<Sample>& pending : narrowed)
  {
    for (std::size_t k = 0; k + 1 < pending.size(); ++k)
    {
      const double cliff = collect(whole, pending[k], pending[k + 1], &Sample::whole, floor, unbounded);

      // The accessible curve is flat below the probe's own diameter and starts there. Since that diameter is a
      // point of the sweep, every interval lies on one side of it, and the ones below hold nothing to collect.
      //
      // What the accessible curve may lose over the interval is held to what the whole void lost there, which
      // is an inequality rather than a safeguard: the accessible region is part of the void, so the volume it
      // loses between two diameters is part of the volume the void loses, and in cubic Angstrom the one is at
      // most the other. The two are normalised by different volumes, whence the ratio.
      if (pending[k].diameter >= probeDiameter)
      {
        collect(reachable, pending[k], pending[k + 1], &Sample::reachable, floor, cliff * allowance);
      }
    }
  }

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

  curve.numberOfEvaluations = evaluations;

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - begin;
  curve.seconds = elapsed.count();
  return curve;
}
