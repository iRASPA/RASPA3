module;

module exact_pore_size_distribution;

import std;

import double3;
import simulationbox;
import voronoi_accessibility;
import exact_surface_patches;
import exact_union_volume;
import exact_solvent_excluded;

namespace
{

// One diameter, and the two things known there: the fraction of the void with a pore size at least this, and
// the derivative of that fraction. They are computed by routes with nothing in common, which is what makes the
// comparison below worth making.
struct Sample
{
  double diameter{0.0};
  double cumulative{0.0};
  double distribution{0.0};
};

// The volume the derivative accounts for between two samples, by the trapezium rule, and the volume actually
// lost there. A smooth stretch has the two agreeing to the truncation error of the rule; a cliff does not.
double unaccounted(const Sample& left, const Sample& right)
{
  double lost = left.cumulative - right.cumulative;
  double accounted = 0.5 * (left.distribution + right.distribution) * (right.diameter - left.diameter);
  return lost - accounted;
}

}  // namespace

PoreSizeDistributionCurve exactPoreSizeDistribution(const std::function<VoronoiAccessibility(double)>& build,
                                                    double cellVolume, double maximumDiameter,
                                                    std::size_t numberOfBins, std::size_t subdivisions,
                                                    std::size_t refinements)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  PoreSizeDistributionCurve curve;
  curve.cellVolume = cellVolume;

  const std::size_t bins = std::max<std::size_t>(1, numberOfBins);
  const double step = maximumDiameter / static_cast<double>(bins);

  // The void volume, which the distribution is normalised by. It is the pore volume at zero probe radius,
  // where the excluded surface is the surface of the bare atoms and the sweep is the void fraction's own.
  {
    VoronoiAccessibility bare = build(0.0);
    curve.voidVolume = cellVolume - unionOfBallsVolume(bare, subdivisions);
    ++curve.numberOfEvaluations;
  }
  const double scale = (curve.voidVolume > 0.0) ? 1.0 / curve.voidVolume : 0.0;

  // What is known at a diameter, and the row of the report that goes with it. The distribution is asked for
  // per unit diameter where the closed form gives the derivative in the probe radius, hence the halving.
  auto evaluate = [&](double diameter, PoreSizeDistributionPoint* row) -> Sample
  {
    const double probeRadius = 0.5 * diameter;
    VoronoiAccessibility accessibility = build(probeRadius);
    SolventExcludedGeometry geometry = solventExcludedGeometry(accessibility, probeRadius, subdivisions);
    ++curve.numberOfEvaluations;

    Sample sample;
    sample.diameter = diameter;
    sample.cumulative = geometry.poreVolume * scale;
    sample.distribution = 0.5 * geometry.distribution * scale;

    if (row != nullptr)
    {
      row->diameter = diameter;
      row->poreVolume = geometry.poreVolume;
      row->cumulative = sample.cumulative;
      row->distribution = sample.distribution;
      row->accessible = 0.5 * geometry.accessibleDistribution * scale;
      row->inaccessible = 0.5 * geometry.inaccessibleDistribution * scale;
      row->undecided = 0.5 * geometry.undecidedDistribution * scale;
      row->numberOfArcs = geometry.numberOfArcs;
      row->cuspedArcs = geometry.cuspedArcs;
      row->numberOfVertices = geometry.numberOfVertices;
      row->clippedVertices = geometry.clippedVertices;
      row->degenerateVertices = geometry.degenerateVertices;
    }
    return sample;
  };

  // At a vanishing probe the excluded surface is the surface of the bare atoms, which is convex everywhere it
  // is exposed, so nothing of it is reentrant and the derivative is zero. The whole void is open, so the
  // cumulative is one. Neither is evaluated; both are what the definitions say.
  Sample origin;
  origin.diameter = 0.0;
  origin.cumulative = 1.0;
  origin.distribution = 0.0;

  curve.points.reserve(bins);
  std::vector<Sample> grid;
  grid.reserve(bins + 1);
  grid.push_back(origin);
  for (std::size_t bin = 0; bin < bins; ++bin)
  {
    curve.points.emplace_back();
    grid.push_back(evaluate(step * (static_cast<double>(bin) + 0.5), &curve.points.back()));
  }

  // A spike has to be worth something to be worth looking for, and below this it would be lost in the
  // trapezium error anyway.
  const double floor = 1.0e-5;
  double lastEnd = -1.0;

  for (std::size_t i = 0; i + 1 < grid.size(); ++i)
  {
    Sample left = grid[i];
    Sample right = grid[i + 1];
    double excess = unaccounted(left, right);
    if (excess <= floor)
    {
      curve.integral += 0.5 * (left.distribution + right.distribution) * (right.diameter - left.diameter);
      continue;
    }

    // Bisect towards it. The half that keeps an excess holds the cliff; a half that loses it was a corner of
    // the continuous part, and its own trapezium is then good enough to add and be done with.
    std::vector<Sample> pending = {left, right};
    for (std::size_t level = 0; level < refinements; ++level)
    {
      std::vector<Sample> next;
      next.reserve(pending.size() + 2);
      bool refined = false;
      for (std::size_t k = 0; k + 1 < pending.size(); ++k)
      {
        double gap = unaccounted(pending[k], pending[k + 1]);
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

    for (std::size_t k = 0; k + 1 < pending.size(); ++k)
    {
      double gap = unaccounted(pending[k], pending[k + 1]);
      curve.integral +=
          0.5 * (pending[k].distribution + pending[k + 1].distribution) * (pending[k + 1].diameter - pending[k].diameter);
      if (gap <= floor) continue;

      PoreSizeSpike spike;
      spike.diameter = 0.5 * (pending[k].diameter + pending[k + 1].diameter);
      spike.weight = gap;
      spike.bracket = pending[k + 1].diameter - pending[k].diameter;

      // A cliff the bisection has straddled leaves an excess on either side of itself, and those two are
      // one cliff rather than two. Anything within sixteen brackets of the last is taken to be the same
      // one, which is a stated resolution rather than a tolerance: two pore sizes that near each other are
      // the same pore size, and two families that are really distinct are further apart than this by
      // orders of magnitude.
      if (!curve.spikes.empty() &&
          (lastEnd == pending[k].diameter ||
           spike.diameter - curve.spikes.back().diameter < 16.0 * (spike.bracket + curve.spikes.back().bracket)))
      {
        PoreSizeSpike& previous = curve.spikes.back();
        double total = previous.weight + spike.weight;
        previous.diameter = (previous.diameter * previous.weight + spike.diameter * spike.weight) / total;
        previous.bracket += spike.bracket;
        previous.weight = total;
      }
      else
      {
        curve.spikes.push_back(spike);
      }
      lastEnd = pending[k + 1].diameter;
      curve.singularWeight += gap;
    }
  }

  curve.truncatedWeight = std::max(0.0, grid.back().cumulative);
  curve.largestDiameter = curve.spikes.empty() ? 0.0 : curve.spikes.back().diameter;

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - begin;
  curve.seconds = elapsed.count();
  return curve;
}
