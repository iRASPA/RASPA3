module;

module mc_backend;

import std;

import double3;
import randomnumbers;
import sampled_structure;
import sampling_backend;

namespace
{
SampledWay straightWayOf(const SampledStructure &structure, const double3 &origin, const double3 &displacement)
{
  SegmentBottleneck bottleneck = structure.segmentBottleneck(origin, displacement);

  double length = displacement.length();
  double3 direction = length > 0.0 ? (1.0 / length) * displacement : double3(0.0, 0.0, 1.0);

  return SampledWay{.radius = bottleneck.radius, .position = bottleneck.position, .direction = direction};
}

SampledWay widestWayOf(const SampledStructure &structure, const double3 &origin, const double3 &displacement,
                       std::size_t depth, RandomNumber &random)
{
  SampledWay straight = straightWayOf(structure, origin, displacement);
  if (depth == 0) return straight;

  double length = displacement.length();
  if (length <= 0.0) return straight;

  double3 lifted = straight.position;
  double liftedRadius = straight.radius;

  double step = std::max(std::abs(straight.radius), 0.1);
  const double shrink = std::pow(1.0e-3, 1.0 / static_cast<double>(samplingRefinementSteps));

  for (std::size_t attempt = 0; attempt < samplingRefinementSteps; ++attempt)
  {
    double3 across = random.randomVectorOnUnitSphere();
    across -= double3::dot(across, straight.direction) * straight.direction;

    double norm = across.length();
    if (norm <= 0.0) continue;

    double3 trial = lifted + (step / norm) * across;
    double clearance = structure.clearance(trial);

    if (clearance > liftedRadius)
    {
      lifted = trial;
      liftedRadius = clearance;
    }

    step *= shrink;
  }

  if (liftedRadius <= straight.radius) return straight;

  SampledWay first = widestWayOf(structure, origin, lifted - origin, depth - 1, random);
  SampledWay second = widestWayOf(structure, lifted, origin + displacement - lifted, depth - 1, random);

  // The narrowest point of the bent way is the narrower of the two halves' own.
  const SampledWay &bent = first.radius <= second.radius ? first : second;

  return bent.radius > straight.radius ? bent : straight;
}
}  // namespace

SamplingBackend samplingBackendCPU()
{
  SamplingBackend backend;
  backend.name = "cpu";

  backend.clearances = [](const SampledStructure &structure, std::span<const double3> positions,
                          std::span<double> into)
  {
#pragma omp parallel for schedule(static)
    for (std::int64_t index = 0; index < static_cast<std::int64_t>(positions.size()); ++index)
    {
      std::size_t i = static_cast<std::size_t>(index);
      into[i] = structure.clearance(positions[i]);
    }
  };

  backend.straightWays = [](const SampledStructure &structure, std::span<const double3> origins,
                            std::span<const double3> displacements, std::span<SampledWay> into)
  {
#pragma omp parallel for schedule(static)
    for (std::int64_t index = 0; index < static_cast<std::int64_t>(origins.size()); ++index)
    {
      std::size_t i = static_cast<std::size_t>(index);
      into[i] = straightWayOf(structure, origins[i], displacements[i]);
    }
  };

  backend.walksUphill = [](const SampledStructure &structure, std::span<const SampledPeak> starts,
                           std::span<const std::uint32_t> seeds, std::size_t steps, std::span<SampledPeak> into)
  {
    const double finalStepFraction = 1.0e-4;
    const double shrink = steps > 0 ? std::pow(finalStepFraction, 1.0 / static_cast<double>(steps)) : 1.0;

#pragma omp parallel for schedule(dynamic, 8)
    for (std::int64_t index = 0; index < static_cast<std::int64_t>(starts.size()); ++index)
    {
      std::size_t i = static_cast<std::size_t>(index);

      RandomNumber walk{static_cast<std::size_t>(seeds[i])};

      double3 best = starts[i].position;
      double bestRadius = starts[i].radius;
      double step = std::max(bestRadius, 0.1);

      for (std::size_t attempt = 0; attempt < steps; ++attempt)
      {
        double3 trial = best + step * walk.randomVectorOnUnitSphere();
        double clearance = structure.clearance(trial);

        if (clearance > bestRadius)
        {
          best = trial;
          bestRadius = clearance;
        }

        step *= shrink;
      }

      into[i] = SampledPeak{.radius = bestRadius, .position = best};
    }
  };

  backend.widestWays = [](const SampledStructure &structure, std::span<const double3> origins,
                          std::span<const double3> displacements, std::span<const std::uint32_t> seeds,
                          std::size_t depth, std::span<SampledWay> into)
  {
#pragma omp parallel for schedule(dynamic, 1)
    for (std::int64_t index = 0; index < static_cast<std::int64_t>(origins.size()); ++index)
    {
      std::size_t i = static_cast<std::size_t>(index);

      RandomNumber bending{static_cast<std::size_t>(seeds[i])};
      into[i] = widestWayOf(structure, origins[i], displacements[i], depth, bending);
    }
  };

  return backend;
}
