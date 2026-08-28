module;

module energy_shared_well_surface;

import std;

import uint3;
import double3;
import double3x3;
import unit_cell;
import crystal;
import units;
import structure_parallel;
import surface_curvature;
import energy_shared_isosurface;

// The measurement itself, in the units the field is held in. Nothing here knows what a Kelvin is: the report
// and the driver next door do that, and the reason they are a separate translation unit is written up there.

// The eight nodes around a point and how much of each of them it takes.
namespace
{
struct VoxelStencil
{
  std::int64_t i{0}, j{0}, k{0};
  double u{0.0}, v{0.0}, w{0.0};
};

VoxelStencil stencilOf(uint3 gridSize, double3 fractional)
{
  const double alongX = fractional.x * static_cast<double>(gridSize.x);
  const double alongY = fractional.y * static_cast<double>(gridSize.y);
  const double alongZ = fractional.z * static_cast<double>(gridSize.z);

  const double flooredX = std::floor(alongX);
  const double flooredY = std::floor(alongY);
  const double flooredZ = std::floor(alongZ);

  return VoxelStencil{static_cast<std::int64_t>(flooredX), static_cast<std::int64_t>(flooredY),
                      static_cast<std::int64_t>(flooredZ), alongX - flooredX, alongY - flooredY,
                      alongZ - flooredZ};
}

// The weight of corner `corner` of the stencil, the bits of the corner index selecting which side of each axis
// it sits on.
double weightOfCorner(const VoxelStencil &stencil, int corner)
{
  const double alongX = (corner & 1) ? stencil.u : 1.0 - stencil.u;
  const double alongY = (corner & 2) ? stencil.v : 1.0 - stencil.v;
  const double alongZ = (corner & 4) ? stencil.w : 1.0 - stencil.w;

  return alongX * alongY * alongZ;
}

void addCurvatureAreas(CurvatureAreas &into, const CurvatureAreas &from)
{
  into.convex += from.convex;
  into.saddle += from.saddle;
  into.concave += from.concave;
  into.flat += from.flat;
  into.unresolved += from.unresolved;

  into.numberOfConvex += from.numberOfConvex;
  into.numberOfSaddle += from.numberOfSaddle;
  into.numberOfConcave += from.numberOfConcave;
  into.numberOfFlat += from.numberOfFlat;
  into.numberOfUnresolved += from.numberOfUnresolved;
  into.numberOfDegenerate += from.numberOfDegenerate;

  into.integratedMeanCurvature += from.integratedMeanCurvature;
  into.integratedGaussianCurvature += from.integratedGaussianCurvature;
}

// What one worker gathers over its own stretch of triangles, so that nothing is written to twice and the
// reduction can be done in worker order.
struct WellAccumulator
{
  double zeroArea{0.0};
  double unmappedZeroArea{0.0};
  double area{0.0};
  double weightedArea{0.0};
  double foldedArea{0.0};
  double sharedArea{0.0};
  double sharedWeightedArea{0.0};

  // Sums of area times the quantity, divided by the area at the end.
  double walkMoment{0.0};
  double depthMoment{0.0};
  double deepest{0.0};
  double gaussianIntegral{0.0};

  std::size_t triangles{0};
  std::size_t folded{0};
  std::size_t shared{0};
  std::size_t unmapped{0};
  std::size_t rejected{0};
  std::size_t atWall{0};

  CurvatureAreas curvature;
  CurvatureAreas weightedCurvature;
};
}  // namespace

bool PeriodicFieldSampler::usable() const
{
  const std::size_t needed = static_cast<std::size_t>(this->gridSize.x) * static_cast<std::size_t>(this->gridSize.y) *
                             static_cast<std::size_t>(this->gridSize.z);

  return needed > 0 && this->field.size() >= needed;
}

double PeriodicFieldSampler::nodeValue(std::int64_t i, std::int64_t j, std::int64_t k) const
{
  const std::int64_t nx = static_cast<std::int64_t>(this->gridSize.x);
  const std::int64_t ny = static_cast<std::int64_t>(this->gridSize.y);
  const std::int64_t nz = static_cast<std::int64_t>(this->gridSize.z);

  const std::int64_t a = ((i % nx) + nx) % nx;
  const std::int64_t b = ((j % ny) + ny) % ny;
  const std::int64_t c = ((k % nz) + nz) % nz;

  return static_cast<double>(this->field[static_cast<std::size_t>((c * ny + b) * nx + a)]);
}

double3 PeriodicFieldSampler::nodeGradient(std::int64_t i, std::int64_t j, std::int64_t k) const
{
  return double3(0.5 * (this->nodeValue(i + 1, j, k) - this->nodeValue(i - 1, j, k)),
                 0.5 * (this->nodeValue(i, j + 1, k) - this->nodeValue(i, j - 1, k)),
                 0.5 * (this->nodeValue(i, j, k + 1) - this->nodeValue(i, j, k - 1)));
}

double PeriodicFieldSampler::value(double3 fractional) const
{
  if (!this->usable()) return 0.0;

  const VoxelStencil stencil = stencilOf(this->gridSize, fractional);

  double result = 0.0;
  for (int corner = 0; corner < 8; ++corner)
  {
    result += weightOfCorner(stencil, corner) * this->nodeValue(stencil.i + (corner & 1 ? 1 : 0),
                                                                stencil.j + (corner & 2 ? 1 : 0),
                                                                stencil.k + (corner & 4 ? 1 : 0));
  }

  return result;
}

double3 PeriodicFieldSampler::gradient(double3 fractional) const
{
  if (!this->usable()) return double3(0.0, 0.0, 0.0);

  const VoxelStencil stencil = stencilOf(this->gridSize, fractional);

  double3 result(0.0, 0.0, 0.0);
  for (int corner = 0; corner < 8; ++corner)
  {
    result = result + this->nodeGradient(stencil.i + (corner & 1 ? 1 : 0), stencil.j + (corner & 2 ? 1 : 0),
                                         stencil.k + (corner & 4 ? 1 : 0)) *
                          weightOfCorner(stencil, corner);
  }

  return result;
}

WellWalk walkToWell(const PeriodicFieldSampler &sampler, const UnitCell &unitCell, double3 start, double3 normal,
                    double step, double longestWalk, double isoValue)
{
  WellWalk walk;
  if (!sampler.usable() || !(step > 0.0) || !(longestWalk > 0.0)) return walk;

  const double3x3 inverseCell = unitCell.inverseCell;

  // The derivative of the field along the ray, which is the thing that has to vanish. The gradient comes back
  // held per grid step, so it is carried to a derivative with respect to position before being projected: on an
  // oblique cell, projecting the raw index-space gradient onto a Cartesian direction is not the same thing and
  // would tilt the stopping point.
  auto slopeAt = [&](double distance)
  {
    const double3 fractional = inverseCell * (start + normal * distance);
    const double3 cartesian =
        cartesianGradientOfGridGradient(inverseCell, sampler.gridSize, sampler.gradient(fractional));

    return double3::dot(cartesian, normal);
  };

  const double slopeAtWall = slopeAt(0.0);
  if (!std::isfinite(slopeAtWall)) return walk;

  // The energy is already rising on the way out, so there is nowhere to walk to and the well is the wall.
  if (slopeAtWall >= 0.0)
  {
    walk.found = true;
    walk.atWall = true;
    walk.depth = sampler.value(inverseCell * start);
    return walk;
  }

  double behind = 0.0;
  for (double distance = step; distance <= longestWalk; distance += step)
  {
    const double slope = slopeAt(distance);
    if (!std::isfinite(slope)) return walk;

    if (slope >= 0.0)
    {
      // Bracketed between `behind` and here. Bisection is enough and is what the interpolation supports: the
      // slope is continuous because the gradient is interpolated from the node differences rather than
      // differentiated out of the trilinear form, so it has a genuine root in the bracket rather than a jump
      // through zero at a voxel face.
      double low = behind;
      double high = distance;
      for (int iteration = 0; iteration < 40 && high - low > 1.0e-5; ++iteration)
      {
        const double middle = 0.5 * (low + high);
        if (slopeAt(middle) < 0.0)
        {
          low = middle;
        }
        else
        {
          high = middle;
        }
      }

      walk.distance = 0.5 * (low + high);
      walk.depth = sampler.value(inverseCell * (start + normal * walk.distance));
      walk.found = true;

      // Whose well is it. Carrying on past the minimum tells a well this wall owns from one it shares with a
      // wall facing it: solid immediately beyond means the minimum sits between the two and will be found again
      // from the other side, while a crest followed by a descent means the far wall has a well of its own.
      //
      // The threshold is a fraction of this well's own depth rather than a fixed energy, so that it scales with
      // whatever the probe and the framework happen to make the landscape. It has to be above the iso-value by
      // something: past its well a ray in a wide pore climbs towards zero and levels off there, and beyond the
      // cutoff it reaches zero exactly, so a test of merely reaching the iso-value would call every wide pore
      // narrow.
      const double threshold = isoValue + 0.05 * std::max(0.0, isoValue - walk.depth);
      for (double beyond = walk.distance + step; beyond <= longestWalk; beyond += step)
      {
        const double3 fractional = inverseCell * (start + normal * beyond);
        if (sampler.value(fractional) > threshold)
        {
          walk.shared = true;
          break;
        }
        if (slopeAt(beyond) <= 0.0) break;
      }

      return walk;
    }

    behind = distance;
  }

  // The energy fell the whole way without turning back up: no wall within reach on this normal.
  return walk;
}

WellSurface wellSurfaceOfField(const Crystal &framework, std::span<const float> field, uint3 gridSize,
                               std::span<const double3> corners, double isoValue, double thermalEnergy,
                               double longestWalk)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  WellSurface result;
  result.isoValue = isoValue;
  result.thermalEnergy = thermalEnergy;
  result.longestWalk = longestWalk;
  result.curvatureBands = curvatureBandsForGrid(framework.unitCell.cell, gridSize);

  PeriodicFieldSampler sampler{gridSize, field};
  if (!sampler.usable() || corners.size() < 3) return result;

  const double3x3 cell = framework.unitCell.cell;
  const double3x3 inverseCell = framework.unitCell.inverseCell;

  // The rays are walked in strides of one voxel. Structure finer than that is not in the field to be found, and
  // a stride shorter than a voxel only costs samples; a stride longer could step over a well and out the far
  // side of it.
  const double step = std::min({(cell * double3(1.0 / static_cast<double>(gridSize.x), 0.0, 0.0)).length(),
                                (cell * double3(0.0, 1.0 / static_cast<double>(gridSize.y), 0.0)).length(),
                                (cell * double3(0.0, 0.0, 1.0 / static_cast<double>(gridSize.z))).length()});
  result.step = step;

  // The zero surface's own triangles are still marching-cubes triangles and are held to the same size guard as
  // anywhere else. The *mapped* ones are not: an offset outward expands a convex patch by (1 + t k)^2 and there
  // is no reason for the result to fit inside a voxel, so applying the guard to them would throw away exactly
  // the convex caps this is meant to measure.
  const double largestPlausible = largestPlausibleTriangleArea(cell, gridSize);

  const double beta = (thermalEnergy > 0.0) ? 1.0 / thermalEnergy : 0.0;

  const std::size_t numberOfTriangles = corners.size() / 3;
  const std::size_t workers = workersAvailable();
  std::vector<WellAccumulator> ofWorker(workers);

  forEachBlock(numberOfTriangles, workers,
               [&](std::size_t worker, std::size_t begin, std::size_t end)
               {
                 WellAccumulator local;

                 for (std::size_t triangle = begin; triangle < end; ++triangle)
                 {
                   const std::size_t first = 3 * triangle;

                   const std::array<double3, 3> wall{cell * corners[first], cell * corners[first + 1],
                                                     cell * corners[first + 2]};

                   const double3 wallFace = double3::cross(wall[1] - wall[0], wall[2] - wall[0]);
                   const double wallArea = 0.5 * wallFace.length();
                   if (!std::isfinite(wallArea) || wallArea >= largestPlausible)
                   {
                     ++local.rejected;
                     continue;
                   }
                   local.zeroArea += wallArea;

                   std::array<double3, 3> normals{};
                   std::array<double3, 3> well{};
                   std::array<double, 3> depth{};
                   std::array<double, 3> distance{};
                   bool mapped = true;
                   std::size_t atWall = 0;
                   std::size_t shared = 0;

                   for (std::size_t vertex = 0; vertex < 3 && mapped; ++vertex)
                   {
                     normals[vertex] = cartesianNormalOfGridGradient(
                         inverseCell, gridSize, sampler.gradient(corners[first + vertex]), FieldSense::GrowsIntoSolid);

                     if (!(normals[vertex].length() > 0.0))
                     {
                       mapped = false;
                       break;
                     }

                     const WellWalk walk = walkToWell(sampler, framework.unitCell, wall[vertex], normals[vertex],
                                                      step, longestWalk, isoValue);
                     if (!walk.found)
                     {
                       mapped = false;
                       break;
                     }

                     if (walk.atWall) ++atWall;
                     if (walk.shared) ++shared;
                     distance[vertex] = walk.distance;
                     depth[vertex] = walk.depth;
                     well[vertex] = wall[vertex] + normals[vertex] * walk.distance;
                   }

                   const double3 wellFace =
                       mapped ? double3::cross(well[1] - well[0], well[2] - well[0]) : double3(0.0, 0.0, 0.0);
                   const double wellArea = 0.5 * wellFace.length();

                   if (!mapped || !std::isfinite(wellArea))
                   {
                     local.unmappedZeroArea += wallArea;
                     ++local.unmapped;
                     continue;
                   }

                   local.atWall += atWall;

                   // The weight of the triangle is the mean of the weights of its corners, which is the
                   // quadrature that goes with a mesh whose fields live at the vertices. Taking the weight of
                   // the mean depth instead would understate it, exp being convex.
                   double weight = 0.0;
                   for (std::size_t vertex = 0; vertex < 3; ++vertex)
                   {
                     const double exponent = std::min(-depth[vertex] * beta, 700.0);
                     weight += (beta > 0.0) ? std::exp(exponent) : 1.0;
                   }
                   weight /= 3.0;

                   // The patch has turned over: the offset has carried it past the focal point of one of its
                   // principal directions, through a line and out the other side.
                   const bool folded = double3::dot(wellFace, wallFace) < 0.0;
                   if (folded)
                   {
                     local.foldedArea += wellArea;
                     ++local.folded;
                   }

                   // A triangle can straddle the boundary between a wall's own well and one it shares, so it
                   // contributes the share of itself whose corners were on the shared side rather than being
                   // counted whole either way.
                   if (shared > 0)
                   {
                     const double fraction = static_cast<double>(shared) / 3.0;
                     local.sharedArea += wellArea * fraction;
                     local.sharedWeightedArea += wellArea * weight * fraction;
                     ++local.shared;
                   }

                   local.area += wellArea;
                   local.weightedArea += wellArea * weight;
                   local.walkMoment += wellArea * (distance[0] + distance[1] + distance[2]) / 3.0;
                   local.depthMoment += wellArea * (depth[0] + depth[1] + depth[2]) / 3.0;
                   local.deepest = std::min({local.deepest, depth[0], depth[1], depth[2]});
                   ++local.triangles;

                   const TriangleCurvature curvature = triangleCurvature(well, normals);

                   // The Gaussian integral wants the whole surface, folds included, because it is a topological
                   // count and counting needs the surface closed.
                   if (curvature.resolved && std::isfinite(curvature.kappa1) && std::isfinite(curvature.kappa2))
                   {
                     local.gaussianIntegral += wellArea * curvature.gaussian();
                   }

                   // The split, and the mean-curvature integral with it, want the part where shape still means
                   // something. An offset carries a curvature k to k/(1 + t k), which keeps its sign until
                   // 1 + t k turns negative --- and losing that is exactly what folding is. So a folded patch
                   // comes back with its curvatures reversed, and a concave corner the offset has overshot would
                   // be counted here as a convex cap. Past its focal point the offset is not a surface a
                   // molecule sits on. Its area is still in `area` and in `foldedArea`.
                   if (!folded)
                   {
                     local.curvature.add(wellArea, curvature, result.curvatureBands);
                     local.weightedCurvature.add(wellArea * weight, curvature, result.curvatureBands);
                   }
                 }

                 ofWorker[worker] = local;
               });

  for (const WellAccumulator &local : ofWorker)
  {
    result.zeroArea += local.zeroArea;
    result.unmappedZeroArea += local.unmappedZeroArea;
    result.area += local.area;
    result.weightedArea += local.weightedArea;
    result.foldedArea += local.foldedArea;
    result.sharedArea += local.sharedArea;
    result.sharedWeightedArea += local.sharedWeightedArea;
    result.meanWalk += local.walkMoment;
    result.meanDepth += local.depthMoment;
    result.deepestWell = std::min(result.deepestWell, local.deepest);
    result.integratedGaussianCurvature += local.gaussianIntegral;
    result.numberOfTriangles += local.triangles;
    result.numberOfFoldedTriangles += local.folded;
    result.numberOfSharedTriangles += local.shared;
    result.numberOfUnmappedTriangles += local.unmapped;
    result.numberOfRejectedTriangles += local.rejected;
    result.numberOfVerticesAtWall += local.atWall;
    addCurvatureAreas(result.curvature, local.curvature);
    addCurvatureAreas(result.weightedCurvature, local.weightedCurvature);
  }

  if (result.area > 0.0)
  {
    result.meanWalk /= result.area;
    result.meanDepth /= result.area;
  }

  if (framework.mass > 0.0)
  {
    const double perGram = Units::Angstrom * Units::Angstrom * Units::AvogadroConstant / framework.mass;
    result.gravimetricArea = result.deduplicatedArea() * perGram;
    result.gravimetricWeightedArea = result.deduplicatedWeightedArea() * perGram;
  }
  if (framework.unitCell.volume > 0.0)
  {
    result.volumetricArea = 1.0e4 * result.deduplicatedArea() / framework.unitCell.volume;
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  result.seconds = elapsed.count();

  return result;
}
