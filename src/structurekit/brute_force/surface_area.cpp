module;

module brute_force_surface_area;

import std;

import double3;
import double3x3;
import unit_cell;
import randomnumbers;
import brute_force_structure;
import brute_force_voxels;

namespace
{
// A stream of directions of its own per atom, so that the answer does not depend on how the atoms were
// shared out among threads.
double3 directionFrom(RandomNumber &random) { return random.randomVectorOnUnitSphere(); }
}  // namespace

BruteForceSurfaceArea BruteForceSurfaceArea::compute(const BruteForceStructure &structure,
                                                     const BruteForceVoxels &voxels,
                                                     std::size_t samplesPerAtom, bool keepExposedPoints)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  BruteForceSurfaceArea self;
  self.samplesPerAtom = samplesPerAtom;

  std::size_t numberOfAtoms = structure.size();
  self.areaOfAtom.assign(numberOfAtoms, 0.0);

  std::vector<double> accessibleOfAtom(numberOfAtoms, 0.0);
  std::vector<double> inaccessibleOfAtom(numberOfAtoms, 0.0);
  std::vector<double> undecidedOfAtom(numberOfAtoms, 0.0);
  std::vector<double> varianceOfAtom(numberOfAtoms, 0.0);
  std::vector<std::size_t> exposedOfAtom(numberOfAtoms, 0);

  std::vector<std::vector<double3>> pointsOfAtom(keepExposedPoints ? numberOfAtoms : 0);

  // How far off the surface to step before asking which piece of the void a point faces. The surface point
  // itself sits on the boundary of the union, where a grid labelled by its voxel centres has it inside an
  // atom as often as not; a little way out along the normal it is unambiguously in the void, and still far
  // short of anywhere it could be mistaken for a different pore.
  double stepOff = 0.5 * std::max({voxels.spacing.x, voxels.spacing.y, voxels.spacing.z});

#pragma omp parallel for schedule(dynamic, 1)
  for (std::int64_t index = 0; index < static_cast<std::int64_t>(numberOfAtoms); ++index)
  {
    std::size_t atom = static_cast<std::size_t>(index);

    double radius = structure.radii[atom];
    double sphereArea = 4.0 * std::numbers::pi * radius * radius;

    RandomNumber random{atom};

    std::size_t exposed = 0;
    std::size_t accessible = 0;
    std::size_t inaccessible = 0;
    std::size_t undecided = 0;

    for (std::size_t sample = 0; sample < samplesPerAtom; ++sample)
    {
      double3 direction = directionFrom(random);
      double3 point = structure.positions[atom] + radius * direction;

      // The point sits on its own atom's sphere, so that atom contributes exactly zero to the clearance
      // there and the minimum over all of them is zero when nothing else reaches it and negative when
      // something does. No separate occlusion test is needed, and none of the bookkeeping about which atom
      // to skip: the point is exposed when the clearance has not gone negative.
      if (structure.clearance(point) < -1.0e-9) continue;
      ++exposed;

      if (keepExposedPoints) pointsOfAtom[atom].push_back(point);

      double3 outside = point + stepOff * direction;
      std::int32_t region = voxels.regionNear(structure, outside);

      if (region < 0)
        ++undecided;
      else if (voxels.regions[static_cast<std::size_t>(region)].percolates)
        ++accessible;
      else
        ++inaccessible;
    }

    double perSample = sphereArea / static_cast<double>(samplesPerAtom);

    exposedOfAtom[atom] = exposed;
    self.areaOfAtom[atom] = static_cast<double>(exposed) * perSample;
    accessibleOfAtom[atom] = static_cast<double>(accessible) * perSample;
    inaccessibleOfAtom[atom] = static_cast<double>(inaccessible) * perSample;
    undecidedOfAtom[atom] = static_cast<double>(undecided) * perSample;

    // The count of exposed directions is binomial, so its variance is n p (1 - p) and the area's is that
    // scaled by the area each sample stands for.
    double fraction = static_cast<double>(exposed) / static_cast<double>(samplesPerAtom);
    varianceOfAtom[atom] =
        sphereArea * sphereArea * fraction * (1.0 - fraction) / static_cast<double>(samplesPerAtom);
  }

  double variance = 0.0;
  for (std::size_t atom = 0; atom < numberOfAtoms; ++atom)
  {
    self.totalArea += self.areaOfAtom[atom];
    self.accessibleArea += accessibleOfAtom[atom];
    self.inaccessibleArea += inaccessibleOfAtom[atom];
    self.undecidedArea += undecidedOfAtom[atom];
    self.numberOfExposedSamples += exposedOfAtom[atom];
    variance += varianceOfAtom[atom];
  }
  self.totalAreaError = std::sqrt(variance);

  if (keepExposedPoints)
  {
    for (const std::vector<double3> &points : pointsOfAtom)
    {
      self.exposedPoints.insert(self.exposedPoints.end(), points.begin(), points.end());
    }
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - begin;
  self.seconds = timing.count();

  return self;
}
