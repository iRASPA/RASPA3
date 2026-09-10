#include <gtest/gtest.h>

#define CL_TARGET_OPENCL_VERSION 120
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_well_surface;
import energy_shared_blocking_mask;
import energy_shared_bet_surface_area;
import energy_shared_electrostatic_potential_grid;
import energy_shared_well_field;
import energy_isosurface;
import energy_shared_isosurface;
import energy_shared_energy_backend;
import energy_well_field;
import energy_opencl_well_field;
import energy_opencl_lewiner_isosurface;
import surface_curvature;
import opencl;

// The Apollonius well surface, on landscapes whose wells can be written down.
//
// A Lennard-Jones pair bottoms out at 2^(1/6) sigma, so everything about this construction on an isolated
// atom is known in closed form: the contact distance is zero there, the energy is -epsilon, and the sheet
// is the sphere of that radius.

namespace
{
constexpr double sigma = 3.0;
constexpr double epsilon = 100.0;
constexpr double ceiling = 1.0e7;

PairInteractions oneType(double size, double strength)
{
  PairInteractions interactions;
  interactions.numberOfTypes = 1;
  interactions.names = {"X"};
  interactions.charges = {0.0};
  interactions.parameters = {PairParameters{size, strength, 0.0}};
  interactions.cutOffVDW = 20.0;
  interactions.cutOffCoulomb = 20.0;
  return interactions;
}

Crystal oneAtom(double edge)
{
  Crystal framework;
  framework.name = "one-atom";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 1.0;
  framework.atoms.push_back(CrystalAtom{double3(0.5 * edge, 0.5 * edge, 0.5 * edge), 0.0, 0});
  framework.fractionalPositions.push_back(double3(0.5, 0.5, 0.5));
  return framework;
}

LinearProbe probeOf(const PairInteractions &interactions, const std::string &name)
{
  return LinearProbe::singleSite(interactions, name).value();
}

void addAtom(Crystal &framework, double x, double y, double z)
{
  framework.atoms.push_back(CrystalAtom{double3(x, y, z), 0.0, 0});
  framework.fractionalPositions.push_back(double3(x / framework.unitCell.lengthA, y / framework.unitCell.lengthB,
                                                 z / framework.unitCell.lengthC));
}

double packingLength(const WellSurface &surface)
{
  double sum = 0.0;
  for (const FilamentVoxel &voxel : surface.filamentVoxels) sum += voxel.length;
  return sum;
}

double packingArea(const WellSurface &surface)
{
  double sum = 0.0;
  for (const FilamentVoxel &voxel : surface.filamentVoxels) sum += voxel.area;
  return sum;
}

WellField blankField(uint3 gridSize, const UnitCell &cell)
{
  WellField field;
  field.gridSize = gridSize;
  field.unitCell = cell;
  const std::size_t n =
      static_cast<std::size_t>(gridSize.x) * static_cast<std::size_t>(gridSize.y) * static_cast<std::size_t>(gridSize.z);
  field.energy.assign(n, 10.0f);
  field.distance.assign(n, 1.0f);
  field.reliability.assign(n, 1.0f);
  field.numberOfOrientations = 1;
  field.ceiling = ceiling;
  return field;
}

void markFilament(WellField &field, std::size_t i)
{
  field.energy[i] = -100.0f;
  field.distance[i] = -0.5f;
  field.reliability[i] = 0.0f;
}

Crystal emptyBox(const UnitCell &cell)
{
  Crystal framework;
  framework.name = "box";
  framework.unitCell = cell;
  framework.mass = 1.0;
  return framework;
}

std::size_t voxelOf(uint3 gridSize, double3 fractional)
{
  auto wrap = [](double x, std::size_t n)
  {
    x -= std::floor(x);
    std::size_t i = static_cast<std::size_t>(std::floor(x * static_cast<double>(n) + 1.0e-12)) % n;
    return i;
  };
  return (wrap(fractional.z, gridSize.z) * gridSize.y + wrap(fractional.y, gridSize.y)) * gridSize.x +
         wrap(fractional.x, gridSize.x);
}
}  // namespace


// The well of an isolated Lennard-Jones atom sits at 2^(1/6) sigma, one epsilon deep, and the contact
// distance is zero there.
TEST(energy_well_surface, well_field_of_an_isolated_atom)
{
  const double edge = 16.0;
  const uint3 gridSize{64, 64, 64};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = oneAtom(edge);

  WellField field =
      computeWellField(interactions, framework, probeOf(interactions, "X"), gridSize, 1, 0.0, {}, 1000.0, ceiling);
  ASSERT_EQ(field.numberOfVoxels(), 64u * 64u * 64u);

  const double rmin = wellContactPrefactor * sigma;
  const double3 centre(0.5, 0.5, 0.5);
  const double3 radial = double3(rmin / edge, 0.0, 0.0);
  const std::size_t voxel = voxelOf(gridSize, centre + radial);

  EXPECT_NEAR(static_cast<double>(field.distance[voxel]), 0.0, 0.15);
  EXPECT_NEAR(static_cast<double>(field.energy[voxel]), -epsilon, 0.15 * epsilon);
  EXPECT_GT(static_cast<double>(field.reliability[voxel]), 0.8);
}


// A blocking sphere overwrites the energy inside it with the depth ramp, so a level set can follow the
// sphere instead of stepping along the grid planes.
TEST(energy_well_surface, blocking_sphere_ramps_the_energy)
{
  const double edge = 10.0;
  const uint3 gridSize{32, 32, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = oneAtom(edge);

  BlockingSphere sphere;
  sphere.centerFractional = double3(0.5, 0.5, 0.5);
  sphere.radius = 2.0;

  WellField field = computeWellField(interactions, framework, probeOf(interactions, "X"), gridSize, 1, 0.0,
                                     std::span(&sphere, 1), 1000.0, ceiling);

  const std::size_t centre = voxelOf(gridSize, double3(0.5, 0.5, 0.5));
  EXPECT_GT(static_cast<double>(field.energy[centre]), 1000.0);
  EXPECT_LT(static_cast<double>(field.distance[centre]), 0.0);
}


TEST(energy_shared_blocking_mask, apply_to_a_flat_field)
{
  UnitCell unitCell(10.0, 10.0, 10.0);
  uint3 gridSize{8, 8, 8};
  std::vector<float> energy(8 * 8 * 8, -10.0f);

  BlockingSphere sphere;
  sphere.centerFractional = double3(0.5, 0.5, 0.5);
  sphere.radius = 2.0;
  applyBlockingSpheresToEnergy(energy, gridSize, unitCell, std::span(&sphere, 1), 1000.0, 1.0e7);

  const std::size_t centre = voxelOf(gridSize, double3(0.5, 0.5, 0.5));
  EXPECT_NEAR(static_cast<double>(energy[centre]), 2000.0, 1.0);
  EXPECT_NEAR(static_cast<double>(energy[0]), -10.0, 1.0e-6);
}


// The sheet of an isolated atom is the sphere of radius 2^(1/6) sigma. Marching cubes on a modest grid is
// a few percent off the closed form; the test asks for that, not for round-off.
TEST(energy_well_surface, isolated_atom_sheet_is_the_contact_sphere)
{
  const double edge = 16.0;
  const uint3 gridSize{96, 96, 96};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = oneAtom(edge);

  WellField field =
      computeWellField(interactions, framework, probeOf(interactions, "X"), gridSize, 1, 0.0, {}, 1000.0, ceiling);
  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.5 * epsilon, {});

  ASSERT_GT(surface.numberOfTriangles, 100u);
  const double rmin = wellContactPrefactor * sigma;
  const double expected = 4.0 * std::numbers::pi * rmin * rmin;
  EXPECT_NEAR(surface.area, expected, 0.08 * expected);
  EXPECT_NEAR(surface.meanDepth, -epsilon, 0.05 * epsilon);
  EXPECT_EQ(surface.filamentVolume, 0.0);
}


// Every edge of a closed periodic mesh appears twice, opposite ways. Quantizing the vertices the way the
// refinement welds them makes the pairing well-defined.
TEST(energy_well_surface, isolated_atom_sheet_is_watertight)
{
  const double edge = 16.0;
  const uint3 gridSize{64, 64, 64};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = oneAtom(edge);

  WellField field =
      computeWellField(interactions, framework, probeOf(interactions, "X"), gridSize, 1, 0.0, {}, 1000.0, ceiling);
  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});
  ASSERT_GT(surface.numberOfTriangles, 50u);

  std::vector<float> combined(field.numberOfVoxels());
  const float iso = effectiveTrimIsovalue(field, 0.0);
  for (std::size_t i = 0; i < field.numberOfVoxels(); ++i)
  {
    combined[i] = std::max(-field.distance[i], 0.001f * (field.energy[i] - iso));
  }
  std::vector<double3> corners = EnergyIsosurface::trianglesOfIsosurface(combined, gridSize, 0.0);

  auto quantize = [](double x) -> std::int32_t
  {
    double wrapped = x - std::floor(x);
    std::int32_t r = static_cast<std::int32_t>(std::rint(wrapped * 1048576.0));
    return r == 1048576 ? 0 : r;
  };

  struct Edge
  {
    std::int32_t ax, ay, az, bx, by, bz;
    bool operator==(const Edge &) const = default;
  };
  struct EdgeHash
  {
    std::size_t operator()(const Edge &e) const
    {
      std::size_t h = static_cast<std::size_t>(e.ax);
      h = h * 6364136223846793005ull + static_cast<std::size_t>(e.ay);
      h = h * 6364136223846793005ull + static_cast<std::size_t>(e.az);
      h = h * 6364136223846793005ull + static_cast<std::size_t>(e.bx);
      h = h * 6364136223846793005ull + static_cast<std::size_t>(e.by);
      h = h * 6364136223846793005ull + static_cast<std::size_t>(e.bz);
      return h;
    }
  };

  std::unordered_map<Edge, int, EdgeHash> counts;
  auto add = [&](double3 a, double3 b)
  {
    Edge forward{quantize(a.x), quantize(a.y), quantize(a.z), quantize(b.x), quantize(b.y), quantize(b.z)};
    Edge reverse{forward.bx, forward.by, forward.bz, forward.ax, forward.ay, forward.az};
    ++counts[forward];
    --counts[reverse];
  };

  for (std::size_t i = 0; i + 2 < corners.size(); i += 3)
  {
    add(corners[i], corners[i + 1]);
    add(corners[i + 1], corners[i + 2]);
    add(corners[i + 2], corners[i]);
  }

  std::size_t unpaired = 0;
    for (const auto &[key, count] : counts)
  {
    if (count != 0) ++unpaired;
  }
  EXPECT_EQ(unpaired, 0u);
}


// Two atoms closer than a contact diameter have d < 0 and cancelling wall directions on the mid-plane:
// the well has merged off the contact sheet onto a filament.
TEST(energy_well_surface, a_narrow_slit_has_a_filament)
{
  const double rmin = wellContactPrefactor * sigma;
  const double separation = 1.85 * rmin;
  PairInteractions interactions = oneType(sigma, epsilon);

  Crystal framework;
  framework.name = "dimer";
  framework.unitCell = UnitCell(20.0, 20.0, 20.0);
  framework.mass = 2.0;
  framework.atoms.push_back(CrystalAtom{double3(10.0 - 0.5 * separation, 10.0, 10.0), 0.0, 0});
  framework.atoms.push_back(CrystalAtom{double3(10.0 + 0.5 * separation, 10.0, 10.0), 0.0, 0});
  framework.fractionalPositions.push_back(double3((10.0 - 0.5 * separation) / 20.0, 0.5, 0.5));
  framework.fractionalPositions.push_back(double3((10.0 + 0.5 * separation) / 20.0, 0.5, 0.5));

  const uint3 gridSize{64, 64, 64};
  WellField field =
      computeWellField(interactions, framework, probeOf(interactions, "X"), gridSize, 1, 0.0, {}, 1000.0, ceiling);

  const std::size_t mid = voxelOf(gridSize, double3(0.5, 0.5, 0.5));
  EXPECT_LT(static_cast<double>(field.distance[mid]), 0.0);
  EXPECT_LT(static_cast<double>(field.reliability[mid]), wellFilamentReliabilityThreshold);
}


// A periodic pancake in the well-field (reliability, distance and energy all interior in a slab spanning
// yz): packing is the midplane, half the overlay, not a cylinder of that slab and not a 1-D file.
TEST(energy_well_surface, a_periodic_slit_packs_on_the_midplane)
{
  const double lx = 12.0, ly = 12.0, lz = 12.0;
  const uint3 gridSize{32, 32, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = emptyBox(UnitCell(lx, ly, lz));
  addAtom(framework, 0.5, 0.5, 0.5);

  WellField field = blankField(gridSize, framework.unitCell);
  field.probe = probeOf(interactions, "X");
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 14; i <= 18; ++i)
      {
        markFilament(field, field.voxelIndex(i, j, k));
      }
    }
  }

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});

  const double expected = ly * lz;
  EXPECT_GT(surface.filamentVolume, 0.0);
  EXPECT_NEAR(surface.filamentMedialArea, expected, 0.20 * expected);
  EXPECT_NEAR(packingArea(surface), surface.filamentMedialArea, 1.0e-9);
  EXPECT_LT(surface.filamentRidgeLength, 0.20 * ly);
  EXPECT_LT(packingLength(surface), 0.20 * ly);
}


// A periodic square tube along z: packing is the cell edge, not a midplane, and not A²/(4πV) of the
// voxelized overlay.
TEST(energy_well_surface, a_periodic_tube_packs_as_a_file)
{
  const double lx = 12.0, ly = 12.0, lz = 16.0;
  const uint3 gridSize{32, 32, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = emptyBox(UnitCell(lx, ly, lz));
  addAtom(framework, 0.5, 0.5, 0.5);

  WellField field = blankField(gridSize, framework.unitCell);
  field.probe = probeOf(interactions, "X");
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 14; j <= 18; ++j)
    {
      for (std::size_t i = 14; i <= 18; ++i)
      {
        markFilament(field, field.voxelIndex(i, j, k));
      }
    }
  }

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});

  EXPECT_GT(surface.filamentVolume, 0.0);
  EXPECT_NEAR(surface.filamentRidgeLength, lz, 0.20 * lz);
  EXPECT_NEAR(packingLength(surface), surface.filamentRidgeLength, 1.0e-9);
  EXPECT_LT(surface.filamentMedialArea, 0.20 * lx * ly);
  EXPECT_LT(packingArea(surface), 0.20 * lx * ly);
}


// A flattened wrapping ribbon still wraps one axis: it is a file, not a midplane of the ribbon.
TEST(energy_well_surface, a_flattened_tube_packs_as_a_file)
{
  const double lx = 12.0, ly = 12.0, lz = 16.0;
  const uint3 gridSize{32, 32, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = emptyBox(UnitCell(lx, ly, lz));
  addAtom(framework, 0.5, 0.5, 0.5);

  WellField field = blankField(gridSize, framework.unitCell);
  field.probe = probeOf(interactions, "X");
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 15; j <= 16; ++j)
    {
      for (std::size_t i = 12; i <= 20; ++i)
      {
        markFilament(field, field.voxelIndex(i, j, k));
      }
    }
  }

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});

  EXPECT_GT(surface.filamentVolume, 0.0);
  EXPECT_NEAR(surface.filamentRidgeLength, lz, 0.20 * lz);
  EXPECT_NEAR(packingLength(surface), surface.filamentRidgeLength, 1.0e-9);
  EXPECT_LT(surface.filamentMedialArea, 0.20 * lx * ly);
}


// A channel shorter than ~1.5 v_L^{1/3} holds one molecule per file per cell, not a fractional continuum.
TEST(energy_well_surface, a_short_periodic_tube_packs_one_site_per_cell)
{
  const double lx = 12.0, ly = 12.0, lz = 5.256;
  const uint3 gridSize{32, 32, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = emptyBox(UnitCell(lx, ly, lz));
  addAtom(framework, 0.5, 0.5, 0.5);

  WellField field = blankField(gridSize, framework.unitCell);
  field.probe = probeOf(interactions, "X");
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 14; j <= 18; ++j)
    {
      for (std::size_t i = 14; i <= 18; ++i)
      {
        markFilament(field, field.voxelIndex(i, j, k));
      }
    }
  }

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});

  EXPECT_GT(surface.filamentVolume, 0.0);
  EXPECT_NEAR(surface.filamentRidgeLength, nitrogenLinearSpacing, 0.20 * nitrogenLinearSpacing);
  EXPECT_NEAR(packingLength(surface), surface.filamentRidgeLength, 1.0e-9);
}


// Energy pinches a wrapping channel into disconnected pockets. The medial tube still wraps, so packing
// is one file, not the leftover midplane of the pockets.
TEST(energy_well_surface, a_pinched_tube_still_packs_as_a_file)
{
  const double lx = 12.0, ly = 12.0, lz = 16.0;
  const uint3 gridSize{32, 32, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = emptyBox(UnitCell(lx, ly, lz));
  addAtom(framework, 0.5, 0.5, 0.5);

  WellField field = blankField(gridSize, framework.unitCell);
  field.probe = probeOf(interactions, "X");
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 15; j <= 16; ++j)
    {
      for (std::size_t i = 15; i <= 16; ++i)
      {
        const std::size_t v = field.voxelIndex(i, j, k);
        field.distance[v] = -0.5f;
        field.reliability[v] = 0.0f;
        field.energy[v] = 10.0f;
      }
    }
  }
  for (std::size_t k : {4u, 5u, 6u, 20u, 21u, 22u})
  {
    for (std::size_t j = 15; j <= 16; ++j)
    {
      for (std::size_t i = 15; i <= 16; ++i)
      {
        field.energy[field.voxelIndex(i, j, k)] = -100.0f;
      }
    }
  }

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});

  EXPECT_GT(surface.filamentVolume, 0.0);
  EXPECT_NEAR(surface.filamentRidgeLength, lz, 0.20 * lz);
  EXPECT_NEAR(packingLength(surface), surface.filamentRidgeLength, 1.0e-9);
  EXPECT_LT(surface.filamentMedialArea, 0.20 * lx * ly);
}


// Two disjoint tubes along z: packing is two files of length lz, not the graph of either overlay.
TEST(energy_well_surface, two_parallel_tubes_pack_as_two_files)
{
  const double lx = 16.0, ly = 12.0, lz = 16.0;
  const uint3 gridSize{32, 32, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = emptyBox(UnitCell(lx, ly, lz));
  addAtom(framework, 0.5, 0.5, 0.5);

  WellField field = blankField(gridSize, framework.unitCell);
  field.probe = probeOf(interactions, "X");
  auto markTube = [&](int i0, int i1, int j0, int j1)
  {
    for (std::size_t k = 0; k < gridSize.z; ++k)
    {
      for (int j = j0; j <= j1; ++j)
      {
        for (int i = i0; i <= i1; ++i)
        {
          markFilament(field, field.voxelIndex(static_cast<std::size_t>(i), static_cast<std::size_t>(j), k));
        }
      }
    }
  };
  markTube(6, 10, 14, 18);
  markTube(21, 25, 14, 18);

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});

  EXPECT_GT(surface.filamentVolume, 0.0);
  EXPECT_NEAR(surface.filamentRidgeLength, 2.0 * lz, 0.20 * 2.0 * lz);
  EXPECT_NEAR(packingLength(surface), surface.filamentRidgeLength, 1.0e-9);
}


// Reliability noise punches holes in the overlay; the hop-count ridge of that swiss cheese is the whole
// blob. The energy void is still a solid 5×5 tube, and the Gelb-Gubbins inscribed radius of U = iso is the
// channel axis, so packing stays one file of length lz.
TEST(energy_well_surface, a_holey_filament_still_packs_as_a_file)
{
  const double lx = 12.0, ly = 12.0, lz = 16.0;
  const uint3 gridSize{32, 32, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = emptyBox(UnitCell(lx, ly, lz));
  addAtom(framework, 0.5, 0.5, 0.5);

  WellField field = blankField(gridSize, framework.unitCell);
  field.probe = probeOf(interactions, "X");
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 14; j <= 18; ++j)
    {
      for (std::size_t i = 14; i <= 18; ++i)
      {
        markFilament(field, field.voxelIndex(i, j, k));
        // Punch holes off the axis so the overlay DT is a plateau. The energy void stays the 5×5
        // tube, whose Gelb-Gubbins inscribed radius is the unique centre file.
        if ((i != 16 || j != 16) && (i % 2 == 0) && (j % 2 == 0))
        {
          field.reliability[field.voxelIndex(i, j, k)] = 1.0f;
        }
      }
    }
  }

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});

  EXPECT_GT(surface.filamentVolume, 0.0);
  EXPECT_NEAR(surface.filamentRidgeLength, lz, 0.20 * lz);
  EXPECT_NEAR(packingLength(surface), surface.filamentRidgeLength, 1.0e-9);
  EXPECT_LT(surface.filamentMedialArea, 0.20 * lx * ly);
}


TEST(energy_well_surface, a_blocking_sphere_removes_the_sheet_it_covers)
{
  const double edge = 16.0;
  const uint3 gridSize{64, 64, 64};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = oneAtom(edge);

  WellSurface open = wellSurfaceOfField(
      framework, interactions,
      computeWellField(interactions, framework, probeOf(interactions, "X"), gridSize, 1, 0.0, {}, 1000.0, ceiling),
      0.0, 0.001, 0.0, {});

  BlockingSphere sphere;
  sphere.centerFractional = double3(0.5, 0.5, 0.5);
  sphere.radius = edge;
  WellSurface blocked = wellSurfaceOfField(
      framework, interactions,
      computeWellField(interactions, framework, probeOf(interactions, "X"), gridSize, 1, 0.0,
                       std::span(&sphere, 1), 1000.0, ceiling),
      0.0, 0.001, 0.0, std::span(&sphere, 1));

  EXPECT_GT(open.area, 50.0);
  EXPECT_LT(blocked.area, 0.25 * open.area);
}


// Two half-strength sites on the same centre are the single site sent the long way round: the pair energies
// add back to one epsilon and the contact radius is the same sigma, so the whole field has to come back
// unchanged. This is the invariance the orientational machinery is held to.
TEST(energy_well_surface, a_zero_length_dimer_is_the_single_site)
{
  const double edge = 16.0;
  const uint3 gridSize{32, 32, 32};

  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "D"};
  interactions.charges = {0.0, 0.0};
  interactions.parameters = {PairParameters{sigma, epsilon, 0.0}, PairParameters{sigma, 0.5 * epsilon, 0.0},
                             PairParameters{sigma, 0.5 * epsilon, 0.0}, PairParameters{sigma, 0.5 * epsilon, 0.0}};
  interactions.cutOffVDW = 20.0;
  interactions.cutOffCoulomb = 20.0;

  Crystal framework = oneAtom(edge);

  LinearProbe dimer;
  dimer.name = "dimer";
  dimer.headTailSymmetric = true;
  dimer.sites.push_back(LinearProbe::Site{1, 0.0, 0.0, "D"});
  dimer.sites.push_back(LinearProbe::Site{1, 0.0, 0.0, "D"});

  WellField one =
      computeWellField(interactions, framework, probeOf(interactions, "X"), gridSize, 1, 0.0, {}, 1000.0, ceiling);
  WellField two = computeWellField(interactions, framework, dimer, gridSize, 64, 0.5 * epsilon, {}, 1000.0, ceiling);

  // No length, so one direction stands for the sphere and the free energy is the energy itself.
  ASSERT_EQ(two.numberOfOrientations, 1u);
  ASSERT_EQ(one.numberOfVoxels(), two.numberOfVoxels());

  double worstEnergy = 0.0;
  double worstDistance = 0.0;
  double worstReliability = 0.0;
  for (std::size_t i = 0; i < one.numberOfVoxels(); ++i)
  {
    double a = static_cast<double>(one.energy[i]);
    double b = static_cast<double>(two.energy[i]);
    worstEnergy = std::max(worstEnergy, std::abs(a - b) / std::max(1.0, std::abs(a)));
    worstDistance = std::max(worstDistance, std::abs(static_cast<double>(one.distance[i] - two.distance[i])));
    worstReliability =
        std::max(worstReliability, std::abs(static_cast<double>(one.reliability[i] - two.reliability[i])));
  }
  EXPECT_LT(worstEnergy, 1.0e-4);
  EXPECT_LT(worstDistance, 1.0e-6);
  EXPECT_LT(worstReliability, 1.0e-6);
}


// A dimer against an isolated atom fits best broadside on: both sites then sit at sqrt(R² + l²/4) from the
// atom, so contact falls where that reaches 2^(1/6) sigma --- for Lennard-Jones the contact radius and the
// well radius are one number --- and the well floor at kT = 0 is two epsilon deep on the sphere of radius
// sqrt(rmin² - l²/4). The refined sheet has to be that smaller sphere, twice as deep as the single site's.
TEST(energy_well_surface, a_dimer_touches_broadside_on)
{
  const double edge = 16.0;
  const uint3 gridSize{64, 64, 64};
  const double bond = 2.4;
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = oneAtom(edge);

  LinearProbe dimer;
  dimer.name = "dimer";
  dimer.headTailSymmetric = true;
  dimer.sites.push_back(LinearProbe::Site{0, -0.5 * bond, 0.0, "X"});
  dimer.sites.push_back(LinearProbe::Site{0, 0.5 * bond, 0.0, "X"});

  WellField field = computeWellField(interactions, framework, dimer, gridSize, 128, 0.0, {}, 1000.0, ceiling);
  ASSERT_EQ(field.numberOfOrientations, 128u);

  // Both sites in their own minimum at once, which no single orientation of a stiffer landscape can beat.
  EXPECT_NEAR(field.deepestEnergy(), -2.0 * epsilon, 0.1 * epsilon);

  const double rmin = wellContactPrefactor * sigma;
  const double radius = std::sqrt(rmin * rmin - 0.25 * bond * bond);

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});
  ASSERT_GT(surface.numberOfTriangles, 100u);
  const double expected = 4.0 * std::numbers::pi * radius * radius;
  EXPECT_NEAR(surface.area, expected, 0.10 * expected);
  EXPECT_NEAR(surface.meanDepth, -2.0 * epsilon, 0.1 * epsilon);
  EXPECT_EQ(surface.filamentVolume, 0.0);
}


// A charged probe handed the framework's potential feels it, and by how much is checkable without the
// engine's Ewald builder: hand the field a potential whose two halves are known exactly. The far half is a
// cosine wave, which trilinear interpolation reproduces exactly on a grid node; the split is put at an
// alpha so small that the erfc near half is the bare Coulomb sum, and the cell is large enough that only
// the neutral +1/-1 pair is inside the cutoff. The shift between the charged and uncharged fields is then
// q phi(x) + f q (1/r+ - 1/r-) to the last float digit, and it only comes out right if the smooth lookup,
// the pair walk, and the conversion factor all agree. The probe's charge rides on a site with no dispersion
// of its own --- the massless centre site of the N2 model --- so this also checks that such a site is kept
// when there is a potential to act on, and the contact distance must not move at all: charges change how
// deep the well is, not where the wall stands.
TEST(energy_well_surface, a_charged_probe_feels_the_framework)
{
  const double edge = 40.0;
  const uint3 gridSize{40, 40, 40};
  const double conversionFactor = 167101.0;  // any nonzero factor serves; the field takes it as an argument
  const double waveHeight = 250.0;
  const double probeCharge = 0.5;

  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "Q"};
  interactions.charges = {0.0, probeCharge};
  interactions.parameters = {PairParameters{sigma, epsilon, 0.0}, PairParameters{0.0, 0.0, 0.0},
                             PairParameters{0.0, 0.0, 0.0}, PairParameters{0.0, 0.0, 0.0}};
  interactions.cutOffVDW = 12.0;
  interactions.cutOffCoulomb = 12.0;

  Crystal framework;
  framework.name = "dipole";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 2.0;
  framework.atoms.push_back(CrystalAtom{double3(21.5, 20.0, 20.0), 1.0, 0});
  framework.atoms.push_back(CrystalAtom{double3(18.5, 20.0, 20.0), -1.0, 0});
  framework.fractionalPositions.push_back(double3(21.5 / edge, 0.5, 0.5));
  framework.fractionalPositions.push_back(double3(18.5 / edge, 0.5, 0.5));

  // A dispersion site and, on the same spot, a bare charge site: the neutral N2 pattern.
  LinearProbe probe;
  probe.name = "Q";
  probe.headTailSymmetric = true;
  probe.sites.push_back(LinearProbe::Site{0, 0.0, 0.0, "X"});
  probe.sites.push_back(LinearProbe::Site{1, 0.0, probeCharge, "Q"});
  ASSERT_TRUE(probe.isCharged());

  ElectrostaticPotentialGrid potential;
  potential.gridSize = gridSize;
  potential.unitCell = framework.unitCell;
  potential.alpha = 1.0e-6;
  potential.cutOff = interactions.cutOffCoulomb;
  potential.smoothPotential.resize(40u * 40u * 40u);
  for (std::size_t k = 0; k < 40; ++k)
  {
    for (std::size_t j = 0; j < 40; ++j)
    {
      for (std::size_t i = 0; i < 40; ++i)
      {
        double wave = waveHeight * std::cos(2.0 * std::numbers::pi * static_cast<double>(i) / 40.0);
        potential.smoothPotential[(k * 40 + j) * 40 + i] = static_cast<float>(wave);
      }
    }
  }

  WellField without = computeWellField(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling);
  WellField with = computeWellField(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling,
                                    &potential, conversionFactor);

  EXPECT_TRUE(without.chargesIgnored);
  EXPECT_TRUE(with.chargesIncluded);

  // The probe voxel at (26, 20, 20): 4.5 Å from the positive atom, 7.5 Å from the negative one, and every
  // periodic image outside the cutoff.
  const std::size_t voxel = voxelOf(gridSize, double3(26.0 / edge, 0.5, 0.5));
  const double smooth = probeCharge * waveHeight * std::cos(2.0 * std::numbers::pi * 0.65);
  const double near = conversionFactor * probeCharge * (1.0 / 4.5 - 1.0 / 7.5);
  const double shift = static_cast<double>(with.energy[voxel]) - static_cast<double>(without.energy[voxel]);
  EXPECT_NEAR(shift, smooth + near, 1.0e-3 * std::abs(smooth + near));

  for (std::size_t i = 0; i < with.numberOfVoxels(); ++i)
  {
    ASSERT_EQ(with.distance[i], without.distance[i]);
  }

  // The same sum at a real split, where the screening is doing something. The near half is read off a table
  // rather than from the library's error function --- it is the innermost thing the walk does and the exact
  // one costs more than everything around it --- so this holds the table to the function it stands for, at
  // the two separations above, through the field itself.
  ElectrostaticPotentialGrid screened = potential;
  screened.alpha = 0.3;
  std::ranges::fill(screened.smoothPotential, 0.0f);

  WellField damped = computeWellField(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling,
                                      &screened, conversionFactor);
  const double exactly = conversionFactor * probeCharge *
                         (std::erfc(0.3 * 4.5) / 4.5 - std::erfc(0.3 * 7.5) / 7.5);
  const double tabulated =
      static_cast<double>(damped.energy[voxel]) - static_cast<double>(without.energy[voxel]);
  EXPECT_NEAR(tabulated, exactly, 1.0e-5 * std::abs(exactly));
}


// TraPPE N2 puts the compensating charge on a massless site with no Lennard-Jones of its own. Place a
// dispersing N just inside σ of an opposite framework charge: Coulomb is gated (otherwise the dummy's
// 1/r through the atom is a well of tens of thousands of kelvin). The contact-distance field must not
// move; only electrostatics are refused. The dummy's own 1 Å floor does not fire here.
TEST(energy_well_surface, a_bare_charge_cannot_see_through_the_core)
{
  const double edge = 40.0;
  const uint3 gridSize{40, 40, 40};
  const double conversionFactor = 167101.0;
  const double coreSigma = 3.10;
  const double bond = 0.55;

  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "Q"};
  interactions.charges = {0.0, 1.0};
  interactions.parameters = {PairParameters{coreSigma, epsilon, 0.0}, PairParameters{0.0, 0.0, 0.0},
                             PairParameters{0.0, 0.0, 0.0}, PairParameters{0.0, 0.0, 0.0}};
  interactions.cutOffVDW = 12.0;
  interactions.cutOffCoulomb = 12.0;

  Crystal framework;
  framework.name = "core";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 1.0;
  framework.atoms.push_back(CrystalAtom{double3(20.0, 20.0, 20.0), -1.0, 0});
  framework.fractionalPositions.push_back(double3(0.5, 0.5, 0.5));

  LinearProbe probe;
  probe.name = "N2";
  probe.headTailSymmetric = true;
  probe.sites.push_back(LinearProbe::Site{0, -bond, 0.0, "X"});
  probe.sites.push_back(LinearProbe::Site{1, 0.0, 1.0, "Q"});
  probe.sites.push_back(LinearProbe::Site{0, bond, 0.0, "X"});

  ElectrostaticPotentialGrid potential;
  potential.gridSize = gridSize;
  potential.unitCell = framework.unitCell;
  potential.alpha = 1.0e-6;
  potential.cutOff = interactions.cutOffCoulomb;
  potential.smoothPotential.assign(40u * 40u * 40u, 1.0e5f);

  WellField without = computeWellField(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling);
  WellField with = computeWellField(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling,
                                    &potential, conversionFactor);

  const std::size_t voxel = voxelOf(gridSize, double3(23.0 / edge, 0.5, 0.5));
  const double uncharged = static_cast<double>(without.energy[voxel]);
  const double charged = static_cast<double>(with.energy[voxel]);
  const double leaked = conversionFactor / 3.0;

  EXPECT_GT(uncharged, -2.0 * epsilon);
  EXPECT_NEAR(charged, uncharged, 1.0e-2 * std::max(1.0, std::abs(uncharged)));
  EXPECT_GT(charged, uncharged - 0.05 * leaked);
  EXPECT_EQ(with.distance[voxel], without.distance[voxel]);
}


// Soft TraPPE: rmin is the well, Coulomb is gated at r < σ. One N inside σ while N_com is not:
// skip electrostatics (otherwise q_COM · q_O / r is a fake monopole) but keep the repulsive LJ —
// do not ceiling the pose. Between σ and rmin the inner well shoulder is legal, quadrupole on.
TEST(energy_well_surface, an_overlapping_n2_orientation_is_not_a_coulomb_well)
{
  const double edge = 40.0;
  const uint3 gridSize{40, 40, 40};
  const double conversionFactor = 167101.0;
  const double nSigma = 3.30;
  const double bond = 0.55;

  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "Q"};
  interactions.charges = {0.0, 0.810};
  interactions.parameters = {PairParameters{nSigma, epsilon, 0.0}, PairParameters{0.0, 0.0, 0.0},
                             PairParameters{0.0, 0.0, 0.0}, PairParameters{0.0, 0.0, 0.0}};
  interactions.cutOffVDW = 12.0;
  interactions.cutOffCoulomb = 12.0;

  Crystal framework;
  framework.name = "oxygen";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 1.0;
  framework.atoms.push_back(CrystalAtom{double3(20.0, 20.0, 20.0), -1.025, 0});
  framework.fractionalPositions.push_back(double3(0.5, 0.5, 0.5));

  LinearProbe probe;
  probe.name = "N2";
  probe.headTailSymmetric = true;
  probe.sites.push_back(LinearProbe::Site{0, -bond, -0.405, "X"});
  probe.sites.push_back(LinearProbe::Site{1, 0.0, 0.810, "Q"});
  probe.sites.push_back(LinearProbe::Site{0, bond, -0.405, "X"});

  ElectrostaticPotentialGrid potential;
  potential.gridSize = gridSize;
  potential.unitCell = framework.unitCell;
  potential.alpha = 1.0e-6;
  potential.cutOff = interactions.cutOffCoulomb;
  potential.smoothPotential.assign(40u * 40u * 40u, 0.0f);

  WellField field = computeWellField(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling,
                                     &potential, conversionFactor);

  // Default axis is z. COM 3 Å along z: closer N at 2.45 Å < σ, N_com at 3 Å >> 1 Å.
  // Coulomb is skipped; repulsive LJ is kept (not ceilinged, not a monopole well).
  const std::size_t overlapping = voxelOf(gridSize, double3(0.5, 0.5, 23.0 / edge));
  const double overlapEnergy = static_cast<double>(field.energy[overlapping]);
  const double leaked = conversionFactor * 0.810 * 1.025 / 3.0;
  EXPECT_GT(overlapEnergy, 0.0);
  EXPECT_LT(overlapEnergy, ceiling * 0.5);
  EXPECT_GT(overlapEnergy, -0.05 * leaked);

  // σ < r < rmin: inner LJ well, full quadrupole. GCMC samples this; it is not a wall.
  const std::size_t shoulder = voxelOf(gridSize, double3(0.5, 0.5, 24.0 / edge));
  const double shoulderEnergy = static_cast<double>(field.energy[shoulder]);
  EXPECT_LT(shoulderEnergy, ceiling * 0.5);
  EXPECT_GT(shoulderEnergy, -5000.0);

  const std::size_t clear = voxelOf(gridSize, double3(0.5, 0.5, 26.0 / edge));
  const double clearEnergy = static_cast<double>(field.energy[clear]);
  EXPECT_LT(clearEnergy, ceiling * 0.5);
  EXPECT_GT(clearEnergy, -5000.0);
}


// T-shaped contact: both N just outside rmin, COM closer to the oxygen than either N (~3.70 Å).
// VDW does not reject. Full quadrupole electrostatics must stay on — not the old rmin core on
// N_com, which left a −70 kT monopole. A 1 Å dummy floor does not fire here.
TEST(energy_well_surface, a_t_shaped_n2_contact_is_a_quadrupole_not_a_monopole)
{
  const double edge = 40.0;
  const uint3 gridSize{40, 400, 40};
  const double conversionFactor = 167101.0;
  const double nSigma = 3.30;
  const double bond = 0.55;

  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "Q"};
  interactions.charges = {0.0, 0.810};
  interactions.parameters = {PairParameters{nSigma, epsilon, 0.0}, PairParameters{0.0, 0.0, 0.0},
                             PairParameters{0.0, 0.0, 0.0}, PairParameters{0.0, 0.0, 0.0}};
  interactions.cutOffVDW = 12.0;
  interactions.cutOffCoulomb = 12.0;

  Crystal framework;
  framework.name = "oxygen";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 1.0;
  framework.atoms.push_back(CrystalAtom{double3(20.0, 20.0, 20.0), -1.025, 0});
  framework.fractionalPositions.push_back(double3(0.5, 0.5, 0.5));

  LinearProbe probe;
  probe.name = "N2";
  probe.headTailSymmetric = true;
  probe.sites.push_back(LinearProbe::Site{0, -bond, -0.405, "X"});
  probe.sites.push_back(LinearProbe::Site{1, 0.0, 0.810, "Q"});
  probe.sites.push_back(LinearProbe::Site{0, bond, -0.405, "X"});

  ElectrostaticPotentialGrid potential;
  potential.gridSize = gridSize;
  potential.unitCell = framework.unitCell;
  potential.alpha = 1.0e-6;
  potential.cutOff = interactions.cutOffCoulomb;
  potential.smoothPotential.assign(static_cast<std::size_t>(gridSize.x) * gridSize.y * gridSize.z, 0.0f);

  WellField field = computeWellField(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling,
                                     &potential, conversionFactor);

  // Default axis is z. COM 3.7 Å along y: N at ~3.74 Å > rmin ≈ 3.70 Å, COM at 3.7 Å >> 1 Å.
  const std::size_t voxel = voxelOf(gridSize, double3(0.5, 23.7 / edge, 0.5));
  const double energy = static_cast<double>(field.energy[voxel]);
  EXPECT_GT(energy, -3000.0);
  EXPECT_LT(energy, 500.0);
  EXPECT_LT(energy, ceiling * 0.5);
}


// The double-float helpers live in the OpenCL prelude and are compiled without the host's -ffast-math,
// so the sum is run on the device: 1 + 1e-8, a hundred thousand times. A plain float stays at 1.
TEST(energy_well_surface, double_float_accumulates_an_alternating_sum)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value())
  {
    GTEST_SKIP() << "no OpenCL device";
  }

  std::string source = std::string(WellFieldOpenCL::doubleFloatSource) + R"foo(
__kernel void Accumulate(__global float *out)
{
  df a = (df)(0.0f, 0.0f);
  a = df_add_f(a, 1.0f);
  for(int i = 0; i < 100000; i++) a = df_add_f(a, 1.0e-8f);
  out[0] = df_hi(a);
}
)foo";
  const char *ptr = source.c_str();
  cl_int err = CL_SUCCESS;
  cl_program program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &ptr, nullptr, &err);
  ASSERT_EQ(err, CL_SUCCESS);
  err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    clReleaseProgram(program);
    FAIL() << "double-float prelude failed to build";
  }
  cl_kernel kernel = clCreateKernel(program, "Accumulate", &err);
  ASSERT_EQ(err, CL_SUCCESS);
  cl_mem out = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_WRITE_ONLY, sizeof(cl_float), nullptr, &err);
  ASSERT_EQ(err, CL_SUCCESS);
  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &out);
  ASSERT_EQ(err, CL_SUCCESS);
  std::size_t one = 1;
  err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 1, nullptr, &one, nullptr, 0, nullptr, nullptr);
  ASSERT_EQ(err, CL_SUCCESS);
  cl_float recovered = 0.0f;
  err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), out, CL_TRUE, 0, sizeof(cl_float), &recovered, 0, nullptr,
                            nullptr);
  ASSERT_EQ(err, CL_SUCCESS);
  clReleaseMemObject(out);
  clReleaseKernel(kernel);
  clReleaseProgram(program);

  EXPECT_NEAR(static_cast<double>(recovered), 1.001, 2.0e-6);
}


TEST(energy_well_surface, gpu_well_field_agrees_with_cpu)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value())
  {
    GTEST_SKIP() << "no OpenCL device";
  }

  const double edge = 16.0;
  const uint3 gridSize{32, 32, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = oneAtom(edge);
  LinearProbe probe = probeOf(interactions, "X");

  WellField cpu = WellFieldCPU::compute(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling);
  WellField gpu = WellFieldOpenCL::compute(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling);

  ASSERT_EQ(cpu.numberOfVoxels(), gpu.numberOfVoxels());
  double worstEnergy = 0.0;
  double worstDistance = 0.0;
  double worstReliability = 0.0;
  for (std::size_t i = 0; i < cpu.numberOfVoxels(); ++i)
  {
    double a = static_cast<double>(cpu.energy[i]);
    double b = static_cast<double>(gpu.energy[i]);
    worstEnergy = std::max(worstEnergy, std::abs(a - b) / std::max(1.0, std::abs(a)));
    worstDistance = std::max(worstDistance, std::abs(static_cast<double>(cpu.distance[i] - gpu.distance[i])));
    worstReliability =
        std::max(worstReliability, std::abs(static_cast<double>(cpu.reliability[i] - gpu.reliability[i])));
  }
  EXPECT_LT(worstEnergy, 2.0e-3);
  EXPECT_LT(worstDistance, 0.05);
  EXPECT_LT(worstReliability, 0.05);

  WellSurface cpuSurface = wellSurfaceOfField(framework, interactions, cpu, 0.0, 0.001, 0.0, {});
  EnergyBackend gpuBackend;
  gpuBackend.isosurfaceTriangles = [](std::span<const float> field, uint3 size, double isoValue)
  { return trianglesOfLewinerIsosurface(field, size, isoValue); };
  gpuBackend.refineWellVertices = [](std::vector<double3> &corners, std::vector<double> &refined,
                                     const PairInteractions &pairs, const Crystal &crystal,
                                     const NeighbourhoodParameters &parameters,
                                     std::span<const BlockingSphere> spheres, double trim)
  { WellFieldOpenCL::refineVertices(corners, refined, pairs, crystal, parameters, spheres, trim); };
  WellSurface gpuSurface = wellSurfaceOfField(framework, interactions, gpu, 0.0, 0.001, 0.0, {}, &gpuBackend);
  ASSERT_GT(cpuSurface.area, 10.0);
  EXPECT_NEAR(gpuSurface.area, cpuSurface.area, 0.02 * cpuSurface.area);
}

TEST(energy_well_surface, gpu_lewiner_agrees_with_cpu)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value())
  {
    GTEST_SKIP() << "no OpenCL device";
  }

  // A sphere plus a weaker offset blob so some cubes are the ambiguous Lewiner cases,
  // not only the single-cut Bourke ones.
  const uint3 gridSize{32, 32, 32};
  const std::size_t n = static_cast<std::size_t>(gridSize.x) * gridSize.y * gridSize.z;
  std::vector<float> field(n);
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < gridSize.x; ++i)
      {
        auto at = [&](double cx, double cy, double cz)
        {
          double x = (static_cast<double>(i) + 0.5) / static_cast<double>(gridSize.x) - cx;
          double y = (static_cast<double>(j) + 0.5) / static_cast<double>(gridSize.y) - cy;
          double z = (static_cast<double>(k) + 0.5) / static_cast<double>(gridSize.z) - cz;
          return x * x + y * y + z * z;
        };
        field[(k * gridSize.y + j) * gridSize.x + i] =
            static_cast<float>(std::min(at(0.5, 0.5, 0.5) - 0.12, 1.4 * at(0.62, 0.48, 0.51) - 0.05));
      }
    }
  }

  std::vector<double3> cpuGradients;
  std::vector<double3> gpuGradients;
  std::vector<double3> cpu = EnergyIsosurface::trianglesOfIsosurface(field, gridSize, 0.0, &cpuGradients);
  std::vector<double3> gpu = trianglesOfLewinerIsosurface(field, gridSize, 0.0, &gpuGradients);

  Crystal box;
  box.unitCell = UnitCell(10.0, 10.0, 10.0);
  IsosurfaceArea cpuArea = EnergyIsosurface::areaOfIsosurface(box, field, gridSize, 0.0);
  IsosurfaceArea gpuArea = accumulateTriangleAreas(box.unitCell.cell, gridSize, gpu, gpuGradients,
                                                   FieldSense::GrowsIntoSolid);

  ASSERT_GT(cpuArea.numberOfTriangles, 100uz);
  EXPECT_NEAR(static_cast<double>(gpu.size() / 3), static_cast<double>(cpu.size() / 3),
              0.02 * static_cast<double>(cpu.size() / 3));
  EXPECT_NEAR(gpuArea.area, cpuArea.area, 0.01 * cpuArea.area);
}
