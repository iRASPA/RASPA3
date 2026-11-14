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


// The same sphere with the atom on the cell's corner, so that the sheet crosses all three periodic faces. The
// triangles straddling a face have their corners wrapped into the cell by the refinement, and an edge taken
// between two wrapped corners without the nearest image runs the length of the cell: such a triangle was thrown
// out as implausibly large, and the sheet lost its area along every face.
TEST(energy_well_surface, isolated_atom_on_the_cell_corner_has_the_same_sheet)
{
  const double edge = 16.0;
  const uint3 gridSize{96, 96, 96};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = oneAtom(edge);
  framework.atoms.front().position = double3(0.0, 0.0, 0.0);
  framework.fractionalPositions.front() = double3(0.0, 0.0, 0.0);

  WellField field =
      computeWellField(interactions, framework, probeOf(interactions, "X"), gridSize, 1, 0.0, {}, 1000.0, ceiling);
  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.5 * epsilon, {});

  ASSERT_GT(surface.numberOfTriangles, 100u);
  const double rmin = wellContactPrefactor * sigma;
  const double expected = 4.0 * std::numbers::pi * rmin * rmin;
  EXPECT_NEAR(surface.area, expected, 0.08 * expected);
  EXPECT_EQ(surface.numberOfRejectedTriangles, 0u) << "no triangle of a sphere is implausibly large";
  EXPECT_NEAR(surface.meanDepth, -epsilon, 0.05 * epsilon);
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

// A straight tube along [110] of a cubic cell wraps a and b, but along one translation (a + b): a file with that
// period, not the slit the two axis flags made of it (no length at all) in a cell whose c is the shortest edge.
TEST(energy_well_surface, a_diagonal_tube_packs_as_a_file)
{
  const double l = 12.0;
  const uint3 gridSize{32, 32, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = emptyBox(UnitCell(l, l, l));
  addAtom(framework, 0.5, 0.5, 0.5);

  WellField field = blankField(gridSize, framework.unitCell);
  field.probe = probeOf(interactions, "X");
  for (std::size_t k = 14; k <= 18; ++k)
  {
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < gridSize.x; ++i)
      {
        // within two voxels of the line i - j = 0 (mod 32)
        const int offset = (static_cast<int>(i) - static_cast<int>(j) + 48) % 32 - 16;
        if (std::abs(offset) <= 2) markFilament(field, field.voxelIndex(i, j, k));
      }
    }
  }

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});

  const double period = std::sqrt(2.0) * l;
  EXPECT_GT(surface.filamentVolume, 0.0);
  EXPECT_NEAR(surface.filamentRidgeLength, period, 0.20 * period);
  EXPECT_LT(surface.filamentMedialArea, 0.20 * l * l);
}

// An inclined (111) slab wraps all three axes along two translations: a slit, packed on its midplane of area
// sqrt(3) l^2 per cell, not the files the three axis flags made of it.
TEST(energy_well_surface, an_inclined_slit_packs_on_the_midplane)
{
  const double l = 12.0;
  const uint3 gridSize{32, 32, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = emptyBox(UnitCell(l, l, l));
  addAtom(framework, 0.5, 0.5, 0.5);

  WellField field = blankField(gridSize, framework.unitCell);
  field.probe = probeOf(interactions, "X");
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < gridSize.x; ++i)
      {
        // within a slab about the plane i + j + k = 16 (mod 32), a few voxels thick
        const int offset = static_cast<int>((i + j + k + 16) % 32) - 16;
        if (std::abs(offset) <= 3) markFilament(field, field.voxelIndex(i, j, k));
      }
    }
  }

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});

  const double expected = std::sqrt(3.0) * l * l;
  EXPECT_GT(surface.filamentVolume, 0.0);
  EXPECT_GT(surface.filamentMedialArea, 0.5 * expected);
  EXPECT_LT(surface.filamentRidgeLength, 0.20 * l);
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


// Files that run side by side fill a cross-section, which has two dimensions: how many there may be is
// set by its area and not by the shorter of its two edges. Four tubes along z stand in a row along x, in a
// cell that is long that way and only six Ångström the other. Held to what fits along the six Ångström,
// which is two, half of them went uncounted.
TEST(energy_well_surface, tubes_in_a_row_across_a_flat_cell_are_all_counted)
{
  const double lx = 24.0, ly = 6.0, lz = 16.0;
  const uint3 gridSize{64, 16, 32};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = emptyBox(UnitCell(lx, ly, lz));
  addAtom(framework, 0.5, 0.5, 0.5);

  WellField field = blankField(gridSize, framework.unitCell);
  field.probe = probeOf(interactions, "X");
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 6; j <= 10; ++j)
    {
      for (std::size_t tube = 0; tube < 4; ++tube)
      {
        for (std::size_t i = 6 + 16 * tube; i <= 10 + 16 * tube; ++i)
        {
          markFilament(field, field.voxelIndex(i, j, k));
        }
      }
    }
  }

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});

  // four files of four sites each, the sixteen Ångström holding four
  const double expected = 4.0 * 4.0 * nitrogenLinearSpacing;
  EXPECT_GT(surface.filamentVolume, 0.0);
  EXPECT_NEAR(surface.filamentRidgeLength, expected, 1.0e-9);
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


// The device keeps the energy of every orientation, as the CPU does, so that the molecular density functional
// can be fed from either. Held at the ceiling orientation by orientation, and agreeing with the CPU's double
// precision to what a float carries; the Helmholtz average and the rest of the field are checked beside it.
TEST(energy_well_surface, the_device_keeps_the_energy_of_every_orientation)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  const double edge = 16.0;
  const uint3 gridSize{24, 24, 24};
  const double bond = 2.4;
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = oneAtom(edge);

  LinearProbe dimer;
  dimer.name = "dimer";
  dimer.headTailSymmetric = true;
  dimer.sites.push_back(LinearProbe::Site{0, -0.5 * bond, 0.0, "X"});
  dimer.sites.push_back(LinearProbe::Site{0, 0.5 * bond, 0.0, "X"});
  const std::size_t orientations = 16;

  WellField cpu = WellFieldCPU::compute(interactions, framework, dimer, gridSize, orientations, 0.5 * epsilon, {},
                                        1000.0, ceiling);
  WellField gpu = WellFieldOpenCL::compute(interactions, framework, dimer, gridSize, orientations, 0.5 * epsilon, {},
                                           1000.0, ceiling);

  ASSERT_EQ(cpu.numberOfOrientations, orientations);
  ASSERT_EQ(gpu.numberOfOrientations, orientations);
  ASSERT_EQ(cpu.orientationEnergy.size(), cpu.numberOfVoxels() * orientations);
  ASSERT_EQ(gpu.orientationEnergy.size(), cpu.orientationEnergy.size()) << "the device kept no orientations";

  double worstOrientation = 0.0;
  double worstEnergy = 0.0;
  for (std::size_t i = 0; i < cpu.orientationEnergy.size(); ++i)
  {
    double a = static_cast<double>(cpu.orientationEnergy[i]);
    double b = static_cast<double>(gpu.orientationEnergy[i]);
    // Measured against the well depth: near the zero crossing of a pair sum the repulsion and the attraction
    // cancel to a thousandth of either, which is what a float carries, and a value near nought is no yardstick.
    worstOrientation = std::max(worstOrientation, std::abs(a - b) / std::max(epsilon, std::abs(a)));
  }
  for (std::size_t i = 0; i < cpu.numberOfVoxels(); ++i)
  {
    double a = static_cast<double>(cpu.energy[i]);
    double b = static_cast<double>(gpu.energy[i]);
    worstEnergy = std::max(worstEnergy, std::abs(a - b) / std::max(1.0, std::abs(a)));
  }
  EXPECT_LT(worstOrientation, 1.0e-4);
  EXPECT_LT(worstEnergy, 1.0e-4);

  // The Helmholtz average of the kept orientations is the field's own energy, so the two columns of the same
  // field agree with each other as well.
  const double beta = 1.0 / (0.5 * epsilon);
  double worstAverage = 0.0;
  for (std::size_t voxel = 0; voxel < gpu.numberOfVoxels(); ++voxel)
  {
    double least = std::numeric_limits<double>::max();
    for (std::size_t o = 0; o < orientations; ++o)
      least = std::min(least, static_cast<double>(gpu.orientationEnergy[voxel * orientations + o]));
    double sum = 0.0;
    for (std::size_t o = 0; o < orientations; ++o)
      sum += std::exp(-beta * (static_cast<double>(gpu.orientationEnergy[voxel * orientations + o]) - least));
    double average = least - (0.5 * epsilon) * std::log(sum / static_cast<double>(orientations));
    average = std::min(average, ceiling);
    worstAverage = std::max(worstAverage, std::abs(average - static_cast<double>(gpu.energy[voxel])) /
                                              std::max(1.0, std::abs(average)));
  }
  EXPECT_LT(worstAverage, 1.0e-4);
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


// The charged probes of the tests below, and what they are run on. The nitrogen model is TraPPE's in shape:
// two dispersing sites that carry charge, and the charge that balances them on a massless site between
// them with no Lennard-Jones of its own.
namespace
{
constexpr double coulombConversion = 167101.0;
constexpr double nitrogenSigma = 3.30;
constexpr double nitrogenBond = 0.55;
constexpr double nitrogenCharge = -0.405;
constexpr double oxygenCharge = -1.025;

PairInteractions nitrogenInteractions()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "Q"};
  interactions.charges = {0.0, -2.0 * nitrogenCharge};
  interactions.parameters = {PairParameters{nitrogenSigma, epsilon, 0.0}, PairParameters{0.0, 0.0, 0.0},
                             PairParameters{0.0, 0.0, 0.0}, PairParameters{0.0, 0.0, 0.0}};
  interactions.cutOffVDW = 12.0;
  interactions.cutOffCoulomb = 12.0;
  return interactions;
}

LinearProbe nitrogenProbe()
{
  LinearProbe probe;
  probe.name = "N2";
  probe.headTailSymmetric = true;
  probe.sites.push_back(LinearProbe::Site{0, -nitrogenBond, nitrogenCharge, "X"});
  probe.sites.push_back(LinearProbe::Site{1, 0.0, -2.0 * nitrogenCharge, "Q"});
  probe.sites.push_back(LinearProbe::Site{0, nitrogenBond, nitrogenCharge, "X"});
  return probe;
}

// One atom of the given charge at the middle of a cube 40 Å on a side.
Crystal oneChargedAtom(double charge)
{
  Crystal framework = oneAtom(40.0);
  framework.atoms[0].charge = charge;
  return framework;
}

// A potential with all of the electrostatics in its near half: no screening to speak of, and a smooth half
// that is the same everywhere.
ElectrostaticPotentialGrid unscreenedPotential(const Crystal &framework, uint3 gridSize, float smooth)
{
  ElectrostaticPotentialGrid potential;
  potential.gridSize = gridSize;
  potential.unitCell = framework.unitCell;
  potential.alpha = 1.0e-6;
  potential.cutOff = 12.0;
  const std::size_t points = static_cast<std::size_t>(gridSize.x) * static_cast<std::size_t>(gridSize.y) *
                             static_cast<std::size_t>(gridSize.z);
  potential.smoothPotential.assign(points, smooth);
  return potential;
}

// The field of one orientation, along z, from the host and from the device where there is one.
std::vector<std::pair<std::string, WellField>> onHostAndDevice(const PairInteractions &interactions,
                                                               const Crystal &framework, const LinearProbe &probe,
                                                               uint3 gridSize,
                                                               const ElectrostaticPotentialGrid &potential)
{
  std::vector<std::pair<std::string, WellField>> fields;
  fields.emplace_back("host", WellFieldCPU::compute(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0,
                                                    ceiling, &potential, coulombConversion));
  OpenCL::initialize();
  if (OpenCL::clContext.has_value())
  {
    fields.emplace_back("device", WellFieldOpenCL::compute(interactions, framework, probe, gridSize, 1, 0.0, {},
                                                           1000.0, ceiling, &potential, coulombConversion));
  }
  return fields;
}

double lennardJones(double size, double strength, double r)
{
  const double ratio6 = std::pow(size / r, 6);
  return 4.0 * strength * ratio6 * (ratio6 - 1.0);
}
}  // namespace


// Sigma is where a pair's dispersion crosses nought, not a wall: a site a little inside it costs a few
// epsilon, and a pose like that is one a simulation visits. So the energy there is the model's energy, the
// dispersion and the electrostatics together, and it runs through sigma without a step. Refusing the
// electrostatics from sigma inwards, as this field once did, takes the whole electrostatic energy of the
// orientation out at once: several hundred kelvin here, with a barrier turned into a well end on and a well
// cut off broadside.
//
// A nitrogen is brought up to one charged atom along its own axis and across it, a tenth of an Ångström at a
// time, from well outside sigma to deep inside it, and the field is held to the sum over its three sites.
TEST(energy_well_surface, the_energy_of_a_charged_probe_runs_through_sigma_without_a_step)
{
  const PairInteractions interactions = nitrogenInteractions();
  const Crystal framework = oneChargedAtom(oxygenCharge);
  const LinearProbe probe = nitrogenProbe();

  auto coulomb = [](double charge, double r) { return coulombConversion * charge * oxygenCharge / r; };

  for (bool endOn : {true, false})
  {
    // The axis of the probe is z. The fine direction of the grid is the one it is brought up along.
    const uint3 gridSize = endOn ? uint3{40, 40, 400} : uint3{400, 40, 40};
    const ElectrostaticPotentialGrid potential = unscreenedPotential(framework, gridSize, 0.0f);

    for (const auto &[backend, field] : onHostAndDevice(interactions, framework, probe, gridSize, potential))
    {
      bool wentInsideSigma = false;
      bool stayedOutsideSigma = false;
      for (std::size_t step = 28; step <= 50; ++step)
      {
        const double distance = 0.1 * static_cast<double>(step);
        const double nearer = endOn ? distance - nitrogenBond
                                    : std::sqrt(distance * distance + nitrogenBond * nitrogenBond);
        const double further = endOn ? distance + nitrogenBond : nearer;
        wentInsideSigma = wentInsideSigma || nearer < 0.9 * nitrogenSigma;
        stayedOutsideSigma = stayedOutsideSigma || nearer > 1.1 * nitrogenSigma;

        const double dispersion = lennardJones(nitrogenSigma, epsilon, nearer) +
                                  lennardJones(nitrogenSigma, epsilon, further);
        const double electrostatics = coulomb(nitrogenCharge, nearer) + coulomb(nitrogenCharge, further) +
                                      coulomb(-2.0 * nitrogenCharge, distance);
        const double expected = dispersion + electrostatics;

        const std::size_t voxel = endOn ? ((200 + step) * gridSize.y + 20) * gridSize.x + 20
                                        : (20 * gridSize.y + 20) * gridSize.x + 200 + step;
        EXPECT_NEAR(static_cast<double>(field.energy[voxel]), expected, 0.5 + 1.0e-3 * std::abs(expected))
            << backend << (endOn ? ", end on, " : ", broadside, ") << distance << " Å off, nearer site at "
            << nearer / nitrogenSigma << " sigma, electrostatics " << electrostatics << " K";
      }
      EXPECT_TRUE(wentInsideSigma && stayedOutsideSigma);
    }
  }
}


// The same, with the smooth half of the potential in it and a probe that has a net charge, so that there is
// no mistaking whether the electrostatics were counted: one bare charge between two uncharged dispersing
// sites, both of them a little inside sigma of an atom of the opposite charge. The near half and the smooth
// half are both there in full, and the contact distance has not moved, since charges change how deep a well
// is and not where the wall stands.
TEST(energy_well_surface, a_pose_inside_sigma_keeps_its_electrostatics)
{
  const double edge = 40.0;
  const uint3 gridSize{40, 40, 40};
  const double coreSigma = 3.10;
  const double bond = 0.55;
  const float smooth = 1.0e5f;

  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "Q"};
  interactions.charges = {0.0, 1.0};
  interactions.parameters = {PairParameters{coreSigma, epsilon, 0.0}, PairParameters{0.0, 0.0, 0.0},
                             PairParameters{0.0, 0.0, 0.0}, PairParameters{0.0, 0.0, 0.0}};
  interactions.cutOffVDW = 12.0;
  interactions.cutOffCoulomb = 12.0;

  const Crystal framework = oneChargedAtom(-1.0);

  LinearProbe probe;
  probe.name = "N2";
  probe.headTailSymmetric = true;
  probe.sites.push_back(LinearProbe::Site{0, -bond, 0.0, "X"});
  probe.sites.push_back(LinearProbe::Site{1, 0.0, 1.0, "Q"});
  probe.sites.push_back(LinearProbe::Site{0, bond, 0.0, "X"});

  const ElectrostaticPotentialGrid potential = unscreenedPotential(framework, gridSize, smooth);

  WellField without = computeWellField(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling);

  // The centre 3 Å off along x, across the axis: both dispersing sites at 3.05 Å, inside a sigma of 3.10.
  const std::size_t voxel = voxelOf(gridSize, double3(23.0 / edge, 0.5, 0.5));
  const double uncharged = static_cast<double>(without.energy[voxel]);
  const double electrostatics = -coulombConversion / 3.0 + static_cast<double>(smooth);
  EXPECT_NEAR(uncharged, 2.0 * lennardJones(coreSigma, epsilon, std::sqrt(9.0 + bond * bond)), 1.0e-2);

  for (const auto &[backend, with] : onHostAndDevice(interactions, framework, probe, gridSize, potential))
  {
    const double charged = static_cast<double>(with.energy[voxel]);
    EXPECT_NEAR(charged - uncharged, electrostatics, 1.0e-3 * std::abs(electrostatics)) << backend;
    EXPECT_NEAR(with.distance[voxel], without.distance[voxel], backend == "host" ? 0.0f : 1.0e-4f) << backend;
  }
}


// A nitrogen pressed end on into an oxygen, the nearer site at three quarters of sigma, is high on the
// repulsive wall and nothing else: twelve thousand kelvin of dispersion, and the electrostatics of the whole
// quadrupole on top of it, which against a negative atom end on are repulsive as well. It is not at the
// ceiling, which is for overlaps, and it is not a well, which is what the charge at the centre would make
// of it if it were counted and the two it balances were not.
TEST(energy_well_surface, an_overlapping_n2_orientation_is_not_a_coulomb_well)
{
  const double edge = 40.0;
  const uint3 gridSize{40, 40, 40};

  const PairInteractions interactions = nitrogenInteractions();
  const Crystal framework = oneChargedAtom(oxygenCharge);
  const LinearProbe probe = nitrogenProbe();
  const ElectrostaticPotentialGrid potential = unscreenedPotential(framework, gridSize, 0.0f);

  WellField field = computeWellField(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling,
                                     &potential, coulombConversion);

  // The axis is z. The centre 3 Å off along it: the nearer site at 2.45 Å, the centre well outside its floor.
  const std::size_t overlapping = voxelOf(gridSize, double3(0.5, 0.5, 23.0 / edge));
  const double overlapEnergy = static_cast<double>(field.energy[overlapping]);
  const double dispersion = lennardJones(nitrogenSigma, epsilon, 3.0 - nitrogenBond) +
                            lennardJones(nitrogenSigma, epsilon, 3.0 + nitrogenBond);
  const double quadrupole = coulombConversion * oxygenCharge * nitrogenCharge *
                            (1.0 / (3.0 - nitrogenBond) + 1.0 / (3.0 + nitrogenBond) - 2.0 / 3.0);
  const double monopole = coulombConversion * oxygenCharge * -2.0 * nitrogenCharge / 3.0;
  ASSERT_GT(quadrupole, 0.0);
  ASSERT_LT(monopole, -40000.0);
  EXPECT_NEAR(overlapEnergy, dispersion + quadrupole, 1.0e-3 * (dispersion + quadrupole));
  EXPECT_GT(overlapEnergy, dispersion);
  EXPECT_LT(overlapEnergy, ceiling * 0.5);

  // Between sigma and the bottom of the well, and beyond it: the full quadrupole, and nothing drastic.
  const std::size_t shoulder = voxelOf(gridSize, double3(0.5, 0.5, 24.0 / edge));
  const double shoulderEnergy = static_cast<double>(field.energy[shoulder]);
  EXPECT_LT(shoulderEnergy, ceiling * 0.5);
  EXPECT_GT(shoulderEnergy, -5000.0);

  const std::size_t clear = voxelOf(gridSize, double3(0.5, 0.5, 26.0 / edge));
  const double clearEnergy = static_cast<double>(field.energy[clear]);
  EXPECT_LT(clearEnergy, ceiling * 0.5);
  EXPECT_GT(clearEnergy, -5000.0);
}


// Where the electrostatics are refused is at an overlap, a pair whose dispersion has reached the ceiling,
// and what the refusal is for is this. The dispersion of such a pair is held at the ceiling and can go no
// higher, while the Coulomb term of two opposite charges goes on down without limit as they close, so that
// left in it digs a hole under the ceiling deeper than any well in the field.
//
// Two dispersing sites that carry charge, and no site between them to be caught by its own floor. One of
// them is put on top of an atom of the opposite charge.
TEST(energy_well_surface, a_charge_pressed_onto_an_atom_does_not_dig_under_the_ceiling)
{
  const uint3 gridSize{40, 40, 80};
  const double siteCharge = -0.5;
  const double bond = 0.5;

  PairInteractions interactions = oneType(sigma, epsilon);
  interactions.cutOffVDW = 12.0;
  interactions.cutOffCoulomb = 12.0;

  const Crystal framework = oneChargedAtom(1.0);

  LinearProbe probe;
  probe.name = "dumbbell";
  probe.headTailSymmetric = true;
  probe.sites.push_back(LinearProbe::Site{0, -bond, siteCharge, "X"});
  probe.sites.push_back(LinearProbe::Site{0, bond, siteCharge, "X"});

  const ElectrostaticPotentialGrid potential = unscreenedPotential(framework, gridSize, 0.0f);

  // What the pair would come to at the separation the screened term is held off nought by.
  const double hole = coulombConversion * siteCharge / 1.0e-4;
  ASSERT_LT(hole, -10.0 * ceiling);

  // The deepest the field can honestly be: both sites at the bottom of their wells, and closer to the atom
  // than either can get without paying for it.
  const double deepestHonestWell = -2.0 * epsilon + 2.0 * coulombConversion * siteCharge / (0.8 * sigma);

  for (const auto &[backend, field] : onHostAndDevice(interactions, framework, probe, gridSize, potential))
  {
    // Half an Ångström along the axis from the atom, which is where one of the sites is.
    const std::size_t onTop = (41 * gridSize.y + 20) * gridSize.x + 20;
    EXPECT_FLOAT_EQ(field.energy[onTop], static_cast<float>(ceiling)) << backend;
    EXPECT_GT(field.deepestEnergy(), deepestHonestWell) << backend;
    EXPECT_LT(field.deepestEnergy(), -epsilon) << backend;
  }
}


// A site may have dispersion, and still have none against one of the atoms: a framework atom that is given
// a charge and no Lennard-Jones, as the silicon of a zeolite often is, with the oxygens around it left to
// keep everything off. Against that atom the site is a bare charge, with no wall of its own to keep its 1/r
// off the nucleus, so there it has the one Ångström floor a massless charge site has everywhere. Outside
// the floor the pair is its Coulomb term and nothing else, whatever size the mixing rules gave it.
TEST(energy_well_surface, a_site_with_only_charge_against_an_atom_is_held_off_its_nucleus)
{
  const uint3 gridSize{40, 40, 80};
  const double siteCharge = -0.5;
  const double bond = 0.5;

  // X, the probe; A, an atom it has dispersion with; B, an atom it has a size with and no strength.
  PairInteractions interactions;
  interactions.numberOfTypes = 3;
  interactions.names = {"X", "A", "B"};
  interactions.charges = {0.0, 0.0, 1.0};
  interactions.parameters.assign(9, PairParameters{sigma, 0.0, 0.0});
  interactions.parameters[0 * 3 + 1] = PairParameters{sigma, epsilon, 0.0};
  interactions.parameters[1 * 3 + 0] = PairParameters{sigma, epsilon, 0.0};
  interactions.cutOffVDW = 12.0;
  interactions.cutOffCoulomb = 12.0;

  // B at the middle of the cell, A ten Ångström above it.
  Crystal framework = oneChargedAtom(1.0);
  framework.atoms[0].type = 2;
  framework.atoms.push_back(CrystalAtom{double3(20.0, 20.0, 30.0), 0.0, 1});
  framework.fractionalPositions.push_back(double3(0.5, 0.5, 0.75));

  LinearProbe probe;
  probe.name = "dumbbell";
  probe.headTailSymmetric = true;
  probe.sites.push_back(LinearProbe::Site{0, -bond, siteCharge, "X"});
  probe.sites.push_back(LinearProbe::Site{0, bond, siteCharge, "X"});

  const ElectrostaticPotentialGrid potential = unscreenedPotential(framework, gridSize, 0.0f);

  for (const auto &[backend, field] : onHostAndDevice(interactions, framework, probe, gridSize, potential))
  {
    // One site on top of B: inside the floor, and at the ceiling for it.
    const std::size_t onTop = (41 * gridSize.y + 20) * gridSize.x + 20;
    EXPECT_FLOAT_EQ(field.energy[onTop], static_cast<float>(ceiling)) << backend;

    // The centre 2.5 Å above B: the sites 2 and 3 Å from it, inside the size of the pair and outside the
    // floor, and 8 and 7 Å below A.
    const std::size_t outside = (45 * gridSize.y + 20) * gridSize.x + 20;
    const double expected = coulombConversion * siteCharge * (1.0 / 2.0 + 1.0 / 3.0) +
                            lennardJones(sigma, epsilon, 8.0) + lennardJones(sigma, epsilon, 7.0);
    EXPECT_NEAR(static_cast<double>(field.energy[outside]), expected, 1.0e-3 * std::abs(expected)) << backend;

    // And nowhere is there a hole: the deepest point is no deeper than both sites a floor away from B.
    EXPECT_GT(field.deepestEnergy(), 2.0 * coulombConversion * siteCharge / wellDummyCoreRadius - 2.0 * epsilon)
        << backend;
  }
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

// The two extractors cut the same cubes the same ways, and a cube cut one of the ambiguous ways is given a vertex of
// its own inside it. Its normal has to be the same from both, as the normals on the edges are: the mean of the field's
// gradients on the cube's edges. The processor one took the mean of the unit normals there instead, and on a field
// folded tightly enough to have such vertices the two were fourteen degrees apart on average, on every one of them.
TEST(energy_well_surface, gpu_lewiner_normals_agree_with_cpu_inside_a_cube)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value())
  {
    GTEST_SKIP() << "no OpenCL device";
  }

  // Two dozen plane waves under an envelope that makes the field several times steeper on one side of the cell, so
  // that the gradients on the edges of one cube are of very different lengths.
  const std::size_t n = 24;
  const uint3 gridSize{24, 24, 24};
  struct Wave
  {
    int a, b, c;
    double phase, amplitude;
  };
  std::vector<Wave> waves;
  std::uint64_t state = 0x9E3779B97F4A7C15ull;
  auto next = [&]()
  {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    return static_cast<double>(state >> 11) / 9007199254740992.0;
  };
  for (int w = 0; w < 24; ++w)
  {
    const int a = static_cast<int>(next() * 11.0) - 5;
    const int b = static_cast<int>(next() * 11.0) - 5;
    const int c = static_cast<int>(next() * 11.0) - 5;
    const double phase = 2.0 * std::numbers::pi * next();
    waves.push_back(Wave{a, b, c, phase, 0.2 + next()});
  }
  std::vector<float> field(n * n * n);
  for (std::size_t k = 0; k < n; ++k)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        double value = 0.0;
        for (const Wave &wave : waves)
        {
          value += wave.amplitude * std::cos(2.0 * std::numbers::pi *
                                                 (wave.a * static_cast<double>(i) + wave.b * static_cast<double>(j) +
                                                  wave.c * static_cast<double>(k)) /
                                                 static_cast<double>(n) +
                                             wave.phase);
        }
        const double envelope = std::exp(1.5 * std::sin(2.0 * std::numbers::pi * static_cast<double>(i) / 24.0));
        field[(k * n + j) * n + i] = static_cast<float>(value * envelope);
      }
    }
  }

  std::vector<double3> cpuGradients;
  std::vector<double3> gpuGradients;
  std::vector<double3> cpu = EnergyIsosurface::trianglesOfIsosurface(field, gridSize, 0.0, &cpuGradients);
  std::vector<double3> gpu = trianglesOfLewinerIsosurface(field, gridSize, 0.0, &gpuGradients);
  ASSERT_EQ(gpu.size(), cpu.size()) << "the two have to be one mesh for its vertices to be compared";
  ASSERT_EQ(gpuGradients.size(), cpuGradients.size());

  std::size_t inside = 0;
  for (std::size_t v = 0; v < cpu.size(); ++v)
  {
    ASSERT_LT((cpu[v] - gpu[v]).length(), 1.0e-4) << "vertex " << v;
    int onGridPlanes = 0;
    for (double coordinate : {cpu[v].x, cpu[v].y, cpu[v].z})
    {
      const double steps = coordinate * static_cast<double>(n);
      if (std::abs(steps - std::round(steps)) < 1.0e-9) ++onGridPlanes;
    }
    if (onGridPlanes < 2) ++inside;

    const double alike =
        double3::dot(cpuGradients[v], gpuGradients[v]) / (cpuGradients[v].length() * gpuGradients[v].length());
    EXPECT_GT(alike, std::cos(1.0e-3)) << (onGridPlanes < 2 ? "inside a cube" : "on an edge") << ", vertex " << v;
  }
  EXPECT_GT(inside, 50uz) << "the surface has to have vertices inside its cubes for this to have asked anything";

  // And with the normals one, the shape of the surface is read one way: the same triangles in each column.
  const UnitCell cell(10.0, 10.0, 10.0);
  const CurvatureAreas onProcessor =
      accumulateTriangleAreas(cell.cell, gridSize, cpu, cpuGradients, FieldSense::GrowsIntoSolid).curvature;
  const CurvatureAreas onGPU =
      accumulateTriangleAreas(cell.cell, gridSize, gpu, gpuGradients, FieldSense::GrowsIntoSolid).curvature;
  EXPECT_NEAR(onGPU.convex, onProcessor.convex, 1.0e-3 * onProcessor.total());
  EXPECT_NEAR(onGPU.saddle, onProcessor.saddle, 1.0e-3 * onProcessor.total());
  EXPECT_NEAR(onGPU.concave, onProcessor.concave, 1.0e-3 * onProcessor.total());
  EXPECT_NEAR(onGPU.unresolved, onProcessor.unresolved, 1.0e-3 * onProcessor.total());
}

// A field with no well has no triangles and no deepest well: it started at +infinity and was reset where
// std::isfinite said so, which the fast-math build folds to true, so the report gave a deepest well of inf.
TEST(energy_well_surface, a_surface_without_triangles_has_no_deepest_well)
{
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = emptyBox(UnitCell(12.0, 12.0, 12.0));
  addAtom(framework, 0.5, 0.5, 0.5);
  WellField field = blankField(uint3{16, 16, 16}, framework.unitCell);
  field.probe = probeOf(interactions, "X");

  WellSurface surface = wellSurfaceOfField(framework, interactions, field, 0.0, 0.001, 0.0, {});
  ASSERT_EQ(surface.numberOfTriangles, 0uz);
  EXPECT_EQ(std::to_string(surface.deepestWell), std::to_string(0.0));
}

// More orientations than the device kernel holds (128): the GPU route refused the field ('--gpu --bet
// --orientations 256' failed outright) instead of computing it on the processor.
TEST(energy_well_surface, the_device_route_computes_a_field_with_more_orientations_than_it_holds)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  const uint3 gridSize{8, 8, 8};
  PairInteractions interactions = oneType(sigma, epsilon);
  Crystal framework = oneAtom(16.0);
  LinearProbe dimer;
  dimer.name = "dimer";
  dimer.headTailSymmetric = false;
  dimer.sites.push_back(LinearProbe::Site{0, -1.2, 0.0, "X"});
  dimer.sites.push_back(LinearProbe::Site{0, 1.2, 0.0, "X"});
  const std::size_t orientations = 256;

  WellField cpu =
      WellFieldCPU::compute(interactions, framework, dimer, gridSize, orientations, 0.5 * epsilon, {}, 1000.0, ceiling);
  WellField gpu;
  ASSERT_NO_THROW(gpu = WellFieldOpenCL::compute(interactions, framework, dimer, gridSize, orientations, 0.5 * epsilon,
                                                 {}, 1000.0, ceiling));
  ASSERT_GT(cpu.numberOfOrientations, 128uz);
  EXPECT_EQ(gpu.numberOfOrientations, cpu.numberOfOrientations);
  ASSERT_EQ(gpu.energy.size(), cpu.energy.size());
  for (std::size_t i = 0; i < cpu.energy.size(); ++i) EXPECT_EQ(gpu.energy[i], cpu.energy[i]);
}
