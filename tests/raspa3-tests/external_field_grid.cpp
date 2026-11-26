#include <gtest/gtest.h>

import std;

import double3;
import uint3;
import units;
import atom;
import forcefield;
import forcefield_settings;
import component;
import simulationbox;
import system;
import running_energy;
import interpolation_energy_grid;
import interactions_external_field;
import triquintic_derivatives_external_field;

namespace
{
ForceField makeForceField()
{
  return ForceField({{"O", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot,
                    1.0, 1.0, 1.0, false, false, false);
}
}  // namespace

TEST(EXTERNAL_FIELD_GRID, analytic_muller_brown_grid_matches_the_field_with_a_shifted_origin)
{
  // An analytic field is tabulated on the grid size of the settings (no CUBE file is read), and the tabulated
  // values (and derivatives) must be those of the Muller-Brown surface shifted by the potential-energy-surface
  // origin, as evaluated without the grid.
  for (ForceFieldSettings::InterpolationScheme scheme :
       {ForceFieldSettings::InterpolationScheme::Tricubic, ForceFieldSettings::InterpolationScheme::Triquintic})
  {
    ForceField forceField = makeForceField();
    forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::MullerBrown;
    forceField.potentialEnergySurfaceOrigin = double3(-1.5, -0.5, 0.0);
    forceField.externalFieldGridFileName = "this_cube_file_does_not_exist.cube";
    forceField.settings.useExternalFieldGrid = true;
    forceField.settings.interpolationScheme = scheme;
    forceField.settings.numberOfExternalFieldGridPoints = uint3(105, 105, 3);
    forceField.settings.numberOfGridTestPoints = 0;
    const std::uint16_t type = static_cast<std::uint16_t>(*forceField.findPseudoAtom("O"));
    Component atomic(forceField, "atomic", 64.0, 1.0e6, 0.02,
                     {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, type, 0, false, false)}, {}, {}, 5, 21);
    System system(forceField, SimulationBox(2.6, 2.6, 2.6), true, 300.0, 1e4, 1.0, {}, {atomic}, {{}}, {0}, 5);

    std::ostringstream stream;
    system.createExternalFieldInterpolationGrid(stream, 0);
    ASSERT_TRUE(system.externalFieldInterpolationGrid.has_value());

    for (const double3& pos : {double3(1.2, 0.6, 1.0), double3(2.0, 0.55, 0.3), double3(0.9, 1.9, 2.1),
                               double3(0.53, 2.31, 1.7), double3(1.77, 1.03, 0.8)})
    {
      const double analytic = Interactions::analyticExternalFieldEnergy(forceField, pos);
      EXPECT_NEAR(system.externalFieldInterpolationGrid->interpolate(pos), analytic, 1.0e-5 * std::abs(analytic));
    }
  }
}

TEST(EXTERNAL_FIELD_GRID, cube_file_skips_the_atom_lines_and_stores_x_fastest)
{
  // A CUBE file lists its voxels with z fastest after the atom lines; the grid stores x fastest (std::layout_left)
  // and holds only the energy at each node, so a Hermite scheme (which needs the derivatives too) is refused.
  const uint3 dims(3, 4, 5);
  const std::filesystem::path filename =
      std::filesystem::temp_directory_path() / "raspa3_external_field_grid_test.cube";
  {
    std::ofstream file(filename);
    std::print(file, "comment\ncomment\n2 0.0 0.0 0.0\n");
    std::print(file, "-3 5.0 0.0 0.0\n-4 0.0 3.333333333333 0.0\n-5 0.0 0.0 2.5\n");
    std::print(file, "8 0.0 1.0 1.0 1.0\n1 0.0 2.0 2.0 2.0\n");
    for (std::size_t i = 0; i < dims.x; ++i)
    {
      for (std::size_t j = 0; j < dims.y; ++j)
      {
        for (std::size_t k = 0; k < dims.z; ++k)
        {
          std::print(file, "{}\n", 100.0 * static_cast<double>(i) + 10.0 * static_cast<double>(j) + k);
        }
      }
    }
  }

  const std::vector<double> grid = InterpolationEnergyGrid::parseExternalFieldGridCube(filename.string());
  ASSERT_EQ(grid.size(), dims.x * dims.y * dims.z);
  for (std::size_t i = 0; i < dims.x; ++i)
  {
    for (std::size_t j = 0; j < dims.y; ++j)
    {
      for (std::size_t k = 0; k < dims.z; ++k)
      {
        EXPECT_DOUBLE_EQ(grid[i + dims.x * (j + dims.y * k)],
                         (100.0 * static_cast<double>(i) + 10.0 * static_cast<double>(j) + k) * Units::KelvinToEnergy);
      }
    }
  }

  ForceField forceField = makeForceField();
  forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::GridFile;
  forceField.externalFieldGridFileName = filename.string();
  const SimulationBox box(10.0, 10.0, 10.0);
  std::ostringstream stream;

  InterpolationEnergyGrid polynomial(box, double3(), dims, ForceFieldSettings::InterpolationScheme::Polynomial);
  polynomial.makeExternalFieldInterpolationGrid(stream, forceField, box);
  ASSERT_EQ(polynomial.data.size(), dims.x * dims.y * dims.z);
  const double3 node = box.cell * double3(1.0 / 2.0, 2.0 / 3.0, 3.0 / 4.0);
  EXPECT_NEAR(polynomial.interpolate(node), 123.0 * Units::KelvinToEnergy, 1.0e-8);

  InterpolationEnergyGrid tricubic(box, double3(), dims, ForceFieldSettings::InterpolationScheme::Tricubic);
  EXPECT_THROW(tricubic.makeExternalFieldInterpolationGrid(stream, forceField, box), std::runtime_error);

  std::filesystem::remove(filename);
}

// The Muller-Brown surface has its three minima at the published points and values (Muller and Brown, Theor.
// Chim. Acta 53, 75 (1979)), and its fourth term is a confining +15 exp(0.7 dx^2 + 0.6 dx dy + 0.7 dy^2): with
// A_4 = -15 and a_4 = -0.7 the minima were off by up to 73 and the surface fell to minus infinity along y.
TEST(EXTERNAL_FIELD_GRID, muller_brown_has_the_published_minima_and_is_bounded_below)
{
  ForceField forceField = makeForceField();
  forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::MullerBrown;
  forceField.potentialEnergySurfaceOrigin = double3(0.0, 0.0, 0.0);

  const std::array<std::pair<double3, double>, 3> minima{{{double3(-0.558224, 1.441726, 0.0), -146.6995},
                                                          {double3(0.623499, 0.028038, 0.0), -108.1667},
                                                          {double3(-0.050011, 0.466694, 0.0), -80.7678}}};
  for (const auto& [position, energy] : minima)
  {
    EXPECT_NEAR(Interactions::analyticExternalFieldEnergy(forceField, position), energy, 1.0e-3);
    const auto derivatives = Potentials::potentialTriquinticMullerBrown(position);
    EXPECT_NEAR(std::get<0>(derivatives), energy, 1.0e-3);
    EXPECT_NEAR(std::get<1>(derivatives)[0], 0.0, 1.0e-2);
    EXPECT_NEAR(std::get<1>(derivatives)[1], 0.0, 1.0e-2);
  }

  for (double3 far : {double3(0.0, 5.0, 0.0), double3(0.0, -5.0, 0.0), double3(5.0, 5.0, 0.0), double3(-5.0, 0.0, 0.0)})
  {
    EXPECT_GT(Interactions::analyticExternalFieldEnergy(forceField, far), -150.0);
  }
}

// The grid of an analytic surface holds it over the cell only: the surface is no periodic function, and the energy of
// an atom beyond the cell is that of the surface itself and not of its wrapped image (which differed by hundreds of
// Kelvin for Muller-Brown, and jumped at every face of the cell).
TEST(EXTERNAL_FIELD_GRID, analytic_grid_is_not_wrapped_beyond_the_cell)
{
  ForceField forceField = makeForceField();
  forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::MullerBrown;
  forceField.potentialEnergySurfaceOrigin = double3(-1.5, -0.5, 0.0);
  forceField.externalFieldGridFileName = "this_cube_file_does_not_exist.cube";
  forceField.settings.useExternalFieldGrid = true;
  forceField.settings.interpolationScheme = ForceFieldSettings::InterpolationScheme::Tricubic;
  forceField.settings.numberOfExternalFieldGridPoints = uint3(105, 105, 3);
  forceField.settings.numberOfGridTestPoints = 0;
  const std::uint16_t type = static_cast<std::uint16_t>(*forceField.findPseudoAtom("O"));
  Component atomic(forceField, "atomic", 64.0, 1.0e6, 0.02, {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, type, 0, false, false)},
                   {}, {}, 5, 21);
  System system(forceField, SimulationBox(2.6, 2.6, 2.6), true, 300.0, 1e4, 1.0, {}, {atomic}, {{}}, {0}, 5);
  std::ostringstream stream;
  system.createExternalFieldInterpolationGrid(stream, 0);
  ASSERT_TRUE(system.externalFieldInterpolationGrid.has_value());

  for (const double3& pos : {double3(1.2, 0.6, 1.0), double3(2.9, 0.6, 1.0), double3(-0.4, 1.9, 1.0)})
  {
    const std::vector<Atom> atoms{Atom(pos, 0.0, 1.0, 0, type, 0, false, false)};
    RunningEnergy energy{};
    Interactions::computeExternalFieldEnergy(true, forceField, system.simulationBox, atoms, energy,
                                             system.externalFieldInterpolationGrid);
    const double analytic = Interactions::analyticExternalFieldEnergy(forceField, pos);
    EXPECT_NEAR(energy.externalFieldVDW, analytic, 1.0e-5 * std::abs(analytic));
  }
}

// The first voxel of a CUBE file lies at the origin in its header: a file whose grid starts at (1, 2, 3) holds the
// field at a node there, and not at the node of the cell at the same fractional position.
TEST(EXTERNAL_FIELD_GRID, cube_file_grid_starts_at_its_origin)
{
  const uint3 dims(3, 4, 5);
  const std::filesystem::path filename =
      std::filesystem::temp_directory_path() / "raspa3_external_field_grid_origin_test.cube";
  const double3 origin(1.0, 2.0, 3.0);
  {
    std::ofstream file(filename);
    std::print(file, "comment\ncomment\n0 {} {} {}\n", origin.x, origin.y, origin.z);
    std::print(file, "-3 5.0 0.0 0.0\n-4 0.0 3.333333333333 0.0\n-5 0.0 0.0 2.5\n");
    for (std::size_t i = 0; i < dims.x; ++i)
    {
      for (std::size_t j = 0; j < dims.y; ++j)
      {
        for (std::size_t k = 0; k < dims.z; ++k)
        {
          std::print(file, "{}\n", 100.0 * static_cast<double>(i) + 10.0 * static_cast<double>(j) + k);
        }
      }
    }
  }
  // the negative counts give the steps (and so the origin) in Angstrom
  EXPECT_NEAR((InterpolationEnergyGrid::parseExternalFieldGridOrigin(filename.string()) - origin).length(), 0.0,
              1.0e-12);

  ForceField forceField = makeForceField();
  forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::GridFile;
  forceField.externalFieldGridFileName = filename.string();
  const SimulationBox box(10.0, 10.0, 10.0);
  std::ostringstream stream;
  InterpolationEnergyGrid polynomial(box, double3(), dims, ForceFieldSettings::InterpolationScheme::Polynomial);
  polynomial.makeExternalFieldInterpolationGrid(stream, forceField, box);
  const double3 node = box.cell * double3(1.0 / 2.0, 2.0 / 3.0, 3.0 / 4.0) + origin;
  EXPECT_NEAR(polynomial.interpolate(node), 123.0 * Units::KelvinToEnergy, 1.0e-8);

  std::filesystem::remove(filename);
}

// The overlap sentinel is a mark in the energy slot of a vertex inside a framework atom and stays when the grid is
// scaled; every other entry is scaled. The derivatives of a vertex just outside an atom are stored as computed and
// can exceed the criterion: tested entry by entry against it, a positive one stayed while the negative ones were
// scaled, and the interpolant no longer was the scaled one.
TEST(EXTERNAL_FIELD_GRID, scaled_grid_keeps_only_the_overlap_sentinel)
{
  const double criterion = 1.0e5;
  InterpolationEnergyGrid grid(SimulationBox(10.0, 10.0, 10.0), double3(), uint3(2, 2, 2),
                               ForceFieldSettings::InterpolationScheme::Tricubic);
  grid.data[0] = 100.0;
  grid.data[1] = 3.0 * criterion;
  grid.data[2] = -3.0 * criterion;
  grid.data[8] = 2.0 * criterion;
  grid.data[9] = 0.5 * criterion;

  grid.scaleEnergy(0.5, criterion);
  EXPECT_DOUBLE_EQ(grid.data[0], 50.0);
  EXPECT_DOUBLE_EQ(grid.data[1], 1.5 * criterion);
  EXPECT_DOUBLE_EQ(grid.data[2], -1.5 * criterion);
  EXPECT_DOUBLE_EQ(grid.data[8], 2.0 * criterion);
  EXPECT_DOUBLE_EQ(grid.data[9], 0.25 * criterion);
}
