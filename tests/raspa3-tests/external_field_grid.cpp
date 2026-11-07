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
import interpolation_energy_grid;
import interactions_external_field;

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
