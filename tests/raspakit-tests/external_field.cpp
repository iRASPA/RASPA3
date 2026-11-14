#include <gtest/gtest.h>

import std;

import double3;
import uint3;
import atom;
import forcefield;
import forcefield_settings;
import simulationbox;
import running_energy;
import interpolation_energy_grid;
import interactions_external_field;
import cbmc_external_energy;
import units;
import system;
import component;
import energy_status;

namespace
{
ForceField makeForceField(ForceField::PotentialEnergySurfaceType type)
{
  ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, false, false, false);
  forceField.potentialEnergySurfaceType = type;
  return forceField;
}

Atom makeAtom(double3 position, double scalingVDW = 1.0, std::uint8_t groupId = 0)
{
  return Atom(position, 0.0, scalingVDW, 0.0, 0u, 0u, 0u, groupId, groupId != 0);
}

double totalEnergy(const ForceField& forceField, const SimulationBox& box,
                   const std::optional<InterpolationEnergyGrid>& grid, std::span<const Atom> atoms)
{
  RunningEnergy energy;
  Interactions::computeExternalFieldEnergy(true, forceField, box, atoms, energy, grid);
  return energy.externalFieldVDW;
}

std::optional<RunningEnergy> cbmcEnergy(const ForceField& forceField, const SimulationBox& box,
                                        const std::optional<InterpolationEnergyGrid>& grid, std::span<const Atom> atoms)
{
  return CBMC::computeExternalFieldEnergy(true, forceField, box, grid, 8.0, 8.0, atoms);
}
}  // namespace

// A rectangular channel with width != height: the total energy, the move energy difference, and the CBMC trial
// energy must agree on which positions are inside. The move difference and the CBMC energy used to transpose the
// width and the height for 'RectangleZ' relative to the total energy, so the sampled channel differed from the
// one the running total and the drift check assumed.
TEST(external_field, rectangular_channel_total_difference_and_cbmc_agree)
{
  const SimulationBox box(20.0, 20.0, 20.0);
  ForceField forceField = makeForceField(ForceField::PotentialEnergySurfaceType::RectangleZ);
  forceField.externalFieldGeometryParameters.x = 3.0;  // width
  forceField.externalFieldGeometryParameters.y = 6.0;  // height
  const std::optional<InterpolationEnergyGrid> noGrid = std::nullopt;
  const double sentinel = 10.0 * forceField.energyOverlapCriteria;

  const double3 centre(10.0, 10.0, 5.0);
  // (offset from the channel axis, inside?) for the documented convention: along z the width is the extent along
  // y and the height the extent along x
  const std::vector<std::pair<double3, bool>> probes = {
      {double3(0.0, 0.0, 0.0), true},   {double3(5.0, 0.0, 0.0), true},   {double3(5.5, 0.0, 7.0), true},
      {double3(6.5, 0.0, 0.0), false},  {double3(0.0, 2.9, 0.0), true},   {double3(0.0, 3.5, 0.0), false},
      {double3(-5.0, 2.0, 12.0), true}, {double3(-5.0, -3.5, 0.0), false}};

  for (const auto& [offset, inside] : probes)
  {
    const std::vector<Atom> atom{makeAtom(centre + offset)};

    EXPECT_DOUBLE_EQ(totalEnergy(forceField, box, noGrid, atom), inside ? 0.0 : sentinel)
        << "offset " << offset.x << " " << offset.y << " " << offset.z;

    // trial position: an overlap when outside
    std::optional<RunningEnergy> asNew =
        Interactions::computeExternalFieldEnergyDifference(true, forceField, box, noGrid, atom, {});
    EXPECT_EQ(asNew.has_value(), inside);
    if (asNew.has_value()) EXPECT_DOUBLE_EQ(asNew->externalFieldVDW, 0.0);

    // existing position: carries exactly the sentinel the total energy assigned to it
    std::optional<RunningEnergy> asOld =
        Interactions::computeExternalFieldEnergyDifference(true, forceField, box, noGrid, {}, atom);
    ASSERT_TRUE(asOld.has_value());
    EXPECT_DOUBLE_EQ(asOld->externalFieldVDW, inside ? 0.0 : -sentinel);

    EXPECT_EQ(cbmcEnergy(forceField, box, noGrid, atom).has_value(), inside);
  }
}

// A molecule that starts outside a confining cylinder (e.g. from a restart file) is not trapped there: the old
// position carries the sentinel energy, so the move that brings the molecule back in is accepted and the running
// energy stays consistent with the total energy. Trial positions outside remain an overlap.
TEST(external_field, molecule_outside_a_cylinder_can_move_back_in)
{
  const SimulationBox box(20.0, 20.0, 20.0);
  ForceField forceField = makeForceField(ForceField::PotentialEnergySurfaceType::CylinderZ);
  forceField.externalFieldGeometryParameters.x = 4.0;  // radius
  const std::optional<InterpolationEnergyGrid> noGrid = std::nullopt;
  const double sentinel = 10.0 * forceField.energyOverlapCriteria;

  const std::vector<Atom> outside{makeAtom(double3(16.0, 10.0, 3.0))};
  const std::vector<Atom> inside{makeAtom(double3(11.0, 10.0, 3.0))};
  const std::vector<Atom> alsoOutside{makeAtom(double3(10.0, 16.5, 3.0))};

  EXPECT_DOUBLE_EQ(totalEnergy(forceField, box, noGrid, outside), sentinel);
  EXPECT_DOUBLE_EQ(totalEnergy(forceField, box, noGrid, inside), 0.0);

  std::optional<RunningEnergy> backIn =
      Interactions::computeExternalFieldEnergyDifference(true, forceField, box, noGrid, inside, outside);
  ASSERT_TRUE(backIn.has_value());
  EXPECT_DOUBLE_EQ(backIn->externalFieldVDW, -sentinel);

  EXPECT_FALSE(Interactions::computeExternalFieldEnergyDifference(true, forceField, box, noGrid, alsoOutside, outside)
                   .has_value());
  EXPECT_FALSE(
      Interactions::computeExternalFieldEnergyDifference(true, forceField, box, noGrid, outside, inside).has_value());
  EXPECT_FALSE(cbmcEnergy(forceField, box, noGrid, outside).has_value());
}

// The confining geometry is periodic: an atom one cell (or several) away from the home cell, as a molecule that
// drifted or a chain bead overhanging the edge of the cell, sits in its own image of the channel. The total energy,
// the move difference and the CBMC energy tested the unwrapped position against the channel of the home cell and
// reported a hard-wall overlap.
TEST(external_field, confining_geometry_is_periodic)
{
  const SimulationBox box(20.0, 20.0, 20.0);
  const std::optional<InterpolationEnergyGrid> noGrid = std::nullopt;
  for (const ForceField::PotentialEnergySurfaceType type :
       {ForceField::PotentialEnergySurfaceType::CylinderZ, ForceField::PotentialEnergySurfaceType::RectangleX})
  {
    ForceField forceField = makeForceField(type);
    forceField.externalFieldGeometryParameters.x = 4.0;
    forceField.externalFieldGeometryParameters.y = 4.0;
    for (const double3& position :
         {double3(31.0, 10.5, 12.0), double3(-9.0, 9.0, -11.0), double3(11.0, 50.0, 30.5), double3(-30.0, -29.5, 9.0)})
    {
      const std::vector<Atom> atom{makeAtom(position)};
      EXPECT_DOUBLE_EQ(totalEnergy(forceField, box, noGrid, atom), 0.0) << position.x << " " << position.y;
      EXPECT_TRUE(
          Interactions::computeExternalFieldEnergyDifference(true, forceField, box, noGrid, atom, {}).has_value());
      EXPECT_TRUE(cbmcEnergy(forceField, box, noGrid, atom).has_value());
    }
    // and an image of a position outside the channel is outside as well
    const std::vector<Atom> outside{makeAtom(double3(36.5, 36.5, 30.0))};
    EXPECT_GT(totalEnergy(forceField, box, noGrid, outside), 0.0);
    EXPECT_FALSE(cbmcEnergy(forceField, box, noGrid, outside).has_value());
  }
}

// A tabulated field couples linearly to the fractional molecules, U = scalingVDW * U_ext, in the total energy,
// the move difference, and the CBMC trial energy alike, and all three report dU/d(scalingVDW) = U_ext for the
// dU/dlambda estimators (which used to be zero, and the old position used to be read unscaled).
TEST(external_field, tabulated_field_scales_with_lambda_and_carries_dudlambda)
{
  const SimulationBox box(10.0, 10.0, 10.0);
  ForceField forceField = makeForceField(ForceField::PotentialEnergySurfaceType::ThirdOrderPolynomialTestFunction);

  std::optional<InterpolationEnergyGrid> grid = InterpolationEnergyGrid(
      box, double3(0.0, 0.0, 0.0), uint3{16, 16, 16}, ForceFieldSettings::InterpolationScheme::Tricubic);
  std::stringstream stream;
  grid->makeExternalFieldInterpolationGrid(stream, forceField, box);

  const double3 position(3.3, 7.1, 4.6);
  const double fieldValue = grid->interpolate(position);
  ASSERT_GT(std::abs(fieldValue), 1.0);

  const std::uint8_t group = 1;
  const std::vector<Atom> integer{makeAtom(position)};
  const std::vector<Atom> fractionalLow{makeAtom(position, 0.3, group)};
  const std::vector<Atom> fractionalHigh{makeAtom(position, 0.7, group)};

  RunningEnergy total;
  Interactions::computeExternalFieldEnergy(true, forceField, box, fractionalLow, total, grid);
  EXPECT_NEAR(total.externalFieldVDW, 0.3 * fieldValue, 1e-10);
  EXPECT_NEAR(total.dudlambdaVDW[group - 1], fieldValue, 1e-10);
  EXPECT_NEAR(totalEnergy(forceField, box, grid, integer), fieldValue, 1e-10);

  // the lambda change of a fractional molecule at a fixed position
  std::optional<RunningEnergy> lambdaChange =
      Interactions::computeExternalFieldEnergyDifference(true, forceField, box, grid, fractionalHigh, fractionalLow);
  ASSERT_TRUE(lambdaChange.has_value());
  EXPECT_NEAR(lambdaChange->externalFieldVDW, 0.4 * fieldValue, 1e-10);
  EXPECT_NEAR(lambdaChange->dudlambdaVDW[group - 1], 0.0, 1e-10);

  // removing the fractional molecule: minus its scaled energy and minus its dU/dlambda factor
  std::optional<RunningEnergy> removal =
      Interactions::computeExternalFieldEnergyDifference(true, forceField, box, grid, {}, fractionalLow);
  ASSERT_TRUE(removal.has_value());
  EXPECT_NEAR(removal->externalFieldVDW, -0.3 * fieldValue, 1e-10);
  EXPECT_NEAR(removal->dudlambdaVDW[group - 1], -fieldValue, 1e-10);

  std::optional<RunningEnergy> trial = cbmcEnergy(forceField, box, grid, fractionalLow);
  ASSERT_TRUE(trial.has_value());
  EXPECT_NEAR(trial->externalFieldVDW, 0.3 * fieldValue, 1e-10);
  EXPECT_NEAR(trial->dudlambdaVDW[group - 1], fieldValue, 1e-10);
}

namespace
{
// Runs in an empty directory, since the grid routines read and write files relative to the working directory.
class ScopedTemporaryDirectory
{
 public:
  explicit ScopedTemporaryDirectory(const std::string& name)
      : original(std::filesystem::current_path()), path(std::filesystem::temp_directory_path() / name)
  {
    std::filesystem::remove_all(path);
    std::filesystem::create_directories(path);
    std::filesystem::current_path(path);
  }
  ~ScopedTemporaryDirectory()
  {
    std::filesystem::current_path(original);
    std::filesystem::remove_all(path);
  }

 private:
  std::filesystem::path original;
  std::filesystem::path path;
};

std::vector<double> readCubeValues(const std::string& fileName)
{
  std::ifstream stream(fileName);
  std::string line;
  // two comment lines, the origin, three grid vectors, and the single atom
  for (std::size_t i = 0; i < 7; ++i) std::getline(stream, line);
  std::vector<double> values;
  double value{};
  while (stream >> value) values.push_back(value);
  return values;
}
}  // namespace

// The analytic field, the CBMC energy and the interpolation grid are the same surface, shifted by the same origin.
TEST(external_field, analytic_cbmc_and_grid_agree_on_every_test_function)
{
  const SimulationBox box(4.0, 4.0, 4.0);
  const std::optional<InterpolationEnergyGrid> noGrid = std::nullopt;
  const std::vector<Atom> atom{makeAtom(double3(1.31, 2.13, 0.72))};

  for (const ForceField::PotentialEnergySurfaceType type :
       {ForceField::PotentialEnergySurfaceType::SecondOrderPolynomialTestFunction,
        ForceField::PotentialEnergySurfaceType::ThirdOrderPolynomialTestFunction,
        ForceField::PotentialEnergySurfaceType::FourthOrderPolynomialTestFunction,
        ForceField::PotentialEnergySurfaceType::FifthOrderPolynomialTestFunction,
        ForceField::PotentialEnergySurfaceType::SixthOrderPolynomialTestFunction,
        ForceField::PotentialEnergySurfaceType::ExponentialNonPolynomialTestFunction})
  {
    for (const double3 origin : {double3(0.0, 0.0, 0.0), double3(0.5, -0.25, 0.75)})
    {
      ForceField forceField = makeForceField(type);
      forceField.potentialEnergySurfaceOrigin = origin;

      std::optional<InterpolationEnergyGrid> grid = InterpolationEnergyGrid(
          box, origin, uint3{65, 65, 65}, ForceFieldSettings::InterpolationScheme::Triquintic);
      std::stringstream stream;
      grid->makeExternalFieldInterpolationGrid(stream, forceField, box);

      const double fromGrid = totalEnergy(forceField, box, grid, atom);
      const double tolerance = 1e-5 * std::max(1.0, std::abs(fromGrid));
      ASSERT_NE(fromGrid, 0.0) << "type " << std::to_underlying(type);
      EXPECT_NEAR(totalEnergy(forceField, box, noGrid, atom), fromGrid, tolerance)
          << "type " << std::to_underlying(type) << " origin " << origin.x;

      std::optional<RunningEnergy> trial = cbmcEnergy(forceField, box, noGrid, atom);
      ASSERT_TRUE(trial.has_value());
      EXPECT_NEAR(trial->externalFieldVDW, fromGrid, tolerance)
          << "type " << std::to_underlying(type) << " origin " << origin.x;
    }
  }
}

// One energy per grid vertex in cube order (x slowest, z fastest), not the Hermite derivatives stored with it.
TEST(external_field, cube_output_holds_one_energy_per_grid_vertex)
{
  const ScopedTemporaryDirectory directory("raspa3_external_field_cube_output");
  const SimulationBox box(8.0, 8.0, 8.0);
  ForceField forceField = makeForceField(ForceField::PotentialEnergySurfaceType::ThirdOrderPolynomialTestFunction);

  for (const ForceFieldSettings::InterpolationScheme scheme : {ForceFieldSettings::InterpolationScheme::Polynomial,
                                                               ForceFieldSettings::InterpolationScheme::Tricubic,
                                                               ForceFieldSettings::InterpolationScheme::Triquintic})
  {
    const uint3 points{3, 4, 5};
    InterpolationEnergyGrid grid(box, double3(0.0, 0.0, 0.0), points, scheme);
    std::stringstream stream;
    grid.makeExternalFieldInterpolationGrid(stream, forceField, box);
    grid.writeOutput(0, box, forceField);

    const std::vector<double> values = readCubeValues("interpolation_grids/grid_external_field.s0.cube");
    ASSERT_EQ(values.size(), points.x * points.y * points.z);

    std::size_t index{0};
    for (std::size_t i = 0; i < points.x; ++i)
    {
      for (std::size_t j = 0; j < points.y; ++j)
      {
        for (std::size_t k = 0; k < points.z; ++k)
        {
          const double x = 8.0 * static_cast<double>(i) / 2.0;
          const double y = 8.0 * static_cast<double>(j) / 3.0;
          const double z = 8.0 * static_cast<double>(k) / 4.0;
          const double expected = (x * x * x + 2.0 * x * x * y * z * z - y * y * z + 1.0) * Units::EnergyToKelvin;
          EXPECT_NEAR(values[index], expected, 1e-8 * std::max(1.0, std::abs(expected)));
          ++index;
        }
      }
    }
  }
}

// Asking for the grid output without a grid writes nothing; a system without a framework can build and write a grid.
TEST(external_field, grid_setup_without_framework_or_without_grid)
{
  const ScopedTemporaryDirectory directory("raspa3_external_field_grid_setup");
  const SimulationBox box(8.0, 8.0, 8.0);
  ForceField forceField = makeForceField(ForceField::PotentialEnergySurfaceType::ThirdOrderPolynomialTestFunction);
  forceField.settings.writeExternalFieldInterpolationGrid = true;
  forceField.settings.interpolationScheme = ForceFieldSettings::InterpolationScheme::Triquintic;
  forceField.settings.numberOfGridTestPoints = 10;
  const std::string output = "interpolation_grids/grid_external_field.s0.cube";

  std::ostringstream log;
  forceField.settings.useExternalFieldGrid = false;
  System withoutGrid = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {}, {}, {}, 5);
  withoutGrid.hasExternalField = true;
  withoutGrid.createExternalFieldInterpolationGrid(log, 0);
  EXPECT_FALSE(withoutGrid.externalFieldInterpolationGrid.has_value());
  EXPECT_FALSE(std::filesystem::exists(output));

  // an analytic field is tabulated on the requested grid (a field read from a cube file takes the file's
  // dimensions; this test used to expect the file's dimensions for the analytic field too)
  forceField.settings.numberOfExternalFieldGridPoints = uint3{5, 5, 5};
  forceField.settings.useExternalFieldGrid = true;
  System withGrid = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {}, {}, {}, 5);
  withGrid.hasExternalField = true;
  withGrid.createExternalFieldInterpolationGrid(log, 0);
  ASSERT_TRUE(withGrid.externalFieldInterpolationGrid.has_value());
  EXPECT_EQ(readCubeValues(output).size(), 125uz);
}

// A Gaussian CUBE file gives its step vectors in Bohr when the voxel counts are positive (what Gaussian, CP2K, ORCA,
// Multiwfn and Psi4 write) and in Angstrom when they are negative. The reader took only the negative form and
// refused every file in the default, positive one; a mix of the two signs is still malformed.
TEST(external_field, cube_files_with_bohr_or_angstrom_steps_are_read)
{
  const ScopedTemporaryDirectory directory("raspa3_external_field_cube_units");
  const auto write = [](const std::string& name, int sign, int signZ)
  {
    std::ofstream file(name);
    std::print(file, "comment\ncomment\n    1    0.000000    0.000000    0.000000\n");
    std::print(file, "{:5d}    1.000000    0.000000    0.000000\n", sign * 2);
    std::print(file, "{:5d}    0.000000    1.000000    0.000000\n", sign * 3);
    std::print(file, "{:5d}    0.000000    0.000000    1.000000\n", signZ * 4);
    std::print(file, "    8    0.000000    0.000000    0.000000    0.000000\n");
    for (int i = 0; i < 24; ++i) std::print(file, "{} ", static_cast<double>(i));
    std::print(file, "\n");
  };
  for (const int sign : {1, -1})
  {
    const std::string name = sign > 0 ? "bohr.cube" : "angstrom.cube";
    write(name, sign, sign);
    const uint3 dimensions = InterpolationEnergyGrid::parseExternalFieldGridDimensions(name);
    EXPECT_EQ(dimensions.x, 2u);
    EXPECT_EQ(dimensions.y, 3u);
    EXPECT_EQ(dimensions.z, 4u);
    EXPECT_EQ(InterpolationEnergyGrid::parseExternalFieldGridCube(name).size(), 24uz) << name;
  }
  write("mixed.cube", 1, -1);
  EXPECT_THROW(InterpolationEnergyGrid::parseExternalFieldGridDimensions("mixed.cube"), std::runtime_error);
}

// The energy status of the pressure computation (the source of the averaged energies) carries the external-field
// energy of each component; before, it was left at zero and missing from the total.
TEST(external_field, energy_status_holds_the_external_field_energy)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::SecondOrderPolynomialTestFunction;
  Component methane = Component::makeMethane(forceField, 0);
  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), true, 300.0, 1e4, 1.0, {}, {methane}, {}, {3}, 5);
  for (std::size_t i = 0; i < 3; ++i)
  {
    const double3 position(4.0 + 6.0 * static_cast<double>(i), 7.0, 11.0);
    system.moleculeData[i].centerOfMassPosition = position;
    system.spanOfMoleculeAtoms()[i].position = position;
  }

  const RunningEnergy running = system.computeTotalEnergies();
  ASSERT_GT(std::abs(running.externalFieldVDW), 1.0);
  EnergyStatus status = system.computeMolecularPressure().first;
  EXPECT_NEAR(status.externalFieldComponentEnergy(0, 0).VanDerWaals.energy, running.externalFieldVDW,
              1e-9 * std::abs(running.externalFieldVDW));
  EXPECT_NEAR(status.totalEnergy.energy, running.potentialEnergy(), 1e-6 * std::abs(running.potentialEnergy()));
}

// A confining geometry is not tabulated, even with the default 'UseExternalFieldGrid': the grid of a hard wall was
// zero everywhere and, read first by the energy routines, let molecules sit outside the channel.
TEST(external_field, confining_geometry_is_not_tabulated)
{
  const SimulationBox box(20.0, 20.0, 20.0);
  ForceField forceField = makeForceField(ForceField::PotentialEnergySurfaceType::CylinderZ);
  forceField.externalFieldGeometryParameters.x = 4.0;
  forceField.settings.numberOfExternalFieldGridPoints = uint3{5, 5, 5};
  ASSERT_TRUE(forceField.settings.useExternalFieldGrid);
  System system = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {}, {}, {}, 5);
  system.hasExternalField = true;
  std::ostringstream log;
  system.createExternalFieldInterpolationGrid(log, 0);
  EXPECT_FALSE(system.externalFieldInterpolationGrid.has_value());

  const std::vector<Atom> outside{makeAtom(double3(16.0, 10.0, 3.0))};
  EXPECT_GT(totalEnergy(forceField, box, system.externalFieldInterpolationGrid, outside), 0.0);
}

// A grid of one cell along an axis: the six-point polynomial stencil reaches two cells beyond it on either side, and
// adding the period once left the index at -1 (read as SIZE_MAX, out of bounds). Wrapped modulo the period, a
// constant field interpolates to that constant.
TEST(external_field, polynomial_stencil_wraps_a_grid_of_one_cell)
{
  const SimulationBox box(10.0, 10.0, 10.0);
  InterpolationEnergyGrid grid(box, double3(), uint3(2, 2, 2), ForceFieldSettings::InterpolationScheme::Polynomial);
  grid.data.assign(8, 5.0);
  for (const double3& position : {double3(0.0, 0.0, 0.0), double3(3.0, 7.0, 1.0), double3(9.9, 0.1, 5.0)})
  {
    EXPECT_NEAR(grid.interpolate(position), 5.0, 1e-12) << position.x << " " << position.y << " " << position.z;
  }
  EXPECT_THROW(
      InterpolationEnergyGrid(box, double3(), uint3(1, 4, 4), ForceFieldSettings::InterpolationScheme::Polynomial),
      std::runtime_error);
}

// A periodic CUBE file (as CP2K writes it) holds N voxels at s = i/N a step L/N apart, without the copy of s = 0 at
// s = 1. The grid used to take the N planes as spanning the cell end to end, so its nodes sat at s = i/(N-1) and the
// field was stretched over the cell. RASPA's own files (N points a step L/(N-1) apart) keep their counts.
TEST(external_field, periodic_cube_file_gets_the_plane_at_the_cell_boundary)
{
  const ScopedTemporaryDirectory directory("raspa3_external_field_periodic_cube");
  const SimulationBox box(10.0, 10.0, 10.0);
  const auto write = [](const std::string& name, double step)
  {
    std::ofstream file(name);
    std::print(file, "comment\ncomment\n    0    0.000000    0.000000    0.000000\n");
    std::print(file, "   -4    {}    0.000000    0.000000\n", step);
    std::print(file, "   -4    0.000000    {}    0.000000\n", step);
    std::print(file, "   -4    0.000000    0.000000    {}\n", step);
    // x slowest, z fastest; the value is that of the x plane
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 16; ++j) std::print(file, "{} ", 10.0 * static_cast<double>(i));
    std::print(file, "\n");
  };

  write("periodic.cube", 2.5);
  const uint3 dimensions = InterpolationEnergyGrid::externalFieldGridDimensions("periodic.cube", box);
  EXPECT_EQ(dimensions.x, 5u);
  EXPECT_EQ(dimensions.y, 5u);
  EXPECT_EQ(dimensions.z, 5u);

  ForceField forceField = makeForceField(ForceField::PotentialEnergySurfaceType::GridFile);
  forceField.externalFieldGridFileName = "periodic.cube";
  InterpolationEnergyGrid grid(box, double3(), dimensions, ForceFieldSettings::InterpolationScheme::Polynomial);
  std::ostringstream log;
  grid.makeExternalFieldInterpolationGrid(log, forceField, box);
  for (int i = 0; i < 4; ++i)
  {
    const double3 node(2.5 * static_cast<double>(i), 5.0, 5.0);
    EXPECT_NEAR(grid.interpolate(node), 10.0 * static_cast<double>(i) * Units::KelvinToEnergy, 1e-9) << i;
  }

  write("end_to_end.cube", 10.0 / 3.0);
  const uint3 own = InterpolationEnergyGrid::externalFieldGridDimensions("end_to_end.cube", box);
  EXPECT_EQ(own.x, 4u);
  EXPECT_EQ(own.z, 4u);
}
