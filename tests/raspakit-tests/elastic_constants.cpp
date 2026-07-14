#include <gtest/gtest.h>

import std;

import atom;
import bond_potential;
import component;
import double3;
import elastic_constants;
import forcefield;
import framework;
import generalized_hessian;
import int3;
import minimization;
import minimization_cell_layout;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import minimization_generalized_coordinates;
import minimization_options;
import simulationbox;
import system;
import units;

namespace
{
System makeHarmonicDimer(double pressure = 0.0)
{
  ForceField forceField = ForceField::makeZeoliteForceField(10.0, true, false, true);
  const SimulationBox box(10.0, 10.0, 10.0);
  const Atom atomA({0.44, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true);
  const Atom atomB({0.56, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "harmonic-dimer", box, 1, {atomA, atomB}, {atomA, atomB}, int3(1, 1, 1));
  framework.rigid = false;
  framework.intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {2000.0, 1.2})};
  framework.intraMolecularImageShifts.bonds = {{{int3(0, 0, 0), int3(0, 0, 0)}}};

  System system(forceField, box, false, 300.0, pressure, 1.0, framework, {}, {}, {}, 5);
  system.spanOfFrameworkAtoms()[0].position = {4.4, 5.0, 5.0};
  system.spanOfFrameworkAtoms()[1].position = {5.6, 5.0, 5.0};
  return system;
}

System makeHarmonicTetrahedron()
{
  ForceField forceField = ForceField::makeZeoliteForceField(10.0, true, false, true);
  const SimulationBox box(12.0, 12.0, 12.0);
  constexpr double edge = 1.2;
  const std::vector<double3> positions = {
      {4.0, 4.0, 4.0},
      {4.0 + edge, 4.0, 4.0},
      {4.0 + 0.5 * edge, 4.0 + 0.5 * std::sqrt(3.0) * edge, 4.0},
      {4.0 + 0.5 * edge, 4.0 + std::sqrt(3.0) * edge / 6.0, 4.0 + std::sqrt(2.0 / 3.0) * edge}};
  std::vector<Atom> atoms;
  for (const double3& position : positions)
  {
    atoms.emplace_back(position, 0.0, 1.0, 0, 0, 0, 0, true);
  }
  Framework framework(forceField, "harmonic-tetrahedron", box, 1, atoms, atoms, int3(1, 1, 1));
  framework.rigid = false;
  constexpr std::array<std::array<std::size_t, 2>, 6> pairs = {
      std::array<std::size_t, 2>{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  for (const auto& pair : pairs)
  {
    framework.intraMolecularPotentials.bonds.emplace_back(pair, BondType::Harmonic, std::vector<double>{2000.0, edge});
    framework.intraMolecularImageShifts.bonds.push_back({int3(0, 0, 0), int3(0, 0, 0)});
  }

  System system(forceField, box, false, 300.0, 0.0, 1.0, framework, {}, {}, {}, 5);
  for (std::size_t atom = 0; atom < positions.size(); ++atom)
  {
    system.spanOfFrameworkAtoms()[atom].position = positions[atom];
  }
  return system;
}

double energyAtLogarithmicStrain(System system, const std::array<double, 6>& strain)
{
  system.cellMinimizationType = CellMinimizationType::Regular;
  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, system.spanOfFrameworkAtoms().size(), 6);
  std::vector<double> displacement(layout.numDofs(), 0.0);
  for (std::size_t coordinate = 0; coordinate < strain.size(); ++coordinate)
  {
    displacement[*layout.cellDof(coordinate)] = strain[coordinate];
  }
  applyGeneralizedDisplacement(system, layout, displacement);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs(), 0.0);
  DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = true};
  DerivativeResults derivatives{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout, capabilities, derivatives);
  return derivatives.energy;
}
}  // namespace

TEST(elastic_constants, relaxed_harmonic_dimer_removes_clamped_axial_stiffness)
{
  System system = makeHarmonicDimer();
  const ElasticConstantsResult result = computeElasticConstants(system);

  EXPECT_GT(result.born[0], 0.0);
  EXPECT_NEAR(result.born[0], 2000.0 * Units::KelvinToEnergy * 1.2 * 1.2 / 1000.0, 1.0e-12);
  EXPECT_GT(result.relaxation[0], 0.0);
  EXPECT_NEAR(result.stiffness[0], 0.0, 1.0e-9 * result.born[0]);
  EXPECT_FALSE(result.complianceAvailable);
  EXPECT_GE(result.discardedInternalModes, 3u);
  for (std::size_t row = 0; row < 6; ++row)
  {
    for (std::size_t column = 0; column < 6; ++column)
    {
      EXPECT_NEAR(result.stiffness[row * 6 + column], result.stiffness[column * 6 + row], 1.0e-12);
    }
  }
}

TEST(elastic_constants, applies_raspa2_hydrostatic_pressure_correction_once)
{
  constexpr double pressurePa = 2.0e5;
  System zeroPressureSystem = makeHarmonicDimer();
  System pressuredSystem = makeHarmonicDimer(pressurePa);
  const ElasticConstantsResult zeroPressure = computeElasticConstants(zeroPressureSystem);
  const ElasticConstantsResult pressured = computeElasticConstants(pressuredSystem);
  const double pressureInternal = pressuredSystem.pressure;

  for (std::size_t direction = 0; direction < 3; ++direction)
  {
    EXPECT_DOUBLE_EQ(pressured.pressureCorrection[direction * 6 + direction], -pressureInternal);
  }
  for (std::size_t shear = 3; shear < 6; ++shear)
  {
    EXPECT_DOUBLE_EQ(pressured.pressureCorrection[shear * 6 + shear], -0.5 * pressureInternal);
  }
  for (std::size_t index = 0; index < 36; ++index)
  {
    EXPECT_NEAR(pressured.born[index], zeroPressure.born[index],
                1.0e-12 * std::max(1.0, std::abs(zeroPressure.born[index])));
    EXPECT_NEAR(pressured.relaxation[index], zeroPressure.relaxation[index],
                1.0e-12 * std::max(1.0, std::abs(zeroPressure.relaxation[index])));
  }
}

TEST(elastic_constants, minimization_driver_computes_requested_postprocessing)
{
  MinimizationOptions options{};
  options.maximumNumberOfSteps = 2;
  options.computeElasticConstants = true;
  Minimization minimization(options, {makeHarmonicDimer()}, false);

  minimization.setup();
  minimization.runPhase();

  ASSERT_EQ(minimization.results.size(), 1u);
  EXPECT_TRUE(minimization.results[0].converged);
  ASSERT_TRUE(minimization.results[0].elasticConstants.has_value());
  EXPECT_GT(minimization.results[0].elasticConstants->born[0], 0.0);
}

TEST(elastic_constants, all_six_born_strains_match_energy_finite_differences)
{
  System system = makeHarmonicTetrahedron();
  const ElasticConstantsResult result = computeElasticConstants(system);
  constexpr std::array<std::size_t, 6> voigtToRegular = {0, 3, 5, 4, 2, 1};
  constexpr double delta = 2.0e-5;
  const double referenceEnergy = energyAtLogarithmicStrain(system, {});
  const double volume = system.simulationBox.volume;

  for (std::size_t row = 0; row < 6; ++row)
  {
    for (std::size_t column = row; column < 6; ++column)
    {
      double numerical{};
      if (row == column)
      {
        std::array<double, 6> plus{};
        std::array<double, 6> minus{};
        plus[voigtToRegular[row]] = delta;
        minus[voigtToRegular[row]] = -delta;
        numerical = (energyAtLogarithmicStrain(system, plus) - 2.0 * referenceEnergy +
                     energyAtLogarithmicStrain(system, minus)) /
                    (delta * delta * volume);
      }
      else
      {
        std::array<double, 6> pp{};
        std::array<double, 6> pm{};
        std::array<double, 6> mp{};
        std::array<double, 6> mm{};
        pp[voigtToRegular[row]] = pp[voigtToRegular[column]] = delta;
        pm[voigtToRegular[row]] = delta;
        pm[voigtToRegular[column]] = -delta;
        mp[voigtToRegular[row]] = -delta;
        mp[voigtToRegular[column]] = delta;
        mm[voigtToRegular[row]] = mm[voigtToRegular[column]] = -delta;
        numerical = (energyAtLogarithmicStrain(system, pp) - energyAtLogarithmicStrain(system, pm) -
                     energyAtLogarithmicStrain(system, mp) + energyAtLogarithmicStrain(system, mm)) /
                    (4.0 * delta * delta * volume);
      }
      const double analytic = result.born[row * 6 + column];
      const double scale = std::max({1.0, std::abs(analytic), std::abs(numerical)});
      EXPECT_NEAR(analytic, numerical, 2.0e-5 * scale) << "row=" << row << " column=" << column;
    }
  }
}
