#include <gtest/gtest.h>

import std;

import atom;
import archive;
import bond_potential;
import component;
import double3;
import double3x3;
import elastic_constants;
import thermobarostat;
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
import property_elastic_constants_fluctuation;
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

// a rest length other than the edge leaves the tetrahedron under stress
System makeHarmonicTetrahedron(std::optional<double> restLength = std::nullopt)
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
    framework.intraMolecularPotentials.bonds.emplace_back(pair, BondType::Harmonic,
                                                          std::vector<double>{2000.0, restLength.value_or(edge)});
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

// The Barron-Klein stress-strain coefficients at the hydrostatic stress -P from the infinitesimal-strain tensor: +P
// between different normal strains, nothing on the normal diagonal and -P/2 on the shear diagonal (the RASPA2
// correction put -P on the normal diagonal and nothing off it, leaving every normal entry P too low).
TEST(elastic_constants, applies_the_barron_klein_pressure_correction_once)
{
  constexpr double pressurePa = 2.0e5;
  System zeroPressureSystem = makeHarmonicDimer();
  System pressuredSystem = makeHarmonicDimer(pressurePa);
  const ElasticConstantsResult zeroPressure = computeElasticConstants(zeroPressureSystem);
  const ElasticConstantsResult pressured = computeElasticConstants(pressuredSystem);
  const double pressureInternal = pressuredSystem.pressure;

  for (std::size_t row = 0; row < 3; ++row)
  {
    for (std::size_t column = 0; column < 3; ++column)
    {
      EXPECT_DOUBLE_EQ(pressured.pressureCorrection[row * 6 + column], row == column ? 0.0 : pressureInternal);
    }
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

TEST(elastic_constants, affine_born_helper_matches_static_born_term)
{
  const System system = makeHarmonicTetrahedron();
  const std::array<double, 36> affine = computeAffineBornTensor(system);
  const ElasticConstantsResult staticResult = computeElasticConstants(system);
  for (std::size_t i = 0; i < affine.size(); ++i)
    EXPECT_NEAR(affine[i], staticResult.born[i], 1.0e-12 * std::max(1.0, std::abs(staticResult.born[i])));
}

TEST(elastic_constants_fluctuation, synthetic_samples_recover_covariance_and_kinetic_terms)
{
  std::array<double, 6> plus{};
  std::array<double, 6> minus{};
  plus[0] = 1.0;
  minus[0] = -1.0;
  std::array<double, 36> born{};
  born[0] = 10.0;

  ElasticFluctuationTerms terms(plus, {}, born, 2.0, 3.0, 4.0);
  terms += ElasticFluctuationTerms(minus, {}, born, 2.0, 3.0, 4.0);
  const ElasticFluctuationData result = (terms / 2.0).compositeProperty();

  EXPECT_DOUBLE_EQ(result.fluctuation[0], -6.0);
  EXPECT_NEAR(result.kinetic[0], 4.0 / 3.0, 1.0e-14);
  EXPECT_NEAR(result.kinetic[3 * 6 + 3], 2.0 / 3.0, 1.0e-14);
  EXPECT_NEAR(result.stiffness[0], 10.0 - 6.0 + 4.0 / 3.0, 1.0e-14);
  EXPECT_DOUBLE_EQ(result.stiffness[1], result.stiffness[6]);
}

TEST(elastic_constants_fluctuation, harmonic_canonical_quadrature_recovers_known_stiffness)
{
  constexpr double bornStiffness = 10.0;
  constexpr double expectedIsothermalStiffness = 5.0;
  const double stressAmplitude = std::sqrt(bornStiffness - expectedIsothermalStiffness);
  std::array<double, 6> plus{};
  std::array<double, 6> minus{};
  plus[0] = stressAmplitude;
  minus[0] = -stressAmplitude;
  std::array<double, 36> born{};
  born[0] = bornStiffness;

  ElasticFluctuationTerms terms(plus, {}, born, 1.0, 1.0, 0.0);
  terms += ElasticFluctuationTerms(minus, {}, born, 1.0, 1.0, 0.0);
  const ElasticFluctuationData result = (terms / 2.0).compositeProperty();

  EXPECT_NEAR(result.stiffness[0], expectedIsothermalStiffness, 1.0e-14);
}

TEST(elastic_constants_fluctuation, kinetic_virial_uses_flexible_framework_velocities)
{
  System system = makeHarmonicDimer();
  ASSERT_EQ(system.spanOfFrameworkDynamics().size(), 2u);
  system.forceField.pseudoAtoms[0].mass = 1.0;
  system.spanOfFrameworkDynamics()[0].velocity = {1.0, 2.0, 3.0};
  system.spanOfFrameworkDynamics()[1].velocity = {};
  const double3x3 kinetic = computeMolecularKineticVirial(system, BarostatCoupling::Atomic);
  EXPECT_DOUBLE_EQ(kinetic.ax, 1.0);
  EXPECT_DOUBLE_EQ(kinetic.by, 4.0);
  EXPECT_DOUBLE_EQ(kinetic.cz, 9.0);
  EXPECT_DOUBLE_EQ(kinetic.ay, 2.0);
  EXPECT_DOUBLE_EQ(kinetic.bz, 6.0);
}

TEST(elastic_constants_fluctuation, archive_round_trip_preserves_accumulator)
{
  PropertyElasticConstantsFluctuation original(3);
  std::array<double, 6> stress{};
  stress[0] = 2.0;
  std::array<double, 36> born{};
  born[0] = 7.0;
  original.addSample(1, ElasticFluctuationTerms(stress, {}, born, 1.5, 8.0, 3.0), 1.0);
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "raspa3_elastic_fluctuation_round_trip.bin";
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  PropertyElasticConstantsFluctuation restored;
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);
  ASSERT_EQ(restored.numberOfBlocks, 3u);
  EXPECT_DOUBLE_EQ(restored.bookKeeping[1].first.configurationalStress[0], 2.0);
  EXPECT_DOUBLE_EQ(restored.bookKeeping[1].first.born[0], 7.0);
}

namespace
{
System makeMixedFrameworkSystem()
{
  ForceField forceField({{"C", true, 12.0, 0.0, 0.0, 6, true}}, {{80.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot,
                        5.0, 5.0, 5.0, true, false, false);
  const SimulationBox box(14.0, 14.0, 14.0);
  std::vector<Atom> atoms;
  for (const double3& position : {double3(4.0, 7.0, 7.0), double3(5.4, 7.3, 7.0), double3(6.8, 6.8, 7.0)})
  {
    atoms.push_back(Atom(position / 14.0, 0.0, 1.0, 0, 0, 0, 0, true));
  }
  Framework framework(forceField, "mixed-chain", box, 1, atoms, atoms, int3(1, 1, 1));

  const std::filesystem::path definition = std::filesystem::temp_directory_path() / "raspa3-elastic-mixed-framework.json";
  std::ofstream(definition) << R"({"Type": "Flexible",
                                   "Groups": [{"Type": "Fixed", "Atoms": [0]}, {"Type": "Flexible", "Atoms": [1, 2]}],
                                   "Intra14VanDerWaalsScalingValue": 0.0,
                                   "Intra14ChargeChargeScalingValue": 0.0})";
  framework.readFrameworkDefinition(forceField, definition.string());
  std::filesystem::remove(definition);
  EXPECT_TRUE(framework.isMixed());
  return System(forceField, box, false, 300.0, 1.0e4, 1.0, framework, {}, {}, {}, 5);
}

template <typename Function>
void expectUnsupportedMixedFramework(Function&& function)
{
  try
  {
    function();
    FAIL() << "a mixed framework must be rejected";
  }
  catch (const std::runtime_error& error)
  {
    EXPECT_NE(std::string_view(error.what()).find("Elastic constants"), std::string_view::npos) << error.what();
  }
}
}  // namespace

// The cell derivatives behind the elastic constants are not available for a mixed framework: a clear error, not an
// inconsistent degree-of-freedom layout.
TEST(elastic_constants, mixed_framework_is_rejected)
{
  const System system = makeMixedFrameworkSystem();
  expectUnsupportedMixedFramework([&] { (void)computeElasticConstants(system, 1.0e-8); });
  expectUnsupportedMixedFramework([&] { (void)computeAffineBornTensor(system); });
}

// The Green-Lagrange Born term of the stress-fluctuation formula is the second derivative of the energy in the strain
// eta = (F^T F - I) / 2 (here F = diag(sqrt(1 + 2 eta))); under stress it differs from the infinitesimal-strain term.
TEST(elastic_constants, green_lagrange_born_matches_energy_finite_differences)
{
  const System system = makeHarmonicTetrahedron(1.35);
  const std::array<double, 36> greenLagrange = computeGreenLagrangeBornTensor(system);
  const std::array<double, 36> infinitesimal = computeAffineBornTensor(system);
  ASSERT_GT(std::abs(greenLagrange[0] - infinitesimal[0]), 1.0e-3 * std::abs(infinitesimal[0]));

  const double volume = system.simulationBox.volume;
  const auto energyAt = [&](double etaX, double etaY)
  {
    std::array<double, 6> logStrain{};
    logStrain[0] = 0.5 * std::log(1.0 + 2.0 * etaX);  // regular coordinate of xx
    logStrain[3] = 0.5 * std::log(1.0 + 2.0 * etaY);  // regular coordinate of yy
    return energyAtLogarithmicStrain(system, logStrain);
  };
  constexpr double delta = 1.0e-4;
  const double xx =
      (energyAt(delta, 0.0) - 2.0 * energyAt(0.0, 0.0) + energyAt(-delta, 0.0)) / (delta * delta * volume);
  const double xy =
      (energyAt(delta, delta) - energyAt(delta, -delta) - energyAt(-delta, delta) + energyAt(-delta, -delta)) /
      (4.0 * delta * delta * volume);
  EXPECT_NEAR(greenLagrange[0], xx, 1.0e-4 * std::max(1.0, std::abs(xx)));
  EXPECT_NEAR(greenLagrange[1], xy, 1.0e-4 * std::max(1.0, std::abs(xy)));
}

// An ideal gas: no Born or fluctuation term, the kinetic Green-Lagrange term 2P on the normal and P on the shear
// diagonal, and the Barron-Klein coefficients P on every normal entry and no shear stiffness -- bulk modulus P,
// shear modulus zero (the former correction gave P/3 and P/2).
TEST(elastic_constants_fluctuation, ideal_gas_has_bulk_modulus_p_and_no_shear_stiffness)
{
  constexpr double beta = 2.0;
  constexpr double volume = 3.0;
  constexpr double entities = 4.0;
  const double pressure = entities / (beta * volume);
  const std::array<double, 6> kineticStress = {pressure, pressure, pressure, 0.0, 0.0, 0.0};
  const ElasticFluctuationTerms terms({}, kineticStress, {}, beta, volume, entities);
  const ElasticFluctuationData result = terms.compositeProperty();
  for (std::size_t row = 0; row < 3; ++row)
  {
    for (std::size_t column = 0; column < 3; ++column)
    {
      EXPECT_NEAR(result.tangentStiffness[row * 6 + column], pressure, 1.0e-14) << row << " " << column;
    }
  }
  for (std::size_t shear = 3; shear < 6; ++shear)
  {
    EXPECT_NEAR(result.tangentStiffness[shear * 6 + shear], 0.0, 1.0e-14) << shear;
  }
}
