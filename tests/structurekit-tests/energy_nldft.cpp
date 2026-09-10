#include <gtest/gtest.h>

import std;

import uint3;
import unit_cell;
import energy_shared_bet_surface_area;
import energy_shared_nldft;

// The classical DFT isotherm on the energy grid: bulk coexistence, the uniform limit, the Henry limit and
// pore condensation, each on a grid small enough to run in seconds.

// The model N2 fluid at 77 K is subcritical with a liquid density near the real 0.0174 /Å ³ and a
// saturation pressure within an order of magnitude of the experimental atmosphere.
TEST(energy_nldft, bulk_coexistence_is_liquid_nitrogen_like)
{
  NLDFTOptions options;
  NLDFTBulk bulk = nldftBulkCoexistence(options);

  EXPECT_GT(bulk.saturationPressure, 1.0e4);
  EXPECT_LT(bulk.saturationPressure, 1.0e6);
  EXPECT_GT(bulk.liquidDensity, 0.012);
  EXPECT_LT(bulk.liquidDensity, 0.030);
  EXPECT_LT(bulk.gasDensity, 1.0e-3);
  EXPECT_LT(bulk.meanFieldIntegral, 0.0);
}


// An empty periodic box: the converged profile is the uniform bulk gas, which checks the whole FFT
// machinery at once --- if any weight or the attraction were wrong, the uniform state would not be
// self-consistent with the same bulk equation of state that set the reservoir.
TEST(energy_nldft, an_empty_box_holds_the_bulk_gas)
{
  NLDFTOptions options;
  options.pressurePoints = 2;
  options.xLow = 0.1;
  options.xHigh = 0.3;

  uint3 gridSize{32, 32, 32};
  UnitCell cell(20.0, 20.0, 20.0);
  std::vector<float> energy(32 * 32 * 32, 0.0f);

  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options);
  ASSERT_EQ(result.isotherm.size(), 2u);
  EXPECT_EQ(result.unconverged, 0u);

  for (const IsothermPoint &point : result.isotherm)
  {
    double expected = nldftGasDensity(options, result.bulk, point.relativePressure) * cell.volume;
    EXPECT_NEAR(point.moleculesPerCell, expected, 0.02 * expected);
  }
}


// The x -> 0 limit is the exact Boltzmann average of the field: rho_b integral exp(-U/kT).
TEST(energy_nldft, the_henry_limit_is_the_boltzmann_average)
{
  NLDFTOptions options;
  options.pressurePoints = 1;
  options.xLow = 1.0e-5;
  options.xHigh = 1.0e-5;

  uint3 gridSize{32, 32, 32};
  UnitCell cell(20.0, 20.0, 20.0);
  const std::size_t N = 32 * 32 * 32;
  std::vector<float> energy(N, 0.0f);
  for (std::size_t i = 0; i < N / 4; ++i) energy[i] = -150.0f;  // a shallow attractive quarter

  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options);
  ASSERT_EQ(result.isotherm.size(), 1u);

  double rhoBulk = nldftGasDensity(options, result.bulk, 1.0e-5);
  double boltzmann = 0.0;
  for (float u : energy) boltzmann += std::exp(-static_cast<double>(u) / options.temperature);
  double expected = rhoBulk * boltzmann * cell.volume / static_cast<double>(N);
  EXPECT_NEAR(result.isotherm.front().moleculesPerCell, expected, 0.05 * expected);
}


// A deep half-box condenses to a liquid-like density well below x = 1 while the empty half stays gaseous:
// pore filling emerges from the field, with no site model underneath.
TEST(energy_nldft, a_deep_region_condenses_before_saturation)
{
  NLDFTOptions options;
  options.pressurePoints = 3;
  options.xLow = 0.05;
  options.xHigh = 0.3;

  uint3 gridSize{32, 32, 32};
  UnitCell cell(20.0, 20.0, 20.0);
  const std::size_t N = 32 * 32 * 32;
  std::vector<float> energy(N, 0.0f);
  for (std::size_t i = 0; i < N / 2; ++i) energy[i] = -700.0f;

  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options);
  ASSERT_GE(result.isotherm.size(), 3u);

  double halfVolume = 0.5 * cell.volume;
  double liquidLike = result.bulk.liquidDensity * halfVolume;
  EXPECT_GT(result.isotherm.back().moleculesPerCell, 0.5 * liquidLike);
  EXPECT_LT(result.isotherm.back().moleculesPerCell, 2.5 * result.bulk.liquidDensity * cell.volume);

  double previous = 0.0;
  for (const IsothermPoint &point : result.isotherm)
  {
    EXPECT_GE(point.moleculesPerCell, previous * (1.0 - 1.0e-6));
    previous = point.moleculesPerCell;
  }
}


// An empty box seeded at the bulk liquid, held at the saturation chemical potential, is a stationary
// point of the functional: the FFT weights, the White Bear excess and the mean-field tail must all
// conspire to match Carnahan-Starling plus a rho.
TEST(energy_nldft, an_empty_box_holds_the_bulk_liquid)
{
  NLDFTOptions options;
  options.pressurePoints = 1;
  options.xLow = 1.0;
  options.xHigh = 1.0;
  options.seedLiquid = true;

  uint3 gridSize{32, 32, 32};
  UnitCell cell(20.0, 20.0, 20.0);
  std::vector<float> energy(32 * 32 * 32, 0.0f);

  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options);
  ASSERT_EQ(result.isotherm.size(), 1u);
  EXPECT_EQ(result.unconverged, 0u);
  double expected = result.bulk.liquidDensity * cell.volume;
  EXPECT_NEAR(result.isotherm.front().moleculesPerCell, expected, 0.05 * expected);
}


// A deep slab (a FER-like slit, a quarter of the cell) fills to a liquid-like film well below x = 1
// without overpacking the whole cell, and the isotherm is monotone.
TEST(energy_nldft, a_slit_fills_to_a_liquid_film_without_overpacking)
{
  NLDFTOptions options;
  options.pressurePoints = 6;
  options.xLow = 0.01;
  options.xHigh = 0.3;

  uint3 gridSize{32, 32, 32};
  UnitCell cell(20.0, 20.0, 20.0);
  const std::size_t N = 32 * 32 * 32;
  std::vector<float> energy(N, 0.0f);
  for (std::size_t k = 0; k < 8; ++k)
  {
    for (std::size_t j = 0; j < 32; ++j)
    {
      for (std::size_t i = 0; i < 32; ++i) energy[(k * 32 + j) * 32 + i] = -700.0f;
    }
  }

  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options);
  ASSERT_FALSE(result.isotherm.empty());

  double wellVolume = 0.25 * cell.volume;
  double film = result.bulk.liquidDensity * wellVolume;
  EXPECT_GT(result.isotherm.back().moleculesPerCell, 0.3 * film);
  EXPECT_LT(result.isotherm.back().moleculesPerCell, 2.0 * result.bulk.liquidDensity * cell.volume);

  double previous = 0.0;
  for (const IsothermPoint &point : result.isotherm)
  {
    EXPECT_GE(point.moleculesPerCell, previous * (1.0 - 1.0e-6));
    previous = point.moleculesPerCell;
  }
}


// Wells deeper than the floor are clipped so a Coulomb hole cannot set the Henry coefficient.
TEST(energy_nldft, wells_deeper_than_the_floor_are_clipped)
{
  const double temperature = 77.355;
  const float ceiling = 1.0e10f;
  std::vector<float> energy{0.0f, -20000.0f, -5000.0f, -100.0f};

  NLDFTExternalFieldStats stats =
      nldftMaskAndClipExternalPotential(energy, temperature, ceiling, nldftWellFloorInKT);

  EXPECT_EQ(stats.clippedVoxels, 2u);
  EXPECT_EQ(energy[0], 0.0f);
  EXPECT_NEAR(energy[1], -nldftWellFloorInKT * temperature, 1.0e-3);
  EXPECT_NEAR(energy[2], -nldftWellFloorInKT * temperature, 1.0e-3);
  EXPECT_EQ(energy[3], -100.0f);
}


// An empty box with several orientations, all U = 0: ρ_n = (1/M) Σ exp(...) reduces to the bulk gas.
TEST(energy_nldft, an_empty_box_with_orientations_holds_the_bulk_gas)
{
  NLDFTOptions options;
  options.pressurePoints = 2;
  options.xLow = 0.1;
  options.xHigh = 0.3;

  uint3 gridSize{32, 32, 32};
  UnitCell cell(20.0, 20.0, 20.0);
  const std::size_t N = 32 * 32 * 32;
  const std::size_t M = 8;
  std::vector<float> energy(N, 0.0f);
  std::vector<float> orientations(N * M, 0.0f);

  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options, orientations, M);
  ASSERT_EQ(result.isotherm.size(), 2u);
  EXPECT_EQ(result.unconverged, 0u);
  EXPECT_TRUE(result.molecular);
  EXPECT_EQ(result.numberOfOrientations, M);

  for (const IsothermPoint &point : result.isotherm)
  {
    double expected = nldftGasDensity(options, result.bulk, point.relativePressure) * cell.volume;
    EXPECT_NEAR(point.moleculesPerCell, expected, 0.02 * expected);
  }
}


// Molecular Henry is the orientation average, not the Helmholtz PMF stuffed into spherical ρ:
// ρ_n = ρ_b (1/M) Σ_ω exp(-βU_ω).
TEST(energy_nldft, molecular_henry_is_the_orientation_average)
{
  NLDFTOptions options;
  options.pressurePoints = 1;
  options.xLow = 1.0e-5;
  options.xHigh = 1.0e-5;

  uint3 gridSize{32, 32, 32};
  UnitCell cell(20.0, 20.0, 20.0);
  const std::size_t N = 32 * 32 * 32;
  const std::size_t M = 4;
  std::vector<float> energy(N, 0.0f);
  std::vector<float> orientations(N * M, 0.0f);
  for (std::size_t i = 0; i < N / 4; ++i)
  {
    orientations[i * M + 0] = -150.0f;  // one attractive orientation in a quarter of the box
  }

  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options, orientations, M);
  ASSERT_EQ(result.isotherm.size(), 1u);
  EXPECT_TRUE(result.molecular);

  double rhoBulk = nldftGasDensity(options, result.bulk, 1.0e-5);
  double boltzmann = 0.0;
  for (std::size_t i = 0; i < N; ++i)
  {
    double sum = 0.0;
    for (std::size_t o = 0; o < M; ++o)
      sum += std::exp(-static_cast<double>(orientations[i * M + o]) / options.temperature);
    boltzmann += sum / static_cast<double>(M);
  }
  double expected = rhoBulk * boltzmann * cell.volume / static_cast<double>(N);
  EXPECT_NEAR(result.isotherm.front().moleculesPerCell, expected, 0.05 * expected);
}


TEST(energy_nldft, a_fused_dumbbell_is_subcritical_at_nitrogen_temperature)
{
  NLDFTOptions options;
  options.dumbbellSiteOffsets = {-0.55, 0.55};
  options.dumbbellSiteDiameter = 3.306;
  NLDFTBulk bulk = nldftBulkCoexistence(options);

  EXPECT_TRUE(bulk.experimentalSaturation);
  EXPECT_NEAR(bulk.saturationPressure, nitrogenSaturationPressure, 1.0);
  EXPECT_GT(bulk.liquidDensity, 0.005);
  EXPECT_LT(bulk.liquidDensity, 0.030);
  EXPECT_LT(bulk.gasDensity, 1.0e-3);
}


// Two hard-sphere sites ±0.55 Å (TraPPE N2) in an empty box: the uniform gas is still the bulk of the
// same functional. Packing is the fused dumbbell (η = ρ V_union), not two independent spheres.
TEST(energy_nldft, an_empty_box_with_a_dumbbell_holds_the_bulk_gas)
{
  NLDFTOptions options;
  options.pressurePoints = 2;
  options.xLow = 0.1;
  options.xHigh = 0.3;
  options.dumbbellSiteOffsets = {-0.55, 0.55};
  options.dumbbellSiteDiameter = 3.306;

  uint3 gridSize{32, 32, 32};
  UnitCell cell(20.0, 20.0, 20.0);
  const std::size_t N = 32 * 32 * 32;
  const std::size_t M = 8;
  std::vector<float> energy(N, 0.0f);
  std::vector<float> orientations(N * M, 0.0f);

  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options, orientations, M);
  ASSERT_EQ(result.isotherm.size(), 2u);
  EXPECT_EQ(result.unconverged, 0u);
  EXPECT_TRUE(result.bulk.experimentalSaturation);
  EXPECT_NEAR(result.bulk.saturationPressure, nitrogenSaturationPressure, 1.0);
  EXPECT_GT(result.bulk.liquidDensity, 0.005);
  EXPECT_LT(result.bulk.liquidDensity, 0.030);

  for (const IsothermPoint &point : result.isotherm)
  {
    double expectedGas = nldftGasDensity(options, result.bulk, point.relativePressure) * cell.volume;
    EXPECT_NEAR(point.moleculesPerCell, expectedGas, 0.03 * expectedGas);
  }
}


// TraPPE N–N 12-6 mean field at ±0.55 Å, COM FMT. Empty box: the uniform gas of that functional.
TEST(energy_nldft, an_empty_box_with_trappe_site_attraction_holds_the_bulk_gas)
{
  NLDFTOptions options;
  options.pressurePoints = 2;
  options.xLow = 0.1;
  options.xHigh = 0.3;
  options.attractiveSiteOffsets = {-0.55, 0.55};
  options.attractiveSigma = 3.306;
  options.attractiveEpsilon = 36.0;

  uint3 gridSize{32, 32, 32};
  UnitCell cell(20.0, 20.0, 20.0);
  const std::size_t N = 32 * 32 * 32;
  const std::size_t M = 8;
  std::vector<float> energy(N, 0.0f);
  std::vector<float> orientations(N * M, 0.0f);

  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options, orientations, M);
  ASSERT_EQ(result.isotherm.size(), 2u);
  EXPECT_EQ(result.unconverged, 0u);
  EXPECT_LT(result.bulk.meanFieldIntegral, 0.0);
  EXPECT_TRUE(result.molecular);

  for (const IsothermPoint &point : result.isotherm)
  {
    double expectedGas = nldftGasDensity(options, result.bulk, point.relativePressure) * cell.volume;
    EXPECT_NEAR(point.moleculesPerCell, expectedGas, 0.03 * expectedGas);
  }
}


// CLI molecular path: fused-dumbbell FMT at TraPPE σ plus site-site LJ attraction together.
TEST(energy_nldft, an_empty_box_with_trappe_dumbbell_and_site_attraction_holds_the_bulk_gas)
{
  NLDFTOptions options;
  options.pressurePoints = 2;
  options.xLow = 0.1;
  options.xHigh = 0.3;
  options.attractiveSiteOffsets = {-0.55, 0.55};
  options.attractiveSigma = 3.306;
  options.attractiveEpsilon = 36.0;
  options.dumbbellSiteOffsets = {-0.55, 0.55};
  options.dumbbellSiteDiameter = 3.306;
  options.hardSphereDiameter = 3.306;
  options.sigma = 3.306;
  options.cutoff = 5.0 * 3.306;

  uint3 gridSize{32, 32, 32};
  UnitCell cell(20.0, 20.0, 20.0);
  const std::size_t N = 32 * 32 * 32;
  const std::size_t M = 8;
  std::vector<float> energy(N, 0.0f);
  std::vector<float> orientations(N * M, 0.0f);

  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options, orientations, M);
  ASSERT_EQ(result.isotherm.size(), 2u);
  EXPECT_EQ(result.unconverged, 0u);
  EXPECT_TRUE(result.bulk.experimentalSaturation);
  EXPECT_LT(result.bulk.meanFieldIntegral, 0.0);

  for (const IsothermPoint &point : result.isotherm)
  {
    double expectedGas = nldftGasDensity(options, result.bulk, point.relativePressure) * cell.volume;
    EXPECT_NEAR(point.moleculesPerCell, expectedGas, 0.03 * expectedGas);
  }
}


TEST(energy_nldft, trappe_site_attraction_mean_field_scales_as_site_count_squared)
{
  NLDFTOptions one;
  one.attractiveSiteOffsets = {0.0};
  one.attractiveSigma = 3.306;
  one.attractiveEpsilon = 36.0;
  NLDFTOptions two = one;
  two.attractiveSiteOffsets = {-0.55, 0.55};

  NLDFTBulk a1 = nldftBulkCoexistence(one);
  NLDFTBulk a2 = nldftBulkCoexistence(two);
  EXPECT_LT(a1.meanFieldIntegral, 0.0);
  EXPECT_NEAR(a2.meanFieldIntegral, 4.0 * a1.meanFieldIntegral, 1.0e-6);
}


// Strong grid Henry (FER-like): the window must open well below the old fixed x=1e-5 floor so the
// Type I rise is inside the solve, and still end at x=1 (experimental / isotherm P0).
TEST(energy_nldft, henry_pressure_window_matches_wham_policy)
{
  NLDFTOptions options;
  const double cellVolume = 1954.0;
  const double packing = cellVolume / nitrogenLiquidVolume;
  const double henryCoefficient = 236.0;  // /cell/Pa, TraPPE Widom scale
  const double henryAtP0 = henryCoefficient * nitrogenSaturationPressure;

  nldftApplyHenryPressureWindow(options, henryAtP0, cellVolume, nitrogenSaturationPressure);

  const double expectedPLow =
      std::min(nldftPressureWindowFloorPa, nldftPressureWindowFillFraction * packing / henryCoefficient);
  EXPECT_NEAR(options.xLow, expectedPLow / nitrogenSaturationPressure, 1.0e-12);
  EXPECT_NEAR(options.xHigh, 1.0, 1.0e-12);
  EXPECT_LT(options.xLow, 1.0e-5);

  // Weak binder: floor at 1 Pa, so xLow = 1/P0.
  NLDFTOptions weak;
  nldftApplyHenryPressureWindow(weak, 0.01 * nitrogenSaturationPressure, cellVolume, nitrogenSaturationPressure);
  EXPECT_NEAR(weak.xLow, nldftPressureWindowFloorPa / nitrogenSaturationPressure, 1.0e-12);
  EXPECT_NEAR(weak.xHigh, 1.0, 1.0e-12);
}
