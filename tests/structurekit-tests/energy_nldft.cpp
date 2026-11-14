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


// The nitrogen dumbbell with site-site attraction, in an empty box at one and three atmospheres. The sites
// are deposited scaled so that their packing is the fused molecule's, and the attraction has to count them
// unscaled, or the pore fluid is short of a quarter of its attraction: at these densities that shifts the box
// away from the bulk gas by a few parts in a hundred, where the two otherwise agree to a part in a thousand.
// The liquid cannot be held this way: for the molecular models the relative pressure is measured against the
// experimental saturation pressure, not the model's own, so the model's liquid is not at the reservoir's
// chemical potential at x = 1.
TEST(energy_nldft, a_molecular_empty_box_holds_the_bulk_gas)
{
  NLDFTOptions options;
  options.pressurePoints = 2;
  options.xLow = 1.0;
  options.xHigh = 3.0;
  options.dumbbellSiteOffsets = {-0.55, 0.55};
  options.dumbbellSiteDiameter = 3.31;
  options.attractiveSiteOffsets = {-0.55, 0.55};
  options.attractiveSigma = 3.31;
  options.attractiveEpsilon = 36.0;

  uint3 gridSize{32, 32, 32};
  UnitCell cell(20.0, 20.0, 20.0);
  const std::size_t orientations = 8;
  std::vector<float> energy(32 * 32 * 32, 0.0f);
  std::vector<float> orientationEnergy(32 * 32 * 32 * orientations, 0.0f);

  NLDFTGridIsotherm result =
      nldftIsothermOnGrid(energy, gridSize, cell, options, orientationEnergy, orientations);
  ASSERT_EQ(result.isotherm.size(), 2u);
  EXPECT_EQ(result.unconverged, 0u);
  for (const IsothermPoint &point : result.isotherm)
  {
    double expected = nldftGasDensity(options, result.bulk, point.relativePressure) * cell.volume;
    EXPECT_NEAR(point.moleculesPerCell, expected, 0.01 * expected) << "x = " << point.relativePressure;
  }
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


namespace
{

// A narrow tube of deep well behind walls no molecule can enter, which fills in a step on the way up and has to be
// held on the way down. Nothing but the functional limits what it holds.
struct DeepTube
{
  static constexpr std::size_t n = 24;
  static constexpr std::size_t width = 5;
  uint3 gridSize{n, n, n};
  UnitCell cell{16.0, 16.0, 16.0};
  std::vector<float> energy = std::vector<float>(n * n * n, 1.0e5f);

  DeepTube()
  {
    for (std::size_t k = 0; k < n; ++k)
      for (std::size_t j = 0; j < width; ++j)
        for (std::size_t i = 0; i < width; ++i) energy[(k * n + j) * n + i] = -900.0f;
  }
};

NLDFTOptions tubeOptions()
{
  NLDFTOptions options;
  options.pressurePoints = 8;
  options.xLow = 1.0e-6;
  options.xHigh = 0.9;
  return options;
}

}  // namespace


// Orientations that all feel the same field are one spherical density, and a molecular density of them has to be that
// density: the same loading at every pressure of the sweep up and of the sweep down, through the filling of the tube.
// The two get there over different unknowns, the density at a voxel and what each orientation holds of it, and by the
// same functional. The iteration this replaced held a molecular density under a cap on the packing instead, and the
// cap was well under what the tube holds.
TEST(energy_nldft, orientations_that_are_all_alike_are_the_spherical_density)
{
  DeepTube tube;
  NLDFTOptions options = tubeOptions();
  const std::size_t M = 2;
  std::vector<float> orientations(tube.energy.size() * M);
  for (std::size_t i = 0; i < tube.energy.size(); ++i)
  {
    for (std::size_t o = 0; o < M; ++o) orientations[i * M + o] = tube.energy[i];
  }
  NLDFTGridIsotherm spherical = nldftIsothermOnGrid(tube.energy, tube.gridSize, tube.cell, options);
  NLDFTGridIsotherm molecular = nldftIsothermOnGrid(tube.energy, tube.gridSize, tube.cell, options, orientations, M);

  EXPECT_EQ(spherical.unconverged, 0uz);
  EXPECT_EQ(molecular.unconverged, 0uz);
  EXPECT_EQ(molecular.pointsWithoutASolution, 0uz);
  ASSERT_EQ(molecular.adsorptionIsotherm.size(), spherical.adsorptionIsotherm.size());
  ASSERT_GT(spherical.adsorptionIsotherm.back().moleculesPerCell,
            10.0 * spherical.adsorptionIsotherm.front().moleculesPerCell)
      << "the tube has to fill for this to have asked anything";

  for (std::size_t p = 0; p < spherical.adsorptionIsotherm.size(); ++p)
  {
    const NLDFTBranchPoint &up = spherical.adsorptionIsotherm[p];
    const NLDFTBranchPoint &down = spherical.desorptionIsotherm[p];
    EXPECT_NEAR(molecular.adsorptionIsotherm[p].moleculesPerCell, up.moleculesPerCell, 1.0e-3 * up.moleculesPerCell)
        << "on the way up at x = " << up.relativePressure;
    EXPECT_NEAR(molecular.desorptionIsotherm[p].moleculesPerCell, down.moleculesPerCell, 1.0e-3 * down.moleculesPerCell)
        << "on the way down at x = " << down.relativePressure;
  }
}

namespace
{

// The relative pressures of a molecular model are against the experimental saturation pressure, and the hard-sphere
// and attractive parts of each can be switched on apart, so each of these is its own functional.
void withDumbbell(NLDFTOptions &options)
{
  options.dumbbellSiteOffsets = {-0.55, 0.55};
  options.dumbbellSiteDiameter = 3.306;
  options.hardSphereDiameter = 3.306;
  options.sigma = 3.306;
  options.cutoff = 5.0 * 3.306;
}

void withSiteAttraction(NLDFTOptions &options)
{
  options.attractiveSiteOffsets = {-0.55, 0.55};
  options.attractiveSigma = 3.306;
  options.attractiveEpsilon = 36.0;
  options.hardSphereDiameter = 3.306;
  options.sigma = 3.306;
  options.cutoff = 5.0 * 3.306;
}

struct Model
{
  const char *name;
  bool molecular;
  bool dumbbell;
  bool siteAttraction;
};

const std::vector<Model> everyModel = {{"spherical", false, false, false},
                                       {"orientations alone", true, false, false},
                                       {"dumbbell", true, true, false},
                                       {"site attraction", true, false, true},
                                       {"dumbbell and site attraction", true, true, true}};

NLDFTGridIsotherm sweepOf(const Model &model, NLDFTOptions options, std::span<const float> energy, uint3 gridSize,
                          const UnitCell &cell, std::span<const float> orientations, std::size_t M)
{
  if (model.dumbbell) withDumbbell(options);
  if (model.siteAttraction) withSiteAttraction(options);
  return model.molecular ? nldftIsothermOnGrid(energy, gridSize, cell, options, orientations, M)
                         : nldftIsothermOnGrid(energy, gridSize, cell, options);
}

}  // namespace


// An empty box is uniform fluid, whose grand potential is minus its pressure times the volume, whichever model it is
// and whatever the iteration made of getting there. The pressure is the one the reservoir was set to.
TEST(energy_nldft, the_grand_potential_of_an_empty_box_is_minus_the_pressure_times_the_volume)
{
  const std::size_t N = 24 * 24 * 24;
  const std::size_t M = 8;
  UnitCell cell(20.0, 20.0, 20.0);
  std::vector<float> energy(N, 0.0f);
  std::vector<float> orientations(N * M, 0.0f);

  NLDFTOptions options;
  options.pressurePoints = 2;
  options.xLow = 0.1;
  options.xHigh = 0.3;
  for (const Model &model : everyModel)
  {
    NLDFTGridIsotherm result = sweepOf(model, options, energy, uint3{24, 24, 24}, cell, orientations, M);
    ASSERT_EQ(result.adsorptionIsotherm.size(), 2uz) << model.name;
    for (const NLDFTBranchPoint &point : result.adsorptionIsotherm)
    {
      ASSERT_TRUE(point.stationary) << model.name;
      const double pressure = point.relativePressure * result.bulk.saturationPressure / 1.380649e7;
      const double expected = -pressure * cell.volume / options.temperature;
      EXPECT_NEAR(point.grandPotential, expected, 5.0e-6 * std::abs(expected)) << model.name;
    }
  }
}

// At saturation the bulk's gas and liquid are at one pressure, which is what saturation means, with nothing else in
// common: about 175 times the molecules in this box. So they have one grand potential.
TEST(energy_nldft, gas_and_liquid_tie_at_saturation)
{
  const std::size_t N = 24 * 24 * 24;
  UnitCell cell(20.0, 20.0, 20.0);
  std::vector<float> energy(N, 0.0f);

  NLDFTOptions options;
  options.pressurePoints = 1;
  options.xLow = 1.0;
  options.xHigh = 1.0;
  NLDFTGridIsotherm gas = nldftIsothermOnGrid(energy, uint3{24, 24, 24}, cell, options);
  options.seedLiquid = true;
  NLDFTGridIsotherm liquid = nldftIsothermOnGrid(energy, uint3{24, 24, 24}, cell, options);

  ASSERT_TRUE(gas.adsorptionIsotherm.front().stationary);
  ASSERT_TRUE(liquid.adsorptionIsotherm.front().stationary);
  ASSERT_GT(liquid.adsorptionIsotherm.front().moleculesPerCell, 100.0 * gas.adsorptionIsotherm.front().moleculesPerCell);
  const double omega = gas.adsorptionIsotherm.front().grandPotential;
  EXPECT_NEAR(liquid.adsorptionIsotherm.front().grandPotential, omega, 1.0e-5 * std::abs(omega));
}

// d Omega / d mu = -N is what it takes for Omega and the update to be one functional, and here it is asked of profiles
// that are not uniform, in a well, with every term of the functional at work: the difference of Omega across a small
// step in pressure against the integral of the loading over it, which owes nothing to how Omega was added up.
TEST(energy_nldft, the_grand_potential_falls_with_the_chemical_potential_at_the_rate_of_the_loading)
{
  const std::size_t n = 24;
  const std::size_t N = n * n * n;
  const std::size_t M = 8;
  UnitCell cell(20.0, 20.0, 20.0);
  std::vector<float> energy(N);
  std::vector<float> orientations(N * M);
  for (std::size_t i = 0; i < N; ++i)
  {
    const float value = i / (n * n) < 8 ? -300.0f : 0.0f;
    energy[i] = value;
    for (std::size_t o = 0; o < M; ++o) orientations[i * M + o] = value * (1.0f + 0.05f * static_cast<float>(o));
  }

  NLDFTOptions options;
  options.pressurePoints = 5;
  options.xLow = 1.0e-3;
  options.xHigh = 1.2e-3;  // the trapezium is good to (d mu)^2 / 12 on a loading exponential in mu, which is 2e-4 here
  options.tolerance = 1.0e-8;
  for (const Model &model : everyModel)
  {
    NLDFTGridIsotherm result = sweepOf(model, options, energy, uint3{24, 24, 24}, cell, orientations, M);
    const std::vector<NLDFTBranchPoint> &branch = result.adsorptionIsotherm;
    ASSERT_EQ(branch.size(), 5uz) << model.name;
    for (std::size_t p = 0; p + 1 < branch.size(); ++p)
    {
      ASSERT_TRUE(branch[p].stationary && branch[p + 1].stationary) << model.name;
      const double step = branch[p + 1].chemicalPotential - branch[p].chemicalPotential;
      const double expected = -0.5 * (branch[p].moleculesPerCell + branch[p + 1].moleculesPerCell) * step;
      EXPECT_NEAR(branch[p + 1].grandPotential - branch[p].grandPotential, expected, 2.0e-3 * std::abs(expected))
          << model.name << " from x = " << branch[p].relativePressure;
    }
  }
}

// A slab of deep well in a cell it does not fill, which condenses with a loop: the sweep up stays empty past the point
// where filling would pay and the sweep down stays full below it. The equilibrium isotherm is neither. It is the
// branch of lower grand potential at each pressure, so inside the loop it leaves the sweep up at the crossing.
TEST(energy_nldft, the_isotherm_is_the_branch_of_lower_grand_potential_across_a_loop)
{
  const std::size_t n = 24;
  const std::size_t N = n * n * n;
  std::vector<float> energy(N, 0.0f);
  for (std::size_t i = 0; i < N / n * 6; ++i) energy[i] = -500.0f;

  NLDFTOptions options;
  options.pressurePoints = 10;
  options.xLow = 1.0e-3;
  options.xHigh = 0.9;
  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, uint3{n, n, n}, UnitCell(16.0, 16.0, 16.0), options);

  ASSERT_EQ(result.adsorptionIsotherm.size(), result.isotherm.size());
  ASSERT_EQ(result.desorptionIsotherm.size(), result.isotherm.size());

  std::size_t stayedOnTheWayUp = 0;  // inside a loop, the filled state metastable
  std::size_t tookTheWayDown = 0;    // inside a loop, the empty state metastable
  double firstFilled = 0.0;
  double previous = 0.0;
  for (std::size_t p = 0; p < result.isotherm.size(); ++p)
  {
    const NLDFTBranchPoint &up = result.adsorptionIsotherm[p];
    const NLDFTBranchPoint &down = result.desorptionIsotherm[p];
    EXPECT_DOUBLE_EQ(up.relativePressure, result.isotherm[p].relativePressure);
    EXPECT_DOUBLE_EQ(down.relativePressure, result.isotherm[p].relativePressure);

    const double loading = result.isotherm[p].moleculesPerCell;
    EXPECT_GE(loading, previous * 0.99) << "at x = " << up.relativePressure;
    previous = loading;

    // Only two profiles that solve the equation are two states to choose between.
    if (!up.stationary || !down.stationary) continue;
    if (down.moleculesPerCell < 1.05 * up.moleculesPerCell) continue;

    // A real loop: two states, and the lower grand potential is the one reported.
    if (down.grandPotential < up.grandPotential)
    {
      EXPECT_DOUBLE_EQ(loading, down.moleculesPerCell) << "at x = " << up.relativePressure;
      if (tookTheWayDown == 0) firstFilled = up.relativePressure;
      ++tookTheWayDown;
    }
    else
    {
      EXPECT_DOUBLE_EQ(loading, up.moleculesPerCell) << "at x = " << up.relativePressure;
      ++stayedOnTheWayUp;
    }
  }

  // Both halves of the loop have to be there for this to have asked anything.
  EXPECT_GT(stayedOnTheWayUp, 0uz);
  EXPECT_GT(tookTheWayDown, 0uz);
  EXPECT_DOUBLE_EQ(result.equilibriumTransition, firstFilled);
}

// A molecular profile is worked out from the density alone, and a sweep carries nothing from one pressure to the next
// but the profile it starts from. So where there is one state to find, a sweep that stops at a pressure on the way and
// one that does not end on the same loading and the same grand potential. The iteration this replaced carried the
// fields of its last step as well, and a point that had diverged could leave its own for the next.
TEST(energy_nldft, a_molecular_profile_owes_nothing_to_the_pressures_before_it)
{
  uint3 gridSize{16, 16, 16};
  const std::size_t n = 16;
  const std::size_t N = n * n * n;
  const std::size_t M = 4;
  UnitCell cell(14.0, 14.0, 14.0);

  std::vector<float> energy(N, 0.0f);
  std::vector<float> orientations(N * M, 0.0f);
  for (std::size_t i = 0; i < N; ++i)
  {
    const float value = i / (n * n) < 4 ? -400.0f : 0.0f;
    energy[i] = value;
    for (std::size_t o = 0; o < M; ++o) orientations[i * M + o] = value * (1.0f + 0.1f * static_cast<float>(o));
  }

  auto sweep = [&](int points)
  {
    NLDFTOptions options;
    options.pressurePoints = points;
    options.xLow = 1.0e-3 * std::pow(600.0, 3.0 / 5.0);
    options.xHigh = 0.6;
    options.attractiveSiteOffsets = {-0.55, 0.55};
    options.attractiveSigma = 3.306;
    options.attractiveEpsilon = 36.0;
    options.dumbbellSiteOffsets = {-0.55, 0.55};
    options.dumbbellSiteDiameter = 3.306;
    options.hardSphereDiameter = 3.306;
    options.sigma = 3.306;
    options.cutoff = 5.0;
    return nldftIsothermOnGrid(energy, gridSize, cell, options, orientations, M);
  };

  NLDFTGridIsotherm three = sweep(3);
  NLDFTGridIsotherm two = sweep(2);
  ASSERT_EQ(three.adsorptionIsotherm.size(), 3uz);
  ASSERT_EQ(two.adsorptionIsotherm.size(), 2uz);
  EXPECT_EQ(three.unconverged, 0uz);
  EXPECT_EQ(two.unconverged, 0uz);

  const NLDFTBranchPoint &stopped = three.adsorptionIsotherm.back();
  const NLDFTBranchPoint &direct = two.adsorptionIsotherm.back();
  ASSERT_TRUE(stopped.stationary);
  ASSERT_TRUE(direct.stationary);
  EXPECT_GT(direct.moleculesPerCell, 1.5 * two.adsorptionIsotherm.front().moleculesPerCell)
      << "the two ends of the sweep have to be different profiles";
  EXPECT_NEAR(stopped.moleculesPerCell, direct.moleculesPerCell, 1.0e-3 * direct.moleculesPerCell);
  EXPECT_NEAR(stopped.grandPotential, direct.grandPotential, 1.0e-4 * std::abs(direct.grandPotential));
}

// A pocket too small for two: every point of it overlaps every other, so the functional is at its zero-dimensional
// limit, which is where a density overlapping itself has to cost without bound. There the exact answer is known, with
// no functional in it: it holds nobody or one molecule, so N = z Q / (1 + z Q), for Q the integral of exp(-U/kT) over
// it. With the attraction off and the gas dilute, z is the density of the gas. The loading runs up to one molecule and
// stops there however hard the reservoir pushes, with n3 exponentially close to one and the equation still solved. The
// form of the third term that goes negative past |n2v|^2 = n2^2 / 3 has no such limit and runs on past one.
TEST(energy_nldft, a_pocket_that_fits_one_molecule_holds_one)
{
  const std::size_t n = 32;
  UnitCell cell(12.0, 12.0, 12.0);
  std::vector<float> energy(n * n * n, 1.0e6f);
  const double spacing = cell.lengthA / static_cast<double>(n);
  double boltzmannIntegral = 0.0;
  const float well = -900.0f;
  NLDFTOptions options;
  options.epsilon = 1.0e-6;
  options.pressurePoints = 7;
  options.xLow = 1.0e-4;
  options.xHigh = 1.0;
  const double voxel = cell.volume / static_cast<double>(n * n * n);
  for (std::size_t k = 0; k < n; ++k)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        const double dx = (static_cast<double>(i) - 16.0) * spacing;
        const double dy = (static_cast<double>(j) - 16.0) * spacing;
        const double dz = (static_cast<double>(k) - 16.0) * spacing;
        if (dx * dx + dy * dy + dz * dz < 1.0)
        {
          energy[(k * n + j) * n + i] = well;
          boltzmannIntegral += std::exp(-static_cast<double>(well) / options.temperature) * voxel;
        }
      }
    }
  }

  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, uint3{32, 32, 32}, cell, options);

  EXPECT_EQ(result.pointsWithoutASolution, 0uz);
  ASSERT_FALSE(result.adsorptionIsotherm.empty());
  for (const NLDFTBranchPoint &point : result.adsorptionIsotherm)
  {
    const double zQ = nldftGasDensity(options, result.bulk, point.relativePressure) * boltzmannIntegral;
    EXPECT_TRUE(point.stationary) << "at x = " << point.relativePressure;
    EXPECT_LT(point.moleculesPerCell, 1.0) << "at x = " << point.relativePressure;
    EXPECT_NEAR(point.moleculesPerCell, zQ / (1.0 + zQ), 0.02) << "at x = " << point.relativePressure;
  }
  EXPECT_GT(result.adsorptionIsotherm.back().moleculesPerCell, 0.9)
      << "the pocket has to come close to full for this to ask anything";
}

namespace
{
// A cubic lattice of round cages joined by windows along the axes, with the deep attractive skin of a zeolite wall and
// a floor at the depth a framework reaches. This is the shape that stopped the old iteration at a few molecules per
// cage.
struct CageLattice
{
  std::size_t n;
  UnitCell cell{16.0, 16.0, 16.0};
  std::vector<float> energy;

  explicit CageLattice(std::size_t points) : n(points), energy(points * points * points)
  {
    const double spacing = cell.lengthA / static_cast<double>(n);
    for (std::size_t k = 0; k < n; ++k)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        for (std::size_t i = 0; i < n; ++i)
        {
          auto wrapped = [&](std::size_t index)
          {
            const double value = static_cast<double>(index) * spacing;
            return value > 8.0 ? value - 16.0 : value;
          };
          const double x = wrapped(i);
          const double y = wrapped(j);
          const double z = wrapped(k);
          const double radial = std::sqrt(x * x + y * y + z * z);
          const double windowAxis =
              std::min({std::sqrt(y * y + z * z), std::sqrt(x * x + z * z), std::sqrt(x * x + y * y)});
          const double clearance = std::max(5.5 - radial, 1.8 - windowAxis);
          // the distance into the open space sets the energy, a well just inside the wall and a smooth rise out of it
          // into the wall, as a real field has: a step would leave the open volume to the grid
          const double value = clearance > 0.0 ? -850.0 * std::exp(-clearance * clearance / 1.5) - 100.0
                                               : -950.0 + 4000.0 * clearance * clearance;
          energy[(k * n + j) * n + i] = static_cast<float>(std::min(value, 1.0e6));
        }
      }
    }
  }

  uint3 gridSize() const
  {
    return uint3{static_cast<std::uint32_t>(n), static_cast<std::uint32_t>(n), static_cast<std::uint32_t>(n)};
  }
};
}  // namespace

// Every pressure from the first molecules up through the fill is a solution, and an isotherm of them behaves as one:
// the loading never falls and neither does anything but the grand potential. Most of this curve, on a field like this,
// used to be profiles that had stopped moving without solving anything.
TEST(energy_nldft, a_lattice_of_cages_solves_at_every_pressure_through_the_fill)
{
  CageLattice lattice(32);
  NLDFTOptions options;
  options.pressurePoints = 10;
  options.xLow = 1.0e-5;
  options.xHigh = 0.1;
  NLDFTGridIsotherm result = nldftIsothermOnGrid(lattice.energy, lattice.gridSize(), lattice.cell, options);

  EXPECT_EQ(result.pointsWithoutASolution, 0uz);
  EXPECT_EQ(result.unconverged, 0uz);
  const std::vector<NLDFTBranchPoint> &branch = result.adsorptionIsotherm;
  ASSERT_GE(branch.size(), 10uz);
  EXPECT_GT(branch.back().moleculesPerCell, 10.0 * branch.front().moleculesPerCell)
      << "a pore that never filled proves little";

  for (std::size_t p = 0; p < branch.size(); ++p)
  {
    EXPECT_TRUE(branch[p].stationary) << "at x = " << branch[p].relativePressure;
    if (p == 0) continue;
    EXPECT_GE(branch[p].moleculesPerCell, branch[p - 1].moleculesPerCell) << "at x = " << branch[p].relativePressure;
    EXPECT_LT(branch[p].grandPotential, branch[p - 1].grandPotential) << "at x = " << branch[p].relativePressure;
  }
}

// The descent ends on the residual, and it also ends where the grand potential and the loading have stopped moving,
// which is there for a pore packed close to solid: the last of these pressures, where the residual is still coming
// down hundreds of steps after neither of the two shows it. With no residual let through, every point here ends that
// way and no other. Where it ends has to be where the residual would have ended it, the loading to a few parts in ten
// thousand and the grand potential, which is at its minimum there, far closer. A stop that only looked steady, on a
// profile that could not move, would be neither.
TEST(energy_nldft, a_descent_that_ends_on_a_steady_grand_potential_ends_at_the_solution)
{
  CageLattice lattice(32);
  NLDFTOptions options;
  options.pressurePoints = 3;
  options.xLow = 1.0e-3;
  options.xHigh = 0.1;

  NLDFTOptions steady = options;
  steady.tolerance = 0.0;  // which no residual is under
  steady.steadyTolerance = 1.0e-5;
  NLDFTOptions solved = options;
  solved.tolerance = 1.0e-6;

  NLDFTGridIsotherm ended = nldftIsothermOnGrid(lattice.energy, lattice.gridSize(), lattice.cell, steady);
  NLDFTGridIsotherm reference = nldftIsothermOnGrid(lattice.energy, lattice.gridSize(), lattice.cell, solved);
  ASSERT_EQ(ended.adsorptionIsotherm.size(), 3uz);
  ASSERT_EQ(reference.adsorptionIsotherm.size(), 3uz);
  EXPECT_EQ(ended.unconverged, 0uz) << "a profile that has gone steady is one that converged";
  ASSERT_EQ(reference.unconverged, 0uz);
  EXPECT_GT(reference.adsorptionIsotherm.back().moleculesPerCell, 30.0)
      << "the cage has to be packed for the tail to be there";

  for (std::size_t p = 0; p < 3; ++p)
  {
    const NLDFTBranchPoint &point = ended.adsorptionIsotherm[p];
    const NLDFTBranchPoint &exact = reference.adsorptionIsotherm[p];
    ASSERT_TRUE(exact.stationary);
    ASSERT_LT(exact.residual, solved.tolerance);
    EXPECT_TRUE(point.stationary) << "at x = " << point.relativePressure;
    EXPECT_NEAR(point.moleculesPerCell, exact.moleculesPerCell, 5.0e-4 * exact.moleculesPerCell)
        << "at x = " << point.relativePressure;
    EXPECT_NEAR(point.grandPotential, exact.grandPotential, 1.0e-5 * std::abs(exact.grandPotential))
        << "at x = " << point.relativePressure;
  }
}

// The functional is a continuum one and the grid only a way to evaluate it, so the loading has to settle as the grid
// is refined, and it does to the size of the smoothing of the kernels. A kernel wrong at any wavevector shows here as a
// gap that does not close, however good each result looks on its own grid.
TEST(energy_nldft, the_loading_settles_as_the_grid_is_refined)
{
  auto top = [](std::size_t points)
  {
    CageLattice lattice(points);
    NLDFTOptions options;
    options.pressurePoints = 2;
    options.xLow = 3.0e-4;
    options.xHigh = 3.0e-3;
    options.targetSpacing = 0.0;
    NLDFTGridIsotherm result = nldftIsothermOnGrid(lattice.energy, lattice.gridSize(), lattice.cell, options);
    EXPECT_EQ(result.pointsWithoutASolution, 0uz);
    return result.isotherm.back().moleculesPerCell;
  };
  const double coarse = top(32);
  const double fine = top(48);
  EXPECT_NEAR(coarse, fine, 0.03 * fine);
}

// A cell is any cell: the weights are laid out through the inverse of its matrix, and a skew cell holds the same
// uniform gas as a square one. It would not if a sample landed anywhere but where the cell says it is.
TEST(energy_nldft, a_skew_empty_box_holds_the_bulk_gas)
{
  NLDFTOptions options;
  options.pressurePoints = 2;
  options.xLow = 0.1;
  options.xHigh = 0.3;

  const double degrees = std::numbers::pi / 180.0;
  UnitCell cell(20.0, 22.0, 24.0, 75.0 * degrees, 80.0 * degrees, 65.0 * degrees);
  std::vector<float> energy(32 * 32 * 32, 0.0f);
  NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, uint3{32, 32, 32}, cell, options);
  ASSERT_EQ(result.isotherm.size(), 2uz);
  for (const IsothermPoint &point : result.isotherm)
  {
    const double expected = nldftGasDensity(options, result.bulk, point.relativePressure) * cell.volume;
    EXPECT_NEAR(point.moleculesPerCell, expected, 1.0e-3 * expected);
  }
}


// A field finer than the density needs is brought down to the spacing it does need, a block of voxels to
// one, and each block is given the free energy -kT ln <exp(-U/kT)> of its voxels, so that the Boltzmann
// weight of the cell is what it was: a well one fine voxel wide is drawn wider and no deeper in total, and a
// block that is wall all through is wall. Here every block is of eight voxels. Half the box is wall, one
// voxel in eight of the rest is a well, and the gas is thin enough that the loading is the Henry limit.
TEST(energy_nldft, a_field_brought_down_to_a_coarser_grid_keeps_its_boltzmann_weight)
{
  NLDFTOptions options;
  options.pressurePoints = 1;
  options.xLow = 1.0e-5;
  options.xHigh = 1.0e-5;
  options.targetSpacing = 1.3;  // two voxels of 0.625 Å to a block along every axis

  const std::size_t n = 32;
  const uint3 gridSize{32, 32, 32};
  const UnitCell cell(20.0, 20.0, 20.0);
  const std::size_t N = n * n * n;
  const std::size_t M = 4;
  const float wall = 1.0e10f;
  const float well = -150.0f;

  std::vector<float> energy(N, 0.0f);
  std::vector<float> orientations(N * M, 0.0f);
  for (std::size_t k = 0; k < n; ++k)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        const std::size_t voxel = (k * n + j) * n + i;
        if (k < n / 2)
        {
          energy[voxel] = wall;
          for (std::size_t o = 0; o < M; ++o) orientations[voxel * M + o] = wall;
        }
        else if (i % 2 == 1 && j % 2 == 0 && k % 2 == 1)
        {
          // one corner of its block, and for the molecule one way round only
          energy[voxel] = well;
          orientations[voxel * M + 1] = well;
        }
      }
    }
  }

  auto boltzmannOf = [&](std::span<const float> field)
  {
    double sum = 0.0;
    for (float u : field) sum += std::exp(-static_cast<double>(u) / options.temperature);
    return sum / static_cast<double>(field.size());
  };

  NLDFTGridIsotherm spherical = nldftIsothermOnGrid(energy, gridSize, cell, options);
  ASSERT_EQ(spherical.isotherm.size(), 1u);
  EXPECT_EQ(spherical.pointsWithoutASolution, 0u);
  const double rhoBulk = nldftGasDensity(options, spherical.bulk, 1.0e-5);
  const double expectedSpherical = rhoBulk * cell.volume * boltzmannOf(energy);
  EXPECT_NEAR(spherical.isotherm.front().moleculesPerCell, expectedSpherical, 0.01 * expectedSpherical);

  NLDFTGridIsotherm molecular = nldftIsothermOnGrid(energy, gridSize, cell, options, orientations, M);
  ASSERT_EQ(molecular.isotherm.size(), 1u);
  EXPECT_TRUE(molecular.molecular);
  EXPECT_EQ(molecular.pointsWithoutASolution, 0u);
  const double expectedMolecular = rhoBulk * cell.volume * boltzmannOf(orientations);
  EXPECT_NEAR(molecular.isotherm.front().moleculesPerCell, expectedMolecular, 0.01 * expectedMolecular);

  // the wells are worth something here, or this has asked nothing
  EXPECT_GT(boltzmannOf(energy), 1.2 * 0.5);
}


// The free energy of a block is a number however far down one of its states lies. An orientation energy is
// not held to any floor, and exp(-U/kT) of one more than 709 kT down is past what a double holds: summed as
// it stood the block's weight was infinite and its free energy minus infinity, and no pressure could be
// solved for. Worked out from the lowest energy of the block down, the weights are all between nought and
// one. A state like that is a site that holds its one molecule at any pressure at all, and how far past
// 709 kT it lies changes nothing.
TEST(energy_nldft, a_state_past_what_an_exponential_holds_is_still_a_site_for_one_molecule)
{
  const std::size_t n = 32;
  const uint3 gridSize{32, 32, 32};
  const UnitCell cell(20.0, 20.0, 20.0);
  const std::size_t N = n * n * n;
  const std::size_t M = 4;
  const std::size_t voxel = (16 * n + 16) * n + 16;

  auto loadingsWith = [&](double depthInKT)
  {
    NLDFTOptions options;
    options.pressurePoints = 2;
    options.xLow = 1.0e-5;
    options.xHigh = 1.0e-3;
    options.targetSpacing = 1.3;  // two voxels to a block along every axis

    std::vector<float> energy(N, 0.0f);
    std::vector<float> orientations(N * M, 0.0f);
    orientations[voxel * M + 0] = static_cast<float>(-depthInKT * options.temperature);

    NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options, orientations, M);
    EXPECT_EQ(result.pointsWithoutASolution, 0u) << depthInKT << " kT down";
    EXPECT_EQ(result.unconverged, 0u) << depthInKT << " kT down";
    std::vector<double> loadings;
    for (const IsothermPoint &point : result.isotherm) loadings.push_back(point.moleculesPerCell);
    return loadings;
  };

  ASSERT_LT(std::exp(705.0), std::numeric_limits<double>::max());
  const std::vector<double> within = loadingsWith(705.0);
  const std::vector<double> past = loadingsWith(720.0);
  const std::vector<double> farPast = loadingsWith(5000.0);

  ASSERT_EQ(within.size(), 2u);
  ASSERT_EQ(past.size(), 2u);
  ASSERT_EQ(farPast.size(), 2u);
  for (std::size_t point = 0; point < 2; ++point)
  {
    // one molecule, to what the descent is stopped at
    EXPECT_NEAR(within[point], 1.0, 0.01);
    EXPECT_NEAR(past[point], within[point], 1.0e-3);
    EXPECT_NEAR(farPast[point], within[point], 1.0e-3);
  }
}

// The spherical solver (one orientation) holds its pore against the spherical bulk, whatever site offsets the
// options carry: the bulk used to be the multi-site fluid while the pore solved the one-site functional.
TEST(energy_nldft, the_spherical_solver_ignores_site_offsets)
{
  NLDFTOptions spherical;
  spherical.pressurePoints = 2;
  spherical.xLow = 0.1;
  spherical.xHigh = 0.3;
  NLDFTOptions withSites = spherical;
  withSites.dumbbellSiteOffsets = {-0.55, 0.55};
  withSites.dumbbellSiteDiameter = 3.306;

  uint3 gridSize{8, 8, 8};
  UnitCell cell(20.0, 20.0, 20.0);
  std::vector<float> energy(8 * 8 * 8, 0.0f);
  const NLDFTGridIsotherm reference = nldftIsothermOnGrid(energy, gridSize, cell, spherical);
  const NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, withSites);
  EXPECT_FALSE(result.molecular);
  EXPECT_DOUBLE_EQ(result.bulk.liquidDensity, reference.bulk.liquidDensity);
  EXPECT_DOUBLE_EQ(result.bulk.saturationPressure, reference.bulk.saturationPressure);
  ASSERT_EQ(result.isotherm.size(), reference.isotherm.size());
  for (std::size_t i = 0; i < result.isotherm.size(); ++i)
    EXPECT_NEAR(result.isotherm[i].moleculesPerCell, reference.isotherm[i].moleculesPerCell,
                1e-9 * std::max(1.0, reference.isotherm[i].moleculesPerCell));
}
