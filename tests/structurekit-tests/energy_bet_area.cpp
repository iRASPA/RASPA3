#include <gtest/gtest.h>

import std;

import energy_shared_well_surface;
import energy_shared_bet_surface_area;

// The BET fit, on samples whose isotherm is known in closed form.

namespace
{
SheetPatch patchOf(double area, double energy)
{
  SheetPatch patch;
  patch.area = area;
  patch.energy = {energy, energy, energy};
  return patch;
}
}  // namespace


// Homogeneous sites, Henry matched to n_m C. The BET plot is exactly linear and the fit recovers n_m and C.
TEST(energy_bet_area, homogeneous_sites_recover_capacity_and_C)
{
  const double nMono = 10.0;
  const double C = 80.0;
  const double thermalEnergy = 1.0;
  const double heatOfLiquefaction = 0.0;
  const double energy = -thermalEnergy * std::log(C);

  std::vector<SheetPatch> patches;
  patches.push_back(patchOf(nMono * nitrogenCrossSection, energy));

  const double henry = nMono * C;
  BETSurfaceArea bet =
      BETSurfaceArea::fromSamples(patches, {}, henry, thermalEnergy, heatOfLiquefaction, 1.0, 1000.0);

  EXPECT_NEAR(bet.henryAnchor, 1.0, 1.0e-6);
  EXPECT_NEAR(bet.monolayerCapacity, nMono, 0.05 * nMono);
  EXPECT_NEAR(bet.cConstant, C, 0.08 * C);
  EXPECT_GT(bet.rSquared, 0.999);
  EXPECT_NEAR(bet.filamentCapacity, 0.0, 1.0e-12);
  EXPECT_NEAR(bet.sheetCapacity, nMono, 1.0e-12);
  EXPECT_NEAR(bet.filamentFileCapacity, 0.0, 1.0e-12);
  EXPECT_NEAR(bet.filamentMidplaneCapacity, 0.0, 1.0e-12);
}


// A deep filament sitting next to a moderate sheet inflates the fitted n_m: the micropore filling is
// already saturated in the BET window and acts as a constant offset.
TEST(energy_bet_area, filament_inflates_the_monolayer_capacity)
{
  const double nSheet = 10.0;
  const double C = 80.0;
  const double thermalEnergy = 1.0;
  const double energy = -thermalEnergy * std::log(C);

  std::vector<SheetPatch> patches;
  patches.push_back(patchOf(nSheet * nitrogenCrossSection, energy));

  FilamentVoxel voxel;
  voxel.volume = 5.0 * nitrogenLiquidVolume;
  voxel.energy = -thermalEnergy * std::log(1.0e4);

  const double henrySheet = nSheet * C;
  const double henryFilament = (voxel.volume / nitrogenLiquidVolume) * 1.0e4;
  BETSurfaceArea withFilament = BETSurfaceArea::fromSamples(patches, std::span(&voxel, 1), henrySheet + henryFilament,
                                                            thermalEnergy, 0.0, 1.0, 1000.0);
  BETSurfaceArea sheetOnly =
      BETSurfaceArea::fromSamples(patches, {}, henrySheet, thermalEnergy, 0.0, 1.0, 1000.0);

  EXPECT_GT(withFilament.monolayerCapacity, sheetOnly.monolayerCapacity);
  EXPECT_NEAR(withFilament.filamentCapacity, 5.0, 1.0e-9);
  EXPECT_NEAR(withFilament.filamentBlobCapacity, 5.0, 1.0e-9);
  EXPECT_NEAR(withFilament.filamentFileCapacity, 0.0, 1.0e-12);
  EXPECT_NEAR(withFilament.filamentMidplaneCapacity, 0.0, 1.0e-12);
}


// A merged-well filament is a line of molecules. Length / v_L^{1/3} is how many fit; using the volume of
// the tube of centres would undercount by the factor by which that tube is thinner than a liquid N2.
TEST(energy_bet_area, filament_length_is_a_one_dimensional_file)
{
  const double thermalEnergy = 1.0;
  const double C = 1.0e4;
  const double length = 8.72;
  const double expectedCapacity = length / nitrogenLinearSpacing;

  FilamentVoxel voxel;
  voxel.volume = 3.33;
  voxel.length = length;
  voxel.energy = -thermalEnergy * std::log(C);

  const double henry = expectedCapacity * C;
  BETSurfaceArea bet =
      BETSurfaceArea::fromSamples({}, std::span(&voxel, 1), henry, thermalEnergy, 0.0, 721.002, 642.88);

  EXPECT_NEAR(bet.filamentCapacity, expectedCapacity, 1.0e-9);
  EXPECT_NEAR(bet.filamentFileCapacity, expectedCapacity, 1.0e-9);
  EXPECT_NEAR(bet.filamentMidplaneCapacity, 0.0, 1.0e-12);
  EXPECT_NEAR(bet.monolayerCapacity, expectedCapacity, 0.15 * expectedCapacity);
  EXPECT_GT(bet.cConstant, 100.0);
}


// A 1-D file whose Henry comes from the whole cell, as in BIK: many voxels of mixed depth must not
// fragment the Langmuir site, or the coldest one steals the Henry match and n_m collapses. One site
// has C = Henry / n_fil, the isotherm saturates at the filling, and the fit recovers that capacity.
TEST(energy_bet_area, a_filament_file_recovers_the_filling_as_n_m)
{
  const double thermalEnergy = 1.0;
  const double nFil = 2.2763;
  const double henry = 78.73;
  const double mass = 721.002;
  const double volume = 642.88;

  FilamentVoxel deep;
  deep.volume = 0.1;
  deep.length = 0.1 * nFil * nitrogenLinearSpacing;
  deep.energy = -20.0 * thermalEnergy;
  FilamentVoxel shallow;
  shallow.volume = 0.9;
  shallow.length = 0.9 * nFil * nitrogenLinearSpacing;
  shallow.energy = -0.1 * thermalEnergy;
  FilamentVoxel voxels[2] = {deep, shallow};

  BETSurfaceArea bet =
      BETSurfaceArea::fromSamples({}, voxels, henry, thermalEnergy, 0.0, mass, volume);

  EXPECT_NEAR(bet.filamentCapacity, nFil, 1.0e-6);
  EXPECT_NEAR(bet.monolayerCapacity, nFil, 0.15 * nFil);
  const double expectedArea = nFil * nitrogenCrossSection * 6.0221419947e3 / mass;
  EXPECT_NEAR(bet.gravimetricArea, expectedArea, 0.15 * expectedArea);
}


// A micropore-dominated sample, after Walton and Snurr (JACS 129, 8552 (2007)): the pores fill well below
// the conventional 0.05, so n(1-x) stops rising there and the window has to follow the filling step down.
// The fit over that window still recovers the total capacity, which is their central result.
TEST(energy_bet_area, micropore_filling_pulls_the_window_below_the_standard_range)
{
  const double thermalEnergy = 1.0;
  const double sheetCapacity = 10.0;
  const double sheetC = 2000.0;
  const double filamentCapacity = 10.0;
  const double filamentC = 1.0e5;

  std::vector<SheetPatch> patches;
  patches.push_back(patchOf(sheetCapacity * nitrogenCrossSection, -thermalEnergy * std::log(sheetC)));

  FilamentVoxel voxel;
  voxel.volume = filamentCapacity * nitrogenLiquidVolume;
  voxel.energy = -thermalEnergy * std::log(filamentC);

  const double henry = sheetCapacity * sheetC + filamentCapacity * filamentC;
  BETSurfaceArea bet =
      BETSurfaceArea::fromSamples(patches, std::span(&voxel, 1), henry, thermalEnergy, 0.0, 1.0, 1000.0);

  EXPECT_LT(bet.windowHigh, 0.05);
  EXPECT_GT(bet.windowHigh, 0.0);
  EXPECT_GT(bet.cConstant, 0.0);
  EXPECT_NEAR(bet.monolayerCapacity, sheetCapacity + filamentCapacity, 0.15 * (sheetCapacity + filamentCapacity));
}


// An ultramicropore, after Bae, Yazaydin and Snurr (Langmuir 26, 5475 (2010)): the pore volume holds only
// one monolayer, so the multilayer stack terminates, the isotherm turns type I and saturates at the filling
// step, and the consistent window drops orders of magnitude below the textbook 0.05-0.30 (their MFI window
// is 5e-5 to 1e-2). The fit over that low window still recovers the capacity.
TEST(energy_bet_area, finite_room_terminates_the_stack_and_the_window_follows_the_filling_step)
{
  const double thermalEnergy = 1.0;
  const double nMono = 10.0;
  const double C = 1.0e4;
  const double energy = -thermalEnergy * std::log(C);

  std::vector<SheetPatch> patches;
  patches.push_back(patchOf(nMono * nitrogenCrossSection, energy));

  const double henry = nMono * C;
  BETSurfaceArea unlimited =
      BETSurfaceArea::fromSamples(patches, {}, henry, thermalEnergy, 0.0, 1.0, 1000.0);
  BETSurfaceArea oneLayer =
      BETSurfaceArea::fromSamples(patches, {}, henry, thermalEnergy, 0.0, 1.0, 1000.0, 1.0);

  // With unlimited room the multilayer keeps n(1-x) rising and the window stays conventional.
  EXPECT_GT(unlimited.windowHigh, 0.05);

  // With room for a single layer the sample saturates at the filling step and the window follows it down.
  EXPECT_GT(oneLayer.windowHigh, 0.0);
  EXPECT_LT(oneLayer.windowHigh, 0.05);
  EXPECT_LT(oneLayer.windowLow, 0.01);
  EXPECT_GT(oneLayer.cConstant, 100.0);
  EXPECT_NEAR(oneLayer.monolayerCapacity, nMono, 0.15 * nMono);
}


// A 2-D merged well (a slit midplane) packs as area / 16.2, not as a 1-D file along the cylinder
// reconstruction of that pancake.
TEST(energy_bet_area, a_planar_filament_packs_as_area)
{
  const double thermalEnergy = 1.0;
  const double C = 1.0e4;
  const double nPlane = 8.0;

  FilamentVoxel voxel;
  voxel.volume = 1.0;
  voxel.area = nPlane * nitrogenCrossSection;
  voxel.energy = -thermalEnergy * std::log(C);

  const double henry = nPlane * C;
  BETSurfaceArea bet =
      BETSurfaceArea::fromSamples({}, std::span(&voxel, 1), henry, thermalEnergy, 0.0, 2050.664, 1954.33);

  EXPECT_NEAR(bet.filamentCapacity, nPlane, 1.0e-9);
  EXPECT_NEAR(bet.filamentMidplaneCapacity, nPlane, 1.0e-9);
  EXPECT_NEAR(bet.filamentFileCapacity, 0.0, 1.0e-12);
  EXPECT_NEAR(bet.monolayerCapacity, nPlane, 0.15 * nPlane);
  EXPECT_GT(bet.cConstant, 100.0);
}


// A deep 2-D core next to a moderate sheet must not steal the Henry match: the filament uses the
// sheet's C, both fill, and n_m is the sum.
TEST(energy_bet_area, a_deep_planar_filament_does_not_steal_the_sheet)
{
  const double thermalEnergy = 1.0;
  const double nSheet = 10.0;
  const double nPlane = 5.0;
  const double sheetC = 80.0;
  const double energy = -thermalEnergy * std::log(sheetC);

  std::vector<SheetPatch> patches;
  patches.push_back(patchOf(nSheet * nitrogenCrossSection, energy));

  FilamentVoxel voxel;
  voxel.volume = 1.0;
  voxel.area = nPlane * nitrogenCrossSection;
  voxel.energy = -20.0 * thermalEnergy;

  const double henry = (nSheet + nPlane) * sheetC;
  BETSurfaceArea bet =
      BETSurfaceArea::fromSamples(patches, std::span(&voxel, 1), henry, thermalEnergy, 0.0, 1.0, 1000.0);

  EXPECT_NEAR(bet.filamentCapacity, nPlane, 1.0e-9);
  EXPECT_NEAR(bet.filamentMidplaneCapacity, nPlane, 1.0e-9);
  EXPECT_NEAR(bet.sheetCapacity, nSheet, 1.0e-9);
  EXPECT_NEAR(bet.monolayerCapacity, nSheet + nPlane, 0.15 * (nSheet + nPlane));
}


// The lattice-gas isotherm on geometric sites: Henry-anchored at x -> 0, and a well of depth U condenses
// at x = exp(U/kT)/f, filling to its capacity.
TEST(energy_bet_area, lattice_isotherm_is_henry_anchored_and_condenses_at_the_well_depth)
{
  const double kT = 1.0;
  const double qL = 8.0;  // W = -16, far below the mean-field critical point |W| = 4 kT
  const double capacity = 5.0;
  const double energy = -2.0;

  LatticeSite site{capacity, energy};
  // Pass the model's own bare Henry slope as the anchor, so the prefactor must come out as 1.
  const double henry = capacity * std::exp((-qL - energy) / kT);
  LatticeIsotherm model = latticeGasIsotherm(std::span(&site, 1), henry, kT, qL);

  EXPECT_NEAR(model.henryPrefactor, 1.0, 0.02);
  ASSERT_FALSE(model.isotherm.empty());
  EXPECT_NEAR(model.saturationCapacity, capacity, 1.0e-12);

  // Below the step the isotherm is the Henry line.
  const IsothermPoint &front = model.isotherm.front();
  EXPECT_NEAR(front.moleculesPerCell, henry * front.relativePressure, 0.05 * henry * front.relativePressure);

  // The step sits at x = exp(-2); by x = 0.35 the site is full.
  for (const IsothermPoint &point : model.isotherm)
  {
    if (point.relativePressure < 0.10) EXPECT_LT(point.moleculesPerCell, 0.05 * capacity);
    if (point.relativePressure > 0.20) EXPECT_GT(point.moleculesPerCell, 0.95 * capacity);
  }
}


// The lattice gas on a graph, where the attraction is between sites. A site that holds several molecules is
// that many molecules in one place, alike in everything and no neighbours of each other, and that is all it
// is: the isotherm of a graph with such sites in it is the isotherm of the graph with each of them written
// out as its molecules. That holds only if a site counts for its capacity everywhere a molecule is counted,
// in the loading, in the field its neighbours feel and in the grand potential that picks between the dilute
// branch and the filled one. The capacities here are unequal, so that leaving them out of the grand
// potential is not just a scale on it, and the wells are of different depths.
namespace
{
struct LatticeGraph
{
  std::vector<LatticeSite> sites;
  std::vector<std::vector<std::uint32_t>> neighbours;
};

// Groups of sites, every site a neighbour of every site outside its own group and of none inside it. With
// `merged` each group is one site holding as many molecules as the group has; without it the group is
// written out, a site to a molecule.
LatticeGraph groupsOf(std::span<const std::size_t> sizes, std::span<const double> energies, bool merged)
{
  LatticeGraph graph;
  std::vector<std::size_t> groupOfSite;
  for (std::size_t group = 0; group < sizes.size(); ++group)
  {
    if (merged)
    {
      graph.sites.push_back(LatticeSite{static_cast<double>(sizes[group]), energies[group]});
      groupOfSite.push_back(group);
      continue;
    }
    for (std::size_t member = 0; member < sizes[group]; ++member)
    {
      graph.sites.push_back(LatticeSite{1.0, energies[group]});
      groupOfSite.push_back(group);
    }
  }
  graph.neighbours.resize(graph.sites.size());
  for (std::size_t i = 0; i < graph.sites.size(); ++i)
  {
    for (std::size_t j = 0; j < graph.sites.size(); ++j)
    {
      if (groupOfSite[i] != groupOfSite[j]) graph.neighbours[i].push_back(static_cast<std::uint32_t>(j));
    }
  }
  return graph;
}
}  // namespace

TEST(energy_bet_area, a_site_that_holds_several_molecules_is_that_many_sites_of_one)
{
  const double kT = 1.0;
  const double qL = 8.0;
  const std::vector<std::size_t> sizes{3, 1, 2, 1, 3, 2};
  const std::vector<double> energies{-4.0, -7.0, -5.0, -6.5, -4.5, -5.5};

  double henry = 0.0;
  for (std::size_t group = 0; group < sizes.size(); ++group)
  {
    henry += static_cast<double>(sizes[group]) * std::exp((-qL - energies[group]) / kT);
  }

  const LatticeGraph merged = groupsOf(sizes, energies, true);
  const LatticeGraph writtenOut = groupsOf(sizes, energies, false);
  ASSERT_EQ(merged.sites.size(), 6uz);
  ASSERT_EQ(writtenOut.sites.size(), 12uz);

  const LatticeIsotherm one = latticeGasIsotherm(merged.sites, henry, kT, qL, merged.neighbours);
  const LatticeIsotherm other = latticeGasIsotherm(writtenOut.sites, henry, kT, qL, writtenOut.neighbours);

  EXPECT_NEAR(one.saturationCapacity, 12.0, 1.0e-12);
  EXPECT_NEAR(other.saturationCapacity, 12.0, 1.0e-12);
  ASSERT_EQ(one.isotherm.size(), other.isotherm.size());
  ASSERT_FALSE(one.isotherm.empty());

  // There is a step in it, the graph being as well connected as a liquid: this is asked across a
  // transition, where the grand potential has a choice to make.
  double largestStep = 0.0;
  for (std::size_t k = 1; k < other.isotherm.size(); ++k)
  {
    largestStep = std::max(largestStep, other.isotherm[k].moleculesPerCell - other.isotherm[k - 1].moleculesPerCell);
  }
  EXPECT_GT(largestStep, 3.0);

  for (std::size_t k = 0; k < one.isotherm.size(); ++k)
  {
    EXPECT_DOUBLE_EQ(one.isotherm[k].relativePressure, other.isotherm[k].relativePressure);
    EXPECT_NEAR(one.isotherm[k].moleculesPerCell, other.isotherm[k].moleculesPerCell, 1.0e-6)
        << "at x = " << one.isotherm[k].relativePressure;
  }
}


// A site with no room in it holds nothing, and is as good as not there: it adds nothing to the loading, and
// nothing to the field of its neighbours or to the grand potential either, however deep its well.
TEST(energy_bet_area, a_site_with_no_room_in_it_changes_nothing_on_the_graph)
{
  const double kT = 1.0;
  const double qL = 8.0;
  const std::vector<std::size_t> sizes{1, 1, 1, 1, 1, 1, 1, 1};
  const std::vector<double> energies{-4.0, -7.0, -5.0, -6.5, -4.5, -5.5, -6.0, -5.0};

  double henry = 0.0;
  for (double energy : energies) henry += std::exp((-qL - energy) / kT);

  const LatticeGraph without = groupsOf(sizes, energies, true);

  // The same, and a site of no capacity in a deep well that is a neighbour of them all.
  LatticeGraph with = without;
  const std::uint32_t empty = static_cast<std::uint32_t>(with.sites.size());
  with.sites.push_back(LatticeSite{0.0, -9.0});
  with.neighbours.emplace_back();
  for (std::uint32_t i = 0; i < empty; ++i)
  {
    with.neighbours[i].push_back(empty);
    with.neighbours[empty].push_back(i);
  }

  const LatticeIsotherm one = latticeGasIsotherm(without.sites, henry, kT, qL, without.neighbours);
  const LatticeIsotherm other = latticeGasIsotherm(with.sites, henry, kT, qL, with.neighbours);

  ASSERT_EQ(one.isotherm.size(), other.isotherm.size());
  ASSERT_FALSE(one.isotherm.empty());
  EXPECT_GT(one.isotherm.back().moleculesPerCell, 7.0);  // filled by the top of the range, so across a step
  for (std::size_t k = 0; k < one.isotherm.size(); ++k)
  {
    EXPECT_NEAR(one.isotherm[k].moleculesPerCell, other.isotherm[k].moleculesPerCell, 1.0e-6)
        << "at x = " << one.isotherm[k].relativePressure;
  }
}


// A heterogeneous micropore (a spread of well depths) through the full pipeline: the isotherm is type I,
// the Rouquerol cap lands at the filling step, and the fitted n_m reads near the geometric capacity.
TEST(energy_bet_area, lattice_micropore_reads_the_filling_capacity_as_the_monolayer)
{
  const double kT = 1.0;
  const double qL = 8.0;
  const double mass = 1000.0;
  const double cellVolume = 1000.0;

  std::vector<LatticeSite> sites;
  double henryRaw = 0.0;
  double totalCapacity = 0.0;
  for (int i = 0; i < 20; ++i)
  {
    // Depths spread uniformly over [-9, -3] kT: condensation from x ~ 1e-4 up to x ~ 0.05.
    LatticeSite site;
    site.capacity = 0.5;
    site.energy = -9.0 + 6.0 * static_cast<double>(i) / 19.0;
    sites.push_back(site);
    totalCapacity += site.capacity;
    henryRaw += site.capacity * std::exp((-qL - site.energy) / kT);
  }
  LatticeIsotherm model = latticeGasIsotherm(sites, henryRaw, kT, qL);

  double previous = 0.0;
  for (const IsothermPoint &point : model.isotherm)
  {
    EXPECT_GE(point.moleculesPerCell, previous * (1.0 - 1.0e-9));
    previous = point.moleculesPerCell;
  }
  EXPECT_NEAR(model.saturationCapacity, totalCapacity, 1.0e-12);

  BETSurfaceArea bet = BETSurfaceArea::fromIsotherm(std::move(model.isotherm), mass, cellVolume, nitrogenCrossSection,
                                                    nitrogenLiquidVolume);

  EXPECT_LT(bet.windowHigh, 0.10);
  EXPECT_NEAR(bet.monolayerCapacity, totalCapacity, 0.35 * totalCapacity);
}


// Fixed-window refit keeps the Rouquerol bounds from the full search and tracks a loading shift.
TEST(energy_bet_area, fixed_window_refit_keeps_window_and_tracks_loading)
{
  const double mass = 100.0;
  const double cellVolume = 1000.0;
  const double nMono = 20.0;
  const double C = 100.0;

  auto classicalBET = [&](double x) {
    const double denom = (1.0 - x) * (1.0 - x + C * x);
    return nMono * C * x / denom;
  };

  std::vector<IsothermPoint> isotherm;
  for (int i = 0; i < 120; ++i)
  {
    const double logX = std::log(1.0e-6) + (std::log(0.2) - std::log(1.0e-6)) * static_cast<double>(i) / 119.0;
    const double x = std::exp(logX);
    isotherm.push_back(IsothermPoint{x, classicalBET(x)});
  }

  const BETSurfaceArea full =
      BETSurfaceArea::fromIsotherm(isotherm, mass, cellVolume, nitrogenCrossSection, nitrogenLiquidVolume);
  ASSERT_FALSE(full.plateauReading);
  ASSERT_GT(full.monolayerCapacity, 0.0);
  ASSERT_GT(full.windowHigh, full.windowLow);

  const BETSurfaceArea sameWindow = BETSurfaceArea::fromIsothermFixedWindow(
      isotherm, mass, cellVolume, full.windowLow, full.windowHigh, nitrogenCrossSection, nitrogenLiquidVolume);
  EXPECT_DOUBLE_EQ(sameWindow.windowLow, full.windowLow);
  EXPECT_DOUBLE_EQ(sameWindow.windowHigh, full.windowHigh);
  EXPECT_NEAR(sameWindow.monolayerCapacity, full.monolayerCapacity, 1.0e-6 * full.monolayerCapacity);
  EXPECT_NEAR(sameWindow.gravimetricArea, full.gravimetricArea, 1.0e-6 * full.gravimetricArea);

  std::vector<IsothermPoint> shifted = isotherm;
  for (IsothermPoint& point : shifted)
  {
    point.moleculesPerCell *= 1.02;
  }
  const BETSurfaceArea shiftedFit = BETSurfaceArea::fromIsothermFixedWindow(
      shifted, mass, cellVolume, full.windowLow, full.windowHigh, nitrogenCrossSection, nitrogenLiquidVolume);
  EXPECT_DOUBLE_EQ(shiftedFit.windowLow, full.windowLow);
  EXPECT_DOUBLE_EQ(shiftedFit.windowHigh, full.windowHigh);
  EXPECT_GT(shiftedFit.gravimetricArea, full.gravimetricArea);
  EXPECT_NEAR(shiftedFit.gravimetricArea / full.gravimetricArea, 1.02, 0.01);
}


namespace
{
// A micropore that fills and a second, weaker kind of site that goes on taking molecules after it: two Langmuir terms,
// not a BET line anywhere, so the area read off it depends on the window, and n(1-x) has its maximum inside the range
// a window is taken from, flat as a maximum is.
double fillingThenRising(double x) { return 14.0 * 1.0e4 * x / (1.0 + 1.0e4 * x) + 6.0 * 30.0 * x / (1.0 + 30.0 * x); }

// The same isotherm at `points` pressures evenly apart in ln x from 1e-5 up, begun `shift` of a step later.
std::vector<IsothermPoint> tabulated(std::size_t points, double shift)
{
  const double step = std::log(1.0 / 1.0e-5) / static_cast<double>(points - 1);
  std::vector<IsothermPoint> isotherm;
  for (std::size_t index = 0; index < points; ++index)
  {
    const double x = 1.0e-5 * std::exp(step * (static_cast<double>(index) + shift));
    if (x > 1.0) break;
    isotherm.push_back(IsothermPoint{x, fillingThenRising(x)});
  }
  return isotherm;
}
}  // namespace

// An isotherm is a curve and the pressures it is tabulated at are a choice, so the area read off it should not turn on
// the choice. Read off the curve it does not: the ceiling of the window is where n(1-x) has its maximum, wherever the
// pressures fall around it, and the area is the same to a part in ten thousand over grids begun a tenth of a step
// apart and grids of other sizes. Read off the points, the ceiling is one of the points, and moves with them.
TEST(energy_bet_area, a_window_read_off_the_curve_does_not_turn_on_where_the_pressures_fall)
{
  const double mass = 5768.0;
  const double cellVolume = 5332.0;

  double ceiling = 0.0;
  double highest = 0.0;
  for (double x = 0.01; x < 0.3; x += 1.0e-5)
  {
    if (fillingThenRising(x) * (1.0 - x) > highest)
    {
      highest = fillingThenRising(x) * (1.0 - x);
      ceiling = x;
    }
  }
  ASSERT_GT(ceiling, 0.02);
  ASSERT_LT(ceiling, 0.2);

  std::vector<std::vector<IsothermPoint>> tabulations;
  for (int tenth = 0; tenth < 10; ++tenth) tabulations.push_back(tabulated(80, 0.1 * static_cast<double>(tenth)));
  for (std::size_t points : {60uz, 100uz, 120uz}) tabulations.push_back(tabulated(points, 0.0));

  const BETSurfaceArea first = BETSurfaceArea::fromIsotherm(tabulations.front(), mass, cellVolume, nitrogenCrossSection,
                                                            nitrogenLiquidVolume, WindowReading::Curve);
  ASSERT_GT(first.gravimetricArea, 0.0);
  double lowestCeilingOffPoints = 1.0;
  double highestCeilingOffPoints = 0.0;
  for (const std::vector<IsothermPoint> &isotherm : tabulations)
  {
    const BETSurfaceArea offCurve = BETSurfaceArea::fromIsotherm(isotherm, mass, cellVolume, nitrogenCrossSection,
                                                               nitrogenLiquidVolume, WindowReading::Curve);
    EXPECT_FALSE(offCurve.plateauReading);
    EXPECT_NEAR(offCurve.windowHigh, ceiling, 0.01 * ceiling) << isotherm.size() << " points";
    EXPECT_NEAR(offCurve.gravimetricArea, first.gravimetricArea, 2.0e-4 * first.gravimetricArea)
        << isotherm.size() << " points from x = " << isotherm.front().relativePressure;

    const BETSurfaceArea offPoints =
        BETSurfaceArea::fromIsotherm(isotherm, mass, cellVolume, nitrogenCrossSection, nitrogenLiquidVolume);
    lowestCeilingOffPoints = std::min(lowestCeilingOffPoints, offPoints.windowHigh);
    highestCeilingOffPoints = std::max(highestCeilingOffPoints, offPoints.windowHigh);
  }
  EXPECT_GT(highestCeilingOffPoints, 1.05 * lowestCeilingOffPoints)
      << "the points have to move the other reading's ceiling for this to have asked anything";
}

// On an isotherm that is a BET line the curve reading finds the line: the monolayer and C it was drawn with, to what
// a straight line in ln x between neighbouring points leaves of them. And an isotherm that is not in order of
// pressure has no curve through it, and is read off its points.
TEST(energy_bet_area, the_curve_reading_recovers_a_classical_isotherm)
{
  const double nMono = 20.0;
  const double C = 100.0;
  std::vector<IsothermPoint> isotherm;
  for (int i = 0; i < 80; ++i)
  {
    const double x = std::exp(std::log(1.0e-5) + (std::log(0.9) - std::log(1.0e-5)) * static_cast<double>(i) / 79.0);
    isotherm.push_back(IsothermPoint{x, nMono * C * x / ((1.0 - x) * (1.0 - x + C * x))});
  }
  const BETSurfaceArea fit = BETSurfaceArea::fromIsotherm(isotherm, 100.0, 1000.0, nitrogenCrossSection,
                                                          nitrogenLiquidVolume, WindowReading::Curve);
  EXPECT_FALSE(fit.plateauReading);
  EXPECT_NEAR(fit.monolayerCapacity, nMono, 5.0e-4 * nMono);
  EXPECT_NEAR(fit.cConstant, C, 2.0e-3 * C);
  EXPECT_LE(fit.windowHigh, 0.30);

  std::vector<IsothermPoint> backwards(isotherm.rbegin(), isotherm.rend());
  const BETSurfaceArea offPoints =
      BETSurfaceArea::fromIsotherm(backwards, 100.0, 1000.0, nitrogenCrossSection, nitrogenLiquidVolume);
  const BETSurfaceArea asked = BETSurfaceArea::fromIsotherm(backwards, 100.0, 1000.0, nitrogenCrossSection,
                                                            nitrogenLiquidVolume, WindowReading::Curve);
  EXPECT_DOUBLE_EQ(asked.monolayerCapacity, offPoints.monolayerCapacity);
  EXPECT_DOUBLE_EQ(asked.windowHigh, offPoints.windowHigh);
}


// Finite-layer (n-layer) BET recovers n_m, C and n from a synthetic isotherm.
TEST(energy_bet_area, finite_layer_bet_recovers_synthetic_isotherm)
{
  const double mass = 100.0;
  const double cellVolume = 1000.0;
  const double nMono = 12.0;
  const double C = 80.0;
  const double nLayers = 4.0;

  std::vector<IsothermPoint> isotherm;
  for (int i = 0; i < 160; ++i)
  {
    const double logX = std::log(1.0e-5) + (std::log(0.85) - std::log(1.0e-5)) * static_cast<double>(i) / 159.0;
    const double x = std::exp(logX);
    isotherm.push_back(IsothermPoint{x, nMono * finiteLayerBETOccupancy(x, C, nLayers)});
  }

  const FiniteLayerBETFit fit =
      fitFiniteLayerBET(isotherm, mass, cellVolume, nitrogenCrossSection, nitrogenLiquidVolume);
  ASSERT_TRUE(fit.ok);
  EXPECT_NEAR(fit.numberOfLayers, nLayers, 0.5);
  EXPECT_NEAR(fit.monolayerCapacity, nMono, 0.05 * nMono);
  EXPECT_NEAR(fit.cConstant, C, 0.15 * C);
  EXPECT_GT(fit.rSquared, 0.999);
  EXPECT_NEAR(fit.gravimetricArea, nMono * nitrogenCrossSection * 6.0221419947e3 / mass, 1.0);

  const FiniteLayerBETFit frozen = fitFiniteLayerBETFixedLayers(
      isotherm, mass, cellVolume, fit.numberOfLayers, fit.windowLow, fit.windowHigh, nitrogenCrossSection,
      nitrogenLiquidVolume);
  ASSERT_TRUE(frozen.ok);
  EXPECT_DOUBLE_EQ(frozen.numberOfLayers, fit.numberOfLayers);
  EXPECT_NEAR(frozen.monolayerCapacity, fit.monolayerCapacity, 1.0e-6 * fit.monolayerCapacity);
  EXPECT_NEAR(frozen.gravimetricArea, fit.gravimetricArea, 1.0e-6 * fit.gravimetricArea);
}


// Geometric packing for a 1-D filament: L / λ sites, which is what an ultramicroporous channel needs
// (volume packing of the filament tube undercounts by orders of magnitude).
TEST(energy_bet_area, geometric_sites_pack_a_filament_by_length)
{
  FilamentVoxel voxel;
  voxel.volume = 3.0;  // tiny tube of centres
  voxel.length = 2.0 * nitrogenLinearSpacing;
  voxel.energy = -5.0;

  GeometricSites geometry = geometricAdsorptionSites({}, std::span(&voxel, 1));
  EXPECT_NEAR(geometry.filamentFileCapacity, 2.0, 1.0e-12);
  EXPECT_NEAR(geometry.filamentCapacity, 2.0, 1.0e-12);
  ASSERT_EQ(geometry.sites.size(), 1u);
  EXPECT_NEAR(geometry.sites.front().capacity, 2.0, 1.0e-12);
}
