#include <benchmark/benchmark.h>

#include <fstream>
#include <ios>
#include <iostream>
#include <print>

import int3;
import randomnumbers;
import forcefield;
import framework;
import component;
import system;
import atom;
import factory;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_framework_molecule_grid;
import interactions_external_field;
import interactions_ewald;
import energy_status;
import interpolation_energy_grid;
import running_energy;

// InterpolationEnergyGrid::makeInterpolationGrid
// ----------------------------------------------------------------------------------------------------

template <ForceField::InterpolationScheme Scheme, ForceField::InterpolationGridType GridType>
static void BM_EnergyGridCreationMFI(benchmark::State& state)
{
  int N = static_cast<int>(state.range(0));

  ForceField forceField =
      TestFactories::makeDefaultFF(12.0, true, false, GridType == ForceField::InterpolationGridType::EwaldReal);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(1, 1, 1));

  int3 numberOfGridPoints{N, N, N};
  InterpolationEnergyGrid grid(f.simulationBox, numberOfGridPoints, Scheme);

  std::stringstream stream;
  for (auto _ : state)
  {
    grid.makeInterpolationGrid(stream, GridType, forceField, f, 12.0, 5);
  }

  constexpr const char* schemeName = Scheme == ForceField::InterpolationScheme::Tricubic ? "Tricubic" : "Triquintic";
  constexpr const char* interactionName =
      GridType == ForceField::InterpolationGridType::LennardJones ? "VDW" : "EwaldReal";
  state.SetLabel(std::format("{} {}x{}x{} {}", schemeName, N, N, N, interactionName));
}

template <ForceField::InterpolationScheme Scheme, ForceField::InterpolationGridType GridType>
static void BM_EnergyGridCreationCHA(benchmark::State& state)
{
  int N = static_cast<int>(state.range(0));

  ForceField forceField =
      TestFactories::makeDefaultFF(12.0, true, false, GridType == ForceField::InterpolationGridType::EwaldReal);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));

  int3 numberOfGridPoints{N, N, N};
  InterpolationEnergyGrid grid(f.simulationBox, numberOfGridPoints, Scheme);

  std::stringstream stream;
  for (auto _ : state)
  {
    grid.makeInterpolationGrid(stream, GridType, forceField, f, 12.0, 5);
  }

  constexpr const char* schemeName = Scheme == ForceField::InterpolationScheme::Tricubic ? "Tricubic" : "Triquintic";
  constexpr const char* interactionName =
      GridType == ForceField::InterpolationGridType::LennardJones ? "VDW" : "EwaldReal";
  state.SetLabel(std::format("{} {}x{}x{} {}", schemeName, N, N, N, interactionName));
}

// Interactions::computeFrameworkMoleculeEnergy
// ----------------------------------------------------------------------------------------------------
static const InterpolationEnergyGrid& getGrid(ForceField::InterpolationScheme scheme,
                                              ForceField::InterpolationGridType type)
{
  // Make sure grids are only computed once by making them static within this benchmarking routine
  static std::optional<InterpolationEnergyGrid> cache[4];
  static std::once_flag flags[4];
  
  size_t scheme_idx = (scheme == ForceField::InterpolationScheme::Tricubic) ? 0 : 1;
  size_t type_idx = (type == ForceField::InterpolationGridType::LennardJones) ? 0 : 2;
  size_t idx = scheme_idx + type_idx;
  
  std::call_once(flags[idx],
    [&]
    {
        ForceField ff = TestFactories::makeDefaultFF(12.0, true, false,
            type == ForceField::InterpolationGridType::EwaldReal);
            Framework f = TestFactories::makeCHA(ff, int3(1, 1, 1));
            
            int3 numberOfGridPoints{64, 64, 64};
            cache[idx].emplace(f.simulationBox, numberOfGridPoints, scheme);
            
            std::stringstream null_stream;
                   cache[idx]->makeInterpolationGrid(null_stream, type, ff, f, 12.0,
                                                     type == ForceField::InterpolationGridType::EwaldReal ? 5 : 2);
                 });

  return *cache[idx];
}

static void BM_ComputeFrameworkEnergyVDW(benchmark::State& state)
{
  size_t N = static_cast<size_t>(state.range(0));
  size_t mode = static_cast<size_t>(state.range(1));
  int uc = static_cast<int>(state.range(2));

  RandomNumber rng(std::nullopt);
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, false);

  Component c = TestFactories::makeMethane(forceField, 0);
  Framework f = TestFactories::makeCHA(forceField, int3(uc, uc, uc));

  std::vector<Atom> atoms(N);
  for (auto& atom : atoms)
  {
    atom.position = f.simulationBox.randomPosition(rng);
  }

  std::optional<InterpolationEnergyGrid> grid;
  if (mode == 1)
  {
    grid = getGrid(ForceField::InterpolationScheme::Tricubic, ForceField::InterpolationGridType::LennardJones);
  }
  else if (mode == 2)
  {
    grid = getGrid(ForceField::InterpolationScheme::Triquintic, ForceField::InterpolationGridType::LennardJones);
  }

  for (auto _ : state)
  {
    auto _ = Interactions::computeFrameworkMoleculeEnergy(forceField, f.simulationBox, {grid}, {f}, f.atoms, atoms);
  }
  const char* names[]{"full", "tricubic", "triquintic"};
  state.SetLabel(std::format("{} Energy VDW, {} particles, {}x{}x{} uc", names[mode], N, uc, uc, uc));
}

static void BM_ComputeFrameworkEnergyEwald(benchmark::State& state)
{
  size_t N = static_cast<size_t>(state.range(0));
  size_t mode = static_cast<size_t>(state.range(1));
  int uc = static_cast<int>(state.range(2));

  RandomNumber rng(std::nullopt);
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  // Turn off cut off vdw interactions. dist is calculated for ewald anyway, only overhead is "if (r < rc)"
  forceField.cutOffFrameworkVDW = 0.0;

  Component c = TestFactories::makeIon(forceField, 0, "Na", 5, 1.0);
  Framework f = TestFactories::makeCHA(forceField, int3(uc, uc, uc));

  std::vector<Atom> atoms(N);
  for (auto& atom : atoms)
  {
    atom.position = f.simulationBox.randomPosition(rng);
  }

  std::optional<InterpolationEnergyGrid> grid;
  if (mode == 1)
  {
    grid = getGrid(ForceField::InterpolationScheme::Tricubic, ForceField::InterpolationGridType::EwaldReal);
  }
  else if (mode == 2)
  {
    grid = getGrid(ForceField::InterpolationScheme::Triquintic, ForceField::InterpolationGridType::EwaldReal);
  }

  for (auto _ : state)
  {
    auto _ = Interactions::computeFrameworkMoleculeEnergy(forceField, f.simulationBox, {grid}, {f}, f.atoms, atoms);
  }

  const char* names[]{"full", "tricubic", "triquintic"};
  state.SetLabel(std::format("{} Energy Ewald, {} particles, {}x{}x{} uc", names[mode], N, uc, uc, uc));
}

// Interactions::computeFrameworkMoleculeGradient
// ----------------------------------------------------------------------------------------------------
static void BM_ComputeFrameworkGradientVDW(benchmark::State& state)
{
  size_t N = static_cast<size_t>(state.range(0));
  size_t mode = static_cast<size_t>(state.range(1));
  int uc = static_cast<int>(state.range(2));

  RandomNumber rng(std::nullopt);
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, false);

  Component c = TestFactories::makeMethane(forceField, 0);
  Framework f = TestFactories::makeCHA(forceField, int3(uc, uc, uc));

  std::vector<Atom> atoms(N);
  for (auto& atom : atoms)
  {
    atom.position = f.simulationBox.randomPosition(rng);
  }

  std::optional<InterpolationEnergyGrid> grid;
  if (mode == 1)
  {
    grid = getGrid(ForceField::InterpolationScheme::Tricubic, ForceField::InterpolationGridType::LennardJones);
  }
  else if (mode == 2)
  {
    grid = getGrid(ForceField::InterpolationScheme::Triquintic, ForceField::InterpolationGridType::LennardJones);
  }

  for (auto _ : state)
  {
    auto _ = Interactions::computeFrameworkMoleculeGradient(forceField, f.simulationBox, f.atoms, atoms, {grid});
  }

  const char* names[]{"full", "tricubic", "triquintic"};
  state.SetLabel(std::format("{} Grad VDW, {} particles, {}x{}x{} uc", names[mode], N, uc, uc, uc));
}

static void BM_ComputeFrameworkGradientEwald(benchmark::State& state)
{
  size_t N = static_cast<size_t>(state.range(0));
  size_t mode = static_cast<size_t>(state.range(1));
  int uc = static_cast<int>(state.range(2));

  RandomNumber rng(std::nullopt);
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  // Turn off cut off vdw interactions. dist is calculated for ewald anyway, only overhead is "if (r < rc)"
  forceField.cutOffFrameworkVDW = 0.0;

  Component c = TestFactories::makeIon(forceField, 0, "Na", 5, 1.0);
  Framework f = TestFactories::makeCHA(forceField, int3(uc, uc, uc));

  std::vector<Atom> atoms(N);
  for (auto& atom : atoms)
  {
    atom.position = f.simulationBox.randomPosition(rng);
  }

  std::optional<InterpolationEnergyGrid> grid;
  if (mode == 1)
  {
    grid = getGrid(ForceField::InterpolationScheme::Tricubic, ForceField::InterpolationGridType::EwaldReal);
  }
  else if (mode == 2)
  {
    grid = getGrid(ForceField::InterpolationScheme::Triquintic, ForceField::InterpolationGridType::EwaldReal);
  }

  for (auto _ : state)
  {
    auto _ = Interactions::computeFrameworkMoleculeGradient(forceField, f.simulationBox, f.atoms, atoms, {grid});
  }

  const char* names[]{"full", "tricubic", "triquintic"};
  state.SetLabel(std::format("{} Grad Ewald, {} particles, {}x{}x{} uc", names[mode], N, uc, uc, uc));
}

// Run benchmarks
// ----------------------------------------------------------------------------------------------------

// Grid creation
// BENCHMARK_TEMPLATE(BM_EnergyGridCreation, ForceField::InterpolationScheme::Tricubic,
//                    ForceField::InterpolationGridType::LennardJones)
//     ->RangeMultiplier(2)
//     ->Range(1, 64);

// BENCHMARK_TEMPLATE(BM_EnergyGridCreation, ForceField::InterpolationScheme::Tricubic,
//                    ForceField::InterpolationGridType::EwaldReal)
//     ->RangeMultiplier(2)
//     ->Range(1, 64);

// BENCHMARK_TEMPLATE(BM_EnergyGridCreation, ForceField::InterpolationScheme::Triquintic,
//                    ForceField::InterpolationGridType::LennardJones)
//     ->RangeMultiplier(2)
//     ->Range(1, 64);

// BENCHMARK_TEMPLATE(BM_EnergyGridCreation, ForceField::InterpolationScheme::Triquintic,
//                    ForceField::InterpolationGridType::EwaldReal)
//     ->RangeMultiplier(2)
//     ->Range(1, 64);

// Energy computation
BENCHMARK(BM_ComputeFrameworkEnergyVDW)
    ->ArgsProduct({benchmark::CreateRange(10, 1000000, 10), {0, 1, 2}, {1, 2, 3, 4, 5, 6}});

BENCHMARK(BM_ComputeFrameworkEnergyEwald)
    ->ArgsProduct({benchmark::CreateRange(10, 1000000, 10), {0, 1, 2}, {1, 2, 3, 4, 5, 6}});

// Gradient computation
BENCHMARK(BM_ComputeFrameworkGradientVDW)
    ->ArgsProduct({benchmark::CreateRange(10, 1000000, 10), {0, 1, 2}, {1, 2, 3, 4, 5, 6}});

BENCHMARK(BM_ComputeFrameworkGradientEwald)
    ->ArgsProduct({benchmark::CreateRange(10, 1000000, 10), {0, 1, 2}, {1, 2, 3, 4, 5, 6}});
