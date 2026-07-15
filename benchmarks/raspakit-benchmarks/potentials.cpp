#include <benchmark/benchmark.h>

#include <cstddef>
#include <random>
#include <vector>

import forcefield;
import vdwparameters;
import pseudo_atom;
import potential_pair_derivatives;
import potential_pair_vdw;

// Micro-benchmark for the Lennard-Jones hot path of Potentials::potentialVDW<0>.
// Loops over a pre-generated set of squared distances and type pairs, mimicking the
// access pattern of the interaction loops.

static ForceField makeBenchmarkForceField()
{
  std::vector<PseudoAtom> pseudoAtoms{
      PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
      PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
      PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
      PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
      PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
  };
  std::vector<VDWParameters> parameters{VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72),
                                        VDWParameters(29.933, 2.745), VDWParameters(85.671, 3.017)};

  return ForceField(pseudoAtoms, parameters, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false,
                    true);
}

static void BM_PotentialVDWEnergyLennardJones(benchmark::State& state)
{
  ForceField forceField = makeBenchmarkForceField();

  constexpr std::size_t N = 4096;
  std::mt19937_64 gen(12345);
  std::uniform_real_distribution<double> dist(0.8, 12.0);
  std::uniform_int_distribution<std::size_t> typeDist(0, 4);

  std::vector<double> rrs(N);
  std::vector<std::size_t> typeAs(N);
  std::vector<std::size_t> typeBs(N);
  for (std::size_t i = 0; i < N; ++i)
  {
    double r = dist(gen);
    rrs[i] = r * r;
    typeAs[i] = typeDist(gen);
    typeBs[i] = typeDist(gen);
  }

  for (auto _ : state)
  {
    double energySum = 0.0;
    double dudlambdaSum = 0.0;
    for (std::size_t i = 0; i < N; ++i)
    {
      Potentials::PairDerivatives<0> ef =
          Potentials::potentialVDW<0>(forceField, 1.0, 1.0, rrs[i], typeAs[i], typeBs[i]);
      energySum += ef.energy;
      dudlambdaSum += ef.dUdlambda;
    }
    benchmark::DoNotOptimize(energySum);
    benchmark::DoNotOptimize(dudlambdaSum);
  }
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations() * N));
}

static void BM_PotentialVDWEnergyLennardJonesFractional(benchmark::State& state)
{
  ForceField forceField = makeBenchmarkForceField();

  constexpr std::size_t N = 4096;
  std::mt19937_64 gen(12345);
  std::uniform_real_distribution<double> dist(0.8, 12.0);
  std::uniform_int_distribution<std::size_t> typeDist(0, 4);

  std::vector<double> rrs(N);
  std::vector<std::size_t> typeAs(N);
  std::vector<std::size_t> typeBs(N);
  for (std::size_t i = 0; i < N; ++i)
  {
    double r = dist(gen);
    rrs[i] = r * r;
    typeAs[i] = typeDist(gen);
    typeBs[i] = typeDist(gen);
  }

  for (auto _ : state)
  {
    double energySum = 0.0;
    double dudlambdaSum = 0.0;
    for (std::size_t i = 0; i < N; ++i)
    {
      Potentials::PairDerivatives<0> ef =
          Potentials::potentialVDW<0>(forceField, 0.65, 1.0, rrs[i], typeAs[i], typeBs[i]);
      energySum += ef.energy;
      dudlambdaSum += ef.dUdlambda;
    }
    benchmark::DoNotOptimize(energySum);
    benchmark::DoNotOptimize(dudlambdaSum);
  }
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations() * N));
}

BENCHMARK(BM_PotentialVDWEnergyLennardJones);
BENCHMARK(BM_PotentialVDWEnergyLennardJonesFractional);
