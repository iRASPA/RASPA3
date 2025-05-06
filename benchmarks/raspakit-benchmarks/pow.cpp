#include <benchmark/benchmark.h>

#include <cmath>
#include <string>

static void BM_LJPow(benchmark::State& state)
{
  double epsilon = 1.0;
  double sigma = 1.0;
  for (auto _ : state)
  {
    double dummy = 0.0;
    for (size_t i = 0; i < 10000; i++)
    {
      double r = i * 0.01 + 1.0;
      double scaledR = sigma / r;
      dummy += 4.0 * epsilon * (std::pow(scaledR, 12) - std::pow(scaledR, 6));
    }
    benchmark::DoNotOptimize(dummy);
  }
}

static void BM_LJSmart(benchmark::State& state)
{
  double epsilon = 1.0;
  double sigma = 1.0;
  for (auto _ : state)
  {
    double dummy = 0.0;
    for (size_t i = 0; i < 10000; i++)
    {
      double r = i * 0.01 + 1.0;
      double scaledR = sigma / r;
      double r3 = scaledR * scaledR * scaledR;
      double r6 = r3 * r3;
      dummy += 4.0 * epsilon * (r6 * (r6 - 1.0));
    }
    benchmark::DoNotOptimize(dummy);
  }
}

// Register the function as a benchmark
BENCHMARK(BM_LJPow);
BENCHMARK(BM_LJSmart);
