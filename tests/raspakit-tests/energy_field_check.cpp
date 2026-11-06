#include <stdlib.h>

#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import energy_backend;
import energy_shared_energy_backend;
import energy_shared_probe_energy_grid;
import energy_shared_energy_barrier;

// The check of a field against the same sum taken again on the processor, which the barrier route writes
// when RASPA_VERIFY_ENERGY_FIELD is set. A check is worth what its reference is worth: the reference has to
// be the sum the builders are meant to make, shift and all, or it reports a field that is right as wrong.

namespace
{

constexpr double edge = 14.0;
constexpr double sigma = 3.0;
constexpr double epsilon = 100.0;
constexpr double cutOff = 6.0;

// A truncated potential, shifted or not so that it reaches the cutoff at nought.
PairInteractions probeAndAtom(bool shifted)
{
  const double ratio6 = std::pow(sigma / cutOff, 6);
  const double shift = shifted ? 4.0 * epsilon * ratio6 * (ratio6 - 1.0) : 0.0;
  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "O"};
  interactions.charges = {0.0, 0.0};
  interactions.parameters.assign(4, PairParameters{sigma, epsilon, shift});
  interactions.cutOffVDW = cutOff;
  interactions.cutOffCoulomb = cutOff;
  return interactions;
}

// A few atoms scattered through a cell wide enough for the cutoff.
Crystal scatteredAtoms(const std::string &name)
{
  Crystal framework;
  framework.name = name;
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 64.0;
  for (const double3 &fractional :
       {double3(0.30, 0.20, 0.25), double3(0.70, 0.55, 0.40), double3(0.15, 0.80, 0.75), double3(0.60, 0.35, 0.90)})
  {
    framework.atoms.push_back(CrystalAtom{edge * fractional, 0.0, 1});
    framework.fractionalPositions.push_back(fractional);
  }
  return framework;
}

struct VerifyingTheField
{
  VerifyingTheField()
  {
#ifdef _WIN32
    _putenv_s("RASPA_VERIFY_ENERGY_FIELD", "1");
#else
    setenv("RASPA_VERIFY_ENERGY_FIELD", "1", 1);
#endif
  }
  ~VerifyingTheField()
  {
#ifdef _WIN32
    _putenv_s("RASPA_VERIFY_ENERGY_FIELD", "");
#else
    unsetenv("RASPA_VERIFY_ENERGY_FIELD");
#endif
  }
};

struct CheckedPoint
{
  double kernel{0.0};
  double reference{0.0};
  double difference{0.0};
};

// The barrier of the processor's field, and what the check of that field says.
std::vector<CheckedPoint> checkOf(const PairInteractions &interactions, const Crystal &framework)
{
  const ProbeEnergyGrid grid = cpuEnergyBackend().probeEnergyGrid(interactions, framework, "X", uint3{16, 16, 16});
  {
    VerifyingTheField verifying;
    EnergyBarrier barrier;
    barrier.run(interactions, framework, grid, 300.0);
  }

  std::vector<CheckedPoint> points;
  const std::string check = framework.name + ".energyfield.check.txt";
  std::ifstream written(check);
  for (std::string line; std::getline(written, line);)
  {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream values(line);
    std::size_t i = 0, j = 0, k = 0;
    CheckedPoint point;
    values >> i >> j >> k >> point.kernel >> point.reference >> point.difference;
    points.push_back(point);
  }
  written.close();
  std::filesystem::remove(check);
  std::filesystem::remove(framework.name + ".grid.barrier." + grid.backend + ".txt");
  return points;
}

}  // namespace


TEST(energy_field_check, the_reference_is_shifted_as_the_field_is)
{
  for (bool shifted : {false, true})
  {
    const PairInteractions interactions = probeAndAtom(shifted);
    const std::vector<CheckedPoint> points = checkOf(interactions, scatteredAtoms("energy_field_check_test"));
    ASSERT_GE(points.size(), 8uz);

    std::size_t withinReach = 0;
    for (const CheckedPoint &point : points)
    {
      if (std::fabs(point.reference) > 1.0e-3) ++withinReach;
      EXPECT_LT(point.difference, 1.0e-5)
          << (shifted ? "shifted" : "unshifted") << ": field " << point.kernel << ", reference " << point.reference;
    }
    EXPECT_GE(withinReach, 4uz) << "the points checked have to be within reach of an atom to say anything";
  }
}
