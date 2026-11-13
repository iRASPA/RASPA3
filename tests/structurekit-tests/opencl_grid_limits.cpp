#include <gtest/gtest.h>

#define CL_TARGET_OPENCL_VERSION 120
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import opencl;
import energy_opencl_surface_area;
import energy_opencl_void_fraction;
import energy_opencl_probe_energy_grid;

// The device routes of the energy grids without a device, and with a grid larger than their kernels can index.

namespace
{

constexpr double edge = 10.0;

PairInteractions probeAndAtom()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "O"};
  interactions.charges = {0.0, 0.0};
  interactions.parameters = {PairParameters{3.0, 100.0, 0.0}, PairParameters{3.0, 100.0, 0.0},
                             PairParameters{3.0, 100.0, 0.0}, PairParameters{3.0, 100.0, 0.0}};
  interactions.cutOffVDW = 4.0;
  interactions.cutOffCoulomb = 4.0;
  return interactions;
}

Crystal oneAtom()
{
  Crystal framework;
  framework.name = "opencl_grid_limits_test";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 16.0;
  framework.atoms.push_back(CrystalAtom{double3(5.0, 5.0, 5.0), 0.0, 1});
  framework.fractionalPositions.push_back(double3(0.5, 0.5, 0.5));
  return framework;
}

// Takes the device away for the scope, as on a machine without one, and gives it back afterwards.
class WithoutDevice
{
 public:
  WithoutDevice() : context_(OpenCL::clContext), device_(OpenCL::clDeviceId), queue_(OpenCL::clCommandQueue)
  {
    OpenCL::clContext.reset();
    OpenCL::clDeviceId.reset();
    OpenCL::clCommandQueue.reset();
  }
  ~WithoutDevice()
  {
    OpenCL::clContext = context_;
    OpenCL::clDeviceId = device_;
    OpenCL::clCommandQueue = queue_;
  }

 private:
  std::optional<cl_context> context_;
  std::optional<cl_device_id> device_;
  std::optional<cl_command_queue> queue_;
};

template <typename Call>
void expectRefusal(Call call, std::string_view message)
{
  try
  {
    call();
    ADD_FAILURE() << "no exception, expected one saying '" << message << "'";
  }
  catch (const std::runtime_error& error)
  {
    EXPECT_NE(std::string(error.what()).find(message), std::string::npos) << error.what();
  }
}

}  // namespace

// Without a device the surface-area route had a work-group size of nought to divide by, and the void-fraction
// route handed handles it never made to the device and let go of them again when it went.
TEST(opencl_grid_limits, the_energy_routes_refuse_to_run_without_a_device)
{
  OpenCL::initialize();
  const PairInteractions interactions = probeAndAtom();
  const Crystal framework = oneAtom();
  WithoutDevice withoutDevice;

  {
    EnergyOpenCLSurfaceArea route;
    expectRefusal([&] { route.run(interactions, framework, 0.0, "X", uint3{8, 8, 8}); }, "no OpenCL device");
  }
  {
    EnergyOpenCLVoidFraction route;
    expectRefusal([&] { route.run(interactions, framework, "X", uint3{8, 8, 8}); }, "no OpenCL device");
  }
}

// The kernels index the grid with a 32-bit int, and the count of points was a product of three 32-bit values
// that wrapped round at 2^32: a grid of 2048 x 2048 x 1024 points was given no room at all.
TEST(opencl_grid_limits, a_grid_larger_than_the_kernels_can_index_is_refused)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  expectRefusal([] { ProbeEnergyGridOpenCL::compute(probeAndAtom(), oneAtom(), "X", uint3{2048, 2048, 1024}); },
                "more than the device kernels can index");
}
