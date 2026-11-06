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

import opencl;

// What a routine makes on the device is let go of when the routine is left, thrown out of or not. A handle
// cannot be asked whether it has been released once it has, so each is given a second holder here before it
// is handed over: when the routine has gone, the count of holders has to be back to that one.

namespace
{

const char *nothingMuch = "__kernel void Nothing(__global float *out) { out[get_global_id(0)] = 1.0f; }";

struct Made
{
  cl_program program{nullptr};
  cl_kernel kernel{nullptr};
  cl_mem first{nullptr};
  cl_mem second{nullptr};
};

// Builds a program, a kernel and two buffers, keeps a hold of its own on each, and hands them to `handles`.
Made makeAndHandOver(OpenCL::Handles &handles)
{
  Made made;
  cl_int err = CL_SUCCESS;

  made.program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &nothingMuch, nullptr, &err);
  EXPECT_EQ(err, CL_SUCCESS);
  EXPECT_EQ(clBuildProgram(made.program, 0, nullptr, nullptr, nullptr, nullptr), CL_SUCCESS);
  clRetainProgram(made.program);
  handles.program = made.program;

  made.kernel = clCreateKernel(made.program, "Nothing", &err);
  EXPECT_EQ(err, CL_SUCCESS);
  clRetainKernel(made.kernel);
  handles.kernel = made.kernel;

  made.first = handles.keep(clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, 64, nullptr, &err));
  EXPECT_EQ(err, CL_SUCCESS);
  clRetainMemObject(made.first);
  made.second = handles.keep(clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, 64, nullptr, &err));
  EXPECT_EQ(err, CL_SUCCESS);
  clRetainMemObject(made.second);
  return made;
}

cl_uint holdersOf(cl_mem buffer)
{
  cl_uint count = 0;
  clGetMemObjectInfo(buffer, CL_MEM_REFERENCE_COUNT, sizeof(count), &count, nullptr);
  return count;
}

cl_uint holdersOf(cl_kernel kernel)
{
  cl_uint count = 0;
  clGetKernelInfo(kernel, CL_KERNEL_REFERENCE_COUNT, sizeof(count), &count, nullptr);
  return count;
}

// Only this test's own hold is left on each, and it is given up.
void expectOnlyOursLeftAndLetGo(const Made &made)
{
  EXPECT_EQ(holdersOf(made.first), 1u);
  EXPECT_EQ(holdersOf(made.second), 1u);
  EXPECT_EQ(holdersOf(made.kernel), 1u);

  // The program is held by its kernel as well as by us until the kernel goes.
  clReleaseKernel(made.kernel);
  cl_uint programHolders = 0;
  clGetProgramInfo(made.program, CL_PROGRAM_REFERENCE_COUNT, sizeof(programHolders), &programHolders, nullptr);
  EXPECT_EQ(programHolders, 1u);

  clReleaseMemObject(made.first);
  clReleaseMemObject(made.second);
  clReleaseProgram(made.program);
}

}  // namespace

TEST(opencl_handles, what_was_made_is_let_go_of_at_the_end_of_the_routine)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  Made made;
  {
    OpenCL::Handles handles;
    made = makeAndHandOver(handles);
    EXPECT_EQ(holdersOf(made.first), 2u);
    EXPECT_EQ(holdersOf(made.kernel), 2u);
  }
  expectOnlyOursLeftAndLetGo(made);
}

TEST(opencl_handles, what_was_made_is_let_go_of_when_the_routine_is_thrown_out_of)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  Made made;
  auto routine = [&]()
  {
    OpenCL::Handles handles;
    made = makeAndHandOver(handles);
    throw std::runtime_error("a call on the device failed");
  };
  EXPECT_THROW(routine(), std::runtime_error);
  expectOnlyOursLeftAndLetGo(made);
}

TEST(opencl_handles, a_handle_that_was_never_made_is_left_alone)
{
  // A buffer that failed to be made comes back as nothing, and is kept with the rest all the same.
  OpenCL::Handles handles;
  EXPECT_EQ(handles.keep(nullptr), nullptr);
}
