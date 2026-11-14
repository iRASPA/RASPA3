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
import mc_opencl_void_fraction;

namespace
{
std::uint32_t referenceCount(cl_program program)
{
  cl_uint count = 0;
  EXPECT_EQ(clGetProgramInfo(program, CL_PROGRAM_REFERENCE_COUNT, sizeof(count), &count, nullptr), CL_SUCCESS);
  return count;
}

std::uint32_t referenceCount(cl_kernel kernel)
{
  cl_uint count = 0;
  EXPECT_EQ(clGetKernelInfo(kernel, CL_KERNEL_REFERENCE_COUNT, sizeof(count), &count, nullptr), CL_SUCCESS);
  return count;
}
}  // namespace

// The device report lists the six preferred vector widths under their own names (LONG was printed as FLOAT,
// FLOAT as DOUBLE, DOUBLE not at all).
TEST(opencl_device, device_report_lists_every_preferred_vector_width)
{
  OpenCL::initialize();
  if (!OpenCL::clDeviceId.has_value()) GTEST_SKIP() << "no OpenCL device";

  const std::string report = OpenCL::printBestOpenCLDevice();
  const std::size_t header = report.find("CL_DEVICE_PREFERRED_VECTOR_WIDTH\n");
  ASSERT_NE(header, std::string::npos) << report;
  const std::size_t begin = header + std::string("CL_DEVICE_PREFERRED_VECTOR_WIDTH\n").size();
  const std::string line = report.substr(begin, report.find('\n', begin) - begin);

  // the line itself: "CHAR a, SHORT b, INT c, LONG d, FLOAT e, DOUBLE f" with the values the device reports
  const char* names[6] = {"CHAR", "SHORT", "INT", "LONG", "FLOAT", "DOUBLE"};
  const cl_device_info queries[6] = {CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,  CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,
                                     CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT,   CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,
                                     CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE};
  std::string expected;
  for (std::size_t k = 0; k < 6; ++k)
  {
    cl_uint width = 0;
    ASSERT_EQ(clGetDeviceInfo(OpenCL::clDeviceId.value(), queries[k], sizeof(width), &width, nullptr), CL_SUCCESS);
    expected += std::format("{}{} {}", k == 0 ? "" : ", ", names[k], width);
  }
  EXPECT_EQ(line, expected);
}

// The void-fraction helper initializes OpenCL itself, owns its program and kernel, and releases them once.
TEST(opencl_device, void_fraction_helper_owns_and_releases_its_program_and_kernel)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value() || !OpenCL::clDeviceId.has_value()) GTEST_SKIP() << "no OpenCL device";

  cl_program program = nullptr;
  cl_kernel kernel = nullptr;
  std::uint32_t programReferences = 0;
  std::uint32_t kernelReferences = 0;
  {
    MC_OpenCL_VoidFraction helper{};
    ASSERT_NE(helper.voidFractionProgram, nullptr);
    ASSERT_NE(helper.voidFractionKernel, nullptr);
    EXPECT_GT(helper.voidFractionWorkGroupSize, 0uz);

    // hold a reference of our own to observe the release of the helper's reference (the kernel may hold an
    // implementation-defined reference to its program, so only the drop is checked)
    program = helper.voidFractionProgram;
    kernel = helper.voidFractionKernel;
    ASSERT_EQ(clRetainProgram(program), CL_SUCCESS);
    ASSERT_EQ(clRetainKernel(kernel), CL_SUCCESS);
    programReferences = referenceCount(program);
    kernelReferences = referenceCount(kernel);
    EXPECT_GE(programReferences, 2u);
    EXPECT_EQ(kernelReferences, 2u);
  }
  // the kernel's release also drops the reference it held on the program
  EXPECT_EQ(referenceCount(kernel), kernelReferences - 1u);
  EXPECT_LE(referenceCount(program), programReferences - 1u);
  EXPECT_GE(referenceCount(program), 1u);
  clReleaseKernel(kernel);
  clReleaseProgram(program);

  // constructing and destroying repeatedly leaks nothing and stays valid
  for (int repeat = 0; repeat < 3; ++repeat)
  {
    MC_OpenCL_VoidFraction helper{};
    EXPECT_NE(helper.voidFractionKernel, nullptr);
  }
}
