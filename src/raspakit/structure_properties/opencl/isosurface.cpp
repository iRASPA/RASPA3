module;
     
#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>  
#include <chrono> 
#include <complex>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <exception>
#include <format>
#include <iostream>
#include <fstream>
#include <istream>
#include <map>                                      
#include <ostream>                                  
#include <print>                                    
#include <source_location>
#include <sstream>
#include <tuple>          
#include <utility>
#include <vector>
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
  #include <OpenCL/cl.h>
#elif _WIN32
  #include <CL/cl.h>
#else
  #include <CL/opencl.h>
#endif
#endif
  
module isosurface;

import opencl;
import float4;
import double4;

Isosurface::Isosurface()
{
  if(OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;
    const char* shaderSourceCode = Isosurface::marchingCubesKernel.c_str();
    program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &shaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("Isosurface: OpenCL clCreateProgramWithSource failed at {} line {}\n", __FILE__, __LINE__));}

    // Build the program executable
    err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
      throw std::runtime_error(std::format("Isorface: OpenCL Failed to build program at {}1 (line {} error: {})\n", __FILE__, __LINE__, buffer));
    }

    constructHPLevelKernel = clCreateKernel(program, "constructHPLevel", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateProgramWithSource failed at {} line {}\n", __FILE__, __LINE__));}
    err = clGetKernelWorkGroupInfo(constructHPLevelKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &constructHPLevelKernelWorkGroupSize, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));}

    classifyCubesKernel = clCreateKernel(program, "classifyCubes", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateProgramWithSource failed at {} line {}\n", __FILE__, __LINE__));}
    err = clGetKernelWorkGroupInfo(constructHPLevelKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &classifyCubesKernelWorkGroupSize, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));}

    traverseHPKernel[4] = clCreateKernel(program, "traverseHP16", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));}
    err = clGetKernelWorkGroupInfo(traverseHPKernel[4], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &traverseHPKernelWorkGroupSize[4], nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("cOpenCL lGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));}
    traverseHPKernel[5] = clCreateKernel(program, "traverseHP32", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));}
    err = clGetKernelWorkGroupInfo(traverseHPKernel[5], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &traverseHPKernelWorkGroupSize[5], nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));}
    traverseHPKernel[6] = clCreateKernel(program, "traverseHP64", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));}
    err = clGetKernelWorkGroupInfo(traverseHPKernel[6], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &traverseHPKernelWorkGroupSize[6], nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));}
    traverseHPKernel[7] = clCreateKernel(program, "traverseHP128", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));}
    err = clGetKernelWorkGroupInfo(traverseHPKernel[7], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &traverseHPKernelWorkGroupSize[7], nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));}
    traverseHPKernel[8] = clCreateKernel(program, "traverseHP256", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));}
    err = clGetKernelWorkGroupInfo(traverseHPKernel[8], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &traverseHPKernelWorkGroupSize[8], nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));}
    traverseHPKernel[9] = clCreateKernel(program, "traverseHP512", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));}
    err = clGetKernelWorkGroupInfo(traverseHPKernel[9], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &traverseHPKernelWorkGroupSize[9], nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));}
  }
};

Isosurface::~Isosurface()
{
  if(OpenCL::clContext.has_value())
  {
    clReleaseKernel(traverseHPKernel[9]);
    clReleaseKernel(traverseHPKernel[8]);
    clReleaseKernel(traverseHPKernel[7]);
    clReleaseKernel(traverseHPKernel[6]);
    clReleaseKernel(traverseHPKernel[5]);
    clReleaseKernel(traverseHPKernel[4]);
    clReleaseKernel(classifyCubesKernel);
    clReleaseKernel(constructHPLevelKernel);
    clReleaseProgram(program);
  }
}



std::vector<float4> Isosurface::computeIsosurface(int3 dimensions, std::vector<cl_float> *voxels, double isoValue) noexcept(false)
{
  cl_int err;

  size_t largestSize = static_cast<size_t>(std::max({dimensions.x,dimensions.y,dimensions.z}));
  size_t powerOfTwo = 1uz;
  while(largestSize > static_cast<size_t>(pow(2, powerOfTwo)))
  {
    powerOfTwo += 1;
  }
  size_t size = static_cast<size_t>(std::pow(2uz,powerOfTwo)); // the encompassing size to use as textures (size 16,32,64,128,256, and 512 are supported).

  size_t bufferSize = size;
  std::vector<cl_mem> images;
  std::vector<cl_mem> buffers;

  for(size_t i=1; i< powerOfTwo; i++)
  {
    cl_image_format imageFormat{};
    switch(i)
    {
    case 1:
      imageFormat = cl_image_format{CL_RGBA, CL_UNSIGNED_INT8};
      break;
    case 2:
      imageFormat = cl_image_format{CL_R, CL_UNSIGNED_INT8};
      break;
    case 3:
      imageFormat = cl_image_format{CL_R, CL_UNSIGNED_INT16};
      break;
    case 4:
      imageFormat = cl_image_format{CL_R, CL_UNSIGNED_INT16};
      break;
    default:
      imageFormat = cl_image_format{CL_R, CL_UNSIGNED_INT32};
      break;
    }

    cl_image_desc imageDescriptor = cl_image_desc{CL_MEM_OBJECT_IMAGE3D, bufferSize, bufferSize, bufferSize, 0, 0, 0, 0, 0, nullptr};
    images.push_back(clCreateImage(OpenCL::clContext.value(), CL_MEM_READ_WRITE, &imageFormat, &imageDescriptor, nullptr, &err));
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateImage failed at {} line {}\n", __FILE__, __LINE__));}
    bufferSize /= 2;
  }

  // last one should always be CL_UNSIGNED_INT32 (because the sum is read back as integers)
  cl_image_format imageFormat = cl_image_format{CL_R, CL_UNSIGNED_INT32};
  cl_image_desc imageDescriptor = cl_image_desc{CL_MEM_OBJECT_IMAGE3D, bufferSize, bufferSize, bufferSize, 0, 0, 0, 0, 0, nullptr};
  images.push_back(clCreateImage(OpenCL::clContext.value(), CL_MEM_READ_WRITE, &imageFormat, &imageDescriptor, nullptr, &err));
  if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateImage failed at {} line {}\n", __FILE__, __LINE__));}


  // Transfer dataset to device
  imageFormat = cl_image_format{CL_R, CL_FLOAT};
  imageDescriptor = cl_image_desc{CL_MEM_OBJECT_IMAGE3D, size, size, size, 0, 0, 0, 0, 0, nullptr};
  cl_mem rawData = clCreateImage(OpenCL::clContext.value(), CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, &imageFormat, &imageDescriptor, voxels->data(), nullptr);
  if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateImage failed at {} line {}\n", __FILE__, __LINE__));}


  // update scalar field
  //===================================================================================================================================

  cl_float clIsoValue = cl_float(isoValue);
  cl_int4 clDimensions = {{dimensions.x,dimensions.y,dimensions.z,1}};

  clSetKernelArg(classifyCubesKernel,  0, sizeof(cl_mem), &images[0]);
  clSetKernelArg(classifyCubesKernel,  1, sizeof(cl_mem), &rawData);
  clSetKernelArg(classifyCubesKernel,  2, sizeof(cl_int4), &clDimensions);
  clSetKernelArg(classifyCubesKernel,  3, sizeof(cl_float), &clIsoValue);


  clGetKernelWorkGroupInfo(classifyCubesKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &classifyCubesKernelWorkGroupSize, nullptr);

  // set work-item dimensions
  size_t global_work_size[3] = {size, size, size};

  err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), classifyCubesKernel, 3, nullptr, global_work_size, nullptr, 0, nullptr, nullptr);
  if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel classifyCubesKernel failed at {} line {}\n", __FILE__, __LINE__));}

  // histoPyramidConstruction
  //===================================================================================================================================

  // Run base to first level
  clSetKernelArg(constructHPLevelKernel,  0, sizeof(cl_mem), &images[0]);
  clSetKernelArg(constructHPLevelKernel,  1, sizeof(cl_mem), &images[1]);

  size_t global_work_size2[3] = {size/2, size/2, size/2};

  err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), constructHPLevelKernel, 3, nullptr, global_work_size2, nullptr, 0, nullptr, nullptr);
  if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel constructHPLevelKernel failed at {} line {}\n", __FILE__, __LINE__));}


  size_t previous = size / 2;
  // Run level 2 to top level
  for(size_t i = 1; i < static_cast<size_t>(ceil(log2(double(size))) - 1); i++)
  {
    clSetKernelArg(constructHPLevelKernel,  0, sizeof(cl_mem), &images[i]);
    clSetKernelArg(constructHPLevelKernel,  1, sizeof(cl_mem), &images[i+1]);

    previous /= 2;
    size_t global_work_size3[3] = {previous, previous, previous};
    err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), constructHPLevelKernel, 3, nullptr, global_work_size3, nullptr, 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel constructHPLevelKernel failed at {} line {}\n", __FILE__, __LINE__));}
  }

  // Read top of histoPyramid an use this size to allocate VBO below
  //===================================================================================================================================

  cl_int sum[8] = {0,0,0,0,0,0,0,0};
  size_t origin[3] = {0,0,0};
  size_t region[3] = {2,2,2};

  err = clEnqueueReadImage(OpenCL::clCommandQueue.value(), images[images.size()-1], CL_FALSE, origin, region, 0, 0, &sum, 0, nullptr, nullptr);
  if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueReadImage failed at {} line {}\n", __FILE__, __LINE__));}

  clFinish(OpenCL::clCommandQueue.value());

  cl_int sum2 = sum[0] + sum[1] + sum[2] + sum[3] + sum[4] + sum[5] + sum[6] + sum[7];
  size_t numberOfTriangles = size_t(sum2);


  std::cout << "NUmber of triangles: " << numberOfTriangles << std::endl;

  // get the results and convert them to an OpenGL Vertex buffer object
  //===================================================================================================================================

  std::vector<float4> triangleData(numberOfTriangles * 3 * 3);

  if(numberOfTriangles == 0)
  {
    return triangleData;
  }

  if(numberOfTriangles>0)
  {
    cl_int clSum = cl_int(numberOfTriangles);

    // Increase the global_work_size so that it is divideable by the workgroup-size
    size_t workGroupSize = traverseHPKernelWorkGroupSize[powerOfTwo];
    size_t local_work_size[1] = {size_t(workGroupSize)};
    size_t globalWorkSize = numberOfTriangles + workGroupSize - (numberOfTriangles - workGroupSize * (numberOfTriangles / workGroupSize));
    size_t global_work_size4[1] = {size_t(globalWorkSize)};

    for(size_t j=0; j<images.size(); j++)
    {
      clSetKernelArg(traverseHPKernel[powerOfTwo], static_cast<cl_uint>(j), sizeof(cl_mem), &images[j]);
    }
    clSetKernelArg(traverseHPKernel[powerOfTwo], static_cast<cl_uint>(images.size()), sizeof(cl_mem), &rawData);

    size_t i = images.size() + 1;
    cl_mem  VBOBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, triangleData.size()*sizeof(cl_float4), nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed at {} line {}\n", __FILE__, __LINE__));}
    clSetKernelArg(traverseHPKernel[powerOfTwo], static_cast<cl_uint>(i), sizeof(cl_mem), &VBOBuffer);
    clSetKernelArg(traverseHPKernel[powerOfTwo], static_cast<cl_uint>(i+1), sizeof(cl_int4), &clDimensions);
    clSetKernelArg(traverseHPKernel[powerOfTwo], static_cast<cl_uint>(i+2), sizeof(cl_float), &clIsoValue);
    clSetKernelArg(traverseHPKernel[powerOfTwo], static_cast<cl_uint>(i+3), sizeof(cl_int), &clSum);

    // Run a NDRange kernel over this buffer which traverses back to the base level
    err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), traverseHPKernel[powerOfTwo], 1, nullptr, global_work_size4, local_work_size, 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel traverseHPKernel failed at {} line {}\n", __FILE__, __LINE__));}

    clFinish(OpenCL::clCommandQueue.value());

    err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), VBOBuffer, CL_TRUE, 0, triangleData.size() * sizeof(float4), triangleData.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueReadBuffer failed at {} line {}\n", __FILE__, __LINE__));}

    clFinish(OpenCL::clCommandQueue.value());

    clReleaseMemObject(VBOBuffer);
  }

  clReleaseMemObject(rawData);

  for(cl_mem &image: images)
  {
    clReleaseMemObject(image);
  }

  for(cl_mem &buffer:  buffers)
  {
    clReleaseMemObject(buffer);
  }

  return triangleData;
}



std::string Isosurface::marchingCubesKernel = std::string(R"foo(

#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_NEAREST;


// Cube description:
//         7 ________ 6           _____6__             ________
//         /|       /|         7/|       /|          /|       /|
//       /  |     /  |        /  |     /5 |        /  6     /  |
//   4 /_______ /    |      /__4____ /    10     /_______3/    |
//    |     |  |5    |     |    11  |     |     |     |  |   2 |
//    |    3|__|_____|2    |     |__|__2__|     | 4   |__|_____|
//    |    /   |    /      8   3/   9    /      |    /   |    /
//    |  /     |  /        |  /     |  /1       |  /     5  /
//    |/_______|/          |/___0___|/          |/_1_____|/
//   0          1        0          1
//        Nodes                Borders               Faces


__constant int4 cubeOffsets[8] =
{
  {0, 0, 0, 0},
  {1, 0, 0, 0},
  {0, 0, 1, 0},
  {1, 0, 1, 0},
  {0, 1, 0, 0},
  {1, 1, 0, 0},
  {0, 1, 1, 0},
  {1, 1, 1, 0}
};

__constant char  offsets3[72] =
{
  // 0
  0,0,0,
  1,0,0,
  // 1
  1,0,0,
  1,0,1,
  // 2
  1,0,1,
  0,0,1,
  // 3
  0,0,1,
  0,0,0,
  // 4
  0,1,0,
  1,1,0,
  // 5
  1,1,0,
  1,1,1,
  // 6
  1,1,1,
  0,1,1,
  // 7
  0,1,1,
  0,1,0,
  // 8
  0,0,0,
  0,1,0,
  // 9
  1,0,0,
  1,1,0,
  // 10
  1,0,1,
  1,1,1,
  // 11
  0,0,1,
  0,1,1
};


// Look up table for the number of triangles produced for each of the 256 cases. There at most 5 triangular facets necessary.
__constant uchar numberOfTriangles[256] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 2, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 3, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 3, 2, 3, 3, 2, 3, 4, 4, 3, 3, 4, 4, 3, 4, 5, 5, 2, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 3, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 4, 2, 3, 3, 4, 3, 4, 2, 3, 3, 4, 4, 5, 4, 5, 3, 2, 3, 4, 4, 3, 4, 5, 3, 2, 4, 5, 5, 4, 5, 2, 4, 1, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 3, 2, 3, 3, 4, 3, 4, 4, 5, 3, 2, 4, 3, 4, 3, 5, 2, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 4, 3, 4, 4, 3, 4, 5, 5, 4, 4, 3, 5, 2, 5, 4, 2, 1, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 2, 3, 3, 2, 3, 4, 4, 5, 4, 5, 5, 2, 4, 3, 5, 4, 3, 2, 4, 1, 3, 4, 4, 5, 4, 5, 3, 4, 4, 5, 5, 2, 3, 4, 2, 1, 2, 3, 3, 2, 3, 4, 2, 1, 3, 2, 4, 1, 2, 1, 1, 0};

)foo") +

std::string(R"foo(

// The last part of the algorithm involves forming the correct facets from the positions that the isosurface intersects the edges of the grid cell.
// Again a table (by Cory Gene Bloyd) is used which this time uses the same cubeindex but allows the vertex sequence to be looked up for as many triangular
// facets are necessary to represent the isosurface within the grid cell.
__constant int triTable[4096] =
{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1,
  3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1,
  3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1,
  3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1,
  9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1,
  9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1,
  2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1,
  8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1,
  9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1,
  4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1,
  3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1,
  1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1,
  4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1,
  4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1,
  9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1,
  5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1,
  2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1,
  9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1,
  0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1,
  2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1,
  10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1,
  4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1,
  5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1,
  5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1,
  9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1,
  0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1,
  1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1,
  10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1,
  8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1,
  2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1,
  7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1,
  9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1,
  2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1,
  11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1,
  9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1,
  5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1,
  11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1,
  11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1,
  1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1,
  9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1,
  5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1,
  2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1,
  0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1,
  5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1,
  6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1,
  3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1,
  6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1,
  5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1,
  1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1,
  10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1,
  6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1,
  8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1,
  7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1,
  3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1,
  5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1,
  0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1,
  9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1,
  8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1,
  5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1,
  0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1,
  6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1,
  10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1,
  10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1,
  8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1,
  1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1,
  3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1,
  0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1,
  10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1,
  3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1,
  6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1,
  9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1,
  8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1,
  3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1,
  6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1,
  0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1,
  10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1,
  10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1,
  2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1,
  7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1,
  7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1,
  2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1,
  1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1,
  11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1,
  8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1,
  0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1,
  7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1,
  10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1,
  2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1,
  6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1,
  7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1,
  2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1,
  1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1,
  10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1,
  10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1,
  0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1,
  7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1,
  6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1,
  8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1,
  9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1,
  6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1,
  4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1,
  10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1,
  8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1,
  0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1,
  1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1,
  8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1,
  10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1,
  4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1,
  10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1,
  5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1,
  11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1,
  9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1,
  6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1,
  7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1,
  3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1,
  7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1,
  9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1,
  3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1,
  6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1,
  9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1,
  1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1,
  4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1,
  7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1,
  6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1,
  3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1,
  0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1,
  6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1,
  0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1,
  11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1,
  6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1,
  5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1,
  9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1,
  1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1,
  1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1,
  10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1,
  0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1,
  5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1,
  10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1,
  11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1,
  9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1,
  7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1,
  2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1,
  8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1,
  9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1,
  9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1,
  1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1,
  9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1,
  9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1,
  5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1,
  0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1,
  10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1,
  2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1,
  0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1,
  0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1,
  9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1,
  5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1,
  3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1,
  5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1,
  8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1,
  0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1,
  9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1,
  1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1,
  3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1,
  4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1,
  9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1,
  11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1,
  11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1,
  2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1,
  9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1,
  3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1,
  1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1,
  4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1,
  4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1,
  0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1,
  3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1,
  3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1,
  0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1,
  9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1,
  1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

)foo") +

std::string(R"foo(
__kernel void constructHPLevel(
                               __read_only image3d_t readHistoPyramid,
                               __write_only image3d_t writeHistoPyramid
                               ) {

  int4 writePos = {get_global_id(0), get_global_id(1), get_global_id(2), 0};
  int4 readPos = writePos*2;
  int writeValue = read_imagei(readHistoPyramid, sampler, readPos).x + // 0
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[1])).x + // 1
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[2])).x + // 2
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[3])).x + // 3
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[4])).x + // 4
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[5])).x + // 5
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[6])).x + // 6
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[7])).x;  // 7

  write_imagei(writeHistoPyramid, writePos, writeValue);
}

int4 scanHPLevel(int target, __read_only image3d_t hp, int4 current) {

  int8 neighbors = {
    read_imagei(hp, sampler, current).x,
    read_imagei(hp, sampler, (current + cubeOffsets[1])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[2])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[3])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[4])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[5])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[6])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[7])).x
  };

  int acc = current.s3 + neighbors.s0;
  int8 cmp;
  cmp.s0 = acc <= target;
  acc += neighbors.s1;
  cmp.s1 = acc <= target;
  acc += neighbors.s2;
  cmp.s2 = acc <= target;
  acc += neighbors.s3;
  cmp.s3 = acc <= target;
  acc += neighbors.s4;
  cmp.s4 = acc <= target;
  acc += neighbors.s5;
  cmp.s5 = acc <= target;
  acc += neighbors.s6;
  cmp.s6 = acc <= target;
  cmp.s7 = 0;

  current += cubeOffsets[(cmp.s0+cmp.s1+cmp.s2+cmp.s3+cmp.s4+cmp.s5+cmp.s6+cmp.s7)];
  current.s0 = current.s0*2;
  current.s1 = current.s1*2;
  current.s2 = current.s2*2;
  current.s3 = current.s3 +
  cmp.s0*neighbors.s0 +
  cmp.s1*neighbors.s1 +
  cmp.s2*neighbors.s2 +
  cmp.s3*neighbors.s3 +
  cmp.s4*neighbors.s4 +
  cmp.s5*neighbors.s5 +
  cmp.s6*neighbors.s6 +
  cmp.s7*neighbors.s7;
  return current;
}
)foo") +

std::string(R"foo(
__kernel void traverseHP16(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // compute normal
    const float4 forwardDifference0 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 forwardDifference1 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = forwardDifference0 + (forwardDifference1 - forwardDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}
)foo") +

std::string(R"foo(
__kernel void traverseHP32(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t hp4,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp4, cubePosition);
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // compute normal
    const float4 forwardDifference0 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 forwardDifference1 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = forwardDifference0 + (forwardDifference1 - forwardDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}
)foo") +

std::string(R"foo(
__kernel void traverseHP64(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t hp4,
                         __read_only image3d_t hp5,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp5, cubePosition);
  cubePosition = scanHPLevel(target, hp4, cubePosition);
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // compute normal
    const float4 forwardDifference0 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 forwardDifference1 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = forwardDifference0 + (forwardDifference1 - forwardDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}

)foo") +

std::string(R"foo(
__kernel void traverseHP128(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t hp4,
                         __read_only image3d_t hp5,
                         __read_only image3d_t hp6,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp6, cubePosition);
  cubePosition = scanHPLevel(target, hp5, cubePosition);
  cubePosition = scanHPLevel(target, hp4, cubePosition);
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // compute normal
    const float4 forwardDifference0 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 forwardDifference1 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = forwardDifference0 + (forwardDifference1 - forwardDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}
)foo") +

std::string(R"foo(
__kernel void traverseHP256(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t hp4,
                         __read_only image3d_t hp5,
                         __read_only image3d_t hp6,
                         __read_only image3d_t hp7,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp7, cubePosition);
  cubePosition = scanHPLevel(target, hp6, cubePosition);
  cubePosition = scanHPLevel(target, hp5, cubePosition);
  cubePosition = scanHPLevel(target, hp4, cubePosition);
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // compute normal
    const float4 forwardDifference0 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 forwardDifference1 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = forwardDifference0 + (forwardDifference1 - forwardDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}
)foo") +

std::string(R"foo(
__kernel void traverseHP512(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t hp4,
                         __read_only image3d_t hp5,
                         __read_only image3d_t hp6,
                         __read_only image3d_t hp7,
                         __read_only image3d_t hp8,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp8, cubePosition);
  cubePosition = scanHPLevel(target, hp7, cubePosition);
  cubePosition = scanHPLevel(target, hp6, cubePosition);
  cubePosition = scanHPLevel(target, hp5, cubePosition);
  cubePosition = scanHPLevel(target, hp4, cubePosition);
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // compute normal
    const float4 forwardDifference0 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 forwardDifference1 = (float4)(
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(-read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x+
                                                        read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = forwardDifference0 + (forwardDifference1 - forwardDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}
)foo") +

std::string(R"foo(
// The first part of the algorithm uses a table (edgeTable) which maps the vertices under the isosurface to the intersecting edges.
// An 8 bit index is formed where each bit corresponds to a vertex.
__kernel void classifyCubes(__write_only image3d_t histoPyramid,
                            __read_only image3d_t rawData,
                            __private int4 dimensions,
                            __private float isolevel)
{
  int4 pos = {get_global_id(0), get_global_id(1), get_global_id(2), 0};

  if(any(pos>=dimensions))
  {
    write_imageui(histoPyramid, pos, (uint4)(0, 0, 0, 0));
    return;
  }

  // Find cube class nr
  const float first = read_imagef(rawData, sampler, pos).x;
  const uchar cubeindex =
       ((first > isolevel)) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[1]) % dimensions).x > isolevel) << 1) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[3]) % dimensions).x > isolevel) << 2) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[2]) % dimensions).x > isolevel) << 3) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[4]) % dimensions).x > isolevel) << 4) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[5]) % dimensions).x > isolevel) << 5) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[7]) % dimensions).x > isolevel) << 6) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[6]) % dimensions).x > isolevel) << 7);

  // Store number of triangles
  write_imageui(histoPyramid, pos, (uint4)(numberOfTriangles[cubeindex], cubeindex, first, 0));
}
)foo");
