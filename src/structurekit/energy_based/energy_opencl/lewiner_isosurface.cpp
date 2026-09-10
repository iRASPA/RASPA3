module;

#define CL_TARGET_OPENCL_VERSION 120
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module energy_opencl_lewiner_isosurface;

import std;

import opencl;
import uint3;
import double3;
import marching_cubes;

namespace
{
void check(cl_int err, const char *what, int line)
{
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(
        std::format("LewinerIsosurface: OpenCL {} failed at {} (err {})\n", what, line, err));
  }
}

void emitChars(std::string &out, const char *name, const signed char *data, std::size_t n)
{
  out += std::format("__constant char {}[{}] = {{", name, n);
  for (std::size_t i = 0; i < n; ++i)
  {
    if (i != 0) out += ',';
    if (i % 20 == 0) out += '\n';
    out += std::to_string(static_cast<int>(data[i]));
  }
  out += "};\n";
}

template <typename T>
void emitTable(std::string &out, const char *name, const T &table)
{
  emitChars(out, name, reinterpret_cast<const signed char *>(&table), sizeof(table));
}

std::string tableSource()
{
  std::string out;
  emitTable(out, "lewinerCases", MarchingCubes::cases);
  emitTable(out, "tiling1", MarchingCubes::tiling1);
  emitTable(out, "tiling2", MarchingCubes::tiling2);
  emitTable(out, "test3", MarchingCubes::test3);
  emitTable(out, "tiling3_1", MarchingCubes::tiling3_1);
  emitTable(out, "tiling3_2", MarchingCubes::tiling3_2);
  emitTable(out, "test4", MarchingCubes::test4);
  emitTable(out, "tiling4_1", MarchingCubes::tiling4_1);
  emitTable(out, "tiling4_2", MarchingCubes::tiling4_2);
  emitTable(out, "tiling5", MarchingCubes::tiling5);
  emitTable(out, "test6", MarchingCubes::test6);
  emitTable(out, "tiling6_1_1", MarchingCubes::tiling6_1_1);
  emitTable(out, "tiling6_1_2", MarchingCubes::tiling6_1_2);
  emitTable(out, "tiling6_2", MarchingCubes::tiling6_2);
  emitTable(out, "test7", MarchingCubes::test7);
  emitTable(out, "tiling7_1", MarchingCubes::tiling7_1);
  emitTable(out, "tiling7_2", MarchingCubes::tiling7_2);
  emitTable(out, "tiling7_3", MarchingCubes::tiling7_3);
  emitTable(out, "tiling7_4_1", MarchingCubes::tiling7_4_1);
  emitTable(out, "tiling7_4_2", MarchingCubes::tiling7_4_2);
  emitTable(out, "tiling8", MarchingCubes::tiling8);
  emitTable(out, "tiling9", MarchingCubes::tiling9);
  emitTable(out, "test10", MarchingCubes::test10);
  emitTable(out, "tiling10_1_1", MarchingCubes::tiling10_1_1);
  emitTable(out, "tiling10_1_1_", MarchingCubes::tiling10_1_1_);
  emitTable(out, "tiling10_1_2", MarchingCubes::tiling10_1_2);
  emitTable(out, "tiling10_2", MarchingCubes::tiling10_2);
  emitTable(out, "tiling10_2_", MarchingCubes::tiling10_2_);
  emitTable(out, "tiling11", MarchingCubes::tiling11);
  emitTable(out, "test12", MarchingCubes::test12);
  emitTable(out, "tiling12_1_1", MarchingCubes::tiling12_1_1);
  emitTable(out, "tiling12_1_1_", MarchingCubes::tiling12_1_1_);
  emitTable(out, "tiling12_1_2", MarchingCubes::tiling12_1_2);
  emitTable(out, "tiling12_2", MarchingCubes::tiling12_2);
  emitTable(out, "tiling12_2_", MarchingCubes::tiling12_2_);
  emitTable(out, "test13", MarchingCubes::test13);
  emitTable(out, "subconfig13", MarchingCubes::subconfig13);
  emitTable(out, "tiling13_1", MarchingCubes::tiling13_1);
  emitTable(out, "tiling13_1_", MarchingCubes::tiling13_1_);
  emitTable(out, "tiling13_2", MarchingCubes::tiling13_2);
  emitTable(out, "tiling13_2_", MarchingCubes::tiling13_2_);
  emitTable(out, "tiling13_3", MarchingCubes::tiling13_3);
  emitTable(out, "tiling13_3_", MarchingCubes::tiling13_3_);
  emitTable(out, "tiling13_4", MarchingCubes::tiling13_4);
  emitTable(out, "tiling13_5_1", MarchingCubes::tiling13_5_1);
  emitTable(out, "tiling13_5_2", MarchingCubes::tiling13_5_2);
  emitTable(out, "tiling14", MarchingCubes::tiling14);
  return out;
}

const char *kernelSource = R"CL(

int wrap_index(int a, int n)
{
  int r = a % n;
  return r < 0 ? r + n : r;
}

float sample_field(__global const float *field, int x, int y, int z, int nx, int ny, int nz)
{
  return field[(wrap_index(z, nz) * ny + wrap_index(y, ny)) * nx + wrap_index(x, nx)];
}

int test_face(char face, __private const float *cube)
{
  const int corner_lookup[24] = {
    0, 4, 5, 1,
    1, 5, 6, 2,
    2, 6, 7, 3,
    3, 7, 4, 0,
    0, 3, 2, 1,
    4, 7, 6, 5
  };
  int idx = (face < 0 ? -face : face) - 1;
  float A = cube[corner_lookup[idx * 4 + 0]];
  float B = cube[corner_lookup[idx * 4 + 1]];
  float C = cube[corner_lookup[idx * 4 + 2]];
  float D = cube[corner_lookup[idx * 4 + 3]];
  float det = A * C - B * D;
  if (fabs(det) < 1.0e-7f) return face >= 0;
  return (float)(face) * A * det >= 0.0f;
}

int test_interior(char s, char caseId, char config, char subconfig, __private const float *cube)
{
  float t, At = 0.0f, Bt = 0.0f, Ct = 0.0f, Dt = 0.0f, a, b;
  int test = 0;
  int edge = -1;

  if (caseId == 4 || caseId == 10)
  {
    a = (cube[4] - cube[0]) * (cube[6] - cube[2]) - (cube[7] - cube[3]) * (cube[5] - cube[1]);
    b = cube[2] * (cube[4] - cube[0]) + cube[0] * (cube[6] - cube[2]) - cube[1] * (cube[7] - cube[3]) -
        cube[3] * (cube[5] - cube[1]);
    t = -b / (2.0f * a);
    if (t < 0.0f || t > 1.0f) return s > 0;
    At = cube[0] + (cube[4] - cube[0]) * t;
    Bt = cube[3] + (cube[7] - cube[3]) * t;
    Ct = cube[2] + (cube[6] - cube[2]) * t;
    Dt = cube[1] + (cube[5] - cube[1]) * t;
  }
  else
  {
    if (caseId == 6) edge = test6[config * 3 + 2];
    else if (caseId == 7) edge = test7[config * 5 + 4];
    else if (caseId == 12) edge = test12[config * 4 + 3];
    else if (caseId == 13) edge = tiling13_5_1[(config * 4 + subconfig) * 18];
    else return s > 0;

    if (edge == 0)
    {
      t = cube[0] / (cube[0] - cube[1]);
      At = 0.0f;
      Bt = cube[3] + (cube[2] - cube[3]) * t;
      Ct = cube[7] + (cube[6] - cube[7]) * t;
      Dt = cube[4] + (cube[5] - cube[4]) * t;
    }
    else if (edge == 1)
    {
      t = cube[1] / (cube[1] - cube[2]);
      At = 0.0f;
      Bt = cube[0] + (cube[3] - cube[0]) * t;
      Ct = cube[4] + (cube[7] - cube[4]) * t;
      Dt = cube[5] + (cube[6] - cube[5]) * t;
    }
    else if (edge == 2)
    {
      t = cube[2] / (cube[2] - cube[3]);
      At = 0.0f;
      Bt = cube[1] + (cube[0] - cube[1]) * t;
      Ct = cube[5] + (cube[4] - cube[5]) * t;
      Dt = cube[6] + (cube[7] - cube[6]) * t;
    }
    else if (edge == 3)
    {
      t = cube[3] / (cube[3] - cube[0]);
      At = 0.0f;
      Bt = cube[2] + (cube[1] - cube[2]) * t;
      Ct = cube[6] + (cube[5] - cube[6]) * t;
      Dt = cube[7] + (cube[4] - cube[7]) * t;
    }
    else if (edge == 4)
    {
      t = cube[4] / (cube[4] - cube[5]);
      At = 0.0f;
      Bt = cube[7] + (cube[6] - cube[7]) * t;
      Ct = cube[3] + (cube[2] - cube[3]) * t;
      Dt = cube[0] + (cube[1] - cube[0]) * t;
    }
    else if (edge == 5)
    {
      t = cube[5] / (cube[5] - cube[6]);
      At = 0.0f;
      Bt = cube[4] + (cube[7] - cube[4]) * t;
      Ct = cube[0] + (cube[3] - cube[0]) * t;
      Dt = cube[1] + (cube[2] - cube[1]) * t;
    }
    else if (edge == 6)
    {
      t = cube[6] / (cube[6] - cube[7]);
      At = 0.0f;
      Bt = cube[5] + (cube[4] - cube[5]) * t;
      Ct = cube[1] + (cube[0] - cube[1]) * t;
      Dt = cube[2] + (cube[3] - cube[2]) * t;
    }
    else if (edge == 7)
    {
      t = cube[7] / (cube[7] - cube[4]);
      At = 0.0f;
      Bt = cube[6] + (cube[5] - cube[6]) * t;
      Ct = cube[2] + (cube[1] - cube[2]) * t;
      Dt = cube[3] + (cube[0] - cube[3]) * t;
    }
    else if (edge == 8)
    {
      t = cube[0] / (cube[0] - cube[4]);
      At = 0.0f;
      Bt = cube[3] + (cube[7] - cube[3]) * t;
      Ct = cube[2] + (cube[6] - cube[2]) * t;
      Dt = cube[1] + (cube[5] - cube[1]) * t;
    }
    else if (edge == 9)
    {
      t = cube[1] / (cube[1] - cube[5]);
      At = 0.0f;
      Bt = cube[0] + (cube[4] - cube[0]) * t;
      Ct = cube[3] + (cube[7] - cube[3]) * t;
      Dt = cube[2] + (cube[6] - cube[2]) * t;
    }
    else if (edge == 10)
    {
      t = cube[2] / (cube[2] - cube[6]);
      At = 0.0f;
      Bt = cube[1] + (cube[5] - cube[1]) * t;
      Ct = cube[0] + (cube[4] - cube[0]) * t;
      Dt = cube[3] + (cube[7] - cube[3]) * t;
    }
    else if (edge == 11)
    {
      t = cube[3] / (cube[3] - cube[7]);
      At = 0.0f;
      Bt = cube[2] + (cube[6] - cube[2]) * t;
      Ct = cube[1] + (cube[5] - cube[1]) * t;
      Dt = cube[0] + (cube[4] - cube[0]) * t;
    }
    else return s > 0;
  }

  if (At >= 0.0f) test += 1;
  if (Bt >= 0.0f) test += 2;
  if (Ct >= 0.0f) test += 4;
  if (Dt >= 0.0f) test += 8;

  if (test == 5)
  {
    if (At * Ct - Bt * Dt < 1.0e-7f) return s > 0;
    return s < 0;
  }
  if (test == 10)
  {
    if (At * Ct - Bt * Dt >= 1.0e-7f) return s > 0;
    return s < 0;
  }
  if (test == 7 || test == 11 || test == 13 || test == 14 || test == 15) return s < 0;
  return s > 0;
}

int fill_edges(__private char *edges, __constant char *table, int base, int n)
{
  int count = 3 * n;
  for (int i = 0; i < count; ++i) edges[i] = table[base + i];
  return n;
}

int tessellate(__private const float *cube, __private char *edges)
{
  int lut = 0;
  for (int p = 0; p < 8; ++p)
  {
    if (cube[p] > 0.0f) lut += 1 << p;
  }

  char caseId = lewinerCases[lut * 2];
  char config = lewinerCases[lut * 2 + 1];
  char subconfig = 0;

  if (caseId == 0) return 0;
  if (caseId == 1) return fill_edges(edges, tiling1, config * 3, 1);
  if (caseId == 2) return fill_edges(edges, tiling2, config * 6, 2);
  if (caseId == 3)
  {
    if (test_face(test3[config], cube)) return fill_edges(edges, tiling3_2, config * 12, 4);
    return fill_edges(edges, tiling3_1, config * 6, 2);
  }
  if (caseId == 4)
  {
    if (test_interior(test4[config], caseId, config, subconfig, cube))
      return fill_edges(edges, tiling4_1, config * 6, 2);
    return fill_edges(edges, tiling4_2, config * 18, 6);
  }
  if (caseId == 5) return fill_edges(edges, tiling5, config * 9, 3);
  if (caseId == 6)
  {
    if (test_face(test6[config * 3], cube)) return fill_edges(edges, tiling6_2, config * 15, 5);
    if (test_interior(test6[config * 3 + 1], caseId, config, subconfig, cube))
      return fill_edges(edges, tiling6_1_1, config * 9, 3);
    return fill_edges(edges, tiling6_1_2, config * 27, 9);
  }
  if (caseId == 7)
  {
    if (test_face(test7[config * 5], cube)) subconfig += 1;
    if (test_face(test7[config * 5 + 1], cube)) subconfig += 2;
    if (test_face(test7[config * 5 + 2], cube)) subconfig += 4;
    if (subconfig == 0) return fill_edges(edges, tiling7_1, config * 9, 3);
    if (subconfig == 1) return fill_edges(edges, tiling7_2, (config * 3 + 0) * 15, 5);
    if (subconfig == 2) return fill_edges(edges, tiling7_2, (config * 3 + 1) * 15, 5);
    if (subconfig == 3) return fill_edges(edges, tiling7_3, (config * 3 + 0) * 27, 9);
    if (subconfig == 4) return fill_edges(edges, tiling7_2, (config * 3 + 2) * 15, 5);
    if (subconfig == 5) return fill_edges(edges, tiling7_3, (config * 3 + 1) * 27, 9);
    if (subconfig == 6) return fill_edges(edges, tiling7_3, (config * 3 + 2) * 27, 9);
    if (test_interior(test7[config * 5 + 3], caseId, config, subconfig, cube))
      return fill_edges(edges, tiling7_4_2, config * 27, 9);
    return fill_edges(edges, tiling7_4_1, config * 15, 5);
  }
  if (caseId == 8) return fill_edges(edges, tiling8, config * 6, 2);
  if (caseId == 9) return fill_edges(edges, tiling9, config * 12, 4);
  if (caseId == 10)
  {
    if (test_face(test10[config * 3], cube))
    {
      if (test_face(test10[config * 3 + 1], cube)) return fill_edges(edges, tiling10_1_1_, config * 12, 4);
      return fill_edges(edges, tiling10_2, config * 24, 8);
    }
    if (test_face(test10[config * 3 + 1], cube)) return fill_edges(edges, tiling10_2_, config * 24, 8);
    if (test_interior(test10[config * 3 + 2], caseId, config, subconfig, cube))
      return fill_edges(edges, tiling10_1_1, config * 12, 4);
    return fill_edges(edges, tiling10_1_2, config * 24, 8);
  }
  if (caseId == 11) return fill_edges(edges, tiling11, config * 12, 4);
  if (caseId == 12)
  {
    if (test_face(test12[config * 4], cube))
    {
      if (test_face(test12[config * 4 + 1], cube)) return fill_edges(edges, tiling12_1_1_, config * 12, 4);
      return fill_edges(edges, tiling12_2, config * 24, 8);
    }
    if (test_face(test12[config * 4 + 1], cube)) return fill_edges(edges, tiling12_2_, config * 24, 8);
    if (test_interior(test12[config * 4 + 2], caseId, config, subconfig, cube))
      return fill_edges(edges, tiling12_1_1, config * 12, 4);
    return fill_edges(edges, tiling12_1_2, config * 24, 8);
  }
  if (caseId == 13)
  {
    if (test_face(test13[config * 7], cube)) subconfig += 1;
    if (test_face(test13[config * 7 + 1], cube)) subconfig += 2;
    if (test_face(test13[config * 7 + 2], cube)) subconfig += 4;
    if (test_face(test13[config * 7 + 3], cube)) subconfig += 8;
    if (test_face(test13[config * 7 + 4], cube)) subconfig += 16;
    if (test_face(test13[config * 7 + 5], cube)) subconfig += 32;
    char kind = subconfig13[subconfig];
    if (kind == 0) return fill_edges(edges, tiling13_1, config * 12, 4);
    if (kind >= 1 && kind <= 6) return fill_edges(edges, tiling13_2, (config * 6 + (kind - 1)) * 18, 6);
    if (kind >= 7 && kind <= 18) return fill_edges(edges, tiling13_3, (config * 12 + (kind - 7)) * 30, 10);
    if (kind >= 19 && kind <= 22) return fill_edges(edges, tiling13_4, (config * 4 + (kind - 19)) * 36, 12);
    if (kind >= 23 && kind <= 26)
    {
      char faceSub = (char)(kind - 23);
      if (test_interior(test13[config * 7 + 6], caseId, config, faceSub, cube))
        return fill_edges(edges, tiling13_5_1, (config * 4 + faceSub) * 18, 6);
      return fill_edges(edges, tiling13_5_2, (config * 4 + faceSub) * 30, 10);
    }
    if (kind >= 27 && kind <= 38) return fill_edges(edges, tiling13_3_, (config * 12 + (kind - 27)) * 30, 10);
    if (kind >= 39 && kind <= 44) return fill_edges(edges, tiling13_2_, (config * 6 + (kind - 39)) * 18, 6);
    if (kind == 45) return fill_edges(edges, tiling13_1_, config * 12, 4);
    return 0;
  }
  if (caseId == 14) return fill_edges(edges, tiling14, config * 12, 4);
  return 0;
}

void load_cube(__global const float *field, int i, int j, int k, int nx, int ny, int nz, float iso,
               __private float *cube)
{
  const int ox[8] = {0, 1, 1, 0, 0, 1, 1, 0};
  const int oy[8] = {0, 0, 1, 1, 0, 0, 1, 1};
  const int oz[8] = {0, 0, 0, 0, 1, 1, 1, 1};
  for (int p = 0; p < 8; ++p)
  {
    cube[p] = sample_field(field, i + ox[p], j + oy[p], k + oz[p], nx, ny, nz) - iso;
    if (fabs(cube[p]) < 1.0e-7f) cube[p] = 1.0e-7f;
  }
}

__constant int edgeA[12] = {0, 1, 3, 0, 4, 5, 7, 4, 0, 1, 2, 3};
__constant int edgeB[12] = {1, 2, 2, 3, 5, 6, 6, 7, 4, 5, 6, 7};
__constant int edgeD0[36] = {
  0,0,0,  1,0,0,  0,1,0,  0,0,0,
  0,0,1,  1,0,1,  0,1,1,  0,0,1,
  0,0,0,  1,0,0,  1,1,0,  0,1,0
};
__constant int edgeD1[36] = {
  1,0,0,  1,1,0,  1,1,0,  0,1,0,
  1,0,1,  1,1,1,  1,1,1,  0,1,1,
  0,0,1,  1,0,1,  1,1,1,  0,1,1
};

float3 grid_grad(__global const float *field, int x, int y, int z, int nx, int ny, int nz)
{
  return (float3)(
      0.5f * (sample_field(field, x + 1, y, z, nx, ny, nz) - sample_field(field, x - 1, y, z, nx, ny, nz)),
      0.5f * (sample_field(field, x, y + 1, z, nx, ny, nz) - sample_field(field, x, y - 1, z, nx, ny, nz)),
      0.5f * (sample_field(field, x, y, z + 1, nx, ny, nz) - sample_field(field, x, y, z - 1, nx, ny, nz)));
}

void interpolate_edge(__global const float *field, int i, int j, int k, int nx, int ny, int nz,
                      __private const float *cube, int edge, __private float3 *pos, __private float3 *nrm)
{
  int a = edgeA[edge];
  int b = edgeB[edge];
  float u = cube[a] / (cube[a] - cube[b]);
  float3 d0 = (float3)((float)edgeD0[edge * 3], (float)edgeD0[edge * 3 + 1], (float)edgeD0[edge * 3 + 2]);
  float3 d1 = (float3)((float)edgeD1[edge * 3], (float)edgeD1[edge * 3 + 1], (float)edgeD1[edge * 3 + 2]);
  *pos = (float3)((float)i, (float)j, (float)k) + d0 + (d1 - d0) * u;

  int ax = i + edgeD0[edge * 3];
  int ay = j + edgeD0[edge * 3 + 1];
  int az = k + edgeD0[edge * 3 + 2];
  int bx = i + edgeD1[edge * 3];
  int by = j + edgeD1[edge * 3 + 1];
  int bz = k + edgeD1[edge * 3 + 2];
  *nrm = (1.0f - u) * grid_grad(field, ax, ay, az, nx, ny, nz) +
         u * grid_grad(field, bx, by, bz, nx, ny, nz);
}

void interior_vertex(__global const float *field, int i, int j, int k, int nx, int ny, int nz,
                     __private const float *cube, __private float3 *pos, __private float3 *nrm)
{
  float3 p = (float3)(0.0f);
  float3 n = (float3)(0.0f);
  float weight = 0.0f;
  for (int e = 0; e < 12; ++e)
  {
    if (cube[edgeA[e]] * cube[edgeB[e]] > 0.0f) continue;
    float3 ep, en;
    interpolate_edge(field, i, j, k, nx, ny, nz, cube, e, &ep, &en);
    p += ep;
    n += en;
    weight += 1.0f;
  }
  if (weight > 0.0f)
  {
    p /= weight;
    n /= weight;
  }
  *pos = p;
  *nrm = n;
}

__kernel void classifyLewiner(__global const float *field,
                              __global uint *counts,
                              const int nx, const int ny, const int nz,
                              const float iso,
                              const uint nCubes)
{
  uint gid = get_global_id(0);
  if (gid >= nCubes)
  {
    counts[gid] = 0;
    return;
  }
  int i = (int)(gid % (uint)nx);
  int j = (int)((gid / (uint)nx) % (uint)ny);
  int k = (int)(gid / ((uint)nx * (uint)ny));
  float cube[8];
  char edges[36];
  load_cube(field, i, j, k, nx, ny, nz, iso, cube);
  counts[gid] = (uint)tessellate(cube, edges);
}

#define SCAN_MAX 256

__kernel void exclusiveScanBlocks(__global const uint *input,
                                  __global uint *output,
                                  __global uint *blockSums,
                                  const uint n)
{
  __local uint temp[SCAN_MAX];
  uint lid = get_local_id(0);
  uint gid = get_global_id(0);
  uint ls = get_local_size(0);
  uint val = (gid < n) ? input[gid] : 0;
  temp[lid] = val;
  barrier(CLK_LOCAL_MEM_FENCE);

  for (uint offset = 1; offset < ls; offset <<= 1)
  {
    uint t = temp[lid];
    if (lid >= offset) t += temp[lid - offset];
    barrier(CLK_LOCAL_MEM_FENCE);
    temp[lid] = t;
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if (gid < n) output[gid] = temp[lid] - val;
  if (lid == ls - 1) blockSums[get_group_id(0)] = temp[lid];
}

__kernel void addBlockExclusive(__global uint *output,
                                __global const uint *blockExclusive,
                                const uint n)
{
  uint gid = get_global_id(0);
  if (gid < n) output[gid] += blockExclusive[get_group_id(0)];
}

__kernel void emitLewiner(__global const float *field,
                          __global const uint *offsets,
                          __global float *vertices,
                          const int nx, const int ny, const int nz,
                          const float iso,
                          const uint nCubes)
{
  uint gid = get_global_id(0);
  if (gid >= nCubes) return;

  int i = (int)(gid % (uint)nx);
  int j = (int)((gid / (uint)nx) % (uint)ny);
  int k = (int)(gid / ((uint)nx * (uint)ny));
  float cube[8];
  char edges[36];
  load_cube(field, i, j, k, nx, ny, nz, iso, cube);
  int nTriangles = tessellate(cube, edges);
  if (nTriangles == 0) return;

  float3 interiorPos, interiorNrm;
  int haveInterior = 0;
  uint base = offsets[gid];
  float invx = 1.0f / (float)nx;
  float invy = 1.0f / (float)ny;
  float invz = 1.0f / (float)nz;

  for (int t = 0; t < nTriangles; ++t)
  {
    for (int v = 0; v < 3; ++v)
    {
      int edge = edges[t * 3 + v];
      float3 pos, nrm;
      if (edge == 12)
      {
        if (!haveInterior)
        {
          interior_vertex(field, i, j, k, nx, ny, nz, cube, &interiorPos, &interiorNrm);
          haveInterior = 1;
        }
        pos = interiorPos;
        nrm = interiorNrm;
      }
      else
      {
        interpolate_edge(field, i, j, k, nx, ny, nz, cube, edge, &pos, &nrm);
      }
      float4 scaled = (float4)(pos.x * invx, pos.y * invy, pos.z * invz, 1.0f);
      float4 gradient = (float4)(nrm.x, nrm.y, nrm.z, 0.0f);
      uint slot = (base + (uint)t) * 6 + (uint)v * 2;
      vstore4(scaled, slot, vertices);
      vstore4(gradient, slot + 1, vertices);
    }
  }
}

)CL";

struct LewinerGpu
{
  cl_context context{};
  cl_program program{};
  cl_kernel classify{};
  cl_kernel scan{};
  cl_kernel addBlocks{};
  cl_kernel emit{};
  bool ready{false};

  void release()
  {
    if (!ready) return;
    clReleaseKernel(emit);
    clReleaseKernel(addBlocks);
    clReleaseKernel(scan);
    clReleaseKernel(classify);
    clReleaseProgram(program);
    ready = false;
  }

  LewinerGpu() = default;
  LewinerGpu(const LewinerGpu &) = delete;
  LewinerGpu &operator=(const LewinerGpu &) = delete;
  LewinerGpu(LewinerGpu &&other) noexcept { *this = std::move(other); }
  LewinerGpu &operator=(LewinerGpu &&other) noexcept
  {
    if (this == &other) return *this;
    release();
    context = other.context;
    program = other.program;
    classify = other.classify;
    scan = other.scan;
    addBlocks = other.addBlocks;
    emit = other.emit;
    ready = other.ready;
    other.ready = false;
    return *this;
  }
  ~LewinerGpu() { release(); }
};

LewinerGpu buildGpu()
{
  if (!OpenCL::clContext.has_value() || !OpenCL::clDeviceId.has_value() || !OpenCL::clCommandQueue.has_value())
  {
    throw std::runtime_error("LewinerIsosurface: no OpenCL device found\n");
  }

  std::string source = tableSource() + kernelSource;
  const char *ptr = source.c_str();
  cl_int err = CL_SUCCESS;
  LewinerGpu gpu;
  gpu.program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &ptr, nullptr, &err);
  check(err, "clCreateProgramWithSource", __LINE__);

  err = clBuildProgram(gpu.program, 0, nullptr, nullptr, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    std::size_t length = 0;
    clGetProgramBuildInfo(gpu.program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, 0, nullptr, &length);
    std::string log(length, '\0');
    clGetProgramBuildInfo(gpu.program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, length, log.data(),
                          nullptr);
    clReleaseProgram(gpu.program);
    throw std::runtime_error(std::format("LewinerIsosurface: OpenCL failed to build program (error: {})\n", log));
  }

  gpu.classify = clCreateKernel(gpu.program, "classifyLewiner", &err);
  check(err, "clCreateKernel classifyLewiner", __LINE__);
  gpu.scan = clCreateKernel(gpu.program, "exclusiveScanBlocks", &err);
  check(err, "clCreateKernel exclusiveScanBlocks", __LINE__);
  gpu.addBlocks = clCreateKernel(gpu.program, "addBlockExclusive", &err);
  check(err, "clCreateKernel addBlockExclusive", __LINE__);
  gpu.emit = clCreateKernel(gpu.program, "emitLewiner", &err);
  check(err, "clCreateKernel emitLewiner", __LINE__);
  gpu.context = OpenCL::clContext.value();
  gpu.ready = true;
  return gpu;
}

LewinerGpu &device()
{
  static LewinerGpu gpu;
  if (!gpu.ready || gpu.context != OpenCL::clContext.value())
  {
    gpu.release();
    gpu = buildGpu();
  }
  return gpu;
}

cl_mem buffer(cl_mem_flags flags, std::size_t bytes, const void *host)
{
  cl_mem mem = OpenCL::createBuffer(flags, bytes);
  if (host != nullptr && bytes > 0)
  {
    OpenCL::writeBuffer(mem, bytes, host);
  }
  return mem;
}

void exclusiveScan(LewinerGpu &gpu, cl_mem counts, cl_mem offsets, std::size_t nCubes, cl_uint &total)
{
  std::size_t maxGroup = 256;
  check(clGetKernelWorkGroupInfo(gpu.scan, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(std::size_t),
                                 &maxGroup, nullptr),
        "clGetKernelWorkGroupInfo", __LINE__);
  std::size_t group = 256;
  while (group > maxGroup) group /= 2;
  if (group < 1) group = 1;

  std::size_t nBlocks = (nCubes + group - 1) / group;
  std::size_t global = nBlocks * group;
  cl_uint n = static_cast<cl_uint>(nCubes);

  cl_mem blockSums = buffer(CL_MEM_READ_WRITE, nBlocks * sizeof(cl_uint), nullptr);

  check(clSetKernelArg(gpu.scan, 0, sizeof(cl_mem), &counts), "clSetKernelArg scan", __LINE__);
  check(clSetKernelArg(gpu.scan, 1, sizeof(cl_mem), &offsets), "clSetKernelArg scan", __LINE__);
  check(clSetKernelArg(gpu.scan, 2, sizeof(cl_mem), &blockSums), "clSetKernelArg scan", __LINE__);
  check(clSetKernelArg(gpu.scan, 3, sizeof(cl_uint), &n), "clSetKernelArg scan", __LINE__);
  check(clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), gpu.scan, 1, nullptr, &global, &group, 0, nullptr,
                               nullptr),
        "clEnqueueNDRangeKernel scan", __LINE__);

  std::vector<cl_uint> sums(nBlocks, 0);
  check(clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), blockSums, CL_TRUE, 0, nBlocks * sizeof(cl_uint),
                            sums.data(), 0, nullptr, nullptr),
        "clEnqueueReadBuffer blockSums", __LINE__);

  std::vector<cl_uint> exclusive(nBlocks, 0);
  cl_uint running = 0;
  for (std::size_t i = 0; i < nBlocks; ++i)
  {
    exclusive[i] = running;
    running += sums[i];
  }
  total = running;
  if (total == 0)
  {
    clReleaseMemObject(blockSums);
    return;
  }

  check(clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), blockSums, CL_TRUE, 0, nBlocks * sizeof(cl_uint),
                             exclusive.data(), 0, nullptr, nullptr),
        "clEnqueueWriteBuffer blockExclusive", __LINE__);
  check(clSetKernelArg(gpu.addBlocks, 0, sizeof(cl_mem), &offsets), "clSetKernelArg addBlocks", __LINE__);
  check(clSetKernelArg(gpu.addBlocks, 1, sizeof(cl_mem), &blockSums), "clSetKernelArg addBlocks", __LINE__);
  check(clSetKernelArg(gpu.addBlocks, 2, sizeof(cl_uint), &n), "clSetKernelArg addBlocks", __LINE__);
  check(clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), gpu.addBlocks, 1, nullptr, &global, &group, 0, nullptr,
                               nullptr),
        "clEnqueueNDRangeKernel addBlocks", __LINE__);
  clReleaseMemObject(blockSums);
}

}  // namespace

std::vector<double3> trianglesOfLewinerIsosurface(std::span<const float> field, uint3 gridSize, double isoValue,
                                                  std::vector<double3> *gradients)
{
  const std::size_t nx = gridSize.x;
  const std::size_t ny = gridSize.y;
  const std::size_t nz = gridSize.z;
  const std::size_t nVoxels = nx * ny * nz;
  if (field.size() < nVoxels)
  {
    throw std::runtime_error(std::format(
        "LewinerIsosurface: the field has {} values, too few for a {} x {} x {} grid\n", field.size(), nx, ny, nz));
  }
  if (nVoxels == 0) return {};

  LewinerGpu &gpu = device();
  const std::size_t nCubes = nVoxels;
  const cl_uint nCubesU = static_cast<cl_uint>(nCubes);
  const cl_int cnx = static_cast<cl_int>(nx);
  const cl_int cny = static_cast<cl_int>(ny);
  const cl_int cnz = static_cast<cl_int>(nz);
  const cl_float iso = static_cast<cl_float>(isoValue);

  std::size_t maxGroup = 256;
  check(clGetKernelWorkGroupInfo(gpu.classify, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                 sizeof(std::size_t), &maxGroup, nullptr),
        "clGetKernelWorkGroupInfo classify", __LINE__);
  std::size_t group = std::min(maxGroup, std::size_t{256});
  std::size_t global = ((nCubes + group - 1) / group) * group;

  cl_mem fieldBuf = buffer(CL_MEM_READ_ONLY, nVoxels * sizeof(cl_float), field.data());
  cl_mem counts = buffer(CL_MEM_READ_WRITE, global * sizeof(cl_uint), nullptr);
  cl_mem offsets = buffer(CL_MEM_READ_WRITE, global * sizeof(cl_uint), nullptr);

  check(clSetKernelArg(gpu.classify, 0, sizeof(cl_mem), &fieldBuf), "clSetKernelArg classify", __LINE__);
  check(clSetKernelArg(gpu.classify, 1, sizeof(cl_mem), &counts), "clSetKernelArg classify", __LINE__);
  check(clSetKernelArg(gpu.classify, 2, sizeof(cl_int), &cnx), "clSetKernelArg classify", __LINE__);
  check(clSetKernelArg(gpu.classify, 3, sizeof(cl_int), &cny), "clSetKernelArg classify", __LINE__);
  check(clSetKernelArg(gpu.classify, 4, sizeof(cl_int), &cnz), "clSetKernelArg classify", __LINE__);
  check(clSetKernelArg(gpu.classify, 5, sizeof(cl_float), &iso), "clSetKernelArg classify", __LINE__);
  check(clSetKernelArg(gpu.classify, 6, sizeof(cl_uint), &nCubesU), "clSetKernelArg classify", __LINE__);
  check(clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), gpu.classify, 1, nullptr, &global, &group, 0, nullptr,
                               nullptr),
        "clEnqueueNDRangeKernel classify", __LINE__);

  cl_uint total = 0;
  exclusiveScan(gpu, counts, offsets, nCubes, total);
  if (total == 0)
  {
    clReleaseMemObject(offsets);
    clReleaseMemObject(counts);
    clReleaseMemObject(fieldBuf);
    if (gradients != nullptr) gradients->clear();
    return {};
  }

  cl_mem vertices = buffer(CL_MEM_WRITE_ONLY, static_cast<std::size_t>(total) * 6 * sizeof(cl_float4), nullptr);
  check(clSetKernelArg(gpu.emit, 0, sizeof(cl_mem), &fieldBuf), "clSetKernelArg emit", __LINE__);
  check(clSetKernelArg(gpu.emit, 1, sizeof(cl_mem), &offsets), "clSetKernelArg emit", __LINE__);
  check(clSetKernelArg(gpu.emit, 2, sizeof(cl_mem), &vertices), "clSetKernelArg emit", __LINE__);
  check(clSetKernelArg(gpu.emit, 3, sizeof(cl_int), &cnx), "clSetKernelArg emit", __LINE__);
  check(clSetKernelArg(gpu.emit, 4, sizeof(cl_int), &cny), "clSetKernelArg emit", __LINE__);
  check(clSetKernelArg(gpu.emit, 5, sizeof(cl_int), &cnz), "clSetKernelArg emit", __LINE__);
  check(clSetKernelArg(gpu.emit, 6, sizeof(cl_float), &iso), "clSetKernelArg emit", __LINE__);
  check(clSetKernelArg(gpu.emit, 7, sizeof(cl_uint), &nCubesU), "clSetKernelArg emit", __LINE__);
  check(clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), gpu.emit, 1, nullptr, &global, &group, 0, nullptr,
                               nullptr),
        "clEnqueueNDRangeKernel emit", __LINE__);

  std::vector<float> packed(static_cast<std::size_t>(total) * 6 * 4);
  OpenCL::readBuffer(vertices, packed.size() * sizeof(float), packed.data());

  clReleaseMemObject(vertices);
  clReleaseMemObject(offsets);
  clReleaseMemObject(counts);
  clReleaseMemObject(fieldBuf);

  std::vector<double3> corners;
  corners.reserve(3 * static_cast<std::size_t>(total));
  if (gradients != nullptr)
  {
    gradients->clear();
    gradients->reserve(3 * static_cast<std::size_t>(total));
  }

  for (std::size_t t = 0; t < static_cast<std::size_t>(total); ++t)
  {
    for (std::size_t v = 0; v < 3; ++v)
    {
      const std::size_t base = (t * 6 + v * 2) * 4;
      corners.emplace_back(static_cast<double>(packed[base]), static_cast<double>(packed[base + 1]),
                           static_cast<double>(packed[base + 2]));
      if (gradients != nullptr)
      {
        gradients->emplace_back(static_cast<double>(packed[base + 4]), static_cast<double>(packed[base + 5]),
                                static_cast<double>(packed[base + 6]));
      }
    }
  }
  return corners;
}
