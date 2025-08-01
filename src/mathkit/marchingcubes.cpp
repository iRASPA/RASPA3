/**
 * @file    MarchingCubes.cpp
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.2
 * @date    12/08/2002
 *
 * @brief   MarchingCubes Algorithm with modifications by Samuel H. Wilks <sw463@cam.ac.uk>
 */
//________________________________________________

module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#endif

module marching_cubes;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import double3;

//_____________________________________________________________________________
// print cube for debug
void print_cube([[maybe_unused]] double *cube)
{
  // Rcpp::Rcout << "\t";
  //   for(int i=0; i < 8; ++i) {
  //     Rcpp::Rcout << cube[i] << " ";
  //   }
  //   Rcpp::Rcout << "\n";
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// Constructor
MarchingCubes::MarchingCubes(const int size_x /*= -1*/, const int size_y /*= -1*/, const int size_z /*= -1*/)
    :  //-----------------------------------------------------------------------------
      _originalMC(false),
      _size_x(size_x),
      _size_y(size_y),
      _size_z(size_z)
{
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// main algorithm
void MarchingCubes::run(real iso)
//-----------------------------------------------------------------------------
{
  compute_intersection_points(iso);

  for (_k = 0; _k < _size_z - 1; _k++)
    for (_j = 0; _j < _size_y - 1; _j++)
      for (_i = 0; _i < _size_x - 1; _i++)
      {
        double cube[8];
        _lut_entry = 0;
        for (int p = 0; p < 8; ++p)
        {
          cube[p] = get_data(int3{_i + ((p ^ (p >> 1)) & 1), _j + ((p >> 1) & 1), _k + ((p >> 2) & 1)}) - iso;
          if (std::abs(cube[p]) < std::numeric_limits<double>::epsilon())
          {
            cube[p] = std::numeric_limits<double>::epsilon();
          }
          if (cube[p] > 0) _lut_entry += 1 << p;
        }

        process_cube(cube);
      }
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// init temporary structures (must set sizes before call)
void MarchingCubes::init_temps()
//-----------------------------------------------------------------------------
{
  _data.resize(static_cast<std::size_t>(_size_x * _size_y * _size_z));
  _x_verts.resize(static_cast<std::size_t>(_size_x * _size_y * _size_z));
  _y_verts.resize(static_cast<std::size_t>(_size_x * _size_y * _size_z));
  _z_verts.resize(static_cast<std::size_t>(_size_x * _size_y * _size_z));

  std::memset(_x_verts.data(), -1, _x_verts.size() * sizeof(int));
  std::memset(_y_verts.data(), -1, _y_verts.size() * sizeof(int));
  std::memset(_z_verts.data(), -1, _z_verts.size() * sizeof(int));
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// init all structures (must set sizes before call)
void MarchingCubes::init_all()
//-----------------------------------------------------------------------------
{
  init_temps();
}
//_____________________________________________________________________________

//_____________________________________________________________________________
//_____________________________________________________________________________

//_____________________________________________________________________________
// Compute the intersection points
void MarchingCubes::compute_intersection_points(real iso)
//-----------------------------------------------------------------------------
{
  for (_k = 0; _k < _size_z; ++_k)
  {
    for (_j = 0; _j < _size_y; ++_j)
    {
      for (_i = 0; _i < _size_x; ++_i)
      {
        auto grid_coord = int3{_i, _j, _k};

        double cube[8];
        cube[0] = get_data(grid_coord) - iso;
        if (_i < _size_x - 1)
          cube[1] = get_data(int3{_i + 1, _j, _k}) - iso;
        else
          cube[1] = cube[0];

        if (_j < _size_y - 1)
          cube[3] = get_data(int3{_i, _j + 1, _k}) - iso;
        else
          cube[3] = cube[0];

        if (_k < _size_z - 1)
          cube[4] = get_data(int3{_i, _j, _k + 1}) - iso;
        else
          cube[4] = cube[0];

        if (std::abs(cube[0]) < std::numeric_limits<double>::epsilon())
          cube[0] = std::numeric_limits<double>::epsilon();
        if (std::abs(cube[1]) < std::numeric_limits<double>::epsilon())
          cube[1] = std::numeric_limits<double>::epsilon();
        if (std::abs(cube[3]) < std::numeric_limits<double>::epsilon())
          cube[3] = std::numeric_limits<double>::epsilon();
        if (std::abs(cube[4]) < std::numeric_limits<double>::epsilon())
          cube[4] = std::numeric_limits<double>::epsilon();

        if (cube[0] < 0)
        {
          if (cube[1] > 0) set_x_vert(add_vertex(grid_coord, int3{1, 0, 0}, 1, cube), _i, _j, _k);
          if (cube[3] > 0) set_y_vert(add_vertex(grid_coord, int3{0, 1, 0}, 3, cube), _i, _j, _k);
          if (cube[4] > 0) set_z_vert(add_vertex(grid_coord, int3{0, 0, 1}, 4, cube), _i, _j, _k);
        }
        else
        {
          if (cube[1] < 0) set_x_vert(add_vertex(grid_coord, int3{1, 0, 0}, 1, cube), _i, _j, _k);
          if (cube[3] < 0) set_y_vert(add_vertex(grid_coord, int3{0, 1, 0}, 3, cube), _i, _j, _k);
          if (cube[4] < 0) set_z_vert(add_vertex(grid_coord, int3{0, 0, 1}, 4, cube), _i, _j, _k);
        }
      }
    }
  }
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// tests if the components of the tesselation of the cube should be connected by the interior of an ambiguous face
// Test a face
// if face>0 return true if the face contains a part of the surface
bool MarchingCubes::test_face(schar face, double *cube)
{
  static int corner_lookup[6][4] = {{0, 4, 5, 1}, {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}, {0, 3, 2, 1}, {4, 7, 6, 5}};

  int idx = std::abs(face) - 1;
  auto corners = corner_lookup[idx];
  auto A = cube[corners[0]];
  auto B = cube[corners[1]];
  auto C = cube[corners[2]];
  auto D = cube[corners[3]];

  if (std::abs(A * C - B * D) < std::numeric_limits<double>::epsilon())
  {
    return face >= 0;
  }

  // face and A invert signs
  return face * A * (A * C - B * D) >= 0.0;
}

//_____________________________________________________________________________
// Test the interior of a cube
// if s == 7, return true  if the interior is empty
// if s ==-7, return false if the interior is empty
bool MarchingCubes::test_interior(schar s, double *cube)
//-----------------------------------------------------------------------------
{
  real t, At = 0, Bt = 0, Ct = 0, Dt = 0, a, b;
  char test = 0;
  char edge = -1;  // reference edge of the triangulation

  switch (_case)
  {
    case 4:
    case 10:
      a = (cube[4] - cube[0]) * (cube[6] - cube[2]) - (cube[7] - cube[3]) * (cube[5] - cube[1]);
      b = cube[2] * (cube[4] - cube[0]) + cube[0] * (cube[6] - cube[2]) - cube[1] * (cube[7] - cube[3]) -
          cube[3] * (cube[5] - cube[1]);
      t = -b / (2 * a);
      if (t < 0 || t > 1) return s > 0;

      At = cube[0] + (cube[4] - cube[0]) * t;
      Bt = cube[3] + (cube[7] - cube[3]) * t;
      Ct = cube[2] + (cube[6] - cube[2]) * t;
      Dt = cube[1] + (cube[5] - cube[1]) * t;
      break;

    case 6:
    case 7:
    case 12:
    case 13:
      switch (_case)
      {
        case 6:
          edge = test6[_config][2];
          break;
        case 7:
          edge = test7[_config][4];
          break;
        case 12:
          edge = test12[_config][3];
          break;
        case 13:
          edge = tiling13_5_1[_config][_subconfig][0];
          break;
      }
      switch (edge)
      {
        case 0:
          t = cube[0] / (cube[0] - cube[1]);
          At = 0;
          Bt = cube[3] + (cube[2] - cube[3]) * t;
          Ct = cube[7] + (cube[6] - cube[7]) * t;
          Dt = cube[4] + (cube[5] - cube[4]) * t;
          break;
        case 1:
          t = cube[1] / (cube[1] - cube[2]);
          At = 0;
          Bt = cube[0] + (cube[3] - cube[0]) * t;
          Ct = cube[4] + (cube[7] - cube[4]) * t;
          Dt = cube[5] + (cube[6] - cube[5]) * t;
          break;
        case 2:
          t = cube[2] / (cube[2] - cube[3]);
          At = 0;
          Bt = cube[1] + (cube[0] - cube[1]) * t;
          Ct = cube[5] + (cube[4] - cube[5]) * t;
          Dt = cube[6] + (cube[7] - cube[6]) * t;
          break;
        case 3:
          t = cube[3] / (cube[3] - cube[0]);
          At = 0;
          Bt = cube[2] + (cube[1] - cube[2]) * t;
          Ct = cube[6] + (cube[5] - cube[6]) * t;
          Dt = cube[7] + (cube[4] - cube[7]) * t;
          break;
        case 4:
          t = cube[4] / (cube[4] - cube[5]);
          At = 0;
          Bt = cube[7] + (cube[6] - cube[7]) * t;
          Ct = cube[3] + (cube[2] - cube[3]) * t;
          Dt = cube[0] + (cube[1] - cube[0]) * t;
          break;
        case 5:
          t = cube[5] / (cube[5] - cube[6]);
          At = 0;
          Bt = cube[4] + (cube[7] - cube[4]) * t;
          Ct = cube[0] + (cube[3] - cube[0]) * t;
          Dt = cube[1] + (cube[2] - cube[1]) * t;
          break;
        case 6:
          t = cube[6] / (cube[6] - cube[7]);
          At = 0;
          Bt = cube[5] + (cube[4] - cube[5]) * t;
          Ct = cube[1] + (cube[0] - cube[1]) * t;
          Dt = cube[2] + (cube[3] - cube[2]) * t;
          break;
        case 7:
          t = cube[7] / (cube[7] - cube[4]);
          At = 0;
          Bt = cube[6] + (cube[5] - cube[6]) * t;
          Ct = cube[2] + (cube[1] - cube[2]) * t;
          Dt = cube[3] + (cube[0] - cube[3]) * t;
          break;
        case 8:
          t = cube[0] / (cube[0] - cube[4]);
          At = 0;
          Bt = cube[3] + (cube[7] - cube[3]) * t;
          Ct = cube[2] + (cube[6] - cube[2]) * t;
          Dt = cube[1] + (cube[5] - cube[1]) * t;
          break;
        case 9:
          t = cube[1] / (cube[1] - cube[5]);
          At = 0;
          Bt = cube[0] + (cube[4] - cube[0]) * t;
          Ct = cube[3] + (cube[7] - cube[3]) * t;
          Dt = cube[2] + (cube[6] - cube[2]) * t;
          break;
        case 10:
          t = cube[2] / (cube[2] - cube[6]);
          At = 0;
          Bt = cube[1] + (cube[5] - cube[1]) * t;
          Ct = cube[0] + (cube[4] - cube[0]) * t;
          Dt = cube[3] + (cube[7] - cube[3]) * t;
          break;
        case 11:
          t = cube[3] / (cube[3] - cube[7]);
          At = 0;
          Bt = cube[2] + (cube[6] - cube[2]) * t;
          Ct = cube[1] + (cube[5] - cube[1]) * t;
          Dt = cube[0] + (cube[4] - cube[0]) * t;
          break;
        default:
          std::cout << " Invalid edge " << edge << "\n";
          print_cube(cube);
          break;
      }
      break;

    default:
      std::cout << " Invalid ambiguous case " << _case << "\n";
      print_cube(cube);
      break;
  }

  if (At >= 0) test++;
  if (Bt >= 0) test += 2;
  if (Ct >= 0) test += 4;
  if (Dt >= 0) test += 8;
  switch (test)
  {
    case 0:
      return s > 0;
    case 1:
      return s > 0;
    case 2:
      return s > 0;
    case 3:
      return s > 0;
    case 4:
      return s > 0;
    case 5:
      if (At * Ct - Bt * Dt < std::numeric_limits<double>::epsilon()) return s > 0;
      break;
    case 6:
      return s > 0;
    case 7:
      return s < 0;
    case 8:
      return s > 0;
    case 9:
      return s > 0;
    case 10:
      if (At * Ct - Bt * Dt >= std::numeric_limits<double>::epsilon()) return s > 0;
      break;
    case 11:
      return s < 0;
    case 12:
      return s > 0;
    case 13:
      return s < 0;
    case 14:
      return s < 0;
    case 15:
      return s < 0;
  }

  return s < 0;
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// Process a unit cube
void MarchingCubes::process_cube(double *cube)
//-----------------------------------------------------------------------------
{
  if (_originalMC)
  {
    char nt = 0;
    while (casesClassic[_lut_entry][3 * nt] != -1) nt++;
    add_triangle(casesClassic[_lut_entry], nt);
    return;
  }

  int v12 = -1;
  _case = static_cast<uchar>(cases[_lut_entry][0]);
  _config = static_cast<uchar>(cases[_lut_entry][1]);
  _subconfig = 0;

  switch (_case)
  {
    case 0:
      break;

    case 1:
      add_triangle(tiling1[_config], 1);
      break;

    case 2:
      add_triangle(tiling2[_config], 2);
      break;

    case 3:
      if (test_face(test3[_config], cube))
        add_triangle(tiling3_2[_config], 4);  // 3.2
      else
        add_triangle(tiling3_1[_config], 2);  // 3.1
      break;

    case 4:
      if (test_interior(test4[_config], cube))
        add_triangle(tiling4_1[_config], 2);  // 4.1.1
      else
        add_triangle(tiling4_2[_config], 6);  // 4.1.2
      break;

    case 5:
      add_triangle(tiling5[_config], 3);
      break;

    case 6:
      if (test_face(test6[_config][0], cube))
        add_triangle(tiling6_2[_config], 5);  // 6.2
      else
      {
        if (test_interior(test6[_config][1], cube))
          add_triangle(tiling6_1_1[_config], 3);  // 6.1.1
        else
        {
          v12 = add_c_vertex();
          add_triangle(tiling6_1_2[_config], 9, v12);  // 6.1.2
        }
      }
      break;

    case 7:
      if (test_face(test7[_config][0], cube)) _subconfig += 1;
      if (test_face(test7[_config][1], cube)) _subconfig += 2;
      if (test_face(test7[_config][2], cube)) _subconfig += 4;
      switch (_subconfig)
      {
        case 0:
          add_triangle(tiling7_1[_config], 3);
          break;
        case 1:
          add_triangle(tiling7_2[_config][0], 5);
          break;
        case 2:
          add_triangle(tiling7_2[_config][1], 5);
          break;
        case 3:
          v12 = add_c_vertex();
          add_triangle(tiling7_3[_config][0], 9, v12);
          break;
        case 4:
          add_triangle(tiling7_2[_config][2], 5);
          break;
        case 5:
          v12 = add_c_vertex();
          add_triangle(tiling7_3[_config][1], 9, v12);
          break;
        case 6:
          v12 = add_c_vertex();
          add_triangle(tiling7_3[_config][2], 9, v12);
          break;
        case 7:
          if (test_interior(test7[_config][3], cube))
            add_triangle(tiling7_4_2[_config], 9);
          else
            add_triangle(tiling7_4_1[_config], 5);
          break;
      };
      break;

    case 8:
      add_triangle(tiling8[_config], 2);
      break;

    case 9:
      add_triangle(tiling9[_config], 4);
      break;

    case 10:
      if (test_face(test10[_config][0], cube))
      {
        if (test_face(test10[_config][1], cube))
          add_triangle(tiling10_1_1_[_config], 4);  // 10.1.1
        else
        {
          v12 = add_c_vertex();
          add_triangle(tiling10_2[_config], 8, v12);  // 10.2
        }
      }
      else
      {
        if (test_face(test10[_config][1], cube))
        {
          v12 = add_c_vertex();
          add_triangle(tiling10_2_[_config], 8, v12);  // 10.2
        }
        else
        {
          if (test_interior(test10[_config][2], cube))
            add_triangle(tiling10_1_1[_config], 4);  // 10.1.1
          else
            add_triangle(tiling10_1_2[_config], 8);  // 10.1.2
        }
      }
      break;

    case 11:
      add_triangle(tiling11[_config], 4);
      break;

    case 12:
      if (test_face(test12[_config][0], cube))
      {
        if (test_face(test12[_config][1], cube))
          add_triangle(tiling12_1_1_[_config], 4);  // 12.1.1
        else
        {
          v12 = add_c_vertex();
          add_triangle(tiling12_2[_config], 8, v12);  // 12.2
        }
      }
      else
      {
        if (test_face(test12[_config][1], cube))
        {
          v12 = add_c_vertex();
          add_triangle(tiling12_2_[_config], 8, v12);  // 12.2
        }
        else
        {
          if (test_interior(test12[_config][2], cube))
            add_triangle(tiling12_1_1[_config], 4);  // 12.1.1
          else
            add_triangle(tiling12_1_2[_config], 8);  // 12.1.2
        }
      }
      break;

    case 13:
      if (test_face(test13[_config][0], cube)) _subconfig += 1;
      if (test_face(test13[_config][1], cube)) _subconfig += 2;
      if (test_face(test13[_config][2], cube)) _subconfig += 4;
      if (test_face(test13[_config][3], cube)) _subconfig += 8;
      if (test_face(test13[_config][4], cube)) _subconfig += 16;
      if (test_face(test13[_config][5], cube)) _subconfig += 32;
      switch (subconfig13[_subconfig])
      {
        case 0: /* 13.1 */
          add_triangle(tiling13_1[_config], 4);
          break;

        case 1: /* 13.2 */
          add_triangle(tiling13_2[_config][0], 6);
          break;
        case 2: /* 13.2 */
          add_triangle(tiling13_2[_config][1], 6);
          break;
        case 3: /* 13.2 */
          add_triangle(tiling13_2[_config][2], 6);
          break;
        case 4: /* 13.2 */
          add_triangle(tiling13_2[_config][3], 6);
          break;
        case 5: /* 13.2 */
          add_triangle(tiling13_2[_config][4], 6);
          break;
        case 6: /* 13.2 */
          add_triangle(tiling13_2[_config][5], 6);
          break;

        case 7: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][0], 10, v12);
          break;
        case 8: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][1], 10, v12);
          break;
        case 9: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][2], 10, v12);
          break;
        case 10: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][3], 10, v12);
          break;
        case 11: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][4], 10, v12);
          break;
        case 12: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][5], 10, v12);
          break;
        case 13: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][6], 10, v12);
          break;
        case 14: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][7], 10, v12);
          break;
        case 15: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][8], 10, v12);
          break;
        case 16: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][9], 10, v12);
          break;
        case 17: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][10], 10, v12);
          break;
        case 18: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3[_config][11], 10, v12);
          break;

        case 19: /* 13.4 */
          v12 = add_c_vertex();
          add_triangle(tiling13_4[_config][0], 12, v12);
          break;
        case 20: /* 13.4 */
          v12 = add_c_vertex();
          add_triangle(tiling13_4[_config][1], 12, v12);
          break;
        case 21: /* 13.4 */
          v12 = add_c_vertex();
          add_triangle(tiling13_4[_config][2], 12, v12);
          break;
        case 22: /* 13.4 */
          v12 = add_c_vertex();
          add_triangle(tiling13_4[_config][3], 12, v12);
          break;

        case 23: /* 13.5 */
          _subconfig = 0;
          if (test_interior(test13[_config][6], cube))
            add_triangle(tiling13_5_1[_config][0], 6);
          else
            add_triangle(tiling13_5_2[_config][0], 10);
          break;
        case 24: /* 13.5 */
          _subconfig = 1;
          if (test_interior(test13[_config][6], cube))
            add_triangle(tiling13_5_1[_config][1], 6);
          else
            add_triangle(tiling13_5_2[_config][1], 10);
          break;
        case 25: /* 13.5 */
          _subconfig = 2;
          if (test_interior(test13[_config][6], cube))
            add_triangle(tiling13_5_1[_config][2], 6);
          else
            add_triangle(tiling13_5_2[_config][2], 10);
          break;
        case 26: /* 13.5 */
          _subconfig = 3;
          if (test_interior(test13[_config][6], cube))
            add_triangle(tiling13_5_1[_config][3], 6);
          else
            add_triangle(tiling13_5_2[_config][3], 10);
          break;

        case 27: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][0], 10, v12);
          break;
        case 28: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][1], 10, v12);
          break;
        case 29: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][2], 10, v12);
          break;
        case 30: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][3], 10, v12);
          break;
        case 31: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][4], 10, v12);
          break;
        case 32: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][5], 10, v12);
          break;
        case 33: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][6], 10, v12);
          break;
        case 34: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][7], 10, v12);
          break;
        case 35: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][8], 10, v12);
          break;
        case 36: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][9], 10, v12);
          break;
        case 37: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][10], 10, v12);
          break;
        case 38: /* 13.3 */
          v12 = add_c_vertex();
          add_triangle(tiling13_3_[_config][11], 10, v12);
          break;

        case 39: /* 13.2 */
          add_triangle(tiling13_2_[_config][0], 6);
          break;
        case 40: /* 13.2 */
          add_triangle(tiling13_2_[_config][1], 6);
          break;
        case 41: /* 13.2 */
          add_triangle(tiling13_2_[_config][2], 6);
          break;
        case 42: /* 13.2 */
          add_triangle(tiling13_2_[_config][3], 6);
          break;
        case 43: /* 13.2 */
          add_triangle(tiling13_2_[_config][4], 6);
          break;
        case 44: /* 13.2 */
          add_triangle(tiling13_2_[_config][5], 6);
          break;

        case 45: /* 13.1 */
          add_triangle(tiling13_1_[_config], 4);
          break;

        default:
          std::cout << "Marching Cubes: Impossible case 13?\n";
          print_cube(cube);
      }
      break;

    case 14:
      add_triangle(tiling14[_config], 4);
      break;
  };
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// Adding triangles
void MarchingCubes::add_triangle(const schar *trig, char n, int v12)
{
  int i = 0;
  while (i < 3 * n)
  {
    int tv[3];

    for (int t = 0; t < 3; ++t, ++i)
    {
      switch (trig[i])
      {
        case 0:
          tv[t] = get_x_vert(_i, _j, _k);
          break;
        case 1:
          tv[t] = get_y_vert(_i + 1, _j, _k);
          break;
        case 2:
          tv[t] = get_x_vert(_i, _j + 1, _k);
          break;
        case 3:
          tv[t] = get_y_vert(_i, _j, _k);
          break;
        case 4:
          tv[t] = get_x_vert(_i, _j, _k + 1);
          break;
        case 5:
          tv[t] = get_y_vert(_i + 1, _j, _k + 1);
          break;
        case 6:
          tv[t] = get_x_vert(_i, _j + 1, _k + 1);
          break;
        case 7:
          tv[t] = get_y_vert(_i, _j, _k + 1);
          break;
        case 8:
          tv[t] = get_z_vert(_i, _j, _k);
          break;
        case 9:
          tv[t] = get_z_vert(_i + 1, _j, _k);
          break;
        case 10:
          tv[t] = get_z_vert(_i + 1, _j + 1, _k);
          break;
        case 11:
          tv[t] = get_z_vert(_i, _j + 1, _k);
          break;
        case 12:
          tv[t] = v12;
          break;
        default:
          break;
      }

      // left operand is garbage
      // if( tv[t] == -1 ) {
      //  std::cout << "Marching Cubes: invalid triangle " << (ntrigs() + 1) << "\n";
      //    //print_cube() ;
      //}
    }

    _triangles.push_back(Triangle{tv[0], tv[1], tv[2]});
  }
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// Calculating gradient

real MarchingCubes::get_x_grad(const int i, const int j, const int k) const
//-----------------------------------------------------------------------------
{
  if (i > 0)
  {
    if (i < _size_x - 1)
    {
      return (get_data(int3{i + 1, j, k}) - get_data(int3{i - 1, j, k})) / 2;
    }
    else
    {
      return get_data(int3{i, j, k}) - get_data(int3{i - 1, j, k});
    }
  }
  else
  {
    return get_data(int3{i + 1, j, k}) - get_data(int3{i, j, k});
  }
}
//-----------------------------------------------------------------------------

real MarchingCubes::get_y_grad(const int i, const int j, const int k) const
//-----------------------------------------------------------------------------
{
  if (j > 0)
  {
    if (j < _size_y - 1)
    {
      return (get_data(int3{i, j + 1, k}) - get_data(int3{i, j - 1, k})) / 2;
    }
    else
    {
      return get_data(int3{i, j, k}) - get_data(int3{i, j - 1, k});
    }
  }
  else
  {
    return get_data(int3{i, j + 1, k}) - get_data(int3{i, j, k});
  }
}
//-----------------------------------------------------------------------------

real MarchingCubes::get_z_grad(const int i, const int j, const int k) const
//-----------------------------------------------------------------------------
{
  if (k > 0)
  {
    if (k < _size_z - 1)
    {
      return (get_data(int3{i, j, k + 1}) - get_data(int3{i, j, k - 1})) / 2;
    }
    else
    {
      return get_data(int3{i, j, k}) - get_data(int3{i, j, k - 1});
    }
  }
  else
  {
    return get_data(int3{i, j, k + 1}) - get_data(int3{i, j, k});
  }
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// Adding vertices

int MarchingCubes::add_vertex(const int3 &grid_coord, const int3 &dir, int corner, double *cube)
{
  double u = cube[0] / (cube[0] - cube[corner]);

  double grid_coord_x = grid_coord.x;
  double grid_coord_y = grid_coord.y;
  double grid_coord_z = grid_coord.z;

  double dir_x = dir.x;
  double dir_y = dir.y;
  double dir_z = dir.z;

  double3 pos = double3{grid_coord_x, grid_coord_y, grid_coord_z} + double3{dir_x, dir_y, dir_z} * u;

  int3 grid_coord2 = grid_coord + dir;

  int grid_coord2_x = grid_coord2.x;
  int grid_coord2_y = grid_coord2.y;
  int grid_coord2_z = grid_coord2.z;

  double nx = (1.0 - u) * get_x_grad(static_cast<int>(grid_coord_x), static_cast<int>(grid_coord_y),
                                     static_cast<int>(grid_coord_z)) +
              u * get_x_grad(grid_coord2_x, grid_coord2_y, grid_coord2_z);
  double ny = (1.0 - u) * get_y_grad(static_cast<int>(grid_coord_x), static_cast<int>(grid_coord_y),
                                     static_cast<int>(grid_coord_z)) +
              u * get_y_grad(grid_coord2_x, grid_coord2_y, grid_coord2_z);
  double nz = (1.0 - u) * get_z_grad(static_cast<int>(grid_coord_x), static_cast<int>(grid_coord_y),
                                     static_cast<int>(grid_coord_z)) +
              u * get_z_grad(grid_coord2_x, grid_coord2_y, grid_coord2_z);

  // double3 n = arma::normalise(double3{nx, ny, nz});
  double3 n = double3{nx, ny, nz}.normalized();

  double pos_x = pos.x;
  double pos_y = pos.y;
  double pos_z = pos.z;

  double n_x = n.x;
  double n_y = n.y;
  double n_z = n.z;

  _vertices.push_back(Vertex{pos_x, pos_y, pos_z, n_x, n_y, n_z});
  return static_cast<int>(_vertices.size() - 1);
}

int MarchingCubes::add_c_vertex()
//-----------------------------------------------------------------------------
{
  auto u = double{0.0};
  auto pos = double3(0, 0, 0);
  auto n = double3(0, 0, 0);

  // Computes the average of the intersection points of the cube
  // x-face
  for (auto t : {0, 1})
  {
    for (auto s : {0, 1})
    {
      auto vid = get_x_vert(_i, _j + s, _k + t);
      if (vid != -1)
      {
        ++u;
        const Vertex &v = _vertices[static_cast<std::size_t>(vid)];
        pos += double3{v.x, v.y, v.z};
        n += double3{v.nx, v.ny, v.nz};
      }
    }
  }

  // y-face
  for (auto t : {0, 1})
  {
    for (auto s : {0, 1})
    {
      auto vid = get_y_vert(_i + t, _j, _k + s);
      if (vid != -1)
      {
        ++u;
        const Vertex &v = _vertices[static_cast<std::size_t>(vid)];
        pos += double3{v.x, v.y, v.z};
        n += double3{v.nx, v.ny, v.nz};
      }
    }
  }

  // z-face
  for (auto t : {0, 1})
  {
    for (auto s : {0, 1})
    {
      auto vid = get_z_vert(_i + s, _j + t, _k);
      if (vid != -1)
      {
        ++u;
        const Vertex &v = _vertices[static_cast<std::size_t>(vid)];
        pos += double3{v.x, v.y, v.z};
        n += double3{v.nx, v.ny, v.nz};
      }
    }
  }

  pos.x *= 1.0 / u;
  pos.y *= 1.0 / u;
  pos.z *= 1.0 / u;
  // n = arma::normalise(n);
  n = n.normalized();

  double pos_x = pos.x;
  double pos_y = pos.y;
  double pos_z = pos.z;

  double n_x = n.x;
  double n_y = n.y;
  double n_z = n.z;

  _vertices.push_back(Vertex{pos_x, pos_y, pos_z, n_x, n_y, n_z});
  return static_cast<int>(_vertices.size() - 1);
}
//_____________________________________________________________________________
//

//_____________________________________________________________________________
/**
 * \brief case mapping
 * For each of the possible vertex states listed in this table there is a
 * specific triangulation of the edge intersection points.  The table lists
 * all of them in the form of 0-5 edge triples with the list terminated by
 * the invalid value -1.  For example: case[3] list the 2 triangles
 * formed when cube[0] and cube[1] are inside of the surface, but the rest of
 * the cube is not.
 *
 * Cube description:
 *         7 ________ 6           _____6__             ________
 *         /|       /|         7/|       /|          /|       /|
 *       /  |     /  |        /  |     /5 |        /  6     /  |
 *   4 /_______ /    |      /__4____ /    10     /_______3/    |
 *    |     |  |5    |     |    11  |     |     |     |  |   2 |
 *    |    3|__|_____|2    |     |__|__2__|     | 4   |__|_____|
 *    |    /   |    /      8   3/   9    /      |    /   |    /
 *    |  /     |  /        |  /     |  /1       |  /     5  /
 *    |/_______|/          |/___0___|/          |/_1_____|/
 *   0          1        0          1
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::cases[256][2] = {
    /*   0:                          */ {0, -1},
    /*   1: 0,                       */ {1, 0},
    /*   2:    1,                    */ {1, 1},
    /*   3: 0, 1,                    */ {2, 0},
    /*   4:       2,                 */ {1, 2},
    /*   5: 0,    2,                 */ {3, 0},
    /*   6:    1, 2,                 */ {2, 3},
    /*   7: 0, 1, 2,                 */ {5, 0},
    /*   8:          3,              */ {1, 3},
    /*   9: 0,       3,              */ {2, 1},
    /*  10:    1,    3,              */ {3, 3},
    /*  11: 0, 1,    3,              */ {5, 1},
    /*  12:       2, 3,              */ {2, 5},
    /*  13: 0,    2, 3,              */ {5, 4},
    /*  14:    1, 2, 3,              */ {5, 9},
    /*  15: 0, 1, 2, 3,              */ {8, 0},
    /*  16:             4,           */ {1, 4},
    /*  17: 0,          4,           */ {2, 2},
    /*  18:    1,       4,           */ {3, 4},
    /*  19: 0, 1,       4,           */ {5, 2},
    /*  20:       2,    4,           */ {4, 2},
    /*  21: 0,    2,    4,           */ {6, 2},
    /*  22:    1, 2,    4,           */ {6, 9},
    /*  23: 0, 1, 2,    4,           */ {11, 0},
    /*  24:          3, 4,           */ {3, 8},
    /*  25: 0,       3, 4,           */ {5, 5},
    /*  26:    1,    3, 4,           */ {7, 3},
    /*  27: 0, 1,    3, 4,           */ {9, 1},
    /*  28:       2, 3, 4,           */ {6, 16},
    /*  29: 0,    2, 3, 4,           */ {14, 3},
    /*  30:    1, 2, 3, 4,           */ {12, 12},
    /*  31: 0, 1, 2, 3, 4,           */ {5, 24},
    /*  32:                5,        */ {1, 5},
    /*  33: 0,             5,        */ {3, 1},
    /*  34:    1,          5,        */ {2, 4},
    /*  35: 0, 1,          5,        */ {5, 3},
    /*  36:       2,       5,        */ {3, 6},
    /*  37: 0,    2,       5,        */ {7, 0},
    /*  38:    1, 2,       5,        */ {5, 10},
    /*  39: 0, 1, 2,       5,        */ {9, 0},
    /*  40:          3,    5,        */ {4, 3},
    /*  41: 0,       3,    5,        */ {6, 4},
    /*  42:    1,    3,    5,        */ {6, 11},
    /*  43: 0, 1,    3,    5,        */ {14, 1},
    /*  44:       2, 3,    5,        */ {6, 17},
    /*  45: 0,    2, 3,    5,        */ {12, 4},
    /*  46:    1, 2, 3,    5,        */ {11, 6},
    /*  47: 0, 1, 2, 3,    5,        */ {5, 25},
    /*  48:             4, 5,        */ {2, 8},
    /*  49: 0,          4, 5,        */ {5, 7},
    /*  50:    1,       4, 5,        */ {5, 12},
    /*  51: 0, 1,       4, 5,        */ {8, 1},
    /*  52:       2,    4, 5,        */ {6, 18},
    /*  53: 0,    2,    4, 5,        */ {12, 5},
    /*  54:    1, 2,    4, 5,        */ {14, 7},
    /*  55: 0, 1, 2,    4, 5,        */ {5, 28},
    /*  56:          3, 4, 5,        */ {6, 21},
    /*  57: 0,       3, 4, 5,        */ {11, 4},
    /*  58:    1,    3, 4, 5,        */ {12, 15},
    /*  59: 0, 1,    3, 4, 5,        */ {5, 30},
    /*  60:       2, 3, 4, 5,        */ {10, 5},
    /*  61: 0,    2, 3, 4, 5,        */ {6, 32},
    /*  62:    1, 2, 3, 4, 5,        */ {6, 39},
    /*  63: 0, 1, 2, 3, 4, 5,        */ {2, 12},
    /*  64:                   6,     */ {1, 6},
    /*  65: 0,                6,     */ {4, 0},
    /*  66:    1,             6,     */ {3, 5},
    /*  67: 0, 1,             6,     */ {6, 0},
    /*  68:       2,          6,     */ {2, 6},
    /*  69: 0,    2,          6,     */ {6, 3},
    /*  70:    1, 2,          6,     */ {5, 11},
    /*  71: 0, 1, 2,          6,     */ {14, 0},
    /*  72:          3,       6,     */ {3, 9},
    /*  73: 0,       3,       6,     */ {6, 5},
    /*  74:    1,    3,       6,     */ {7, 4},
    /*  75: 0, 1,    3,       6,     */ {12, 1},
    /*  76:       2, 3,       6,     */ {5, 14},
    /*  77: 0,    2, 3,       6,     */ {11, 3},
    /*  78:    1, 2, 3,       6,     */ {9, 4},
    /*  79: 0, 1, 2, 3,       6,     */ {5, 26},
    /*  80:             4,    6,     */ {3, 10},
    /*  81: 0,          4,    6,     */ {6, 6},
    /*  82:    1,       4,    6,     */ {7, 5},
    /*  83: 0, 1,       4,    6,     */ {12, 2},
    /*  84:       2,    4,    6,     */ {6, 19},
    /*  85: 0,    2,    4,    6,     */ {10, 1},
    /*  86:    1, 2,    4,    6,     */ {12, 13},
    /*  87: 0, 1, 2,    4,    6,     */ {6, 24},
    /*  88:          3, 4,    6,     */ {7, 7},
    /*  89: 0,       3, 4,    6,     */ {12, 9},
    /*  90:    1,    3, 4,    6,     */ {13, 1},
    /*  91: 0, 1,    3, 4,    6,     */ {7, 9},
    /*  92:       2, 3, 4,    6,     */ {12, 20},
    /*  93: 0,    2, 3, 4,    6,     */ {6, 33},
    /*  94:    1, 2, 3, 4,    6,     */ {7, 13},
    /*  95: 0, 1, 2, 3, 4,    6,     */ {3, 12},
    /*  96:                5, 6,     */ {2, 10},
    /*  97: 0,             5, 6,     */ {6, 7},
    /*  98:    1,          5, 6,     */ {5, 13},
    /*  99: 0, 1,          5, 6,     */ {11, 2},
    /* 100:       2,       5, 6,     */ {5, 16},
    /* 101: 0,    2,       5, 6,     */ {12, 7},
    /* 102:    1, 2,       5, 6,     */ {8, 3},
    /* 103: 0, 1, 2,       5, 6,     */ {5, 29},
    /* 104:          3,    5, 6,     */ {6, 22},
    /* 105: 0,       3,    5, 6,     */ {10, 2},
    /* 106:    1,    3,    5, 6,     */ {12, 17},
    /* 107: 0, 1,    3,    5, 6,     */ {6, 27},
    /* 108:       2, 3,    5, 6,     */ {14, 9},
    /* 109: 0,    2, 3,    5, 6,     */ {6, 34},
    /* 110:    1, 2, 3,    5, 6,     */ {5, 39},
    /* 111: 0, 1, 2, 3,    5, 6,     */ {2, 14},
    /* 112:             4, 5, 6,     */ {5, 20},
    /* 113: 0,          4, 5, 6,     */ {14, 5},
    /* 114:    1,       4, 5, 6,     */ {9, 5},
    /* 115: 0, 1,       4, 5, 6,     */ {5, 32},
    /* 116:       2,    4, 5, 6,     */ {11, 10},
    /* 117: 0,    2,    4, 5, 6,     */ {6, 35},
    /* 118:    1, 2,    4, 5, 6,     */ {5, 41},
    /* 119: 0, 1, 2,    4, 5, 6,     */ {2, 16},
    /* 120:          3, 4, 5, 6,     */ {12, 23},
    /* 121: 0,       3, 4, 5, 6,     */ {6, 37},
    /* 122:    1,    3, 4, 5, 6,     */ {7, 14},
    /* 123: 0, 1,    3, 4, 5, 6,     */ {3, 16},
    /* 124:       2, 3, 4, 5, 6,     */ {6, 46},
    /* 125: 0,    2, 3, 4, 5, 6,     */ {4, 6},
    /* 126:    1, 2, 3, 4, 5, 6,     */ {3, 21},
    /* 127: 0, 1, 2, 3, 4, 5, 6,     */ {1, 8},
    /* 128:                      7,  */ {1, 7},
    /* 129: 0,                   7,  */ {3, 2},
    /* 130:    1,                7,  */ {4, 1},
    /* 131: 0, 1,                7,  */ {6, 1},
    /* 132:       2,             7,  */ {3, 7},
    /* 133: 0,    2,             7,  */ {7, 1},
    /* 134:    1, 2,             7,  */ {6, 10},
    /* 135: 0, 1, 2,             7,  */ {12, 0},
    /* 136:          3,          7,  */ {2, 7},
    /* 137: 0,       3,          7,  */ {5, 6},
    /* 138:    1,    3,          7,  */ {6, 12},
    /* 139: 0, 1,    3,          7,  */ {11, 1},
    /* 140:       2, 3,          7,  */ {5, 15},
    /* 141: 0,    2, 3,          7,  */ {9, 2},
    /* 142:    1, 2, 3,          7,  */ {14, 6},
    /* 143: 0, 1, 2, 3,          7,  */ {5, 27},
    /* 144:             4,       7,  */ {2, 9},
    /* 145: 0,          4,       7,  */ {5, 8},
    /* 146:    1,       4,       7,  */ {6, 13},
    /* 147: 0, 1,       4,       7,  */ {14, 2},
    /* 148:       2,    4,       7,  */ {6, 20},
    /* 149: 0,    2,    4,       7,  */ {12, 6},
    /* 150:    1, 2,    4,       7,  */ {10, 3},
    /* 151: 0, 1, 2,    4,       7,  */ {6, 25},
    /* 152:          3, 4,       7,  */ {5, 18},
    /* 153: 0,       3, 4,       7,  */ {8, 2},
    /* 154:    1,    3, 4,       7,  */ {12, 16},
    /* 155: 0, 1,    3, 4,       7,  */ {5, 31},
    /* 156:       2, 3, 4,       7,  */ {11, 9},
    /* 157: 0,    2, 3, 4,       7,  */ {5, 34},
    /* 158:    1, 2, 3, 4,       7,  */ {6, 40},
    /* 159: 0, 1, 2, 3, 4,       7,  */ {2, 13},
    /* 160:                5,    7,  */ {3, 11},
    /* 161: 0,             5,    7,  */ {7, 2},
    /* 162:    1,          5,    7,  */ {6, 14},
    /* 163: 0, 1,          5,    7,  */ {12, 3},
    /* 164:       2,       5,    7,  */ {7, 6},
    /* 165: 0,    2,       5,    7,  */ {13, 0},
    /* 166:    1, 2,       5,    7,  */ {12, 14},
    /* 167: 0, 1, 2,       5,    7,  */ {7, 8},
    /* 168:          3,    5,    7,  */ {6, 23},
    /* 169: 0,       3,    5,    7,  */ {12, 10},
    /* 170:    1,    3,    5,    7,  */ {10, 4},
    /* 171: 0, 1,    3,    5,    7,  */ {6, 28},
    /* 172:       2, 3,    5,    7,  */ {12, 21},
    /* 173: 0,    2, 3,    5,    7,  */ {7, 10},
    /* 174:    1, 2, 3,    5,    7,  */ {6, 41},
    /* 175: 0, 1, 2, 3,    5,    7,  */ {3, 13},
    /* 176:             4, 5,    7,  */ {5, 21},
    /* 177: 0,          4, 5,    7,  */ {9, 3},
    /* 178:    1,       4, 5,    7,  */ {11, 8},
    /* 179: 0, 1,       4, 5,    7,  */ {5, 33},
    /* 180:       2,    4, 5,    7,  */ {12, 22},
    /* 181: 0,    2,    4, 5,    7,  */ {7, 11},
    /* 182:    1, 2,    4, 5,    7,  */ {6, 42},
    /* 183: 0, 1, 2,    4, 5,    7,  */ {3, 14},
    /* 184:          3, 4, 5,    7,  */ {14, 11},
    /* 185: 0,       3, 4, 5,    7,  */ {5, 36},
    /* 186:    1,    3, 4, 5,    7,  */ {6, 44},
    /* 187: 0, 1,    3, 4, 5,    7,  */ {2, 17},
    /* 188:       2, 3, 4, 5,    7,  */ {6, 47},
    /* 189: 0,    2, 3, 4, 5,    7,  */ {3, 18},
    /* 190:    1, 2, 3, 4, 5,    7,  */ {4, 7},
    /* 191: 0, 1, 2, 3, 4, 5,    7,  */ {1, 9},
    /* 192:                   6, 7,  */ {2, 11},
    /* 193: 0,                6, 7,  */ {6, 8},
    /* 194:    1,             6, 7,  */ {6, 15},
    /* 195: 0, 1,             6, 7,  */ {10, 0},
    /* 196:       2,          6, 7,  */ {5, 17},
    /* 197: 0,    2,          6, 7,  */ {12, 8},
    /* 198:    1, 2,          6, 7,  */ {11, 7},
    /* 199: 0, 1, 2,          6, 7,  */ {6, 26},
    /* 200:          3,       6, 7,  */ {5, 19},
    /* 201: 0,       3,       6, 7,  */ {14, 4},
    /* 202:    1,    3,       6, 7,  */ {12, 18},
    /* 203: 0, 1,    3,       6, 7,  */ {6, 29},
    /* 204:       2, 3,       6, 7,  */ {8, 4},
    /* 205: 0,    2, 3,       6, 7,  */ {5, 35},
    /* 206:    1, 2, 3,       6, 7,  */ {5, 40},
    /* 207: 0, 1, 2, 3,       6, 7,  */ {2, 15},
    /* 208:             4,    6, 7,  */ {5, 22},
    /* 209: 0,          4,    6, 7,  */ {11, 5},
    /* 210:    1,       4,    6, 7,  */ {12, 19},
    /* 211: 0, 1,       4,    6, 7,  */ {6, 30},
    /* 212:       2,    4,    6, 7,  */ {14, 10},
    /* 213: 0,    2,    4,    6, 7,  */ {6, 36},
    /* 214:    1, 2,    4,    6, 7,  */ {6, 43},
    /* 215: 0, 1, 2,    4,    6, 7,  */ {4, 4},
    /* 216:          3, 4,    6, 7,  */ {9, 7},
    /* 217: 0,       3, 4,    6, 7,  */ {5, 37},
    /* 218:    1,    3, 4,    6, 7,  */ {7, 15},
    /* 219: 0, 1,    3, 4,    6, 7,  */ {3, 17},
    /* 220:       2, 3, 4,    6, 7,  */ {5, 44},
    /* 221: 0,    2, 3, 4,    6, 7,  */ {2, 19},
    /* 222:    1, 2, 3, 4,    6, 7,  */ {3, 22},
    /* 223: 0, 1, 2, 3, 4,    6, 7,  */ {1, 10},
    /* 224:                5, 6, 7,  */ {5, 23},
    /* 225: 0,             5, 6, 7,  */ {12, 11},
    /* 226:    1,          5, 6, 7,  */ {14, 8},
    /* 227: 0, 1,          5, 6, 7,  */ {6, 31},
    /* 228:       2,       5, 6, 7,  */ {9, 6},
    /* 229: 0,    2,       5, 6, 7,  */ {7, 12},
    /* 230:    1, 2,       5, 6, 7,  */ {5, 42},
    /* 231: 0, 1, 2,       5, 6, 7,  */ {3, 15},
    /* 232:          3,    5, 6, 7,  */ {11, 11},
    /* 233: 0,       3,    5, 6, 7,  */ {6, 38},
    /* 234:    1,    3,    5, 6, 7,  */ {6, 45},
    /* 235: 0, 1,    3,    5, 6, 7,  */ {4, 5},
    /* 236:       2, 3,    5, 6, 7,  */ {5, 45},
    /* 237: 0,    2, 3,    5, 6, 7,  */ {3, 19},
    /* 238:    1, 2, 3,    5, 6, 7,  */ {2, 21},
    /* 239: 0, 1, 2, 3,    5, 6, 7,  */ {1, 11},
    /* 240:             4, 5, 6, 7,  */ {8, 5},
    /* 241: 0,          4, 5, 6, 7,  */ {5, 38},
    /* 242:    1,       4, 5, 6, 7,  */ {5, 43},
    /* 243: 0, 1,       4, 5, 6, 7,  */ {2, 18},
    /* 244:       2,    4, 5, 6, 7,  */ {5, 46},
    /* 245: 0,    2,    4, 5, 6, 7,  */ {3, 20},
    /* 246:    1, 2,    4, 5, 6, 7,  */ {2, 22},
    /* 247: 0, 1, 2,    4, 5, 6, 7,  */ {1, 12},
    /* 248:          3, 4, 5, 6, 7,  */ {5, 47},
    /* 249: 0,       3, 4, 5, 6, 7,  */ {2, 20},
    /* 250:    1,    3, 4, 5, 6, 7,  */ {3, 23},
    /* 251: 0, 1,    3, 4, 5, 6, 7,  */ {1, 13},
    /* 252:       2, 3, 4, 5, 6, 7,  */ {2, 23},
    /* 253: 0,    2, 3, 4, 5, 6, 7,  */ {1, 14},
    /* 254:    1, 2, 3, 4, 5, 6, 7,  */ {1, 15},
    /* 255: 0, 1, 2, 3, 4, 5, 6, 7,  */ {0, -1}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief tiling table for case 1
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling1[16][3] = {
    /*   1: 0,                       */ {0, 8, 3},
    /*   2:    1,                    */ {0, 1, 9},
    /*   4:       2,                 */ {1, 2, 10},
    /*   8:          3,              */ {3, 11, 2},
    /*  16:             4,           */ {4, 7, 8},
    /*  32:                5,        */ {9, 5, 4},
    /*  64:                   6,     */ {10, 6, 5},
    /* 128:                      7,  */ {7, 6, 11},
    /* 127: 0, 1, 2, 3, 4, 5, 6,     */ {7, 11, 6},
    /* 191: 0, 1, 2, 3, 4, 5,    7,  */ {10, 5, 6},
    /* 223: 0, 1, 2, 3, 4,    6, 7,  */ {9, 4, 5},
    /* 239: 0, 1, 2, 3,    5, 6, 7,  */ {4, 8, 7},
    /* 247: 0, 1, 2,    4, 5, 6, 7,  */ {3, 2, 11},
    /* 251: 0, 1,    3, 4, 5, 6, 7,  */ {1, 10, 2},
    /* 253: 0,    2, 3, 4, 5, 6, 7,  */ {0, 9, 1},
    /* 254:    1, 2, 3, 4, 5, 6, 7,  */ {0, 3, 8}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief tiling table for case 2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling2[24][6] = {
    /*   3: 0, 1,                    */ {1, 8, 3, 9, 8, 1},
    /*   9: 0,       3,              */ {0, 11, 2, 8, 11, 0},
    /*  17: 0,          4,           */ {4, 3, 0, 7, 3, 4},
    /*   6:    1, 2,                 */ {9, 2, 10, 0, 2, 9},
    /*  34:    1,          5,        */ {0, 5, 4, 1, 5, 0},
    /*  12:       2, 3,              */ {3, 10, 1, 11, 10, 3},
    /*  68:       2,          6,     */ {1, 6, 5, 2, 6, 1},
    /* 136:          3,          7,  */ {7, 2, 3, 6, 2, 7},
    /*  48:             4, 5,        */ {9, 7, 8, 5, 7, 9},
    /* 144:             4,       7,  */ {6, 8, 4, 11, 8, 6},
    /*  96:                5, 6,     */ {10, 4, 9, 6, 4, 10},
    /* 192:                   6, 7,  */ {11, 5, 10, 7, 5, 11},
    /*  63: 0, 1, 2, 3, 4, 5,        */ {11, 10, 5, 7, 11, 5},
    /* 159: 0, 1, 2, 3, 4,       7,  */ {10, 9, 4, 6, 10, 4},
    /* 111: 0, 1, 2, 3,    5, 6,     */ {6, 4, 8, 11, 6, 8},
    /* 207: 0, 1, 2, 3,       6, 7,  */ {9, 8, 7, 5, 9, 7},
    /* 119: 0, 1, 2,    4, 5, 6,     */ {7, 3, 2, 6, 7, 2},
    /* 187: 0, 1,    3, 4, 5,    7,  */ {1, 5, 6, 2, 1, 6},
    /* 243: 0, 1,       4, 5, 6, 7,  */ {3, 1, 10, 11, 3, 10},
    /* 221: 0,    2, 3, 4,    6, 7,  */ {0, 4, 5, 1, 0, 5},
    /* 249: 0,       3, 4, 5, 6, 7,  */ {9, 10, 2, 0, 9, 2},
    /* 238:    1, 2, 3,    5, 6, 7,  */ {4, 0, 3, 7, 4, 3},
    /* 246:    1, 2,    4, 5, 6, 7,  */ {0, 2, 11, 8, 0, 11},
    /* 252:       2, 3, 4, 5, 6, 7,  */ {1, 3, 8, 9, 1, 8}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief test table for case 3
 * One face to test
 * When the test on the specified face is positive : 4 first triangles
 * When the test on the specified face is negative : 2 last triangles
 *
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::test3[24] = {
    /*   5: 0,    2,                 */ 5,
    /*  33: 0,             5,        */ 1,
    /* 129: 0,                   7,  */ 4,
    /*  10:    1,    3,              */ 5,
    /*  18:    1,       4,           */ 1,
    /*  66:    1,             6,     */ 2,
    /*  36:       2,       5,        */ 2,
    /* 132:       2,             7,  */ 3,
    /*  24:          3, 4,           */ 4,
    /*  72:          3,       6,     */ 3,
    /*  80:             4,    6,     */ 6,
    /* 160:                5,    7,  */ 6,
    /*  95: 0, 1, 2, 3, 4,    6,     */ -6,
    /* 175: 0, 1, 2, 3,    5,    7,  */ -6,
    /* 183: 0, 1, 2,    4, 5,    7,  */ -3,
    /* 231: 0, 1, 2,       5, 6, 7,  */ -4,
    /* 123: 0, 1,    3, 4, 5, 6,     */ -3,
    /* 219: 0, 1,    3, 4,    6, 7,  */ -2,
    /* 189: 0,    2, 3, 4, 5,    7,  */ -2,
    /* 237: 0,    2, 3,    5, 6, 7,  */ -1,
    /* 245: 0,    2,    4, 5, 6, 7,  */ -5,
    /* 126:    1, 2, 3, 4, 5, 6,     */ -4,
    /* 222:    1, 2, 3, 4,    6, 7,  */ -1,
    /* 250:    1,    3, 4, 5, 6, 7,  */ -5};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 3.1
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling3_1[24][6] = {
    /*   5: 0,    2,                 */ {0, 8, 3, 1, 2, 10},
    /*  33: 0,             5,        */ {9, 5, 4, 0, 8, 3},
    /* 129: 0,                   7,  */ {3, 0, 8, 11, 7, 6},
    /*  10:    1,    3,              */ {1, 9, 0, 2, 3, 11},
    /*  18:    1,       4,           */ {0, 1, 9, 8, 4, 7},
    /*  66:    1,             6,     */ {9, 0, 1, 5, 10, 6},
    /*  36:       2,       5,        */ {1, 2, 10, 9, 5, 4},
    /* 132:       2,             7,  */ {10, 1, 2, 6, 11, 7},
    /*  24:          3, 4,           */ {8, 4, 7, 3, 11, 2},
    /*  72:          3,       6,     */ {2, 3, 11, 10, 6, 5},
    /*  80:             4,    6,     */ {5, 10, 6, 4, 7, 8},
    /* 160:                5,    7,  */ {4, 9, 5, 7, 6, 11},
    /*  95: 0, 1, 2, 3, 4,    6,     */ {5, 9, 4, 11, 6, 7},
    /* 175: 0, 1, 2, 3,    5,    7,  */ {6, 10, 5, 8, 7, 4},
    /* 183: 0, 1, 2,    4, 5,    7,  */ {11, 3, 2, 5, 6, 10},
    /* 231: 0, 1, 2,       5, 6, 7,  */ {7, 4, 8, 2, 11, 3},
    /* 123: 0, 1,    3, 4, 5, 6,     */ {2, 1, 10, 7, 11, 6},
    /* 219: 0, 1,    3, 4,    6, 7,  */ {10, 2, 1, 4, 5, 9},
    /* 189: 0,    2, 3, 4, 5,    7,  */ {1, 0, 9, 6, 10, 5},
    /* 237: 0,    2, 3,    5, 6, 7,  */ {9, 1, 0, 7, 4, 8},
    /* 245: 0,    2,    4, 5, 6, 7,  */ {0, 9, 1, 11, 3, 2},
    /* 126:    1, 2, 3, 4, 5, 6,     */ {8, 0, 3, 6, 7, 11},
    /* 222:    1, 2, 3, 4,    6, 7,  */ {4, 5, 9, 3, 8, 0},
    /* 250:    1,    3, 4, 5, 6, 7,  */ {3, 8, 0, 10, 2, 1}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 3.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling3_2[24][12] = {
    /*   5: 0,    2,                 */ {10, 3, 2, 10, 8, 3, 10, 1, 0, 8, 10, 0},
    /*  33: 0,             5,        */ {3, 4, 8, 3, 5, 4, 3, 0, 9, 5, 3, 9},
    /* 129: 0,                   7,  */ {6, 8, 7, 6, 0, 8, 6, 11, 3, 0, 6, 3},
    /*  10:    1,    3,              */ {11, 0, 3, 11, 9, 0, 11, 2, 1, 9, 11, 1},
    /*  18:    1,       4,           */ {7, 9, 4, 7, 1, 9, 7, 8, 0, 1, 7, 0},
    /*  66:    1,             6,     */ {6, 1, 10, 6, 0, 1, 9, 0, 6, 9, 6, 5},
    /*  36:       2,       5,        */ {4, 10, 5, 4, 2, 10, 4, 9, 1, 2, 4, 1},
    /* 132:       2,             7,  */ {7, 2, 11, 7, 1, 2, 7, 6, 10, 1, 7, 10},
    /*  24:          3, 4,           */ {2, 7, 11, 2, 4, 7, 2, 3, 8, 4, 2, 8},
    /*  72:          3,       6,     */ {5, 11, 6, 5, 3, 11, 5, 10, 2, 3, 5, 2},
    /*  80:             4,    6,     */ {8, 6, 7, 8, 10, 6, 8, 4, 5, 10, 8, 5},
    /* 160:                5,    7,  */ {11, 5, 6, 11, 9, 5, 11, 7, 4, 9, 11, 4},
    /*  95: 0, 1, 2, 3, 4,    6,     */ {6, 5, 11, 5, 9, 11, 4, 7, 11, 4, 11, 9},
    /* 175: 0, 1, 2, 3,    5,    7,  */ {7, 6, 8, 6, 10, 8, 5, 4, 8, 5, 8, 10},
    /* 183: 0, 1, 2,    4, 5,    7,  */ {6, 11, 5, 11, 3, 5, 2, 10, 5, 2, 5, 3},
    /* 231: 0, 1, 2,       5, 6, 7,  */ {11, 7, 2, 7, 4, 2, 8, 3, 2, 8, 2, 4},
    /* 123: 0, 1,    3, 4, 5, 6,     */ {11, 2, 7, 2, 1, 7, 10, 6, 7, 10, 7, 1},
    /* 219: 0, 1,    3, 4,    6, 7,  */ {5, 10, 4, 10, 2, 4, 1, 9, 4, 1, 4, 2},
    /* 189: 0,    2, 3, 4, 5,    7,  */ {10, 1, 6, 1, 0, 6, 6, 0, 9, 5, 6, 9},
    /* 237: 0,    2, 3,    5, 6, 7,  */ {4, 9, 7, 9, 1, 7, 0, 8, 7, 0, 7, 1},
    /* 245: 0,    2,    4, 5, 6, 7,  */ {3, 0, 11, 0, 9, 11, 1, 2, 11, 1, 11, 9},
    /* 126:    1, 2, 3, 4, 5, 6,     */ {7, 8, 6, 8, 0, 6, 3, 11, 6, 3, 6, 0},
    /* 222:    1, 2, 3, 4,    6, 7,  */ {8, 4, 3, 4, 5, 3, 9, 0, 3, 9, 3, 5},
    /* 250:    1,    3, 4, 5, 6, 7,  */ {2, 3, 10, 3, 8, 10, 0, 1, 10, 0, 10, 8}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief test table for case 4
 * Interior to test
 * When the test on the interior is negative : 2 first triangles
 * When the test on the interior is positive : 6 last triangles
 *
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::test4[8] = {
    /*  65: 0,                6,     */ 7,
    /* 130:    1,                7,  */ 7,
    /*  20:       2,    4,           */ 7,
    /*  40:          3,    5,        */ 7,
    /* 215: 0, 1, 2,    4,    6, 7,  */ -7,
    /* 235: 0, 1,    3,    5, 6, 7,  */ -7,
    /* 125: 0,    2, 3, 4, 5, 6,     */ -7,
    /* 190:    1, 2, 3, 4, 5,    7,  */ -7};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 4.1
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling4_1[8][6] = {
    /*  65: 0,                6,     */ {0, 8, 3, 5, 10, 6},
    /* 130:    1,                7,  */ {0, 1, 9, 11, 7, 6},
    /*  20:       2,    4,           */ {1, 2, 10, 8, 4, 7},
    /*  40:          3,    5,        */ {9, 5, 4, 2, 3, 11},
    /* 215: 0, 1, 2,    4,    6, 7,  */ {4, 5, 9, 11, 3, 2},
    /* 235: 0, 1,    3,    5, 6, 7,  */ {10, 2, 1, 7, 4, 8},
    /* 125: 0,    2, 3, 4, 5, 6,     */ {9, 1, 0, 6, 7, 11},
    /* 190:    1, 2, 3, 4, 5,    7,  */ {3, 8, 0, 6, 10, 5}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 4.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling4_2[8][18] = {
    /*  65: 0,                6,     */ {8, 5, 0, 5, 8, 6, 3, 6, 8, 6, 3, 10, 0, 10, 3, 10, 0, 5},
    /* 130:    1,                7,  */ {9, 6, 1, 6, 9, 7, 0, 7, 9, 7, 0, 11, 1, 11, 0, 11, 1, 6},
    /*  20:       2,    4,           */ {10, 7, 2, 7, 10, 4, 1, 4, 10, 4, 1, 8, 2, 8, 1, 8, 2, 7},
    /*  40:          3,    5,        */ {11, 4, 3, 4, 11, 5, 2, 5, 11, 5, 2, 9, 3, 9, 2, 9, 3, 4},
    /* 215: 0, 1, 2,    4,    6, 7,  */ {3, 4, 11, 5, 11, 4, 11, 5, 2, 9, 2, 5, 2, 9, 3, 4, 3, 9},
    /* 235: 0, 1,    3,    5, 6, 7,  */ {2, 7, 10, 4, 10, 7, 10, 4, 1, 8, 1, 4, 1, 8, 2, 7, 2, 8},
    /* 125: 0,    2, 3, 4, 5, 6,     */ {1, 6, 9, 7, 9, 6, 9, 7, 0, 11, 0, 7, 0, 11, 1, 6, 1, 11},
    /* 190:    1, 2, 3, 4, 5,    7,  */ {0, 5, 8, 6, 8, 5, 8, 6, 3, 10, 3, 6, 3, 10, 0, 5, 0, 10}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief tiling table for case 5
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling5[48][9] = {
    /*   7: 0, 1, 2,                 */ {2, 8, 3, 2, 10, 8, 10, 9, 8},
    /*  11: 0, 1,    3,              */ {1, 11, 2, 1, 9, 11, 9, 8, 11},
    /*  19: 0, 1,       4,           */ {4, 1, 9, 4, 7, 1, 7, 3, 1},
    /*  35: 0, 1,          5,        */ {8, 5, 4, 8, 3, 5, 3, 1, 5},
    /*  13: 0,    2, 3,              */ {0, 10, 1, 0, 8, 10, 8, 11, 10},
    /*  25: 0,       3, 4,           */ {11, 4, 7, 11, 2, 4, 2, 0, 4},
    /* 137: 0,       3,          7,  */ {7, 0, 8, 7, 6, 0, 6, 2, 0},
    /*  49: 0,          4, 5,        */ {9, 3, 0, 9, 5, 3, 5, 7, 3},
    /* 145: 0,          4,       7,  */ {3, 6, 11, 3, 0, 6, 0, 4, 6},
    /*  14:    1, 2, 3,              */ {3, 9, 0, 3, 11, 9, 11, 10, 9},
    /*  38:    1, 2,       5,        */ {5, 2, 10, 5, 4, 2, 4, 0, 2},
    /*  70:    1, 2,          6,     */ {9, 6, 5, 9, 0, 6, 0, 2, 6},
    /*  50:    1,       4, 5,        */ {0, 7, 8, 0, 1, 7, 1, 5, 7},
    /*  98:    1,          5, 6,     */ {10, 0, 1, 10, 6, 0, 6, 4, 0},
    /*  76:       2, 3,       6,     */ {6, 3, 11, 6, 5, 3, 5, 1, 3},
    /* 140:       2, 3,          7,  */ {10, 7, 6, 10, 1, 7, 1, 3, 7},
    /* 100:       2,       5, 6,     */ {1, 4, 9, 1, 2, 4, 2, 6, 4},
    /* 196:       2,          6, 7,  */ {11, 1, 2, 11, 7, 1, 7, 5, 1},
    /* 152:          3, 4,       7,  */ {8, 2, 3, 8, 4, 2, 4, 6, 2},
    /* 200:          3,       6, 7,  */ {2, 5, 10, 2, 3, 5, 3, 7, 5},
    /* 112:             4, 5, 6,     */ {7, 10, 6, 7, 8, 10, 8, 9, 10},
    /* 176:             4, 5,    7,  */ {6, 9, 5, 6, 11, 9, 11, 8, 9},
    /* 208:             4,    6, 7,  */ {5, 8, 4, 5, 10, 8, 10, 11, 8},
    /* 224:                5, 6, 7,  */ {4, 11, 7, 4, 9, 11, 9, 10, 11},
    /*  31: 0, 1, 2, 3, 4,           */ {4, 7, 11, 4, 11, 9, 9, 11, 10},
    /*  47: 0, 1, 2, 3,    5,        */ {5, 4, 8, 5, 8, 10, 10, 8, 11},
    /*  79: 0, 1, 2, 3,       6,     */ {6, 5, 9, 6, 9, 11, 11, 9, 8},
    /* 143: 0, 1, 2, 3,          7,  */ {7, 6, 10, 7, 10, 8, 8, 10, 9},
    /*  55: 0, 1, 2,    4, 5,        */ {2, 10, 5, 2, 5, 3, 3, 5, 7},
    /* 103: 0, 1, 2,       5, 6,     */ {8, 3, 2, 8, 2, 4, 4, 2, 6},
    /*  59: 0, 1,    3, 4, 5,        */ {11, 2, 1, 11, 1, 7, 7, 1, 5},
    /* 155: 0, 1,    3, 4,       7,  */ {1, 9, 4, 1, 4, 2, 2, 4, 6},
    /* 115: 0, 1,       4, 5, 6,     */ {10, 6, 7, 10, 7, 1, 1, 7, 3},
    /* 179: 0, 1,       4, 5,    7,  */ {6, 11, 3, 6, 3, 5, 5, 3, 1},
    /* 157: 0,    2, 3, 4,       7,  */ {10, 1, 0, 10, 0, 6, 6, 0, 4},
    /* 205: 0,    2, 3,       6, 7,  */ {0, 8, 7, 0, 7, 1, 1, 7, 5},
    /* 185: 0,       3, 4, 5,    7,  */ {9, 5, 6, 9, 6, 0, 0, 6, 2},
    /* 217: 0,       3, 4,    6, 7,  */ {5, 10, 2, 5, 2, 4, 4, 2, 0},
    /* 241: 0,          4, 5, 6, 7,  */ {3, 0, 9, 3, 9, 11, 11, 9, 10},
    /* 110:    1, 2, 3,    5, 6,     */ {3, 11, 6, 3, 6, 0, 0, 6, 4},
    /* 206:    1, 2, 3,       6, 7,  */ {9, 0, 3, 9, 3, 5, 5, 3, 7},
    /* 118:    1, 2,    4, 5, 6,     */ {7, 8, 0, 7, 0, 6, 6, 0, 2},
    /* 230:    1, 2,       5, 6, 7,  */ {11, 7, 4, 11, 4, 2, 2, 4, 0},
    /* 242:    1,       4, 5, 6, 7,  */ {0, 1, 10, 0, 10, 8, 8, 10, 11},
    /* 220:       2, 3, 4,    6, 7,  */ {8, 4, 5, 8, 5, 3, 3, 5, 1},
    /* 236:       2, 3,    5, 6, 7,  */ {4, 9, 1, 4, 1, 7, 7, 1, 3},
    /* 244:       2,    4, 5, 6, 7,  */ {1, 2, 11, 1, 11, 9, 9, 11, 8},
    /* 248:          3, 4, 5, 6, 7,  */ {2, 3, 8, 2, 8, 10, 10, 8, 9}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief test table for case 6
 * 1 face to test + eventually the interior
 * When the test on the specified face is positive : 5 first triangles
 * When the test on the specified face is negative :
 * - if the test on the interior is negative : 3 middle triangles
 * - if the test on the interior is positive : 8 last triangles
 * The support edge for the interior test is marked as the 3rd column.
 *
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::test6[48][3] = {
    /*  67: 0, 1,             6,     */ {2, 7, 10},
    /* 131: 0, 1,                7,  */ {4, 7, 11},
    /*  21: 0,    2,    4,           */ {5, 7, 1},
    /*  69: 0,    2,          6,     */ {5, 7, 3},
    /*  41: 0,       3,    5,        */ {1, 7, 9},
    /*  73: 0,       3,       6,     */ {3, 7, 10},
    /*  81: 0,          4,    6,     */ {6, 7, 5},
    /*  97: 0,             5, 6,     */ {1, 7, 8},
    /* 193: 0,                6, 7,  */ {4, 7, 8},
    /*  22:    1, 2,    4,           */ {1, 7, 8},
    /* 134:    1, 2,             7,  */ {3, 7, 11},
    /*  42:    1,    3,    5,        */ {5, 7, 2},
    /* 138:    1,    3,          7,  */ {5, 7, 0},
    /* 146:    1,       4,       7,  */ {1, 7, 9},
    /* 162:    1,          5,    7,  */ {6, 7, 6},
    /* 194:    1,             6, 7,  */ {2, 7, 9},
    /*  28:       2, 3, 4,           */ {4, 7, 8},
    /*  44:       2, 3,    5,        */ {2, 7, 9},
    /*  52:       2,    4, 5,        */ {2, 7, 10},
    /*  84:       2,    4,    6,     */ {6, 7, 7},
    /* 148:       2,    4,       7,  */ {3, 7, 10},
    /*  56:          3, 4, 5,        */ {4, 7, 11},
    /* 104:          3,    5, 6,     */ {3, 7, 11},
    /* 168:          3,    5,    7,  */ {6, 7, 4},
    /*  87: 0, 1, 2,    4,    6,     */ {-6, -7, 4},
    /* 151: 0, 1, 2,    4,       7,  */ {-3, -7, 11},
    /* 199: 0, 1, 2,          6, 7,  */ {-4, -7, 11},
    /* 107: 0, 1,    3,    5, 6,     */ {-3, -7, 10},
    /* 171: 0, 1,    3,    5,    7,  */ {-6, -7, 7},
    /* 203: 0, 1,    3,       6, 7,  */ {-2, -7, 10},
    /* 211: 0, 1,       4,    6, 7,  */ {-2, -7, 9},
    /* 227: 0, 1,          5, 6, 7,  */ {-4, -7, 8},
    /*  61: 0,    2, 3, 4, 5,        */ {-2, -7, 9},
    /*  93: 0,    2, 3, 4,    6,     */ {-6, -7, 6},
    /* 109: 0,    2, 3,    5, 6,     */ {-1, -7, 9},
    /* 117: 0,    2,    4, 5, 6,     */ {-5, -7, 0},
    /* 213: 0,    2,    4,    6, 7,  */ {-5, -7, 2},
    /* 121: 0,       3, 4, 5, 6,     */ {-3, -7, 11},
    /* 233: 0,       3,    5, 6, 7,  */ {-1, -7, 8},
    /*  62:    1, 2, 3, 4, 5,        */ {-4, -7, 8},
    /* 158:    1, 2, 3, 4,       7,  */ {-1, -7, 8},
    /* 174:    1, 2, 3,    5,    7,  */ {-6, -7, 5},
    /* 182:    1, 2,    4, 5,    7,  */ {-3, -7, 10},
    /* 214:    1, 2,    4,    6, 7,  */ {-1, -7, 9},
    /* 186:    1,    3, 4, 5,    7,  */ {-5, -7, 3},
    /* 234:    1,    3,    5, 6, 7,  */ {-5, -7, 1},
    /* 124:       2, 3, 4, 5, 6,     */ {-4, -7, 11},
    /* 188:       2, 3, 4, 5,    7,  */ {-2, -7, 10}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 6.1.1
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling6_1_1[48][9] = {
    /*  67: 0, 1,             6,     */ {6, 5, 10, 3, 1, 8, 9, 8, 1},
    /* 131: 0, 1,                7,  */ {11, 7, 6, 9, 3, 1, 3, 9, 8},
    /*  21: 0,    2,    4,           */ {1, 2, 10, 7, 0, 4, 0, 7, 3},
    /*  69: 0,    2,          6,     */ {3, 0, 8, 5, 2, 6, 2, 5, 1},
    /*  41: 0,       3,    5,        */ {5, 4, 9, 2, 0, 11, 8, 11, 0},
    /*  73: 0,       3,       6,     */ {10, 6, 5, 8, 2, 0, 2, 8, 11},
    /*  81: 0,          4,    6,     */ {10, 6, 5, 0, 4, 3, 7, 3, 4},
    /*  97: 0,             5, 6,     */ {3, 0, 8, 6, 4, 10, 9, 10, 4},
    /* 193: 0,                6, 7,  */ {8, 3, 0, 10, 7, 5, 7, 10, 11},
    /*  22:    1, 2,    4,           */ {8, 4, 7, 10, 0, 2, 0, 10, 9},
    /* 134:    1, 2,             7,  */ {7, 6, 11, 0, 2, 9, 10, 9, 2},
    /*  42:    1,    3,    5,        */ {2, 3, 11, 4, 1, 5, 1, 4, 0},
    /* 138:    1,    3,          7,  */ {0, 1, 9, 6, 3, 7, 3, 6, 2},
    /* 146:    1,       4,       7,  */ {9, 0, 1, 11, 4, 6, 4, 11, 8},
    /* 162:    1,          5,    7,  */ {11, 7, 6, 1, 5, 0, 4, 0, 5},
    /* 194:    1,             6, 7,  */ {0, 1, 9, 7, 5, 11, 10, 11, 5},
    /*  28:       2, 3, 4,           */ {4, 7, 8, 1, 3, 10, 11, 10, 3},
    /*  44:       2, 3,    5,        */ {9, 5, 4, 11, 1, 3, 1, 11, 10},
    /*  52:       2,    4, 5,        */ {10, 1, 2, 8, 5, 7, 5, 8, 9},
    /*  84:       2,    4,    6,     */ {8, 4, 7, 2, 6, 1, 5, 1, 6},
    /* 148:       2,    4,       7,  */ {1, 2, 10, 4, 6, 8, 11, 8, 6},
    /*  56:          3, 4, 5,        */ {2, 3, 11, 5, 7, 9, 8, 9, 7},
    /* 104:          3,    5, 6,     */ {11, 2, 3, 9, 6, 4, 6, 9, 10},
    /* 168:          3,    5,    7,  */ {9, 5, 4, 3, 7, 2, 6, 2, 7},
    /*  87: 0, 1, 2,    4,    6,     */ {4, 5, 9, 2, 7, 3, 7, 2, 6},
    /* 151: 0, 1, 2,    4,       7,  */ {3, 2, 11, 4, 6, 9, 10, 9, 6},
    /* 199: 0, 1, 2,          6, 7,  */ {11, 3, 2, 9, 7, 5, 7, 9, 8},
    /* 107: 0, 1,    3,    5, 6,     */ {10, 2, 1, 8, 6, 4, 6, 8, 11},
    /* 171: 0, 1,    3,    5,    7,  */ {7, 4, 8, 1, 6, 2, 6, 1, 5},
    /* 203: 0, 1,    3,       6, 7,  */ {2, 1, 10, 7, 5, 8, 9, 8, 5},
    /* 211: 0, 1,       4,    6, 7,  */ {4, 5, 9, 3, 1, 11, 10, 11, 1},
    /* 227: 0, 1,          5, 6, 7,  */ {8, 7, 4, 10, 3, 1, 3, 10, 11},
    /*  61: 0,    2, 3, 4, 5,        */ {9, 1, 0, 11, 5, 7, 5, 11, 10},
    /*  93: 0,    2, 3, 4,    6,     */ {6, 7, 11, 0, 5, 1, 5, 0, 4},
    /* 109: 0,    2, 3,    5, 6,     */ {1, 0, 9, 6, 4, 11, 8, 11, 4},
    /* 117: 0,    2,    4, 5, 6,     */ {9, 1, 0, 7, 3, 6, 2, 6, 3},
    /* 213: 0,    2,    4,    6, 7,  */ {11, 3, 2, 5, 1, 4, 0, 4, 1},
    /* 121: 0,       3, 4, 5, 6,     */ {11, 6, 7, 9, 2, 0, 2, 9, 10},
    /* 233: 0,       3,    5, 6, 7,  */ {7, 4, 8, 2, 0, 10, 9, 10, 0},
    /*  62:    1, 2, 3, 4, 5,        */ {0, 3, 8, 5, 7, 10, 11, 10, 7},
    /* 158:    1, 2, 3, 4,       7,  */ {8, 0, 3, 10, 4, 6, 4, 10, 9},
    /* 174:    1, 2, 3,    5,    7,  */ {5, 6, 10, 3, 4, 0, 4, 3, 7},
    /* 182:    1, 2,    4, 5,    7,  */ {5, 6, 10, 0, 2, 8, 11, 8, 2},
    /* 214:    1, 2,    4,    6, 7,  */ {9, 4, 5, 11, 0, 2, 0, 11, 8},
    /* 186:    1,    3, 4, 5,    7,  */ {8, 0, 3, 6, 2, 5, 1, 5, 2},
    /* 234:    1,    3,    5, 6, 7,  */ {10, 2, 1, 4, 0, 7, 3, 7, 0},
    /* 124:       2, 3, 4, 5, 6,     */ {6, 7, 11, 1, 3, 9, 8, 9, 3},
    /* 188:       2, 3, 4, 5,    7,  */ {10, 5, 6, 8, 1, 3, 1, 8, 9}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 6.1.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling6_1_2[48][27] = {
    /*  67: 0, 1,             6,     */ {1, 12, 3, 12, 10, 3, 6, 3, 10, 3,  6,  8, 5, 8,
                                         6, 8,  5, 12, 12, 9, 8, 1, 9,  12, 12, 5, 10},
    /* 131: 0, 1,                7,  */ {1, 12, 3, 1, 11, 12, 11, 1,  6, 9, 6,  1, 6, 9,
                                         7, 12, 7, 9, 9,  8,  12, 12, 8, 3, 11, 7, 12},
    /*  21: 0,    2,    4,           */ {4, 12, 0, 4, 1, 12, 1,  4,  10, 7, 10, 4, 10, 7,
                                         2, 12, 2, 7, 7, 3,  12, 12, 3,  0, 1,  2, 12},
    /*  69: 0,    2,          6,     */ {6, 12, 2, 6, 3, 12, 3,  6,  8, 5, 8, 6, 8, 5,
                                         0, 12, 0, 5, 5, 1,  12, 12, 1, 2, 3, 0, 12},
    /*  41: 0,       3,    5,        */ {0, 12, 2, 12, 9,  2, 5,  2, 9, 2,  5,  11, 4, 11,
                                         5, 11, 4, 12, 12, 8, 11, 0, 8, 12, 12, 4,  9},
    /*  73: 0,       3,       6,     */ {0, 12, 2, 0, 10, 12, 10, 0,  5,  8, 5,  0, 5, 8,
                                         6, 12, 6, 8, 8,  11, 12, 12, 11, 2, 10, 6, 12},
    /*  81: 0,          4,    6,     */ {4,  12, 0, 12, 5,  0, 10, 0, 5, 0,  10, 3, 6, 3,
                                         10, 3,  6, 12, 12, 7, 3,  4, 7, 12, 12, 6, 5},
    /*  97: 0,             5, 6,     */ {4, 12, 6, 12, 8,  6, 3,  6, 8, 6,  3,  10, 0, 10,
                                         3, 10, 0, 12, 12, 9, 10, 4, 9, 12, 12, 0,  8},
    /* 193: 0,                6, 7,  */ {5, 12, 7, 5,  8,  12, 8,  5,  0,  10, 0, 5, 0, 10,
                                         3, 12, 3, 10, 10, 11, 12, 12, 11, 7,  8, 3, 12},
    /*  22:    1, 2,    4,           */ {2, 12, 0, 2,  8,  12, 8,  2,  7, 10, 7, 2, 7, 10,
                                         4, 12, 4, 10, 10, 9,  12, 12, 9, 0,  8, 4, 12},
    /* 134:    1, 2,             7,  */ {2, 12, 0, 12, 11, 0,  7, 0, 11, 0,  7,  9, 6, 9,
                                         7, 9,  6, 12, 12, 10, 9, 2, 10, 12, 12, 6, 11},
    /*  42:    1,    3,    5,        */ {5, 12, 1, 5, 2, 12, 2,  5,  11, 4, 11, 5, 11, 4,
                                         3, 12, 3, 4, 4, 0,  12, 12, 0,  1, 2,  3, 12},
    /* 138:    1,    3,          7,  */ {7, 12, 3, 7, 0, 12, 0,  7,  9, 6, 9, 7, 9, 6,
                                         1, 12, 1, 6, 6, 2,  12, 12, 2, 3, 0, 1, 12},
    /* 146:    1,       4,       7,  */ {6, 12, 4, 6,  9,  12, 9,  6,  1, 11, 1, 6, 1, 11,
                                         0, 12, 0, 11, 11, 8,  12, 12, 8, 4,  9, 0, 12},
    /* 162:    1,          5,    7,  */ {5,  12, 1, 12, 6,  1, 11, 1, 6, 1,  11, 0, 7, 0,
                                         11, 0,  7, 12, 12, 4, 0,  5, 4, 12, 12, 7, 6},
    /* 194:    1,             6, 7,  */ {5, 12, 7, 12, 9,  7,  0,  7, 9,  7,  0,  11, 1, 11,
                                         0, 11, 1, 12, 12, 10, 11, 5, 10, 12, 12, 1,  9},
    /*  28:       2, 3, 4,           */ {3, 12, 1, 12, 8,  1,  4,  1, 8,  1,  4,  10, 7, 10,
                                         4, 10, 7, 12, 12, 11, 10, 3, 11, 12, 12, 7,  8},
    /*  44:       2, 3,    5,        */ {3, 12, 1, 3,  9,  12, 9,  3,  4,  11, 4, 3, 4, 11,
                                         5, 12, 5, 11, 11, 10, 12, 12, 10, 1,  9, 5, 12},
    /*  52:       2,    4, 5,        */ {7, 12, 5, 7, 10, 12, 10, 7,  2, 8, 2,  7, 2, 8,
                                         1, 12, 1, 8, 8,  9,  12, 12, 9, 5, 10, 1, 12},
    /*  84:       2,    4,    6,     */ {6, 12, 2, 12, 7,  2, 8, 2, 7, 2,  8,  1, 4, 1,
                                         8, 1,  4, 12, 12, 5, 1, 6, 5, 12, 12, 4, 7},
    /* 148:       2,    4,       7,  */ {6, 12, 4, 12, 10, 4,  1, 4, 10, 4,  1,  8, 2, 8,
                                         1, 8,  2, 12, 12, 11, 8, 6, 11, 12, 12, 2, 10},
    /*  56:          3, 4, 5,        */ {7, 12, 5, 12, 11, 5, 2, 5, 11, 5,  2,  9, 3, 9,
                                         2, 9,  3, 12, 12, 8, 9, 7, 8,  12, 12, 3, 11},
    /* 104:          3,    5, 6,     */ {4, 12, 6, 4, 11, 12, 11, 4,  3,  9, 3,  4, 3, 9,
                                         2, 12, 2, 9, 9,  10, 12, 12, 10, 6, 11, 2, 12},
    /* 168:          3,    5,    7,  */ {7, 12, 3, 12, 4,  3, 9, 3, 4, 3,  9,  2, 5, 2,
                                         9, 2,  5, 12, 12, 6, 2, 7, 6, 12, 12, 5, 4},
    /*  87: 0, 1, 2,    4,    6,     */ {3, 12, 7, 3, 4, 12, 4,  3,  9, 2, 9, 3, 9, 2,
                                         5, 12, 5, 2, 2, 6,  12, 12, 6, 7, 4, 5, 12},
    /* 151: 0, 1, 2,    4,       7,  */ {6, 12, 4, 12, 11, 4,  3, 4, 11, 4,  3,  9, 2, 9,
                                         3, 9,  2, 12, 12, 10, 9, 6, 10, 12, 12, 2, 11},
    /* 199: 0, 1, 2,          6, 7,  */ {5, 12, 7, 5, 11, 12, 11, 5,  2, 9, 2,  5, 2, 9,
                                         3, 12, 3, 9, 9,  8,  12, 12, 8, 7, 11, 3, 12},
    /* 107: 0, 1,    3,    5, 6,     */ {4, 12, 6, 4, 10, 12, 10, 4,  1,  8, 1,  4, 1, 8,
                                         2, 12, 2, 8, 8,  11, 12, 12, 11, 6, 10, 2, 12},
    /* 171: 0, 1,    3,    5,    7,  */ {2, 12, 6, 2, 7, 12, 7,  2,  8, 1, 8, 2, 8, 1,
                                         4, 12, 4, 1, 1, 5,  12, 12, 5, 6, 7, 4, 12},
    /* 203: 0, 1,    3,       6, 7,  */ {5, 12, 7, 12, 10, 7, 2, 7, 10, 7,  2,  8, 1, 8,
                                         2, 8,  1, 12, 12, 9, 8, 5, 9,  12, 12, 1, 10},
    /* 211: 0, 1,       4,    6, 7,  */ {1, 12, 3, 12, 9,  3,  4,  3, 9,  3,  4,  11, 5, 11,
                                         4, 11, 5, 12, 12, 10, 11, 1, 10, 12, 12, 5,  9},
    /* 227: 0, 1,          5, 6, 7,  */ {1, 12, 3, 1,  8,  12, 8,  1,  4,  10, 4, 1, 4, 10,
                                         7, 12, 7, 10, 10, 11, 12, 12, 11, 3,  8, 7, 12},
    /*  61: 0,    2, 3, 4, 5,        */ {7, 12, 5, 7,  9,  12, 9,  7,  0,  11, 0, 7, 0, 11,
                                         1, 12, 1, 11, 11, 10, 12, 12, 10, 5,  9, 1, 12},
    /*  93: 0,    2, 3, 4,    6,     */ {1, 12, 5, 1, 6, 12, 6,  1,  11, 0, 11, 1, 11, 0,
                                         7, 12, 7, 0, 0, 4,  12, 12, 4,  5, 6,  7, 12},
    /* 109: 0,    2, 3,    5, 6,     */ {4, 12, 6, 12, 9,  6, 1,  6, 9, 6,  1,  11, 0, 11,
                                         1, 11, 0, 12, 12, 8, 11, 4, 8, 12, 12, 0,  9},
    /* 117: 0,    2,    4, 5, 6,     */ {3, 12, 7, 12, 0,  7, 9, 7, 0, 7,  9,  6, 1, 6,
                                         9, 6,  1, 12, 12, 2, 6, 3, 2, 12, 12, 1, 0},
    /* 213: 0,    2,    4,    6, 7,  */ {1,  12, 5, 12, 2,  5, 11, 5, 2, 5,  11, 4, 3, 4,
                                         11, 4,  3, 12, 12, 0, 4,  1, 0, 12, 12, 3, 2},
    /* 121: 0,       3, 4, 5, 6,     */ {0, 12, 2, 0, 11, 12, 11, 0,  7,  9, 7,  0, 7, 9,
                                         6, 12, 6, 9, 9,  10, 12, 12, 10, 2, 11, 6, 12},
    /* 233: 0,       3,    5, 6, 7,  */ {0, 12, 2, 12, 8,  2, 7,  2, 8, 2,  7,  10, 4, 10,
                                         7, 10, 4, 12, 12, 9, 10, 0, 9, 12, 12, 4,  8},
    /*  62:    1, 2, 3, 4, 5,        */ {7, 12, 5, 12, 8,  5,  0,  5, 8,  5,  0,  10, 3, 10,
                                         0, 10, 3, 12, 12, 11, 10, 7, 11, 12, 12, 3,  8},
    /* 158:    1, 2, 3, 4,       7,  */ {6, 12, 4, 6,  8,  12, 8,  6,  3, 10, 3, 6, 3, 10,
                                         0, 12, 0, 10, 10, 9,  12, 12, 9, 4,  8, 0, 12},
    /* 174:    1, 2, 3,    5,    7,  */ {0, 12, 4, 0, 5, 12, 5,  0,  10, 3, 10, 0, 10, 3,
                                         6, 12, 6, 3, 3, 7,  12, 12, 7,  4, 5,  6, 12},
    /* 182:    1, 2,    4, 5,    7,  */ {2, 12, 0, 12, 10, 0,  5, 0, 10, 0,  5,  8, 6, 8,
                                         5, 8,  6, 12, 12, 11, 8, 2, 11, 12, 12, 6, 10},
    /* 214:    1, 2,    4,    6, 7,  */ {2, 12, 0, 2,  9,  12, 9,  2,  5, 11, 5, 2, 5, 11,
                                         4, 12, 4, 11, 11, 8,  12, 12, 8, 0,  9, 4, 12},
    /* 186:    1,    3, 4, 5,    7,  */ {2, 12, 6, 12, 3,  6, 8, 6, 3, 6,  8,  5, 0, 5,
                                         8, 5,  0, 12, 12, 1, 5, 2, 1, 12, 12, 0, 3},
    /* 234:    1,    3,    5, 6, 7,  */ {0,  12, 4, 12, 1,  4, 10, 4, 1, 4,  10, 7, 2, 7,
                                         10, 7,  2, 12, 12, 3, 7,  0, 3, 12, 12, 2, 1},
    /* 124:       2, 3, 4, 5, 6,     */ {3, 12, 1, 12, 11, 1, 6, 1, 11, 1,  6,  9, 7, 9,
                                         6, 9,  7, 12, 12, 8, 9, 3, 8,  12, 12, 7, 11},
    /* 188:       2, 3, 4, 5,    7,  */ {3, 12, 1, 3, 10, 12, 10, 3,  6, 8, 6,  3, 6, 8,
                                         5, 12, 5, 8, 8,  9,  12, 12, 9, 1, 10, 5, 12},
};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 6.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling6_2[48][15] = {
    /*  67: 0, 1,             6,     */ {1, 10, 3, 6, 3, 10, 3, 6, 8, 5, 8, 6, 8, 5, 9},
    /* 131: 0, 1,                7,  */ {1, 11, 3, 11, 1, 6, 9, 6, 1, 6, 9, 7, 8, 7, 9},
    /*  21: 0,    2,    4,           */ {4, 1, 0, 1, 4, 10, 7, 10, 4, 10, 7, 2, 3, 2, 7},
    /*  69: 0,    2,          6,     */ {6, 3, 2, 3, 6, 8, 5, 8, 6, 8, 5, 0, 1, 0, 5},
    /*  41: 0,       3,    5,        */ {0, 9, 2, 5, 2, 9, 2, 5, 11, 4, 11, 5, 11, 4, 8},
    /*  73: 0,       3,       6,     */ {0, 10, 2, 10, 0, 5, 8, 5, 0, 5, 8, 6, 11, 6, 8},
    /*  81: 0,          4,    6,     */ {4, 5, 0, 10, 0, 5, 0, 10, 3, 6, 3, 10, 3, 6, 7},
    /*  97: 0,             5, 6,     */ {4, 8, 6, 3, 6, 8, 6, 3, 10, 0, 10, 3, 10, 0, 9},
    /* 193: 0,                6, 7,  */ {5, 8, 7, 8, 5, 0, 10, 0, 5, 0, 10, 3, 11, 3, 10},
    /*  22:    1, 2,    4,           */ {2, 8, 0, 8, 2, 7, 10, 7, 2, 7, 10, 4, 9, 4, 10},
    /* 134:    1, 2,             7,  */ {2, 11, 0, 7, 0, 11, 0, 7, 9, 6, 9, 7, 9, 6, 10},
    /*  42:    1,    3,    5,        */ {5, 2, 1, 2, 5, 11, 4, 11, 5, 11, 4, 3, 0, 3, 4},
    /* 138:    1,    3,          7,  */ {7, 0, 3, 0, 7, 9, 6, 9, 7, 9, 6, 1, 2, 1, 6},
    /* 146:    1,       4,       7,  */ {6, 9, 4, 9, 6, 1, 11, 1, 6, 1, 11, 0, 8, 0, 11},
    /* 162:    1,          5,    7,  */ {5, 6, 1, 11, 1, 6, 1, 11, 0, 7, 0, 11, 0, 7, 4},
    /* 194:    1,             6, 7,  */ {5, 9, 7, 0, 7, 9, 7, 0, 11, 1, 11, 0, 11, 1, 10},
    /*  28:       2, 3, 4,           */ {3, 8, 1, 4, 1, 8, 1, 4, 10, 7, 10, 4, 10, 7, 11},
    /*  44:       2, 3,    5,        */ {3, 9, 1, 9, 3, 4, 11, 4, 3, 4, 11, 5, 10, 5, 11},
    /*  52:       2,    4, 5,        */ {7, 10, 5, 10, 7, 2, 8, 2, 7, 2, 8, 1, 9, 1, 8},
    /*  84:       2,    4,    6,     */ {6, 7, 2, 8, 2, 7, 2, 8, 1, 4, 1, 8, 1, 4, 5},
    /* 148:       2,    4,       7,  */ {6, 10, 4, 1, 4, 10, 4, 1, 8, 2, 8, 1, 8, 2, 11},
    /*  56:          3, 4, 5,        */ {7, 11, 5, 2, 5, 11, 5, 2, 9, 3, 9, 2, 9, 3, 8},
    /* 104:          3,    5, 6,     */ {4, 11, 6, 11, 4, 3, 9, 3, 4, 3, 9, 2, 10, 2, 9},
    /* 168:          3,    5,    7,  */ {7, 4, 3, 9, 3, 4, 3, 9, 2, 5, 2, 9, 2, 5, 6},
    /*  87: 0, 1, 2,    4,    6,     */ {3, 4, 7, 4, 3, 9, 2, 9, 3, 9, 2, 5, 6, 5, 2},
    /* 151: 0, 1, 2,    4,       7,  */ {6, 11, 4, 3, 4, 11, 4, 3, 9, 2, 9, 3, 9, 2, 10},
    /* 199: 0, 1, 2,          6, 7,  */ {5, 11, 7, 11, 5, 2, 9, 2, 5, 2, 9, 3, 8, 3, 9},
    /* 107: 0, 1,    3,    5, 6,     */ {4, 10, 6, 10, 4, 1, 8, 1, 4, 1, 8, 2, 11, 2, 8},
    /* 171: 0, 1,    3,    5,    7,  */ {2, 7, 6, 7, 2, 8, 1, 8, 2, 8, 1, 4, 5, 4, 1},
    /* 203: 0, 1,    3,       6, 7,  */ {5, 10, 7, 2, 7, 10, 7, 2, 8, 1, 8, 2, 8, 1, 9},
    /* 211: 0, 1,       4,    6, 7,  */ {1, 9, 3, 4, 3, 9, 3, 4, 11, 5, 11, 4, 11, 5, 10},
    /* 227: 0, 1,          5, 6, 7,  */ {1, 8, 3, 8, 1, 4, 10, 4, 1, 4, 10, 7, 11, 7, 10},
    /*  61: 0,    2, 3, 4, 5,        */ {7, 9, 5, 9, 7, 0, 11, 0, 7, 0, 11, 1, 10, 1, 11},
    /*  93: 0,    2, 3, 4,    6,     */ {1, 6, 5, 6, 1, 11, 0, 11, 1, 11, 0, 7, 4, 7, 0},
    /* 109: 0,    2, 3,    5, 6,     */ {4, 9, 6, 1, 6, 9, 6, 1, 11, 0, 11, 1, 11, 0, 8},
    /* 117: 0,    2,    4, 5, 6,     */ {3, 0, 7, 9, 7, 0, 7, 9, 6, 1, 6, 9, 6, 1, 2},
    /* 213: 0,    2,    4,    6, 7,  */ {1, 2, 5, 11, 5, 2, 5, 11, 4, 3, 4, 11, 4, 3, 0},
    /* 121: 0,       3, 4, 5, 6,     */ {0, 11, 2, 11, 0, 7, 9, 7, 0, 7, 9, 6, 10, 6, 9},
    /* 233: 0,       3,    5, 6, 7,  */ {0, 8, 2, 7, 2, 8, 2, 7, 10, 4, 10, 7, 10, 4, 9},
    /*  62:    1, 2, 3, 4, 5,        */ {7, 8, 5, 0, 5, 8, 5, 0, 10, 3, 10, 0, 10, 3, 11},
    /* 158:    1, 2, 3, 4,       7,  */ {6, 8, 4, 8, 6, 3, 10, 3, 6, 3, 10, 0, 9, 0, 10},
    /* 174:    1, 2, 3,    5,    7,  */ {0, 5, 4, 5, 0, 10, 3, 10, 0, 10, 3, 6, 7, 6, 3},
    /* 182:    1, 2,    4, 5,    7,  */ {2, 10, 0, 5, 0, 10, 0, 5, 8, 6, 8, 5, 8, 6, 11},
    /* 214:    1, 2,    4,    6, 7,  */ {2, 9, 0, 9, 2, 5, 11, 5, 2, 5, 11, 4, 8, 4, 11},
    /* 186:    1,    3, 4, 5,    7,  */ {2, 3, 6, 8, 6, 3, 6, 8, 5, 0, 5, 8, 5, 0, 1},
    /* 234:    1,    3,    5, 6, 7,  */ {0, 1, 4, 10, 4, 1, 4, 10, 7, 2, 7, 10, 7, 2, 3},
    /* 124:       2, 3, 4, 5, 6,     */ {3, 11, 1, 6, 1, 11, 1, 6, 9, 7, 9, 6, 9, 7, 8},
    /* 188:       2, 3, 4, 5,    7,  */ {3, 10, 1, 10, 3, 6, 8, 6, 3, 6, 8, 5, 9, 5, 8}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief test table for case 7
 * 3 faces to test + eventually the interior
 * When the tests on the 3 specified faces are positive :
 * - if the test on the interior is positive : 5 first triangles
 * - if the test on the interior is negative : 9 next triangles
 * When the tests on the first  and the second specified faces are positive : 9 next triangles
 * When the tests on the first  and the third  specified faces are positive : 9 next triangles
 * When the tests on the second and the third  specified faces are positive : 9 next triangles
 * When the test on the first  specified face is positive : 5 next triangles
 * When the test on the second specified face is positive : 5 next triangles
 * When the test on the third  specified face is positive : 5 next triangles
 * When the tests on the 3 specified faces are negative : 3 last triangles
 * The support edge for the interior test is marked as the 5th column.
 *
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::test7[16][5] = {
    /*  37: 0,    2,       5,        */ {1, 2, 5, 7, 1},
    /* 133: 0,    2,             7,  */ {3, 4, 5, 7, 3},
    /* 161: 0,             5,    7,  */ {4, 1, 6, 7, 4},
    /*  26:    1,    3, 4,           */ {4, 1, 5, 7, 0},
    /*  74:    1,    3,       6,     */ {2, 3, 5, 7, 2},
    /*  82:    1,       4,    6,     */ {1, 2, 6, 7, 5},
    /* 164:       2,       5,    7,  */ {2, 3, 6, 7, 6},
    /*  88:          3, 4,    6,     */ {3, 4, 6, 7, 7},
    /* 167: 0, 1, 2,       5,    7,  */ {-3, -4, -6, -7, 7},
    /*  91: 0, 1,    3, 4,    6,     */ {-2, -3, -6, -7, 6},
    /* 173: 0,    2, 3,    5,    7,  */ {-1, -2, -6, -7, 5},
    /* 181: 0,    2,    4, 5,    7,  */ {-2, -3, -5, -7, 2},
    /* 229: 0,    2,       5, 6, 7,  */ {-4, -1, -5, -7, 0},
    /*  94:    1, 2, 3, 4,    6,     */ {-4, -1, -6, -7, 4},
    /* 122:    1,    3, 4, 5, 6,     */ {-3, -4, -5, -7, 3},
    /* 218:    1,    3, 4,    6, 7,  */ {-1, -2, -5, -7, 1}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 7.1
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling7_1[16][9] = {
    /*  37: 0,    2,       5,        */ {9, 5, 4, 10, 1, 2, 8, 3, 0},
    /* 133: 0,    2,             7,  */ {11, 7, 6, 8, 3, 0, 10, 1, 2},
    /* 161: 0,             5,    7,  */ {3, 0, 8, 5, 4, 9, 7, 6, 11},
    /*  26:    1,    3, 4,           */ {8, 4, 7, 9, 0, 1, 11, 2, 3},
    /*  74:    1,    3,       6,     */ {10, 6, 5, 11, 2, 3, 9, 0, 1},
    /*  82:    1,       4,    6,     */ {0, 1, 9, 6, 5, 10, 4, 7, 8},
    /* 164:       2,       5,    7,  */ {1, 2, 10, 7, 6, 11, 5, 4, 9},
    /*  88:          3, 4,    6,     */ {2, 3, 11, 4, 7, 8, 6, 5, 10},
    /* 167: 0, 1, 2,       5,    7,  */ {11, 3, 2, 8, 7, 4, 10, 5, 6},
    /*  91: 0, 1,    3, 4,    6,     */ {10, 2, 1, 11, 6, 7, 9, 4, 5},
    /* 173: 0,    2, 3,    5,    7,  */ {9, 1, 0, 10, 5, 6, 8, 7, 4},
    /* 181: 0,    2,    4, 5,    7,  */ {5, 6, 10, 3, 2, 11, 1, 0, 9},
    /* 229: 0,    2,       5, 6, 7,  */ {7, 4, 8, 1, 0, 9, 3, 2, 11},
    /*  94:    1, 2, 3, 4,    6,     */ {8, 0, 3, 9, 4, 5, 11, 6, 7},
    /* 122:    1,    3, 4, 5, 6,     */ {6, 7, 11, 0, 3, 8, 2, 1, 10},
    /* 218:    1,    3, 4,    6, 7,  */ {4, 5, 9, 2, 1, 10, 0, 3, 8}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 7.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling7_2[16][3][15] = {
    /*  37: 0,    2,       5,        */ {/* 1,0 */ {1, 2, 10, 3, 4, 8, 4, 3, 5, 0, 5, 3, 5, 0, 9},
                                         /* 0,1 */ {3, 0, 8, 9, 1, 4, 2, 4, 1, 4, 2, 5, 10, 5, 2},
                                         /* 1,1 */ {9, 5, 4, 0, 10, 1, 10, 0, 8, 10, 8, 2, 3, 2, 8}},
    /* 133: 0,    2,             7,  */
    {/* 1,0 */ {3, 0, 8, 1, 6, 10, 6, 1, 7, 2, 7, 1, 7, 2, 11},
     /* 0,1 */ {1, 2, 10, 11, 3, 6, 0, 6, 3, 6, 0, 7, 8, 7, 0},
     /* 1,1 */ {11, 7, 6, 2, 8, 3, 8, 2, 10, 8, 10, 0, 1, 0, 10}},
    /* 161: 0,             5,    7,  */
    {/* 1,0 */ {9, 5, 4, 11, 3, 6, 0, 6, 3, 6, 0, 7, 8, 7, 0},
     /* 0,1 */ {11, 7, 6, 3, 4, 8, 4, 3, 5, 0, 5, 3, 5, 0, 9},
     /* 1,1 */ {3, 0, 8, 4, 9, 7, 11, 7, 9, 5, 11, 9, 11, 5, 6}},
    /*  26:    1,    3, 4,           */
    {/* 1,0 */ {0, 1, 9, 2, 7, 11, 7, 2, 4, 3, 4, 2, 4, 3, 8},
     /* 0,1 */ {2, 3, 11, 8, 0, 7, 1, 7, 0, 7, 1, 4, 9, 4, 1},
     /* 1,1 */ {8, 4, 7, 3, 9, 0, 9, 3, 11, 9, 11, 1, 2, 1, 11}},
    /*  74:    1,    3,       6,     */
    {/* 1,0 */ {2, 3, 11, 0, 5, 9, 5, 0, 6, 1, 6, 0, 6, 1, 10},
     /* 0,1 */ {0, 1, 9, 10, 2, 5, 3, 5, 2, 5, 3, 6, 11, 6, 3},
     /* 1,1 */ {6, 5, 10, 1, 11, 2, 11, 1, 9, 11, 9, 3, 0, 3, 9}},
    /*  82:    1,       4,    6,     */
    {/* 1,0 */ {6, 5, 10, 8, 0, 7, 1, 7, 0, 7, 1, 4, 9, 4, 1},
     /* 0,1 */ {8, 4, 7, 0, 5, 9, 5, 0, 6, 1, 6, 0, 6, 1, 10},
     /* 1,1 */ {0, 1, 9, 5, 10, 4, 8, 4, 10, 6, 8, 10, 8, 6, 7}},
    /* 164:       2,       5,    7,  */
    {/* 1,0 */ {11, 7, 6, 9, 1, 4, 2, 4, 1, 4, 2, 5, 10, 5, 2},
     /* 0,1 */ {9, 5, 4, 1, 6, 10, 6, 1, 7, 2, 7, 1, 7, 2, 11},
     /* 1,1 */ {1, 2, 10, 6, 11, 5, 9, 5, 11, 7, 9, 11, 9, 7, 4}},
    /*  88:          3, 4,    6,     */
    {/* 1,0 */ {8, 4, 7, 10, 2, 5, 3, 5, 2, 5, 3, 6, 11, 6, 3},
     /* 0,1 */ {6, 5, 10, 2, 7, 11, 7, 2, 4, 3, 4, 2, 4, 3, 8},
     /* 1,1 */ {2, 3, 11, 7, 8, 6, 10, 6, 8, 4, 10, 8, 10, 4, 5}},
    /* 167: 0, 1, 2,       5,    7,  */
    {/* 1,0 */ {7, 4, 8, 5, 2, 10, 2, 5, 3, 6, 3, 5, 3, 6, 11},
     /* 0,1 */ {10, 5, 6, 11, 7, 2, 4, 2, 7, 2, 4, 3, 8, 3, 4},
     /* 1,1 */ {11, 3, 2, 6, 8, 7, 8, 6, 10, 8, 10, 4, 5, 4, 10}},
    /*  91: 0, 1,    3, 4,    6,     */
    {/* 1,0 */ {6, 7, 11, 4, 1, 9, 1, 4, 2, 5, 2, 4, 2, 5, 10},
     /* 0,1 */ {4, 5, 9, 10, 6, 1, 7, 1, 6, 1, 7, 2, 11, 2, 7},
     /* 1,1 */ {10, 2, 1, 5, 11, 6, 11, 5, 9, 11, 9, 7, 4, 7, 9}},
    /* 173: 0,    2, 3,    5,    7,  */
    {/* 1,0 */ {10, 5, 6, 7, 0, 8, 0, 7, 1, 4, 1, 7, 1, 4, 9},
     /* 0,1 */ {7, 4, 8, 9, 5, 0, 6, 0, 5, 0, 6, 1, 10, 1, 6},
     /* 1,1 */ {9, 1, 0, 4, 10, 5, 10, 4, 8, 10, 8, 6, 7, 6, 8}},
    /* 181: 0,    2,    4, 5,    7,  */
    {/* 1,0 */ {11, 3, 2, 9, 5, 0, 6, 0, 5, 0, 6, 1, 10, 1, 6},
     /* 0,1 */ {9, 1, 0, 5, 2, 10, 2, 5, 3, 6, 3, 5, 3, 6, 11},
     /* 1,1 */ {10, 5, 6, 2, 11, 1, 9, 1, 11, 3, 9, 11, 9, 3, 0}},
    /* 229: 0,    2,       5, 6, 7,  */
    {/* 1,0 */ {9, 1, 0, 11, 7, 2, 4, 2, 7, 2, 4, 3, 8, 3, 4},
     /* 0,1 */ {11, 3, 2, 7, 0, 8, 0, 7, 1, 4, 1, 7, 1, 4, 9},
     /* 1,1 */ {7, 4, 8, 0, 9, 3, 11, 3, 9, 1, 11, 9, 11, 1, 2}},
    /*  94:    1, 2, 3, 4,    6,     */
    {/* 1,0 */ {4, 5, 9, 6, 3, 11, 3, 6, 0, 7, 0, 6, 0, 7, 8},
     /* 0,1 */ {6, 7, 11, 8, 4, 3, 5, 3, 4, 3, 5, 0, 9, 0, 5},
     /* 1,1 */ {8, 0, 3, 7, 9, 4, 9, 7, 11, 9, 11, 5, 6, 5, 11}},
    /* 122:    1,    3, 4, 5, 6,     */
    {/* 1,0 */ {8, 0, 3, 10, 6, 1, 7, 1, 6, 1, 7, 2, 11, 2, 7},
     /* 0,1 */ {10, 2, 1, 6, 3, 11, 3, 6, 0, 7, 0, 6, 0, 7, 8},
     /* 1,1 */ {6, 7, 11, 3, 8, 2, 10, 2, 8, 0, 10, 8, 10, 0, 1}},
    /* 218:    1,    3, 4,    6, 7,  */
    {/* 1,0 */ {10, 2, 1, 8, 4, 3, 5, 3, 4, 3, 5, 0, 9, 0, 5},
     /* 0,1 */ {8, 0, 3, 4, 1, 9, 1, 4, 2, 5, 2, 4, 2, 5, 10},
     /* 1,1 */ {4, 5, 9, 1, 10, 0, 8, 0, 10, 2, 8, 10, 8, 2, 3}}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 7.3
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling7_3[16][3][27] = {
    /*  37: 0,    2,       5,        */ {
        /* 1,0 */ {12, 2, 10, 12, 10, 5, 12, 5, 4, 12, 4, 8, 12, 8, 3, 12, 3, 0, 12, 0, 9, 12, 9, 1, 12, 1, 2},
        /* 0,1 */ {12, 5, 4, 12, 4, 8, 12, 8, 3, 12, 3, 2, 12, 2, 10, 12, 10, 1, 12, 1, 0, 12, 0, 9, 12, 9, 5},
        /* 1,1 */ {5, 4, 12, 10, 5, 12, 2, 10, 12, 3, 2, 12, 8, 3, 12, 0, 8, 12, 1, 0, 12, 9, 1, 12, 4, 9, 12}},
    /* 133: 0,    2,             7,  */
    {/* 1,0 */ {12, 0, 8, 12, 8, 7, 12, 7, 6, 12, 6, 10, 12, 10, 1, 12, 1, 2, 12, 2, 11, 12, 11, 3, 12, 3, 0},
     /* 0,1 */ {12, 7, 6, 12, 6, 10, 12, 10, 1, 12, 1, 0, 12, 0, 8, 12, 8, 3, 12, 3, 2, 12, 2, 11, 12, 11, 7},
     /* 1,1 */ {7, 6, 12, 8, 7, 12, 0, 8, 12, 1, 0, 12, 10, 1, 12, 2, 10, 12, 3, 2, 12, 11, 3, 12, 6, 11, 12}},
    /* 161: 0,             5,    7,  */
    {/* 1,0 */ {9, 5, 12, 0, 9, 12, 3, 0, 12, 11, 3, 12, 6, 11, 12, 7, 6, 12, 8, 7, 12, 4, 8, 12, 5, 4, 12},
     /* 0,1 */ {3, 0, 12, 11, 3, 12, 6, 11, 12, 5, 6, 12, 9, 5, 12, 4, 9, 12, 7, 4, 12, 8, 7, 12, 0, 8, 12},
     /* 1,1 */ {12, 3, 0, 12, 0, 9, 12, 9, 5, 12, 5, 6, 12, 6, 11, 12, 11, 7, 12, 7, 4, 12, 4, 8, 12, 8, 3}},
    /*  26:    1,    3, 4,           */
    {/* 1,0 */ {12, 1, 9, 12, 9, 4, 12, 4, 7, 12, 7, 11, 12, 11, 2, 12, 2, 3, 12, 3, 8, 12, 8, 0, 12, 0, 1},
     /* 0,1 */ {12, 4, 7, 12, 7, 11, 12, 11, 2, 12, 2, 1, 12, 1, 9, 12, 9, 0, 12, 0, 3, 12, 3, 8, 12, 8, 4},
     /* 1,1 */ {4, 7, 12, 9, 4, 12, 1, 9, 12, 2, 1, 12, 11, 2, 12, 3, 11, 12, 0, 3, 12, 8, 0, 12, 7, 8, 12}},
    /*  74:    1,    3,       6,     */
    {/* 1,0 */ {12, 3, 11, 12, 11, 6, 12, 6, 5, 12, 5, 9, 12, 9, 0, 12, 0, 1, 12, 1, 10, 12, 10, 2, 12, 2, 3},
     /* 0,1 */ {12, 6, 5, 12, 5, 9, 12, 9, 0, 12, 0, 3, 12, 3, 11, 12, 11, 2, 12, 2, 1, 12, 1, 10, 12, 10, 6},
     /* 1,1 */ {6, 5, 12, 11, 6, 12, 3, 11, 12, 0, 3, 12, 9, 0, 12, 1, 9, 12, 2, 1, 12, 10, 2, 12, 5, 10, 12}},
    /*  82:    1,       4,    6,     */
    {/* 1,0 */ {10, 6, 12, 1, 10, 12, 0, 1, 12, 8, 0, 12, 7, 8, 12, 4, 7, 12, 9, 4, 12, 5, 9, 12, 6, 5, 12},
     /* 0,1 */ {0, 1, 12, 8, 0, 12, 7, 8, 12, 6, 7, 12, 10, 6, 12, 5, 10, 12, 4, 5, 12, 9, 4, 12, 1, 9, 12},
     /* 1,1 */ {12, 0, 1, 12, 1, 10, 12, 10, 6, 12, 6, 7, 12, 7, 8, 12, 8, 4, 12, 4, 5, 12, 5, 9, 12, 9, 0}},
    /* 164:       2,       5,    7,  */
    {/* 1,0 */ {11, 7, 12, 2, 11, 12, 1, 2, 12, 9, 1, 12, 4, 9, 12, 5, 4, 12, 10, 5, 12, 6, 10, 12, 7, 6, 12},
     /* 0,1 */ {1, 2, 12, 9, 1, 12, 4, 9, 12, 7, 4, 12, 11, 7, 12, 6, 11, 12, 5, 6, 12, 10, 5, 12, 2, 10, 12},
     /* 1,1 */ {12, 1, 2, 12, 2, 11, 12, 11, 7, 12, 7, 4, 12, 4, 9, 12, 9, 5, 12, 5, 6, 12, 6, 10, 12, 10, 1}},
    /*  88:          3, 4,    6,     */
    {/* 1,0 */ {8, 4, 12, 3, 8, 12, 2, 3, 12, 10, 2, 12, 5, 10, 12, 6, 5, 12, 11, 6, 12, 7, 11, 12, 4, 7, 12},
     /* 0,1 */ {2, 3, 12, 10, 2, 12, 5, 10, 12, 4, 5, 12, 8, 4, 12, 7, 8, 12, 6, 7, 12, 11, 6, 12, 3, 11, 12},
     /* 1,1 */ {12, 2, 3, 12, 3, 8, 12, 8, 4, 12, 4, 5, 12, 5, 10, 12, 10, 6, 12, 6, 7, 12, 7, 11, 12, 11, 2}},
    /* 167: 0, 1, 2,       5,    7,  */
    {/* 1,0 */ {12, 4, 8, 12, 8, 3, 12, 3, 2, 12, 2, 10, 12, 10, 5, 12, 5, 6, 12, 6, 11, 12, 11, 7, 12, 7, 4},
     /* 0,1 */ {12, 3, 2, 12, 2, 10, 12, 10, 5, 12, 5, 4, 12, 4, 8, 12, 8, 7, 12, 7, 6, 12, 6, 11, 12, 11, 3},
     /* 1,1 */ {3, 2, 12, 8, 3, 12, 4, 8, 12, 5, 4, 12, 10, 5, 12, 6, 10, 12, 7, 6, 12, 11, 7, 12, 2, 11, 12}},
    /*  91: 0, 1,    3, 4,    6,     */
    {/* 1,0 */ {12, 7, 11, 12, 11, 2, 12, 2, 1, 12, 1, 9, 12, 9, 4, 12, 4, 5, 12, 5, 10, 12, 10, 6, 12, 6, 7},
     /* 0,1 */ {12, 2, 1, 12, 1, 9, 12, 9, 4, 12, 4, 7, 12, 7, 11, 12, 11, 6, 12, 6, 5, 12, 5, 10, 12, 10, 2},
     /* 1,1 */ {2, 1, 12, 11, 2, 12, 7, 11, 12, 4, 7, 12, 9, 4, 12, 5, 9, 12, 6, 5, 12, 10, 6, 12, 1, 10, 12}},
    /* 173: 0,    2, 3,    5,    7,  */
    {/* 1,0 */ {12, 6, 10, 12, 10, 1, 12, 1, 0, 12, 0, 8, 12, 8, 7, 12, 7, 4, 12, 4, 9, 12, 9, 5, 12, 5, 6},
     /* 0,1 */ {12, 1, 0, 12, 0, 8, 12, 8, 7, 12, 7, 6, 12, 6, 10, 12, 10, 5, 12, 5, 4, 12, 4, 9, 12, 9, 1},
     /* 1,1 */ {1, 0, 12, 10, 1, 12, 6, 10, 12, 7, 6, 12, 8, 7, 12, 4, 8, 12, 5, 4, 12, 9, 5, 12, 0, 9, 12}},
    /* 181: 0,    2,    4, 5,    7,  */
    {/* 1,0 */ {11, 3, 12, 6, 11, 12, 5, 6, 12, 9, 5, 12, 0, 9, 12, 1, 0, 12, 10, 1, 12, 2, 10, 12, 3, 2, 12},
     /* 0,1 */ {5, 6, 12, 9, 5, 12, 0, 9, 12, 3, 0, 12, 11, 3, 12, 2, 11, 12, 1, 2, 12, 10, 1, 12, 6, 10, 12},
     /* 1,1 */ {12, 5, 6, 12, 6, 11, 12, 11, 3, 12, 3, 0, 12, 0, 9, 12, 9, 1, 12, 1, 2, 12, 2, 10, 12, 10, 5}},
    /* 229: 0,    2,       5, 6, 7,  */
    {/* 1,0 */ {9, 1, 12, 4, 9, 12, 7, 4, 12, 11, 7, 12, 2, 11, 12, 3, 2, 12, 8, 3, 12, 0, 8, 12, 1, 0, 12},
     /* 0,1 */ {7, 4, 12, 11, 7, 12, 2, 11, 12, 1, 2, 12, 9, 1, 12, 0, 9, 12, 3, 0, 12, 8, 3, 12, 4, 8, 12},
     /* 1,1 */ {12, 7, 4, 12, 4, 9, 12, 9, 1, 12, 1, 2, 12, 2, 11, 12, 11, 3, 12, 3, 0, 12, 0, 8, 12, 8, 7}},
    /*  94:    1, 2, 3, 4,    6,     */
    {/* 1,0 */ {12, 5, 9, 12, 9, 0, 12, 0, 3, 12, 3, 11, 12, 11, 6, 12, 6, 7, 12, 7, 8, 12, 8, 4, 12, 4, 5},
     /* 0,1 */ {12, 0, 3, 12, 3, 11, 12, 11, 6, 12, 6, 5, 12, 5, 9, 12, 9, 4, 12, 4, 7, 12, 7, 8, 12, 8, 0},
     /* 1,1 */ {0, 3, 12, 9, 0, 12, 5, 9, 12, 6, 5, 12, 11, 6, 12, 7, 11, 12, 4, 7, 12, 8, 4, 12, 3, 8, 12}},
    /* 122:    1,    3, 4, 5, 6,     */
    {/* 1,0 */ {8, 0, 12, 7, 8, 12, 6, 7, 12, 10, 6, 12, 1, 10, 12, 2, 1, 12, 11, 2, 12, 3, 11, 12, 0, 3, 12},
     /* 0,1 */ {6, 7, 12, 10, 6, 12, 1, 10, 12, 0, 1, 12, 8, 0, 12, 3, 8, 12, 2, 3, 12, 11, 2, 12, 7, 11, 12},
     /* 1,1 */ {12, 6, 7, 12, 7, 8, 12, 8, 0, 12, 0, 1, 12, 1, 10, 12, 10, 2, 12, 2, 3, 12, 3, 11, 12, 11, 6}},
    /* 218:    1,    3, 4,    6, 7,  */
    {/* 1,0 */ {10, 2, 12, 5, 10, 12, 4, 5, 12, 8, 4, 12, 3, 8, 12, 0, 3, 12, 9, 0, 12, 1, 9, 12, 2, 1, 12},
     /* 0,1 */ {4, 5, 12, 8, 4, 12, 3, 8, 12, 2, 3, 12, 10, 2, 12, 1, 10, 12, 0, 1, 12, 9, 0, 12, 5, 9, 12},
     /* 1,1 */ {12, 4, 5, 12, 5, 10, 12, 10, 2, 12, 2, 3, 12, 3, 8, 12, 8, 0, 12, 0, 1, 12, 1, 9, 12, 9, 4}}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 7.4.1
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling7_4_1[16][15] = {
    /*  37: 0,    2,       5,        */ {3, 4, 8, 4, 3, 10, 2, 10, 3, 4, 10, 5, 9, 1, 0},
    /* 133: 0,    2,             7,  */ {1, 6, 10, 6, 1, 8, 0, 8, 1, 6, 8, 7, 11, 3, 2},
    /* 161: 0,             5,    7,  */ {11, 3, 6, 9, 6, 3, 6, 9, 5, 0, 9, 3, 7, 4, 8},
    /*  26:    1,    3, 4,           */ {2, 7, 11, 7, 2, 9, 1, 9, 2, 7, 9, 4, 8, 0, 3},
    /*  74:    1,    3,       6,     */ {0, 5, 9, 5, 0, 11, 3, 11, 0, 5, 11, 6, 10, 2, 1},
    /*  82:    1,       4,    6,     */ {8, 0, 7, 10, 7, 0, 7, 10, 6, 1, 10, 0, 4, 5, 9},
    /* 164:       2,       5,    7,  */ {9, 1, 4, 11, 4, 1, 4, 11, 7, 2, 11, 1, 5, 6, 10},
    /*  88:          3, 4,    6,     */ {10, 2, 5, 8, 5, 2, 5, 8, 4, 3, 8, 2, 6, 7, 11},
    /* 167: 0, 1, 2,       5,    7,  */ {5, 2, 10, 2, 5, 8, 4, 8, 5, 2, 8, 3, 11, 7, 6},
    /*  91: 0, 1,    3, 4,    6,     */ {4, 1, 9, 1, 4, 11, 7, 11, 4, 1, 11, 2, 10, 6, 5},
    /* 173: 0,    2, 3,    5,    7,  */ {7, 0, 8, 0, 7, 10, 6, 10, 7, 0, 10, 1, 9, 5, 4},
    /* 181: 0,    2,    4, 5,    7,  */ {9, 5, 0, 11, 0, 5, 0, 11, 3, 6, 11, 5, 1, 2, 10},
    /* 229: 0,    2,       5, 6, 7,  */ {11, 7, 2, 9, 2, 7, 2, 9, 1, 4, 9, 7, 3, 0, 8},
    /*  94:    1, 2, 3, 4,    6,     */ {6, 3, 11, 3, 6, 9, 5, 9, 6, 3, 9, 0, 8, 4, 7},
    /* 122:    1,    3, 4, 5, 6,     */ {10, 6, 1, 8, 1, 6, 1, 8, 0, 7, 8, 6, 2, 3, 11},
    /* 218:    1,    3, 4,    6, 7,  */ {8, 4, 3, 10, 3, 4, 3, 10, 2, 5, 10, 4, 0, 1, 9}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 7.4.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling7_4_2[16][27] = {
    /*  37: 0,    2,       5,        */ {9, 4, 8, 4, 9, 5, 10, 5, 9, 1, 10, 9, 10, 1,
                                         2, 0, 2, 1, 2, 0, 3,  8, 3, 0, 9,  8, 0},
    /* 133: 0,    2,             7,  */ {11, 6, 10, 6, 11, 7, 8, 7,  11, 3, 8,  11, 8, 3,
                                         0,  2, 0,  3, 0,  2, 1, 10, 1,  2, 11, 10, 2},
    /* 161: 0,             5,    7,  */ {11, 3, 8, 0, 8, 3, 8, 0, 9, 8,  9, 4,  5, 4,
                                         9,  4, 5, 7, 6, 7, 5, 7, 6, 11, 7, 11, 8},
    /*  26:    1,    3, 4,           */ {8, 7, 11, 7, 8, 4, 9, 4,  8, 0, 9, 8,  9, 0,
                                         1, 3, 1,  0, 1, 3, 2, 11, 2, 3, 8, 11, 3},
    /*  74:    1,    3,       6,     */ {10, 5, 9, 5, 10, 6, 11, 6, 10, 2, 11, 10, 11, 2,
                                         3,  1, 3, 2, 3,  1, 0,  9, 0,  1, 10, 9,  1},
    /*  82:    1,       4,    6,     */ {8,  0, 9, 1, 9, 0, 9, 1, 10, 9, 10, 5, 6, 5,
                                         10, 5, 6, 4, 7, 4, 6, 4, 7,  8, 4,  8, 9},
    /* 164:       2,       5,    7,  */ {9,  1, 10, 2, 10, 1, 10, 2, 11, 10, 11, 6, 7, 6,
                                         11, 6, 7,  5, 4,  5, 7,  5, 4,  9,  5,  9, 10},
    /*  88:          3, 4,    6,     */ {10, 2, 11, 3, 11, 2, 11, 3, 8, 11, 8, 7,  4, 7,
                                         8,  7, 4,  6, 5,  6, 4,  6, 5, 10, 6, 10, 11},
    /* 167: 0, 1, 2,       5,    7,  */ {11, 2, 10, 2, 11, 3, 8, 3,  11, 7, 8,  11, 8, 7,
                                         4,  6, 4,  7, 4,  6, 5, 10, 5,  6, 11, 10, 6},
    /*  91: 0, 1,    3, 4,    6,     */ {10, 1, 9, 1, 10, 2, 11, 2, 10, 6, 11, 10, 11, 6,
                                         7,  5, 7, 6, 7,  5, 4,  9, 4,  5, 10, 9,  5},
    /* 173: 0,    2, 3,    5,    7,  */ {9, 0, 8, 0, 9, 1, 10, 1, 9, 5, 10, 9, 10, 5,
                                         6, 4, 6, 5, 6, 4, 7,  8, 7, 4, 9,  8, 4},
    /* 181: 0,    2,    4, 5,    7,  */ {9,  5, 10, 6, 10, 5, 10, 6, 11, 10, 11, 2, 3, 2,
                                         11, 2, 3,  1, 0,  1, 3,  1, 0,  9,  1,  9, 10},
    /* 229: 0,    2,       5, 6, 7,  */ {11, 7, 8, 4, 8, 7, 8, 4, 9, 8,  9, 0,  1, 0,
                                         9,  0, 1, 3, 2, 3, 1, 3, 2, 11, 3, 11, 8},
    /*  94:    1, 2, 3, 4,    6,     */ {8, 3, 11, 3, 8, 0, 9, 0,  8, 4, 9, 8,  9, 4,
                                         5, 7, 5,  4, 5, 7, 6, 11, 6, 7, 8, 11, 7},
    /* 122:    1,    3, 4, 5, 6,     */ {10, 6, 11, 7, 11, 6, 11, 7, 8, 11, 8, 3,  0, 3,
                                         8,  3, 0,  2, 1,  2, 0,  2, 1, 10, 2, 10, 11},
    /* 218:    1,    3, 4,    6, 7,  */ {8,  4, 9, 5, 9, 4, 9, 5, 10, 9, 10, 1, 2, 1,
                                         10, 1, 2, 0, 3, 0, 2, 0, 3,  8, 0,  8, 9}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief tiling table for case 8
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling8[6][6] = {
    /*  15: 0, 1, 2, 3,              */ {9, 8, 10, 10, 8, 11},
    /*  51: 0, 1,       4, 5,        */ {1, 5, 3, 3, 5, 7},
    /* 153: 0,       3, 4,       7,  */ {0, 4, 2, 4, 6, 2},
    /* 102:    1, 2,       5, 6,     */ {0, 2, 4, 4, 2, 6},
    /* 204:       2, 3,       6, 7,  */ {1, 3, 5, 3, 7, 5},
    /* 240:             4, 5, 6, 7,  */ {9, 10, 8, 10, 11, 8}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief tiling table for case 9
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling9[8][12] = {
    /*  39: 0, 1, 2,       5,        */ {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8},
    /*  27: 0, 1,    3, 4,           */ {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1},
    /* 141: 0,    2, 3,          7,  */ {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8},
    /* 177: 0,          4, 5,    7,  */ {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5},
    /*  78:    1, 2, 3,       6,     */ {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9},
    /* 114:    1,       4, 5, 6,     */ {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0},
    /* 228:       2,       5, 6, 7,  */ {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2},
    /* 216:          3, 4,    6, 7,  */ {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief test table for case 10
 * 2 faces to test + eventually the interior
 * When the tests on both specified faces are positive : 4 middle triangles (1)
 * When the test on the first  specified face is positive : 8 first triangles
 * When the test on the second specified face is positive : 8 next triangles
 * When the tests on both specified faces are negative :
 * - if the test on the interior is negative : 4 middle triangles
 * - if the test on the interior is positive : 8 last triangles
 *
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::test10[6][3] = {
    /* 195: 0, 1,             6, 7,  */ {2, 4, 7},
    /*  85: 0,    2,    4,    6,     */ {5, 6, 7},
    /* 105: 0,       3,    5, 6,     */ {1, 3, 7},
    /* 150:    1, 2,    4,       7,  */ {1, 3, 7},
    /* 170:    1,    3,    5,    7,  */ {5, 6, 7},
    /*  60:       2, 3, 4, 5,        */ {2, 4, 7}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 10.1.1
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling10_1_1[6][12] = {
    /* 195: 0, 1,             6, 7,  */ {5, 10, 7, 11, 7, 10, 8, 1, 9, 1, 8, 3},
    /*  85: 0,    2,    4,    6,     */ {1, 2, 5, 6, 5, 2, 4, 3, 0, 3, 4, 7},
    /* 105: 0,       3,    5, 6,     */ {11, 0, 8, 0, 11, 2, 4, 9, 6, 10, 6, 9},
    /* 150:    1, 2,    4,       7,  */ {9, 0, 10, 2, 10, 0, 6, 8, 4, 8, 6, 11},
    /* 170:    1,    3,    5,    7,  */ {7, 2, 3, 2, 7, 6, 0, 1, 4, 5, 4, 1},
    /*  60:       2, 3, 4, 5,        */ {7, 9, 5, 9, 7, 8, 10, 1, 11, 3, 11, 1}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 10.1.1 inverted
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling10_1_1_[6][12] = {
    /* 195: 0, 1,             6, 7,  */ {5, 9, 7, 8, 7, 9, 11, 1, 10, 1, 11, 3},
    /*  85: 0,    2,    4,    6,     */ {3, 2, 7, 6, 7, 2, 4, 1, 0, 1, 4, 5},
    /* 105: 0,       3,    5, 6,     */ {10, 0, 9, 0, 10, 2, 4, 8, 6, 11, 6, 8},
    /* 150:    1, 2,    4,       7,  */ {8, 0, 11, 2, 11, 0, 6, 9, 4, 9, 6, 10},
    /* 170:    1,    3,    5,    7,  */ {5, 2, 1, 2, 5, 6, 0, 3, 4, 7, 4, 3},
    /*  60:       2, 3, 4, 5,        */ {7, 10, 5, 10, 7, 11, 9, 1, 8, 3, 8, 1}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 10.1.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling10_1_2[6][24] = {
    /* 195: 0, 1,             6, 7,  */ {3, 11, 7, 3, 7, 8, 9, 8, 7, 5, 9, 7, 9, 5, 10, 9, 10, 1, 3, 1, 10, 11, 3, 10},
    /*  85: 0,    2,    4,    6,     */ {7, 6, 5, 7, 5, 4, 0, 4, 5, 1, 0, 5, 0, 1, 2, 0, 2, 3, 7, 3, 2, 6, 7, 2},
    /* 105: 0,       3,    5, 6,     */ {11, 2, 10, 6, 11, 10, 11, 6, 4,  11, 4,  8,
                                         0,  8, 4,  9, 0,  4,  0,  9, 10, 0,  10, 2},
    /* 150:    1, 2,    4,       7,  */ {11, 2, 10, 11, 10, 6, 4,  6, 10, 9, 4,  10,
                                         4,  9, 0,  4,  0,  8, 11, 8, 0,  2, 11, 0},
    /* 170:    1,    3,    5,    7,  */ {7, 6, 5, 4, 7, 5, 7, 4, 0, 7, 0, 3, 2, 3, 0, 1, 2, 0, 2, 1, 5, 2, 5, 6},
    /*  60:       2, 3, 4, 5,        */ {7, 8, 3, 11, 7, 3, 7, 11, 10, 7, 10, 5, 9, 5, 10, 1, 9, 10, 9, 1, 3, 9, 3, 8}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 10.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling10_2[6][24] = {
    /* 195: 0, 1,             6, 7,  */ {12, 5, 9,  12, 9,  8,  12, 8,  3, 12, 3, 1,
                                         12, 1, 10, 12, 10, 11, 12, 11, 7, 12, 7, 5},
    /*  85: 0,    2,    4,    6,     */ {12, 1, 0, 12, 0, 4, 12, 4, 7, 12, 7, 3,
                                         12, 3, 2, 12, 2, 6, 12, 6, 5, 12, 5, 1},
    /* 105: 0,       3,    5, 6,     */ {4, 8, 12, 6, 4, 12, 10, 6, 12, 9, 10, 12,
                                         0, 9, 12, 2, 0, 12, 11, 2, 12, 8, 11, 12},
    /* 150:    1, 2,    4,       7,  */ {12, 9, 4, 12, 4, 6, 12, 6, 11, 12, 11, 8,
                                         12, 8, 0, 12, 0, 2, 12, 2, 10, 12, 10, 9},
    /* 170:    1,    3,    5,    7,  */ {0, 3, 12, 4, 0, 12, 5, 4, 12, 1, 5, 12,
                                         2, 1, 12, 6, 2, 12, 7, 6, 12, 3, 7, 12},
    /*  60:       2, 3, 4, 5,        */ {10, 5, 12, 11, 10, 12, 3, 11, 12, 1, 3, 12,
                                         9,  1, 12, 8,  9,  12, 7, 8,  12, 5, 7, 12}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 10.2 inverted
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling10_2_[6][24] = {
    /* 195: 0, 1,             6, 7,  */ {8,  7, 12, 9,  8,  12, 1, 9,  12, 3, 1, 12,
                                         11, 3, 12, 10, 11, 12, 5, 10, 12, 7, 5, 12},
    /*  85: 0,    2,    4,    6,     */ {4, 5, 12, 0, 4, 12, 3, 0, 12, 7, 3, 12,
                                         6, 7, 12, 2, 6, 12, 1, 2, 12, 5, 1, 12},
    /* 105: 0,       3,    5, 6,     */ {12, 11, 6, 12, 6, 4, 12, 4, 9, 12, 9, 10,
                                         12, 10, 2, 12, 2, 0, 12, 0, 8, 12, 8, 11},
    /* 150:    1, 2,    4,       7,  */ {6, 10, 12, 4, 6, 12, 8, 4, 12, 11, 8, 12,
                                         2, 11, 12, 0, 2, 12, 9, 0, 12, 10, 9, 12},
    /* 170:    1,    3,    5,    7,  */ {12, 7, 4, 12, 4, 0, 12, 0, 1, 12, 1, 5,
                                         12, 5, 6, 12, 6, 2, 12, 2, 3, 12, 3, 7},
    /*  60:       2, 3, 4, 5,        */ {12, 7, 11, 12, 11, 10, 12, 10, 1, 12, 1, 3,
                                         12, 3, 8,  12, 8,  9,  12, 9,  5, 12, 5, 7}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief tiling table for case 11
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling11[12][12] = {
    /*  23: 0, 1, 2,    4,           */ {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4},
    /* 139: 0, 1,    3,          7,  */ {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6},
    /*  99: 0, 1,          5, 6,     */ {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10},
    /*  77: 0,    2, 3,       6,     */ {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6},
    /*  57: 0,       3, 4, 5,        */ {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11},
    /* 209: 0,          4,    6, 7,  */ {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0},
    /*  46:    1, 2, 3,    5,        */ {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3},
    /* 198:    1, 2,          6, 7,  */ {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7},
    /* 178:    1,       4, 5,    7,  */ {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11},
    /* 156:       2, 3, 4,       7,  */ {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1},
    /* 116:       2,    4, 5, 6,     */ {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7},
    /* 232:          3,    5, 6, 7,  */ {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief test table for case 12
 * 2 faces to test + eventually the interior
 * When the tests on both specified faces are positive : 4 middle triangles (1)
 * When the test on the first  specified face is positive : 8 first triangles
 * When the test on the second specified face is positive : 8 next triangles
 * When the tests on both specified faces are negative :
 * - if the test on the interior is negative : 4 middle triangles
 * - if the test on the interior is positive : 8 last triangles
 * The support edge for the interior test is marked as the 4th column.
 *
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::test12[24][4] = {
    /* 135: 0, 1, 2,             7,  */ {4, 3, 7, 11},
    /*  75: 0, 1,    3,       6,     */ {3, 2, 7, 10},
    /*  83: 0, 1,       4,    6,     */ {2, 6, 7, 5},
    /* 163: 0, 1,          5,    7,  */ {6, 4, 7, 7},
    /*  45: 0,    2, 3,    5,        */ {2, 1, 7, 9},
    /*  53: 0,    2,    4, 5,        */ {5, 2, 7, 1},
    /* 149: 0,    2,    4,       7,  */ {5, 3, 7, 2},
    /* 101: 0,    2,       5, 6,     */ {5, 1, 7, 0},
    /* 197: 0,    2,          6, 7,  */ {5, 4, 7, 3},
    /*  89: 0,       3, 4,    6,     */ {6, 3, 7, 6},
    /* 169: 0,       3,    5,    7,  */ {1, 6, 7, 4},
    /* 225: 0,             5, 6, 7,  */ {1, 4, 7, 8},
    /*  30:    1, 2, 3, 4,           */ {4, 1, 7, 8},
    /*  86:    1, 2,    4,    6,     */ {6, 1, 7, 4},
    /* 166:    1, 2,       5,    7,  */ {3, 6, 7, 6},
    /*  58:    1,    3, 4, 5,        */ {4, 5, 7, 3},
    /* 154:    1,    3, 4,       7,  */ {1, 5, 7, 0},
    /* 106:    1,    3,    5, 6,     */ {3, 5, 7, 2},
    /* 202:    1,    3,       6, 7,  */ {2, 5, 7, 1},
    /* 210:    1,       4,    6, 7,  */ {1, 2, 7, 9},
    /*  92:       2, 3, 4,    6,     */ {4, 6, 7, 7},
    /* 172:       2, 3,    5,    7,  */ {6, 2, 7, 5},
    /* 180:       2,    4, 5,    7,  */ {2, 3, 7, 10},
    /* 120:          3, 4, 5, 6,     */ {3, 4, 7, 11}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 12.1.1
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling12_1_1[24][12] = {
    /* 135: 0, 1, 2,             7,  */ {7, 6, 11, 10, 3, 2, 3, 10, 8, 9, 8, 10},
    /*  75: 0, 1,    3,       6,     */ {6, 5, 10, 9, 2, 1, 2, 9, 11, 8, 11, 9},
    /*  83: 0, 1,       4,    6,     */ {10, 6, 5, 7, 9, 4, 9, 7, 1, 3, 1, 7},
    /* 163: 0, 1,          5,    7,  */ {7, 6, 11, 4, 8, 5, 3, 5, 8, 5, 3, 1},
    /*  45: 0,    2, 3,    5,        */ {5, 4, 9, 8, 1, 0, 1, 8, 10, 11, 10, 8},
    /*  53: 0,    2,    4, 5,        */ {1, 2, 10, 0, 9, 3, 5, 3, 9, 3, 5, 7},
    /* 149: 0,    2,    4,       7,  */ {10, 1, 2, 0, 11, 3, 11, 0, 6, 4, 6, 0},
    /* 101: 0,    2,       5, 6,     */ {8, 3, 0, 2, 9, 1, 9, 2, 4, 6, 4, 2},
    /* 197: 0,    2,          6, 7,  */ {3, 0, 8, 2, 11, 1, 7, 1, 11, 1, 7, 5},
    /*  89: 0,       3, 4,    6,     */ {6, 5, 10, 7, 11, 4, 2, 4, 11, 4, 2, 0},
    /* 169: 0,       3,    5,    7,  */ {9, 5, 4, 6, 8, 7, 8, 6, 0, 2, 0, 6},
    /* 225: 0,             5, 6, 7,  */ {8, 3, 0, 7, 4, 11, 9, 11, 4, 11, 9, 10},
    /*  30:    1, 2, 3, 4,           */ {4, 7, 8, 11, 0, 3, 0, 11, 9, 10, 9, 11},
    /*  86:    1, 2,    4,    6,     */ {4, 7, 8, 5, 9, 6, 0, 6, 9, 6, 0, 2},
    /* 166:    1, 2,       5,    7,  */ {11, 7, 6, 4, 10, 5, 10, 4, 2, 0, 2, 4},
    /*  58:    1,    3, 4, 5,        */ {11, 2, 3, 1, 8, 0, 8, 1, 7, 5, 7, 1},
    /* 154:    1,    3, 4,       7,  */ {0, 1, 9, 3, 8, 2, 4, 2, 8, 2, 4, 6},
    /* 106:    1,    3,    5, 6,     */ {2, 3, 11, 1, 10, 0, 6, 0, 10, 0, 6, 4},
    /* 202:    1,    3,       6, 7,  */ {9, 0, 1, 3, 10, 2, 10, 3, 5, 7, 5, 3},
    /* 210:    1,       4,    6, 7,  */ {9, 0, 1, 4, 5, 8, 10, 8, 5, 8, 10, 11},
    /*  92:       2, 3, 4,    6,     */ {8, 4, 7, 5, 11, 6, 11, 5, 3, 1, 3, 5},
    /* 172:       2, 3,    5,    7,  */ {5, 4, 9, 6, 10, 7, 1, 7, 10, 7, 1, 3},
    /* 180:       2,    4, 5,    7,  */ {10, 1, 2, 5, 6, 9, 11, 9, 6, 9, 11, 8},
    /* 120:          3, 4, 5, 6,     */ {11, 2, 3, 6, 7, 10, 8, 10, 7, 10, 8, 9}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 12.1.1 inverted
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling12_1_1_[24][12] = {
    /* 135: 0, 1, 2,             7,  */ {3, 2, 11, 10, 7, 6, 7, 10, 8, 9, 8, 10},
    /*  75: 0, 1,    3,       6,     */ {2, 1, 10, 9, 6, 5, 6, 9, 11, 8, 11, 9},
    /*  83: 0, 1,       4,    6,     */ {9, 4, 5, 7, 10, 6, 10, 7, 1, 3, 1, 7},
    /* 163: 0, 1,          5,    7,  */ {7, 4, 8, 6, 11, 5, 3, 5, 11, 5, 3, 1},
    /*  45: 0,    2, 3,    5,        */ {1, 0, 9, 8, 5, 4, 5, 8, 10, 11, 10, 8},
    /*  53: 0,    2,    4, 5,        */ {1, 0, 9, 2, 10, 3, 5, 3, 10, 3, 5, 7},
    /* 149: 0,    2,    4,       7,  */ {11, 3, 2, 0, 10, 1, 10, 0, 6, 4, 6, 0},
    /* 101: 0,    2,       5, 6,     */ {9, 1, 0, 2, 8, 3, 8, 2, 4, 6, 4, 2},
    /* 197: 0,    2,          6, 7,  */ {3, 2, 11, 0, 8, 1, 7, 1, 8, 1, 7, 5},
    /*  89: 0,       3, 4,    6,     */ {6, 7, 11, 5, 10, 4, 2, 4, 10, 4, 2, 0},
    /* 169: 0,       3,    5,    7,  */ {8, 7, 4, 6, 9, 5, 9, 6, 0, 2, 0, 6},
    /* 225: 0,             5, 6, 7,  */ {8, 7, 4, 3, 0, 11, 9, 11, 0, 11, 9, 10},
    /*  30:    1, 2, 3, 4,           */ {0, 3, 8, 11, 4, 7, 4, 11, 9, 10, 9, 11},
    /*  86:    1, 2,    4,    6,     */ {4, 5, 9, 7, 8, 6, 0, 6, 8, 6, 0, 2},
    /* 166:    1, 2,       5,    7,  */ {10, 5, 6, 4, 11, 7, 11, 4, 2, 0, 2, 4},
    /*  58:    1,    3, 4, 5,        */ {8, 0, 3, 1, 11, 2, 11, 1, 7, 5, 7, 1},
    /* 154:    1,    3, 4,       7,  */ {0, 3, 8, 1, 9, 2, 4, 2, 9, 2, 4, 6},
    /* 106:    1,    3,    5, 6,     */ {2, 1, 10, 3, 11, 0, 6, 0, 11, 0, 6, 4},
    /* 202:    1,    3,       6, 7,  */ {10, 2, 1, 3, 9, 0, 9, 3, 5, 7, 5, 3},
    /* 210:    1,       4,    6, 7,  */ {9, 4, 5, 0, 1, 8, 10, 8, 1, 8, 10, 11},
    /*  92:       2, 3, 4,    6,     */ {11, 6, 7, 5, 8, 4, 8, 5, 3, 1, 3, 5},
    /* 172:       2, 3,    5,    7,  */ {5, 6, 10, 4, 9, 7, 1, 7, 9, 7, 1, 3},
    /* 180:       2,    4, 5,    7,  */ {10, 5, 6, 1, 2, 9, 11, 9, 2, 9, 11, 8},
    /* 120:          3, 4, 5, 6,     */ {11, 6, 7, 2, 3, 10, 8, 10, 3, 10, 8, 9}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 12.1.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling12_1_2[24][24] = {
    /* 135: 0, 1, 2,             7,  */ {7, 3, 11, 3, 7, 8, 9, 8, 7, 6, 9, 7, 9, 6, 10, 2, 10, 6, 11, 2, 6, 2, 11, 3},
    /*  75: 0, 1,    3,       6,     */ {6, 2, 10, 2, 6, 11, 8, 11, 6, 5, 8, 6, 8, 5, 9, 1, 9, 5, 10, 1, 5, 1, 10, 2},
    /*  83: 0, 1,       4,    6,     */ {10, 9, 5, 9, 10, 1, 3, 1, 10, 6, 3, 10, 3, 6, 7, 4, 7, 6, 5, 4, 6, 4, 5, 9},
    /* 163: 0, 1,          5,    7,  */ {7, 8, 11, 3, 11, 8, 11, 3, 1, 11, 1, 6, 5, 6, 1, 6, 5, 4, 6, 4, 7, 8, 7, 4},
    /*  45: 0,    2, 3,    5,        */ {5, 1, 9, 1, 5, 10, 11, 10, 5, 4, 11, 5, 11, 4, 8, 0, 8, 4, 9, 0, 4, 0, 9, 1},
    /*  53: 0,    2,    4, 5,        */ {1, 9, 10, 5, 10, 9, 10, 5, 7, 10, 7, 2, 3, 2, 7, 2, 3, 0, 2, 0, 1, 9, 1, 0},
    /* 149: 0,    2,    4,       7,  */ {10, 11, 2, 11, 10, 6, 4, 6, 10, 1, 4, 10, 4, 1, 0, 3, 0, 1, 2, 3, 1, 3, 2, 11},
    /* 101: 0,    2,       5, 6,     */ {8, 9, 0, 9, 8, 4, 6, 4, 8, 3, 6, 8, 6, 3, 2, 1, 2, 3, 0, 1, 3, 1, 0, 9},
    /* 197: 0,    2,          6, 7,  */ {3, 11, 8, 7, 8, 11, 8, 7, 5, 8, 5, 0, 1, 0, 5, 0, 1, 2, 0, 2, 3, 11, 3, 2},
    /*  89: 0,       3, 4,    6,     */ {6, 11, 10, 2, 10, 11, 10, 2, 0, 10, 0, 5, 4, 5, 0, 5, 4, 7, 5, 7, 6, 11, 6, 7},
    /* 169: 0,       3,    5,    7,  */ {9, 8, 4, 8, 9, 0, 2, 0, 9, 5, 2, 9, 2, 5, 6, 7, 6, 5, 4, 7, 5, 7, 4, 8},
    /* 225: 0,             5, 6, 7,  */ {8, 4, 0, 9, 0, 4, 0, 9, 10, 0, 10, 3, 11, 3, 10, 3, 11, 7, 3, 7, 8, 4, 8, 7},
    /*  30:    1, 2, 3, 4,           */ {4, 0, 8, 0, 4, 9, 10, 9, 4, 7, 10, 4, 10, 7, 11, 3, 11, 7, 8, 3, 7, 3, 8, 0},
    /*  86:    1, 2,    4,    6,     */ {4, 9, 8, 0, 8, 9, 8, 0, 2, 8, 2, 7, 6, 7, 2, 7, 6, 5, 7, 5, 4, 9, 4, 5},
    /* 166:    1, 2,       5,    7,  */ {11, 10, 6, 10, 11, 2, 0, 2, 11, 7, 0, 11, 0, 7, 4, 5, 4, 7, 6, 5, 7, 5, 6, 10},
    /*  58:    1,    3, 4, 5,        */ {11, 8, 3, 8, 11, 7, 5, 7, 11, 2, 5, 11, 5, 2, 1, 0, 1, 2, 3, 0, 2, 0, 3, 8},
    /* 154:    1,    3, 4,       7,  */ {0, 8, 9, 4, 9, 8, 9, 4, 6, 9, 6, 1, 2, 1, 6, 1, 2, 3, 1, 3, 0, 8, 0, 3},
    /* 106:    1,    3,    5, 6,     */ {2, 10, 11, 6, 11, 10, 11, 6, 4, 11, 4, 3, 0, 3, 4, 3, 0, 1, 3, 1, 2, 10, 2, 1},
    /* 202:    1,    3,       6, 7,  */ {9, 10, 1, 10, 9, 5, 7, 5, 9, 0, 7, 9, 7, 0, 3, 2, 3, 0, 1, 2, 0, 2, 1, 10},
    /* 210:    1,       4,    6, 7,  */ {9, 5, 1, 10, 1, 5, 1, 10, 11, 1, 11, 0, 8, 0, 11, 0, 8, 4, 0, 4, 9, 5, 9, 4},
    /*  92:       2, 3, 4,    6,     */ {8, 11, 7, 11, 8, 3, 1, 3, 8, 4, 1, 8, 1, 4, 5, 6, 5, 4, 7, 6, 4, 6, 7, 11},
    /* 172:       2, 3,    5,    7,  */ {5, 10, 9, 1, 9, 10, 9, 1, 3, 9, 3, 4, 7, 4, 3, 4, 7, 6, 4, 6, 5, 10, 5, 6},
    /* 180:       2,    4, 5,    7,  */ {10, 6, 2, 11, 2, 6, 2, 11, 8, 2, 8, 1, 9, 1, 8, 1, 9, 5, 1, 5, 10, 6, 10, 5},
    /* 120:          3, 4, 5, 6,     */ {11, 7, 3, 8, 3, 7, 3, 8, 9, 3, 9, 2, 10, 2, 9, 2, 10, 6, 2, 6, 11, 7, 11, 6}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 12.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling12_2[24][24] = {
    /* 135: 0, 1, 2,             7,  */ {9,  8, 12, 10, 9,  12, 2, 10, 12, 3, 2, 12,
                                         11, 3, 12, 6,  11, 12, 7, 6,  12, 8, 7, 12},
    /*  75: 0, 1,    3,       6,     */ {8,  11, 12, 9, 8,  12, 1, 9, 12, 2,  1, 12,
                                         10, 2,  12, 5, 10, 12, 6, 5, 12, 11, 6, 12},
    /*  83: 0, 1,       4,    6,     */ {3, 1, 12, 7, 3, 12, 4,  7, 12, 9, 4,  12,
                                         5, 9, 12, 6, 5, 12, 10, 6, 12, 1, 10, 12},
    /* 163: 0, 1,          5,    7,  */ {12, 3,  1, 12, 1, 5, 12, 5, 6, 12, 6, 11,
                                         12, 11, 7, 12, 7, 4, 12, 4, 8, 12, 8, 3},
    /*  45: 0,    2, 3,    5,        */ {11, 10, 12, 8, 11, 12, 0, 8, 12, 1,  0, 12,
                                         9,  1,  12, 4, 9,  12, 5, 4, 12, 10, 5, 12},
    /*  53: 0,    2,    4, 5,        */ {12, 5,  7, 12, 7, 3, 12, 3, 2, 12, 2, 10,
                                         12, 10, 1, 12, 1, 0, 12, 0, 9, 12, 9, 5},
    /* 149: 0,    2,    4,       7,  */ {4, 6,  12, 0, 4, 12, 1,  0, 12, 10, 1,  12,
                                         2, 10, 12, 3, 2, 12, 11, 3, 12, 6,  11, 12},
    /* 101: 0,    2,       5, 6,     */ {6, 4, 12, 2, 6, 12, 3, 2, 12, 8, 3, 12,
                                         0, 8, 12, 1, 0, 12, 9, 1, 12, 4, 9, 12},
    /* 197: 0,    2,          6, 7,  */ {12, 7, 5, 12, 5, 1, 12, 1, 0,  12, 0,  8,
                                         12, 8, 3, 12, 3, 2, 12, 2, 11, 12, 11, 7},
    /*  89: 0,       3, 4,    6,     */ {12, 2,  0, 12, 0, 4, 12, 4, 5,  12, 5,  10,
                                         12, 10, 6, 12, 6, 7, 12, 7, 11, 12, 11, 2},
    /* 169: 0,       3,    5,    7,  */ {2, 0, 12, 6, 2, 12, 7, 6, 12, 8, 7, 12,
                                         4, 8, 12, 5, 4, 12, 9, 5, 12, 0, 9, 12},
    /* 225: 0,             5, 6, 7,  */ {12, 9, 10, 12, 10, 11, 12, 11, 7, 12, 7, 4,
                                         12, 4, 8,  12, 8,  3,  12, 3,  0, 12, 0, 9},
    /*  30:    1, 2, 3, 4,           */ {10, 9, 12, 11, 10, 12, 7, 11, 12, 4, 7, 12,
                                         8,  4, 12, 3,  8,  12, 0, 3,  12, 9, 0, 12},
    /*  86:    1, 2,    4,    6,     */ {12, 0, 2, 12, 2, 6, 12, 6, 7, 12, 7, 8,
                                         12, 8, 4, 12, 4, 5, 12, 5, 9, 12, 9, 0},
    /* 166:    1, 2,       5,    7,  */ {0, 2,  12, 4, 0, 12, 5,  4, 12, 10, 5,  12,
                                         6, 10, 12, 7, 6, 12, 11, 7, 12, 2,  11, 12},
    /*  58:    1,    3, 4, 5,        */ {5, 7, 12, 1, 5, 12, 0,  1, 12, 8, 0,  12,
                                         3, 8, 12, 2, 3, 12, 11, 2, 12, 7, 11, 12},
    /* 154:    1,    3, 4,       7,  */ {12, 4, 6, 12, 6, 2, 12, 2, 3, 12, 3, 8,
                                         12, 8, 0, 12, 0, 1, 12, 1, 9, 12, 9, 4},
    /* 106:    1,    3,    5, 6,     */ {12, 6,  4, 12, 4, 0, 12, 0, 1,  12, 1,  10,
                                         12, 10, 2, 12, 2, 3, 12, 3, 11, 12, 11, 6},
    /* 202:    1,    3,       6, 7,  */ {7, 5,  12, 3, 7, 12, 2, 3, 12, 10, 2, 12,
                                         1, 10, 12, 0, 1, 12, 9, 0, 12, 5,  9, 12},
    /* 210:    1,       4,    6, 7,  */ {12, 10, 11, 12, 11, 8, 12, 8, 0, 12, 0, 1,
                                         12, 1,  9,  12, 9,  4, 12, 4, 5, 12, 5, 10},
    /*  92:       2, 3, 4,    6,     */ {1, 3,  12, 5, 1, 12, 6, 5, 12, 11, 6, 12,
                                         7, 11, 12, 4, 7, 12, 8, 4, 12, 3,  8, 12},
    /* 172:       2, 3,    5,    7,  */ {12, 1, 3, 12, 3, 7, 12, 7, 4,  12, 4,  9,
                                         12, 9, 5, 12, 5, 6, 12, 6, 10, 12, 10, 1},
    /* 180:       2,    4, 5,    7,  */ {12, 11, 8,  12, 8,  9, 12, 9, 1, 12, 1, 2,
                                         12, 2,  10, 12, 10, 5, 12, 5, 6, 12, 6, 11},
    /* 120:          3, 4, 5, 6,     */ {12, 8, 9,  12, 9,  10, 12, 10, 2, 12, 2, 3,
                                         12, 3, 11, 12, 11, 6,  12, 6,  7, 12, 7, 8}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 12.2 inverted
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling12_2_[24][24] = {
    /* 135: 0, 1, 2,             7,  */ {12, 2,  11, 12, 11, 7, 12, 7, 6, 12, 6, 10,
                                         12, 10, 9,  12, 9,  8, 12, 8, 3, 12, 3, 2},
    /*  75: 0, 1,    3,       6,     */ {12, 1, 10, 12, 10, 6,  12, 6,  5, 12, 5, 9,
                                         12, 9, 8,  12, 8,  11, 12, 11, 2, 12, 2, 1},
    /*  83: 0, 1,       4,    6,     */ {12, 4, 5, 12, 5, 10, 12, 10, 6, 12, 6, 7,
                                         12, 7, 3, 12, 3, 1,  12, 1,  9, 12, 9, 4},
    /* 163: 0, 1,          5,    7,  */ {7, 6, 12, 8, 7, 12, 4,  8, 12, 5, 4,  12,
                                         1, 5, 12, 3, 1, 12, 11, 3, 12, 6, 11, 12},
    /*  45: 0,    2, 3,    5,        */ {12, 0, 9,  12, 9,  5,  12, 5,  4, 12, 4, 8,
                                         12, 8, 11, 12, 11, 10, 12, 10, 1, 12, 1, 0},
    /*  53: 0,    2,    4, 5,        */ {1, 2, 12, 9, 1, 12, 0,  9, 12, 3, 0,  12,
                                         7, 3, 12, 5, 7, 12, 10, 5, 12, 2, 10, 12},
    /* 149: 0,    2,    4,       7,  */ {12, 1, 2, 12, 2, 11, 12, 11, 3,  12, 3,  0,
                                         12, 0, 4, 12, 4, 6,  12, 6,  10, 12, 10, 1},
    /* 101: 0,    2,       5, 6,     */ {12, 3, 0, 12, 0, 9, 12, 9, 1, 12, 1, 2,
                                         12, 2, 6, 12, 6, 4, 12, 4, 8, 12, 8, 3},
    /* 197: 0,    2,          6, 7,  */ {3, 0, 12, 11, 3, 12, 2, 11, 12, 1, 2, 12,
                                         5, 1, 12, 7,  5, 12, 8, 7,  12, 0, 8, 12},
    /*  89: 0,       3, 4,    6,     */ {6, 5, 12, 11, 6, 12, 7,  11, 12, 4, 7,  12,
                                         0, 4, 12, 2,  0, 12, 10, 2,  12, 5, 10, 12},
    /* 169: 0,       3,    5,    7,  */ {12, 7, 4, 12, 4, 9, 12, 9, 5, 12, 5, 6,
                                         12, 6, 2, 12, 2, 0, 12, 0, 8, 12, 8, 7},
    /* 225: 0,             5, 6, 7,  */ {8,  7,  12, 0, 8,  12, 3, 0, 12, 11, 3, 12,
                                         10, 11, 12, 9, 10, 12, 4, 9, 12, 7,  4, 12},
    /*  30:    1, 2, 3, 4,           */ {12, 7,  8,  12, 8,  0, 12, 0, 3, 12, 3, 11,
                                         12, 11, 10, 12, 10, 9, 12, 9, 4, 12, 4, 7},
    /*  86:    1, 2,    4,    6,     */ {4, 7, 12, 9, 4, 12, 5, 9, 12, 6, 5, 12,
                                         2, 6, 12, 0, 2, 12, 8, 0, 12, 7, 8, 12},
    /* 166:    1, 2,       5,    7,  */ {12, 5, 6, 12, 6, 11, 12, 11, 7,  12, 7,  4,
                                         12, 4, 0, 12, 0, 2,  12, 2,  10, 12, 10, 5},
    /*  58:    1,    3, 4, 5,        */ {12, 0, 3, 12, 3, 11, 12, 11, 2, 12, 2, 1,
                                         12, 1, 5, 12, 5, 7,  12, 7,  8, 12, 8, 0},
    /* 154:    1,    3, 4,       7,  */ {0, 3, 12, 9, 0, 12, 1, 9, 12, 2, 1, 12,
                                         6, 2, 12, 4, 6, 12, 8, 4, 12, 3, 8, 12},
    /* 106:    1,    3,    5, 6,     */ {2, 1, 12, 11, 2, 12, 3,  11, 12, 0, 3,  12,
                                         4, 0, 12, 6,  4, 12, 10, 6,  12, 1, 10, 12},
    /* 202:    1,    3,       6, 7,  */ {12, 2, 1, 12, 1, 9, 12, 9, 0,  12, 0,  3,
                                         12, 3, 7, 12, 7, 5, 12, 5, 10, 12, 10, 2},
    /* 210:    1,       4,    6, 7,  */ {9,  0, 12, 5,  9,  12, 4, 5,  12, 8, 4, 12,
                                         11, 8, 12, 10, 11, 12, 1, 10, 12, 0, 1, 12},
    /*  92:       2, 3, 4,    6,     */ {12, 6, 7, 12, 7, 8, 12, 8, 4,  12, 4,  5,
                                         12, 5, 1, 12, 1, 3, 12, 3, 11, 12, 11, 6},
    /* 172:       2, 3,    5,    7,  */ {5, 4, 12, 10, 5, 12, 6, 10, 12, 7, 6, 12,
                                         3, 7, 12, 1,  3, 12, 9, 1,  12, 4, 9, 12},
    /* 180:       2,    4, 5,    7,  */ {10, 1, 12, 6,  10, 12, 5, 6,  12, 9, 5, 12,
                                         8,  9, 12, 11, 8,  12, 2, 11, 12, 1, 2, 12},
    /* 120:          3, 4, 5, 6,     */ {11, 2,  12, 7, 11, 12, 6, 7, 12, 10, 6, 12,
                                         9,  10, 12, 8, 9,  12, 3, 8, 12, 2,  3, 12}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief test table for case 13
 * All faces are to be tested
 *
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
/* 13: face test */
const schar MarchingCubes::test13[2][7] = {
    /* 165: 0,    2,       5,    7,  */ {1, 2, 3, 4, 5, 6, 7},
    /*  90:    1,    3, 4,    6,     */ {2, 3, 4, 1, 5, 6, 7},
};

//_____________________________________________________________________________
/**
 * \brief subconfiguration table for case 13
 * Hard-coded tests for the subconfiguration determination
 *
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
/* 13: sub configs */
const schar MarchingCubes::subconfig13[64] = {
    /*  0: 0,0,0,0,0,0 */ 0,
    /*  1: 1,0,0,0,0,0 */ 1,
    /*  2: 0,1,0,0,0,0 */ 2,
    /*  3: 1,1,0,0,0,0 */ 7,
    /*  4: 0,0,1,0,0,0 */ 3,
    /*  5: 1,0,1,0,0,0 */ -1,
    /*  6: 0,1,1,0,0,0 */ 11,
    /*  7: 1,1,1,0,0,0 */ -1,
    /*  8: 0,0,0,1,0,0 */ 4,
    /*  9: 1,0,0,1,0,0 */ 8,
    /* 10: 0,1,0,1,0,0 */ -1,
    /* 11: 1,1,0,1,0,0 */ -1,
    /* 12: 0,0,1,1,0,0 */ 14,
    /* 13: 1,0,1,1,0,0 */ -1,
    /* 14: 0,1,1,1,0,0 */ -1,
    /* 15: 1,1,1,1,0,0 */ -1,
    /* 16: 0,0,0,0,1,0 */ 5,
    /* 17: 1,0,0,0,1,0 */ 9,
    /* 18: 0,1,0,0,1,0 */ 12,
    /* 19: 1,1,0,0,1,0 */ 23,
    /* 20: 0,0,1,0,1,0 */ 15,
    /* 21: 1,0,1,0,1,0 */ -1,
    /* 22: 0,1,1,0,1,0 */ 21,
    /* 23: 1,1,1,0,1,0 */ 38,
    /* 24: 0,0,0,1,1,0 */ 17,
    /* 25: 1,0,0,1,1,0 */ 20,
    /* 26: 0,1,0,1,1,0 */ -1,
    /* 27: 1,1,0,1,1,0 */ 36,
    /* 28: 0,0,1,1,1,0 */ 26,
    /* 29: 1,0,1,1,1,0 */ 33,
    /* 30: 0,1,1,1,1,0 */ 30,
    /* 31: 1,1,1,1,1,0 */ 44,
    /* 32: 0,0,0,0,0,1 */ 6,
    /* 33: 1,0,0,0,0,1 */ 10,
    /* 34: 0,1,0,0,0,1 */ 13,
    /* 35: 1,1,0,0,0,1 */ 19,
    /* 36: 0,0,1,0,0,1 */ 16,
    /* 37: 1,0,1,0,0,1 */ -1,
    /* 38: 0,1,1,0,0,1 */ 25,
    /* 39: 1,1,1,0,0,1 */ 37,
    /* 40: 0,0,0,1,0,1 */ 18,
    /* 41: 1,0,0,1,0,1 */ 24,
    /* 42: 0,1,0,1,0,1 */ -1,
    /* 43: 1,1,0,1,0,1 */ 35,
    /* 44: 0,0,1,1,0,1 */ 22,
    /* 45: 1,0,1,1,0,1 */ 32,
    /* 46: 0,1,1,1,0,1 */ 29,
    /* 47: 1,1,1,1,0,1 */ 43,
    /* 48: 0,0,0,0,1,1 */ -1,
    /* 49: 1,0,0,0,1,1 */ -1,
    /* 50: 0,1,0,0,1,1 */ -1,
    /* 51: 1,1,0,0,1,1 */ 34,
    /* 52: 0,0,1,0,1,1 */ -1,
    /* 53: 1,0,1,0,1,1 */ -1,
    /* 54: 0,1,1,0,1,1 */ 28,
    /* 55: 1,1,1,0,1,1 */ 42,
    /* 56: 0,0,0,1,1,1 */ -1,
    /* 57: 1,0,0,1,1,1 */ 31,
    /* 58: 0,1,0,1,1,1 */ -1,
    /* 59: 1,1,0,1,1,1 */ 41,
    /* 60: 0,0,1,1,1,1 */ 27,
    /* 61: 1,0,1,1,1,1 */ 40,
    /* 62: 0,1,1,1,1,1 */ 39,
    /* 63: 1,1,1,1,1,1 */ 45,
};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 13.1
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
/* 13.1 */
const schar MarchingCubes::tiling13_1[2][12] = {
    /* 165: 0,    2,       5,    7,  */ {11, 7, 6, 1, 2, 10, 8, 3, 0, 9, 5, 4},
    /*  90:    1,    3, 4,    6,     */ {8, 4, 7, 2, 3, 11, 9, 0, 1, 10, 6, 5}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 13.1 inverted
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
/* 13.1 */
const schar MarchingCubes::tiling13_1_[2][12] = {
    /* 165: 0,    2,       5,    7,  */ {7, 4, 8, 11, 3, 2, 1, 0, 9, 5, 6, 10},
    /*  90:    1,    3, 4,    6,     */ {6, 7, 11, 10, 2, 1, 0, 3, 8, 4, 5, 9}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 13.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
/* 13.2 */
const schar MarchingCubes::tiling13_2[2][6][18] = {
    /* 165: 0,    2,       5,    7,  */ {/* 1 */ {1, 2, 10, 11, 7, 6, 3, 4, 8, 4, 3, 5, 0, 5, 3, 5, 0, 9},
                                         /* 2 */ {8, 3, 0, 11, 7, 6, 9, 1, 4, 2, 4, 1, 4, 2, 5, 10, 5, 2},
                                         /* 3 */ {9, 5, 4, 8, 3, 0, 1, 6, 10, 6, 1, 7, 2, 7, 1, 7, 2, 11},
                                         /* 4 */ {9, 5, 4, 1, 2, 10, 11, 3, 6, 0, 6, 3, 6, 0, 7, 8, 7, 0},
                                         /* 5 */ {9, 5, 4, 11, 7, 6, 0, 10, 1, 10, 0, 8, 10, 8, 2, 3, 2, 8},
                                         /* 6 */ {1, 2, 10, 3, 0, 8, 4, 9, 7, 11, 7, 9, 5, 11, 9, 11, 5, 6}},
    /*  90:    1,    3, 4,    6,     */ {/* 1 */ {2, 3, 11, 8, 4, 7, 0, 5, 9, 5, 0, 6, 1, 6, 0, 6, 1, 10},
                                         /* 2 */ {9, 0, 1, 8, 4, 7, 10, 2, 5, 3, 5, 2, 5, 3, 6, 11, 6, 3},
                                         /* 3 */ {6, 5, 10, 9, 0, 1, 2, 7, 11, 7, 2, 4, 3, 4, 2, 4, 3, 8},
                                         /* 4 */ {6, 5, 10, 2, 3, 11, 8, 0, 7, 1, 7, 0, 7, 1, 4, 9, 4, 1},
                                         /* 5 */ {6, 5, 10, 8, 4, 7, 1, 11, 2, 11, 1, 9, 11, 9, 3, 0, 3, 9},
                                         /* 6 */ {2, 3, 11, 0, 1, 9, 5, 10, 4, 8, 4, 10, 6, 8, 10, 8, 6, 7}}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 13.2 inverted
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
/* 13.2 */
const schar MarchingCubes::tiling13_2_[2][6][18] = {
    /* 165: 0,    2,       5,    7,  */ {/* 1 */ {10, 5, 6, 11, 3, 2, 7, 0, 8, 0, 7, 1, 4, 1, 7, 1, 4, 9},
                                         /* 2 */ {11, 3, 2, 7, 4, 8, 9, 5, 0, 6, 0, 5, 0, 6, 1, 10, 1, 6},
                                         /* 3 */ {1, 0, 9, 7, 4, 8, 5, 2, 10, 2, 5, 3, 6, 3, 5, 3, 6, 11},
                                         /* 4 */ {10, 5, 6, 1, 0, 9, 11, 7, 2, 4, 2, 7, 2, 4, 3, 8, 3, 4},
                                         /* 5 */ {10, 5, 6, 7, 4, 8, 2, 11, 1, 9, 1, 11, 3, 9, 11, 9, 3, 0},
                                         /* 6 */ {11, 3, 2, 9, 1, 0, 4, 10, 5, 10, 4, 8, 10, 8, 6, 7, 6, 8}},
    /*  90:    1,    3, 4,    6,     */ {/* 1 */ {6, 7, 11, 8, 0, 3, 4, 1, 9, 1, 4, 2, 5, 2, 4, 2, 5, 10},
                                         /* 2 */ {8, 0, 3, 4, 5, 9, 10, 6, 1, 7, 1, 6, 1, 7, 2, 11, 2, 7},
                                         /* 3 */ {2, 1, 10, 4, 5, 9, 6, 3, 11, 3, 6, 0, 7, 0, 6, 0, 7, 8},
                                         /* 4 */ {6, 7, 11, 2, 1, 10, 8, 4, 3, 5, 3, 4, 3, 5, 0, 9, 0, 5},
                                         /* 5 */ {6, 7, 11, 4, 5, 9, 3, 8, 2, 10, 2, 8, 0, 10, 8, 10, 0, 1},
                                         /* 6 */ {8, 0, 3, 10, 2, 1, 5, 11, 6, 11, 5, 9, 11, 9, 7, 4, 7, 9}}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 13.3
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
/* 13.3 */
const schar MarchingCubes::tiling13_3[2][12][30] = {
    /* 165: 0,    2,       5,    7,  */
    {/* 1,2 */ {11, 7, 6, 12, 2, 10, 12, 10, 5, 12, 5, 4, 12, 4, 8, 12, 8, 3, 12, 3, 0, 12, 0, 9, 12, 9, 1, 12, 1, 2},
     /* 1,4 */ {1, 2, 10, 9, 5, 12, 0, 9, 12, 3, 0, 12, 11, 3, 12, 6, 11, 12, 7, 6, 12, 8, 7, 12, 4, 8, 12, 5, 4, 12},
     /* 1,5 */ {11, 7, 6, 12, 5, 4, 12, 4, 8, 12, 8, 3, 12, 3, 2, 12, 2, 10, 12, 10, 1, 12, 1, 0, 12, 0, 9, 12, 9, 5},
     /* 1,6 */ {1, 2, 10, 12, 3, 0, 12, 0, 9, 12, 9, 5, 12, 5, 6, 12, 6, 11, 12, 11, 7, 12, 7, 4, 12, 4, 8, 12, 8, 3},
     /* 2,3 */ {8, 3, 0, 11, 7, 12, 2, 11, 12, 1, 2, 12, 9, 1, 12, 4, 9, 12, 5, 4, 12, 10, 5, 12, 6, 10, 12, 7, 6, 12},
     /* 2,5 */ {11, 7, 6, 5, 4, 12, 10, 5, 12, 2, 10, 12, 3, 2, 12, 8, 3, 12, 0, 8, 12, 1, 0, 12, 9, 1, 12, 4, 9, 12},
     /* 2,6 */ {8, 3, 0, 1, 2, 12, 9, 1, 12, 4, 9, 12, 7, 4, 12, 11, 7, 12, 6, 11, 12, 5, 6, 12, 10, 5, 12, 2, 10, 12},
     /* 3,4 */ {9, 5, 4, 12, 0, 8, 12, 8, 7, 12, 7, 6, 12, 6, 10, 12, 10, 1, 12, 1, 2, 12, 2, 11, 12, 11, 3, 12, 3, 0},
     /* 3,5 */ {9, 5, 4, 12, 7, 6, 12, 6, 10, 12, 10, 1, 12, 1, 0, 12, 0, 8, 12, 8, 3, 12, 3, 2, 12, 2, 11, 12, 11, 7},
     /* 3,6 */ {8, 3, 0, 12, 1, 2, 12, 2, 11, 12, 11, 7, 12, 7, 4, 12, 4, 9, 12, 9, 5, 12, 5, 6, 12, 6, 10, 12, 10, 1},
     /* 4,5 */ {9, 5, 4, 7, 6, 12, 8, 7, 12, 0, 8, 12, 1, 0, 12, 10, 1, 12, 2, 10, 12, 3, 2, 12, 11, 3, 12, 6, 11, 12},
     /* 4,6 */ {1, 2, 10, 3, 0, 12, 11, 3, 12, 6, 11, 12, 5, 6, 12, 9, 5, 12, 4, 9, 12, 7, 4, 12, 8, 7, 12, 0, 8, 12}},
    /*  90:    1,    3, 4,    6,     */
    {/* 1,2 */ {8, 4, 7, 12, 3, 11, 12, 11, 6, 12, 6, 5, 12, 5, 9, 12, 9, 0, 12, 0, 1, 12, 1, 10, 12, 10, 2, 12, 2, 3},
     /* 1,4 */ {2, 3, 11, 10, 6, 12, 1, 10, 12, 0, 1, 12, 8, 0, 12, 7, 8, 12, 4, 7, 12, 9, 4, 12, 5, 9, 12, 6, 5, 12},
     /* 1,5 */ {8, 4, 7, 12, 6, 5, 12, 5, 9, 12, 9, 0, 12, 0, 3, 12, 3, 11, 12, 11, 2, 12, 2, 1, 12, 1, 10, 12, 10, 6},
     /* 1,6 */ {2, 3, 11, 12, 0, 1, 12, 1, 10, 12, 10, 6, 12, 6, 7, 12, 7, 8, 12, 8, 4, 12, 4, 5, 12, 5, 9, 12, 9, 0},
     /* 2,3 */ {0, 1, 9, 8, 4, 12, 3, 8, 12, 2, 3, 12, 10, 2, 12, 5, 10, 12, 6, 5, 12, 11, 6, 12, 7, 11, 12, 4, 7, 12},
     /* 2,5 */ {8, 4, 7, 6, 5, 12, 11, 6, 12, 3, 11, 12, 0, 3, 12, 9, 0, 12, 1, 9, 12, 2, 1, 12, 10, 2, 12, 5, 10, 12},
     /* 2,6 */ {9, 0, 1, 2, 3, 12, 10, 2, 12, 5, 10, 12, 4, 5, 12, 8, 4, 12, 7, 8, 12, 6, 7, 12, 11, 6, 12, 3, 11, 12},
     /* 3,4 */ {6, 5, 10, 12, 1, 9, 12, 9, 4, 12, 4, 7, 12, 7, 11, 12, 11, 2, 12, 2, 3, 12, 3, 8, 12, 8, 0, 12, 0, 1},
     /* 3,5 */ {6, 5, 10, 12, 4, 7, 12, 7, 11, 12, 11, 2, 12, 2, 1, 12, 1, 9, 12, 9, 0, 12, 0, 3, 12, 3, 8, 12, 8, 4},
     /* 3,6 */ {9, 0, 1, 12, 2, 3, 12, 3, 8, 12, 8, 4, 12, 4, 5, 12, 5, 10, 12, 10, 6, 12, 6, 7, 12, 7, 11, 12, 11, 2},
     /* 4,5 */ {6, 5, 10, 4, 7, 12, 9, 4, 12, 1, 9, 12, 2, 1, 12, 11, 2, 12, 3, 11, 12, 0, 3, 12, 8, 0, 12, 7, 8, 12},
     /* 4,6 */ {2, 3, 11, 0, 1, 12, 8, 0, 12, 7, 8, 12, 6, 7, 12, 10, 6, 12, 5, 10, 12, 4, 5, 12, 9, 4, 12, 1, 9, 12}}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 13.3, inverted
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
/* 13.3 */
const schar MarchingCubes::tiling13_3_[2][12][30] = {
    /* 165: 0,    2,       5,    7,  */
    {/* 1,2 */ {3, 2, 11, 8, 7, 12, 0, 8, 12, 1, 0, 12, 10, 1, 12, 6, 10, 12, 5, 6, 12, 9, 5, 12, 4, 9, 12, 7, 4, 12},
     /* 1,4 */ {5, 6, 10, 12, 2, 11, 12, 11, 7, 12, 7, 4, 12, 4, 9, 12, 9, 1, 12, 1, 0, 12, 0, 8, 12, 8, 3, 12, 3, 2},
     /* 1,5 */ {10, 5, 6, 12, 7, 4, 12, 4, 9, 12, 9, 1, 12, 1, 2, 12, 2, 11, 12, 11, 3, 12, 3, 0, 12, 0, 8, 12, 8, 7},
     /* 1,6 */ {11, 3, 2, 12, 1, 0, 12, 0, 8, 12, 8, 7, 12, 7, 6, 12, 6, 10, 12, 10, 5, 12, 5, 4, 12, 4, 9, 12, 9, 1},
     /* 2,3 */ {7, 4, 8, 11, 3, 12, 6, 11, 12, 5, 6, 12, 9, 5, 12, 0, 9, 12, 1, 0, 12, 10, 1, 12, 2, 10, 12, 3, 2, 12},
     /* 2,5 */ {7, 4, 8, 5, 6, 12, 9, 5, 12, 0, 9, 12, 3, 0, 12, 11, 3, 12, 2, 11, 12, 1, 2, 12, 10, 1, 12, 6, 10, 12},
     /* 2,6 */ {11, 3, 2, 1, 0, 12, 10, 1, 12, 6, 10, 12, 7, 6, 12, 8, 7, 12, 4, 8, 12, 5, 4, 12, 9, 5, 12, 0, 9, 12},
     /* 3,4 */ {1, 0, 9, 12, 4, 8, 12, 8, 3, 12, 3, 2, 12, 2, 10, 12, 10, 5, 12, 5, 6, 12, 6, 11, 12, 11, 7, 12, 7, 4},
     /* 3,5 */ {7, 4, 8, 12, 5, 6, 12, 6, 11, 12, 11, 3, 12, 3, 0, 12, 0, 9, 12, 9, 1, 12, 1, 2, 12, 2, 10, 12, 10, 5},
     /* 3,6 */ {1, 0, 9, 12, 3, 2, 12, 2, 10, 12, 10, 5, 12, 5, 4, 12, 4, 8, 12, 8, 7, 12, 7, 6, 12, 6, 11, 12, 11, 3},
     /* 4,5 */ {10, 5, 6, 7, 4, 12, 11, 7, 12, 2, 11, 12, 1, 2, 12, 9, 1, 12, 0, 9, 12, 3, 0, 12, 8, 3, 12, 4, 8, 12},
     /* 4,6 */ {9, 1, 0, 3, 2, 12, 8, 3, 12, 4, 8, 12, 5, 4, 12, 10, 5, 12, 6, 10, 12, 7, 6, 12, 11, 7, 12, 2, 11, 12}},
    /*  90:    1,    3, 4,    6,     */
    {/* 1,2 */ {0, 3, 8, 9, 4, 12, 1, 9, 12, 2, 1, 12, 11, 2, 12, 7, 11, 12, 6, 7, 12, 10, 6, 12, 5, 10, 12, 4, 5, 12},
     /* 1,4 */ {11, 6, 7, 12, 3, 8, 12, 8, 4, 12, 4, 5, 12, 5, 10, 12, 10, 2, 12, 2, 1, 12, 1, 9, 12, 9, 0, 12, 0, 3},
     /* 1,5 */ {6, 7, 11, 12, 4, 5, 12, 5, 10, 12, 10, 2, 12, 2, 3, 12, 3, 8, 12, 8, 0, 12, 0, 1, 12, 1, 9, 12, 9, 4},
     /* 1,6 */ {8, 0, 3, 12, 2, 1, 12, 1, 9, 12, 9, 4, 12, 4, 7, 12, 7, 11, 12, 11, 6, 12, 6, 5, 12, 5, 10, 12, 10, 2},
     /* 2,3 */ {4, 5, 9, 8, 0, 12, 7, 8, 12, 6, 7, 12, 10, 6, 12, 1, 10, 12, 2, 1, 12, 11, 2, 12, 3, 11, 12, 0, 3, 12},
     /* 2,5 */ {4, 5, 9, 6, 7, 12, 10, 6, 12, 1, 10, 12, 0, 1, 12, 8, 0, 12, 3, 8, 12, 2, 3, 12, 11, 2, 12, 7, 11, 12},
     /* 2,6 */ {8, 0, 3, 2, 1, 12, 11, 2, 12, 7, 11, 12, 4, 7, 12, 9, 4, 12, 5, 9, 12, 6, 5, 12, 10, 6, 12, 1, 10, 12},
     /* 3,4 */ {2, 1, 10, 12, 5, 9, 12, 9, 0, 12, 0, 3, 12, 3, 11, 12, 11, 6, 12, 6, 7, 12, 7, 8, 12, 8, 4, 12, 4, 5},
     /* 3,5 */ {4, 5, 9, 12, 6, 7, 12, 7, 8, 12, 8, 0, 12, 0, 1, 12, 1, 10, 12, 10, 2, 12, 2, 3, 12, 3, 11, 12, 11, 6},
     /* 3,6 */ {2, 1, 10, 12, 0, 3, 12, 3, 11, 12, 11, 6, 12, 6, 5, 12, 5, 9, 12, 9, 4, 12, 4, 7, 12, 7, 8, 12, 8, 0},
     /* 4,5 */ {6, 7, 11, 4, 5, 12, 8, 4, 12, 3, 8, 12, 2, 3, 12, 10, 2, 12, 1, 10, 12, 0, 1, 12, 9, 0, 12, 5, 9, 12},
     /* 4,6 */ {10, 2, 1, 0, 3, 12, 9, 0, 12, 5, 9, 12, 6, 5, 12, 11, 6, 12, 7, 11, 12, 4, 7, 12, 8, 4, 12, 3, 8, 12}}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 13.4
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
/* 13.4 */
const schar MarchingCubes::tiling13_4[2][4][36] = {
    /* 165: 0,    2,       5,    7,  */ {/* 1,2,6 */ {12, 2, 10, 12, 10, 5, 12, 5, 6, 12, 6, 11, 12, 11, 7, 12, 7, 4,
                                                      12, 4, 8,  12, 8,  3, 12, 3, 0, 12, 0, 9,  12, 9,  1, 12, 1, 2},
                                         /* 1,4,5 */ {11, 3, 12, 6, 11, 12, 7, 6, 12, 8,  7, 12, 4, 8,  12, 5, 4, 12,
                                                      9,  5, 12, 0, 9,  12, 1, 0, 12, 10, 1, 12, 2, 10, 12, 3, 2, 12},
                                         /* 2,3,5 */ {9,  1, 12, 4, 9,  12, 5, 4, 12, 10, 5, 12, 6, 10, 12, 7, 6, 12,
                                                      11, 7, 12, 2, 11, 12, 3, 2, 12, 8,  3, 12, 0, 8,  12, 1, 0, 12},
                                         /* 3,4,6 */ {12, 0, 8,  12, 8,  7, 12, 7, 4, 12, 4, 9,  12, 9,  5, 12, 5, 6,
                                                      12, 6, 10, 12, 10, 1, 12, 1, 2, 12, 2, 11, 12, 11, 3, 12, 3, 0}},
    /*  90:    1,    3, 4,    6,     */ {/* 1,2,6 */ {12, 3, 11, 12, 11, 6, 12, 6, 7, 12, 7, 8,  12, 8,  4, 12, 4, 5,
                                                      12, 5, 9,  12, 9,  0, 12, 0, 1, 12, 1, 10, 12, 10, 2, 12, 2, 3},
                                         /* 1,4,5 */ {8,  0, 12, 7, 8,  12, 4, 7, 12, 9,  4, 12, 5, 9,  12, 6, 5, 12,
                                                      10, 6, 12, 1, 10, 12, 2, 1, 12, 11, 2, 12, 3, 11, 12, 0, 3, 12},
                                         /* 2,3,5 */ {10, 2, 12, 5, 10, 12, 6, 5, 12, 11, 6, 12, 7, 11, 12, 4, 7, 12,
                                                      8,  4, 12, 3, 8,  12, 0, 3, 12, 9,  0, 12, 1, 9,  12, 2, 1, 12},
                                         /* 3,4,6 */ {12, 1, 9,  12, 9,  4, 12, 4, 5, 12, 5, 10, 12, 10, 6, 12, 6, 7,
                                                      12, 7, 11, 12, 11, 2, 12, 2, 3, 12, 3, 8,  12, 8,  0, 12, 0, 1}}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 13.5.1
 * The support edge for the interior test is marked as the 1st column.
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
/* 13.5.1 */
const schar MarchingCubes::tiling13_5_1[2][4][18] = {
    /* 165: 0,    2,       5,    7,  */ {/* 1,2,5 */ {7, 6, 11, 1, 0, 9, 10, 3, 2, 3, 10, 5, 3, 5, 8, 4, 8, 5},
                                         /* 1,4,6 */ {1, 2, 10, 7, 4, 8, 3, 0, 11, 6, 11, 0, 9, 6, 0, 6, 9, 5},
                                         /* 2,3,6 */ {3, 0, 8, 5, 6, 10, 1, 2, 9, 4, 9, 2, 11, 4, 2, 4, 11, 7},
                                         /* 3,4,5 */ {5, 4, 9, 3, 2, 11, 8, 1, 0, 1, 8, 7, 1, 7, 10, 6, 10, 7}},
    /*  90:    1,    3, 4,    6,     */ {/* 1,2,5 */ {4, 7, 8, 2, 1, 10, 11, 0, 3, 0, 11, 6, 0, 6, 9, 5, 9, 6},
                                         /* 1,4,6 */ {2, 3, 11, 4, 5, 9, 0, 1, 8, 7, 8, 1, 10, 7, 1, 7, 10, 6},
                                         /* 2,3,6 */ {0, 1, 9, 6, 7, 11, 2, 3, 10, 5, 10, 3, 8, 5, 3, 5, 8, 4},
                                         /* 3,4,5 */ {6, 5, 10, 0, 3, 8, 9, 2, 1, 2, 9, 4, 2, 4, 11, 7, 11, 4}}};

//_____________________________________________________________________________
/**
 * \brief tiling table for case 13.5.2
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
/* 13.5.2 */
const schar MarchingCubes::tiling13_5_2[2][4][30] = {
    /* 165: 0,    2,       5,    7,  */ {
        /* 1,2,5 */ {1, 0, 9, 7, 4, 8, 7, 8, 3, 7, 3, 11, 2, 11, 3, 11, 2, 10, 11, 10, 6, 5, 6, 10, 6, 5, 7, 4, 7, 5},
        /* 1,4,6 */ {7, 4, 8, 11, 3, 2, 6, 11, 2, 10, 6, 2, 6, 10, 5, 9, 5, 10, 1, 9, 10, 9, 1, 0, 2, 0, 1, 0, 2, 3},
        /* 2,3,6 */ {5, 6, 10, 9, 1, 0, 4, 9, 0, 8, 4, 0, 4, 8, 7, 11, 7, 8, 3, 11, 8, 11, 3, 2, 0, 2, 3, 2, 0, 1},
        /* 3,4,5 */ {3, 2, 11, 5, 6, 10, 5, 10, 1, 5, 1, 9, 0, 9, 1, 9, 0, 8, 9, 8, 4, 4, 8, 7, 4, 7, 5, 6, 5, 7}},
    /*  90:    1,    3, 4,    6,     */ {
        /* 1,2,5 */ {2, 1, 10, 4, 5, 9, 4, 9, 0, 4, 0, 8, 3, 8, 0, 8, 3, 11, 8, 11, 7, 6, 7, 11, 7, 6, 4, 5, 4, 6},
        /* 1,4,6 */ {4, 5, 9, 8, 0, 3, 7, 8, 3, 11, 7, 3, 7, 11, 6, 10, 6, 11, 2, 10, 11, 10, 2, 1, 3, 1, 2, 1, 3, 0},
        /* 2,3,6 */ {6, 7, 11, 10, 2, 1, 5, 10, 1, 9, 5, 1, 5, 9, 4, 8, 4, 9, 0, 8, 9, 8, 0, 3, 1, 3, 0, 3, 1, 2},
        /* 3,4,5 */ {0, 3, 8, 6, 7, 11, 6, 11, 2, 6, 2, 10, 1, 10, 2, 10, 1, 9, 10, 9, 5, 5, 9, 4, 5, 4, 6, 7, 6, 4}}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief tiling table for case 14
 * For each of the case above, the specific triangulation of the edge
 * intersection points is given.
 * When a case is ambiguous, there is an auxiliary table that contains
 * the face number to test and the tiling table contains the specific
 * triangulations depending on the results
 * A minus sign means to invert the result of the test.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::tiling14[12][12] = {
    /*  71: 0, 1, 2,          6,     */ {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8},
    /*  43: 0, 1,    3,    5,        */ {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5},
    /* 147: 0, 1,       4,       7,  */ {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6},
    /*  29: 0,    2, 3, 4,           */ {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4},
    /* 201: 0,       3,       6, 7,  */ {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5},
    /* 113: 0,          4, 5, 6,     */ {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10},
    /* 142:    1, 2, 3,          7,  */ {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7},
    /*  54:    1, 2,    4, 5,        */ {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2},
    /* 226:    1,          5, 6, 7,  */ {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11},
    /* 108:       2, 3,    5, 6,     */ {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3},
    /* 212:       2,    4,    6, 7,  */ {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8},
    /* 184:          3, 4, 5,    7,  */ {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2}};
//_____________________________________________________________________________

//_____________________________________________________________________________
/**
 * \brief original Marching Cubes implementation
 * For each of the possible vertex states listed in this table there is a
 * specific triangulation of the edge intersection points.  The table lists
 * all of them in the form of 0-5 edge triples with the list terminated by
 * the invalid value -1.  For example: casesClassic[3] list the 2 triangles
 * formed when cube[0] and cube[1] are inside of the surface, but the rest of
 * the cube is not.
 */
//-----------------------------------------------------------------------------
const schar MarchingCubes::casesClassic[256][16] = {
    /*   0:                          */ {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*   1: 0,                       */ {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*   2:    1,                    */ {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*   3: 0, 1,                    */ {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*   4:       2,                 */ {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*   5: 0,    2,                 */ {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*   6:    1, 2,                 */ {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*   7: 0, 1, 2,                 */ {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    /*   8:          3,              */ {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*   9: 0,       3,              */ {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  10:    1,    3,              */ {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  11: 0, 1,    3,              */ {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    /*  12:       2, 3,              */ {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  13: 0,    2, 3,              */ {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    /*  14:    1, 2, 3,              */ {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    /*  15: 0, 1, 2, 3,              */ {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  16:             4,           */ {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  17: 0,          4,           */ {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  18:    1,       4,           */ {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  19: 0, 1,       4,           */ {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    /*  20:       2,    4,           */ {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  21: 0,    2,    4,           */ {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    /*  22:    1, 2,    4,           */ {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    /*  23: 0, 1, 2,    4,           */ {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    /*  24:          3, 4,           */ {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  25: 0,       3, 4,           */ {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    /*  26:    1,    3, 4,           */ {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    /*  27: 0, 1,    3, 4,           */ {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    /*  28:       2, 3, 4,           */ {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    /*  29: 0,    2, 3, 4,           */ {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    /*  30:    1, 2, 3, 4,           */ {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    /*  31: 0, 1, 2, 3, 4,           */ {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    /*  32:                5,        */ {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  33: 0,             5,        */ {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  34:    1,          5,        */ {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  35: 0, 1,          5,        */ {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    /*  36:       2,       5,        */ {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  37: 0,    2,       5,        */ {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    /*  38:    1, 2,       5,        */ {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    /*  39: 0, 1, 2,       5,        */ {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    /*  40:          3,    5,        */ {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  41: 0,       3,    5,        */ {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    /*  42:    1,    3,    5,        */ {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    /*  43: 0, 1,    3,    5,        */ {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    /*  44:       2, 3,    5,        */ {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    /*  45: 0,    2, 3,    5,        */ {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    /*  46:    1, 2, 3,    5,        */ {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    /*  47: 0, 1, 2, 3,    5,        */ {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    /*  48:             4, 5,        */ {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  49: 0,          4, 5,        */ {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    /*  50:    1,       4, 5,        */ {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    /*  51: 0, 1,       4, 5,        */ {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  52:       2,    4, 5,        */ {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    /*  53: 0,    2,    4, 5,        */ {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    /*  54:    1, 2,    4, 5,        */ {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    /*  55: 0, 1, 2,    4, 5,        */ {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    /*  56:          3, 4, 5,        */ {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    /*  57: 0,       3, 4, 5,        */ {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    /*  58:    1,    3, 4, 5,        */ {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    /*  59: 0, 1,    3, 4, 5,        */ {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    /*  60:       2, 3, 4, 5,        */ {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    /*  61: 0,    2, 3, 4, 5,        */ {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    /*  62:    1, 2, 3, 4, 5,        */ {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    /*  63: 0, 1, 2, 3, 4, 5,        */ {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  64:                   6,     */ {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  65: 0,                6,     */ {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  66:    1,             6,     */ {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  67: 0, 1,             6,     */ {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    /*  68:       2,          6,     */ {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  69: 0,    2,          6,     */ {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    /*  70:    1, 2,          6,     */ {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    /*  71: 0, 1, 2,          6,     */ {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    /*  72:          3,       6,     */ {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  73: 0,       3,       6,     */ {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    /*  74:    1,    3,       6,     */ {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    /*  75: 0, 1,    3,       6,     */ {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    /*  76:       2, 3,       6,     */ {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    /*  77: 0,    2, 3,       6,     */ {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    /*  78:    1, 2, 3,       6,     */ {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    /*  79: 0, 1, 2, 3,       6,     */ {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    /*  80:             4,    6,     */ {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  81: 0,          4,    6,     */ {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    /*  82:    1,       4,    6,     */ {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    /*  83: 0, 1,       4,    6,     */ {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    /*  84:       2,    4,    6,     */ {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    /*  85: 0,    2,    4,    6,     */ {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    /*  86:    1, 2,    4,    6,     */ {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    /*  87: 0, 1, 2,    4,    6,     */ {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    /*  88:          3, 4,    6,     */ {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    /*  89: 0,       3, 4,    6,     */ {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    /*  90:    1,    3, 4,    6,     */ {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    /*  91: 0, 1,    3, 4,    6,     */ {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    /*  92:       2, 3, 4,    6,     */ {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    /*  93: 0,    2, 3, 4,    6,     */ {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    /*  94:    1, 2, 3, 4,    6,     */ {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    /*  95: 0, 1, 2, 3, 4,    6,     */ {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    /*  96:                5, 6,     */ {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /*  97: 0,             5, 6,     */ {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    /*  98:    1,          5, 6,     */ {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    /*  99: 0, 1,          5, 6,     */ {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    /* 100:       2,       5, 6,     */ {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    /* 101: 0,    2,       5, 6,     */ {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    /* 102:    1, 2,       5, 6,     */ {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 103: 0, 1, 2,       5, 6,     */ {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    /* 104:          3,    5, 6,     */ {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    /* 105: 0,       3,    5, 6,     */ {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    /* 106:    1,    3,    5, 6,     */ {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    /* 107: 0, 1,    3,    5, 6,     */ {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    /* 108:       2, 3,    5, 6,     */ {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    /* 109: 0,    2, 3,    5, 6,     */ {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    /* 110:    1, 2, 3,    5, 6,     */ {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    /* 111: 0, 1, 2, 3,    5, 6,     */ {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 112:             4, 5, 6,     */ {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    /* 113: 0,          4, 5, 6,     */ {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    /* 114:    1,       4, 5, 6,     */ {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    /* 115: 0, 1,       4, 5, 6,     */ {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    /* 116:       2,    4, 5, 6,     */ {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    /* 117: 0,    2,    4, 5, 6,     */ {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    /* 118:    1, 2,    4, 5, 6,     */ {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    /* 119: 0, 1, 2,    4, 5, 6,     */ {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 120:          3, 4, 5, 6,     */ {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    /* 121: 0,       3, 4, 5, 6,     */ {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    /* 122:    1,    3, 4, 5, 6,     */ {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    /* 123: 0, 1,    3, 4, 5, 6,     */ {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    /* 124:       2, 3, 4, 5, 6,     */ {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    /* 125: 0,    2, 3, 4, 5, 6,     */ {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 126:    1, 2, 3, 4, 5, 6,     */ {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    /* 127: 0, 1, 2, 3, 4, 5, 6,     */ {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 128:                      7,  */ {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 129: 0,                   7,  */ {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 130:    1,                7,  */ {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 131: 0, 1,                7,  */ {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    /* 132:       2,             7,  */ {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 133: 0,    2,             7,  */ {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    /* 134:    1, 2,             7,  */ {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    /* 135: 0, 1, 2,             7,  */ {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    /* 136:          3,          7,  */ {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 137: 0,       3,          7,  */ {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    /* 138:    1,    3,          7,  */ {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    /* 139: 0, 1,    3,          7,  */ {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    /* 140:       2, 3,          7,  */ {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    /* 141: 0,    2, 3,          7,  */ {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    /* 142:    1, 2, 3,          7,  */ {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    /* 143: 0, 1, 2, 3,          7,  */ {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    /* 144:             4,       7,  */ {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 145: 0,          4,       7,  */ {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    /* 146:    1,       4,       7,  */ {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    /* 147: 0, 1,       4,       7,  */ {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    /* 148:       2,    4,       7,  */ {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    /* 149: 0,    2,    4,       7,  */ {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    /* 150:    1, 2,    4,       7,  */ {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    /* 151: 0, 1, 2,    4,       7,  */ {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    /* 152:          3, 4,       7,  */ {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    /* 153: 0,       3, 4,       7,  */ {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 154:    1,    3, 4,       7,  */ {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    /* 155: 0, 1,    3, 4,       7,  */ {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    /* 156:       2, 3, 4,       7,  */ {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    /* 157: 0,    2, 3, 4,       7,  */ {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    /* 158:    1, 2, 3, 4,       7,  */ {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    /* 159: 0, 1, 2, 3, 4,       7,  */ {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 160:                5,    7,  */ {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 161: 0,             5,    7,  */ {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    /* 162:    1,          5,    7,  */ {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    /* 163: 0, 1,          5,    7,  */ {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    /* 164:       2,       5,    7,  */ {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    /* 165: 0,    2,       5,    7,  */ {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    /* 166:    1, 2,       5,    7,  */ {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    /* 167: 0, 1, 2,       5,    7,  */ {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    /* 168:          3,    5,    7,  */ {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    /* 169: 0,       3,    5,    7,  */ {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    /* 170:    1,    3,    5,    7,  */ {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    /* 171: 0, 1,    3,    5,    7,  */ {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    /* 172:       2, 3,    5,    7,  */ {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    /* 173: 0,    2, 3,    5,    7,  */ {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    /* 174:    1, 2, 3,    5,    7,  */ {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    /* 175: 0, 1, 2, 3,    5,    7,  */ {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    /* 176:             4, 5,    7,  */ {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    /* 177: 0,          4, 5,    7,  */ {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    /* 178:    1,       4, 5,    7,  */ {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    /* 179: 0, 1,       4, 5,    7,  */ {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    /* 180:       2,    4, 5,    7,  */ {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    /* 181: 0,    2,    4, 5,    7,  */ {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    /* 182:    1, 2,    4, 5,    7,  */ {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    /* 183: 0, 1, 2,    4, 5,    7,  */ {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    /* 184:          3, 4, 5,    7,  */ {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    /* 185: 0,       3, 4, 5,    7,  */ {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    /* 186:    1,    3, 4, 5,    7,  */ {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    /* 187: 0, 1,    3, 4, 5,    7,  */ {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 188:       2, 3, 4, 5,    7,  */ {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    /* 189: 0,    2, 3, 4, 5,    7,  */ {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    /* 190:    1, 2, 3, 4, 5,    7,  */ {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 191: 0, 1, 2, 3, 4, 5,    7,  */ {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 192:                   6, 7,  */ {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 193: 0,                6, 7,  */ {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    /* 194:    1,             6, 7,  */ {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    /* 195: 0, 1,             6, 7,  */ {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    /* 196:       2,          6, 7,  */ {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    /* 197: 0,    2,          6, 7,  */ {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    /* 198:    1, 2,          6, 7,  */ {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    /* 199: 0, 1, 2,          6, 7,  */ {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    /* 200:          3,       6, 7,  */ {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    /* 201: 0,       3,       6, 7,  */ {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    /* 202:    1,    3,       6, 7,  */ {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    /* 203: 0, 1,    3,       6, 7,  */ {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    /* 204:       2, 3,       6, 7,  */ {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 205: 0,    2, 3,       6, 7,  */ {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    /* 206:    1, 2, 3,       6, 7,  */ {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    /* 207: 0, 1, 2, 3,       6, 7,  */ {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 208:             4,    6, 7,  */ {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    /* 209: 0,          4,    6, 7,  */ {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    /* 210:    1,       4,    6, 7,  */ {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    /* 211: 0, 1,       4,    6, 7,  */ {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    /* 212:       2,    4,    6, 7,  */ {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    /* 213: 0,    2,    4,    6, 7,  */ {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    /* 214:    1, 2,    4,    6, 7,  */ {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    /* 215: 0, 1, 2,    4,    6, 7,  */ {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 216:          3, 4,    6, 7,  */ {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    /* 217: 0,       3, 4,    6, 7,  */ {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    /* 218:    1,    3, 4,    6, 7,  */ {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    /* 219: 0, 1,    3, 4,    6, 7,  */ {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    /* 220:       2, 3, 4,    6, 7,  */ {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    /* 221: 0,    2, 3, 4,    6, 7,  */ {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 222:    1, 2, 3, 4,    6, 7,  */ {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    /* 223: 0, 1, 2, 3, 4,    6, 7,  */ {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 224:                5, 6, 7,  */ {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    /* 225: 0,             5, 6, 7,  */ {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    /* 226:    1,          5, 6, 7,  */ {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    /* 227: 0, 1,          5, 6, 7,  */ {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    /* 228:       2,       5, 6, 7,  */ {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    /* 229: 0,    2,       5, 6, 7,  */ {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    /* 230:    1, 2,       5, 6, 7,  */ {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    /* 231: 0, 1, 2,       5, 6, 7,  */ {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    /* 232:          3,    5, 6, 7,  */ {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    /* 233: 0,       3,    5, 6, 7,  */ {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    /* 234:    1,    3,    5, 6, 7,  */ {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    /* 235: 0, 1,    3,    5, 6, 7,  */ {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 236:       2, 3,    5, 6, 7,  */ {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    /* 237: 0,    2, 3,    5, 6, 7,  */ {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    /* 238:    1, 2, 3,    5, 6, 7,  */ {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 239: 0, 1, 2, 3,    5, 6, 7,  */ {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 240:             4, 5, 6, 7,  */ {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 241: 0,          4, 5, 6, 7,  */ {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    /* 242:    1,       4, 5, 6, 7,  */ {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    /* 243: 0, 1,       4, 5, 6, 7,  */ {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 244:       2,    4, 5, 6, 7,  */ {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    /* 245: 0,    2,    4, 5, 6, 7,  */ {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    /* 246:    1, 2,    4, 5, 6, 7,  */ {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 247: 0, 1, 2,    4, 5, 6, 7,  */ {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 248:          3, 4, 5, 6, 7,  */ {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    /* 249: 0,       3, 4, 5, 6, 7,  */ {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 250:    1,    3, 4, 5, 6, 7,  */ {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    /* 251: 0, 1,    3, 4, 5, 6, 7,  */ {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 252:       2, 3, 4, 5, 6, 7,  */ {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 253: 0,    2, 3, 4, 5, 6, 7,  */ {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 254:    1, 2, 3, 4, 5, 6, 7,  */ {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    /* 255: 0, 1, 2, 3, 4, 5, 6, 7,  */ {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};
//_____________________________________________________________________________
