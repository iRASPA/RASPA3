module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>
#endif

module triquintic_derivatives_external_field;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <array>;
import <tuple>;
import <cstddef>;
#endif

import double3;

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
Potentials::potentialTriquinticSecondOrderPolynomialTestFunction(const double3 &pos)
{
  // Factor out common terms
  const double x2 = pos.x * pos.x;
  const double y2 = pos.y * pos.y;
  const double z2 = pos.z * pos.z;
  const double r2 = x2 + y2 + z2;

  return {r2, {{2.0 * pos.x, 2.0 * pos.y, 2.0 * pos.z}}, {{{{2, 0, 0}}, {{0, 2, 0}}, {{0, 0, 2}}}}, {}, {}, {}, {}};
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
Potentials::potentialTriquinticThirdOrderPolynomialTestFunction(const double3 &pos)
{
  // Factor out common terms
  const double x2 = pos.x * pos.x;
  const double y2 = pos.y * pos.y;
  const double z2 = pos.z * pos.z;
  const double xy = pos.x * pos.y;
  const double xz = pos.x * pos.z;
  const double yz = pos.y * pos.z;
  const double x3 = pos.x * x2;

  return {x3 + 2.0 * x2 * pos.y * z2 - y2 * pos.z + 1.0,
          {{3 * x2 + 4 * xy * z2, -2 * yz + 2 * x2 * z2, -y2 + 4 * x2 * yz}},
          {{{{6 * pos.x + 4 * pos.y * z2, 4 * pos.x * z2, 8 * xy * pos.z}},
            {{4 * pos.x * z2, -2 * pos.z, -2 * pos.y + 4 * x2 * pos.z}},
            {{8 * xy * pos.z, -2 * pos.y + 4 * x2 * pos.z, 4 * x2 * pos.y}}}},
          {{{{{{6, 4 * z2, 8 * yz}}, {{4 * z2, 0, 8 * xz}}, {{8 * yz, 8 * xz, 8 * xy}}}},
            {{{{4 * z2, 0, 8 * xz}}, {{0, 0, -2}}, {{8 * xz, -2, 4 * x2}}}},
            {{{{8 * yz, 8 * xz, 8 * xy}}, {{8 * xz, -2, 4 * x2}}, {{8 * xy, 4 * x2, 0}}}}}},
          {{{{{{{{0, 0, 0}}, {{0, 0, 8 * pos.z}}, {{0, 8 * pos.z, 8 * pos.y}}}},
              {{{{0, 0, 8 * pos.z}}, {{0, 0, 0}}, {{8 * pos.z, 0, 8 * pos.x}}}},
              {{{{0, 8 * pos.z, 8 * pos.y}}, {{8 * pos.z, 0, 8 * pos.x}}, {{8 * pos.y, 8 * pos.x, 0}}}}}},
            {{{{{{0, 0, 8 * pos.z}}, {{0, 0, 0}}, {{8 * pos.z, 0, 8 * pos.x}}}},
              {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
              {{{{8 * pos.z, 0, 8 * pos.x}}, {{0, 0, 0}}, {{8 * pos.x, 0, 0}}}}}},
            {{{{{{0, 8 * pos.z, 8 * pos.y}}, {{8 * pos.z, 0, 8 * pos.x}}, {{8 * pos.y, 8 * pos.x, 0}}}},
              {{{{8 * pos.z, 0, 8 * pos.x}}, {{0, 0, 0}}, {{8 * pos.x, 0, 0}}}},
              {{{{8 * pos.y, 8 * pos.x, 0}}, {{8 * pos.x, 0, 0}}, {{0, 0, 0}}}}}}}},
          {{{{{{{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 8}}}},
                {{{{0, 0, 0}}, {{0, 0, 8}}, {{0, 8, 0}}}}}},
              {{{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 8}}}},
                {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                {{{{0, 0, 8}}, {{0, 0, 0}}, {{8, 0, 0}}}}}},
              {{{{{{0, 0, 0}}, {{0, 0, 8}}, {{0, 8, 0}}}},
                {{{{0, 0, 8}}, {{0, 0, 0}}, {{8, 0, 0}}}},
                {{{{0, 8, 0}}, {{8, 0, 0}}, {{0, 0, 0}}}}}}}},
            {{{{{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 8}}}},
                {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                {{{{0, 0, 8}}, {{0, 0, 0}}, {{8, 0, 0}}}}}},
              {{{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
              {{{{{{0, 0, 8}}, {{0, 0, 0}}, {{8, 0, 0}}}},
                {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                {{{{8, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}},
            {{{{{{{{0, 0, 0}}, {{0, 0, 8}}, {{0, 8, 0}}}},
                {{{{0, 0, 8}}, {{0, 0, 0}}, {{8, 0, 0}}}},
                {{{{0, 8, 0}}, {{8, 0, 0}}, {{0, 0, 0}}}}}},
              {{{{{{0, 0, 8}}, {{0, 0, 0}}, {{8, 0, 0}}}},
                {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                {{{{8, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
              {{{{{{0, 8, 0}}, {{8, 0, 0}}, {{0, 0, 0}}}},
                {{{{8, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}}}},
          {}};
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
Potentials::potentialTriquinticFourthOrderPolynomialTestFunction(const double3 &pos)
{
  // Factor out common terms
  const double x2 = pos.x * pos.x;
  const double y2 = pos.y * pos.y;
  const double z2 = pos.z * pos.z;
  const double xy = pos.x * pos.y;
  const double xz = pos.x * pos.z;
  const double yz = pos.y * pos.z;
  const double r2 = x2 + y2 + z2;
  const double r4 = r2 * r2;

  return {
      r4,
      {{4 * pos.x * r2, 4 * pos.y * r2, 4 * pos.z * r2}},
      {{{{8 * x2 + 4 * r2, 8 * xy, 8 * xz}}, {{8 * xy, 8 * y2 + 4 * r2, 8 * yz}}, {{8 * xz, 8 * yz, 8 * z2 + 4 * r2}}}},
      {{{{{{24 * pos.x, 8 * pos.y, 8 * pos.z}}, {{8 * pos.y, 8 * pos.x, 0}}, {{8 * pos.z, 0, 8 * pos.x}}}},
        {{{{8 * pos.y, 8 * pos.x, 0}}, {{8 * pos.x, 24 * pos.y, 8 * pos.z}}, {{0, 8 * pos.z, 8 * pos.y}}}},
        {{{{8 * pos.z, 0, 8 * pos.x}}, {{0, 8 * pos.z, 8 * pos.y}}, {{8 * pos.x, 8 * pos.y, 24 * pos.z}}}}}},
      {{{{{{{{24, 0, 0}}, {{0, 8, 0}}, {{0, 0, 8}}}},
          {{{{0, 8, 0}}, {{8, 0, 0}}, {{0, 0, 0}}}},
          {{{{0, 0, 8}}, {{0, 0, 0}}, {{8, 0, 0}}}}}},
        {{{{{{0, 8, 0}}, {{8, 0, 0}}, {{0, 0, 0}}}},
          {{{{8, 0, 0}}, {{0, 24, 0}}, {{0, 0, 8}}}},
          {{{{0, 0, 0}}, {{0, 0, 8}}, {{0, 8, 0}}}}}},
        {{{{{{0, 0, 8}}, {{0, 0, 0}}, {{8, 0, 0}}}},
          {{{{0, 0, 0}}, {{0, 0, 8}}, {{0, 8, 0}}}},
          {{{{8, 0, 0}}, {{0, 8, 0}}, {{0, 0, 24}}}}}}}},
      {},
      {}};
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
Potentials::potentialTriquinticFifthOrderPolynomialTestFunction(const double3 &pos)
{
  // Factor out common terms
  const double x2 = pos.x * pos.x;
  const double y2 = pos.y * pos.y;
  const double z2 = pos.z * pos.z;
  const double xy = pos.x * pos.y;
  const double xz = pos.x * pos.z;
  const double yz = pos.y * pos.z;
  const double x3 = pos.x * x2;
  const double y3 = pos.y * y2;
  const double x4 = x2 * x2;

  return {x2 * x3 * y2 * yz + 2.0 * pos.x * y3 + 4.0 * x2 * pos.y * z2,
          {{2 * y3 + 5 * x4 * y2 * yz + 8 * xy * z2, 6 * pos.x * y2 + 3 * x2 * x3 * y2 * pos.z + 4 * x2 * z2,
            x2 * x3 * y3 + 8 * x2 * yz}},
          {{{{20 * x3 * y2 * yz + 8 * pos.y * z2, 6 * y2 + 15 * x4 * y2 * pos.z + 8 * pos.x * z2,
              5 * x4 * y3 + 16 * xy * pos.z}},
            {{6 * y2 + 15 * x4 * y2 * pos.z + 8 * pos.x * z2, 12 * xy + 6 * x4 * xy * pos.z,
              3 * x2 * x3 * y2 + 8 * x2 * pos.z}},
            {{5 * x4 * y3 + 16 * xy * pos.z, 3 * x2 * x3 * y2 + 8 * x2 * pos.z, 8 * x2 * pos.y}}}},
          {{{{{{60 * x2 * y2 * yz, 60 * x3 * y2 * pos.z + 8 * z2, 20 * x3 * y3 + 16 * yz}},
              {{60 * x3 * y2 * pos.z + 8 * z2, 12 * pos.y + 30 * x4 * yz, 15 * x4 * y2 + 16 * xz}},
              {{20 * x3 * y3 + 16 * yz, 15 * x4 * y2 + 16 * xz, 16 * xy}}}},
            {{{{60 * x3 * y2 * pos.z + 8 * z2, 12 * pos.y + 30 * x4 * yz, 15 * x4 * y2 + 16 * xz}},
              {{12 * pos.y + 30 * x4 * yz, 12 * pos.x + 6 * x4 * xz, 6 * x4 * xy}},
              {{15 * x4 * y2 + 16 * xz, 6 * x4 * xy, 8 * x2}}}},
            {{{{20 * x3 * y3 + 16 * yz, 15 * x4 * y2 + 16 * xz, 16 * xy}},
              {{15 * x4 * y2 + 16 * xz, 6 * x4 * xy, 8 * x2}},
              {{16 * xy, 8 * x2, 0}}}}}},
          {{{{{{{{120 * pos.x * y2 * yz, 180 * x2 * y2 * pos.z, 60 * x2 * y3}},
                {{180 * x2 * y2 * pos.z, 120 * x2 * xy * pos.z, 60 * x3 * y2 + 16 * pos.z}},
                {{60 * x2 * y3, 60 * x3 * y2 + 16 * pos.z, 16 * pos.y}}}},
              {{{{180 * x2 * y2 * pos.z, 120 * x2 * xy * pos.z, 60 * x3 * y2 + 16 * pos.z}},
                {{120 * x2 * xy * pos.z, 12 + 30 * x4 * pos.z, 30 * x4 * pos.y}},
                {{60 * x3 * y2 + 16 * pos.z, 30 * x4 * pos.y, 16 * pos.x}}}},
              {{{{60 * x2 * y3, 60 * x3 * y2 + 16 * pos.z, 16 * pos.y}},
                {{60 * x3 * y2 + 16 * pos.z, 30 * x4 * pos.y, 16 * pos.x}},
                {{16 * pos.y, 16 * pos.x, 0}}}}}},
            {{{{{{180 * x2 * y2 * pos.z, 120 * x2 * xy * pos.z, 60 * x3 * y2 + 16 * pos.z}},
                {{120 * x2 * xy * pos.z, 12 + 30 * x4 * pos.z, 30 * x4 * pos.y}},
                {{60 * x3 * y2 + 16 * pos.z, 30 * x4 * pos.y, 16 * pos.x}}}},
              {{{{120 * x2 * xy * pos.z, 12 + 30 * x4 * pos.z, 30 * x4 * pos.y}},
                {{12 + 30 * x4 * pos.z, 0, 6 * x2 * x3}},
                {{30 * x4 * pos.y, 6 * x2 * x3, 0}}}},
              {{{{60 * x3 * y2 + 16 * pos.z, 30 * x4 * pos.y, 16 * pos.x}},
                {{30 * x4 * pos.y, 6 * x2 * x3, 0}},
                {{16 * pos.x, 0, 0}}}}}},
            {{{{{{60 * x2 * y3, 60 * x3 * y2 + 16 * pos.z, 16 * pos.y}},
                {{60 * x3 * y2 + 16 * pos.z, 30 * x4 * pos.y, 16 * pos.x}},
                {{16 * pos.y, 16 * pos.x, 0}}}},
              {{{{60 * x3 * y2 + 16 * pos.z, 30 * x4 * pos.y, 16 * pos.x}},
                {{30 * x4 * pos.y, 6 * x2 * x3, 0}},
                {{16 * pos.x, 0, 0}}}},
              {{{{16 * pos.y, 16 * pos.x, 0}}, {{16 * pos.x, 0, 0}}, {{0, 0, 0}}}}}}}},
          {{{{{{{{{{120 * y2 * yz, 360 * pos.x * y2 * pos.z, 120 * pos.x * y3}},
                  {{360 * pos.x * y2 * pos.z, 360 * x2 * yz, 180 * x2 * y2}},
                  {{120 * pos.x * y3, 180 * x2 * y2, 0}}}},
                {{{{360 * pos.x * y2 * pos.z, 360 * x2 * yz, 180 * x2 * y2}},
                  {{360 * x2 * yz, 120 * x2 * xz, 120 * x2 * xy}},
                  {{180 * x2 * y2, 120 * x2 * xy, 16}}}},
                {{{{120 * pos.x * y3, 180 * x2 * y2, 0}}, {{180 * x2 * y2, 120 * x2 * xy, 16}}, {{0, 16, 0}}}}}},
              {{{{{{360 * pos.x * y2 * pos.z, 360 * x2 * yz, 180 * x2 * y2}},
                  {{360 * x2 * yz, 120 * x2 * xz, 120 * x2 * xy}},
                  {{180 * x2 * y2, 120 * x2 * xy, 16}}}},
                {{{{360 * x2 * yz, 120 * x2 * xz, 120 * x2 * xy}},
                  {{120 * x2 * xz, 0, 30 * x4}},
                  {{120 * x2 * xy, 30 * x4, 0}}}},
                {{{{180 * x2 * y2, 120 * x2 * xy, 16}}, {{120 * x2 * xy, 30 * x4, 0}}, {{16, 0, 0}}}}}},
              {{{{{{120 * pos.x * y3, 180 * x2 * y2, 0}}, {{180 * x2 * y2, 120 * x2 * xy, 16}}, {{0, 16, 0}}}},
                {{{{180 * x2 * y2, 120 * x2 * xy, 16}}, {{120 * x2 * xy, 30 * x4, 0}}, {{16, 0, 0}}}},
                {{{{0, 16, 0}}, {{16, 0, 0}}, {{0, 0, 0}}}}}}}},
            {{{{{{{{360 * pos.x * y2 * pos.z, 360 * x2 * yz, 180 * x2 * y2}},
                  {{360 * x2 * yz, 120 * x2 * xz, 120 * x2 * xy}},
                  {{180 * x2 * y2, 120 * x2 * xy, 16}}}},
                {{{{360 * x2 * yz, 120 * x2 * xz, 120 * x2 * xy}},
                  {{120 * x2 * xz, 0, 30 * x4}},
                  {{120 * x2 * xy, 30 * x4, 0}}}},
                {{{{180 * x2 * y2, 120 * x2 * xy, 16}}, {{120 * x2 * xy, 30 * x4, 0}}, {{16, 0, 0}}}}}},
              {{{{{{360 * x2 * yz, 120 * x2 * xz, 120 * x2 * xy}},
                  {{120 * x2 * xz, 0, 30 * x4}},
                  {{120 * x2 * xy, 30 * x4, 0}}}},
                {{{{120 * x2 * xz, 0, 30 * x4}}, {{0, 0, 0}}, {{30 * x4, 0, 0}}}},
                {{{{120 * x2 * xy, 30 * x4, 0}}, {{30 * x4, 0, 0}}, {{0, 0, 0}}}}}},
              {{{{{{180 * x2 * y2, 120 * x2 * xy, 16}}, {{120 * x2 * xy, 30 * x4, 0}}, {{16, 0, 0}}}},
                {{{{120 * x2 * xy, 30 * x4, 0}}, {{30 * x4, 0, 0}}, {{0, 0, 0}}}},
                {{{{16, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}},
            {{{{{{{{120 * pos.x * y3, 180 * x2 * y2, 0}}, {{180 * x2 * y2, 120 * x2 * xy, 16}}, {{0, 16, 0}}}},
                {{{{180 * x2 * y2, 120 * x2 * xy, 16}}, {{120 * x2 * xy, 30 * x4, 0}}, {{16, 0, 0}}}},
                {{{{0, 16, 0}}, {{16, 0, 0}}, {{0, 0, 0}}}}}},
              {{{{{{180 * x2 * y2, 120 * x2 * xy, 16}}, {{120 * x2 * xy, 30 * x4, 0}}, {{16, 0, 0}}}},
                {{{{120 * x2 * xy, 30 * x4, 0}}, {{30 * x4, 0, 0}}, {{0, 0, 0}}}},
                {{{{16, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
              {{{{{{0, 16, 0}}, {{16, 0, 0}}, {{0, 0, 0}}}},
                {{{{16, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}}}},
          {{{{{{{{{{{{0, 360 * y2 * pos.z, 120 * y3}},
                    {{360 * y2 * pos.z, 720 * xy * pos.z, 360 * pos.x * y2}},
                    {{120 * y3, 360 * pos.x * y2, 0}}}},
                  {{{{360 * y2 * pos.z, 720 * xy * pos.z, 360 * pos.x * y2}},
                    {{720 * xy * pos.z, 360 * x2 * pos.z, 360 * x2 * pos.y}},
                    {{360 * pos.x * y2, 360 * x2 * pos.y, 0}}}},
                  {{{{120 * y3, 360 * pos.x * y2, 0}}, {{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{0, 0, 0}}}}}},
                {{{{{{360 * y2 * pos.z, 720 * xy * pos.z, 360 * pos.x * y2}},
                    {{720 * xy * pos.z, 360 * x2 * pos.z, 360 * x2 * pos.y}},
                    {{360 * pos.x * y2, 360 * x2 * pos.y, 0}}}},
                  {{{{720 * xy * pos.z, 360 * x2 * pos.z, 360 * x2 * pos.y}},
                    {{360 * x2 * pos.z, 0, 120 * x3}},
                    {{360 * x2 * pos.y, 120 * x3, 0}}}},
                  {{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}}}},
                {{{{{{120 * y3, 360 * pos.x * y2, 0}}, {{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{0, 0, 0}}}},
                  {{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}},
              {{{{{{{{360 * y2 * pos.z, 720 * xy * pos.z, 360 * pos.x * y2}},
                    {{720 * xy * pos.z, 360 * x2 * pos.z, 360 * x2 * pos.y}},
                    {{360 * pos.x * y2, 360 * x2 * pos.y, 0}}}},
                  {{{{720 * xy * pos.z, 360 * x2 * pos.z, 360 * x2 * pos.y}},
                    {{360 * x2 * pos.z, 0, 120 * x3}},
                    {{360 * x2 * pos.y, 120 * x3, 0}}}},
                  {{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}}}},
                {{{{{{720 * xy * pos.z, 360 * x2 * pos.z, 360 * x2 * pos.y}},
                    {{360 * x2 * pos.z, 0, 120 * x3}},
                    {{360 * x2 * pos.y, 120 * x3, 0}}}},
                  {{{{360 * x2 * pos.z, 0, 120 * x3}}, {{0, 0, 0}}, {{120 * x3, 0, 0}}}},
                  {{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}},
                  {{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}},
              {{{{{{{{120 * y3, 360 * pos.x * y2, 0}}, {{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{0, 0, 0}}}},
                  {{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}},
                  {{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}}}},
            {{{{{{{{{{360 * y2 * pos.z, 720 * xy * pos.z, 360 * pos.x * y2}},
                    {{720 * xy * pos.z, 360 * x2 * pos.z, 360 * x2 * pos.y}},
                    {{360 * pos.x * y2, 360 * x2 * pos.y, 0}}}},
                  {{{{720 * xy * pos.z, 360 * x2 * pos.z, 360 * x2 * pos.y}},
                    {{360 * x2 * pos.z, 0, 120 * x3}},
                    {{360 * x2 * pos.y, 120 * x3, 0}}}},
                  {{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}}}},
                {{{{{{720 * xy * pos.z, 360 * x2 * pos.z, 360 * x2 * pos.y}},
                    {{360 * x2 * pos.z, 0, 120 * x3}},
                    {{360 * x2 * pos.y, 120 * x3, 0}}}},
                  {{{{360 * x2 * pos.z, 0, 120 * x3}}, {{0, 0, 0}}, {{120 * x3, 0, 0}}}},
                  {{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}},
                  {{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}},
              {{{{{{{{720 * xy * pos.z, 360 * x2 * pos.z, 360 * x2 * pos.y}},
                    {{360 * x2 * pos.z, 0, 120 * x3}},
                    {{360 * x2 * pos.y, 120 * x3, 0}}}},
                  {{{{360 * x2 * pos.z, 0, 120 * x3}}, {{0, 0, 0}}, {{120 * x3, 0, 0}}}},
                  {{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{360 * x2 * pos.z, 0, 120 * x3}}, {{0, 0, 0}}, {{120 * x3, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{120 * x3, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}},
                  {{{{120 * x3, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}},
              {{{{{{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}},
                  {{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}},
                  {{{{120 * x3, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}}}},
            {{{{{{{{{{120 * y3, 360 * pos.x * y2, 0}}, {{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{0, 0, 0}}}},
                  {{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}},
                  {{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}},
              {{{{{{{{360 * pos.x * y2, 360 * x2 * pos.y, 0}}, {{360 * x2 * pos.y, 120 * x3, 0}}, {{0, 0, 0}}}},
                  {{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{360 * x2 * pos.y, 120 * x3, 0}}, {{120 * x3, 0, 0}}, {{0, 0, 0}}}},
                  {{{{120 * x3, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}},
              {{{{{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}}}}}}}}}};
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
Potentials::potentialTriquinticSixthOrderPolynomialTestFunction(const double3 &pos)
{
  // Factor out common terms
  const double x2 = pos.x * pos.x;
  const double y2 = pos.y * pos.y;
  const double z2 = pos.z * pos.z;
  const double xy = pos.x * pos.y;
  const double xz = pos.x * pos.z;
  const double yz = pos.y * pos.z;
  const double x3 = pos.x * x2;
  const double y3 = pos.y * y2;
  const double z3 = pos.z * z2;
  const double r2 = x2 + y2 + z2;
  const double r4 = r2 * r2;
  const double r6 = r4 * r2;

  return {r6,
          {{6 * pos.x * r4, 6 * pos.y * r4, 6 * pos.z * r4}},
          {{{{24 * x2 * r2 + 6 * r2 * r2, 24 * xy * r2, 24 * xz * r2}},
            {{24 * xy * r2, 24 * y2 * r2 + 6 * r2 * r2, 24 * yz * r2}},
            {{24 * xz * r2, 24 * yz * r2, 24 * z2 * r2 + 6 * r2 * r2}}}},
          {{{{{{48 * x3 + 72 * pos.x * r2, 48 * x2 * pos.y + 24 * pos.y * r2, 48 * x2 * pos.z + 24 * pos.z * r2}},
              {{48 * x2 * pos.y + 24 * pos.y * r2, 48 * pos.x * y2 + 24 * pos.x * r2, 48 * xy * pos.z}},
              {{48 * x2 * pos.z + 24 * pos.z * r2, 48 * xy * pos.z, 48 * pos.x * z2 + 24 * pos.x * r2}}}},
            {{{{48 * x2 * pos.y + 24 * pos.y * r2, 48 * pos.x * y2 + 24 * pos.x * r2, 48 * xy * pos.z}},
              {{48 * pos.x * y2 + 24 * pos.x * r2, 48 * y3 + 72 * pos.y * r2, 48 * y2 * pos.z + 24 * pos.z * r2}},
              {{48 * xy * pos.z, 48 * y2 * pos.z + 24 * pos.z * r2, 48 * pos.y * z2 + 24 * pos.y * r2}}}},
            {{{{48 * x2 * pos.z + 24 * pos.z * r2, 48 * xy * pos.z, 48 * pos.x * z2 + 24 * pos.x * r2}},
              {{48 * xy * pos.z, 48 * y2 * pos.z + 24 * pos.z * r2, 48 * pos.y * z2 + 24 * pos.y * r2}},
              {{48 * pos.x * z2 + 24 * pos.x * r2, 48 * pos.y * z2 + 24 * pos.y * r2, 48 * z3 + 72 * pos.z * r2}}}}}},
          {{{{{{{{288 * x2 + 72 * r2, 144 * xy, 144 * xz}},
                {{144 * xy, 48 * x2 + 48 * y2 + 24 * r2, 48 * yz}},
                {{144 * xz, 48 * yz, 48 * x2 + 48 * z2 + 24 * r2}}}},
              {{{{144 * xy, 48 * x2 + 48 * y2 + 24 * r2, 48 * yz}},
                {{48 * x2 + 48 * y2 + 24 * r2, 144 * xy, 48 * xz}},
                {{48 * yz, 48 * xz, 48 * xy}}}},
              {{{{144 * xz, 48 * yz, 48 * x2 + 48 * z2 + 24 * r2}},
                {{48 * yz, 48 * xz, 48 * xy}},
                {{48 * x2 + 48 * z2 + 24 * r2, 48 * xy, 144 * xz}}}}}},
            {{{{{{144 * xy, 48 * x2 + 48 * y2 + 24 * r2, 48 * yz}},
                {{48 * x2 + 48 * y2 + 24 * r2, 144 * xy, 48 * xz}},
                {{48 * yz, 48 * xz, 48 * xy}}}},
              {{{{48 * x2 + 48 * y2 + 24 * r2, 144 * xy, 48 * xz}},
                {{144 * xy, 288 * y2 + 72 * r2, 144 * yz}},
                {{48 * xz, 144 * yz, 48 * y2 + 48 * z2 + 24 * r2}}}},
              {{{{48 * yz, 48 * xz, 48 * xy}},
                {{48 * xz, 144 * yz, 48 * y2 + 48 * z2 + 24 * r2}},
                {{48 * xy, 48 * y2 + 48 * z2 + 24 * r2, 144 * yz}}}}}},
            {{{{{{144 * xz, 48 * yz, 48 * x2 + 48 * z2 + 24 * r2}},
                {{48 * yz, 48 * xz, 48 * xy}},
                {{48 * x2 + 48 * z2 + 24 * r2, 48 * xy, 144 * xz}}}},
              {{{{48 * yz, 48 * xz, 48 * xy}},
                {{48 * xz, 144 * yz, 48 * y2 + 48 * z2 + 24 * r2}},
                {{48 * xy, 48 * y2 + 48 * z2 + 24 * r2, 144 * yz}}}},
              {{{{48 * x2 + 48 * z2 + 24 * r2, 48 * xy, 144 * xz}},
                {{48 * xy, 48 * y2 + 48 * z2 + 24 * r2, 144 * yz}},
                {{144 * xz, 144 * yz, 288 * z2 + 72 * r2}}}}}}}},
          {{{{{{{{{{720 * pos.x, 144 * pos.y, 144 * pos.z}},
                  {{144 * pos.y, 144 * pos.x, 0}},
                  {{144 * pos.z, 0, 144 * pos.x}}}},
                {{{{144 * pos.y, 144 * pos.x, 0}},
                  {{144 * pos.x, 144 * pos.y, 48 * pos.z}},
                  {{0, 48 * pos.z, 48 * pos.y}}}},
                {{{{144 * pos.z, 0, 144 * pos.x}},
                  {{0, 48 * pos.z, 48 * pos.y}},
                  {{144 * pos.x, 48 * pos.y, 144 * pos.z}}}}}},
              {{{{{{144 * pos.y, 144 * pos.x, 0}},
                  {{144 * pos.x, 144 * pos.y, 48 * pos.z}},
                  {{0, 48 * pos.z, 48 * pos.y}}}},
                {{{{144 * pos.x, 144 * pos.y, 48 * pos.z}},
                  {{144 * pos.y, 144 * pos.x, 0}},
                  {{48 * pos.z, 0, 48 * pos.x}}}},
                {{{{0, 48 * pos.z, 48 * pos.y}}, {{48 * pos.z, 0, 48 * pos.x}}, {{48 * pos.y, 48 * pos.x, 0}}}}}},
              {{{{{{144 * pos.z, 0, 144 * pos.x}},
                  {{0, 48 * pos.z, 48 * pos.y}},
                  {{144 * pos.x, 48 * pos.y, 144 * pos.z}}}},
                {{{{0, 48 * pos.z, 48 * pos.y}}, {{48 * pos.z, 0, 48 * pos.x}}, {{48 * pos.y, 48 * pos.x, 0}}}},
                {{{{144 * pos.x, 48 * pos.y, 144 * pos.z}},
                  {{48 * pos.y, 48 * pos.x, 0}},
                  {{144 * pos.z, 0, 144 * pos.x}}}}}}}},
            {{{{{{{{144 * pos.y, 144 * pos.x, 0}},
                  {{144 * pos.x, 144 * pos.y, 48 * pos.z}},
                  {{0, 48 * pos.z, 48 * pos.y}}}},
                {{{{144 * pos.x, 144 * pos.y, 48 * pos.z}},
                  {{144 * pos.y, 144 * pos.x, 0}},
                  {{48 * pos.z, 0, 48 * pos.x}}}},
                {{{{0, 48 * pos.z, 48 * pos.y}}, {{48 * pos.z, 0, 48 * pos.x}}, {{48 * pos.y, 48 * pos.x, 0}}}}}},
              {{{{{{144 * pos.x, 144 * pos.y, 48 * pos.z}},
                  {{144 * pos.y, 144 * pos.x, 0}},
                  {{48 * pos.z, 0, 48 * pos.x}}}},
                {{{{144 * pos.y, 144 * pos.x, 0}},
                  {{144 * pos.x, 720 * pos.y, 144 * pos.z}},
                  {{0, 144 * pos.z, 144 * pos.y}}}},
                {{{{48 * pos.z, 0, 48 * pos.x}},
                  {{0, 144 * pos.z, 144 * pos.y}},
                  {{48 * pos.x, 144 * pos.y, 144 * pos.z}}}}}},
              {{{{{{0, 48 * pos.z, 48 * pos.y}}, {{48 * pos.z, 0, 48 * pos.x}}, {{48 * pos.y, 48 * pos.x, 0}}}},
                {{{{48 * pos.z, 0, 48 * pos.x}},
                  {{0, 144 * pos.z, 144 * pos.y}},
                  {{48 * pos.x, 144 * pos.y, 144 * pos.z}}}},
                {{{{48 * pos.y, 48 * pos.x, 0}},
                  {{48 * pos.x, 144 * pos.y, 144 * pos.z}},
                  {{0, 144 * pos.z, 144 * pos.y}}}}}}}},
            {{{{{{{{144 * pos.z, 0, 144 * pos.x}},
                  {{0, 48 * pos.z, 48 * pos.y}},
                  {{144 * pos.x, 48 * pos.y, 144 * pos.z}}}},
                {{{{0, 48 * pos.z, 48 * pos.y}}, {{48 * pos.z, 0, 48 * pos.x}}, {{48 * pos.y, 48 * pos.x, 0}}}},
                {{{{144 * pos.x, 48 * pos.y, 144 * pos.z}},
                  {{48 * pos.y, 48 * pos.x, 0}},
                  {{144 * pos.z, 0, 144 * pos.x}}}}}},
              {{{{{{0, 48 * pos.z, 48 * pos.y}}, {{48 * pos.z, 0, 48 * pos.x}}, {{48 * pos.y, 48 * pos.x, 0}}}},
                {{{{48 * pos.z, 0, 48 * pos.x}},
                  {{0, 144 * pos.z, 144 * pos.y}},
                  {{48 * pos.x, 144 * pos.y, 144 * pos.z}}}},
                {{{{48 * pos.y, 48 * pos.x, 0}},
                  {{48 * pos.x, 144 * pos.y, 144 * pos.z}},
                  {{0, 144 * pos.z, 144 * pos.y}}}}}},
              {{{{{{144 * pos.x, 48 * pos.y, 144 * pos.z}},
                  {{48 * pos.y, 48 * pos.x, 0}},
                  {{144 * pos.z, 0, 144 * pos.x}}}},
                {{{{48 * pos.y, 48 * pos.x, 0}},
                  {{48 * pos.x, 144 * pos.y, 144 * pos.z}},
                  {{0, 144 * pos.z, 144 * pos.y}}}},
                {{{{144 * pos.z, 0, 144 * pos.x}},
                  {{0, 144 * pos.z, 144 * pos.y}},
                  {{144 * pos.x, 144 * pos.y, 720 * pos.z}}}}}}}}}},
          {{{{{{{{{{{{720, 0, 0}}, {{0, 144, 0}}, {{0, 0, 144}}}},
                  {{{{0, 144, 0}}, {{144, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 144}}, {{0, 0, 0}}, {{144, 0, 0}}}}}},
                {{{{{{0, 144, 0}}, {{144, 0, 0}}, {{0, 0, 0}}}},
                  {{{{144, 0, 0}}, {{0, 144, 0}}, {{0, 0, 48}}}},
                  {{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}}}},
                {{{{{{0, 0, 144}}, {{0, 0, 0}}, {{144, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}},
                  {{{{144, 0, 0}}, {{0, 48, 0}}, {{0, 0, 144}}}}}}}},
              {{{{{{{{0, 144, 0}}, {{144, 0, 0}}, {{0, 0, 0}}}},
                  {{{{144, 0, 0}}, {{0, 144, 0}}, {{0, 0, 48}}}},
                  {{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}}}},
                {{{{{{144, 0, 0}}, {{0, 144, 0}}, {{0, 0, 48}}}},
                  {{{{0, 144, 0}}, {{144, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}}}},
                {{{{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}},
                  {{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}},
                  {{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}}}}}},
              {{{{{{{{0, 0, 144}}, {{0, 0, 0}}, {{144, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}},
                  {{{{144, 0, 0}}, {{0, 48, 0}}, {{0, 0, 144}}}}}},
                {{{{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}},
                  {{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}},
                  {{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{144, 0, 0}}, {{0, 48, 0}}, {{0, 0, 144}}}},
                  {{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 144}}, {{0, 0, 0}}, {{144, 0, 0}}}}}}}}}},
            {{{{{{{{{{0, 144, 0}}, {{144, 0, 0}}, {{0, 0, 0}}}},
                  {{{{144, 0, 0}}, {{0, 144, 0}}, {{0, 0, 48}}}},
                  {{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}}}},
                {{{{{{144, 0, 0}}, {{0, 144, 0}}, {{0, 0, 48}}}},
                  {{{{0, 144, 0}}, {{144, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}}}},
                {{{{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}},
                  {{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}},
                  {{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}}}}}},
              {{{{{{{{144, 0, 0}}, {{0, 144, 0}}, {{0, 0, 48}}}},
                  {{{{0, 144, 0}}, {{144, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}}}},
                {{{{{{0, 144, 0}}, {{144, 0, 0}}, {{0, 0, 0}}}},
                  {{{{144, 0, 0}}, {{0, 720, 0}}, {{0, 0, 144}}}},
                  {{{{0, 0, 0}}, {{0, 0, 144}}, {{0, 144, 0}}}}}},
                {{{{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 144}}, {{0, 144, 0}}}},
                  {{{{48, 0, 0}}, {{0, 144, 0}}, {{0, 0, 144}}}}}}}},
              {{{{{{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}},
                  {{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}},
                  {{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 144}}, {{0, 144, 0}}}},
                  {{{{48, 0, 0}}, {{0, 144, 0}}, {{0, 0, 144}}}}}},
                {{{{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}},
                  {{{{48, 0, 0}}, {{0, 144, 0}}, {{0, 0, 144}}}},
                  {{{{0, 0, 0}}, {{0, 0, 144}}, {{0, 144, 0}}}}}}}}}},
            {{{{{{{{{{0, 0, 144}}, {{0, 0, 0}}, {{144, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}},
                  {{{{144, 0, 0}}, {{0, 48, 0}}, {{0, 0, 144}}}}}},
                {{{{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}},
                  {{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}},
                  {{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{144, 0, 0}}, {{0, 48, 0}}, {{0, 0, 144}}}},
                  {{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 144}}, {{0, 0, 0}}, {{144, 0, 0}}}}}}}},
              {{{{{{{{0, 0, 0}}, {{0, 0, 48}}, {{0, 48, 0}}}},
                  {{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}},
                  {{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}}}},
                {{{{{{0, 0, 48}}, {{0, 0, 0}}, {{48, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 144}}, {{0, 144, 0}}}},
                  {{{{48, 0, 0}}, {{0, 144, 0}}, {{0, 0, 144}}}}}},
                {{{{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}},
                  {{{{48, 0, 0}}, {{0, 144, 0}}, {{0, 0, 144}}}},
                  {{{{0, 0, 0}}, {{0, 0, 144}}, {{0, 144, 0}}}}}}}},
              {{{{{{{{144, 0, 0}}, {{0, 48, 0}}, {{0, 0, 144}}}},
                  {{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}},
                  {{{{0, 0, 144}}, {{0, 0, 0}}, {{144, 0, 0}}}}}},
                {{{{{{0, 48, 0}}, {{48, 0, 0}}, {{0, 0, 0}}}},
                  {{{{48, 0, 0}}, {{0, 144, 0}}, {{0, 0, 144}}}},
                  {{{{0, 0, 0}}, {{0, 0, 144}}, {{0, 144, 0}}}}}},
                {{{{{{0, 0, 144}}, {{0, 0, 0}}, {{144, 0, 0}}}},
                  {{{{0, 0, 0}}, {{0, 0, 144}}, {{0, 144, 0}}}},
                  {{{{144, 0, 0}}, {{0, 144, 0}}, {{0, 0, 720}}}}}}}}}}}}};
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
Potentials::potentialTriquinticExponentialNonPolynomialTestFunction(const double3 &pos)
{
  // Factor out common terms
  const double x2 = pos.x * pos.x;
  const double y2 = pos.y * pos.y;
  const double z2 = pos.z * pos.z;
  const double xy = pos.x * pos.y;
  const double xz = pos.x * pos.z;
  const double yz = pos.y * pos.z;
  const double x3 = pos.x * x2;
  const double y3 = pos.y * y2;
  const double z3 = pos.z * z2;
  const double x4 = x2 * x2;
  const double y4 = y2 * y2;
  const double z4 = z2 * z2;
  const double x6 = x4 * x2;
  const double y6 = y4 * y2;
  const double z6 = z4 * z2;
  const double r2 = x2 + y2 + z2;
  double exp_term = std::exp(-r2);

  const double exp_x2 = exp_term * x2;
  const double exp_y2 = exp_term * y2;
  const double exp_z2 = exp_term * z2;
  const double exp_xy = exp_term * xy;
  const double exp_xz = exp_term * xz;
  const double exp_yz = exp_term * yz;
  const double exp_r2 = exp_term * r2;

  const double exp_x2_r2 = exp_x2 * r2;
  const double exp_y2_r2 = exp_y2 * r2;
  const double exp_z2_r2 = exp_z2 * r2;
  const double exp_xy_r2 = exp_xy * r2;
  const double exp_xz_r2 = exp_xz * r2;
  const double exp_yz_r2 = exp_yz * r2;

  const double exp_x3 = exp_term * x3;
  const double exp_y3 = exp_term * y3;
  const double exp_z3 = exp_term * z3;

  const double exp_x2_xy = exp_x2 * xy;
  const double exp_x2_xz = exp_x2 * xz;
  const double exp_x2_yz = exp_x2 * yz;
  const double exp_y2_yz = exp_y2 * yz;

  const double exp_term_x4 = exp_term * x4;
  const double exp_term_y4 = exp_term * y4;
  const double exp_term_z4 = exp_term * z4;

  const double exp_term_pos_x = exp_term * pos.x;
  const double exp_term_pos_y = exp_term * pos.y;
  const double exp_term_pos_z = exp_term * pos.z;

  const double exp_xy_z2 = exp_xy * z2;
  const double exp_xy_z4 = exp_xy * z4;
  const double exp_x2_xy_z2 = exp_x2_xy * z2;

  const double exp_x3_y2 = exp_x3 * y2;
  const double exp_x3_z2 = exp_x3 * z2;
  const double exp_y3_z2 = exp_y3 * z2;
  const double exp_y3_z3 = exp_y3 * z3;

  const double exp_term_x4_y2 = exp_term_x4 * y2;
  const double exp_term_x4_z2 = exp_term_x4 * z2;
  const double exp_term_y4_z2 = exp_term_y4 * z2;

  const double exp_term_pos_x_y2 = exp_term_pos_x * y2;
  const double exp_term_pos_x_z2 = exp_term_pos_x * z2;
  const double exp_term_pos_y_z2 = exp_term_pos_y * z2;
  const double exp_term_pos_y_z3 = exp_term_pos_y * z3;

  return {
      r2 * exp_term,
      {{2 * exp_term_pos_x - 2 * exp_term_pos_x * r2, 2 * exp_term_pos_y - 2 * exp_term_pos_y * r2,
        2 * exp_term_pos_z - 2 * exp_term_pos_z * r2}},
      {{{{2 * exp_term - 8 * exp_x2 - 2 * exp_r2 + 4 * exp_x2_r2, -8 * exp_xy + 4 * exp_xy_r2,
          -8 * exp_xz + 4 * exp_xz_r2}},
        {{-8 * exp_xy + 4 * exp_xy_r2, 2 * exp_term - 8 * exp_y2 - 2 * exp_r2 + 4 * exp_y2_r2,
          -8 * exp_yz + 4 * exp_yz_r2}},
        {{-8 * exp_xz + 4 * exp_xz_r2, -8 * exp_yz + 4 * exp_yz_r2,
          2 * exp_term - 8 * exp_z2 - 2 * exp_r2 + 4 * exp_z2_r2}}}},
      {{{{{{-24 * exp_term_pos_x + 24 * exp_x3 + 12 * exp_term_pos_x * r2 - 8 * exp_x3 * r2,
            -8 * exp_term_pos_y + 24 * exp_x2 * pos.y + 4 * exp_term_pos_y * r2 - 8 * exp_x2 * pos.y * r2,
            -8 * exp_term_pos_z + 24 * exp_x2 * pos.z + 4 * exp_term_pos_z * r2 - 8 * exp_x2 * pos.z * r2}},
          {{-8 * exp_term_pos_y + 24 * exp_x2 * pos.y + 4 * exp_term_pos_y * r2 - 8 * exp_x2 * pos.y * r2,
            -8 * exp_term_pos_x + 24 * exp_term_pos_x_y2 + 4 * exp_term_pos_x * r2 - 8 * exp_term_pos_x_y2 * r2,
            24 * exp_xy * pos.z - 8 * exp_xy * pos.z * r2}},
          {{-8 * exp_term_pos_z + 24 * exp_x2 * pos.z + 4 * exp_term_pos_z * r2 - 8 * exp_x2 * pos.z * r2,
            24 * exp_xy * pos.z - 8 * exp_xy * pos.z * r2,
            -8 * exp_term_pos_x + 24 * exp_term_pos_x_z2 + 4 * exp_term_pos_x * r2 - 8 * exp_term_pos_x_z2 * r2}}}},
        {{{{-8 * exp_term_pos_y + 24 * exp_x2 * pos.y + 4 * exp_term_pos_y * r2 - 8 * exp_x2 * pos.y * r2,
            -8 * exp_term_pos_x + 24 * exp_term_pos_x_y2 + 4 * exp_term_pos_x * r2 - 8 * exp_term_pos_x_y2 * r2,
            24 * exp_xy * pos.z - 8 * exp_xy * pos.z * r2}},
          {{-8 * exp_term_pos_x + 24 * exp_term_pos_x_y2 + 4 * exp_term_pos_x * r2 - 8 * exp_term_pos_x_y2 * r2,
            -24 * exp_term_pos_y + 24 * exp_y3 + 12 * exp_term_pos_y * r2 - 8 * exp_y3 * r2,
            -8 * exp_term_pos_z + 24 * exp_y2 * pos.z + 4 * exp_term_pos_z * r2 - 8 * exp_y2 * pos.z * r2}},
          {{24 * exp_xy * pos.z - 8 * exp_xy * pos.z * r2,
            -8 * exp_term_pos_z + 24 * exp_y2 * pos.z + 4 * exp_term_pos_z * r2 - 8 * exp_y2 * pos.z * r2,
            -8 * exp_term_pos_y + 24 * exp_term_pos_y_z2 + 4 * exp_term_pos_y * r2 - 8 * exp_term_pos_y_z2 * r2}}}},
        {{{{-8 * exp_term_pos_z + 24 * exp_x2 * pos.z + 4 * exp_term_pos_z * r2 - 8 * exp_x2 * pos.z * r2,
            24 * exp_xy * pos.z - 8 * exp_xy * pos.z * r2,
            -8 * exp_term_pos_x + 24 * exp_term_pos_x_z2 + 4 * exp_term_pos_x * r2 - 8 * exp_term_pos_x_z2 * r2}},
          {{24 * exp_xy * pos.z - 8 * exp_xy * pos.z * r2,
            -8 * exp_term_pos_z + 24 * exp_y2 * pos.z + 4 * exp_term_pos_z * r2 - 8 * exp_y2 * pos.z * r2,
            -8 * exp_term_pos_y + 24 * exp_term_pos_y_z2 + 4 * exp_term_pos_y * r2 - 8 * exp_term_pos_y_z2 * r2}},
          {{-8 * exp_term_pos_x + 24 * exp_term_pos_x_z2 + 4 * exp_term_pos_x * r2 - 8 * exp_term_pos_x_z2 * r2,
            -8 * exp_term_pos_y + 24 * exp_term_pos_y_z2 + 4 * exp_term_pos_y * r2 - 8 * exp_term_pos_y_z2 * r2,
            -24 * exp_term_pos_z + 24 * exp_z3 + 12 * exp_term_pos_z * r2 - 8 * exp_z3 * r2}}}}}},
      {{{{{{{{-24 * exp_term + 144 * exp_x2 - 64 * exp_term_x4 + 12 * exp_r2 - 48 * exp_x2_r2 + 16 * exp_term_x4 * r2,
              72 * exp_xy - 64 * exp_x2_xy - 24 * exp_xy_r2 + 16 * exp_x2_xy * r2,
              72 * exp_xz - 64 * exp_x2_xz - 24 * exp_xz_r2 + 16 * exp_x2_xz * r2}},
            {{72 * exp_xy - 64 * exp_x2_xy - 24 * exp_xy_r2 + 16 * exp_x2_xy * r2,
              -8 * exp_term + 24 * exp_x2 + 24 * exp_y2 - 64 * exp_x2 * y2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_y2_r2 + 16 * exp_x2 * y2 * r2,
              24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2}},
            {{72 * exp_xz - 64 * exp_x2_xz - 24 * exp_xz_r2 + 16 * exp_x2_xz * r2,
              24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2,
              -8 * exp_term + 24 * exp_x2 + 24 * exp_z2 - 64 * exp_x2 * z2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_z2_r2 + 16 * exp_x2 * z2 * r2}}}},
          {{{{72 * exp_xy - 64 * exp_x2_xy - 24 * exp_xy_r2 + 16 * exp_x2_xy * r2,
              -8 * exp_term + 24 * exp_x2 + 24 * exp_y2 - 64 * exp_x2 * y2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_y2_r2 + 16 * exp_x2 * y2 * r2,
              24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2}},
            {{-8 * exp_term + 24 * exp_x2 + 24 * exp_y2 - 64 * exp_x2 * y2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_y2_r2 + 16 * exp_x2 * y2 * r2,
              72 * exp_xy - 64 * exp_term_pos_x * y3 - 24 * exp_xy_r2 + 16 * exp_term_pos_x * y3 * r2,
              24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2}},
            {{24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2,
              24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2,
              24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2}}}},
          {{{{72 * exp_xz - 64 * exp_x2_xz - 24 * exp_xz_r2 + 16 * exp_x2_xz * r2,
              24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2,
              -8 * exp_term + 24 * exp_x2 + 24 * exp_z2 - 64 * exp_x2 * z2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_z2_r2 + 16 * exp_x2 * z2 * r2}},
            {{24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2,
              24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2,
              24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2}},
            {{-8 * exp_term + 24 * exp_x2 + 24 * exp_z2 - 64 * exp_x2 * z2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_z2_r2 + 16 * exp_x2 * z2 * r2,
              24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2,
              72 * exp_xz - 64 * exp_term_pos_x * z3 - 24 * exp_xz_r2 + 16 * exp_term_pos_x * z3 * r2}}}}}},
        {{{{{{72 * exp_xy - 64 * exp_x2_xy - 24 * exp_xy_r2 + 16 * exp_x2_xy * r2,
              -8 * exp_term + 24 * exp_x2 + 24 * exp_y2 - 64 * exp_x2 * y2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_y2_r2 + 16 * exp_x2 * y2 * r2,
              24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2}},
            {{-8 * exp_term + 24 * exp_x2 + 24 * exp_y2 - 64 * exp_x2 * y2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_y2_r2 + 16 * exp_x2 * y2 * r2,
              72 * exp_xy - 64 * exp_term_pos_x * y3 - 24 * exp_xy_r2 + 16 * exp_term_pos_x * y3 * r2,
              24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2}},
            {{24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2,
              24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2,
              24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2}}}},
          {{{{-8 * exp_term + 24 * exp_x2 + 24 * exp_y2 - 64 * exp_x2 * y2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_y2_r2 + 16 * exp_x2 * y2 * r2,
              72 * exp_xy - 64 * exp_term_pos_x * y3 - 24 * exp_xy_r2 + 16 * exp_term_pos_x * y3 * r2,
              24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2}},
            {{72 * exp_xy - 64 * exp_term_pos_x * y3 - 24 * exp_xy_r2 + 16 * exp_term_pos_x * y3 * r2,
              -24 * exp_term + 144 * exp_y2 - 64 * exp_term_y4 + 12 * exp_r2 - 48 * exp_y2_r2 + 16 * exp_term_y4 * r2,
              72 * exp_yz - 64 * exp_y2_yz - 24 * exp_yz_r2 + 16 * exp_y2_yz * r2}},
            {{24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2,
              72 * exp_yz - 64 * exp_y2_yz - 24 * exp_yz_r2 + 16 * exp_y2_yz * r2,
              -8 * exp_term + 24 * exp_y2 + 24 * exp_z2 - 64 * exp_y2 * z2 + 4 * exp_r2 - 8 * exp_y2_r2 -
                  8 * exp_z2_r2 + 16 * exp_y2 * z2 * r2}}}},
          {{{{24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2,
              24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2,
              24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2}},
            {{24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2,
              72 * exp_yz - 64 * exp_y2_yz - 24 * exp_yz_r2 + 16 * exp_y2_yz * r2,
              -8 * exp_term + 24 * exp_y2 + 24 * exp_z2 - 64 * exp_y2 * z2 + 4 * exp_r2 - 8 * exp_y2_r2 -
                  8 * exp_z2_r2 + 16 * exp_y2 * z2 * r2}},
            {{24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2,
              -8 * exp_term + 24 * exp_y2 + 24 * exp_z2 - 64 * exp_y2 * z2 + 4 * exp_r2 - 8 * exp_y2_r2 -
                  8 * exp_z2_r2 + 16 * exp_y2 * z2 * r2,
              72 * exp_yz - 64 * exp_term_pos_y_z3 - 24 * exp_yz_r2 + 16 * exp_term_pos_y_z3 * r2}}}}}},
        {{{{{{72 * exp_xz - 64 * exp_x2_xz - 24 * exp_xz_r2 + 16 * exp_x2_xz * r2,
              24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2,
              -8 * exp_term + 24 * exp_x2 + 24 * exp_z2 - 64 * exp_x2 * z2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_z2_r2 + 16 * exp_x2 * z2 * r2}},
            {{24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2,
              24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2,
              24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2}},
            {{-8 * exp_term + 24 * exp_x2 + 24 * exp_z2 - 64 * exp_x2 * z2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_z2_r2 + 16 * exp_x2 * z2 * r2,
              24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2,
              72 * exp_xz - 64 * exp_term_pos_x * z3 - 24 * exp_xz_r2 + 16 * exp_term_pos_x * z3 * r2}}}},
          {{{{24 * exp_yz - 64 * exp_x2_yz - 8 * exp_yz_r2 + 16 * exp_x2_yz * r2,
              24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2,
              24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2}},
            {{24 * exp_xz - 64 * exp_term_pos_x_y2 * pos.z - 8 * exp_xz_r2 + 16 * exp_term_pos_x_y2 * pos.z * r2,
              72 * exp_yz - 64 * exp_y2_yz - 24 * exp_yz_r2 + 16 * exp_y2_yz * r2,
              -8 * exp_term + 24 * exp_y2 + 24 * exp_z2 - 64 * exp_y2 * z2 + 4 * exp_r2 - 8 * exp_y2_r2 -
                  8 * exp_z2_r2 + 16 * exp_y2 * z2 * r2}},
            {{24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2,
              -8 * exp_term + 24 * exp_y2 + 24 * exp_z2 - 64 * exp_y2 * z2 + 4 * exp_r2 - 8 * exp_y2_r2 -
                  8 * exp_z2_r2 + 16 * exp_y2 * z2 * r2,
              72 * exp_yz - 64 * exp_term_pos_y_z3 - 24 * exp_yz_r2 + 16 * exp_term_pos_y_z3 * r2}}}},
          {{{{-8 * exp_term + 24 * exp_x2 + 24 * exp_z2 - 64 * exp_x2 * z2 + 4 * exp_r2 - 8 * exp_x2_r2 -
                  8 * exp_z2_r2 + 16 * exp_x2 * z2 * r2,
              24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2,
              72 * exp_xz - 64 * exp_term_pos_x * z3 - 24 * exp_xz_r2 + 16 * exp_term_pos_x * z3 * r2}},
            {{24 * exp_xy - 64 * exp_xy_z2 - 8 * exp_xy_r2 + 16 * exp_xy_z2 * r2,
              -8 * exp_term + 24 * exp_y2 + 24 * exp_z2 - 64 * exp_y2 * z2 + 4 * exp_r2 - 8 * exp_y2_r2 -
                  8 * exp_z2_r2 + 16 * exp_y2 * z2 * r2,
              72 * exp_yz - 64 * exp_term_pos_y_z3 - 24 * exp_yz_r2 + 16 * exp_term_pos_y_z3 * r2}},
            {{72 * exp_xz - 64 * exp_term_pos_x * z3 - 24 * exp_xz_r2 + 16 * exp_term_pos_x * z3 * r2,
              72 * exp_yz - 64 * exp_term_pos_y_z3 - 24 * exp_yz_r2 + 16 * exp_term_pos_y_z3 * r2,
              -24 * exp_term + 144 * exp_z2 - 64 * exp_term_z4 + 12 * exp_r2 - 48 * exp_z2_r2 +
                  16 * exp_term_z4 * r2}}}}}}}},
      {{{{{{{{{{360 * exp_term_pos_x - 640 * exp_x3 + 160 * exp_x2 * x3 - 120 * exp_term_pos_x * r2 +
                    160 * exp_x3 * r2 - 32 * exp_x2 * x3 * r2,
                72 * exp_term_pos_y - 384 * exp_x2 * pos.y + 160 * exp_term_x4 * pos.y - 24 * exp_term_pos_y * r2 +
                    96 * exp_x2 * pos.y * r2 - 32 * exp_term_x4 * pos.y * r2,
                72 * exp_term_pos_z - 384 * exp_x2 * pos.z + 160 * exp_term_x4 * pos.z - 24 * exp_term_pos_z * r2 +
                    96 * exp_x2 * pos.z * r2 - 32 * exp_term_x4 * pos.z * r2}},
              {{72 * exp_term_pos_y - 384 * exp_x2 * pos.y + 160 * exp_term_x4 * pos.y - 24 * exp_term_pos_y * r2 +
                    96 * exp_x2 * pos.y * r2 - 32 * exp_term_x4 * pos.y * r2,
                72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_y2 + 160 * exp_x3_y2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_y2 * r2 - 32 * exp_x3_y2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2}},
              {{72 * exp_term_pos_z - 384 * exp_x2 * pos.z + 160 * exp_term_x4 * pos.z - 24 * exp_term_pos_z * r2 +
                    96 * exp_x2 * pos.z * r2 - 32 * exp_term_x4 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_z2 + 160 * exp_x3_z2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_z2 * r2 - 32 * exp_x3_z2 * r2}}}},
            {{{{72 * exp_term_pos_y - 384 * exp_x2 * pos.y + 160 * exp_term_x4 * pos.y - 24 * exp_term_pos_y * r2 +
                    96 * exp_x2 * pos.y * r2 - 32 * exp_term_x4 * pos.y * r2,
                72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_y2 + 160 * exp_x3_y2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_y2 * r2 - 32 * exp_x3_y2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2}},
              {{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_y2 + 160 * exp_x3_y2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_y2 * r2 - 32 * exp_x3_y2 * r2,
                72 * exp_term_pos_y - 192 * exp_x2 * pos.y - 64 * exp_y3 + 160 * exp_x2 * y3 -
                    24 * exp_term_pos_y * r2 + 48 * exp_x2 * pos.y * r2 + 16 * exp_y3 * r2 - 32 * exp_x2 * y3 * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}}}},
            {{{{72 * exp_term_pos_z - 384 * exp_x2 * pos.z + 160 * exp_term_x4 * pos.z - 24 * exp_term_pos_z * r2 +
                    96 * exp_x2 * pos.z * r2 - 32 * exp_term_x4 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_z2 + 160 * exp_x3_z2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_z2 * r2 - 32 * exp_x3_z2 * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}},
              {{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_z2 + 160 * exp_x3_z2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_z2 * r2 - 32 * exp_x3_z2 * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                72 * exp_term_pos_z - 192 * exp_x2 * pos.z - 64 * exp_z3 + 160 * exp_x2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_x2 * pos.z * r2 + 16 * exp_z3 * r2 -
                    32 * exp_x2 * z3 * r2}}}}}},
          {{{{{{72 * exp_term_pos_y - 384 * exp_x2 * pos.y + 160 * exp_term_x4 * pos.y - 24 * exp_term_pos_y * r2 +
                    96 * exp_x2 * pos.y * r2 - 32 * exp_term_x4 * pos.y * r2,
                72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_y2 + 160 * exp_x3_y2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_y2 * r2 - 32 * exp_x3_y2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2}},
              {{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_y2 + 160 * exp_x3_y2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_y2 * r2 - 32 * exp_x3_y2 * r2,
                72 * exp_term_pos_y - 192 * exp_x2 * pos.y - 64 * exp_y3 + 160 * exp_x2 * y3 -
                    24 * exp_term_pos_y * r2 + 48 * exp_x2 * pos.y * r2 + 16 * exp_y3 * r2 - 32 * exp_x2 * y3 * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}}}},
            {{{{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_y2 + 160 * exp_x3_y2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_y2 * r2 - 32 * exp_x3_y2 * r2,
                72 * exp_term_pos_y - 192 * exp_x2 * pos.y - 64 * exp_y3 + 160 * exp_x2 * y3 -
                    24 * exp_term_pos_y * r2 + 48 * exp_x2 * pos.y * r2 + 16 * exp_y3 * r2 - 32 * exp_x2 * y3 * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2}},
              {{72 * exp_term_pos_y - 192 * exp_x2 * pos.y - 64 * exp_y3 + 160 * exp_x2 * y3 -
                    24 * exp_term_pos_y * r2 + 48 * exp_x2 * pos.y * r2 + 16 * exp_y3 * r2 - 32 * exp_x2 * y3 * r2,
                72 * exp_term_pos_x - 384 * exp_term_pos_x_y2 + 160 * exp_term_pos_x * y4 - 24 * exp_term_pos_x * r2 +
                    96 * exp_term_pos_x_y2 * r2 - 32 * exp_term_pos_x * y4 * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2}},
              {{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}}}},
            {{{{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}},
              {{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}},
              {{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}}}}}},
          {{{{{{72 * exp_term_pos_z - 384 * exp_x2 * pos.z + 160 * exp_term_x4 * pos.z - 24 * exp_term_pos_z * r2 +
                    96 * exp_x2 * pos.z * r2 - 32 * exp_term_x4 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_z2 + 160 * exp_x3_z2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_z2 * r2 - 32 * exp_x3_z2 * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}},
              {{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_z2 + 160 * exp_x3_z2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_z2 * r2 - 32 * exp_x3_z2 * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                72 * exp_term_pos_z - 192 * exp_x2 * pos.z - 64 * exp_z3 + 160 * exp_x2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_x2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_x2 * z3 * r2}}}},
            {{{{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}},
              {{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}},
              {{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}}}},
            {{{{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_z2 + 160 * exp_x3_z2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_z2 * r2 - 32 * exp_x3_z2 * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                72 * exp_term_pos_z - 192 * exp_x2 * pos.z - 64 * exp_z3 + 160 * exp_x2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_x2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_x2 * z3 * r2}},
              {{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}},
              {{72 * exp_term_pos_z - 192 * exp_x2 * pos.z - 64 * exp_z3 + 160 * exp_x2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_x2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_x2 * z3 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2,
                72 * exp_term_pos_x - 384 * exp_term_pos_x_z2 + 160 * exp_term_pos_x * z4 - 24 * exp_term_pos_x * r2 +
                    96 * exp_term_pos_x_z2 * r2 - 32 * exp_term_pos_x * z4 * r2}}}}}}}},
        {{{{{{{{72 * exp_term_pos_y - 384 * exp_x2 * pos.y + 160 * exp_term_x4 * pos.y - 24 * exp_term_pos_y * r2 +
                    96 * exp_x2 * pos.y * r2 - 32 * exp_term_x4 * pos.y * r2,
                72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_y2 + 160 * exp_x3_y2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_y2 * r2 - 32 * exp_x3_y2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2}},
              {{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_y2 + 160 * exp_x3_y2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_y2 * r2 - 32 * exp_x3_y2 * r2,
                72 * exp_term_pos_y - 192 * exp_x2 * pos.y - 64 * exp_y3 + 160 * exp_x2 * y3 -
                    24 * exp_term_pos_y * r2 + 48 * exp_x2 * pos.y * r2 + 16 * exp_y3 * r2 - 32 * exp_x2 * y3 * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}}}},
            {{{{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_y2 + 160 * exp_x3_y2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_y2 * r2 - 32 * exp_x3_y2 * r2,
                72 * exp_term_pos_y - 192 * exp_x2 * pos.y - 64 * exp_y3 + 160 * exp_x2 * y3 -
                    24 * exp_term_pos_y * r2 + 48 * exp_x2 * pos.y * r2 + 16 * exp_y3 * r2 - 32 * exp_x2 * y3 * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2}},
              {{72 * exp_term_pos_y - 192 * exp_x2 * pos.y - 64 * exp_y3 + 160 * exp_x2 * y3 -
                    24 * exp_term_pos_y * r2 + 48 * exp_x2 * pos.y * r2 + 16 * exp_y3 * r2 - 32 * exp_x2 * y3 * r2,
                72 * exp_term_pos_x - 384 * exp_term_pos_x_y2 + 160 * exp_term_pos_x * y4 - 24 * exp_term_pos_x * r2 +
                    96 * exp_term_pos_x_y2 * r2 - 32 * exp_term_pos_x * y4 * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2}},
              {{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}}}},
            {{{{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}},
              {{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}},
              {{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}}}}}},
          {{{{{{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_y2 + 160 * exp_x3_y2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_y2 * r2 - 32 * exp_x3_y2 * r2,
                72 * exp_term_pos_y - 192 * exp_x2 * pos.y - 64 * exp_y3 + 160 * exp_x2 * y3 -
                    24 * exp_term_pos_y * r2 + 48 * exp_x2 * pos.y * r2 + 16 * exp_y3 * r2 - 32 * exp_x2 * y3 * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2}},
              {{72 * exp_term_pos_y - 192 * exp_x2 * pos.y - 64 * exp_y3 + 160 * exp_x2 * y3 -
                    24 * exp_term_pos_y * r2 + 48 * exp_x2 * pos.y * r2 + 16 * exp_y3 * r2 - 32 * exp_x2 * y3 * r2,
                72 * exp_term_pos_x - 384 * exp_term_pos_x_y2 + 160 * exp_term_pos_x * y4 - 24 * exp_term_pos_x * r2 +
                    96 * exp_term_pos_x_y2 * r2 - 32 * exp_term_pos_x * y4 * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2}},
              {{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}}}},
            {{{{72 * exp_term_pos_y - 192 * exp_x2 * pos.y - 64 * exp_y3 + 160 * exp_x2 * y3 -
                    24 * exp_term_pos_y * r2 + 48 * exp_x2 * pos.y * r2 + 16 * exp_y3 * r2 - 32 * exp_x2 * y3 * r2,
                72 * exp_term_pos_x - 384 * exp_term_pos_x_y2 + 160 * exp_term_pos_x * y4 - 24 * exp_term_pos_x * r2 +
                    96 * exp_term_pos_x_y2 * r2 - 32 * exp_term_pos_x * y4 * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2}},
              {{72 * exp_term_pos_x - 384 * exp_term_pos_x_y2 + 160 * exp_term_pos_x * y4 - 24 * exp_term_pos_x * r2 +
                    96 * exp_term_pos_x_y2 * r2 - 32 * exp_term_pos_x * y4 * r2,
                360 * exp_term_pos_y - 640 * exp_y3 + 160 * exp_y2 * y3 - 120 * exp_term_pos_y * r2 +
                    160 * exp_y3 * r2 - 32 * exp_y2 * y3 * r2,
                72 * exp_term_pos_z - 384 * exp_y2 * pos.z + 160 * exp_term_y4 * pos.z - 24 * exp_term_pos_z * r2 +
                    96 * exp_y2 * pos.z * r2 - 32 * exp_term_y4 * pos.z * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                72 * exp_term_pos_z - 384 * exp_y2 * pos.z + 160 * exp_term_y4 * pos.z - 24 * exp_term_pos_z * r2 +
                    96 * exp_y2 * pos.z * r2 - 32 * exp_term_y4 * pos.z * r2,
                72 * exp_term_pos_y - 64 * exp_y3 - 192 * exp_term_pos_y_z2 + 160 * exp_y3_z2 -
                    24 * exp_term_pos_y * r2 + 16 * exp_y3 * r2 + 48 * exp_term_pos_y_z2 * r2 - 32 * exp_y3_z2 * r2}}}},
            {{{{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                72 * exp_term_pos_z - 384 * exp_y2 * pos.z + 160 * exp_term_y4 * pos.z - 24 * exp_term_pos_z * r2 +
                    96 * exp_y2 * pos.z * r2 - 32 * exp_term_y4 * pos.z * r2,
                72 * exp_term_pos_y - 64 * exp_y3 - 192 * exp_term_pos_y_z2 + 160 * exp_y3_z2 -
                    24 * exp_term_pos_y * r2 + 16 * exp_y3 * r2 + 48 * exp_term_pos_y_z2 * r2 - 32 * exp_y3_z2 * r2}},
              {{24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                72 * exp_term_pos_y - 64 * exp_y3 - 192 * exp_term_pos_y_z2 + 160 * exp_y3_z2 -
                    24 * exp_term_pos_y * r2 + 16 * exp_y3 * r2 + 48 * exp_term_pos_y_z2 * r2 - 32 * exp_y3_z2 * r2,
                72 * exp_term_pos_z - 192 * exp_y2 * pos.z - 64 * exp_z3 + 160 * exp_y2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_y2 * pos.z * r2 + 16 * exp_z3 * r2 -
                    32 * exp_y2 * z3 * r2}}}}}},
          {{{{{{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}},
              {{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}},
              {{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}}}},
            {{{{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                72 * exp_term_pos_z - 384 * exp_y2 * pos.z + 160 * exp_term_y4 * pos.z - 24 * exp_term_pos_z * r2 +
                    96 * exp_y2 * pos.z * r2 - 32 * exp_term_y4 * pos.z * r2,
                72 * exp_term_pos_y - 64 * exp_y3 - 192 * exp_term_pos_y_z2 + 160 * exp_y3_z2 -
                    24 * exp_term_pos_y * r2 + 16 * exp_y3 * r2 + 48 * exp_term_pos_y_z2 * r2 - 32 * exp_y3_z2 * r2}},
              {{24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                72 * exp_term_pos_y - 64 * exp_y3 - 192 * exp_term_pos_y_z2 + 160 * exp_y3_z2 -
                    24 * exp_term_pos_y * r2 + 16 * exp_y3 * r2 + 48 * exp_term_pos_y_z2 * r2 - 32 * exp_y3_z2 * r2,
                72 * exp_term_pos_z - 192 * exp_y2 * pos.z - 64 * exp_z3 + 160 * exp_y2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_y2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_y2 * z3 * r2}}}},
            {{{{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}},
              {{24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                72 * exp_term_pos_y - 64 * exp_y3 - 192 * exp_term_pos_y_z2 + 160 * exp_y3_z2 -
                    24 * exp_term_pos_y * r2 + 16 * exp_y3 * r2 + 48 * exp_term_pos_y_z2 * r2 - 32 * exp_y3_z2 * r2,
                72 * exp_term_pos_z - 192 * exp_y2 * pos.z - 64 * exp_z3 + 160 * exp_y2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_y2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_y2 * z3 * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2,
                72 * exp_term_pos_z - 192 * exp_y2 * pos.z - 64 * exp_z3 + 160 * exp_y2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_y2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_y2 * z3 * r2,
                72 * exp_term_pos_y - 384 * exp_term_pos_y_z2 + 160 * exp_term_pos_y * z4 - 24 * exp_term_pos_y * r2 +
                    96 * exp_term_pos_y_z2 * r2 - 32 * exp_term_pos_y * z4 * r2}}}}}}}},
        {{{{{{{{72 * exp_term_pos_z - 384 * exp_x2 * pos.z + 160 * exp_term_x4 * pos.z - 24 * exp_term_pos_z * r2 +
                    96 * exp_x2 * pos.z * r2 - 32 * exp_term_x4 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_z2 + 160 * exp_x3_z2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_z2 * r2 - 32 * exp_x3_z2 * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}},
              {{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_z2 + 160 * exp_x3_z2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_z2 * r2 - 32 * exp_x3_z2 * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                72 * exp_term_pos_z - 192 * exp_x2 * pos.z - 64 * exp_z3 + 160 * exp_x2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_x2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_x2 * z3 * r2}}}},
            {{{{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}},
              {{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}},
              {{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}}}},
            {{{{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_z2 + 160 * exp_x3_z2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_z2 * r2 - 32 * exp_x3_z2 * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                72 * exp_term_pos_z - 192 * exp_x2 * pos.z - 64 * exp_z3 + 160 * exp_x2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_x2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_x2 * z3 * r2}},
              {{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}},
              {{72 * exp_term_pos_z - 192 * exp_x2 * pos.z - 64 * exp_z3 + 160 * exp_x2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_x2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_x2 * z3 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2,
                72 * exp_term_pos_x - 384 * exp_term_pos_x_z2 + 160 * exp_term_pos_x * z4 - 24 * exp_term_pos_x * r2 +
                    96 * exp_term_pos_x_z2 * r2 - 32 * exp_term_pos_x * z4 * r2}}}}}},
          {{{{{{-192 * exp_xy * pos.z + 160 * exp_x2_xy * pos.z + 48 * exp_xy * pos.z * r2 -
                    32 * exp_x2_xy * pos.z * r2,
                24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2}},
              {{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}},
              {{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}}}},
            {{{{24 * exp_term_pos_z - 64 * exp_x2 * pos.z - 64 * exp_y2 * pos.z + 160 * exp_x2 * y2 * pos.z -
                    8 * exp_term_pos_z * r2 + 16 * exp_x2 * pos.z * r2 + 16 * exp_y2 * pos.z * r2 -
                    32 * exp_x2 * y2 * pos.z * r2,
                -192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_term_pos_x_y2 * yz + 48 * exp_xy * pos.z * r2 -
                    32 * exp_term_pos_x_y2 * yz * r2,
                72 * exp_term_pos_z - 384 * exp_y2 * pos.z + 160 * exp_term_y4 * pos.z - 24 * exp_term_pos_z * r2 +
                    96 * exp_y2 * pos.z * r2 - 32 * exp_term_y4 * pos.z * r2,
                72 * exp_term_pos_y - 64 * exp_y3 - 192 * exp_term_pos_y_z2 + 160 * exp_y3_z2 -
                    24 * exp_term_pos_y * r2 + 16 * exp_y3 * r2 + 48 * exp_term_pos_y_z2 * r2 - 32 * exp_y3_z2 * r2}},
              {{24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                72 * exp_term_pos_y - 64 * exp_y3 - 192 * exp_term_pos_y_z2 + 160 * exp_y3_z2 -
                    24 * exp_term_pos_y * r2 + 16 * exp_y3 * r2 + 48 * exp_term_pos_y_z2 * r2 - 32 * exp_y3_z2 * r2,
                72 * exp_term_pos_z - 192 * exp_y2 * pos.z - 64 * exp_z3 + 160 * exp_y2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_y2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_y2 * z3 * r2}}}},
            {{{{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}},
              {{24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                72 * exp_term_pos_y - 64 * exp_y3 - 192 * exp_term_pos_y_z2 + 160 * exp_y3_z2 -
                    24 * exp_term_pos_y * r2 + 16 * exp_y3 * r2 + 48 * exp_term_pos_y_z2 * r2 - 32 * exp_y3_z2 * r2,
                72 * exp_term_pos_z - 192 * exp_y2 * pos.z - 64 * exp_z3 + 160 * exp_y2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_y2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_y2 * z3 * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2,
                72 * exp_term_pos_z - 192 * exp_y2 * pos.z - 64 * exp_z3 + 160 * exp_y2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_y2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_y2 * z3 * r2,
                72 * exp_term_pos_y - 384 * exp_term_pos_y_z2 + 160 * exp_term_pos_y * z4 - 24 * exp_term_pos_y * r2 +
                    96 * exp_term_pos_y_z2 * r2 - 32 * exp_term_pos_y * z4 * r2}}}}}},
          {{{{{{72 * exp_term_pos_x - 64 * exp_x3 - 192 * exp_term_pos_x_z2 + 160 * exp_x3_z2 -
                    24 * exp_term_pos_x * r2 + 16 * exp_x3 * r2 + 48 * exp_term_pos_x_z2 * r2 - 32 * exp_x3_z2 * r2,
                24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                72 * exp_term_pos_z - 192 * exp_x2 * pos.z - 64 * exp_z3 + 160 * exp_x2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_x2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_x2 * z3 * r2}},
              {{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}},
              {{72 * exp_term_pos_z - 192 * exp_x2 * pos.z - 64 * exp_z3 + 160 * exp_x2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_x2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_x2 * z3 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2,
                72 * exp_term_pos_x - 384 * exp_term_pos_x_z2 + 160 * exp_term_pos_x * z4 - 24 * exp_term_pos_x * r2 +
                    96 * exp_term_pos_x_z2 * r2 - 32 * exp_term_pos_x * z4 * r2}}}},
            {{{{24 * exp_term_pos_y - 64 * exp_x2 * pos.y - 64 * exp_term_pos_y_z2 + 160 * exp_x2 * pos.y * z2 -
                    8 * exp_term_pos_y * r2 + 16 * exp_x2 * pos.y * r2 + 16 * exp_term_pos_y_z2 * r2 -
                    32 * exp_x2 * pos.y * z2 * r2,
                24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2}},
              {{24 * exp_term_pos_x - 64 * exp_term_pos_x_y2 - 64 * exp_term_pos_x_z2 + 160 * exp_term_pos_x_y2 * z2 -
                    8 * exp_term_pos_x * r2 + 16 * exp_term_pos_x_y2 * r2 + 16 * exp_term_pos_x_z2 * r2 -
                    32 * exp_term_pos_x_y2 * z2 * r2,
                72 * exp_term_pos_y - 64 * exp_y3 - 192 * exp_term_pos_y_z2 + 160 * exp_y3_z2 -
                    24 * exp_term_pos_y * r2 + 16 * exp_y3 * r2 + 48 * exp_term_pos_y_z2 * r2 - 32 * exp_y3_z2 * r2,
                72 * exp_term_pos_z - 192 * exp_y2 * pos.z - 64 * exp_z3 + 160 * exp_y2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_y2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_y2 * z3 * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2,
                72 * exp_term_pos_z - 192 * exp_y2 * pos.z - 64 * exp_z3 + 160 * exp_y2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_y2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_y2 * z3 * r2,
                72 * exp_term_pos_y - 384 * exp_term_pos_y_z2 + 160 * exp_term_pos_y * z4 - 24 * exp_term_pos_y * r2 +
                    96 * exp_term_pos_y_z2 * r2 - 32 * exp_term_pos_y * z4 * r2}}}},
            {{{{72 * exp_term_pos_z - 192 * exp_x2 * pos.z - 64 * exp_z3 + 160 * exp_x2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_x2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_x2 * z3 * r2,
                -192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2,
                72 * exp_term_pos_x - 384 * exp_term_pos_x_z2 + 160 * exp_term_pos_x * z4 - 24 * exp_term_pos_x * r2 +
                    96 * exp_term_pos_x_z2 * r2 - 32 * exp_term_pos_x * z4 * r2}},
              {{-192 * exp_xy * pos.z + 160 * exp_xy * z3 + 48 * exp_xy * pos.z * r2 - 32 * exp_xy * z3 * r2,
                72 * exp_term_pos_z - 192 * exp_y2 * pos.z - 64 * exp_z3 + 160 * exp_y2 * z3 -
                    24 * exp_term_pos_z * r2 + 48 * exp_y2 * pos.z * r2 + 16 * exp_z3 * r2 - 32 * exp_y2 * z3 * r2,
                72 * exp_term_pos_y - 384 * exp_term_pos_y_z2 + 160 * exp_term_pos_y * z4 - 24 * exp_term_pos_y * r2 +
                    96 * exp_term_pos_y_z2 * r2 - 32 * exp_term_pos_y * z4 * r2}},
              {{72 * exp_term_pos_x - 384 * exp_term_pos_x_z2 + 160 * exp_term_pos_x * z4 - 24 * exp_term_pos_x * r2 +
                    96 * exp_term_pos_x_z2 * r2 - 32 * exp_term_pos_x * z4 * r2,
                72 * exp_term_pos_y - 384 * exp_term_pos_y_z2 + 160 * exp_term_pos_y * z4 - 24 * exp_term_pos_y * r2 +
                    96 * exp_term_pos_y_z2 * r2 - 32 * exp_term_pos_y * z4 * r2,
                360 * exp_term_pos_z - 640 * exp_z3 + 160 * exp_z2 * z3 - 120 * exp_term_pos_z * r2 +
                    160 * exp_z3 * r2 - 32 * exp_z2 * z3 * r2}}}}}}}}}},
      {{{{{{{{{{{{360 * exp_term - 2880 * exp_x2 + 2400 * exp_term_x4 - 384 * exp_term * x6 - 120 * exp_r2 +
                      720 * exp_x2_r2 - 480 * exp_term_x4 * r2 + 64 * exp_term * x6 * r2,
                  -960 * exp_xy + 1600 * exp_x2_xy - 384 * exp_term_x4 * xy + 240 * exp_xy_r2 - 320 * exp_x2_xy * r2 +
                      64 * exp_term_x4 * xy * r2,
                  -960 * exp_xz + 1600 * exp_x2_xz - 384 * exp_term_x4 * xz + 240 * exp_xz_r2 - 320 * exp_x2_xz * r2 +
                      64 * exp_term_x4 * xz * r2}},
                {{-960 * exp_xy + 1600 * exp_x2_xy - 384 * exp_term_x4 * xy + 240 * exp_xy_r2 - 320 * exp_x2_xy * r2 +
                      64 * exp_term_x4 * xy * r2,
                  72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2}},
                {{-960 * exp_xz + 1600 * exp_x2_xz - 384 * exp_term_x4 * xz + 240 * exp_xz_r2 - 320 * exp_x2_xz * r2 +
                      64 * exp_term_x4 * xz * r2,
                  -192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2}}}},
              {{{{-960 * exp_xy + 1600 * exp_x2_xy - 384 * exp_term_x4 * xy + 240 * exp_xy_r2 - 320 * exp_x2_xy * r2 +
                      64 * exp_term_x4 * xy * r2,
                  72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2}},
                {{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2}},
                {{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}}}},
              {{{{-960 * exp_xz + 1600 * exp_x2_xz - 384 * exp_term_x4 * xz + 240 * exp_xz_r2 - 320 * exp_x2_xz * r2 +
                      64 * exp_term_x4 * xz * r2,
                  -192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2}},
                {{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  -576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2}}}}}},
            {{{{{{-960 * exp_xy + 1600 * exp_x2_xy - 384 * exp_term_x4 * xy + 240 * exp_xy_r2 - 320 * exp_x2_xy * r2 +
                      64 * exp_term_x4 * xy * r2,
                  72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2}},
                {{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2}},
                {{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}}}},
              {{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2}},
                {{-576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}}}},
              {{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}}}},
            {{{{{{-960 * exp_xz + 1600 * exp_x2_xz - 384 * exp_term_x4 * xz + 240 * exp_xz_r2 - 320 * exp_x2_xz * r2 +
                      64 * exp_term_x4 * xz * r2,
                  -192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2}},
                {{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  -576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2}}}},
              {{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}},
              {{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  -576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{-576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2}}}}}}}},
          {{{{{{{{-960 * exp_xy + 1600 * exp_x2_xy - 384 * exp_term_x4 * xy + 240 * exp_xy_r2 - 320 * exp_x2_xy * r2 +
                      64 * exp_term_x4 * xy * r2,
                  72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2}},
                {{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2}},
                {{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}}}},
              {{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2}},
                {{-576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}}}},
              {{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}}}},
            {{{{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2}},
                {{-576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}}}},
              {{{{-576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2}},
                {{72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -960 * exp_xy + 1600 * exp_term_pos_x * y3 - 384 * exp_term_pos_x_y2 * y3 + 240 * exp_xy_r2 -
                      320 * exp_term_pos_x * y3 * r2 + 64 * exp_term_pos_x_y2 * y3 * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}}}},
              {{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}}}},
            {{{{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}},
              {{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}},
              {{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}}}}}},
          {{{{{{{{-960 * exp_xz + 1600 * exp_x2_xz - 384 * exp_term_x4 * xz + 240 * exp_xz_r2 - 320 * exp_x2_xz * r2 +
                      64 * exp_term_x4 * xz * r2,
                  -192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2}},
                {{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  -576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2}}}},
              {{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}},
              {{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  -576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{-576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2}}}}}},
            {{{{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}},
              {{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}},
              {{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}}}},
            {{{{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  -576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{-576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2}}}},
              {{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}},
              {{{{-576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}},
                {{72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2,
                  -960 * exp_xz + 1600 * exp_term_pos_x * z3 - 384 * exp_term_pos_x_z2 * z3 + 240 * exp_xz_r2 -
                      320 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_z2 * z3 * r2}}}}}}}}}},
        {{{{{{{{{{-960 * exp_xy + 1600 * exp_x2_xy - 384 * exp_term_x4 * xy + 240 * exp_xy_r2 - 320 * exp_x2_xy * r2 +
                      64 * exp_term_x4 * xy * r2,
                  72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2}},
                {{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2}},
                {{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}}}},
              {{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2}},
                {{-576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}}}},
              {{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}}}},
            {{{{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2}},
                {{-576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}}}},
              {{{{-576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2}},
                {{72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -960 * exp_xy + 1600 * exp_term_pos_x * y3 - 384 * exp_term_pos_x_y2 * y3 + 240 * exp_xy_r2 -
                      320 * exp_term_pos_x * y3 * r2 + 64 * exp_term_pos_x_y2 * y3 * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}}}},
              {{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}}}},
            {{{{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}},
              {{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}},
              {{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}}}}}},
          {{{{{{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_y2 + 960 * exp_x2 * y2 -
                      384 * exp_term_x4_y2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_y2_r2 -
                      192 * exp_x2 * y2 * r2 + 64 * exp_term_x4_y2 * r2,
                  -576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2}},
                {{-576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}}}},
              {{{{-576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2}},
                {{72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -960 * exp_xy + 1600 * exp_term_pos_x * y3 - 384 * exp_term_pos_x_y2 * y3 + 240 * exp_xy_r2 -
                      320 * exp_term_pos_x * y3 * r2 + 64 * exp_term_pos_x_y2 * y3 * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}}}},
              {{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}}}},
            {{{{{{-576 * exp_xy + 480 * exp_x2_xy + 480 * exp_term_pos_x * y3 - 384 * exp_x3 * y3 + 144 * exp_xy_r2 -
                      96 * exp_x2_xy * r2 - 96 * exp_term_pos_x * y3 * r2 + 64 * exp_x3 * y3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2}},
                {{72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -960 * exp_xy + 1600 * exp_term_pos_x * y3 - 384 * exp_term_pos_x_y2 * y3 + 240 * exp_xy_r2 -
                      320 * exp_term_pos_x * y3 * r2 + 64 * exp_term_pos_x_y2 * y3 * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}}}},
              {{{{72 * exp_term - 192 * exp_x2 - 384 * exp_y2 + 960 * exp_x2 * y2 + 160 * exp_term_y4 -
                      384 * exp_x2 * y4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_y2_r2 - 192 * exp_x2 * y2 * r2 -
                      32 * exp_term_y4 * r2 + 64 * exp_x2 * y4 * r2,
                  -960 * exp_xy + 1600 * exp_term_pos_x * y3 - 384 * exp_term_pos_x_y2 * y3 + 240 * exp_xy_r2 -
                      320 * exp_term_pos_x * y3 * r2 + 64 * exp_term_pos_x_y2 * y3 * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2}},
                {{-960 * exp_xy + 1600 * exp_term_pos_x * y3 - 384 * exp_term_pos_x_y2 * y3 + 240 * exp_xy_r2 -
                      320 * exp_term_pos_x * y3 * r2 + 64 * exp_term_pos_x_y2 * y3 * r2,
                  360 * exp_term - 2880 * exp_y2 + 2400 * exp_term_y4 - 384 * exp_term * y6 - 120 * exp_r2 +
                      720 * exp_y2_r2 - 480 * exp_term_y4 * r2 + 64 * exp_term * y6 * r2,
                  -960 * exp_yz + 1600 * exp_y2_yz - 384 * exp_term_y4 * yz + 240 * exp_yz_r2 - 320 * exp_y2_yz * r2 +
                      64 * exp_term_y4 * yz * r2}},
                {{-192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -960 * exp_yz + 1600 * exp_y2_yz - 384 * exp_term_y4 * yz + 240 * exp_yz_r2 - 320 * exp_y2_yz * r2 +
                      64 * exp_term_y4 * yz * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2}}}},
              {{{{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{-192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -960 * exp_yz + 1600 * exp_y2_yz - 384 * exp_term_y4 * yz + 240 * exp_yz_r2 - 320 * exp_y2_yz * r2 +
                      64 * exp_term_y4 * yz * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2}},
                {{-192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2}}}}}},
            {{{{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}},
              {{{{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{-192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -960 * exp_yz + 1600 * exp_y2_yz - 384 * exp_term_y4 * yz + 240 * exp_yz_r2 - 320 * exp_y2_yz * r2 +
                      64 * exp_term_y4 * yz * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2}},
                {{-192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2}}}},
              {{{{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2}},
                {{-192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2}}}}}}}},
          {{{{{{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}},
              {{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}},
              {{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}}}},
            {{{{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}},
              {{{{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{-192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -960 * exp_yz + 1600 * exp_y2_yz - 384 * exp_term_y4 * yz + 240 * exp_yz_r2 - 320 * exp_y2_yz * r2 +
                      64 * exp_term_y4 * yz * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2}},
                {{-192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2}}}},
              {{{{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2}},
                {{-192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2}}}}}},
            {{{{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}},
              {{{{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2}},
                {{-192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2}}}},
              {{{{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}},
                {{-192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2}},
                {{-192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2,
                  -960 * exp_yz + 1600 * exp_term_pos_y_z3 - 384 * exp_term_pos_y_z2 * z3 + 240 * exp_yz_r2 -
                      320 * exp_term_pos_y_z3 * r2 + 64 * exp_term_pos_y_z2 * z3 * r2}}}}}}}}}},
        {{{{{{{{{{-960 * exp_xz + 1600 * exp_x2_xz - 384 * exp_term_x4 * xz + 240 * exp_xz_r2 - 320 * exp_x2_xz * r2 +
                      64 * exp_term_x4 * xz * r2,
                  -192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2}},
                {{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  -576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2}}}},
              {{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}},
              {{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  -576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{-576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2}}}}}},
            {{{{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}},
              {{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}},
              {{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}}}},
            {{{{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  -576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{-576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2}}}},
              {{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}},
              {{{{-576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}},
                {{72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2,
                  -960 * exp_xz + 1600 * exp_term_pos_x * z3 - 384 * exp_term_pos_x_z2 * z3 + 240 * exp_xz_r2 -
                      320 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_z2 * z3 * r2}}}}}}}},
          {{{{{{{{-192 * exp_yz + 960 * exp_x2_yz - 384 * exp_term_x4 * yz + 48 * exp_yz_r2 - 192 * exp_x2_yz * r2 +
                      64 * exp_term_x4 * yz * r2,
                  -192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2}},
                {{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}}}},
              {{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}},
              {{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}}}},
            {{{{{{-192 * exp_xz + 160 * exp_x2_xz + 480 * exp_term_pos_x_y2 * pos.z - 384 * exp_x3_y2 * pos.z +
                      48 * exp_xz_r2 - 32 * exp_x2_xz * r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 +
                      64 * exp_x3_y2 * pos.z * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}}}},
              {{{{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_y2_yz - 384 * exp_x2 * y2 * yz + 48 * exp_yz_r2 -
                      96 * exp_x2_yz * r2 - 32 * exp_y2_yz * r2 + 64 * exp_x2 * y2 * yz * r2,
                  -192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2}},
                {{-192 * exp_xz + 960 * exp_term_pos_x_y2 * pos.z - 384 * exp_term_pos_x * y4 * pos.z + 48 * exp_xz_r2 -
                      192 * exp_term_pos_x_y2 * pos.z * r2 + 64 * exp_term_pos_x * y4 * pos.z * r2,
                  -960 * exp_yz + 1600 * exp_y2_yz - 384 * exp_term_y4 * yz + 240 * exp_yz_r2 - 320 * exp_y2_yz * r2 +
                      64 * exp_term_y4 * yz * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2}},
                {{-192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2}}}},
              {{{{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2}},
                {{-192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2}}}}}},
            {{{{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}},
              {{{{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2}},
                {{-192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2}}}},
              {{{{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}},
                {{-192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2}},
                {{-192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2,
                  -960 * exp_yz + 1600 * exp_term_pos_y_z3 - 384 * exp_term_pos_y_z2 * z3 + 240 * exp_yz_r2 -
                      320 * exp_term_pos_y_z3 * r2 + 64 * exp_term_pos_y_z2 * z3 * r2}}}}}}}},
          {{{{{{{{72 * exp_term - 384 * exp_x2 + 160 * exp_term_x4 - 192 * exp_z2 + 960 * exp_x2 * z2 -
                      384 * exp_term_x4_z2 - 24 * exp_r2 + 96 * exp_x2_r2 - 32 * exp_term_x4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_x2 * z2 * r2 + 64 * exp_term_x4_z2 * r2,
                  -192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  -576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{-576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2}}}},
              {{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}},
              {{{{-576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}},
                {{72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2,
                  -960 * exp_xz + 1600 * exp_term_pos_x * z3 - 384 * exp_term_pos_x_z2 * z3 + 240 * exp_xz_r2 -
                      320 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_z2 * z3 * r2}}}}}},
            {{{{{{-192 * exp_xy + 160 * exp_x2_xy + 480 * exp_xy_z2 - 384 * exp_x2_xy_z2 + 48 * exp_xy_r2 -
                      32 * exp_x2_xy * r2 - 96 * exp_xy_z2 * r2 + 64 * exp_x2_xy_z2 * r2,
                  24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2}},
                {{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}}}},
              {{{{24 * exp_term - 64 * exp_x2 - 64 * exp_y2 + 160 * exp_x2 * y2 - 64 * exp_z2 + 160 * exp_x2 * z2 +
                      160 * exp_y2 * z2 - 384 * exp_x2 * y2 * z2 - 8 * exp_r2 + 16 * exp_x2_r2 + 16 * exp_y2_r2 -
                      32 * exp_x2 * y2 * r2 + 16 * exp_z2_r2 - 32 * exp_x2 * z2 * r2 - 32 * exp_y2 * z2 * r2 +
                      64 * exp_x2 * y2 * z2 * r2,
                  -192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2}},
                {{-192 * exp_xy + 160 * exp_term_pos_x * y3 + 480 * exp_xy_z2 - 384 * exp_term_pos_x * y3 * z2 +
                      48 * exp_xy_r2 - 32 * exp_term_pos_x * y3 * r2 - 96 * exp_xy_z2 * r2 +
                      64 * exp_term_pos_x * y3 * z2 * r2,
                  72 * exp_term - 384 * exp_y2 + 160 * exp_term_y4 - 192 * exp_z2 + 960 * exp_y2 * z2 -
                      384 * exp_term_y4_z2 - 24 * exp_r2 + 96 * exp_y2_r2 - 32 * exp_term_y4 * r2 + 48 * exp_z2_r2 -
                      192 * exp_y2 * z2 * r2 + 64 * exp_term_y4_z2 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2}},
                {{-192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2}}}},
              {{{{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}},
                {{-192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2}},
                {{-192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2,
                  -960 * exp_yz + 1600 * exp_term_pos_y_z3 - 384 * exp_term_pos_y_z2 * z3 + 240 * exp_yz_r2 -
                      320 * exp_term_pos_y_z3 * r2 + 64 * exp_term_pos_y_z2 * z3 * r2}}}}}},
            {{{{{{-576 * exp_xz + 480 * exp_x2_xz + 480 * exp_term_pos_x * z3 - 384 * exp_x3 * z3 + 144 * exp_xz_r2 -
                      96 * exp_x2_xz * r2 - 96 * exp_term_pos_x * z3 * r2 + 64 * exp_x3 * z3 * r2,
                  -192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2}},
                {{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}},
                {{72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2,
                  -960 * exp_xz + 1600 * exp_term_pos_x * z3 - 384 * exp_term_pos_x_z2 * z3 + 240 * exp_xz_r2 -
                      320 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_z2 * z3 * r2}}}},
              {{{{-192 * exp_yz + 480 * exp_x2_yz + 160 * exp_term_pos_y_z3 - 384 * exp_x2 * pos.y * z3 +
                      48 * exp_yz_r2 - 96 * exp_x2_yz * r2 - 32 * exp_term_pos_y_z3 * r2 +
                      64 * exp_x2 * pos.y * z3 * r2,
                  -192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2}},
                {{-192 * exp_xz + 480 * exp_term_pos_x_y2 * pos.z + 160 * exp_term_pos_x * z3 -
                      384 * exp_term_pos_x_y2 * z3 + 48 * exp_xz_r2 - 96 * exp_term_pos_x_y2 * pos.z * r2 -
                      32 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_y2 * z3 * r2,
                  -576 * exp_yz + 480 * exp_y2_yz + 480 * exp_term_pos_y_z3 - 384 * exp_y3_z3 + 144 * exp_yz_r2 -
                      96 * exp_y2_yz * r2 - 96 * exp_term_pos_y_z3 * r2 + 64 * exp_y3_z3 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2}},
                {{-192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2,
                  -960 * exp_yz + 1600 * exp_term_pos_y_z3 - 384 * exp_term_pos_y_z2 * z3 + 240 * exp_yz_r2 -
                      320 * exp_term_pos_y_z3 * r2 + 64 * exp_term_pos_y_z2 * z3 * r2}}}},
              {{{{72 * exp_term - 192 * exp_x2 - 384 * exp_z2 + 960 * exp_x2 * z2 + 160 * exp_term_z4 -
                      384 * exp_x2 * z4 - 24 * exp_r2 + 48 * exp_x2_r2 + 96 * exp_z2_r2 - 192 * exp_x2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_x2 * z4 * r2,
                  -192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2,
                  -960 * exp_xz + 1600 * exp_term_pos_x * z3 - 384 * exp_term_pos_x_z2 * z3 + 240 * exp_xz_r2 -
                      320 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_z2 * z3 * r2}},
                {{-192 * exp_xy + 960 * exp_xy_z2 - 384 * exp_xy_z4 + 48 * exp_xy_r2 - 192 * exp_xy_z2 * r2 +
                      64 * exp_xy_z4 * r2,
                  72 * exp_term - 192 * exp_y2 - 384 * exp_z2 + 960 * exp_y2 * z2 + 160 * exp_term_z4 -
                      384 * exp_y2 * z4 - 24 * exp_r2 + 48 * exp_y2_r2 + 96 * exp_z2_r2 - 192 * exp_y2 * z2 * r2 -
                      32 * exp_term_z4 * r2 + 64 * exp_y2 * z4 * r2,
                  -960 * exp_yz + 1600 * exp_term_pos_y_z3 - 384 * exp_term_pos_y_z2 * z3 + 240 * exp_yz_r2 -
                      320 * exp_term_pos_y_z3 * r2 + 64 * exp_term_pos_y_z2 * z3 * r2}},
                {{-960 * exp_xz + 1600 * exp_term_pos_x * z3 - 384 * exp_term_pos_x_z2 * z3 + 240 * exp_xz_r2 -
                      320 * exp_term_pos_x * z3 * r2 + 64 * exp_term_pos_x_z2 * z3 * r2,
                  -960 * exp_yz + 1600 * exp_term_pos_y_z3 - 384 * exp_term_pos_y_z2 * z3 + 240 * exp_yz_r2 -
                      320 * exp_term_pos_y_z3 * r2 + 64 * exp_term_pos_y_z2 * z3 * r2,
                  360 * exp_term - 2880 * exp_z2 + 2400 * exp_term_z4 - 384 * exp_term * z6 - 120 * exp_r2 +
                      720 * exp_z2_r2 - 480 * exp_term_z4 * r2 + 64 * exp_term * z6 * r2}}}}}}}}}}}}};
}
