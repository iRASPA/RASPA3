module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>
#endif

module tricubic_derivatives_external_field;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
Potentials::potentialTricubicSecondOrderPolynomialTestFunction(const double3 &pos)
{
  return {pos.x * pos.x + pos.y * pos.y + pos.z * pos.z,
          {{2.0 * pos.x, 2.0 * pos.y, 2.0 * pos.z}},
          {{{{2, 0, 0}}, {{0, 2, 0}}, {{0, 0, 2}}}},
          {}};
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
Potentials::potentialTricubicThirdOrderPolynomialTestFunction(const double3 &pos)
{
  return {
      pos.x * pos.x * pos.x + 2.0 * pos.x * pos.x * pos.y * pos.z * pos.z - pos.y * pos.y * pos.z + 1.0,
      {{3 * pos.x * pos.x + 4 * pos.x * pos.y * pos.z * pos.z, -2 * pos.y * pos.z + 2 * pos.x * pos.x * pos.z * pos.z,
        -pos.y * pos.y + 4 * pos.x * pos.x * pos.y * pos.z}},
      {{{{6 * pos.x + 4 * pos.y * pos.z * pos.z, 4 * pos.x * pos.z * pos.z, 8 * pos.x * pos.y * pos.z}},
        {{4 * pos.x * pos.z * pos.z, -2 * pos.z, -2 * pos.y + 4 * pos.x * pos.x * pos.z}},
        {{8 * pos.x * pos.y * pos.z, -2 * pos.y + 4 * pos.x * pos.x * pos.z, 4 * pos.x * pos.x * pos.y}}}},
      {{{{{{6, 4 * pos.z * pos.z, 8 * pos.y * pos.z}},
          {{4 * pos.z * pos.z, 0, 8 * pos.x * pos.z}},
          {{8 * pos.y * pos.z, 8 * pos.x * pos.z, 8 * pos.x * pos.y}}}},
        {{{{4 * pos.z * pos.z, 0, 8 * pos.x * pos.z}}, {{0, 0, -2}}, {{8 * pos.x * pos.z, -2, 4 * pos.x * pos.x}}}},
        {{{{8 * pos.y * pos.z, 8 * pos.x * pos.z, 8 * pos.x * pos.y}},
          {{8 * pos.x * pos.z, -2, 4 * pos.x * pos.x}},
          {{8 * pos.x * pos.y, 4 * pos.x * pos.x, 0}}}}}}};
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
Potentials::potentialTricubicFourthOrderPolynomialTestFunction(const double3 &pos)
{
  return {(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
          {{4 * pos.x * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            4 * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            4 * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
          {{{{8 * pos.x * pos.x + 4 * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z), 8 * pos.x * pos.y,
              8 * pos.x * pos.z}},
            {{8 * pos.x * pos.y, 8 * pos.y * pos.y + 4 * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
              8 * pos.y * pos.z}},
            {{8 * pos.x * pos.z, 8 * pos.y * pos.z,
              8 * pos.z * pos.z + 4 * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}}}},
          {{{{{{24 * pos.x, 8 * pos.y, 8 * pos.z}}, {{8 * pos.y, 8 * pos.x, 0}}, {{8 * pos.z, 0, 8 * pos.x}}}},
            {{{{8 * pos.y, 8 * pos.x, 0}}, {{8 * pos.x, 24 * pos.y, 8 * pos.z}}, {{0, 8 * pos.z, 8 * pos.y}}}},
            {{{{8 * pos.z, 0, 8 * pos.x}}, {{0, 8 * pos.z, 8 * pos.y}}, {{8 * pos.x, 8 * pos.y, 24 * pos.z}}}}}}};
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
Potentials::potentialTricubicFifthOrderPolynomialTestFunction(const double3 &pos)
{
  return {
      pos.x * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y * pos.y * pos.z + 2.0 * pos.x * pos.y * pos.y * pos.y +
          4.0 * pos.x * pos.x * pos.y * pos.z * pos.z,
      {{2 * pos.y * pos.y * pos.y + 5 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y * pos.y * pos.z +
            8 * pos.x * pos.y * pos.z * pos.z,
        6 * pos.x * pos.y * pos.y + 3 * pos.x * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y * pos.z +
            4 * pos.x * pos.x * pos.z * pos.z,
        pos.x * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y * pos.y + 8 * pos.x * pos.x * pos.y * pos.z}},
      {{{{20 * pos.x * pos.x * pos.x * pos.y * pos.y * pos.y * pos.z + 8 * pos.y * pos.z * pos.z,
          6 * pos.y * pos.y + 15 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y * pos.z + 8 * pos.x * pos.z * pos.z,
          5 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y * pos.y + 16 * pos.x * pos.y * pos.z}},
        {{6 * pos.y * pos.y + 15 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y * pos.z + 8 * pos.x * pos.z * pos.z,
          12 * pos.x * pos.y + 6 * pos.x * pos.x * pos.x * pos.x * pos.x * pos.y * pos.z,
          3 * pos.x * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y + 8 * pos.x * pos.x * pos.z}},
        {{5 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y * pos.y + 16 * pos.x * pos.y * pos.z,
          3 * pos.x * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y + 8 * pos.x * pos.x * pos.z,
          8 * pos.x * pos.x * pos.y}}}},
      {{{{{{60 * pos.x * pos.x * pos.y * pos.y * pos.y * pos.z,
            60 * pos.x * pos.x * pos.x * pos.y * pos.y * pos.z + 8 * pos.z * pos.z,
            20 * pos.x * pos.x * pos.x * pos.y * pos.y * pos.y + 16 * pos.y * pos.z}},
          {{60 * pos.x * pos.x * pos.x * pos.y * pos.y * pos.z + 8 * pos.z * pos.z,
            12 * pos.y + 30 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.z,
            15 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y + 16 * pos.x * pos.z}},
          {{20 * pos.x * pos.x * pos.x * pos.y * pos.y * pos.y + 16 * pos.y * pos.z,
            15 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y + 16 * pos.x * pos.z, 16 * pos.x * pos.y}}}},
        {{{{60 * pos.x * pos.x * pos.x * pos.y * pos.y * pos.z + 8 * pos.z * pos.z,
            12 * pos.y + 30 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.z,
            15 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y + 16 * pos.x * pos.z}},
          {{12 * pos.y + 30 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.z,
            12 * pos.x + 6 * pos.x * pos.x * pos.x * pos.x * pos.x * pos.z,
            6 * pos.x * pos.x * pos.x * pos.x * pos.x * pos.y}},
          {{15 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y + 16 * pos.x * pos.z,
            6 * pos.x * pos.x * pos.x * pos.x * pos.x * pos.y, 8 * pos.x * pos.x}}}},
        {{{{20 * pos.x * pos.x * pos.x * pos.y * pos.y * pos.y + 16 * pos.y * pos.z,
            15 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y + 16 * pos.x * pos.z, 16 * pos.x * pos.y}},
          {{15 * pos.x * pos.x * pos.x * pos.x * pos.y * pos.y + 16 * pos.x * pos.z,
            6 * pos.x * pos.x * pos.x * pos.x * pos.x * pos.y, 8 * pos.x * pos.x}},
          {{16 * pos.x * pos.y, 8 * pos.x * pos.x, 0}}}}}}};
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
Potentials::potentialTricubicSixthOrderPolynomialTestFunction(const double3 &pos)
{
  return {
      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) *
          (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
      {{6 * pos.x * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
        6 * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
        6 * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
      {{{{24 * pos.x * pos.x * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) +
              6 * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
          24 * pos.x * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
          24 * pos.x * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
        {{24 * pos.x * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
          24 * pos.y * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) +
              6 * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
          24 * pos.y * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
        {{24 * pos.x * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
          24 * pos.y * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
          24 * pos.z * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) +
              6 * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}}}},
      {{{{{{48 * pos.x * pos.x * pos.x + 72 * pos.x * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.x * pos.x * pos.y + 24 * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.x * pos.x * pos.z + 24 * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
          {{48 * pos.x * pos.x * pos.y + 24 * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.x * pos.y * pos.y + 24 * pos.x * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.x * pos.y * pos.z}},
          {{48 * pos.x * pos.x * pos.z + 24 * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.x * pos.y * pos.z,
            48 * pos.x * pos.z * pos.z + 24 * pos.x * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}}}},
        {{{{48 * pos.x * pos.x * pos.y + 24 * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.x * pos.y * pos.y + 24 * pos.x * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.x * pos.y * pos.z}},
          {{48 * pos.x * pos.y * pos.y + 24 * pos.x * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.y * pos.y * pos.y + 72 * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.y * pos.y * pos.z + 24 * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
          {{48 * pos.x * pos.y * pos.z,
            48 * pos.y * pos.y * pos.z + 24 * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.y * pos.z * pos.z + 24 * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}}}},
        {{{{48 * pos.x * pos.x * pos.z + 24 * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.x * pos.y * pos.z,
            48 * pos.x * pos.z * pos.z + 24 * pos.x * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
          {{48 * pos.x * pos.y * pos.z,
            48 * pos.y * pos.y * pos.z + 24 * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.y * pos.z * pos.z + 24 * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
          {{48 * pos.x * pos.z * pos.z + 24 * pos.x * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.y * pos.z * pos.z + 24 * pos.y * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            48 * pos.z * pos.z * pos.z + 72 * pos.z * (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}}}}}}};
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
Potentials::potentialTricubicExponentialNonPolynomialTestFunction(const double3 &pos)
{
  return {(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) * std::exp(-(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)),
          {{2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x -
                2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x *
                    (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y -
                2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y *
                    (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
            2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z -
                2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z *
                    (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
          {{{{2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) -
                  8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x -
                  2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) +
                  4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
              -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y +
                  4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
              -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.z +
                  4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.z *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
            {{-8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y +
                  4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
              2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) -
                  8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.y -
                  2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) +
                  4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.y *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
              -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.z +
                  4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.z *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
            {{-8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.z +
                  4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.z *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
              -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.z +
                  4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.z *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
              2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) -
                  8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z * pos.z -
                  2 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) +
                  4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z * pos.z *
                      (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}}}},
          {{{{{{-24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.x +
                    12 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.x *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.y +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
              {{-8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.y +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.y +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
              {{-8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.z * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.z * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}}}},
            {{{{-8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.y +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.y +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
              {{-8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.y +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.y * pos.y +
                    12 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.y * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.y * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.y * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
              {{24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.y * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.y * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.z * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.z * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}}}},
            {{{{-8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.x * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.z * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.z * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
              {{24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.y * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.y * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.y * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.z * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.z * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}},
              {{-8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.z * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.x * pos.z * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.z * pos.z +
                    4 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.y * pos.z * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z),
                -24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z +
                    24 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z * pos.z * pos.z +
                    12 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) -
                    8 * std::exp(-pos.x * pos.x - pos.y * pos.y - pos.z * pos.z) * pos.z * pos.z * pos.z *
                        (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)}}}}}}};
}
