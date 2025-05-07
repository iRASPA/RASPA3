module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <deque>
#include <future>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <semaphore>
#include <span>
#include <thread>
#include <utility>
#include <vector>
#include <array>
#endif

module interactions_framework_molecule_grid;

#ifndef USE_LEGACY_HEADERS
import <numbers>;
import <optional>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <cmath>;
import <thread>;
import <future>;
import <deque>;
import <semaphore>;
import <atomic>;
import <utility>;
import <limits>;
#endif

import double3;
import double4;
import double3x3;
import double3x3x3;
import energy_status;
import potential_energy_vdw;
import potential_energy_coulomb;
import potential_gradient_vdw;
import potential_gradient_coulomb;
import potential_hessian_vdw;
import potential_hessian_coulomb;
import potential_tricubic_derivative_lj;
import potential_tricubic_derivative_real_ewald;
import potential_triquintic_derivative_lj;
import potential_triquintic_derivative_real_ewald;
import potential_electrostatics;
import simulationbox;
import forcefield;
import atom;
import energy_factor;
import energy_status_inter;
import running_energy;
import units;
import threadpool;
// import threading;
import energy_factor;
import gradient_factor;
import hessian_factor;
import tricubic_derivative_factor;
import triquintic_derivative_factor;
import framework;
import component;

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
Interactions::calculateTricubicDerivativeAtPosition(ForceField::InterpolationGridType interpolationGridType,
                                                    const ForceField &forceField, const SimulationBox &simulationBox,
                                                    double3 posA, size_t typeA, std::span<const Atom> frameworkAtoms)
{
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;

  double energy{0.0};
  std::array<double, 3> first_derivative{};
  std::array<std::array<double, 3>, 3> second_derivative{};
  std::array<std::array<std::array<double, 3>, 3>, 3> third_derivative{};
  Potentials::TricubicDerivativeFactor v{};

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posB = it1->position;
    size_t typeB = static_cast<size_t>(it1->type);
    double chargeB = it1->charge;

    double3 dr = posA - posB;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = std::max(0.1, double3::dot(dr, dr));

    if (rr < cutOffFrameworkVDWSquared)
    {
      switch (interpolationGridType)
      {
        case ForceField::InterpolationGridType::LennardJones:
          v = Potentials::potentialLennardJonesTricubicDerivative(forceField, rr, typeA, typeB);
          break;
        case ForceField::InterpolationGridType::LennardJonesRepulsion:
          v = Potentials::potentialLennardJonesRepulsionTricubicDerivative(forceField, rr, cutOffFrameworkVDWSquared,
                                                                           typeA, typeB);
          break;
        case ForceField::InterpolationGridType::LennardJonesAttraction:
          v = Potentials::potentialLennardJonesAttractionTricubicDerivative(forceField, rr, cutOffFrameworkVDWSquared,
                                                                            typeA, typeB);
          break;
        case ForceField::InterpolationGridType::EwaldReal:
          v = Potentials::potentialRealEwaldTricubicDerivative(forceField, rr, std::sqrt(rr), 1.0, chargeB);
          break;
      }

      energy += v.energy;

      for (size_t i = 0; i != 3; ++i)
      {
        first_derivative[i] += dr[i] * v.firstDerivativeFactor;

        for (size_t j = 0; j != 3; ++j)
        {
          second_derivative[i][j] +=
              v.secondDerivativeFactor * dr[i] * dr[j] + (i == j ? v.firstDerivativeFactor : 0.0);

          for (size_t k = 0; k != 3; ++k)
          {
            third_derivative[i][j][k] +=
                v.thirdDerivativeFactor * dr[i] * dr[j] * dr[k] + (j == k ? v.secondDerivativeFactor * dr[i] : 0.0) +
                (i == k ? v.secondDerivativeFactor * dr[j] : 0.0) + (i == j ? v.secondDerivativeFactor * dr[k] : 0.0);
          }
        }
      }
    }
  }

  return {energy, first_derivative, second_derivative, third_derivative};
}

std::array<double, 8> Interactions::calculateTricubicCartesianAtPosition(
    ForceField::InterpolationGridType interpolationGridType, const ForceField &forceField,
    const SimulationBox &simulationBox, double3 posA, size_t typeA, std::span<const Atom> frameworkAtoms)
{
  auto [energy, first_derivative, second_derivative, third_derivative, fourth_derivative, firth_derivative,
        sixth_derivative] =
      Interactions::calculateTriquinticDerivativeAtPosition(interpolationGridType, forceField, simulationBox, posA,
                                                            typeA, frameworkAtoms);

  return {energy,

          first_derivative[0],
          first_derivative[1],
          first_derivative[2],

          second_derivative[0][1],
          second_derivative[0][2],
          second_derivative[1][2],

          third_derivative[0][1][2]};
}

std::array<double, 8> Interactions::calculateTricubicFractionalAtPosition(
    ForceField::InterpolationGridType interpolationGridType, const ForceField &forceField,
    const SimulationBox &simulationBox, double3 posA, size_t typeA, const SimulationBox &frameworkBox,
    std::span<const Atom> frameworkAtoms)
{
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  double energy{0.0};
  double first_derivative_x{};
  double first_derivative_y{};
  double first_derivative_z{};
  double second_derivative_xy{};
  double second_derivative_xz{};
  double second_derivative_yz{};
  double third_derivative_xyz{};
  Potentials::TricubicDerivativeFactor v{};
  double energy_fractional{0.0};
  std::array<double, 3> first_derivative_fractional{};
  std::array<std::array<double, 3>, 3> second_derivative_fractional{};
  std::array<std::array<std::array<double, 3>, 3>, 3> third_derivative_fractional{};

  double3x3 frameworkCell = frameworkBox.cell;
  switch (frameworkBox.type)
  {
    case SimulationBox::Type::Rectangular:
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        double3 posB = it1->position;
        size_t typeB = static_cast<size_t>(it1->type);
        double chargeB = it1->charge;

        double3 dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        double rr = std::max(0.1, double3::dot(dr, dr));

        if (rr < cutOffFrameworkVDWSquared)
        {
          switch (interpolationGridType)
          {
            case ForceField::InterpolationGridType::LennardJones:
              v = Potentials::potentialLennardJonesTricubicDerivative(forceField, rr, typeA, typeB);
              break;
            case ForceField::InterpolationGridType::LennardJonesRepulsion:
              v = Potentials::potentialLennardJonesRepulsionTricubicDerivative(forceField, rr,
                                                                               cutOffFrameworkVDWSquared, typeA, typeB);
              break;
            case ForceField::InterpolationGridType::LennardJonesAttraction:
              v = Potentials::potentialLennardJonesAttractionTricubicDerivative(
                  forceField, rr, cutOffFrameworkVDWSquared, typeA, typeB);
              break;
            case ForceField::InterpolationGridType::EwaldReal:
              v = Potentials::potentialRealEwaldTricubicDerivative(forceField, rr, std::sqrt(rr), 1.0, chargeB);
              break;
          }

          energy += v.energy;

          first_derivative_x += dr[0] * v.firstDerivativeFactor;
          first_derivative_y += dr[1] * v.firstDerivativeFactor;
          first_derivative_z += dr[2] * v.firstDerivativeFactor;

          second_derivative_xy += v.secondDerivativeFactor * dr[0] * dr[1];
          second_derivative_xz += v.secondDerivativeFactor * dr[0] * dr[2];
          second_derivative_yz += v.secondDerivativeFactor * dr[1] * dr[2];

          third_derivative_xyz += v.thirdDerivativeFactor * dr[0] * dr[1] * dr[2];
        }
      }

      return {energy,
              first_derivative_x * frameworkCell[0][0],
              first_derivative_y * frameworkCell[1][1],
              first_derivative_z * frameworkCell[2][2],
              second_derivative_xy * frameworkCell[0][0] * frameworkCell[1][1],
              second_derivative_xz * frameworkCell[0][0] * frameworkCell[2][2],
              second_derivative_yz * frameworkCell[1][1] * frameworkCell[2][2],
              third_derivative_xyz * frameworkCell[0][0] * frameworkCell[1][1] * frameworkCell[2][2]};
    case SimulationBox::Type::Triclinic:
      auto [value, first_derivative, second_derivative, third_derivative] =
          Interactions::calculateTricubicDerivativeAtPosition(interpolationGridType, forceField, simulationBox, posA,
                                                              typeA, frameworkAtoms);
      energy_fractional = value;
      for (const size_t &p : std::array<size_t, 3>{0, 1, 2})
      {
        first_derivative_fractional[0] += frameworkCell[0][p] * first_derivative[p];
        first_derivative_fractional[1] += frameworkCell[1][p] * first_derivative[p];
        first_derivative_fractional[2] += frameworkCell[2][p] * first_derivative[p];

        for (const size_t &q : std::array<size_t, 3>{0, 1, 2})
        {
          second_derivative_fractional[0][1] += frameworkCell[0][p] * frameworkCell[1][q] * second_derivative[p][q];
          second_derivative_fractional[0][2] += frameworkCell[0][p] * frameworkCell[2][q] * second_derivative[p][q];
          second_derivative_fractional[1][2] += frameworkCell[1][p] * frameworkCell[2][q] * second_derivative[p][q];
          for (const size_t &r : std::array<size_t, 3>{0, 1, 2})
          {
            third_derivative_fractional[0][1][2] +=
                frameworkCell[0][p] * frameworkCell[1][q] * frameworkCell[2][r] * third_derivative[p][q][r];
          }
        }
      }

      return {energy_fractional,

              first_derivative_fractional[0],
              first_derivative_fractional[1],
              first_derivative_fractional[2],

              second_derivative_fractional[0][1],
              second_derivative_fractional[0][2],
              second_derivative_fractional[1][2],

              third_derivative_fractional[0][1][2]};
  }
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
Interactions::calculateTriquinticDerivativeAtPosition(ForceField::InterpolationGridType interpolationGridType,
                                                      const ForceField &forceField, const SimulationBox &simulationBox,
                                                      double3 posA, size_t typeA, std::span<const Atom> frameworkAtoms)
{
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;

  double energy{0.0};
  std::array<double, 3> first_derivative{};
  std::array<std::array<double, 3>, 3> second_derivative{};
  std::array<std::array<std::array<double, 3>, 3>, 3> third_derivative{};
  std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> fourth_derivative{};
  std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3> fifth_derivative{};
  std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3> sixth_derivative{};
  Potentials::TriquinticDerivativeFactor v{};

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posB = it1->position;
    double chargeB = it1->charge;
    size_t typeB = static_cast<size_t>(it1->type);

    double3 dr = posA - posB;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = std::max(0.1, double3::dot(dr, dr));

    if (rr < cutOffFrameworkVDWSquared)
    {
      switch (interpolationGridType)
      {
        case ForceField::InterpolationGridType::LennardJones:
          v = Potentials::potentialLennardJonesTriquinticDerivative(forceField, rr, typeA, typeB);
          break;
        case ForceField::InterpolationGridType::LennardJonesRepulsion:
          v = Potentials::potentialLennardJonesTriquinticDerivative(forceField, rr, typeA, typeB);
          break;
        case ForceField::InterpolationGridType::LennardJonesAttraction:
          v = Potentials::potentialLennardJonesTriquinticDerivative(forceField, rr, typeA, typeB);
          break;
        case ForceField::InterpolationGridType::EwaldReal:
          v = Potentials::potentialRealEwaldTriquinticDerivative(forceField, rr, std::sqrt(rr), 1.0, chargeB);
          break;
      }

      energy += v.energy;

      for (size_t i = 0; i != 3; ++i)
      {
        first_derivative[i] += dr[i] * v.firstDerivativeFactor;

        for (size_t j = 0; j != 3; ++j)
        {
          second_derivative[i][j] +=
              v.secondDerivativeFactor * dr[i] * dr[j] + (i == j ? v.firstDerivativeFactor : 0.0);

          for (size_t k = 0; k != 3; ++k)
          {
            third_derivative[i][j][k] +=
                v.thirdDerivativeFactor * dr[i] * dr[j] * dr[k] + (j == k ? v.secondDerivativeFactor * dr[i] : 0.0) +
                (i == k ? v.secondDerivativeFactor * dr[j] : 0.0) + (i == j ? v.secondDerivativeFactor * dr[k] : 0.0);

            for (size_t l = 0; l != 3; ++l)
            {
              fourth_derivative[i][j][k][l] += v.fourthDerivativeFactor * dr[i] * dr[j] * dr[k] * dr[l] +
                                               (i == j ? v.thirdDerivativeFactor * dr[k] * dr[l] : 0.0) +
                                               (i == k ? v.thirdDerivativeFactor * dr[j] * dr[l] : 0.0) +
                                               (i == l ? v.thirdDerivativeFactor * dr[j] * dr[k] : 0.0) +
                                               (j == k ? v.thirdDerivativeFactor * dr[i] * dr[l] : 0.0) +
                                               (j == l ? v.thirdDerivativeFactor * dr[i] * dr[k] : 0.0) +
                                               (k == l ? v.thirdDerivativeFactor * dr[i] * dr[j] : 0.0) +
                                               (i == j && k == l ? v.secondDerivativeFactor : 0.0) +
                                               (i == l && j == k ? v.secondDerivativeFactor : 0.0) +
                                               (i == k && j == l ? v.secondDerivativeFactor : 0.0);
              for (size_t m = 0; m != 3; ++m)
              {
                fifth_derivative[i][j][k][l][m] += v.fifthDerivativeFactor * dr[i] * dr[j] * dr[k] * dr[l] * dr[m] +
                                                   (i == j ? v.fourthDerivativeFactor * dr[k] * dr[l] * dr[m] : 0.0) +
                                                   (i == k ? v.fourthDerivativeFactor * dr[j] * dr[l] * dr[m] : 0.0) +
                                                   (i == l ? v.fourthDerivativeFactor * dr[j] * dr[k] * dr[m] : 0.0) +
                                                   (i == m ? v.fourthDerivativeFactor * dr[j] * dr[k] * dr[l] : 0.0) +
                                                   (j == k ? v.fourthDerivativeFactor * dr[i] * dr[l] * dr[m] : 0.0) +
                                                   (j == l ? v.fourthDerivativeFactor * dr[i] * dr[k] * dr[m] : 0.0) +
                                                   (j == m ? v.fourthDerivativeFactor * dr[i] * dr[k] * dr[l] : 0.0) +
                                                   (k == l ? v.fourthDerivativeFactor * dr[i] * dr[j] * dr[m] : 0.0) +
                                                   (k == m ? v.fourthDerivativeFactor * dr[i] * dr[j] * dr[l] : 0.0) +
                                                   (l == m ? v.fourthDerivativeFactor * dr[i] * dr[j] * dr[k] : 0.0)

                                                   + (j == k && l == m ? v.thirdDerivativeFactor * dr[i] : 0.0) +
                                                   (j == m && k == l ? v.thirdDerivativeFactor * dr[i] : 0.0) +
                                                   (j == l && k == m ? v.thirdDerivativeFactor * dr[i] : 0.0)

                                                   + (k == l && m == i ? v.thirdDerivativeFactor * dr[j] : 0.0) +
                                                   (k == i && l == m ? v.thirdDerivativeFactor * dr[j] : 0.0) +
                                                   (k == m && l == i ? v.thirdDerivativeFactor * dr[j] : 0.0)

                                                   + (l == m && i == j ? v.thirdDerivativeFactor * dr[k] : 0.0) +
                                                   (l == j && m == i ? v.thirdDerivativeFactor * dr[k] : 0.0) +
                                                   (l == i && m == j ? v.thirdDerivativeFactor * dr[k] : 0.0)

                                                   + (m == i && j == k ? v.thirdDerivativeFactor * dr[l] : 0.0) +
                                                   (m == k && i == j ? v.thirdDerivativeFactor * dr[l] : 0.0) +
                                                   (m == j && i == k ? v.thirdDerivativeFactor * dr[l] : 0.0)

                                                   + (i == j && k == l ? v.thirdDerivativeFactor * dr[m] : 0.0) +
                                                   (i == l && j == k ? v.thirdDerivativeFactor * dr[m] : 0.0) +
                                                   (i == k && j == l ? v.thirdDerivativeFactor * dr[m] : 0.0);

                for (size_t n = 0; n != 3; ++n)
                {
                  sixth_derivative[i][j][k][l][m][n] +=
                      v.sixthDerivativeFactor * dr[i] * dr[j] * dr[k] * dr[l] * dr[m] * dr[n]

                      + (i == j ? (v.fifthDerivativeFactor * dr[k] * dr[l] * dr[m] * dr[n]) : 0.0) +
                      (i == k ? (v.fifthDerivativeFactor * dr[j] * dr[l] * dr[m] * dr[n]) : 0.0) +
                      (i == l ? (v.fifthDerivativeFactor * dr[j] * dr[k] * dr[m] * dr[n]) : 0.0) +
                      (i == m ? (v.fifthDerivativeFactor * dr[j] * dr[k] * dr[l] * dr[n]) : 0.0) +
                      (i == n ? (v.fifthDerivativeFactor * dr[j] * dr[k] * dr[l] * dr[m]) : 0.0) +
                      (j == k ? (v.fifthDerivativeFactor * dr[i] * dr[l] * dr[m] * dr[n]) : 0.0) +
                      (j == l ? (v.fifthDerivativeFactor * dr[i] * dr[k] * dr[m] * dr[n]) : 0.0) +
                      (j == m ? (v.fifthDerivativeFactor * dr[i] * dr[k] * dr[l] * dr[n]) : 0.0) +
                      (j == n ? (v.fifthDerivativeFactor * dr[i] * dr[k] * dr[l] * dr[m]) : 0.0) +
                      (k == l ? (v.fifthDerivativeFactor * dr[i] * dr[j] * dr[m] * dr[n]) : 0.0) +
                      (k == m ? (v.fifthDerivativeFactor * dr[i] * dr[j] * dr[l] * dr[n]) : 0.0) +
                      (k == n ? (v.fifthDerivativeFactor * dr[i] * dr[j] * dr[l] * dr[m]) : 0.0) +
                      (l == m ? (v.fifthDerivativeFactor * dr[i] * dr[j] * dr[k] * dr[n]) : 0.0) +
                      (l == n ? (v.fifthDerivativeFactor * dr[i] * dr[j] * dr[k] * dr[m]) : 0.0) +
                      (m == n ? (v.fifthDerivativeFactor * dr[i] * dr[j] * dr[k] * dr[l]) : 0.0)

                      + (k == l && m == n ? v.fourthDerivativeFactor * dr[i] * dr[j] : 0.0) +
                      (k == m && l == n ? v.fourthDerivativeFactor * dr[i] * dr[j] : 0.0) +
                      (l == m && k == n ? v.fourthDerivativeFactor * dr[i] * dr[j] : 0.0)

                      + (j == l && m == n ? v.fourthDerivativeFactor * dr[i] * dr[k] : 0.0) +
                      (j == m && l == n ? v.fourthDerivativeFactor * dr[i] * dr[k] : 0.0) +
                      (j == n && l == m ? v.fourthDerivativeFactor * dr[i] * dr[k] : 0.0)

                      + (j == k && m == n ? v.fourthDerivativeFactor * dr[i] * dr[l] : 0.0) +
                      (j == m && k == n ? v.fourthDerivativeFactor * dr[i] * dr[l] : 0.0) +
                      (j == n && k == m ? v.fourthDerivativeFactor * dr[i] * dr[l] : 0.0)

                      + (j == k && l == n ? v.fourthDerivativeFactor * dr[i] * dr[m] : 0.0) +
                      (j == l && k == n ? v.fourthDerivativeFactor * dr[i] * dr[m] : 0.0) +
                      (j == n && k == l ? v.fourthDerivativeFactor * dr[i] * dr[m] : 0.0)

                      + (j == k && l == m ? v.fourthDerivativeFactor * dr[i] * dr[n] : 0.0) +
                      (j == l && k == m ? v.fourthDerivativeFactor * dr[i] * dr[n] : 0.0) +
                      (j == m && k == l ? v.fourthDerivativeFactor * dr[i] * dr[n] : 0.0)

                      + (i == l && m == n ? v.fourthDerivativeFactor * dr[j] * dr[k] : 0.0) +
                      (i == m && l == n ? v.fourthDerivativeFactor * dr[j] * dr[k] : 0.0) +
                      (i == n && l == m ? v.fourthDerivativeFactor * dr[j] * dr[k] : 0.0)

                      + (i == k && m == n ? v.fourthDerivativeFactor * dr[j] * dr[l] : 0.0) +
                      (i == m && k == n ? v.fourthDerivativeFactor * dr[j] * dr[l] : 0.0) +
                      (i == n && k == m ? v.fourthDerivativeFactor * dr[j] * dr[l] : 0.0)

                      + (i == k && l == n ? v.fourthDerivativeFactor * dr[j] * dr[m] : 0.0) +
                      (i == l && k == n ? v.fourthDerivativeFactor * dr[j] * dr[m] : 0.0) +
                      (i == n && k == l ? v.fourthDerivativeFactor * dr[j] * dr[m] : 0.0)

                      + (i == k && l == m ? v.fourthDerivativeFactor * dr[j] * dr[n] : 0.0) +
                      (i == l && k == m ? v.fourthDerivativeFactor * dr[j] * dr[n] : 0.0) +
                      (i == m && k == l ? v.fourthDerivativeFactor * dr[j] * dr[n] : 0.0)

                      + (i == j && m == n ? v.fourthDerivativeFactor * dr[k] * dr[l] : 0.0) +
                      (i == m && j == n ? v.fourthDerivativeFactor * dr[k] * dr[l] : 0.0) +
                      (i == n && j == m ? v.fourthDerivativeFactor * dr[k] * dr[l] : 0.0)

                      + (i == j && l == n ? v.fourthDerivativeFactor * dr[k] * dr[m] : 0.0) +
                      (i == l && j == n ? v.fourthDerivativeFactor * dr[k] * dr[m] : 0.0) +
                      (i == n && j == l ? v.fourthDerivativeFactor * dr[k] * dr[m] : 0.0)

                      + (i == j && l == m ? v.fourthDerivativeFactor * dr[k] * dr[n] : 0.0) +
                      (i == l && j == m ? v.fourthDerivativeFactor * dr[k] * dr[n] : 0.0) +
                      (i == m && j == l ? v.fourthDerivativeFactor * dr[k] * dr[n] : 0.0)

                      + (i == j && k == n ? v.fourthDerivativeFactor * dr[l] * dr[m] : 0.0) +
                      (i == k && j == n ? v.fourthDerivativeFactor * dr[l] * dr[m] : 0.0) +
                      (i == n && j == k ? v.fourthDerivativeFactor * dr[l] * dr[m] : 0.0)

                      + (i == j && k == m ? v.fourthDerivativeFactor * dr[l] * dr[n] : 0.0) +
                      (i == k && j == m ? v.fourthDerivativeFactor * dr[l] * dr[n] : 0.0) +
                      (i == m && j == k ? v.fourthDerivativeFactor * dr[l] * dr[n] : 0.0)

                      + (i == j && k == l ? v.fourthDerivativeFactor * dr[m] * dr[n] : 0.0) +
                      (i == k && j == l ? v.fourthDerivativeFactor * dr[m] * dr[n] : 0.0) +
                      (i == l && j == k ? v.fourthDerivativeFactor * dr[m] * dr[n] : 0.0)

                      + (i == j && k == l && m == n ? v.thirdDerivativeFactor : 0.0) +
                      (i == j && k == n && l == m ? v.thirdDerivativeFactor : 0.0) +
                      (i == j && k == m && l == n ? v.thirdDerivativeFactor : 0.0) +
                      (i == k && j == l && m == n ? v.thirdDerivativeFactor : 0.0) +
                      (i == k && j == m && l == n ? v.thirdDerivativeFactor : 0.0) +
                      (i == k && j == n && l == m ? v.thirdDerivativeFactor : 0.0) +
                      (i == l && j == k && m == n ? v.thirdDerivativeFactor : 0.0) +
                      (i == l && j == m && k == n ? v.thirdDerivativeFactor : 0.0) +
                      (i == l && j == n && k == m ? v.thirdDerivativeFactor : 0.0) +
                      (i == m && j == k && l == n ? v.thirdDerivativeFactor : 0.0) +
                      (i == m && j == l && k == n ? v.thirdDerivativeFactor : 0.0) +
                      (i == m && j == n && k == l ? v.thirdDerivativeFactor : 0.0) +
                      (i == n && j == k && l == m ? v.thirdDerivativeFactor : 0.0) +
                      (i == n && j == l && k == m ? v.thirdDerivativeFactor : 0.0) +
                      (i == n && j == m && k == l ? v.thirdDerivativeFactor : 0.0);
                }
              }
            }
          }
        }
      }
    }
  }

  return {energy,           first_derivative, second_derivative, third_derivative, fourth_derivative,
          fifth_derivative, sixth_derivative};
}

std::array<double, 27> Interactions::calculateTriquinticCartesianAtPosition(
    ForceField::InterpolationGridType interpolationGridType, const ForceField &forceField,
    const SimulationBox &simulationBox, double3 posA, size_t typeA, std::span<const Atom> frameworkAtoms)
{
  auto [energy, first_derivative, second_derivative, third_derivative, fourth_derivative, fifth_derivative,
        sixth_derivative] =
      Interactions::calculateTriquinticDerivativeAtPosition(interpolationGridType, forceField, simulationBox, posA,
                                                            typeA, frameworkAtoms);

  return std::array<double, 27>{energy,

                                first_derivative[0],
                                first_derivative[1],
                                first_derivative[2],

                                second_derivative[0][0],
                                second_derivative[0][1],
                                second_derivative[0][2],
                                second_derivative[1][1],
                                second_derivative[1][2],
                                second_derivative[2][2],

                                third_derivative[0][0][1],
                                third_derivative[0][0][2],
                                third_derivative[0][1][1],
                                third_derivative[0][1][2],
                                third_derivative[1][1][2],
                                third_derivative[0][2][2],
                                third_derivative[1][2][2],

                                fourth_derivative[0][0][1][1],
                                fourth_derivative[0][0][2][2],
                                fourth_derivative[1][1][2][2],
                                fourth_derivative[0][0][1][2],
                                fourth_derivative[0][1][1][2],
                                fourth_derivative[0][1][2][2],

                                fifth_derivative[0][0][1][1][2],
                                fifth_derivative[0][0][1][2][2],
                                fifth_derivative[0][1][1][2][2],

                                sixth_derivative[0][0][1][1][2][2]};
}

std::array<double, 27> Interactions::calculateTriquinticFractionalAtPosition(
    ForceField::InterpolationGridType interpolationGridType, const ForceField &forceField,
    const SimulationBox &simulationBox, double3 posA, size_t typeA, const SimulationBox &frameworkBox,
    std::span<const Atom> frameworkAtoms)
{
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  std::array<double, 27> derivatives{};
  double energy_fractional{0.0};
  std::array<double, 3> first_derivative_fractional{};
  std::array<std::array<double, 3>, 3> second_derivative_fractional{};
  std::array<std::array<std::array<double, 3>, 3>, 3> third_derivative_fractional{};
  std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> fourth_derivative_fractional{};
  std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3> fifth_derivative_fractional{};
  std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>
      sixth_derivative_fractional{};
  Potentials::TriquinticDerivativeFactor v{};
  double3x3 frameworkCell = frameworkBox.cell;

  switch (frameworkBox.type)
  {
    case SimulationBox::Type::Rectangular:
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        double3 posB = it1->position;
        double chargeB = it1->charge;
        size_t typeB = static_cast<size_t>(it1->type);

        double3 dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        double rr = std::max(0.1, double3::dot(dr, dr));

        if (rr < cutOffFrameworkVDWSquared)
        {
          switch (interpolationGridType)
          {
            case ForceField::InterpolationGridType::LennardJones:
              v = Potentials::potentialLennardJonesTriquinticDerivative(forceField, rr, typeA, typeB);
              break;
            case ForceField::InterpolationGridType::LennardJonesRepulsion:
              v = Potentials::potentialLennardJonesTriquinticDerivative(forceField, rr, typeA, typeB);
              break;
            case ForceField::InterpolationGridType::LennardJonesAttraction:
              v = Potentials::potentialLennardJonesTriquinticDerivative(forceField, rr, typeA, typeB);
              break;
            case ForceField::InterpolationGridType::EwaldReal:
              v = Potentials::potentialRealEwaldTriquinticDerivative(forceField, rr, std::sqrt(rr), 1.0, chargeB);
              break;
          }

          derivatives[0] += v.energy;

          derivatives[1] += v.firstDerivativeFactor * dr[0];
          derivatives[2] += v.firstDerivativeFactor * dr[1];
          derivatives[3] += v.firstDerivativeFactor * dr[2];

          derivatives[4] += v.secondDerivativeFactor * dr[0] * dr[0] + v.firstDerivativeFactor;
          derivatives[5] += v.secondDerivativeFactor * dr[0] * dr[1];
          derivatives[6] += v.secondDerivativeFactor * dr[0] * dr[2];
          derivatives[7] += v.secondDerivativeFactor * dr[1] * dr[1] + v.firstDerivativeFactor;
          derivatives[8] += v.secondDerivativeFactor * dr[1] * dr[2];
          derivatives[9] += v.secondDerivativeFactor * dr[2] * dr[2] + v.firstDerivativeFactor;

          derivatives[10] += v.thirdDerivativeFactor * dr[0] * dr[0] * dr[1] + v.secondDerivativeFactor * dr[1];
          derivatives[11] += v.thirdDerivativeFactor * dr[0] * dr[0] * dr[2] + v.secondDerivativeFactor * dr[2];
          derivatives[12] += v.thirdDerivativeFactor * dr[0] * dr[1] * dr[1] + v.secondDerivativeFactor * dr[0];
          derivatives[13] += v.thirdDerivativeFactor * dr[0] * dr[1] * dr[2];
          derivatives[14] += v.thirdDerivativeFactor * dr[1] * dr[1] * dr[2] + v.secondDerivativeFactor * dr[2];
          derivatives[15] += v.thirdDerivativeFactor * dr[0] * dr[2] * dr[2] + v.secondDerivativeFactor * dr[0];
          derivatives[16] += v.thirdDerivativeFactor * dr[1] * dr[2] * dr[2] + v.secondDerivativeFactor * dr[1];

          derivatives[17] += v.fourthDerivativeFactor * dr[0] * dr[0] * dr[1] * dr[1] +
                             v.thirdDerivativeFactor * dr[1] * dr[1] + v.thirdDerivativeFactor * dr[0] * dr[0] +
                             v.secondDerivativeFactor;
          derivatives[18] += v.fourthDerivativeFactor * dr[0] * dr[0] * dr[2] * dr[2] +
                             v.thirdDerivativeFactor * dr[2] * dr[2] + v.thirdDerivativeFactor * dr[0] * dr[0] +
                             v.secondDerivativeFactor;
          derivatives[19] += v.fourthDerivativeFactor * dr[1] * dr[1] * dr[2] * dr[2] +
                             v.thirdDerivativeFactor * dr[2] * dr[2] + v.thirdDerivativeFactor * dr[1] * dr[1] +
                             v.secondDerivativeFactor;
          derivatives[20] +=
              v.fourthDerivativeFactor * dr[0] * dr[0] * dr[1] * dr[2] + v.thirdDerivativeFactor * dr[1] * dr[2];
          derivatives[21] +=
              v.fourthDerivativeFactor * dr[0] * dr[1] * dr[1] * dr[2] + v.thirdDerivativeFactor * dr[0] * dr[2];
          derivatives[22] +=
              v.fourthDerivativeFactor * dr[0] * dr[1] * dr[2] * dr[2] + v.thirdDerivativeFactor * dr[0] * dr[1];

          derivatives[23] += v.fifthDerivativeFactor * dr[0] * dr[0] * dr[1] * dr[1] * dr[2] +
                             v.fourthDerivativeFactor * dr[1] * dr[1] * dr[2] +
                             v.fourthDerivativeFactor * dr[0] * dr[0] * dr[2] + v.thirdDerivativeFactor * dr[2];

          derivatives[24] += v.fifthDerivativeFactor * dr[0] * dr[0] * dr[1] * dr[2] * dr[2] +
                             v.fourthDerivativeFactor * dr[1] * dr[2] * dr[2] +
                             v.fourthDerivativeFactor * dr[0] * dr[0] * dr[1] + v.thirdDerivativeFactor * dr[1];

          derivatives[25] += v.fifthDerivativeFactor * dr[0] * dr[1] * dr[1] * dr[2] * dr[2] +
                             v.fourthDerivativeFactor * dr[0] * dr[2] * dr[2] +
                             v.fourthDerivativeFactor * dr[0] * dr[1] * dr[1] + v.thirdDerivativeFactor * dr[0];

          derivatives[26] += v.sixthDerivativeFactor * dr[0] * dr[0] * dr[1] * dr[1] * dr[2] * dr[2] +
                             v.fifthDerivativeFactor * dr[1] * dr[1] * dr[2] * dr[2] +
                             v.fifthDerivativeFactor * dr[0] * dr[0] * dr[2] * dr[2] +
                             v.fifthDerivativeFactor * dr[0] * dr[0] * dr[1] * dr[1] +
                             v.fourthDerivativeFactor * dr[0] * dr[0] + v.fourthDerivativeFactor * dr[1] * dr[1] +
                             v.fourthDerivativeFactor * dr[2] * dr[2];
        }
      }

      return {derivatives[0],

              derivatives[1] * frameworkCell[0][0],
              derivatives[2] * frameworkCell[0][0],
              derivatives[3] * frameworkCell[0][0],

              derivatives[4] * frameworkCell[0][0] * frameworkCell[0][0],
              derivatives[5] * frameworkCell[0][0] * frameworkCell[1][1],
              derivatives[6] * frameworkCell[0][0] * frameworkCell[2][2],
              derivatives[7] * frameworkCell[1][1] * frameworkCell[1][1],
              derivatives[8] * frameworkCell[1][1] * frameworkCell[2][2],
              derivatives[9] * frameworkCell[2][2] * frameworkCell[2][2],

              derivatives[10] * frameworkCell[0][0] * frameworkCell[0][0] * frameworkCell[1][1],
              derivatives[11] * frameworkCell[0][0] * frameworkCell[0][0] * frameworkCell[2][2],
              derivatives[12] * frameworkCell[0][0] * frameworkCell[1][1] * frameworkCell[1][1],
              derivatives[13] * frameworkCell[0][0] * frameworkCell[1][1] * frameworkCell[2][2],
              derivatives[14] * frameworkCell[1][1] * frameworkCell[1][1] * frameworkCell[2][2],
              derivatives[15] * frameworkCell[0][0] * frameworkCell[2][2] * frameworkCell[2][2],
              derivatives[16] * frameworkCell[1][1] * frameworkCell[2][2] * frameworkCell[2][2],

              derivatives[17] * frameworkCell[0][0] * frameworkCell[0][0] * frameworkCell[1][1] * frameworkCell[1][1],
              derivatives[18] * frameworkCell[0][0] * frameworkCell[0][0] * frameworkCell[2][2] * frameworkCell[2][2],
              derivatives[19] * frameworkCell[1][1] * frameworkCell[1][1] * frameworkCell[2][2] * frameworkCell[2][2],
              derivatives[20] * frameworkCell[0][0] * frameworkCell[0][0] * frameworkCell[1][1] * frameworkCell[2][2],
              derivatives[21] * frameworkCell[0][0] * frameworkCell[1][1] * frameworkCell[1][1] * frameworkCell[2][2],
              derivatives[22] * frameworkCell[0][0] * frameworkCell[1][1] * frameworkCell[2][2] * frameworkCell[2][2],

              derivatives[23] * frameworkCell[0][0] * frameworkCell[0][0] * frameworkCell[1][1] * frameworkCell[1][1] *
                  frameworkCell[2][2],
              derivatives[24] * frameworkCell[0][0] * frameworkCell[0][0] * frameworkCell[1][1] * frameworkCell[2][2] *
                  frameworkCell[2][2],
              derivatives[25] * frameworkCell[0][0] * frameworkCell[1][1] * frameworkCell[1][1] * frameworkCell[2][2] *
                  frameworkCell[2][2],

              derivatives[26] * frameworkCell[0][0] * frameworkCell[0][0] * frameworkCell[1][1] * frameworkCell[1][1] *
                  frameworkCell[2][2] * frameworkCell[2][2]};
    case SimulationBox::Type::Triclinic:
      auto [energy, first_derivative, second_derivative, third_derivative, fourth_derivative, fifth_derivative,
            sixth_derivative] =
          Interactions::calculateTriquinticDerivativeAtPosition(interpolationGridType, forceField, simulationBox, posA,
                                                                typeA, frameworkAtoms);

      energy_fractional = energy;
      for (const size_t &p : std::array<size_t, 3>{0, 1, 2})
      {
        first_derivative_fractional[0] += frameworkCell[0][p] * first_derivative[p];
        first_derivative_fractional[1] += frameworkCell[1][p] * first_derivative[p];
        first_derivative_fractional[2] += frameworkCell[2][p] * first_derivative[p];

        for (const size_t &q : std::array<size_t, 3>{0, 1, 2})
        {
          second_derivative_fractional[0][0] += frameworkCell[0][p] * frameworkCell[0][q] * second_derivative[p][q];
          second_derivative_fractional[0][1] += frameworkCell[0][p] * frameworkCell[1][q] * second_derivative[p][q];
          second_derivative_fractional[0][2] += frameworkCell[0][p] * frameworkCell[2][q] * second_derivative[p][q];
          second_derivative_fractional[1][1] += frameworkCell[1][p] * frameworkCell[1][q] * second_derivative[p][q];
          second_derivative_fractional[1][2] += frameworkCell[1][p] * frameworkCell[2][q] * second_derivative[p][q];
          second_derivative_fractional[2][2] += frameworkCell[2][p] * frameworkCell[2][q] * second_derivative[p][q];

          for (const size_t &r : std::array<size_t, 3>{0, 1, 2})
          {
            third_derivative_fractional[0][0][1] +=
                frameworkCell[0][p] * frameworkCell[0][q] * frameworkCell[1][r] * third_derivative[p][q][r];
            third_derivative_fractional[0][0][2] +=
                frameworkCell[0][p] * frameworkCell[0][q] * frameworkCell[2][r] * third_derivative[p][q][r];
            third_derivative_fractional[0][1][1] +=
                frameworkCell[0][p] * frameworkCell[1][q] * frameworkCell[1][r] * third_derivative[p][q][r];
            third_derivative_fractional[0][1][2] +=
                frameworkCell[0][p] * frameworkCell[1][q] * frameworkCell[2][r] * third_derivative[p][q][r];
            third_derivative_fractional[1][1][2] +=
                frameworkCell[1][p] * frameworkCell[1][q] * frameworkCell[2][r] * third_derivative[p][q][r];
            third_derivative_fractional[0][2][2] +=
                frameworkCell[0][p] * frameworkCell[2][q] * frameworkCell[2][r] * third_derivative[p][q][r];
            third_derivative_fractional[1][2][2] +=
                frameworkCell[1][p] * frameworkCell[2][q] * frameworkCell[2][r] * third_derivative[p][q][r];

            for (const size_t &s : std::array<size_t, 3>{0, 1, 2})
            {
              fourth_derivative_fractional[0][0][1][1] += frameworkCell[0][p] * frameworkCell[0][q] *
                                                          frameworkCell[1][r] * frameworkCell[1][s] *
                                                          fourth_derivative[p][q][r][s];
              fourth_derivative_fractional[0][0][2][2] += frameworkCell[0][p] * frameworkCell[0][q] *
                                                          frameworkCell[2][r] * frameworkCell[2][s] *
                                                          fourth_derivative[p][q][r][s];
              fourth_derivative_fractional[1][1][2][2] += frameworkCell[1][p] * frameworkCell[1][q] *
                                                          frameworkCell[2][r] * frameworkCell[2][s] *
                                                          fourth_derivative[p][q][r][s];
              fourth_derivative_fractional[0][0][1][2] += frameworkCell[0][p] * frameworkCell[0][q] *
                                                          frameworkCell[1][r] * frameworkCell[2][s] *
                                                          fourth_derivative[p][q][r][s];
              fourth_derivative_fractional[0][1][1][2] += frameworkCell[0][p] * frameworkCell[1][q] *
                                                          frameworkCell[1][r] * frameworkCell[2][s] *
                                                          fourth_derivative[p][q][r][s];
              fourth_derivative_fractional[0][1][2][2] += frameworkCell[0][p] * frameworkCell[1][q] *
                                                          frameworkCell[2][r] * frameworkCell[2][s] *
                                                          fourth_derivative[p][q][r][s];

              for (const size_t &t : std::array<size_t, 3>{0, 1, 2})
              {
                fifth_derivative_fractional[0][0][1][1][2] += frameworkCell[0][p] * frameworkCell[0][q] *
                                                              frameworkCell[1][r] * frameworkCell[1][s] *
                                                              frameworkCell[2][t] * fifth_derivative[p][q][r][s][t];
                fifth_derivative_fractional[0][0][1][2][2] += frameworkCell[0][p] * frameworkCell[0][q] *
                                                              frameworkCell[1][r] * frameworkCell[2][s] *
                                                              frameworkCell[2][t] * fifth_derivative[p][q][r][s][t];
                fifth_derivative_fractional[0][1][1][2][2] += frameworkCell[0][p] * frameworkCell[1][q] *
                                                              frameworkCell[1][r] * frameworkCell[2][s] *
                                                              frameworkCell[2][t] * fifth_derivative[p][q][r][s][t];

                for (const size_t &u : std::array<size_t, 3>{0, 1, 2})
                {
                  sixth_derivative_fractional[0][0][1][1][2][2] +=
                      frameworkCell[0][p] * frameworkCell[0][q] * frameworkCell[1][r] * frameworkCell[1][s] *
                      frameworkCell[2][t] * frameworkCell[2][u] * sixth_derivative[p][q][r][s][t][u];
                }
              }
            }
          }
        }
      }

      return {energy_fractional,

              first_derivative_fractional[0],
              first_derivative_fractional[1],
              first_derivative_fractional[2],

              second_derivative_fractional[0][0],
              second_derivative_fractional[0][1],
              second_derivative_fractional[0][2],
              second_derivative_fractional[1][1],
              second_derivative_fractional[1][2],
              second_derivative_fractional[2][2],

              third_derivative_fractional[0][0][1],
              third_derivative_fractional[0][0][2],
              third_derivative_fractional[0][1][1],
              third_derivative_fractional[0][1][2],
              third_derivative_fractional[1][1][2],
              third_derivative_fractional[0][2][2],
              third_derivative_fractional[1][2][2],

              fourth_derivative_fractional[0][0][1][1],
              fourth_derivative_fractional[0][0][2][2],
              fourth_derivative_fractional[1][1][2][2],
              fourth_derivative_fractional[0][0][1][2],
              fourth_derivative_fractional[0][1][1][2],
              fourth_derivative_fractional[0][1][2][2],

              fifth_derivative_fractional[0][0][1][1][2],
              fifth_derivative_fractional[0][0][1][2][2],
              fifth_derivative_fractional[0][1][1][2][2],

              sixth_derivative_fractional[0][0][1][1][2][2]};
  };
}
