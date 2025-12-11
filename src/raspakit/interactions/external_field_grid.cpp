module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <optional>
#include <span>
#include <type_traits>
#include <array>
#include <vector>
#endif

module interactions_external_field_grid;

#ifdef USE_STD_IMPORT
import std;
#endif

import int3;
import double3;
import double4;
import double3x3;
import atom;
import simulationbox;
import energy_status;
import energy_status_inter;
import units;
import energy_factor;
import gradient_factor;
import running_energy;
import component;
import forcefield;
import tricubic_derivatives_external_field;
import triquintic_derivatives_external_field;


std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
Interactions::calculateThirdDerivativeAtPositionExternalField(const ForceField &forceField,
                                                              [[maybe_unused]] const SimulationBox &simulationBox,
                                                              double3 pos)
{
  double energy{0.0};
  std::array<double, 3> first_derivative{};
  std::array<std::array<double, 3>, 3> second_derivative{};
  std::array<std::array<std::array<double, 3>, 3>, 3> third_derivative{};

  //if (forceField.hasExternalField)
  {
    switch (forceField.potentialEnergySurfaceType)
    {
      case ForceField::PotentialEnergySurfaceType::None:
        break;
      case ForceField::PotentialEnergySurfaceType::SecondOrderPolynomialTestFunction:
        return Potentials::potentialTricubicSecondOrderPolynomialTestFunction(pos);
      case ForceField::PotentialEnergySurfaceType::ThirdOrderPolynomialTestFunction:
        return Potentials::potentialTricubicThirdOrderPolynomialTestFunction(pos);
      case ForceField::PotentialEnergySurfaceType::FourthOrderPolynomialTestFunction:
        return Potentials::potentialTricubicFourthOrderPolynomialTestFunction(pos);
      case ForceField::PotentialEnergySurfaceType::FifthOrderPolynomialTestFunction:
        return Potentials::potentialTricubicFifthOrderPolynomialTestFunction(pos);
      case ForceField::PotentialEnergySurfaceType::SixthOrderPolynomialTestFunction:
        return Potentials::potentialTricubicSixthOrderPolynomialTestFunction(pos);
      case ForceField::PotentialEnergySurfaceType::ExponentialNonPolynomialTestFunction:
        return Potentials::potentialTricubicExponentialNonPolynomialTestFunction(pos);
      default:
        break;
    }
  }

  return {energy, first_derivative, second_derivative, third_derivative};
}

std::array<double, 8> Interactions::calculateTricubicFractionalAtPositionExternalField(
    const ForceField &forceField, [[maybe_unused]] const SimulationBox &simulationBox, double3 posA)
{
  //if (forceField.hasExternalField)
  {
    double energy_fractional{0.0};
    std::array<double, 3> first_derivative_fractional{};
    std::array<std::array<double, 3>, 3> second_derivative_fractional{};
    std::array<std::array<std::array<double, 3>, 3>, 3> third_derivative_fractional{};

    auto [energy, first_derivative, second_derivative, third_derivative] =
        Interactions::calculateThirdDerivativeAtPositionExternalField(forceField, simulationBox, posA);

    energy_fractional = energy;
    for (const std::size_t &p : std::array<std::size_t, 3>{0, 1, 2})
    {
      first_derivative_fractional[0] += simulationBox.cell[0][p] * first_derivative[p];
      first_derivative_fractional[1] += simulationBox.cell[1][p] * first_derivative[p];
      first_derivative_fractional[2] += simulationBox.cell[2][p] * first_derivative[p];
      for (const std::size_t &q : std::array<std::size_t, 3>{0, 1, 2})
      {
        second_derivative_fractional[0][1] +=
            simulationBox.cell[0][p] * simulationBox.cell[1][q] * second_derivative[p][q];
        second_derivative_fractional[0][2] +=
            simulationBox.cell[0][p] * simulationBox.cell[2][q] * second_derivative[p][q];
        second_derivative_fractional[1][2] +=
            simulationBox.cell[1][p] * simulationBox.cell[2][q] * second_derivative[p][q];
        for (const std::size_t &r : std::array<std::size_t, 3>{0, 1, 2})
        {
          third_derivative_fractional[0][1][2] += simulationBox.cell[0][p] * simulationBox.cell[1][q] *
                                                  simulationBox.cell[2][r] * third_derivative[p][q][r];
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

  return {};
}

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
Interactions::calculateSixthOrderDerivativeAtPositionExternalField(const ForceField &forceField, double3 pos)
{
  //if (forceField.hasExternalField)
  {
    switch (forceField.potentialEnergySurfaceType)
    {
      case ForceField::PotentialEnergySurfaceType::None:
        break;
      case ForceField::PotentialEnergySurfaceType::SecondOrderPolynomialTestFunction:
        return Potentials::potentialTriquinticSecondOrderPolynomialTestFunction(pos);
      case ForceField::PotentialEnergySurfaceType::ThirdOrderPolynomialTestFunction:
        return Potentials::potentialTriquinticThirdOrderPolynomialTestFunction(pos);
      case ForceField::PotentialEnergySurfaceType::FourthOrderPolynomialTestFunction:
        return Potentials::potentialTriquinticFourthOrderPolynomialTestFunction(pos);
      case ForceField::PotentialEnergySurfaceType::FifthOrderPolynomialTestFunction:
        return Potentials::potentialTriquinticFifthOrderPolynomialTestFunction(pos);
      case ForceField::PotentialEnergySurfaceType::SixthOrderPolynomialTestFunction:
        return Potentials::potentialTriquinticSixthOrderPolynomialTestFunction(pos);
      case ForceField::PotentialEnergySurfaceType::ExponentialNonPolynomialTestFunction:
        return Potentials::potentialTriquinticExponentialNonPolynomialTestFunction(pos);
      default:
        break;
    }
  }

  return {};
}

std::array<double, 27> Interactions::calculateTriquinticFractionalAtPositionExternalField(
    const ForceField &forceField, const SimulationBox &simulationBox, double3 posA)
{
  //if (forceField.hasExternalField)
  {
    double energy_fractional{0.0};
    std::array<double, 3> first_derivative_fractional{};
    std::array<std::array<double, 3>, 3> second_derivative_fractional{};
    std::array<std::array<std::array<double, 3>, 3>, 3> third_derivative_fractional{};
    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> fourth_derivative_fractional{};
    std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3> fifth_derivative_fractional{};
    std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>
        sixth_derivative_fractional{};

    auto [energy, first_derivative, second_derivative, third_derivative, fourth_derivative, fifth_derivative,
          sixth_derivative] = Interactions::calculateSixthOrderDerivativeAtPositionExternalField(forceField, posA);

    energy_fractional = energy;
    for (const std::size_t &p : std::array<std::size_t, 3>{0, 1, 2})
    {
      first_derivative_fractional[0] += simulationBox.cell[0][p] * first_derivative[p];
      first_derivative_fractional[1] += simulationBox.cell[1][p] * first_derivative[p];
      first_derivative_fractional[2] += simulationBox.cell[2][p] * first_derivative[p];

      for (const std::size_t &q : std::array<std::size_t, 3>{0, 1, 2})
      {
        second_derivative_fractional[0][0] +=
            simulationBox.cell[0][p] * simulationBox.cell[0][q] * second_derivative[p][q];
        second_derivative_fractional[0][1] +=
            simulationBox.cell[0][p] * simulationBox.cell[1][q] * second_derivative[p][q];
        second_derivative_fractional[0][2] +=
            simulationBox.cell[0][p] * simulationBox.cell[2][q] * second_derivative[p][q];
        second_derivative_fractional[1][1] +=
            simulationBox.cell[1][p] * simulationBox.cell[1][q] * second_derivative[p][q];
        second_derivative_fractional[1][2] +=
            simulationBox.cell[1][p] * simulationBox.cell[2][q] * second_derivative[p][q];
        second_derivative_fractional[2][2] +=
            simulationBox.cell[2][p] * simulationBox.cell[2][q] * second_derivative[p][q];

        for (const std::size_t &r : std::array<std::size_t, 3>{0, 1, 2})
        {
          third_derivative_fractional[0][0][1] += simulationBox.cell[0][p] * simulationBox.cell[0][q] *
                                                  simulationBox.cell[1][r] * third_derivative[p][q][r];
          third_derivative_fractional[0][0][2] += simulationBox.cell[0][p] * simulationBox.cell[0][q] *
                                                  simulationBox.cell[2][r] * third_derivative[p][q][r];
          third_derivative_fractional[0][1][1] += simulationBox.cell[0][p] * simulationBox.cell[1][q] *
                                                  simulationBox.cell[1][r] * third_derivative[p][q][r];
          third_derivative_fractional[0][1][2] += simulationBox.cell[0][p] * simulationBox.cell[1][q] *
                                                  simulationBox.cell[2][r] * third_derivative[p][q][r];
          third_derivative_fractional[1][1][2] += simulationBox.cell[1][p] * simulationBox.cell[1][q] *
                                                  simulationBox.cell[2][r] * third_derivative[p][q][r];
          third_derivative_fractional[0][2][2] += simulationBox.cell[0][p] * simulationBox.cell[2][q] *
                                                  simulationBox.cell[2][r] * third_derivative[p][q][r];
          third_derivative_fractional[1][2][2] += simulationBox.cell[1][p] * simulationBox.cell[2][q] *
                                                  simulationBox.cell[2][r] * third_derivative[p][q][r];

          for (const std::size_t &s : std::array<std::size_t, 3>{0, 1, 2})
          {
            fourth_derivative_fractional[0][0][1][1] += simulationBox.cell[0][p] * simulationBox.cell[0][q] *
                                                        simulationBox.cell[1][r] * simulationBox.cell[1][s] *
                                                        fourth_derivative[p][q][r][s];
            fourth_derivative_fractional[0][0][2][2] += simulationBox.cell[0][p] * simulationBox.cell[0][q] *
                                                        simulationBox.cell[2][r] * simulationBox.cell[2][s] *
                                                        fourth_derivative[p][q][r][s];
            fourth_derivative_fractional[1][1][2][2] += simulationBox.cell[1][p] * simulationBox.cell[1][q] *
                                                        simulationBox.cell[2][r] * simulationBox.cell[2][s] *
                                                        fourth_derivative[p][q][r][s];
            fourth_derivative_fractional[0][0][1][2] += simulationBox.cell[0][p] * simulationBox.cell[0][q] *
                                                        simulationBox.cell[1][r] * simulationBox.cell[2][s] *
                                                        fourth_derivative[p][q][r][s];
            fourth_derivative_fractional[0][1][1][2] += simulationBox.cell[0][p] * simulationBox.cell[1][q] *
                                                        simulationBox.cell[1][r] * simulationBox.cell[2][s] *
                                                        fourth_derivative[p][q][r][s];
            fourth_derivative_fractional[0][1][2][2] += simulationBox.cell[0][p] * simulationBox.cell[1][q] *
                                                        simulationBox.cell[2][r] * simulationBox.cell[2][s] *
                                                        fourth_derivative[p][q][r][s];

            for (const std::size_t &t : std::array<std::size_t, 3>{0, 1, 2})
            {
              fifth_derivative_fractional[0][0][1][1][2] += simulationBox.cell[0][p] * simulationBox.cell[0][q] *
                                                            simulationBox.cell[1][r] * simulationBox.cell[1][s] *
                                                            simulationBox.cell[2][t] * fifth_derivative[p][q][r][s][t];
              fifth_derivative_fractional[0][0][1][2][2] += simulationBox.cell[0][p] * simulationBox.cell[0][q] *
                                                            simulationBox.cell[1][r] * simulationBox.cell[2][s] *
                                                            simulationBox.cell[2][t] * fifth_derivative[p][q][r][s][t];
              fifth_derivative_fractional[0][1][1][2][2] += simulationBox.cell[0][p] * simulationBox.cell[1][q] *
                                                            simulationBox.cell[1][r] * simulationBox.cell[2][s] *
                                                            simulationBox.cell[2][t] * fifth_derivative[p][q][r][s][t];

              for (const std::size_t &u : std::array<std::size_t, 3>{0, 1, 2})
              {
                sixth_derivative_fractional[0][0][1][1][2][2] += simulationBox.cell[0][p] * simulationBox.cell[0][q] *
                                                                 simulationBox.cell[1][r] * simulationBox.cell[1][s] *
                                                                 simulationBox.cell[2][t] * simulationBox.cell[2][u] *
                                                                 sixth_derivative[p][q][r][s][t][u];
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
  }

  return {};
}
