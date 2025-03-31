module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <tuple>
#endif

export module tricubic_derivatives_external_field;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <array>;
import <tuple>;
import <iostream>;
#endif

import double3;
import tricubic_derivative_factor;

export namespace Potentials
{
std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
potentialTricubicSecondOrderPolynomialTestFunction(const double3 &pos);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
potentialTricubicThirdOrderPolynomialTestFunction(const double3 &pos);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
potentialTricubicFourthOrderPolynomialTestFunction(const double3 &pos);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
potentialTricubicFifthOrderPolynomialTestFunction(const double3 &pos);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
potentialTricubicSixthOrderPolynomialTestFunction(const double3 &pos);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
potentialTricubicExponentialNonPolynomialTestFunction(const double3 &pos);
}  // namespace Potentials
