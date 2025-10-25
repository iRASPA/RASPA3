module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#endif

export module triquintic_derivatives_external_field;

#ifdef USE_STD_IMPORT
import std;
#endif

import double3;

export namespace Potentials
{
std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
potentialTriquinticSecondOrderPolynomialTestFunction(const double3 &pos);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
potentialTriquinticThirdOrderPolynomialTestFunction(const double3 &pos);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
potentialTriquinticFourthOrderPolynomialTestFunction(const double3 &pos);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
potentialTriquinticFifthOrderPolynomialTestFunction(const double3 &pos);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
potentialTriquinticSixthOrderPolynomialTestFunction(const double3 &pos);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
potentialTriquinticExponentialNonPolynomialTestFunction(const double3 &pos);
}  // namespace Potentials
