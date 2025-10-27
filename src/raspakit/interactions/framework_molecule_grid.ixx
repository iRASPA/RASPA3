module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module interactions_framework_molecule_grid;

#ifdef USE_STD_IMPORT
import std;
#endif

import double3;
import double3x3;
import double3x3x3;
import atom;
import running_energy;
import energy_status;
import simulationbox;
import energy_factor;
import gradient_factor;
import hessian_factor;
import forcefield;
import framework;
import component;

export namespace Interactions
{
double calculateEnergyAtPosition(ForceField::InterpolationGridType interpolationGridType, const ForceField &forceField,
                                 const SimulationBox &simulationBox, double3 posB, std::size_t typeB,
                                 std::span<const Atom> frameworkAtoms);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
calculateTricubicDerivativeAtPosition(ForceField::InterpolationGridType interpolationGridType,
                                      const ForceField &forceField, const SimulationBox &simulationBox, double3 posB,
                                      std::size_t typeB, std::span<const Atom> frameworkAtoms);

std::array<double, 8> calculateTricubicCartesianAtPosition(ForceField::InterpolationGridType interpolationGridType,
                                                           const ForceField &forceField,
                                                           const SimulationBox &simulationBox, double3 posA,
                                                           std::size_t typeA, std::span<const Atom> frameworkAtoms);

std::array<double, 8> calculateTricubicFractionalAtPosition(ForceField::InterpolationGridType interpolationGridType,
                                                            const ForceField &forceField,
                                                            const SimulationBox &simulationBox, double3 posA,
                                                            std::size_t typeA, const SimulationBox &frameworkBox,
                                                            std::span<const Atom> frameworkAtoms);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
calculateTriquinticDerivativeAtPosition(ForceField::InterpolationGridType interpolationGridType,
                                        const ForceField &forceField, const SimulationBox &simulationBox, double3 posA,
                                        std::size_t typeA, std::span<const Atom> frameworkAtoms);

std::array<double, 27> calculateTriquinticCartesianAtPosition(ForceField::InterpolationGridType interpolationGridType,
                                                              const ForceField &forceField,
                                                              const SimulationBox &simulationBox, double3 posA,
                                                              std::size_t typeA, std::span<const Atom> frameworkAtoms);

std::array<double, 27> calculateTriquinticFractionalAtPosition(ForceField::InterpolationGridType interpolationGridType,
                                                               const ForceField &forceField,
                                                               const SimulationBox &simulationBox, double3 posA,
                                                               std::size_t typeA, const SimulationBox &frameworkBox,
                                                               std::span<const Atom> frameworkAtoms);

};  // namespace Interactions
