module;

export module interactions_external_field_grid;

import std;

import double3;
import double3x3;
import atom;
import running_energy;
import energy_status;
import simulationbox;
import gradient_factor;
import forcefield;
import component;

export namespace Interactions
{
std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
calculateThirdDerivativeAtPositionExternalField(const ForceField &forceField, const SimulationBox &simulationBox,
                                                double3 posA);

std::array<double, 8> calculateTricubicFractionalAtPositionExternalField(const ForceField &forceField,
                                                                         const SimulationBox &simulationBox,
                                                                         double3 posA);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
calculateSixthOrderDerivativeAtPositionExternalField(const ForceField &forceField, double3 pos);

std::array<double, 27> calculateTriquinticFractionalAtPositionExternalField(const ForceField &forceField,
                                                                            const SimulationBox &simulationBox,
                                                                            double3 posA);

}  // namespace Interactions
