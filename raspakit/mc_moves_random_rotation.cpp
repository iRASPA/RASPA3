module;

module mc_moves;

import component;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_factor;
import energy_status;
import energy_status_inter;
import running_energy;
import lambda;
import property_widom;
import averages;

import <complex>;
import <vector>;
import <array>;
import <tuple>;
import <optional>;
import <span>;
import <optional>;
import <tuple>;
import <algorithm>;
import <chrono>;
import <cmath>;
import <iostream>;
import <iomanip>;

std::optional<RunningEnergy> MC_Moves::randomRotationMove(System& system, size_t selectedComponent, std::span<Atom> molecule)
{
    double3 angle{};
    std::array<double3,3> axes{double3(1.0,0.0,0.0), double3(0.0,1.0,0.0) ,double3(0.0,0.0,1.0) };
    double3 maxAngle = system.components[selectedComponent].statistics_RandomRotationMove.maxChange;
    size_t selectedDirection = size_t(3.0 * RandomNumber::Uniform());
    angle[selectedDirection] = maxAngle[selectedDirection] * 2.0 * (RandomNumber::Uniform() - 0.5);
    system.components[selectedComponent].statistics_RandomRotationMove.counts[selectedDirection] += 1;

    size_t startingBead = system.components[selectedComponent].startingBead;
    std::vector<Atom> trialPositions(molecule.size());
    double rotationAngle = angle[selectedDirection];
    double3 rotationAxis = double3(axes[selectedDirection]);
    double3x3 rotationMatrix = double3x3(simd_quatd::fromAxisAngle(rotationAngle, rotationAxis));
    std::transform(molecule.begin(), molecule.end(), trialPositions.begin(),
            [&](Atom a) { a.position = rotationMatrix * (a.position - molecule[startingBead].position) 
                          + molecule[startingBead].position; return a; });
    std::span<Atom> newMolecule{trialPositions.begin(), trialPositions.end()};

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkMolecule = system.computeFrameworkMoleculeEnergyDifference(newMolecule, molecule);
    if (!frameworkMolecule.has_value()) return std::nullopt;

    std::optional<RunningEnergy> interMolecule = system.computeInterMolecularEnergyDifference(newMolecule, molecule);
    if (!interMolecule.has_value()) return std::nullopt;
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    system.components[selectedComponent].cpuTime_RandomRotationMove_NonEwald += (t2 - t1);

    std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
    RunningEnergy ewaldFourierEnergy = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, molecule);
    std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
    system.components[selectedComponent].cpuTime_RandomRotationMove_Ewald += (u2 - u1);
    RunningEnergy energyDifference = frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy;

    system.components[selectedComponent].statistics_RandomRotationMove.constructed[selectedDirection] += 1;

    if (RandomNumber::Uniform() < std::exp(-system.Beta * energyDifference.total()))
    {
        system.components[selectedComponent].statistics_RandomRotationMove.accepted[selectedDirection] += 1;

        system.acceptEwaldMove();
        std::copy(trialPositions.cbegin(), trialPositions.cend(), molecule.begin());

        return energyDifference;
    };
    return std::nullopt;
}