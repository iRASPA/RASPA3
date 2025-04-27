module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <map>
#include <ostream>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>
#endif

module integrators_cputime;

#ifndef USE_LEGACY_HEADERS
import <chrono>;
import <string>;
import <sstream>;
import <format>;
import <exception>;
import <source_location>;
import <fstream>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <algorithm>;
import <print>;
#endif

import double3;
import stringutils;
import archive;
import json;

const std::string IntegratorsCPUTime::writeIntegratorsCPUTimeStatistics() const
{
  std::ostringstream stream;

  if (computeTranslationalKineticEnergy > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Compute Translational Kinetic Energy: {:14f} [s]\n", computeTranslationalKineticEnergy.count());
  }

  if (computeRotationalKineticEnergy > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Compute Rotational Kinetic Energy:    {:14f} [s]\n", computeRotationalKineticEnergy.count());
  }

  if (computeCenterOfMass > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Compute Center of Mass:               {:14f} [s]\n", computeCenterOfMass.count());
  }

  if (computeCenterOfMassVelocity > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Compute Center of Mass Velocity:      {:14f} [s]\n", computeCenterOfMassVelocity.count());
  }

  if (computeLinearMomentum > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Compute Linear Momentum:              {:14f} [s]\n", computeLinearMomentum.count());
  }

  if (scaleVelocities > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Scale Velocities:                     {:14f} [s]\n", scaleVelocities.count());
  }

  if (removeCenterOfMassVelocity > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Remove COM velocity:                  {:14f} [s]\n", removeCenterOfMassVelocity.count());
  }

  if (updatePositions > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Update Positions:                     {:14f} [s]\n", updatePositions.count());
  }

  if (updateVelocities > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Update Velocities:                    {:14f} [s]\n", updateVelocities.count());
  }

  if (createCartesianPositions > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Create Cartesian Positions:           {:14f} [s]\n", createCartesianPositions.count());
  }

  if (noSquishFreeRotorOrderTwo > std::chrono::duration<double>::zero())
  {
    std::print(stream, "No Squish Free Rotor Order Two:       {:14f} [s]\n", noSquishFreeRotorOrderTwo.count());
  }

  if (noSquishRotate > std::chrono::duration<double>::zero())
  {
    std::print(stream, "No Squish Rotate:                     {:14f} [s]\n", noSquishRotate.count());
  }

  if (updateCenterOfMassAndQuaternionVelocities > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Update CoM and Quaternion Velocities: {:14f} [s]\n",
               updateCenterOfMassAndQuaternionVelocities.count());
  }

  if (updateCenterOfMassAndQuaternionGradients > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Update CoM and Quaternion Gradients:  {:14f} [s]\n",
               updateCenterOfMassAndQuaternionGradients.count());
  }

  if (updateGradients > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Update Gradients:                     {:14f} [s]\n", updateGradients.count());
  }

  if (velocityVerlet > std::chrono::duration<double>::zero())
  {
    std::print(stream, "Velocity Verlet:                      {:14f} [s]\n", velocityVerlet.count());
  }

  return stream.str();
}

const std::string IntegratorsCPUTime::writeIntegratorsCPUTimeStatistics(
    std::chrono::duration<double> totalSimulation) const
{
  std::ostringstream stream;

  std::print(stream, "Integrators CPU Time Statistics:\n\n");

  std::print(stream, "Compute Translational Kinetic Energy: {:14f} [s]\n", computeTranslationalKineticEnergy.count());
  std::print(stream, "Compute Rotational Kinetic Energy:    {:14f} [s]\n", computeRotationalKineticEnergy.count());
  std::print(stream, "Compute Center of Mass:               {:14f} [s]\n", computeCenterOfMass.count());
  std::print(stream, "Compute Center of Mass Velocity:      {:14f} [s]\n", computeCenterOfMassVelocity.count());
  std::print(stream, "Compute Linear Momentum:              {:14f} [s]\n", computeLinearMomentum.count());
  std::print(stream, "Scale Velocities:                     {:14f} [s]\n", scaleVelocities.count());
  std::print(stream, "Scale Velocities:                     {:14f} [s]\n", removeCenterOfMassVelocity.count());
  std::print(stream, "Update Positions:                     {:14f} [s]\n", updatePositions.count());
  std::print(stream, "Update Velocities:                    {:14f} [s]\n", updateVelocities.count());
  std::print(stream, "Create Cartesian Positions:           {:14f} [s]\n", createCartesianPositions.count());
  std::print(stream, "No Squish Free Rotor Order Two:       {:14f} [s]\n", noSquishFreeRotorOrderTwo.count());
  std::print(stream, "No Squish Rotate:                     {:14f} [s]\n", noSquishRotate.count());
  std::print(stream, "Update CoM and Quaternion Velocities: {:14f} [s]\n",
             updateCenterOfMassAndQuaternionVelocities.count());
  std::print(stream, "Update CoM and Quaternion Gradients:  {:14f} [s]\n",
             updateCenterOfMassAndQuaternionGradients.count());
  std::print(stream, "Update Gradients:                     {:14f} [s]\n", updateGradients.count());
  std::print(stream, "Velocity Verlet:                      {:14f} [s]\n", velocityVerlet.count());

  std::print(stream, "Total Simulation Time:                {:14f} [s]\n", totalSimulation.count());
  std::print(stream, "Overhead:                             {:14f} [s]\n",
             totalSimulation.count() - velocityVerlet.count());

  return stream.str();
}

const nlohmann::json IntegratorsCPUTime::jsonSystemIntegratorsCPUTimeStatistics() const
{
  nlohmann::json status;

  status["computeTranslationalKineticEnergy"] = computeTranslationalKineticEnergy.count();
  status["computeRotationalKineticEnergy"] = computeRotationalKineticEnergy.count();
  status["computeCenterOfMass"] = computeCenterOfMass.count();
  status["computeCenterOfMassVelocity"] = computeCenterOfMassVelocity.count();
  status["computeLinearMomentum"] = computeLinearMomentum.count();
  status["scaleVelocities"] = scaleVelocities.count();
  status["removeCenterOfMassVelocity"] = removeCenterOfMassVelocity.count();
  status["updatePositions"] = updatePositions.count();
  status["updateVelocities"] = updateVelocities.count();
  status["createCartesianPositions"] = createCartesianPositions.count();
  status["noSquishFreeRotorOrderTwo"] = noSquishFreeRotorOrderTwo.count();
  status["noSquishRotate"] = noSquishRotate.count();
  status["updateCenterOfMassAndQuaternionVelocities"] = updateCenterOfMassAndQuaternionVelocities.count();
  status["updateCenterOfMassAndQuaternionGradients"] = updateCenterOfMassAndQuaternionGradients.count();
  status["updateGradients"] = updateGradients.count();
  status["velocityVerlet"] = velocityVerlet.count();

  return status;
}

const nlohmann::json IntegratorsCPUTime::jsonOverallIntegratorsCPUTimeStatistics(
    std::chrono::duration<double> totalSimulation) const
{
  nlohmann::json status = jsonSystemIntegratorsCPUTimeStatistics();

  status["totalIntegratorsCPUTime"] = total().count();
  status["totalSimulationTime"] = totalSimulation.count();
  status["overhead"] = totalSimulation.count() - total().count();

  return status;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const IntegratorsCPUTime& t)
{
  archive << t.versionNumber;

  archive << t.computeTranslationalKineticEnergy;
  archive << t.computeRotationalKineticEnergy;
  archive << t.computeCenterOfMass;
  archive << t.computeCenterOfMassVelocity;
  archive << t.computeLinearMomentum;
  archive << t.scaleVelocities;
  archive << t.removeCenterOfMassVelocity;
  archive << t.updatePositions;
  archive << t.updateVelocities;
  archive << t.createCartesianPositions;
  archive << t.noSquishFreeRotorOrderTwo;
  archive << t.noSquishRotate;
  archive << t.updateCenterOfMassAndQuaternionVelocities;
  archive << t.updateCenterOfMassAndQuaternionGradients;
  archive << t.updateGradients;
  archive << t.velocityVerlet;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, IntegratorsCPUTime& t)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > t.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'IntegratorsCPUTime' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> t.computeTranslationalKineticEnergy;
  archive >> t.computeRotationalKineticEnergy;
  archive >> t.computeCenterOfMass;
  archive >> t.computeCenterOfMassVelocity;
  archive >> t.computeLinearMomentum;
  archive >> t.scaleVelocities;
  archive >> t.updatePositions;
  archive >> t.updateVelocities;
  archive >> t.removeCenterOfMassVelocity;
  archive >> t.createCartesianPositions;
  archive >> t.noSquishFreeRotorOrderTwo;
  archive >> t.noSquishRotate;
  archive >> t.updateCenterOfMassAndQuaternionVelocities;
  archive >> t.updateCenterOfMassAndQuaternionGradients;
  archive >> t.updateGradients;
  archive >> t.velocityVerlet;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("IntegratorsCPUTime: Error in binary restart\n"));
  }
#endif

  return archive;
}
