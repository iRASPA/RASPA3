module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <cstddef>
#include <fstream>
#include <string>
#endif

export module integrators_cputime;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <chrono>;
import <fstream>;
#endif

import double3;
import archive;
import json;

export struct IntegratorsCPUTime
{
  // Constructor
  IntegratorsCPUTime() {};

  // Version number
  uint64_t versionNumber{1};

  // Timings
  std::chrono::duration<double> computeTranslationalKineticEnergy{0.0};
  std::chrono::duration<double> computeRotationalKineticEnergy{0.0};
  std::chrono::duration<double> computeCenterOfMass{0.0};
  std::chrono::duration<double> computeCenterOfMassVelocity{0.0};
  std::chrono::duration<double> computeLinearMomentum{0.0};
  std::chrono::duration<double> scaleVelocities{0.0};
  std::chrono::duration<double> removeCenterOfMassVelocity{0.0};
  std::chrono::duration<double> updatePositions{0.0};
  std::chrono::duration<double> updateVelocities{0.0};
  std::chrono::duration<double> createCartesianPositions{0.0};
  std::chrono::duration<double> noSquishFreeRotorOrderTwo{0.0};
  std::chrono::duration<double> noSquishRotate{0.0};
  std::chrono::duration<double> updateCenterOfMassAndQuaternionVelocities{0.0};
  std::chrono::duration<double> updateCenterOfMassAndQuaternionGradients{0.0};
  std::chrono::duration<double> updateGradients{0.0};
  std::chrono::duration<double> velocityVerlet{0.0};

  // Total timing
  inline std::chrono::duration<double> total() const
  {
    return computeTranslationalKineticEnergy + computeRotationalKineticEnergy + computeCenterOfMass +
           computeCenterOfMassVelocity + computeLinearMomentum + scaleVelocities + removeCenterOfMassVelocity +
           updatePositions + updateVelocities + createCartesianPositions + noSquishFreeRotorOrderTwo + noSquishRotate +
           updateCenterOfMassAndQuaternionVelocities + updateCenterOfMassAndQuaternionGradients + updateGradients +
           velocityVerlet;
  }

  void clearTimingStatistics()
  {
    computeTranslationalKineticEnergy = std::chrono::duration<double>{0.0};
    computeRotationalKineticEnergy = std::chrono::duration<double>{0.0};
    computeCenterOfMass = std::chrono::duration<double>{0.0};
    computeCenterOfMassVelocity = std::chrono::duration<double>{0.0};
    computeLinearMomentum = std::chrono::duration<double>{0.0};
    scaleVelocities = std::chrono::duration<double>{0.0};
    removeCenterOfMassVelocity = std::chrono::duration<double>{0.0};
    updatePositions = std::chrono::duration<double>{0.0};
    updateVelocities = std::chrono::duration<double>{0.0};
    createCartesianPositions = std::chrono::duration<double>{0.0};
    noSquishFreeRotorOrderTwo = std::chrono::duration<double>{0.0};
    noSquishRotate = std::chrono::duration<double>{0.0};
    updateCenterOfMassAndQuaternionVelocities = std::chrono::duration<double>{0.0};
    updateCenterOfMassAndQuaternionGradients = std::chrono::duration<double>{0.0};
    updateGradients = std::chrono::duration<double>{0.0};
    velocityVerlet = std::chrono::duration<double>{0.0};
  }

  const std::string writeIntegratorsCPUTimeStatistics() const;
  const std::string writeIntegratorsCPUTimeStatistics(std::chrono::duration<double> total) const;

  const nlohmann::json jsonSystemIntegratorsCPUTimeStatistics() const;
  const nlohmann::json jsonOverallIntegratorsCPUTimeStatistics(std::chrono::duration<double> total) const;

  // Default copy constructor
  IntegratorsCPUTime(const IntegratorsCPUTime&) = default;

  // operator=
  inline IntegratorsCPUTime& operator=(const IntegratorsCPUTime& b)
  {
    computeTranslationalKineticEnergy = b.computeTranslationalKineticEnergy;
    computeRotationalKineticEnergy = b.computeRotationalKineticEnergy;
    computeCenterOfMass = b.computeCenterOfMass;
    computeCenterOfMassVelocity = b.computeCenterOfMassVelocity;
    computeLinearMomentum = b.computeLinearMomentum;
    scaleVelocities = b.scaleVelocities;
    updatePositions = b.updatePositions;
    updateVelocities = b.updateVelocities;
    createCartesianPositions = b.createCartesianPositions;
    noSquishFreeRotorOrderTwo = b.noSquishFreeRotorOrderTwo;
    noSquishRotate = b.noSquishRotate;
    updateCenterOfMassAndQuaternionVelocities = b.updateCenterOfMassAndQuaternionVelocities;
    updateCenterOfMassAndQuaternionGradients = b.updateCenterOfMassAndQuaternionGradients;
    updateGradients = b.updateGradients;
    velocityVerlet = b.velocityVerlet;

    return *this;
  }

  // operator+=
  inline IntegratorsCPUTime& operator+=(const IntegratorsCPUTime& b)
  {
    computeTranslationalKineticEnergy += b.computeTranslationalKineticEnergy;
    computeRotationalKineticEnergy += b.computeRotationalKineticEnergy;
    computeCenterOfMass += b.computeCenterOfMass;
    computeCenterOfMassVelocity += b.computeCenterOfMassVelocity;
    computeLinearMomentum += b.computeLinearMomentum;
    scaleVelocities += b.scaleVelocities;
    removeCenterOfMassVelocity += b.removeCenterOfMassVelocity;
    updatePositions += b.updatePositions;
    updateVelocities += b.updateVelocities;
    createCartesianPositions += b.createCartesianPositions;
    noSquishFreeRotorOrderTwo += b.noSquishFreeRotorOrderTwo;
    noSquishRotate += b.noSquishRotate;
    updateCenterOfMassAndQuaternionVelocities += b.updateCenterOfMassAndQuaternionVelocities;
    updateCenterOfMassAndQuaternionGradients += b.updateCenterOfMassAndQuaternionGradients;
    updateGradients += b.updateGradients;
    velocityVerlet += b.velocityVerlet;

    return *this;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const IntegratorsCPUTime& t);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, IntegratorsCPUTime& t);
};

export inline IntegratorsCPUTime operator+(const IntegratorsCPUTime& a, const IntegratorsCPUTime& b)
{
  IntegratorsCPUTime m;

  m.computeTranslationalKineticEnergy = a.computeTranslationalKineticEnergy + b.computeTranslationalKineticEnergy;
  m.computeRotationalKineticEnergy = a.computeRotationalKineticEnergy + b.computeRotationalKineticEnergy;
  m.computeCenterOfMass = a.computeCenterOfMass + b.computeCenterOfMass;
  m.computeCenterOfMassVelocity = a.computeCenterOfMassVelocity + b.computeCenterOfMassVelocity;
  m.computeLinearMomentum = a.computeLinearMomentum + b.computeLinearMomentum;
  m.scaleVelocities = a.scaleVelocities + b.scaleVelocities;
  m.removeCenterOfMassVelocity = a.removeCenterOfMassVelocity + b.removeCenterOfMassVelocity;
  m.updatePositions = a.updatePositions + b.updatePositions;
  m.updateVelocities = a.updateVelocities + b.updateVelocities;
  m.createCartesianPositions = a.createCartesianPositions + b.createCartesianPositions;
  m.noSquishFreeRotorOrderTwo = a.noSquishFreeRotorOrderTwo + b.noSquishFreeRotorOrderTwo;
  m.noSquishRotate = a.noSquishRotate + b.noSquishRotate;
  m.updateCenterOfMassAndQuaternionVelocities =
      a.updateCenterOfMassAndQuaternionVelocities + b.updateCenterOfMassAndQuaternionVelocities;
  m.updateCenterOfMassAndQuaternionGradients =
      a.updateCenterOfMassAndQuaternionGradients + b.updateCenterOfMassAndQuaternionGradients;
  m.updateGradients = a.updateGradients + b.updateGradients;
  m.velocityVerlet = a.velocityVerlet + b.velocityVerlet;

  return m;
}

/**
 * \brief NOTE: integratorsCPUTime is now defined as global and therefore defined "program-wide".
 * In its current implementation it can not be tracked per system, which might be necessary in the future.
 */
export IntegratorsCPUTime integratorsCPUTime;
