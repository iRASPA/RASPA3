module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <deque>
#include <functional>
#include <future>
#include <iostream>
#include <numbers>
#include <optional>
#include <print>
#include <semaphore>
#include <span>
#include <thread>
#include <type_traits>
#include <vector>
#endif

module cbmc_interactions_framework_molecule;

#ifdef USE_STD_IMPORT
import std;
#endif

import energy_status;
import potential_energy_vdw;
import potential_gradient_vdw;
import potential_energy_coulomb;
import potential_gradient_coulomb;
import potential_correction_vdw;
import framework;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import energy_factor;
import gradient_factor;
import energy_status_inter;
import running_energy;
import units;
import threadpool;
import interpolation_energy_grid;

template <ThreadPool::ThreadingType T>
[[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms, double cutOffVDW,
    double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) noexcept;

template <>
[[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy<ThreadPool::ThreadingType::Serial>(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    [[maybe_unused]] const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms, double cutOffVDW,
    double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) noexcept
{
  bool useCharge = forceField.useCharge;
  [[maybe_unused]] const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffVDWSquared = cutOffVDW * cutOffVDW;
  const double cutOffChargeSquared = cutOffCoulomb * cutOffCoulomb;

  RunningEnergy energySum;

  for (int index = 0; const Atom &atom : atoms)
  {
    if (index != skip)
    {
      double3 posB = atom.position;
      std::size_t typeB = static_cast<std::size_t>(atom.type);
      bool groupIdB = static_cast<bool>(atom.groupId);
      bool isFractional = static_cast<bool>(atom.isFractional);
      double scalingVDWB = atom.scalingVDW;
      double scalingCoulombB = atom.scalingCoulomb;
      double chargeB = atom.charge;

      if (interpolationGrids[typeB].has_value() && !isFractional)
      {
        double energy = interpolationGrids[typeB]->interpolate(posB);
        if (energy > overlapCriteria)
        {
          return std::nullopt;
        }
        energySum.frameworkMoleculeVDW += energy;
        if (useCharge)
        {
          energySum.frameworkMoleculeCharge += chargeB * interpolationGrids.back()->interpolate(posB);
        }
      }
      else
      {
        for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
        {
          double3 posA = it1->position;
          std::size_t typeA = static_cast<std::size_t>(it1->type);
          bool groupIdA = static_cast<bool>(it1->groupId);
          double scalingVDWA = it1->scalingVDW;
          double scalingCoulombA = it1->scalingCoulomb;
          double chargeA = it1->charge;

          double3 dr = posA - posB;
          dr = simulationBox.applyPeriodicBoundaryConditions(dr);
          double rr = double3::dot(dr, dr);

          if (rr < cutOffVDWSquared)
          {
            Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
                forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
            if (energyFactor.energy > overlapCriteria)
            {
              return std::nullopt;
            }
            energySum.frameworkMoleculeVDW += energyFactor.energy;
            energySum.dudlambdaVDW += energyFactor.dUdlambda;
          }
          if (useCharge && rr < cutOffChargeSquared)
          {
            double r = std::sqrt(rr);
            Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
                forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

            energySum.frameworkMoleculeCharge += energyFactor.energy;
            energySum.dudlambdaCharge += energyFactor.dUdlambda;
          }
        }
      }
    }
    ++index;
  }
  return energySum;
}

template <>
[[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy<ThreadPool::ThreadingType::ThreadPool>(
    const ForceField &forceField, const SimulationBox &simulationBox,
    [[maybe_unused]] const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    [[maybe_unused]] const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms, double cutOffVDW,
    double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) noexcept
{
  std::atomic_flag cancel;

  auto &pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
  const std::size_t numberOfHelperThreads = pool.getThreadCount();

  std::vector<std::future<RunningEnergy>> threads(numberOfHelperThreads);

  // std::size_t const block_size = frameworkAtoms.size() / numberOfHelperThreads;
  std::size_t const block_size = frameworkAtoms.size() / (numberOfHelperThreads + 1);

  auto task = [skip, cutOffVDW, cutOffCoulomb, atoms, &cancel, &forceField, &simulationBox](
                  std::span<const Atom>::iterator startIterator,
                  std::span<const Atom>::iterator endIterator) -> RunningEnergy
  {
    RunningEnergy energySum;

    bool useCharge = forceField.useCharge;
    const double overlapCriteria = forceField.energyOverlapCriteria;
    const double cutOffVDWSquared = cutOffVDW * cutOffVDW;
    const double cutOffChargeSquared = cutOffCoulomb * cutOffCoulomb;

    for (std::span<const Atom>::iterator it1 = startIterator; it1 != endIterator; ++it1)
    {
      if (cancel.test()) return energySum;
      double3 posA = it1->position;
      std::size_t typeA = static_cast<std::size_t>(it1->type);
      bool groupIdA = static_cast<bool>(it1->groupId);
      double scalingVDWA = it1->scalingVDW;
      double scalingCoulombA = it1->scalingCoulomb;
      double chargeA = it1->charge;

      for (int index = 0; const Atom &atom : atoms)
      {
        if (index != skip)
        {
          double3 posB = atom.position;
          std::size_t typeB = static_cast<std::size_t>(atom.type);
          bool groupIdB = static_cast<bool>(atom.groupId);
          double scalingVDWB = atom.scalingVDW;
          double scalingCoulombB = atom.scalingCoulomb;
          double chargeB = atom.charge;

          double3 dr = posA - posB;
          dr = simulationBox.applyPeriodicBoundaryConditions(dr);
          double rr = double3::dot(dr, dr);

          if (rr < cutOffVDWSquared)
          {
            Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
                forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);

            if (energyFactor.energy > overlapCriteria)
            {
              cancel.test_and_set();
              return energySum;
            }
            energySum.frameworkMoleculeVDW += energyFactor.energy;
            energySum.dudlambdaVDW += energyFactor.dUdlambda;
          }
          if (useCharge && rr < cutOffChargeSquared)
          {
            double r = std::sqrt(rr);
            Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
                forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

            energySum.frameworkMoleculeCharge += energyFactor.energy;
            energySum.dudlambdaCharge += energyFactor.dUdlambda;
          }
        }
        ++index;
      }
    }

    return energySum;
  };
  /*
    std::span<const Atom>::iterator block_start = frameworkAtoms.begin();
    for(std::size_t i = 0 ; i != numberOfHelperThreads - 1; ++i)
    {
      std::span<const Atom>::iterator block_end = block_start;
      std::advance(block_end,block_size);

      threads[i] = pool.enqueue(task, block_start, block_end);
      block_start=block_end;
    }
    threads[numberOfHelperThreads - 1] = pool.enqueue(task, block_start, frameworkAtoms.end());

    RunningEnergy energy{};
    for(std::size_t i = 0; i != numberOfHelperThreads; ++i)
    {
      energy += threads[i].get();
    }
    if(cancel.test()) return std::nullopt;
    */

  std::span<const Atom>::iterator block_start = frameworkAtoms.begin();
  for (std::size_t i = 0; i != numberOfHelperThreads; ++i)
  {
    std::span<const Atom>::iterator block_end = block_start;
    std::advance(block_end, block_size);

    threads[i] = pool.enqueue(task, block_start, block_end);
    block_start = block_end;
  }
  RunningEnergy energy = task(block_start, frameworkAtoms.end());

  for (std::size_t i = 0; i != numberOfHelperThreads; ++i)
  {
    energy += threads[i].get();
  }
  if (cancel.test()) return std::nullopt;

  return energy;
}

template <>
[[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy<ThreadPool::ThreadingType::OpenMP>(
    const ForceField &forceField, const SimulationBox &simulationBox,
    [[maybe_unused]] const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    [[maybe_unused]] const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms, double cutOffVDW,
    double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) noexcept
{
  bool useCharge = forceField.useCharge;
  [[maybe_unused]] const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffVDWSquared = cutOffVDW * cutOffVDW;
  const double cutOffChargeSquared = cutOffCoulomb * cutOffCoulomb;

  std::atomic_flag cancel;
  cancel.clear();

  RunningEnergy energySum;

#pragma omp declare reduction \
   (energy_sum_reduction : RunningEnergy : omp_out += omp_in) \
   initializer(omp_priv = RunningEnergy{})

#pragma omp parallel for reduction(energy_sum_reduction : energySum)
  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    if (!cancel.test())
    {
      double3 posA = it1->position;
      std::size_t typeA = static_cast<std::size_t>(it1->type);
      bool groupIdA = static_cast<bool>(it1->groupId);
      double scalingVDWA = it1->scalingVDW;
      double scalingCoulombA = it1->scalingCoulomb;
      double chargeA = it1->charge;

      for (int index = 0; const Atom &atom : atoms)
      {
        if (index != skip)
        {
          double3 posB = atom.position;
          std::size_t typeB = static_cast<std::size_t>(atom.type);
          bool groupIdB = static_cast<bool>(atom.groupId);
          double scalingVDWB = atom.scalingVDW;
          double scalingCoulombB = atom.scalingCoulomb;
          double chargeB = atom.charge;

          double3 dr = posA - posB;
          dr = simulationBox.applyPeriodicBoundaryConditions(dr);
          double rr = double3::dot(dr, dr);

          if (rr < cutOffVDWSquared)
          {
            Potentials::EnergyFactor energyFactor = Potentials::potentialVDWEnergy(
                forceField, groupIdA, groupIdB, scalingVDWA, scalingVDWB, rr, typeA, typeB);
            if (energyFactor.energy > overlapCriteria)
            {
              cancel.test_and_set();
            }
            energySum.frameworkMoleculeVDW += energyFactor.energy;
            energySum.dudlambdaVDW += energyFactor.dUdlambda;
          }
          if (useCharge && rr < cutOffChargeSquared)
          {
            double r = std::sqrt(rr);
            Potentials::EnergyFactor energyFactor = Potentials::potentialCoulombEnergy(
                forceField, groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

            energySum.frameworkMoleculeCharge += energyFactor.energy;
            energySum.dudlambdaCharge += energyFactor.dUdlambda;
          }
        }
        ++index;
      }
    }
  }
  if (cancel.test()) return std::nullopt;
  return energySum;
}

template <>
[[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy<ThreadPool::ThreadingType::GPU_Offload>(
    [[maybe_unused]] const ForceField &forceField, [[maybe_unused]] const SimulationBox &simulationBox,
    [[maybe_unused]] const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    [[maybe_unused]] const std::optional<Framework> &framework, [[maybe_unused]] std::span<const Atom> frameworkAtoms,
    [[maybe_unused]] double cutOffVDW, [[maybe_unused]] double cutOffCoulomb, [[maybe_unused]] std::span<Atom> atoms,
    [[maybe_unused]] std::make_signed_t<std::size_t> skip) noexcept
{
  return std::nullopt;
}

[[nodiscard]] std::optional<RunningEnergy> CBMC::computeFrameworkMoleculeEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms, double cutOffVDW,
    double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) noexcept
{
  auto &pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
  switch (pool.getThreadingType())
  {
    default:
    case ThreadPool::ThreadingType::Serial:
    {
      return computeFrameworkMoleculeEnergy<ThreadPool::ThreadingType::Serial>(
          forceField, simulationBox, interpolationGrids, framework, frameworkAtoms, cutOffVDW, cutOffCoulomb, atoms,
          skip);
    }
    case ThreadPool::ThreadingType::OpenMP:
    {
      return computeFrameworkMoleculeEnergy<ThreadPool::ThreadingType::OpenMP>(
          forceField, simulationBox, interpolationGrids, framework, frameworkAtoms, cutOffVDW, cutOffCoulomb, atoms,
          skip);
    }
    case ThreadPool::ThreadingType::ThreadPool:
    {
      return computeFrameworkMoleculeEnergy<ThreadPool::ThreadingType::ThreadPool>(
          forceField, simulationBox, interpolationGrids, framework, frameworkAtoms, cutOffVDW, cutOffCoulomb, atoms,
          skip);
    }
    case ThreadPool::ThreadingType::GPU_Offload:
    {
      return computeFrameworkMoleculeEnergy<ThreadPool::ThreadingType::GPU_Offload>(
          forceField, simulationBox, interpolationGrids, framework, frameworkAtoms, cutOffVDW, cutOffCoulomb, atoms,
          skip);
    }
  }
}
