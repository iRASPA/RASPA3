module;

module cbmc_external_energy;

import std;

import atom;
import running_energy;
import framework;
import component;
import double3;
import double4;
import forcefield;
import simulationbox;
import threadpool;
import interpolation_energy_grid;
import interactions_pair_kernel;
import potential_pair_derivatives;
import cbmc_grow_context;
import cbmc_results;

// ---------------------------------------------------------------------------------------------------
// Blocking pockets.
// ---------------------------------------------------------------------------------------------------

bool CBMC::insideBlockedPockets(const std::optional<Framework> &framework, const Component &component,
                                std::span<const Atom> molecule_atoms)
{
  if (framework.has_value())
  {
    for (std::size_t i = 0; i != component.blockingPockets.size(); ++i)
    {
      double radius_squared = component.blockingPockets[i].w * component.blockingPockets[i].w;
      double3 pos =
          framework->simulationBox.cell *
          double3(component.blockingPockets[i].x, component.blockingPockets[i].y, component.blockingPockets[i].z);
      for (const Atom &atom : molecule_atoms)
      {
        double3 dr = atom.position - pos;

        // compute the periodic boundary conditions with the single unit cell of the framework
        dr = framework->simulationBox.applyPeriodicBoundaryConditions(dr);

        double vdwScaling = atom.scalingVDW;
        if (dr.length_squared() < vdwScaling * radius_squared)
        {
          return true;
        }
      }
    }
  }
  return false;
}

// ---------------------------------------------------------------------------------------------------
// External field.
// ---------------------------------------------------------------------------------------------------

[[nodiscard]] std::optional<RunningEnergy> CBMC::computeExternalFieldEnergy(
    bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid, [[maybe_unused]] double cutOffVDW,
    [[maybe_unused]] double cutOffCoulomb, std::span<const Atom> atoms) noexcept
{
  RunningEnergy energySum{};

  if (!hasExternalField) return energySum;

  const double4 externalFieldGeometryParameters = forceField.externalFieldGeometryParameters;

  // A confining geometry (cylinder / rectangle along an axis) contributes zero energy inside and an
  // overlap outside: the distance of 'posA' to the axis through the cell centre along 'axisDirection'.
  const auto insideCylinder = [&](const double3 &posA, const double3 &axisBegin, const double3 &axisEnd)
  {
    const double3 v = (axisEnd - axisBegin).normalized();
    const double3 w = double3::cross(posA - axisBegin, v);
    return w.length_squared() < externalFieldGeometryParameters.x * externalFieldGeometryParameters.x;
  };
  // 'first' and 'second' index the two transverse components (0 = x, 1 = y, 2 = z) compared against the
  // half-widths in 'externalFieldGeometryParameters.x' and '.y'.
  const auto insideRectangle = [&](const double3 &posA, const double3 &axisBegin, const double3 &axisEnd,
                                   std::size_t first, std::size_t second)
  {
    const double3 v = (axisEnd - axisBegin).normalized();
    const double3 w = double3::cross(posA - axisBegin, v);
    return std::abs(w[first]) < externalFieldGeometryParameters.x &&
           std::abs(w[second]) < externalFieldGeometryParameters.y;
  };

  for (const Atom &atom : atoms)
  {
    const std::uint8_t groupIdA = atom.groupId;
    const double scalingVDWA = atom.scalingVDW;
    const double3 posA = atom.position;

    Potentials::PairDerivatives<0> energyFactor{0.0, 0.0};

    if (externalFieldInterpolationGrid.has_value())
    {
      energyFactor.energy = scalingVDWA * externalFieldInterpolationGrid->interpolate(posA);
    }
    else
    {
      switch (forceField.potentialEnergySurfaceType)
      {
        case ForceField::PotentialEnergySurfaceType::None:
          break;
        case ForceField::PotentialEnergySurfaceType::GridFile:
          energyFactor.energy = scalingVDWA * externalFieldInterpolationGrid->interpolate(posA);
          break;
        case ForceField::PotentialEnergySurfaceType::MullerBrown:
        {
          const double4 A{-200.0, -100.0, -170.0, -15.0};
          const double4 a{-1.0, -1.0, -6.5, -0.7};
          const double4 x0{1.0, 0.0, -0.5, -1.0};
          const double4 b{0.0, 0.0, 11.0, 0.6};
          const double4 y0{0.0, 0.5, 1.5, 1.0};
          const double4 c{-10.0, -10.0, -6.5, 0.7};
          energyFactor.energy = 0.0;
          for (std::size_t i = 0; i < 4; ++i)
          {
            energyFactor.energy += A[i] * std::exp(a[i] * (posA.x - x0[i]) * (posA.x - x0[i]) +
                                                   b[i] * (posA.x - x0[i]) * (posA.y - y0[i]) +
                                                   c[i] * (posA.y - y0[i]) * (posA.y - y0[i]));
          }
        }
        break;
        case ForceField::PotentialEnergySurfaceType::ThirdOrderPolynomialTestFunction:
          energyFactor.energy = posA.x * posA.y * posA.z;
          break;
        case ForceField::PotentialEnergySurfaceType::CylinderX:
          if (!insideCylinder(posA, simulationBox.cell * double3{0.0, 0.5, 0.5},
                              simulationBox.cell * double3{1.0, 0.5, 0.5}))
            return std::nullopt;
          break;
        case ForceField::PotentialEnergySurfaceType::CylinderY:
          if (!insideCylinder(posA, simulationBox.cell * double3{0.5, 0.0, 0.5},
                              simulationBox.cell * double3{0.5, 1.0, 0.5}))
            return std::nullopt;
          break;
        case ForceField::PotentialEnergySurfaceType::CylinderZ:
          if (!insideCylinder(posA, simulationBox.cell * double3{0.5, 0.5, 0.0},
                              simulationBox.cell * double3{0.5, 0.5, 1.0}))
            return std::nullopt;
          break;
        case ForceField::PotentialEnergySurfaceType::RectangleX:
          if (!insideRectangle(posA, simulationBox.cell * double3{0.0, 0.5, 0.5},
                               simulationBox.cell * double3{1.0, 0.5, 0.5}, 1, 2))
            return std::nullopt;
          break;
        case ForceField::PotentialEnergySurfaceType::RectangleY:
          if (!insideRectangle(posA, simulationBox.cell * double3{0.5, 0.0, 0.5},
                               simulationBox.cell * double3{0.5, 1.0, 0.5}, 2, 0))
            return std::nullopt;
          break;
        case ForceField::PotentialEnergySurfaceType::RectangleZ:
          if (!insideRectangle(posA, simulationBox.cell * double3{0.5, 0.5, 0.0},
                               simulationBox.cell * double3{0.5, 0.5, 1.0}, 1, 0))
            return std::nullopt;
          break;
        default:
          break;
      }
    }
    energySum.externalFieldVDW += energyFactor.energy;
    energySum.addDudlambdaVDW(groupIdA, 0, 1.0, 1.0, energyFactor.dUdlambda);
  }

  return energySum;
}

// ---------------------------------------------------------------------------------------------------
// Framework-molecule.
// ---------------------------------------------------------------------------------------------------

namespace
{
/// Adds one framework-molecule pair to 'energySum'. False on a hard overlap (the VDW energy exceeds
/// 'energyOverlapCriteria'), in which case 'energySum' is meaningless and the caller abandons the trial.
bool accumulateFrameworkPair(const ForceField &forceField, const SimulationBox &simulationBox,
                             const Atom &frameworkAtom, const Atom &atom, double cutOffVDWSquared,
                             double cutOffChargeSquared, RunningEnergy &energySum)
{
  bool overlap = false;
  Interactions::evaluatePair<0>(
      forceField, simulationBox, frameworkAtom, atom, cutOffVDWSquared, cutOffChargeSquared, forceField.useCharge,
      [&](const Potentials::PairDerivatives<0> &factors, const double3 &)
      {
        if (factors.energy > forceField.energyOverlapCriteria)
        {
          overlap = true;
          return;
        }
        energySum.frameworkMoleculeVDW += factors.energy;
        energySum.addDudlambdaVDW(frameworkAtom.groupId, atom.groupId, frameworkAtom.scalingVDW, atom.scalingVDW,
                                  factors.dUdlambda);
      },
      [&](const Potentials::PairDerivatives<0> &factors, const double3 &)
      {
        if (overlap) return;
        energySum.frameworkMoleculeCharge += factors.energy;
        energySum.addDudlambdaCharge(frameworkAtom.groupId, atom.groupId, frameworkAtom.scalingCoulomb,
                                     atom.scalingCoulomb, factors.dUdlambda);
      });
  return !overlap;
}

/// The energy of one trial atom from its interpolation grid (VDW, and the Coulomb grid when charges
/// are on); false when the grid energy signals an overlap.
bool accumulateGridEnergy(const ForceField &forceField,
                          const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
                          const Atom &atom, RunningEnergy &energySum)
{
  const std::size_t type = static_cast<std::size_t>(atom.type);
  const double energy = interpolationGrids[type]->interpolate(atom.position);
  if (energy > forceField.energyOverlapCriteria) return false;

  energySum.frameworkMoleculeVDW += energy;
  if (forceField.useCharge)
  {
    energySum.frameworkMoleculeCharge += atom.charge * interpolationGrids.back()->interpolate(atom.position);
  }
  return true;
}

[[nodiscard]] std::optional<RunningEnergy> computeFrameworkSerial(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids, std::span<const Atom> frameworkAtoms,
    double cutOffVDW, double cutOffCoulomb, std::span<const Atom> atoms) noexcept
{
  const double cutOffVDWSquared = cutOffVDW * cutOffVDW;
  const double cutOffChargeSquared = cutOffCoulomb * cutOffCoulomb;

  RunningEnergy energySum;
  for (const Atom &atom : atoms)
  {
    // A grid exists only for a framework; a fractional atom is never read from the grid.
    if (interpolationGrids[static_cast<std::size_t>(atom.type)].has_value() && !atom.isFractional)
    {
      if (!accumulateGridEnergy(forceField, interpolationGrids, atom, energySum)) return std::nullopt;
      continue;
    }

    for (const Atom &frameworkAtom : frameworkAtoms)
    {
      if (!accumulateFrameworkPair(forceField, simulationBox, frameworkAtom, atom, cutOffVDWSquared,
                                   cutOffChargeSquared, energySum))
      {
        return std::nullopt;
      }
    }
  }
  return energySum;
}

// Thread-pool variant: the framework atoms are split into blocks, one per helper thread plus the
// calling thread; an overlap in any block cancels the others. Grids are not used here.
[[nodiscard]] std::optional<RunningEnergy> computeFrameworkThreaded(const ForceField &forceField,
                                                                    const SimulationBox &simulationBox,
                                                                    std::span<const Atom> frameworkAtoms,
                                                                    double cutOffVDW, double cutOffCoulomb,
                                                                    std::span<const Atom> atoms) noexcept
{
  std::atomic_flag cancel;

  auto &pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
  const std::size_t numberOfHelperThreads = pool.getThreadCount();

  std::vector<std::future<RunningEnergy>> threads(numberOfHelperThreads);
  const std::size_t blockSize = frameworkAtoms.size() / (numberOfHelperThreads + 1);

  auto task = [cutOffVDW, cutOffCoulomb, atoms, &cancel, &forceField, &simulationBox](
                  std::span<const Atom>::iterator begin, std::span<const Atom>::iterator end) -> RunningEnergy
  {
    const double cutOffVDWSquared = cutOffVDW * cutOffVDW;
    const double cutOffChargeSquared = cutOffCoulomb * cutOffCoulomb;

    RunningEnergy energySum;
    for (std::span<const Atom>::iterator it = begin; it != end; ++it)
    {
      if (cancel.test()) return energySum;
      for (const Atom &atom : atoms)
      {
        if (!accumulateFrameworkPair(forceField, simulationBox, *it, atom, cutOffVDWSquared, cutOffChargeSquared,
                                     energySum))
        {
          cancel.test_and_set();
          return energySum;
        }
      }
    }
    return energySum;
  };

  std::span<const Atom>::iterator blockStart = frameworkAtoms.begin();
  for (std::size_t i = 0; i != numberOfHelperThreads; ++i)
  {
    std::span<const Atom>::iterator blockEnd = blockStart;
    std::advance(blockEnd, blockSize);
    threads[i] = pool.enqueue(task, blockStart, blockEnd);
    blockStart = blockEnd;
  }
  RunningEnergy energy = task(blockStart, frameworkAtoms.end());

  for (std::size_t i = 0; i != numberOfHelperThreads; ++i)
  {
    energy += threads[i].get();
  }
  if (cancel.test()) return std::nullopt;

  return energy;
}
}  // namespace

[[nodiscard]] std::optional<RunningEnergy> CBMC::computeFrameworkMoleculeEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    [[maybe_unused]] const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    double cutOffVDW, double cutOffCoulomb, std::span<const Atom> atoms) noexcept
{
  // No framework, no framework energy. Checked before the grids are touched, so an isolated-molecule
  // context ('System::makeIdealGasGrowContext') may pass an empty grid vector.
  if (frameworkAtoms.empty()) return RunningEnergy{};

  auto &pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
  switch (pool.getThreadingType())
  {
    default:
    case ThreadPool::ThreadingType::Serial:
      return computeFrameworkSerial(forceField, simulationBox, interpolationGrids, frameworkAtoms, cutOffVDW,
                                    cutOffCoulomb, atoms);
    case ThreadPool::ThreadingType::ThreadPool:
      return computeFrameworkThreaded(forceField, simulationBox, frameworkAtoms, cutOffVDW, cutOffCoulomb, atoms);
  }
}

// ---------------------------------------------------------------------------------------------------
// Inter-molecular.
// ---------------------------------------------------------------------------------------------------

[[nodiscard]] std::optional<RunningEnergy> CBMC::computeInterMolecularEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
    double cutOffVDW, double cutOffCoulomb, std::span<const Atom> atoms,
    std::optional<std::size_t> skipBackgroundMolecule) noexcept
{
  const bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffVDWSquared = cutOffVDW * cutOffVDW;
  const double cutOffChargeSquared = cutOffCoulomb * cutOffCoulomb;

  RunningEnergy energySum;
  bool overlap = false;

  for (const Atom &backgroundAtom : moleculeAtoms)
  {
    const std::size_t backgroundMolecule = static_cast<std::size_t>(backgroundAtom.moleculeId);
    if (skipBackgroundMolecule == backgroundMolecule) continue;

    for (const Atom &atom : atoms)
    {
      // Atoms of the same molecule never interact (a molecule regrown under its own id).
      if (backgroundMolecule == static_cast<std::size_t>(atom.moleculeId)) continue;

      Interactions::evaluatePair<0>(
          forceField, simulationBox, backgroundAtom, atom, cutOffVDWSquared, cutOffChargeSquared, useCharge,
          [&](const Potentials::PairDerivatives<0> &factors, const double3 &)
          {
            if (factors.energy > overlapCriteria)
            {
              overlap = true;
              return;
            }
            energySum.moleculeMoleculeVDW += factors.energy;
            energySum.addDudlambdaVDW(backgroundAtom.groupId, atom.groupId, backgroundAtom.scalingVDW, atom.scalingVDW,
                                      factors.dUdlambda);
          },
          [&](const Potentials::PairDerivatives<0> &factors, const double3 &)
          {
            if (overlap) return;
            energySum.moleculeMoleculeCharge += factors.energy;
            energySum.addDudlambdaCharge(backgroundAtom.groupId, atom.groupId, backgroundAtom.scalingCoulomb,
                                         atom.scalingCoulomb, factors.dUdlambda);
          });
      if (overlap) return std::nullopt;
    }
  }

  return energySum;
}

// ---------------------------------------------------------------------------------------------------
// The combined evaluation and the dual cut-off correction.
// ---------------------------------------------------------------------------------------------------

[[nodiscard]] std::vector<CBMC::FirstBeadTrial> CBMC::computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component, std::span<const Atom> trialPositions) noexcept
{
  std::vector<CBMC::FirstBeadTrial> energies{};
  energies.reserve(trialPositions.size());

  // Each first-bead trial is a one-atom trial set.
  for (const Atom &trialPosition : trialPositions)
  {
    std::optional<RunningEnergy> energy = computeExternalNonOverlappingEnergy(context, component, {&trialPosition, 1});
    if (!energy.has_value()) continue;
    energies.push_back({trialPosition, energy.value()});
  }
  return energies;
}

std::optional<RunningEnergy> CBMC::computeExternalNonOverlappingEnergy(const GrowContext &context,
                                                                       const Component &component,
                                                                       std::span<const Atom> trialPositionSet) noexcept
{
  if (CBMC::insideBlockedPockets(context.framework, component, trialPositionSet))
  {
    return std::nullopt;
  }

  std::optional<RunningEnergy> externalFieldEnergy =
      CBMC::computeExternalFieldEnergy(context.hasExternalField, context.forceField, context.simulationBox,
                                       context.externalFieldInterpolationGrid, context.cutOffFrameworkVDW,
                                       context.cutOffCoulomb, trialPositionSet);
  if (!externalFieldEnergy.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkEnergy = CBMC::computeFrameworkMoleculeEnergy(
      context.forceField, context.simulationBox, context.interpolationGrids, context.framework,
      context.frameworkAtoms, context.cutOffFrameworkVDW, context.cutOffCoulomb, trialPositionSet);
  if (!frameworkEnergy.has_value()) return std::nullopt;

  std::optional<RunningEnergy> interEnergy =
      CBMC::computeInterMolecularEnergy(context.forceField, context.simulationBox, context.moleculeAtoms,
                                        context.cutOffMoleculeVDW, context.cutOffCoulomb, trialPositionSet,
                                        context.skipBackgroundMolecule);
  if (!interEnergy.has_value()) return std::nullopt;

  return externalFieldEnergy.value() + interEnergy.value() + frameworkEnergy.value();
}

std::optional<RunningEnergy> CBMC::computeDualCutOffCorrection(const GrowContext &context, const Component &component,
                                                               std::span<const Atom> trialPositionSet) noexcept
{
  // The same configuration and background evaluated at the full and at the inner cut-offs.
  const GrowContext fullCutOffContext = context.withFullCutOffs();
  const GrowContext innerCutOffContext = context.withInnerCutOffs();

  std::optional<RunningEnergy> fullCutOffEnergy =
      CBMC::computeExternalNonOverlappingEnergy(fullCutOffContext, component, trialPositionSet);
  if (!fullCutOffEnergy.has_value()) return std::nullopt;

  std::optional<RunningEnergy> innerCutOffEnergy =
      CBMC::computeExternalNonOverlappingEnergy(innerCutOffContext, component, trialPositionSet);
  if (!innerCutOffEnergy.has_value()) return std::nullopt;

  return fullCutOffEnergy.value() - innerCutOffEnergy.value();
}
