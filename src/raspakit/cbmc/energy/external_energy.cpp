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
import potential_pair_vdw;
import potential_pair_coulomb;
import mc_cell_list;
import interactions_neighbourhood;
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

namespace
{
/// Whether a trial atom reads its framework energy from an interpolation grid (one exists for its type and the
/// atom is not fractional).
inline bool readsFrameworkGrid(const CBMC::GrowContext& context, const Atom& atom)
{
  return context.interpolationGrids[static_cast<std::size_t>(atom.type)].has_value() && !atom.isFractional;
}

/// The trial atoms of the sets that are evaluated by framework pairs (not from a grid), in 'pairAtoms'; false
/// when the batch needs no framework pairs at all.
bool frameworkPairAtoms(const CBMC::GrowContext& context, std::span<const std::span<const Atom>> trialSets,
                        std::vector<Atom>& pairAtoms)
{
  pairAtoms.clear();
  if (!context.frameworkCellList || context.frameworkAtoms.empty()) return false;
  for (std::span<const Atom> set : trialSets)
  {
    for (const Atom& atom : set)
    {
      if (!readsFrameworkGrid(context, atom)) pairAtoms.push_back(atom);
    }
  }
  return !pairAtoms.empty();
}

inline double frameworkGatherCutOff(const CBMC::GrowContext& context)
{
  // gathered at the full cut-offs, so the dual cut-off correction is covered as well
  return std::max(context.forceField.cutOffFrameworkVDW,
                  context.forceField.useCharge ? context.forceField.cutOffCoulomb : 0.0);
}
}  // namespace

[[nodiscard]] std::optional<RunningEnergy> CBMC::computeFrameworkMoleculeEnergy(const GrowContext& context,
                                                                                std::span<const Atom> atoms) noexcept
{
  if (context.frameworkAtoms.empty()) return RunningEnergy{};
  if (!context.frameworkCellList)
  {
    return computeFrameworkMoleculeEnergy(context.forceField, context.simulationBox, context.interpolationGrids,
                                          context.framework, context.frameworkAtoms, context.cutOffFrameworkVDW,
                                          context.cutOffCoulomb, atoms);
  }

  // the atoms with a grid are settled here, the others by pairs against the gathered framework atoms
  RunningEnergy energySum{};
  thread_local std::vector<Atom> pairAtoms{};
  pairAtoms.clear();
  for (const Atom& atom : atoms)
  {
    if (!readsFrameworkGrid(context, atom))
    {
      pairAtoms.push_back(atom);
      continue;
    }
    if (!accumulateGridEnergy(context.forceField, context.interpolationGrids, atom, energySum)) return std::nullopt;
  }
  if (pairAtoms.empty()) return energySum;

  std::optional<RunningEnergy> pairEnergy;
  if (context.frameworkNeighbourhood &&
      context.frameworkNeighbourhood->covers(pairAtoms, context.cutOffFrameworkVDW, context.cutOffCoulomb,
                                             context.forceField.useCharge))
  {
    pairEnergy = Interactions::computeNeighbourhoodEnergy(context.forceField, context.simulationBox,
                                                          *context.frameworkNeighbourhood, pairAtoms,
                                                          context.cutOffFrameworkVDW, context.cutOffCoulomb, 1.0, true);
  }
  else
  {
    // a one-off gather around this set; a single atom with no usable list is cheaper by the direct loop
    const MCCellList& list = *context.frameworkCellList;
    const bool usableList = list.enabled && list.valid && list.numberOfAtoms == context.frameworkAtoms.size();
    if (!usableList && pairAtoms.size() < 2)
    {
      pairEnergy =
          computeFrameworkSerial(context.forceField, context.simulationBox, context.interpolationGrids,
                                 context.frameworkAtoms, context.cutOffFrameworkVDW, context.cutOffCoulomb, pairAtoms);
    }
    else if (pairAtoms.size() == 1)
    {
      // a lone trial atom (a first bead): straight out of the cells, without the gather
      pairEnergy = Interactions::computeCellListAtomEnergy(
          context.forceField, context.simulationBox, list, pairAtoms[0], context.cutOffFrameworkVDW,
          context.cutOffCoulomb, 1.0, true, std::nullopt, Interactions::Neighbourhood::Background::Framework);
    }
    else
    {
      thread_local Interactions::Neighbourhood scratch{};
      const auto [centre, radius] = Interactions::Neighbourhood::boundingSphere(pairAtoms);
      scratch.gatherFramework(context.simulationBox, usableList ? &list : nullptr, context.frameworkAtoms, centre,
                              radius, frameworkGatherCutOff(context));
      pairEnergy =
          Interactions::computeNeighbourhoodEnergy(context.forceField, context.simulationBox, scratch, pairAtoms,
                                                   context.cutOffFrameworkVDW, context.cutOffCoulomb, 1.0, true);
    }
  }
  if (!pairEnergy.has_value()) return std::nullopt;
  return energySum + pairEnergy.value();
}

bool CBMC::gatherFrameworkNeighbourhood(const GrowContext& context, Interactions::Neighbourhood& neighbourhood,
                                        std::span<const std::span<const Atom>> trialSets, double extraRadius) noexcept
{
  thread_local std::vector<Atom> pairAtoms{};
  if (!frameworkPairAtoms(context, trialSets, pairAtoms)) return false;
  const auto [centre, radius] = Interactions::Neighbourhood::boundingSphere(pairAtoms);
  neighbourhood.gatherFramework(context.simulationBox, context.frameworkCellList, context.frameworkAtoms, centre,
                                radius + extraRadius, frameworkGatherCutOff(context));
  return true;
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

[[nodiscard]] std::optional<RunningEnergy> CBMC::computeInterMolecularEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox, const MCCellList &cellList,
    std::span<const Atom> moleculeAtoms, double cutOffVDW, double cutOffCoulomb, std::span<const Atom> atoms,
    std::optional<std::size_t> skipBackgroundMolecule) noexcept
{
  if (atoms.empty()) return RunningEnergy{};

  // A one-off neighbourhood of this trial set: the background within the cut-off of its bounding sphere,
  // gathered through the cells when the list is usable and by a pass over all atoms otherwise, into the
  // thread's scratch buffer and evaluated with the dense kernel. The pass over all atoms only pays off when
  // it is amortised over several trial atoms: a single atom (a first bead) is cheaper by the direct loop.
  const bool usableList = cellList.enabled && cellList.valid && cellList.numberOfAtoms == moleculeAtoms.size();
  if (!usableList && atoms.size() < 2)
  {
    return computeInterMolecularEnergy(forceField, simulationBox, moleculeAtoms, cutOffVDW, cutOffCoulomb, atoms,
                                       skipBackgroundMolecule);
  }
  const std::optional<std::uint32_t> skipped =
      skipBackgroundMolecule.has_value() ? std::optional{static_cast<std::uint32_t>(skipBackgroundMolecule.value())}
                                         : std::nullopt;
  if (atoms.size() == 1)
  {
    // a lone trial atom (a first bead): straight out of the cells, without the gather
    return Interactions::computeCellListAtomEnergy(forceField, simulationBox, cellList, atoms[0], cutOffVDW,
                                                   cutOffCoulomb, 1.0, true, skipped,
                                                   Interactions::Neighbourhood::Background::Molecules);
  }
  thread_local Interactions::Neighbourhood scratch{};
  const auto [centre, radius] = Interactions::Neighbourhood::boundingSphere(atoms);
  scratch.gather(simulationBox, usableList ? &cellList : nullptr, moleculeAtoms, centre, radius,
                 std::max(cutOffVDW, forceField.useCharge ? cutOffCoulomb : 0.0), skipped);

  return Interactions::computeNeighbourhoodEnergy(forceField, simulationBox, scratch, atoms, cutOffVDW, cutOffCoulomb,
                                                  1.0, true);
}

[[nodiscard]] std::optional<RunningEnergy> CBMC::computeInterMolecularEnergy(
    const ForceField& forceField, const SimulationBox& simulationBox, const Interactions::Neighbourhood& neighbourhood,
    double cutOffVDW, double cutOffCoulomb, std::span<const Atom> atoms) noexcept
{
  return Interactions::computeNeighbourhoodEnergy(forceField, simulationBox, neighbourhood, atoms, cutOffVDW,
                                                  cutOffCoulomb, 1.0, true);
}

void CBMC::gatherNeighbourhood(const GrowContext& context, Interactions::Neighbourhood& neighbourhood,
                               std::span<const std::span<const Atom>> trialSets, double extraRadius) noexcept
{
  const auto [centre, baseRadius] = Interactions::Neighbourhood::boundingSphere(trialSets);
  const double radius = baseRadius + extraRadius;
  const std::optional<std::uint32_t> skipped =
      context.skipBackgroundMolecule.has_value()
          ? std::optional{static_cast<std::uint32_t>(context.skipBackgroundMolecule.value())}
          : std::nullopt;
  // gathered at the full cut-offs, so the dual cut-off correction (which evaluates at the full and the
  // inner cut-offs) is covered as well
  const double cutOff = std::max(context.forceField.cutOffMoleculeVDW,
                                 context.forceField.useCharge ? context.forceField.cutOffCoulomb : 0.0);
  neighbourhood.gather(context.simulationBox, context.cellList, context.moleculeAtoms, centre, radius, cutOff, skipped);
}

// ---------------------------------------------------------------------------------------------------
// The combined evaluation and the dual cut-off correction.
// ---------------------------------------------------------------------------------------------------

[[nodiscard]] std::vector<CBMC::FirstBeadTrial> CBMC::computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component, std::span<const Atom> trialPositions) noexcept
{
  std::vector<CBMC::FirstBeadTrial> energies{};
  energies.reserve(trialPositions.size());

  // Each first-bead trial is a one-atom trial set. The trials are spread over the whole box, so they share
  // no local neighbourhood; when the system's cell list (molecule or framework) is not usable (a box too
  // small for one) they instead share one unshifted copy of all atoms, over which the dense kernels apply the
  // minimum image per pair -- cheaper than the atom-by-atom loop, provided the copy is amortised over several
  // trials. A context without a cell list keeps the brute-force loop.
  GrowContext shared(context);
  if (trialPositions.size() >= 2)
  {
    const bool usableList = context.cellList && context.cellList->enabled && context.cellList->valid &&
                            context.cellList->numberOfAtoms == context.moleculeAtoms.size();
    if (context.cellList && !usableList && !context.neighbourhood && !context.moleculeAtoms.empty())
    {
      thread_local Interactions::Neighbourhood all{};
      const std::optional<std::uint32_t> skipped =
          context.skipBackgroundMolecule.has_value()
              ? std::optional{static_cast<std::uint32_t>(context.skipBackgroundMolecule.value())}
              : std::nullopt;
      const double cutOff = std::max(context.forceField.cutOffMoleculeVDW,
                                     context.forceField.useCharge ? context.forceField.cutOffCoulomb : 0.0);
      all.gatherAll(context.moleculeAtoms, cutOff, skipped);
      shared.neighbourhood = &all;
    }

    const bool usableFrameworkList = context.frameworkCellList && context.frameworkCellList->enabled &&
                                     context.frameworkCellList->valid &&
                                     context.frameworkCellList->numberOfAtoms == context.frameworkAtoms.size();
    if (context.frameworkCellList && !usableFrameworkList && !context.frameworkNeighbourhood)
    {
      thread_local std::vector<Atom> pairAtoms{};
      if (frameworkPairAtoms(context, std::span<const std::span<const Atom>>(&trialPositions, 1), pairAtoms))
      {
        thread_local Interactions::Neighbourhood allFramework{};
        allFramework.gatherAllFramework(context.frameworkAtoms, frameworkGatherCutOff(context));
        shared.frameworkNeighbourhood = &allFramework;
      }
    }
  }

  for (const Atom &trialPosition : trialPositions)
  {
    std::optional<RunningEnergy> energy = computeExternalNonOverlappingEnergy(shared, component, {&trialPosition, 1});
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

  std::optional<RunningEnergy> frameworkEnergy = CBMC::computeFrameworkMoleculeEnergy(context, trialPositionSet);
  if (!frameworkEnergy.has_value()) return std::nullopt;

  std::optional<RunningEnergy> interEnergy;
  if (context.neighbourhood && context.neighbourhood->covers(trialPositionSet, context.cutOffMoleculeVDW,
                                                             context.cutOffCoulomb, context.forceField.useCharge))
  {
    interEnergy = CBMC::computeInterMolecularEnergy(context.forceField, context.simulationBox, *context.neighbourhood,
                                                    context.cutOffMoleculeVDW, context.cutOffCoulomb, trialPositionSet);
  }
  else if (context.cellList)
  {
    interEnergy = CBMC::computeInterMolecularEnergy(
        context.forceField, context.simulationBox, *context.cellList, context.moleculeAtoms, context.cutOffMoleculeVDW,
        context.cutOffCoulomb, trialPositionSet, context.skipBackgroundMolecule);
  }
  else
  {
    interEnergy = CBMC::computeInterMolecularEnergy(context.forceField, context.simulationBox, context.moleculeAtoms,
                                                    context.cutOffMoleculeVDW, context.cutOffCoulomb, trialPositionSet,
                                                    context.skipBackgroundMolecule);
  }
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
