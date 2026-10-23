module;

module interactions_framework_molecule;

import std;

import double3;
import double4;
import double3x3;
import double3x3x3;
import energy_status;
import potential_pair_derivatives;
import potential_pair_vdw;
import potential_pair_coulomb;
import potential_tricubic_derivative_lj;
import potential_tricubic_derivative_real_ewald;
import potential_triquintic_derivative_lj;
import potential_electrostatics;
import interactions_pair_kernel;
import simulationbox;
import framework;
import forcefield;
import atom;
import energy_dudlambda;
import energy_status_inter;
import running_energy;
import units;
import threadpool;
// import threading;
import tricubic_derivative_factor;
import triquintic_derivative_factor;
import framework;
import component;
import interpolation_energy_grid;
import mc_cell_list;
import interactions_neighbourhood;

namespace
{
// The grids are far larger than the caches and the atoms of a loop are spread over the whole box, so an
// interpolation waits for memory. Loops over many atoms fetch the grid cells this many atoms ahead.
constexpr std::size_t gridPrefetchDistance = 8;

inline void prefetchGridCell(const ForceField& forceField,
                             const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
                             std::span<const Atom> atoms, std::size_t index)
{
  if (index >= atoms.size()) return;
  const Atom& atom = atoms[index];
  const std::optional<InterpolationEnergyGrid>& grid = interpolationGrids[static_cast<std::size_t>(atom.type)];
  if (!grid.has_value() || atom.isFractional) return;
  grid->prefetch(atom.position);
  if (forceField.useCharge) interpolationGrids.back()->prefetch(atom.position);
}

inline void prefetchGridCells(const ForceField& forceField,
                              const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
                              std::span<const Atom> atoms, std::size_t begin, std::size_t end)
{
  for (std::size_t index = begin; index < end; ++index)
  {
    prefetchGridCell(forceField, interpolationGrids, atoms, index);
  }
}

// The grid energies of one atom (the Van der Waals grid of its type, and the Coulomb grid when charges are on).
inline std::pair<double, double> gridEnergies(const ForceField& forceField,
                                              const std::vector<std::optional<InterpolationEnergyGrid>>& grids,
                                              const Atom& atom)
{
  const InterpolationEnergyGrid& grid = grids[static_cast<std::size_t>(atom.type)].value();
  if (forceField.useCharge)
  {
    const auto [energyVDW, energyCoulomb] = grid.interpolate(grids.back().value(), atom.position);
    return {energyVDW, atom.charge * energyCoulomb};
  }
  return {grid.interpolate(atom.position), 0.0};
}
}  // namespace

RunningEnergy Interactions::computeFrameworkMoleculeEnergy(
    const ForceField& forceField, const SimulationBox& simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energySum{};

  if (!framework.has_value()) return energySum;
  if (moleculeAtoms.empty()) return energySum;

  prefetchGridCells(forceField, interpolationGrids, moleculeAtoms, 0, gridPrefetchDistance);
  for (std::size_t index = 0; index < moleculeAtoms.size(); ++index)
  {
    const Atom& moleculeAtom = moleculeAtoms[index];
    std::size_t typeB = static_cast<std::size_t>(moleculeAtom.type);
    bool isFractional = static_cast<bool>(moleculeAtom.isFractional);

    if (interpolationGrids[typeB].has_value() && !isFractional &&
        forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
    {
      prefetchGridCell(forceField, interpolationGrids, moleculeAtoms, index + gridPrefetchDistance);
      const auto [energyVDW, energyCoulomb] = gridEnergies(forceField, interpolationGrids, moleculeAtom);
      energySum.frameworkMoleculeVDW += energyVDW;
      energySum.frameworkMoleculeCharge += energyCoulomb;
    }
    else
    {
      forEachFrameworkMoleculePair<0>(
          forceField, simulationBox, moleculeAtom, frameworkAtoms,
          [&](std::size_t, const Atom& frameworkAtom, const Potentials::PairDerivatives<0>& factors, const double3&)
          {
            energySum.frameworkMoleculeVDW += factors.energy;
            energySum.addDudlambdaVDW(moleculeAtom.groupId, frameworkAtom.groupId, moleculeAtom.scalingVDW,
                                      frameworkAtom.scalingVDW, factors.dUdlambda);
          },
          [&](std::size_t, const Atom& frameworkAtom, const Potentials::PairDerivatives<0>& factors, const double3&)
          {
            energySum.frameworkMoleculeCharge += factors.energy;
            energySum.addDudlambdaCharge(moleculeAtom.groupId, frameworkAtom.groupId, moleculeAtom.scalingCoulomb,
                                         frameworkAtom.scalingCoulomb, factors.dUdlambda);
          });
    }
  }
  return energySum;
}

namespace
{
// The framework-molecule tail correction of a molecule atom B is sum_A s_A s_B 2 (2 pi / V) C_AB over all framework
// atoms A, which depends on the framework only through the scalingVDW-weighted count of its atoms per type (and
// per dU/dlambda group, for the rare framework atom that carries one). The framework is therefore reduced once to
// per-type weights w_B = 2 (2 pi / V) sum_A N_A C_AB, and every molecule atom costs a lookup instead of a pass over
// all framework atoms.
struct FrameworkTailWeights
{
  std::vector<double> perType;                                                   // w_B
  std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> perGroupType;  // w_B restricted to group g

  FrameworkTailWeights(const ForceField& forceField, const SimulationBox& simulationBox,
                       std::span<const Atom> frameworkAtoms)
  {
    const std::size_t numberOfPseudoAtoms = forceField.numberOfPseudoAtoms;
    std::vector<double> counts(numberOfPseudoAtoms, 0.0);
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> groupCounts{};
    for (const Atom& atom : frameworkAtoms)
    {
      counts[static_cast<std::size_t>(atom.type)] += atom.scalingVDW;
      if (atom.groupId != 0)
      {
        std::vector<double>& group = groupCounts[static_cast<std::size_t>(atom.groupId - 1)];
        if (group.empty()) group.assign(numberOfPseudoAtoms, 0.0);
        group[static_cast<std::size_t>(atom.type)] += atom.scalingVDW;
      }
    }

    const double factor = 2.0 * (2.0 * std::numbers::pi / simulationBox.volume);
    auto weights = [&](const std::vector<double>& typeCounts)
    {
      std::vector<double> w(numberOfPseudoAtoms, 0.0);
      for (std::size_t typeB = 0; typeB < numberOfPseudoAtoms; ++typeB)
      {
        double sum = 0.0;
        for (std::size_t typeA = 0; typeA < numberOfPseudoAtoms; ++typeA)
        {
          sum += typeCounts[typeA] * forceField(typeA, typeB).tailCorrectionEnergy;
        }
        w[typeB] = factor * sum;
      }
      return w;
    };
    perType = weights(counts);
    for (std::size_t group = 0; group < maximumNumberOfDUDlambdaGroups; ++group)
    {
      if (!groupCounts[group].empty()) perGroupType[group] = weights(groupCounts[group]);
    }
  }

  // Adds sign times the tail correction of 'atom' with the whole framework (energy and dU/dlambda).
  void accumulate(RunningEnergy& energySum, const Atom& atom, double sign) const
  {
    const std::size_t typeB = static_cast<std::size_t>(atom.type);
    const double scalingVDWB = atom.scalingVDW;
    const double w = perType[typeB];

    energySum.tail += sign * scalingVDWB * w;
    // the molecule atom's group collects sum_A s_A * temp_AB = w_B; a framework group collects s_B times its share
    if (atom.groupId != 0) energySum.dudlambdaVDW[static_cast<std::size_t>(atom.groupId - 1)] += sign * w;
    for (std::size_t group = 0; group < maximumNumberOfDUDlambdaGroups; ++group)
    {
      if (!perGroupType[group].empty())
      {
        energySum.dudlambdaVDW[group] += sign * scalingVDWB * perGroupType[group][typeB];
      }
    }
  }
};

// Whether any framework-molecule pair involving one of 'atoms' carries a tail correction at all; false lets the
// callers skip the framework entirely (the common case: tail corrections off, or only between adsorbates).
bool hasFrameworkTailCorrection(const ForceField& forceField, std::span<const Atom> atoms)
{
  const std::size_t numberOfPseudoAtoms = forceField.numberOfPseudoAtoms;
  for (const Atom& atom : atoms)
  {
    const std::size_t typeB = static_cast<std::size_t>(atom.type);
    for (std::size_t typeA = 0; typeA < numberOfPseudoAtoms; ++typeA)
    {
      if (forceField(typeA, typeB).tailCorrectionEnergy != 0.0) return true;
    }
  }
  return false;
}
}  // namespace

RunningEnergy Interactions::computeFrameworkMoleculeTailEnergy(const ForceField& forceField,
                                                               const SimulationBox& simulationBox,
                                                               std::span<const Atom> frameworkAtoms,
                                                               std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energySum{};
  if (frameworkAtoms.empty() || moleculeAtoms.empty()) return energySum;
  if (!hasFrameworkTailCorrection(forceField, moleculeAtoms)) return energySum;

  const FrameworkTailWeights weights(forceField, simulationBox, frameworkAtoms);
  for (const Atom& atom : moleculeAtoms)
  {
    weights.accumulate(energySum, atom, 1.0);
  }
  return energySum;
}

RunningEnergy Interactions::computeFrameworkMoleculeTailEnergyAggregated(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const double> effectiveTypeCounts,
    const std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& groupCounts) noexcept
{
  RunningEnergy energySum{};

  if (frameworkAtoms.empty()) return energySum;

  const std::size_t numberOfPseudoAtoms = forceField.numberOfPseudoAtoms;
  const double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  // Effective per-type counts of the framework (framework atoms carry no dU/dlambda group).
  std::vector<double> frameworkTypeCounts(numberOfPseudoAtoms, 0.0);
  for (const Atom& atom : frameworkAtoms)
  {
    frameworkTypeCounts[static_cast<std::size_t>(atom.type)] += atom.scalingVDW;
  }

  for (std::size_t typeB = 0; typeB < numberOfPseudoAtoms; ++typeB)
  {
    // w_b = 2 * (2 pi / V) * sum_a N_framework[a] * C_ab
    double wB = 0.0;
    for (std::size_t typeA = 0; typeA < numberOfPseudoAtoms; ++typeA)
    {
      wB += frameworkTypeCounts[typeA] * forceField(typeA, typeB).tailCorrectionEnergy;
    }
    wB *= 2.0 * preFactor;

    energySum.tail += effectiveTypeCounts[typeB] * wB;

    for (std::size_t group = 0; group < maximumNumberOfDUDlambdaGroups; ++group)
    {
      const double nBg = groupCounts[group][typeB];
      if (nBg != 0.0)
      {
        energySum.dudlambdaVDW[group] += nBg * wB;
      }
    }
  }

  return energySum;
}

// Used in Translation and Rotation
//

[[nodiscard]] std::optional<RunningEnergy> Interactions::computeFrameworkMoleculeEnergyDifference(
    const ForceField& forceField, const SimulationBox& simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  if (!framework.has_value()) return energySum;

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (auto& atom : newatoms)
  {
    double3 posB = atom.position;
    std::size_t typeB = static_cast<std::size_t>(atom.type);
    std::uint8_t groupIdB = atom.groupId;
    bool isFractional = static_cast<bool>(atom.isFractional);
    double scalingVDWB = atom.scalingVDW;
    double scalingCoulombB = atom.scalingCoulomb;
    double chargeB = atom.charge;

    if (interpolationGrids[typeB].has_value() && !isFractional &&
        forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
    {
      const auto [energyVDW, energyCoulomb] = gridEnergies(forceField, interpolationGrids, atom);
      if (energyVDW > overlapCriteria)
      {
        return std::nullopt;
      }
      energySum.frameworkMoleculeVDW += energyVDW;
      energySum.frameworkMoleculeCharge += energyCoulomb;
    }
    else
    {
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        double3 posA = it1->position;
        std::size_t typeA = static_cast<std::size_t>(it1->type);
        std::uint8_t groupIdA = it1->groupId;
        double scalingVDWA = it1->scalingVDW;
        double scalingCoulombA = it1->scalingCoulomb;
        double chargeA = it1->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);
          if (energyFactor.energy > overlapCriteria) return std::nullopt;

          energySum.frameworkMoleculeVDW += energyFactor.energy;
          energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, energyFactor.dUdlambda);
        }
        if (useCharge && rr < cutOffChargeSquared && chargeA * chargeB != 0.0 &&
            scalingCoulombA * scalingCoulombB != 0.0)
        {
          double r = std::sqrt(rr);
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.frameworkMoleculeCharge += energyFactor.energy;
          energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, energyFactor.dUdlambda);
        }
      }
    }
  }

  for (auto& atom : oldatoms)
  {
    double3 posB = atom.position;
    std::size_t typeB = static_cast<std::size_t>(atom.type);
    std::uint8_t groupIdB = atom.groupId;
    bool isFractional = static_cast<bool>(atom.isFractional);
    double scalingVDWB = atom.scalingVDW;
    double scalingCoulombB = atom.scalingCoulomb;
    double chargeB = atom.charge;

    if (interpolationGrids[typeB].has_value() && !isFractional &&
        forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
    {
      const auto [energyVDW, energyCoulomb] = gridEnergies(forceField, interpolationGrids, atom);
      energySum.frameworkMoleculeVDW -= energyVDW;
      energySum.frameworkMoleculeCharge -= energyCoulomb;
    }
    else
    {
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        double3 posA = it1->position;
        std::size_t typeA = static_cast<std::size_t>(it1->type);
        std::uint8_t groupIdA = it1->groupId;
        double scalingVDWA = it1->scalingVDW;
        double scalingCoulombA = it1->scalingCoulomb;
        double chargeA = it1->charge;

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);

        if (rr < cutOffFrameworkVDWSquared)
        {
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);

          energySum.frameworkMoleculeVDW -= energyFactor.energy;
          energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, -energyFactor.dUdlambda);
        }
        if (useCharge && rr < cutOffChargeSquared && chargeA * chargeB != 0.0 &&
            scalingCoulombA * scalingCoulombB != 0.0)
        {
          double r = std::sqrt(rr);
          Potentials::PairDerivatives<0> energyFactor =
              Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

          energySum.frameworkMoleculeCharge -= energyFactor.energy;
          energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, -energyFactor.dUdlambda);
        }
      }
    }
  }

  return std::optional{energySum};
}

[[nodiscard]] std::optional<RunningEnergy> Interactions::computeFrameworkMoleculeEnergyDifference(
    const ForceField& forceField, const SimulationBox& simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<Framework>& framework, const MCCellList& frameworkCellList,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  RunningEnergy energySum{};
  if (!framework.has_value()) return energySum;

  const bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.energyOverlapCriteria;

  // The atoms read from a grid are settled here; the others go through the gathered neighbourhood (same grid
  // rule as the brute-force loop above).
  auto readsGrid = [&](const Atom& atom)
  {
    return interpolationGrids[static_cast<std::size_t>(atom.type)].has_value() && !atom.isFractional &&
           forceField.chargeMethod == ForceField::ChargeMethod::Ewald;
  };
  thread_local std::vector<Atom> pairNew{};
  thread_local std::vector<Atom> pairOld{};
  pairNew.clear();
  pairOld.clear();
  // the grid cells of both sets are fetched up front: a move touches a handful of atoms, each a memory wait
  prefetchGridCells(forceField, interpolationGrids, newatoms, 0, newatoms.size());
  prefetchGridCells(forceField, interpolationGrids, oldatoms, 0, oldatoms.size());
  for (const Atom& atom : newatoms)
  {
    if (!readsGrid(atom))
    {
      pairNew.push_back(atom);
      continue;
    }
    const auto [energyVDW, energyCoulomb] = gridEnergies(forceField, interpolationGrids, atom);
    if (energyVDW > overlapCriteria) return std::nullopt;
    energySum.frameworkMoleculeVDW += energyVDW;
    energySum.frameworkMoleculeCharge += energyCoulomb;
  }
  for (const Atom& atom : oldatoms)
  {
    if (!readsGrid(atom))
    {
      pairOld.push_back(atom);
      continue;
    }
    const auto [energyVDW, energyCoulomb] = gridEnergies(forceField, interpolationGrids, atom);
    energySum.frameworkMoleculeVDW -= energyVDW;
    energySum.frameworkMoleculeCharge -= energyCoulomb;
  }
  if ((pairNew.empty() && pairOld.empty()) || frameworkAtoms.empty()) return energySum;

  // The new and the old positions lie within the move's displacement of each other: the framework atoms within
  // the cut-off of their common bounding sphere are gathered once and both sets are evaluated against them.
  const bool usableList =
      frameworkCellList.enabled && frameworkCellList.valid && frameworkCellList.numberOfAtoms == frameworkAtoms.size();
  thread_local Neighbourhood scratch{};
  const std::array<std::span<const Atom>, 2> sets{std::span<const Atom>(pairNew), std::span<const Atom>(pairOld)};
  const auto [centre, radius] = Neighbourhood::boundingSphere(sets);
  const double cutOff = std::max(forceField.cutOffFrameworkVDW, useCharge ? forceField.cutOffCoulomb : 0.0);
  scratch.gatherFramework(simulationBox, usableList ? &frameworkCellList : nullptr, frameworkAtoms, centre, radius,
                          cutOff);

  std::optional<RunningEnergy> energyNew = computeNeighbourhoodEnergy(
      forceField, simulationBox, scratch, pairNew, forceField.cutOffFrameworkVDW, forceField.cutOffCoulomb, 1.0, true);
  if (!energyNew.has_value()) return std::nullopt;
  std::optional<RunningEnergy> energyOld =
      computeNeighbourhoodEnergy(forceField, simulationBox, scratch, pairOld, forceField.cutOffFrameworkVDW,
                                 forceField.cutOffCoulomb, -1.0, false);

  return std::optional{energySum + energyNew.value() + energyOld.value()};
}

std::optional<RunningEnergy> Interactions::computeFrameworkMoleculeEnergyDifference(
    const ForceField& forceField, const SimulationBox& simulationBox,
    [[maybe_unused]] const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    [[maybe_unused]] const std::optional<Framework>& framework, std::span<const Atom> frameworkAtoms,
    std::span<double3> electricFieldMoleculeNew, std::span<double3> electricFieldMoleculeOld,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posA = it1->position;
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = newatoms.begin(); it2 != newatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(newatoms.begin(), it2));
      double3 posB = it2->position;
      std::size_t typeB = static_cast<std::size_t>(it2->type);
      std::uint8_t groupIdB = it2->groupId;
      double scalingVDWB = it2->scalingVDW;
      double scalingCoulombB = it2->scalingCoulomb;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);
      if (rr < cutOffFrameworkVDWSquared)
      {
        Potentials::PairDerivatives<0> energyFactor =
            Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);
        if (energyFactor.energy > overlapCriteria) return std::nullopt;

        energySum.frameworkMoleculeVDW += energyFactor.energy;
        energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, energyFactor.dUdlambda);
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        Potentials::PairDerivatives<0> energyFactor =
            Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, energyFactor.dUdlambda);

        Potentials::PairDerivatives<1> gradientFactor =
            Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);
        electricFieldMoleculeNew[indexB] += scalingCoulombA * chargeA * gradientFactor.firstDerivativeFactor * dr;
      }
    }

    for (std::span<const Atom>::iterator it2 = oldatoms.begin(); it2 != oldatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(oldatoms.begin(), it2));
      double3 posB = it2->position;
      std::size_t typeB = static_cast<std::size_t>(it2->type);
      std::uint8_t groupIdB = it2->groupId;
      double scalingVDWB = it2->scalingVDW;
      double scalingCoulombB = it2->scalingCoulomb;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffFrameworkVDWSquared)
      {
        Potentials::PairDerivatives<0> energyFactor =
            Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW -= energyFactor.energy;
        energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, -energyFactor.dUdlambda);
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        Potentials::PairDerivatives<0> energyFactor =
            Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge -= energyFactor.energy;
        energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, -energyFactor.dUdlambda);

        Potentials::PairDerivatives<1> gradientFactor =
            Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);
        electricFieldMoleculeOld[indexB] -= scalingCoulombA * chargeA * gradientFactor.firstDerivativeFactor * dr;
      }
    }
  }

  return std::optional{energySum};
}

void Interactions::computeFrameworkMoleculeElectricFieldDifference(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<double3> electricFieldMoleculeNew, std::span<double3> electricFieldMoleculeOld,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posA = it1->position;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = newatoms.begin(); it2 != newatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(newatoms.begin(), it2));
      double3 posB = it2->position;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);

        Potentials::PairDerivatives<1> gradientFactor =
            Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);
        electricFieldMoleculeNew[indexB] += scalingCoulombA * chargeA * gradientFactor.firstDerivativeFactor * dr;
      }
    }

    for (std::span<const Atom>::iterator it2 = oldatoms.begin(); it2 != oldatoms.end(); ++it2)
    {
      std::size_t indexB = static_cast<std::size_t>(std::distance(oldatoms.begin(), it2));
      double3 posB = it2->position;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);

        Potentials::PairDerivatives<1> gradientFactor =
            Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);
        electricFieldMoleculeOld[indexB] -= scalingCoulombA * chargeA * gradientFactor.firstDerivativeFactor * dr;
      }
    }
  }
}

[[nodiscard]] RunningEnergy Interactions::computeFrameworkMoleculeTailEnergyDifference(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept
{
  RunningEnergy energySum{};
  if (frameworkAtoms.empty()) return energySum;
  if (!hasFrameworkTailCorrection(forceField, newatoms) && !hasFrameworkTailCorrection(forceField, oldatoms))
  {
    return energySum;
  }

  const FrameworkTailWeights weights(forceField, simulationBox, frameworkAtoms);
  for (const Atom& atom : newatoms)
  {
    weights.accumulate(energySum, atom, 1.0);
  }
  for (const Atom& atom : oldatoms)
  {
    weights.accumulate(energySum, atom, -1.0);
  }
  return energySum;
}

RunningEnergy Interactions::computeFrameworkMoleculeGradient(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, std::span<AtomDynamics> moleculeDynamics,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::optional<Framework>& framework, std::span<AtomDynamics> frameworkDynamics) noexcept
{
  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  // Any mobile host atom disables energy grids and receives Newton's-third-law gradients.
  const bool mobileFramework =
      framework && framework->hasMobileAtoms() && frameworkDynamics.size() == frameworkAtoms.size();

  if (moleculeAtoms.empty()) return energySum;

  if (!mobileFramework) prefetchGridCells(forceField, interpolationGrids, moleculeAtoms, 0, gridPrefetchDistance);
  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    std::size_t indexA = static_cast<std::size_t>(it1 - moleculeAtoms.begin());
    const Atom& moleculeAtom = *it1;
    std::size_t typeA = static_cast<std::size_t>(moleculeAtom.type);
    bool isFractional = static_cast<bool>(moleculeAtom.isFractional);

    if (!mobileFramework && interpolationGrids[typeA].has_value() && !isFractional &&
        forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
    {
      prefetchGridCell(forceField, interpolationGrids, moleculeAtoms, indexA + gridPrefetchDistance);
      if (useCharge)
      {
        const auto [vdw, realEwald] =
            interpolationGrids[typeA]->interpolateGradient(interpolationGrids.back().value(), moleculeAtom.position);
        energySum.frameworkMoleculeVDW += vdw.first;
        energySum.frameworkMoleculeCharge += moleculeAtom.charge * realEwald.first;
        moleculeDynamics[indexA].gradient += vdw.second + moleculeAtom.charge * realEwald.second;
      }
      else
      {
        auto [energy_vdw, gradient_vdw] = interpolationGrids[typeA]->interpolateGradient(moleculeAtom.position);
        energySum.frameworkMoleculeVDW += energy_vdw;
        moleculeDynamics[indexA].gradient += gradient_vdw;
      }
    }
    else
    {
      forEachFrameworkMoleculePair<1>(
          forceField, simulationBox, moleculeAtom, frameworkAtoms,
          [&](std::size_t indexB, const Atom& frameworkAtom, const Potentials::PairDerivatives<1>& factors,
              const double3& dr)
          {
            energySum.frameworkMoleculeVDW += factors.energy;
            energySum.addDudlambdaVDW(moleculeAtom.groupId, frameworkAtom.groupId, moleculeAtom.scalingVDW,
                                      frameworkAtom.scalingVDW, factors.dUdlambda);

            const double3 f = factors.firstDerivativeFactor * dr;

            moleculeDynamics[indexA].gradient += f;
            if (mobileFramework && !framework->isFixedAtom(indexB)) frameworkDynamics[indexB].gradient -= f;
          },
          [&](std::size_t indexB, const Atom& frameworkAtom, const Potentials::PairDerivatives<1>& factors,
              const double3& dr)
          {
            energySum.frameworkMoleculeCharge += factors.energy;
            energySum.addDudlambdaCharge(moleculeAtom.groupId, frameworkAtom.groupId, moleculeAtom.scalingCoulomb,
                                         frameworkAtom.scalingCoulomb, factors.dUdlambda);

            const double3 g = factors.firstDerivativeFactor * dr;

            moleculeDynamics[indexA].gradient += g;
            if (mobileFramework && !framework->isFixedAtom(indexB)) frameworkDynamics[indexB].gradient -= g;
          });
    }
  }
  return energySum;
}

[[nodiscard]] std::pair<EnergyStatus, double3x3> Interactions::computeFrameworkMoleculeEnergyStrainDerivative(
    const ForceField& forceField, const std::optional<Framework>& framework,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::vector<Component>& components, const SimulationBox& simulationBox, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, std::span<AtomDynamics> moleculeDynamics,
    const PolarizationFieldStrain* polarizationGather) noexcept
{
  double3x3 strainDerivativeTensor;
  EnergyStatus energy(1, framework.has_value() ? 1 : 0, components.size());

  bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;

  if (moleculeAtoms.empty()) return std::make_pair(energy, strainDerivativeTensor);

  // Framework atoms are fixed, so their gradient is never consumed and is not accumulated.
  prefetchGridCells(forceField, interpolationGrids, moleculeAtoms, 0, gridPrefetchDistance);
  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    std::size_t indexA = static_cast<std::size_t>(it1 - moleculeAtoms.begin());
    const double3 posA = it1->position;
    std::size_t compA = static_cast<std::size_t>(it1->componentId);
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    bool isFractional = static_cast<bool>(it1->isFractional);
    double scalingVDWA = it1->scalingVDW;
    double chargeA = it1->charge;

    const bool gatherFieldForAtom =
        polarizationGather != nullptr && useCharge && polarizationGather->polarizability[indexA] != 0.0;
    const double3 sigmaA =
        gatherFieldForAtom ? polarizationGather->centerOfMassOffset[indexA] : double3(0.0, 0.0, 0.0);

    if (interpolationGrids[typeA].has_value() && !isFractional &&
        forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
    {
      prefetchGridCell(forceField, interpolationGrids, moleculeAtoms, indexA + gridPrefetchDistance);
      if (useCharge)
      {
        const auto [vdw, realEwald] =
            interpolationGrids[typeA]->interpolateGradient(interpolationGrids.back().value(), posA);
        energy.frameworkComponentEnergy(0, compA).VanDerWaals += EnergyDuDlambda(vdw.first, 0.0);
        energy.frameworkComponentEnergy(0, compA).CoulombicReal += EnergyDuDlambda(chargeA * realEwald.first, 0.0);
        const double3 f = vdw.second + chargeA * realEwald.second;

        moleculeDynamics[indexA].gradient += f;
        accumulateStrainDerivative(strainDerivativeTensor, f, posA);
      }
      else
      {
        auto [energy_vdw, gradient_vdw] = interpolationGrids[typeA]->interpolateGradient(posA);
        energy.frameworkComponentEnergy(0, compA).VanDerWaals += EnergyDuDlambda(energy_vdw, 0.0);
        const double3 f = gradient_vdw;

        moleculeDynamics[indexA].gradient += f;
        accumulateStrainDerivative(strainDerivativeTensor, f, posA);
      }

      if (gatherFieldForAtom)
      {
        // The electric field is never grid-interpolated: walk the framework charges explicitly for the
        // polarization field (matching computeTotalElectricField), even though the energy came from the grid.
        for (const Atom& frameworkAtom : frameworkAtoms)
        {
          double3 dr = posA - frameworkAtom.position;
          dr = simulationBox.applyPeriodicBoundaryConditions(dr);
          const double rr = double3::dot(dr, dr);
          if (rr >= cutOffChargeSquared) continue;
          const double r = std::sqrt(rr);
          const Potentials::PairDerivatives<2> unitFactors =
              Potentials::potentialCoulomb<2>(forceField, 1.0, 1.0, r, 1.0, 1.0);
          accumulatePolarizationFieldStrain(*polarizationGather, indexA,
                                            frameworkAtom.scalingCoulomb * frameworkAtom.charge, dr, dr - sigmaA,
                                            unitFactors.firstDerivativeFactor, unitFactors.secondDerivativeFactor);
        }
      }
    }
    else
    {
      for (std::span<const Atom>::iterator it2 = frameworkAtoms.begin(); it2 != frameworkAtoms.end(); ++it2)
      {
        std::size_t typeB = static_cast<std::size_t>(it2->type);
        double scalingVDWB = it2->scalingVDW;

        EnergyDuDlambda temp(
            preFactor * scalingVDWA * scalingVDWB * forceField(typeB, typeA).tailCorrectionEnergy, 0.0);
        energy.frameworkComponentEnergy(0, compA).VanDerWaalsTailCorrection += 2.0 * temp;

        const auto accumulateGradientAndStrain = [&](const double3& g, const double3& dr)
        {
          moleculeDynamics[indexA].gradient += g;
          accumulateStrainDerivative(strainDerivativeTensor, g, dr);
        };

        if (!gatherFieldForAtom)
        {
          evaluatePair<1>(
              forceField, simulationBox, *it1, *it2, cutOffFrameworkVDWSquared, cutOffChargeSquared, useCharge,
              [&](const Potentials::PairDerivatives<1>& factors, const double3& dr)
              {
                energy.frameworkComponentEnergy(0, compA).VanDerWaals += EnergyDuDlambda(factors.energy, 0.0);
                accumulateGradientAndStrain(factors.firstDerivativeFactor * dr, dr);
              },
              [&](const Potentials::PairDerivatives<1>& factors, const double3& dr)
              {
                energy.frameworkComponentEnergy(0, compA).CoulombicReal += EnergyDuDlambda(factors.energy, 0.0);
                accumulateGradientAndStrain(factors.firstDerivativeFactor * dr, dr);
              });
        }
        else
        {
          // Fused polarization path: the same pair walk yields both the pair energy/virial and the
          // polarization field with its strain response (scaled by the framework source charge only). The
          // pair term is taken at order 1 so it is identical to the non-polarization branch above (the
          // order-1 Ewald term is tabulated, the order-2 field factors are closed-form).
          double3 dr = posA - it2->position;
          dr = simulationBox.applyPeriodicBoundaryConditions(dr);
          const double rr = double3::dot(dr, dr);

          if (rr < cutOffFrameworkVDWSquared)
          {
            const Potentials::PairDerivatives<1> factors =
                Potentials::potentialVDW<1>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);
            energy.frameworkComponentEnergy(0, compA).VanDerWaals += EnergyDuDlambda(factors.energy, 0.0);
            accumulateGradientAndStrain(factors.firstDerivativeFactor * dr, dr);
          }
          if (rr < cutOffChargeSquared)
          {
            const double r = std::sqrt(rr);
            const Potentials::PairDerivatives<1> pairFactors = Potentials::potentialCoulomb<1>(
                forceField, it1->scalingCoulomb, it2->scalingCoulomb, r, chargeA, it2->charge);
            energy.frameworkComponentEnergy(0, compA).CoulombicReal += EnergyDuDlambda(pairFactors.energy, 0.0);
            accumulateGradientAndStrain(pairFactors.firstDerivativeFactor * dr, dr);

            const Potentials::PairDerivatives<2> unitFactors =
                Potentials::potentialCoulomb<2>(forceField, 1.0, 1.0, r, 1.0, 1.0);
            const double scaledChargeB = it2->scalingCoulomb * it2->charge;
            accumulatePolarizationFieldStrain(*polarizationGather, indexA, scaledChargeB, dr, dr - sigmaA,
                                              unitFactors.firstDerivativeFactor, unitFactors.secondDerivativeFactor);
          }
        }
      }
    }
  }

  return std::make_pair(energy, strainDerivativeTensor);
}

void Interactions::computeFrameworkMoleculeElectrostaticPotential(const ForceField& forceField,
                                                                  const SimulationBox& simulationBox,
                                                                  std::span<double> electricPotentialMolecules,
                                                                  std::span<const Atom> frameworkAtoms,
                                                                  std::span<const Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  bool useCharge = forceField.useCharge;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (!useCharge) return;
  if (moleculeAtoms.empty()) return;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      posB = it2->position;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);

        double potential = Potentials::potentialElectrostatics(forceField, 1.0, r, 1.0);

        std::size_t indexB = static_cast<std::size_t>(std::distance(moleculeAtoms.begin(), it2));
        electricPotentialMolecules[indexB] += scalingCoulombA * chargeA * potential;
      }
    }
  }
}

RunningEnergy Interactions::computeFrameworkMoleculeElectricField(const ForceField& forceField,
                                                                  const SimulationBox& simulationBox,
                                                                  std::span<double3> electricFieldMolecules,
                                                                  std::span<const Atom> frameworkAtoms,
                                                                  std::span<const Atom> moleculeAtoms) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  RunningEnergy energySum{};

  bool useCharge = forceField.useCharge;
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (!useCharge) return energySum;
  if (moleculeAtoms.empty()) return energySum;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      posB = it2->position;
      std::size_t typeB = static_cast<std::size_t>(it2->type);
      std::uint8_t groupIdB = it2->groupId;
      double scalingVDWB = it2->scalingVDW;
      double scalingCoulombB = it2->scalingCoulomb;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffFrameworkVDWSquared)
      {
        Potentials::PairDerivatives<0> energyFactor =
            Potentials::potentialVDW<0>(forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);

        energySum.frameworkMoleculeVDW += energyFactor.energy;
        energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, energyFactor.dUdlambda);
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        Potentials::PairDerivatives<0> energyFactor =
            Potentials::potentialCoulomb<0>(forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, energyFactor.dUdlambda);

        Potentials::PairDerivatives<1> gradientFactor =
            Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);
        std::size_t index = static_cast<std::size_t>(std::distance(moleculeAtoms.begin(), it2));
        electricFieldMolecules[index] += scalingCoulombA * chargeA * gradientFactor.firstDerivativeFactor * dr;
      }
    }
  }

  return energySum;
}

std::tuple<double, double3, double3x3> Interactions::calculateHessianAtPositionVDW(const ForceField& forceField,
                                                                                   const SimulationBox& simulationBox,
                                                                                   double3 posA, std::size_t typeA,
                                                                                   std::span<const Atom> frameworkAtoms)
{
  const double cutOffFrameworkVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;

  double energy{0.0};
  double3 first_derivative{};
  double3x3 second_derivative{};

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posB = it1->position;
    std::size_t typeB = static_cast<std::size_t>(it1->type);
    double scalingB = it1->scalingVDW;

    double3 dr = posA - posB;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = double3::dot(dr, dr);

    if (rr < cutOffFrameworkVDWSquared)
    {
      Potentials::PairDerivatives<2> v = Potentials::potentialVDW<2>(forceField, 1.0, scalingB, rr, typeA, typeB);

      energy += v.energy;

      first_derivative.x += dr.x * v.firstDerivativeFactor;
      first_derivative.y += dr.y * v.firstDerivativeFactor;
      first_derivative.z += dr.z * v.firstDerivativeFactor;

      // add contribution to the second derivatives (Hessian matrix)
      second_derivative.ax += v.secondDerivativeFactor * dr.x * dr.x + v.firstDerivativeFactor;
      second_derivative.ay += v.secondDerivativeFactor * dr.x * dr.y;
      second_derivative.az += v.secondDerivativeFactor * dr.x * dr.z;

      second_derivative.bx += v.secondDerivativeFactor * dr.y * dr.x;
      second_derivative.by += v.secondDerivativeFactor * dr.y * dr.y + v.firstDerivativeFactor;
      second_derivative.bz += v.secondDerivativeFactor * dr.y * dr.z;

      second_derivative.cx += v.secondDerivativeFactor * dr.z * dr.x;
      second_derivative.cy += v.secondDerivativeFactor * dr.z * dr.y;
      second_derivative.cz += v.secondDerivativeFactor * dr.z * dr.z + v.firstDerivativeFactor;
    }
  }

  return {energy, first_derivative, second_derivative};
}

std::tuple<double, double3, double3x3> Interactions::calculateHessianAtPositionCoulomb(
    const ForceField& forceField, const SimulationBox& simulationBox, double3 posA, double chargeA,
    std::span<const Atom> frameworkAtoms)
{
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  double energy{0.0};
  double3 first_derivative{};
  double3x3 second_derivative{};

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posB = it1->position;
    double scalingB = it1->scalingVDW;
    double chargeB = it1->charge;

    double3 dr = posA - posB;
    dr = simulationBox.applyPeriodicBoundaryConditions(dr);
    double rr = double3::dot(dr, dr);

    if (rr < cutOffChargeSquared)
    {
      double r = std::sqrt(rr);
      Potentials::PairDerivatives<2> v = Potentials::potentialCoulomb<2>(forceField, 1.0, scalingB, r, chargeA, chargeB);

      energy += v.energy;

      first_derivative.x += dr.x * v.firstDerivativeFactor;
      first_derivative.y += dr.y * v.firstDerivativeFactor;
      first_derivative.z += dr.z * v.firstDerivativeFactor;

      // add contribution to the second derivatives (Hessian matrix)
      second_derivative.ax += v.secondDerivativeFactor * dr.x * dr.x + v.firstDerivativeFactor;
      second_derivative.ay += v.secondDerivativeFactor * dr.x * dr.y;
      second_derivative.az += v.secondDerivativeFactor * dr.x * dr.z;

      second_derivative.bx += v.secondDerivativeFactor * dr.y * dr.x;
      second_derivative.by += v.secondDerivativeFactor * dr.y * dr.y + v.firstDerivativeFactor;
      second_derivative.bz += v.secondDerivativeFactor * dr.y * dr.z;

      second_derivative.cx += v.secondDerivativeFactor * dr.z * dr.x;
      second_derivative.cy += v.secondDerivativeFactor * dr.z * dr.y;
      second_derivative.cz += v.secondDerivativeFactor * dr.z * dr.z + v.firstDerivativeFactor;
    }
  }

  return {energy, first_derivative, second_derivative};
}
