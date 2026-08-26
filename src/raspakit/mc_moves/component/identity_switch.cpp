module;

module mc_moves_identity_switch;

import std;

import component;
import atom;
import molecule;
import double3;
import simd_quatd;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import cbmc_growth_context;
import cbmc_interactions;
import randomnumbers;
import system;
import running_energy;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_polarization;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::identitySwitchMove(RandomNumber &random, System &system,
                                                          std::size_t selectedComponent)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  const Move::Types move = Move::Types::IdentitySwitchCBMC;

  Component &componentAData = system.components[selectedComponent];
  if (componentAData.identitySwitches.empty())
  {
    return std::nullopt;
  }

  const std::size_t componentA = selectedComponent;
  const std::size_t componentB =
      componentAData.identitySwitches[random.uniform_integer(0, componentAData.identitySwitches.size() - 1)];

  if (componentB >= system.components.size() || componentB == componentA)
  {
    return std::nullopt;
  }

  Component &componentBData = system.components[componentB];
  if (componentBData.type != componentAData.type)
  {
    return std::nullopt;
  }

  // Composition is conserved (one A and one B are exchanged), so picking one molecule of each
  // uniformly gives identical proposal probabilities for the forward and reverse moves.
  if (system.numberOfIntegerMoleculesPerComponent[componentA] == 0 ||
      system.numberOfIntegerMoleculesPerComponent[componentB] == 0)
  {
    return std::nullopt;
  }

  componentAData.mc_moves_statistics.addTrial(move);
  componentBData.mc_moves_statistics.addTrial(move);

  const std::size_t moleculeA = system.randomIntegerMoleculeOfComponent(random, componentA);
  const std::size_t moleculeB = system.randomIntegerMoleculeOfComponent(random, componentB);

  std::span<Atom> atomsA = system.spanOfMolecule(componentA, moleculeA);
  std::span<Atom> atomsB = system.spanOfMolecule(componentB, moleculeB);

  // copies: the old configurations, needed after the molecules have been removed from storage and as
  // the background of the second retrace
  std::vector<Atom> oldA(atomsA.begin(), atomsA.end());
  std::vector<Atom> oldB(atomsB.begin(), atomsB.end());
  const Atom startingBeadA = oldA[componentAData.startingBead];
  const Atom startingBeadB = oldB[componentBData.startingBead];

  // Determine cutoff distances based on whether dual cutoff is used.
  const double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;

  // Distinct trial ids that cannot collide with any existing molecule (global ids are 0..N-1), so
  // that the molecule grown second sees the one grown before it in the background.
  const std::size_t trialIdNewB = system.numberOfMolecules();
  const std::size_t trialIdNewA = system.numberOfMolecules() + 1;

  const std::span<Atom> allMoleculeAtoms = system.spanOfMoleculeAtoms();

  // Molecule-molecule polarization has to write the field change of every surviving molecule back
  // into a storage-aligned array.  Record where each background atom came from while the background
  // is built, so that the scatter does not have to repeat the membership test below.
  const bool trackNeighborPolarization =
      system.forceField.computePolarization && !system.forceField.omitInterPolarization;
  std::vector<std::size_t> backgroundStorageIndex;

  // Background of the CBMC phase: every molecule (fractional ones as well) except the two that
  // exchange identities.  Both exchanged molecules are absent from the storage-backed background,
  // so no molecule needs to be excluded by id and the CBMC helpers are called with their default
  // (no) skip.  The first molecule of the pair is appended in place for the second operation; the
  // reserve below covers that append (the pair is not part of 'backgroundWithoutPair'), so the
  // vector never reallocates.
  std::vector<Atom> backgroundWithoutPair;
  backgroundWithoutPair.reserve(allMoleculeAtoms.size());
  if (trackNeighborPolarization) backgroundStorageIndex.reserve(allMoleculeAtoms.size());
  for (std::size_t component = 0; component < system.components.size(); ++component)
  {
    for (std::size_t molecule = 0; molecule < system.numberOfMoleculesPerComponent[component]; ++molecule)
    {
      if ((component == componentA && molecule == moleculeA) || (component == componentB && molecule == moleculeB))
      {
        continue;
      }
      std::span<Atom> atoms = system.spanOfMolecule(component, molecule);
      backgroundWithoutPair.insert(backgroundWithoutPair.end(), atoms.begin(), atoms.end());

      if (trackNeighborPolarization)
      {
        const std::size_t offset = static_cast<std::size_t>(atoms.data() - allMoleculeAtoms.data());
        for (std::size_t k = 0; k != atoms.size(); ++k) backgroundStorageIndex.push_back(offset + k);
      }
    }
  }
  const std::size_t backgroundWithoutPairSize = backgroundWithoutPair.size();

  auto makeContext = [&](std::span<const Atom> background)
  {
    return CBMC::GrowContext{system.hasExternalField,
                             system.forceField,
                             system.simulationBox,
                             system.interpolationGrids,
                             system.externalFieldInterpolationGrid,
                             system.framework,
                             system.spanOfFrameworkAtoms(),
                             background,
                             system.beta,
                             cutOffFrameworkVDW,
                             cutOffMoleculeVDW,
                             cutOffCoulomb};
  };

  // One side of the exchange: the existing molecule (retraced at its own site) together with the
  // new molecule of the same component, which is grown at the *other* molecule's starting bead.
  struct Exchange
  {
    Component *component;
    std::size_t componentId;
    std::span<Atom> oldAtoms;    // the existing molecule, in system storage
    std::vector<Atom> *oldCopy;  // pristine copy of the same molecule
    Atom growStartingBead;       // starting bead of the other molecule: where the new one is grown
    std::size_t trialMoleculeId;
    std::optional<ChainGrowData> grown;
    ChainRetraceData retraced;
  };

  Exchange exchangeA{&componentAData, componentA, atomsA, &oldA, startingBeadB, trialIdNewA, std::nullopt, {}};
  Exchange exchangeB{&componentBData, componentB, atomsB, &oldB, startingBeadA, trialIdNewB, std::nullopt, {}};

  // Rather than hiding both molecules from all four CBMC operations and adding their mutual
  // interaction back by hand, the pair is grown and retraced in a nested background: the first
  // molecule against 'backgroundWithoutPair' and the second against that background plus the
  // first.  Every intra-pair interaction then sits inside exactly one Rosenbluth weight on each
  // side of the move, at the full cut-offs and with the dual cut-off correction applied to it.
  //
  // The nesting order must be the same in the forward and the reverse move, otherwise the retrace
  // of a configuration no longer reproduces the weight that its growth produced and the Rosenbluth
  // acceptance is not exact.  Ordering by role ("the selected molecule first") does not qualify:
  // the exchange moves each component to the other site, so the reverse move would assign the
  // roles the other way around.  The component index is invariant under the exchange and under
  // which of the two components triggered the move, so it is used as the ordering key.
  Exchange *const ordered[2] = {componentA < componentB ? &exchangeA : &exchangeB,
                                componentA < componentB ? &exchangeB : &exchangeA};

  time_begin = std::chrono::steady_clock::now();
  bool constructed = true;
  for (std::size_t step = 0; step != 2 && constructed; ++step)
  {
    Exchange &exchange = *ordered[step];
    const CBMC::GrowContext growContext = makeContext(backgroundWithoutPair);

    exchange.grown = CBMC::growMoleculeIdentityChangeInsertion(random, growContext, *exchange.component,
                                                               exchange.componentId, exchange.trialMoleculeId,
                                                               exchange.growStartingBead, 1.0, 0, false);
    if (!exchange.grown ||
        system.insideBlockedPockets(*exchange.component, std::span<const Atom>(exchange.grown->atoms)))
    {
      constructed = false;
      break;
    }

    if (system.forceField.useDualCutOff)
    {
      // Dual cut-off scheme: correct the grown configuration from the inner cut-off to the full
      // cut-offs, using the same background as the growth.
      std::optional<RunningEnergy> correction =
          CBMC::computeDualCutOffCorrection(growContext, *exchange.component, exchange.grown->atoms);
      if (!correction.has_value())
      {
        constructed = false;
        break;
      }

      exchange.grown->energies += correction.value();
      exchange.grown->RosenbluthWeight *= std::exp(-system.beta * correction->potentialEnergy());
    }

    if (step == 0)
    {
      backgroundWithoutPair.insert(backgroundWithoutPair.end(), exchange.grown->atoms.begin(),
                                   exchange.grown->atoms.end());
    }
  }
  backgroundWithoutPair.resize(backgroundWithoutPairSize);

  for (std::size_t step = 0; step != 2 && constructed; ++step)
  {
    Exchange &exchange = *ordered[step];
    const CBMC::GrowContext retraceContext = makeContext(backgroundWithoutPair);

    exchange.retraced =
        CBMC::retraceMoleculeIdentityChangeDeletion(random, retraceContext, *exchange.component, exchange.oldAtoms);

    if (system.forceField.useDualCutOff)
    {
      // Dual cut-off scheme: correct the retraced configuration from the inner cut-off to the full
      // cut-offs, using the same background as the retrace.
      std::optional<RunningEnergy> correction =
          CBMC::computeDualCutOffCorrection(retraceContext, *exchange.component, *exchange.oldCopy);
      if (!correction.has_value())
      {
        // an existing configuration should never register as an overlap; reject defensively
        constructed = false;
        break;
      }

      exchange.retraced.energies += correction.value();
      exchange.retraced.RosenbluthWeight *= std::exp(-system.beta * correction->potentialEnergy());
    }

    if (step == 0)
    {
      backgroundWithoutPair.insert(backgroundWithoutPair.end(), exchange.oldCopy->begin(), exchange.oldCopy->end());
    }
  }
  backgroundWithoutPair.resize(backgroundWithoutPairSize);
  time_end = std::chrono::steady_clock::now();
  componentAData.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  if (!constructed)
  {
    return std::nullopt;
  }

  componentAData.mc_moves_statistics.addConstructed(move);
  componentBData.mc_moves_statistics.addConstructed(move);

  // Combined Ewald Fourier difference for the simultaneous change of both molecules
  // (single call so the cross terms between the two exchanged molecules are handled exactly).
  std::vector<Atom> newAtoms;
  newAtoms.reserve(exchangeB.grown->atoms.size() + exchangeA.grown->atoms.size());
  newAtoms.insert(newAtoms.end(), exchangeB.grown->atoms.begin(), exchangeB.grown->atoms.end());
  newAtoms.insert(newAtoms.end(), exchangeA.grown->atoms.begin(), exchangeA.grown->atoms.end());

  std::vector<Atom> oldAtoms;
  oldAtoms.reserve(oldA.size() + oldB.size());
  oldAtoms.insert(oldAtoms.end(), oldA.begin(), oldA.end());
  oldAtoms.insert(oldAtoms.end(), oldB.begin(), oldB.end());

  time_begin = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, newAtoms, oldAtoms, system.netCharge);
  time_end = std::chrono::steady_clock::now();
  componentAData.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Tail corrections cancel exactly: the multiset of atom types in the box is unchanged
  // (one A and one B on both sides of the move), for both molecule-molecule and
  // framework-molecule contributions.  Likewise the Ewald net-charge correction vanishes.

  std::vector<double3> electricFieldNeighborDelta;
  RunningEnergy polarizationDifference;
  std::vector<double3> new_electric_field(newAtoms.size());
  std::vector<double3> old_electric_field(oldAtoms.size());
  if (system.forceField.computePolarization)
  {
    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), new_electric_field,
                                                                  old_electric_field, newAtoms, oldAtoms);

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.trialEik, system.forceField, system.simulationBox, new_electric_field, old_electric_field, newAtoms,
        oldAtoms);

    // Real-space molecule-molecule polarization.  The energies these calls return are discarded: the
    // inter-molecular energy is already in the CBMC weights, only the electric fields are wanted.
    if (trackNeighborPolarization)
    {
      electricFieldNeighborDelta.assign(allMoleculeAtoms.size(), double3(0.0, 0.0, 0.0));

      // Old configuration: system storage is the right neighbor set, since each old molecule is
      // skipped against itself by molecule id while the partner it interacts with is still present.
      [[maybe_unused]] std::optional<RunningEnergy> oldFieldEnergy =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, electricFieldNeighborDelta, std::span<double3>{},
              old_electric_field, allMoleculeAtoms, {}, oldAtoms);

      // New configuration: the two new molecules sit on the sites of the old ones, so storage cannot
      // serve as their neighbor set (they would see the molecules they replace at zero distance).
      // 'backgroundWithoutPair' holds exactly the surviving molecules; its field change is scattered
      // back into the storage-aligned array through the index recorded while it was built.
      std::vector<double3> backgroundDelta(backgroundWithoutPairSize, double3(0.0, 0.0, 0.0));
      [[maybe_unused]] std::optional<RunningEnergy> newFieldEnergy =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, backgroundDelta, new_electric_field, std::span<double3>{},
              backgroundWithoutPair, newAtoms, {});

      for (std::size_t i = 0; i != backgroundWithoutPairSize; ++i)
      {
        electricFieldNeighborDelta[backgroundStorageIndex[i]] += backgroundDelta[i];
      }

      // Field the two new molecules exert on each other: neither is in storage yet, and their trial
      // ids differ, so one pass of the pair against itself gives each of them the field of the other.
      std::vector<double3> mutualField(newAtoms.size(), double3(0.0, 0.0, 0.0));
      std::vector<double3> unusedDelta(newAtoms.size(), double3(0.0, 0.0, 0.0));
      [[maybe_unused]] std::optional<RunningEnergy> mutualFieldEnergy =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, unusedDelta, mutualField, std::span<double3>{}, newAtoms,
              newAtoms, {});
      for (std::size_t i = 0; i != newAtoms.size(); ++i)
      {
        new_electric_field[i] += mutualField[i];
      }

      // The two replaced molecules must not appear in the neighbor sum: their own field change is
      // carried by 'old_electric_field' and 'new_electric_field'.
      const std::size_t offsetA = static_cast<std::size_t>(atomsA.data() - allMoleculeAtoms.data());
      for (std::size_t k = 0; k != atomsA.size(); ++k) electricFieldNeighborDelta[offsetA + k] = double3(0.0, 0.0, 0.0);
      const std::size_t offsetB = static_cast<std::size_t>(atomsB.data() - allMoleculeAtoms.data());
      for (std::size_t k = 0; k != atomsB.size(); ++k) electricFieldNeighborDelta[offsetB + k] = double3(0.0, 0.0, 0.0);
    }

    polarizationDifference = Interactions::computePolarizationEnergyDifference(system.forceField, new_electric_field,
                                                                               old_electric_field, newAtoms, oldAtoms);

    if (trackNeighborPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta, allMoleculeAtoms);
    }
  }

  // Canonical acceptance: composition is conserved, so fugacities, ideal-gas Rosenbluth weights and
  // the N_old/(N_new+1) factors of the semi-grand move all cancel identically.  What remains is the
  // Rosenbluth ratio of the two grows over the two retraces, corrected for the energy terms that the
  // CBMC weights do not contain (Fourier-space Ewald and polarization).  The interaction between the
  // two exchanged molecules needs no separate term: the nested background puts it inside the second
  // grow and the second retrace.
  const double correctionFactor =
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + polarizationDifference.potentialEnergy()));

  const double acceptanceProbability =
      correctionFactor * (exchangeB.grown->RosenbluthWeight * exchangeA.grown->RosenbluthWeight) /
      (exchangeA.retraced.RosenbluthWeight * exchangeB.retraced.RosenbluthWeight);

  if (random.uniform() < acceptanceProbability)
  {
    componentAData.mc_moves_statistics.addAccepted(move);
    componentBData.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::vector<Atom> acceptedAtomsNewB(exchangeB.grown->atoms.begin(), exchangeB.grown->atoms.end());
    for (Atom &atom : acceptedAtomsNewB)
    {
      atom.componentId = static_cast<std::uint8_t>(componentB);
    }
    Molecule acceptedMoleculeNewB = exchangeB.grown->molecule;
    acceptedMoleculeNewB.componentId = componentB;

    std::vector<Atom> acceptedAtomsNewA(exchangeA.grown->atoms.begin(), exchangeA.grown->atoms.end());
    for (Atom &atom : acceptedAtomsNewA)
    {
      atom.componentId = static_cast<std::uint8_t>(componentA);
    }
    Molecule acceptedMoleculeNewA = exchangeA.grown->molecule;
    acceptedMoleculeNewA.componentId = componentA;

    // Apply the field changes on the surrounding molecules while the storage layout is still the one
    // 'electricFieldNeighborDelta' was built against, i.e. before any molecule is removed or added.
    if (trackNeighborPolarization)
    {
      std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
      for (std::size_t i = 0; i < storedElectricField.size(); ++i)
      {
        storedElectricField[i] += electricFieldNeighborDelta[i];
      }
    }

    // Replace A first, then B.  Deleting a molecule of component A shifts atom storage but leaves
    // component B's per-component molecule index unchanged (insertions append at the end of a
    // component's block), so 'moleculeB' still refers to the original molecule; its span, however,
    // must be re-fetched because the underlying atom vector has been modified.
    if (system.forceField.computePolarization)
    {
      std::span<double3> electricFieldNewB = std::span(new_electric_field.begin(), exchangeB.grown->atoms.size());
      std::span<double3> electricFieldNewA =
          std::span(new_electric_field.begin() +
                        static_cast<std::vector<double3>::difference_type>(exchangeB.grown->atoms.size()),
                    exchangeA.grown->atoms.size());

      system.deleteMolecule(componentA, moleculeA, atomsA);
      system.insertMoleculePolarization(componentB, acceptedMoleculeNewB, acceptedAtomsNewB, electricFieldNewB);

      std::span<Atom> refetchedAtomsB = system.spanOfMolecule(componentB, moleculeB);
      system.deleteMolecule(componentB, moleculeB, refetchedAtomsB);
      system.insertMoleculePolarization(componentA, acceptedMoleculeNewA, acceptedAtomsNewA, electricFieldNewA);
    }
    else
    {
      system.deleteMolecule(componentA, moleculeA, atomsA);
      system.insertMolecule(componentB, acceptedMoleculeNewB, acceptedAtomsNewB);

      std::span<Atom> refetchedAtomsB = system.spanOfMolecule(componentB, moleculeB);
      system.deleteMolecule(componentB, moleculeB, refetchedAtomsB);
      system.insertMolecule(componentA, acceptedMoleculeNewA, acceptedAtomsNewA);
    }

    const RunningEnergy energyDifference = (exchangeB.grown->energies + exchangeA.grown->energies) -
                                           (exchangeA.retraced.energies + exchangeB.retraced.energies) +
                                           energyFourierDifference + polarizationDifference;

    return energyDifference;
  }

  return std::nullopt;
}
