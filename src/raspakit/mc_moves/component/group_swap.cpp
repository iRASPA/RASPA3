module;

module mc_moves_group_swap;

import std;

import double3;
import component;
import molecule;
import atom;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import cbmc_interactions;
import randomnumbers;
import system;
import energy_status;
import running_energy;
import forcefield;
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

// Group insertion/deletion move for overall-neutral groups of charged molecules (a generalization
// of the Orkoulas-Panagiotopoulos ion-pair swap to an arbitrary number of satellite molecules).
//
// A group consists of one molecule of the defining ("central") component plus the satellite
// molecules listed in Component::groupComponentIds (repetitions allowed, e.g. Ca2+ with
// groupComponentIds = {Cl-, Cl-}). Insertion places the central molecule at a random position in
// the box and then grows every satellite, one after the other, with its starting bead inside the
// sphere of radius R_max = Component::maximumGroupDistance around the central starting bead. Each
// satellite is grown in a background that already contains the previously grown group members, so
// all intra-group interactions are counted exactly once. Deletion is the exact reverse: a random
// integer central molecule is chosen and, for every satellite slot in the listed order, a partner
// is chosen uniformly among the integer molecules of that slot's component whose starting bead
// lies within R_max of the central bead (sequentially, without replacement, excluding the central
// molecule itself).
//
// Detailed balance. With unordered configurations the grand-canonical weight ratio for adding the
// group is (beta f_0 V) * prod_c (beta f_c V)^{n_c} exp(-beta dU) (f = fugacity, n_c = number of
// satellites of component c, the central component's own satellites included in its n_c). The
// insertion proposal density is 1/V for the central position and 1/V_s per satellite
// (V_s = 4 pi R_max^3 / 3) for uniform-in-sphere placement; the reverse deletion selects the
// central molecule with probability 1/N_0(new) (N_0(new) = all integer molecules of the central
// component in the inserted state) and the satellite of slot t of component c with probability
// 1/(k_c - t), where k_c is the number of in-range integer molecules of component c around the
// central bead in the inserted state (central molecule excluded; every inserted satellite is in
// range by construction). The generation multiplicity n_c! of identical satellites cancels against
// the reverse-selection multiplicity, leaving
//
//   acc(ins) = min(1, [beta f_0 V / N_0(new)] prod_slots [beta f_c V_s b / (k_c - t)]
//                     * prod_i (W_i / W_i^IG) * exp(-beta dU_corr))
//
// and the deletion acceptance is the exact reciprocal evaluated in the current state. The factor
// b is the radial-proposal bias: b = 1 for uniform-in-sphere placement (conventional variant,
// r = R_max cbrt(u)) and b = 3 r^2 / R_max^2 for uniform-in-r placement (CBMC variant,
// r = R_max u). W_i / W_i^IG are the Rosenbluth-weight ratios of the grown/retraced molecules and
// dU_corr collects the Ewald, tail-correction, and polarization energies outside those weights.

struct GroupMemberView
{
  std::size_t componentId;
  std::size_t moleculeId;  // molecule index within the component (deletion only)
  std::span<Atom> atoms;
};

// number of integer molecules of 'componentId' whose starting bead lies within R_max of 'position'
// (minimum image), excluding the molecule 'excludedMolecule' (pass a value >= the number of
// molecules to exclude none)
static std::size_t countInRangeIntegerMolecules(const System& system, std::size_t componentId,
                                                const double3& position, double maximumDistance,
                                                std::size_t excludedMolecule)
{
  const std::size_t startingBead = system.components[componentId].startingBead;
  const double maximumDistanceSquared = maximumDistance * maximumDistance;
  std::size_t count = 0;

  const std::size_t firstIntegerMolecule = system.numberOfFractionalMoleculesPerComponent[componentId];
  for (std::size_t molecule = firstIntegerMolecule; molecule < system.numberOfMoleculesPerComponent[componentId];
       ++molecule)
  {
    if (molecule == excludedMolecule) continue;
    const std::span<const Atom> atoms = system.spanOfMolecule(componentId, molecule);
    const double3 dr = system.simulationBox.applyPeriodicBoundaryConditions(position - atoms[startingBead].position);
    if (dr.length_squared() <= maximumDistanceSquared) ++count;
  }
  return count;
}

static bool validGroupDefinition(const System& system, const Component& centralComponent)
{
  if (centralComponent.groupComponentIds.empty() || !centralComponent.maximumGroupDistance.has_value())
  {
    return false;
  }
  if (centralComponent.maximumGroupDistance.value() <= 0.0)
  {
    return false;
  }
  for (std::size_t satelliteComponentId : centralComponent.groupComponentIds)
  {
    if (satelliteComponentId >= system.components.size())
    {
      return false;
    }
  }
  return true;
}

static std::pair<std::optional<RunningEnergy>, double3> groupInsertion(RandomNumber& random, System& system,
                                                                       std::size_t selectedComponent,
                                                                       const Move::Types move)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Component& centralComponent = system.components[selectedComponent];

  if (!validGroupDefinition(system, centralComponent))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const std::vector<std::size_t>& satelliteComponentIds = centralComponent.groupComponentIds;
  const std::size_t numberOfSatellites = satelliteComponentIds.size();
  const double R_max = centralComponent.maximumGroupDistance.value();
  const double sphereVolume = (4.0 / 3.0) * std::numbers::pi * R_max * R_max * R_max;
  const bool distanceBiased = (move == Move::Types::GroupSwapCBMC);

  centralComponent.mc_moves_statistics.addTrial(move, 0);

  // Determine cutoff distances based on whether dual cutoff is used.
  const double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;

  const CBMC::GrowContext growContextCentral{
      system.hasExternalField,     system.forceField,           system.simulationBox,
      system.interpolationGrids,   system.externalFieldInterpolationGrid,
      system.framework,            system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.beta,                 cutOffFrameworkVDW,          cutOffMoleculeVDW,
      cutOffCoulomb};

  time_begin = std::chrono::steady_clock::now();
  std::optional<ChainGrowData> growDataCentral = CBMC::growMoleculeSwapInsertion(
      random, growContextCentral, centralComponent, selectedComponent, system.numberOfMolecules(), 1.0,
      std::uint8_t{0}, false);
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  centralComponent.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  if (!growDataCentral) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  if (system.insideBlockedPockets(centralComponent, std::span<const Atom>(growDataCentral->atoms)))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  if (system.forceField.useDualCutOff)
  {
    // Dual cut-off scheme: correct the central molecule from the inner cut-off to the full cut-offs.
    std::optional<RunningEnergy> correction =
        CBMC::computeDualCutOffCorrection(growContextCentral, centralComponent, growDataCentral->atoms);
    if (!correction.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    growDataCentral->energies += correction.value();
    growDataCentral->RosenbluthWeight *= std::exp(-system.beta * correction->potentialEnergy());
  }

  const double3 centralPosition = growDataCentral->atoms[centralComponent.startingBead].position;

  // Background used to grow the satellites: existing molecules plus the already grown group
  // members. 'memberBackgroundSizes[i]' is the prefix length of 'background' that group member i
  // was grown against, so that every intra-group interaction is counted exactly once.
  std::vector<Atom> background(system.spanOfMoleculeAtoms().begin(), system.spanOfMoleculeAtoms().end());
  std::vector<std::size_t> memberBackgroundSizes;
  memberBackgroundSizes.reserve(1 + numberOfSatellites);
  memberBackgroundSizes.push_back(background.size());
  background.insert(background.end(), growDataCentral->atoms.begin(), growDataCentral->atoms.end());

  std::vector<ChainGrowData> satelliteGrowData;
  satelliteGrowData.reserve(numberOfSatellites);
  std::vector<double> distanceBiasFactors;
  distanceBiasFactors.reserve(numberOfSatellites);

  for (std::size_t j = 0; j < numberOfSatellites; ++j)
  {
    const std::size_t satelliteComponentId = satelliteComponentIds[j];
    Component& satelliteComponent = system.components[satelliteComponentId];

    const double r = distanceBiased ? R_max * random.uniform() : R_max * std::cbrt(random.uniform());
    const double3 direction = random.UnitSphere();
    const double3 fixedFirstBeadPosition = centralPosition + r * direction;

    memberBackgroundSizes.push_back(background.size());

    const CBMC::GrowContext growContext{
        system.hasExternalField,   system.forceField,             system.simulationBox,
        system.interpolationGrids, system.externalFieldInterpolationGrid,
        system.framework,          system.spanOfFrameworkAtoms(), background,
        system.beta,               cutOffFrameworkVDW,            cutOffMoleculeVDW,
        cutOffCoulomb};

    time_begin = std::chrono::steady_clock::now();
    std::optional<ChainGrowData> growData = CBMC::growMoleculePairSecondSwapInsertion(
        random, growContext, satelliteComponent, satelliteComponentId, system.numberOfMolecules() + 1 + j,
        fixedFirstBeadPosition, 1.0, std::uint8_t{0}, false);
    time_end = std::chrono::steady_clock::now();
    system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
    centralComponent.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

    if (!growData) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    if (system.insideBlockedPockets(satelliteComponent, std::span<const Atom>(growData->atoms)))
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    if (system.forceField.useDualCutOff)
    {
      // Correct the satellite from the inner cut-off to the full cut-offs, using the same
      // background (existing molecules plus previously grown group members) as the growth.
      std::optional<RunningEnergy> correction =
          CBMC::computeDualCutOffCorrection(growContext, satelliteComponent, growData->atoms);
      if (!correction.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

      growData->energies += correction.value();
      growData->RosenbluthWeight *= std::exp(-system.beta * correction->potentialEnergy());
    }

    distanceBiasFactors.push_back(distanceBiased ? 3.0 * r * r / (R_max * R_max) : 1.0);
    background.insert(background.end(), growData->atoms.begin(), growData->atoms.end());
    satelliteGrowData.push_back(std::move(growData.value()));
  }

  centralComponent.mc_moves_statistics.addConstructed(move, 0);

  const std::size_t groupSize = 1 + numberOfSatellites;
  std::vector<GroupMemberView> members;
  members.reserve(groupSize);
  members.push_back({selectedComponent, 0, std::span<Atom>(growDataCentral->atoms)});
  for (std::size_t j = 0; j < numberOfSatellites; ++j)
  {
    members.push_back({satelliteComponentIds[j], 0, std::span<Atom>(satelliteGrowData[j].atoms)});
  }

  // Ewald Fourier energy: sequential differences, each molecule against the state that already
  // contains the previously added group members.
  const auto savedStoredEik = system.storedEik;
  const auto savedTrialEik = system.trialEik;

  time_begin = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifference{};
  double runningNetCharge = system.netCharge;
  for (const GroupMemberView& member : members)
  {
    energyFourierDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
        system.simulationBox, member.atoms, {}, runningNetCharge);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    runningNetCharge += system.components[member.componentId].netCharge;
  }
  system.storedEik = savedStoredEik;
  system.trialEik = savedTrialEik;
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  centralComponent.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Tail corrections: each member against the same background prefix it was grown against, so the
  // intra-group cross terms are counted exactly once.
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy tailEnergyDifference{};
  for (std::size_t i = 0; i < groupSize; ++i)
  {
    const std::span<const Atom> backgroundSlice(background.data(), memberBackgroundSizes[i]);
    tailEnergyDifference += Interactions::computeInterMolecularTailEnergyDifference(
                                system.forceField, system.simulationBox, backgroundSlice, members[i].atoms, {}) +
                            Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(),
                                members[i].atoms, {});
  }
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);
  centralComponent.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  std::vector<std::vector<double3>> newElectricFields(groupSize);
  for (std::size_t i = 0; i < groupSize; ++i)
  {
    newElectricFields[i].assign(members[i].atoms.size(), double3(0.0, 0.0, 0.0));
  }

  std::vector<double3> electricFieldNeighborDelta;
  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    for (std::size_t i = 0; i < groupSize; ++i)
    {
      Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                    system.spanOfFrameworkAtoms(),
                                                                    newElectricFields[i], {}, members[i].atoms, {});

      Interactions::computeEwaldFourierElectricFieldDifference(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
          system.trialEik, system.forceField, system.simulationBox, newElectricFields[i], {}, members[i].atoms, {});
    }

    if (!system.forceField.omitInterPolarization)
    {
      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));

      // Field on the inserted molecules from the existing molecules, and the reciprocal field
      // change on those existing molecules (the inter-molecular energy is already captured by the
      // CBMC growth).
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        [[maybe_unused]] std::optional<RunningEnergy> energy =
            Interactions::computeInterMolecularPolarizationElectricFieldDifference(
                system.forceField, system.simulationBox, electricFieldNeighborDelta, newElectricFields[i],
                std::span<double3>{}, system.spanOfMoleculeAtoms(), members[i].atoms, {});
      }

      // Mutual fields between the group members (their interaction energy is already included
      // through the sequential CBMC growth).
      for (std::size_t i = 0; i < groupSize; ++i)
      {
        for (std::size_t j = i + 1; j < groupSize; ++j)
        {
          std::vector<double3> fieldOnI(members[i].atoms.size());
          std::vector<double3> fieldOnJ(members[j].atoms.size());
          [[maybe_unused]] std::optional<RunningEnergy> energy =
              Interactions::computeInterMolecularPolarizationElectricFieldDifference(
                  system.forceField, system.simulationBox, fieldOnI, fieldOnJ, std::span<double3>{}, members[i].atoms,
                  members[j].atoms, {});
          for (std::size_t k = 0; k < fieldOnI.size(); ++k) newElectricFields[i][k] += fieldOnI[k];
          for (std::size_t k = 0; k < fieldOnJ.size(); ++k) newElectricFields[j][k] += fieldOnJ[k];
        }
      }
    }

    for (std::size_t i = 0; i < groupSize; ++i)
    {
      polarizationDifference +=
          Interactions::computePolarizationEnergyDifference(system.forceField, newElectricFields[i], {},
                                                            members[i].atoms);
    }

    if (!system.forceField.omitInterPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  const double correctionFactorEwald =
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy() +
                               polarizationDifference.potentialEnergy()));

  // The reverse deletion selects the central molecule uniformly among all integer molecules of the
  // central component in the inserted state, which includes any satellites of the same component.
  std::size_t insertedCentralComponentMolecules = 1;
  for (std::size_t satelliteComponentId : satelliteComponentIds)
  {
    if (satelliteComponentId == selectedComponent) ++insertedCentralComponentMolecules;
  }

  const double fugacityCentral =
      centralComponent.molFraction * centralComponent.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double idealGasCentral = centralComponent.idealGasRosenbluthWeight.value_or(1.0);
  const double N_central = double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]);

  double Pacc = correctionFactorEwald *
                (system.beta * fugacityCentral * system.simulationBox.volume /
                 (N_central + double(insertedCentralComponentMolecules))) *
                (growDataCentral->RosenbluthWeight / idealGasCentral);

  // Per-slot reverse-selection candidate counts in the inserted state: the in-range integer
  // molecules of the slot's component (all inserted satellites of that component are in range by
  // construction; the central molecule is excluded), minus the satellites selected in earlier
  // slots of the same component.
  std::vector<std::size_t> totalSatellitesOfComponent(system.components.size(), 0);
  for (std::size_t satelliteComponentId : satelliteComponentIds) ++totalSatellitesOfComponent[satelliteComponentId];

  std::vector<std::optional<std::size_t>> existingInRangeCache(system.components.size());
  std::vector<std::size_t> earlierSlotsOfComponent(system.components.size(), 0);
  for (std::size_t j = 0; j < numberOfSatellites; ++j)
  {
    const std::size_t satelliteComponentId = satelliteComponentIds[j];
    const Component& satelliteComponent = system.components[satelliteComponentId];

    if (!existingInRangeCache[satelliteComponentId].has_value())
    {
      existingInRangeCache[satelliteComponentId] = countInRangeIntegerMolecules(
          system, satelliteComponentId, centralPosition, R_max, std::numeric_limits<std::size_t>::max());
    }

    const std::size_t reverseCandidates = existingInRangeCache[satelliteComponentId].value() +
                                          totalSatellitesOfComponent[satelliteComponentId] -
                                          earlierSlotsOfComponent[satelliteComponentId];
    ++earlierSlotsOfComponent[satelliteComponentId];

    const double fugacitySatellite =
        satelliteComponent.molFraction * satelliteComponent.fugacityCoefficient.value_or(1.0) * system.pressure;
    const double idealGasSatellite = satelliteComponent.idealGasRosenbluthWeight.value_or(1.0);

    Pacc *= (system.beta * fugacitySatellite * sphereVolume * distanceBiasFactors[j] / double(reverseCandidates)) *
            (satelliteGrowData[j].RosenbluthWeight / idealGasSatellite);
  }

  const std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

  if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN >= system.tmmc.maxMacrostate)
  {
    return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
  }

  const std::size_t newN = oldN == std::numeric_limits<std::size_t>::max() ? oldN : oldN + 1;
  const double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    centralComponent.mc_moves_statistics.addAccepted(move, 0);

    // Commit the field changes on the existing molecules before appending the group (and its fields).
    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
      for (std::size_t i = 0; i < storedElectricField.size(); ++i)
      {
        storedElectricField[i] += electricFieldNeighborDelta[i];
      }
    }

    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                               system.storedEik, system.trialEik, system.forceField,
                                               system.simulationBox, members[0].atoms, {});
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    system.insertMoleculePolarization(selectedComponent, growDataCentral->molecule, growDataCentral->atoms,
                                      newElectricFields[0]);

    for (std::size_t j = 0; j < numberOfSatellites; ++j)
    {
      Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                 system.storedEik, system.trialEik, system.forceField,
                                                 system.simulationBox, members[1 + j].atoms, {});
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
      system.insertMoleculePolarization(satelliteComponentIds[j], satelliteGrowData[j].molecule,
                                        satelliteGrowData[j].atoms, newElectricFields[1 + j]);
    }

    RunningEnergy totalEnergyDifference =
        growDataCentral->energies + energyFourierDifference + tailEnergyDifference + polarizationDifference;
    for (const ChainGrowData& growData : satelliteGrowData)
    {
      totalEnergyDifference += growData.energies;
    }

    return {totalEnergyDifference, double3(0.0, 1.0 - Pacc, Pacc)};
  }

  return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}

static std::pair<std::optional<RunningEnergy>, double3> groupDeletion(RandomNumber& random, System& system,
                                                                      std::size_t selectedComponent,
                                                                      std::size_t selectedMolecule,
                                                                      const Move::Types move)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Component& centralComponent = system.components[selectedComponent];

  if (!validGroupDefinition(system, centralComponent))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  const std::vector<std::size_t>& satelliteComponentIds = centralComponent.groupComponentIds;
  const std::size_t numberOfSatellites = satelliteComponentIds.size();
  const double R_max = centralComponent.maximumGroupDistance.value();
  const double sphereVolume = (4.0 / 3.0) * std::numbers::pi * R_max * R_max * R_max;
  const bool distanceBiased = (move == Move::Types::GroupSwapCBMC);

  centralComponent.mc_moves_statistics.addTrial(move, 1);

  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0)
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  std::span<Atom> centralAtoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  const double3 centralPosition = centralAtoms[centralComponent.startingBead].position;

  // Sequential selection of the satellites: for every slot, choose uniformly among the in-range
  // integer molecules of the slot's component, excluding the central molecule and the satellites
  // selected in earlier slots (matching the reverse-count bookkeeping of the insertion move).
  std::vector<std::size_t> selectedSatelliteMolecules;
  selectedSatelliteMolecules.reserve(numberOfSatellites);
  std::vector<std::size_t> slotCandidateCounts;
  slotCandidateCounts.reserve(numberOfSatellites);
  std::vector<double> distanceBiasFactors;
  distanceBiasFactors.reserve(numberOfSatellites);

  const double maximumDistanceSquared = R_max * R_max;
  for (std::size_t j = 0; j < numberOfSatellites; ++j)
  {
    const std::size_t satelliteComponentId = satelliteComponentIds[j];
    const std::size_t startingBead = system.components[satelliteComponentId].startingBead;

    std::vector<std::size_t> candidates;
    const std::size_t firstIntegerMolecule = system.numberOfFractionalMoleculesPerComponent[satelliteComponentId];
    for (std::size_t molecule = firstIntegerMolecule;
         molecule < system.numberOfMoleculesPerComponent[satelliteComponentId]; ++molecule)
    {
      if (satelliteComponentId == selectedComponent && molecule == selectedMolecule) continue;

      bool alreadySelected = false;
      for (std::size_t earlier = 0; earlier < j; ++earlier)
      {
        if (satelliteComponentIds[earlier] == satelliteComponentId && selectedSatelliteMolecules[earlier] == molecule)
        {
          alreadySelected = true;
          break;
        }
      }
      if (alreadySelected) continue;

      const std::span<const Atom> atoms = system.spanOfMolecule(satelliteComponentId, molecule);
      const double3 dr =
          system.simulationBox.applyPeriodicBoundaryConditions(centralPosition - atoms[startingBead].position);
      if (dr.length_squared() <= maximumDistanceSquared)
      {
        candidates.push_back(molecule);
      }
    }

    if (candidates.empty())
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    const std::size_t chosen = candidates[random.uniform_integer(0, candidates.size() - 1)];
    const std::span<const Atom> chosenAtoms = system.spanOfMolecule(satelliteComponentId, chosen);
    const double3 dr = system.simulationBox.applyPeriodicBoundaryConditions(centralPosition -
                                                                            chosenAtoms[startingBead].position);
    const double r = dr.length();
    if (r <= 0.0 || r > R_max)
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    selectedSatelliteMolecules.push_back(chosen);
    slotCandidateCounts.push_back(candidates.size());
    distanceBiasFactors.push_back(distanceBiased ? 3.0 * r * r / (R_max * R_max) : 1.0);
  }

  const std::size_t groupSize = 1 + numberOfSatellites;
  std::vector<GroupMemberView> members;
  members.reserve(groupSize);
  members.push_back({selectedComponent, selectedMolecule, centralAtoms});
  for (std::size_t j = 0; j < numberOfSatellites; ++j)
  {
    members.push_back({satelliteComponentIds[j], selectedSatelliteMolecules[j],
                       system.spanOfMolecule(satelliteComponentIds[j], selectedSatelliteMolecules[j])});
  }

  // Determine cutoff distances based on whether dual cutoff is used.
  const double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;

  // Background without the group, then cumulative prefixes exactly mirroring the insertion order:
  // member i is retraced against 'backgroundWithoutGroup' plus the group members 0..i-1.
  std::vector<Atom> cumulativeBackground;
  cumulativeBackground.reserve(system.spanOfMoleculeAtoms().size());
  for (std::size_t component = 0; component < system.components.size(); ++component)
  {
    // include all molecules (fractional ones as well) except the group that is being deleted
    for (std::size_t molecule = 0; molecule < system.numberOfMoleculesPerComponent[component]; ++molecule)
    {
      bool isGroupMember = false;
      for (const GroupMemberView& member : members)
      {
        if (member.componentId == component && member.moleculeId == molecule)
        {
          isGroupMember = true;
          break;
        }
      }
      if (isGroupMember) continue;

      std::span<Atom> atoms = system.spanOfMolecule(component, molecule);
      cumulativeBackground.insert(cumulativeBackground.end(), atoms.begin(), atoms.end());
    }
  }

  std::vector<std::size_t> memberBackgroundSizes;
  memberBackgroundSizes.reserve(groupSize);
  for (const GroupMemberView& member : members)
  {
    memberBackgroundSizes.push_back(cumulativeBackground.size());
    cumulativeBackground.insert(cumulativeBackground.end(), member.atoms.begin(), member.atoms.end());
  }

  // Retrace the group members, each against the same background prefix it was grown against.
  std::vector<ChainRetraceData> retraceData;
  retraceData.reserve(groupSize);
  for (std::size_t i = 0; i < groupSize; ++i)
  {
    const Component& memberComponent = system.components[members[i].componentId];
    const std::span<const Atom> backgroundSlice(cumulativeBackground.data(), memberBackgroundSizes[i]);
    const CBMC::GrowContext retraceContext{
        system.hasExternalField,   system.forceField,             system.simulationBox,
        system.interpolationGrids, system.externalFieldInterpolationGrid,
        system.framework,          system.spanOfFrameworkAtoms(), backgroundSlice,
        system.beta,               cutOffFrameworkVDW,            cutOffMoleculeVDW,
        cutOffCoulomb};

    std::span<Atom> memberAtoms = system.spanOfMolecule(members[i].componentId, members[i].moleculeId);

    time_begin = std::chrono::steady_clock::now();
    ChainRetraceData retrace =
        (i == 0) ? CBMC::retraceMoleculeSwapDeletion(random, retraceContext, memberComponent, memberAtoms)
                 : CBMC::retraceMoleculePairSecondSwapDeletion(retraceContext, memberComponent, memberAtoms);
    time_end = std::chrono::steady_clock::now();
    system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
    centralComponent.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

    if (system.forceField.useDualCutOff)
    {
      // Dual cut-off scheme: correct the retraced configuration from the inner cut-off to the full
      // cut-offs, using the same background as the retrace.
      std::vector<Atom> memberAtomsCopy(memberAtoms.begin(), memberAtoms.end());
      std::optional<RunningEnergy> correction =
          CBMC::computeDualCutOffCorrection(retraceContext, memberComponent, memberAtomsCopy);
      if (!correction.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

      retrace.energies += correction.value();
      retrace.RosenbluthWeight *= std::exp(-system.beta * correction->potentialEnergy());
    }

    retraceData.push_back(std::move(retrace));
  }

  // Ewald Fourier energy: sequential removal differences, satellites in reverse insertion order,
  // the central molecule last; each difference is evaluated at the net charge before that removal.
  const auto savedStoredEik = system.storedEik;
  const auto savedTrialEik = system.trialEik;

  time_begin = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifference{};
  double runningNetCharge = system.netCharge;
  for (std::size_t i = groupSize; i-- > 0;)
  {
    energyFourierDifference += Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
        system.simulationBox, {}, members[i].atoms, runningNetCharge);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    runningNetCharge -= system.components[members[i].componentId].netCharge;
  }
  system.storedEik = savedStoredEik;
  system.trialEik = savedTrialEik;
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  centralComponent.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Tail corrections: exact mirror of the insertion bookkeeping.
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy tailEnergyDifference{};
  for (std::size_t i = 0; i < groupSize; ++i)
  {
    const std::span<const Atom> backgroundSlice(cumulativeBackground.data(), memberBackgroundSizes[i]);
    tailEnergyDifference += Interactions::computeInterMolecularTailEnergyDifference(
                                system.forceField, system.simulationBox, backgroundSlice, {}, members[i].atoms) +
                            Interactions::computeFrameworkMoleculeTailEnergyDifference(
                                system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), {},
                                members[i].atoms);
  }
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);
  centralComponent.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  // Polarization energy change when the group is removed. The removed molecules lose their own
  // polarization energy (evaluated with their full stored field) and, when inter-molecular
  // polarization is active, the field on all remaining molecules changes because the group
  // disappears. The field change on the remaining molecules is returned in
  // 'electricFieldNeighborDelta' so it can be committed on acceptance.
  std::vector<double3> electricFieldNeighborDelta;
  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    if (!system.forceField.omitInterPolarization)
    {
      for (const GroupMemberView& member : members)
      {
        std::span<double3> storedField = system.spanElectricFieldOld(member.componentId, member.moleculeId);
        polarizationDifference +=
            Interactions::computePolarizationEnergyDifference(system.forceField, {}, storedField, {}, member.atoms);
      }

      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
      for (const GroupMemberView& member : members)
      {
        std::vector<double3> removedField(member.atoms.size());
        [[maybe_unused]] std::optional<RunningEnergy> energy =
            Interactions::computeInterMolecularPolarizationElectricFieldDifference(
                system.forceField, system.simulationBox, electricFieldNeighborDelta, std::span<double3>{},
                removedField, system.spanOfMoleculeAtoms(), {}, member.atoms);
      }

      // The removed molecules must not appear in the neighbor sum; their own change is already
      // accounted for above.
      std::span<Atom> allAtoms = system.spanOfMoleculeAtoms();
      for (const GroupMemberView& member : members)
      {
        const std::size_t offset = static_cast<std::size_t>(member.atoms.data() - allAtoms.data());
        for (std::size_t k = 0; k < member.atoms.size(); ++k)
        {
          electricFieldNeighborDelta[offset + k] = double3(0.0, 0.0, 0.0);
        }
      }

      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
    else
    {
      for (const GroupMemberView& member : members)
      {
        std::vector<double3> oldElectricField(member.atoms.size());
        Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                      system.spanOfFrameworkAtoms(), {},
                                                                      oldElectricField, {}, member.atoms);
        Interactions::computeEwaldFourierElectricFieldDifference(
            system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
            system.trialEik, system.forceField, system.simulationBox, {}, oldElectricField, {}, member.atoms);
        polarizationDifference += Interactions::computePolarizationEnergyDifference(system.forceField, {},
                                                                                    oldElectricField, {},
                                                                                    member.atoms);
      }
    }
  }

  const double correctionFactorEwald =
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy() +
                               polarizationDifference.potentialEnergy()));

  const double fugacityCentral =
      centralComponent.molFraction * centralComponent.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double idealGasCentral = centralComponent.idealGasRosenbluthWeight.value_or(1.0);
  const double N_central = double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]);

  double Pacc = correctionFactorEwald * (N_central / (system.beta * fugacityCentral * system.simulationBox.volume)) *
                (idealGasCentral / retraceData[0].RosenbluthWeight);

  for (std::size_t j = 0; j < numberOfSatellites; ++j)
  {
    const Component& satelliteComponent = system.components[satelliteComponentIds[j]];
    const double fugacitySatellite =
        satelliteComponent.molFraction * satelliteComponent.fugacityCoefficient.value_or(1.0) * system.pressure;
    const double idealGasSatellite = satelliteComponent.idealGasRosenbluthWeight.value_or(1.0);

    Pacc *= (double(slotCandidateCounts[j]) /
             (system.beta * fugacitySatellite * sphereVolume * distanceBiasFactors[j])) *
            (idealGasSatellite / retraceData[1 + j].RosenbluthWeight);
  }

  const std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

  if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN <= system.tmmc.minMacrostate)
  {
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  const std::size_t newN = oldN == 0 ? 0 : oldN - 1;
  const double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

  centralComponent.mc_moves_statistics.addConstructed(move, 1);

  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    centralComponent.mc_moves_statistics.addAccepted(move, 1);

    // Commit the field changes on the remaining molecules (the removed molecules' own entries were zeroed).
    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
      for (std::size_t i = 0; i < storedElectricField.size(); ++i)
      {
        storedElectricField[i] += electricFieldNeighborDelta[i];
      }
    }

    // Delete the group members with, per component, descending molecule indices so that the
    // molecule indices of the members still to be deleted remain valid; the spans are re-acquired
    // before every deletion because each deletion shifts the atom storage.
    std::vector<std::pair<std::size_t, std::size_t>> deletionOrder;
    deletionOrder.reserve(groupSize);
    for (const GroupMemberView& member : members)
    {
      deletionOrder.push_back({member.componentId, member.moleculeId});
    }
    std::sort(deletionOrder.begin(), deletionOrder.end(), std::greater<>());

    for (const auto& [componentId, moleculeId] : deletionOrder)
    {
      std::span<Atom> atoms = system.spanOfMolecule(componentId, moleculeId);
      Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                 system.storedEik, system.trialEik, system.forceField,
                                                 system.simulationBox, {}, std::span<const Atom>(atoms));
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
      system.deleteMolecule(componentId, moleculeId, atoms);
    }

    RunningEnergy totalEnergyDifference = -energyFourierDifference - tailEnergyDifference - polarizationDifference;
    for (const ChainRetraceData& retrace : retraceData)
    {
      totalEnergyDifference += retrace.energies;
    }

    return {totalEnergyDifference, double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
}

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::groupInsertionMove(RandomNumber& random, System& system,
                                                                              std::size_t selectedComponent)
{
  return groupInsertion(random, system, selectedComponent, Move::Types::GroupSwap);
}

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::groupDeletionMove(RandomNumber& random, System& system,
                                                                             std::size_t selectedComponent,
                                                                             std::size_t selectedMolecule)
{
  return groupDeletion(random, system, selectedComponent, selectedMolecule, Move::Types::GroupSwap);
}

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::groupInsertionMoveCBMC(RandomNumber& random, System& system,
                                                                                  std::size_t selectedComponent)
{
  return groupInsertion(random, system, selectedComponent, Move::Types::GroupSwapCBMC);
}

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::groupDeletionMoveCBMC(RandomNumber& random, System& system,
                                                                                 std::size_t selectedComponent,
                                                                                 std::size_t selectedMolecule)
{
  return groupDeletion(random, system, selectedComponent, selectedMolecule, Move::Types::GroupSwapCBMC);
}
