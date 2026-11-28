module;

module reweighted_histogram;

import std;

import stringutils;
import hardware_info;
import archive;
import graceful_shutdown;
import system;
import framework;
import randomnumbers;
import input_reader;
import component;
import averages;
import property_loading;
import units;
import simulationbox;
import forcefield;
import vdwparameters;
import equation_of_states;
import energy_status;
import running_energy;
import atom;
import int3;
import double3;
import double3x3;
import property_lambda_probability_histogram;
import mc_moves;
import mc_moves_cputime;
import mc_moves_statistics;
import mc_moves_parallel_tempering_swap;
import system_preparation;
import json;
import isotherm_bet;
import nitrogen_bet_setup;
import molecule;
import interactions_internal;
import mc_moves_move_types;

// The analysis-property writers (RDFs, density grid, histograms, molecule properties) gate
// themselves on their own 'writeEvery'; a cycle argument of 0 forces the write (used for the
// final flush at the end of the run). The replica id keys the output filenames, so every replica
// writes its own set of files.
static void writeReplicaAnalysisOutputs(System& system, std::size_t replicaId, std::size_t cycle)
{
  if (system.propertyConventionalRadialDistributionFunction.has_value())
  {
    system.propertyConventionalRadialDistributionFunction->writeOutput(
        system.forceField, replicaId, system.simulationBox.volume, system.totalNumberOfPseudoAtoms, cycle);
  }
  if (system.propertyRadialDistributionFunction.has_value())
  {
    system.propertyRadialDistributionFunction->writeOutput(system.forceField, replicaId, system.simulationBox.volume,
                                                           system.totalNumberOfPseudoAtoms, cycle);
  }
  if (system.propertyDensityGrid.has_value())
  {
    system.propertyDensityGrid->writeOutput(replicaId, system.simulationBox, system.forceField, system.framework,
                                            system.components, cycle);
  }
  if (system.averageEnergyHistogram.has_value())
  {
    system.averageEnergyHistogram->writeOutput(replicaId, cycle);
  }
  if (system.averageNumberOfMoleculesHistogram.has_value())
  {
    system.averageNumberOfMoleculesHistogram->writeOutput(replicaId, system.components, cycle);
  }
  if (system.propertyMoleculeProperties.has_value())
  {
    system.propertyMoleculeProperties->writeOutput(replicaId, system.components, cycle);
  }
  // (sampled with the other properties of the replica, and written by the parallel-tempering driver only)
  if (system.propertyEndToEndACF.has_value())
  {
    system.propertyEndToEndACF->writeOutput(replicaId, system.components, cycle);
  }
  if (system.propertyMoleculeShape.has_value())
  {
    system.propertyMoleculeShape->writeOutput(replicaId, system.components, cycle);
  }
  if (system.propertyMoleculeBackbone.has_value())
  {
    system.propertyMoleculeBackbone->writeOutput(replicaId, system.components, cycle);
  }
}

// The template system must be one the reweighting analysis can describe. The sample of a replica is (N, U) with
// N the count of a single adsorbate and U the potential energy, and the reduced potential of a state is
// u_i = beta_i U - N ln(beta_i f_i); everything that this pair or that reduced potential does not carry is
// refused here, for the input reader and for the programmatic constructor alike.
//
// 'reweightsBeyondTheSample': a single state point that is reweighted to another temperature or pressure. Only a
// single state point whose results are written at that very state point needs no reduced potential at all.
static void validateReweightingTemplate(const System& system, std::span<const double> temperatures,
                                        std::size_t numberOfPressures, bool reweightsBeyondTheSample = false)
{
  if (system.components.empty())
  {
    throw std::runtime_error(
        "[ReweightedHistogram]: an adsorbate component is required (the reweighting histograms are collected over "
        "the molecule count of the first adsorbate)\n");
  }
  // the state points are (temperature, pressure) pairs: a chemical potential on a component would be
  // overwritten by the fugacity of each (the input reader refuses the keyword for this driver)
  for (const Component& component : system.components)
  {
    if (component.chemicalPotential.has_value())
    {
      throw std::runtime_error(std::format(
          "[ReweightedHistogram]: component '{}' has a chemical potential; the state points are given by "
          "pressures\n",
          component.name));
    }
  }
  // The replicas sample the grand-canonical Boltzmann distribution of their own state point. A transition-matrix
  // bias, or its macrostate bounds, makes it another distribution, which the reduced potential does not carry.
  if (system.tmmc.doTMMC)
  {
    throw std::runtime_error(
        "[ReweightedHistogram]: transition-matrix Monte Carlo ('UseBiasOnMacrostate', 'doTMMC') can not be combined "
        "with histogram reweighting: the samples have to be drawn from the unbiased grand-canonical ensemble\n");
  }
  // A reaction or an identity change turns molecules of the counted component into another one. The sample has
  // one count and one fugacity, and a reaction's own fractional molecules (their lambda is not looked at when a
  // sample is taken) would be binned at partial coupling as if they were whole.
  {
    auto probability = [&](Move::Types move)
    {
      double total = system.mc_moves_probabilities.getProbability(move);
      for (const Component& component : system.components)
      {
        total += component.mc_moves_probabilities.getProbability(move);
      }
      return total;
    };
    for (const Move::Types move :
         {Move::Types::ReactionCBMC, Move::Types::ReactionConventionalCFCMC, Move::Types::ReactionConventionalCBCFCMC,
          Move::Types::ReactionCFCMC, Move::Types::ReactionCBCFCMC, Move::Types::IdentityChangeCBMC})
    {
      if (probability(move) > 0.0)
      {
        throw std::runtime_error(std::format(
            "[ReweightedHistogram]: the move '{}' is not supported: the (N, U) samples carry the molecule count of "
            "one adsorbate at one fugacity, which a reaction or an identity change does not conserve\n",
            Move::moveNames[std::to_underlying(move)]));
      }
    }
    if (!system.reactions.list.empty())
    {
      throw std::runtime_error(
          "[ReweightedHistogram]: reactions are not supported: the (N, U) samples carry the molecule count of one "
          "adsorbate at one fugacity\n");
    }
  }
  // Across state points the reduced potential carries the molecule count of the adsorbate only: a second
  // swappable component would need its own -N_1 ln(beta f_1), whose fugacity changes along the ladder.
  // The same holds for one state point that is reweighted to another: a programmatic run of a single state with
  // a reweighting range passed these checks, and a binary mixture, or a system without an insertion move, was
  // then reweighted along the pressure with N_0 ln(beta f_0) alone.
  const bool severalStatePoints = temperatures.size() * numberOfPressures > 1uz || reweightsBeyondTheSample;
  const std::size_t numberOfSwappable = static_cast<std::size_t>(
      std::ranges::count_if(system.components, [](const Component& component) { return component.swappable; }));
  if (severalStatePoints && numberOfSwappable > 1uz)
  {
    throw std::runtime_error(std::format(
        "[ReweightedHistogram]: {} swappable components were given; the (N, U) samples record the molecule count of "
        "the first adsorbate only, so a ladder of state points can only be reweighted for a single adsorbate\n",
        numberOfSwappable));
  }
  // A pair or a group swap inserts the counted molecule together with others, at the activity of the pair or
  // the group: the sample has one count and one fugacity. Their fractional molecules are not carried by the
  // configuration swaps either, which refuse every exchange of such a system: a homodimer group
  // ('GroupComponents': [0, 0]) is a single swappable component, passed the checks above, and ran its ladder
  // without one accepted swap.
  for (const Component& component : system.components)
  {
    for (const Move::Types move : {Move::Types::PairSwap, Move::Types::PairSwapCBMC, Move::Types::PairSwapCFCMC,
                                   Move::Types::PairSwapCBCFCMC, Move::Types::GroupSwap, Move::Types::GroupSwapCBMC,
                                   Move::Types::GroupSwapCFCMC, Move::Types::GroupSwapCBCFCMC})
    {
      if (component.mc_moves_probabilities.getProbability(move) > 0.0)
      {
        throw std::runtime_error(std::format(
            "[ReweightedHistogram]: the move '{}' of component '{}' is not supported: the (N, U) samples carry the "
            "molecule count of one adsorbate at its own fugacity, not that of a pair or a group\n",
            Move::moveNames[std::to_underlying(move)], component.name));
      }
    }
  }
  // A fractional molecule pinned at a lambda (thermodynamic integration) is never decoupled, and a sample is
  // taken at lambda = 0 only: the run ended without one.
  for (const Component& component : system.components)
  {
    if (component.fixedLambdaBin.has_value())
    {
      throw std::runtime_error(std::format(
          "[ReweightedHistogram]: component '{}' has a fixed lambda ('LambdaBinIndex', 'ThermodynamicIntegration'); "
          "the (N, U) samples are taken with the fractional molecules decoupled, which a pinned one never is\n",
          component.name));
    }
  }
  for (const VDWParameters& parameter : system.forceField.data)
  {
    if (parameter.type == VDWParameters::Type::FeynmannHibbs)
    {
      throw std::runtime_error(
          std::format("[ReweightedHistogram]: temperature-dependent potentials (Feynman-Hibbs) are not "
                      "compatible with histogram reweighting (a sampled energy U(x) is reused at all "
                      "temperatures)\n"));
    }
  }
  if (system.framework.has_value() && system.framework->hasMobileAtoms())
  {
    throw std::runtime_error(
        "[ReweightedHistogram]: a framework with mobile (flexible) atoms is not supported: the configuration "
        "swaps exchange the adsorbates only, and a replica would keep its own framework with its partner's "
        "energies\n");
  }
  // Volume moves make V a sampled variable: the (N, U) sample does not record it, the reduced potential lacks
  // beta P V - N ln V, and the bulk coexistence analysis would take one replica's final volume for all.
  for (const auto& [move, name] : {std::pair{Move::Types::VolumeChange, "VolumeMoveProbability"},
                                   std::pair{Move::Types::AnisotropicVolumeChange, "AnisotropicVolumeMoveProbability"}})
  {
    if (system.mc_moves_probabilities.getProbability(move) > 0.0)
    {
      throw std::runtime_error(std::format(
          "[ReweightedHistogram]: '{}' is not supported: histogram reweighting is done in the grand-canonical "
          "ensemble at fixed volume\n",
          name));
    }
  }

  // the counted component is the adsorbate, which extra-framework cations of a fixed number may precede
  const Component& component = system.components[system.tmmcComponent()];
  const bool changesTheMoleculeCount = std::ranges::any_of(
      std::array{Move::Types::Swap, Move::Types::SwapCBMC, Move::Types::SwapCFCMC, Move::Types::SwapCBCFCMC,
                 Move::Types::SwapMultiTrialCBMC},
      [&](Move::Types move) { return component.mc_moves_probabilities.getProbability(move) > 0.0; });
  if (severalStatePoints && (!component.swappable || !changesTheMoleculeCount))
  {
    throw std::runtime_error(std::format(
        "[ReweightedHistogram]: component '{}' needs an insertion/deletion move ('SwapProbability', "
        "'CFCMC_SwapProbability' or 'CFCMC_CBMC_SwapProbability'): without one the molecule "
        "count never changes and the pressure ladder has nothing to reweight\n",
        component.name));
  }

  // (A flexible adsorbate across temperatures: the activity it is sampled at contains its intramolecular
  // configuration integral q(beta), which 'Preparation::prepareIdealGasReferences' computes for the sampled and the
  // reweighted temperatures when the replicas are made.)
}

// The numbers of the run itself, for the input reader and the programmatic constructor alike: the block count
// sizes the per-block histograms (zero blocks wrote past the end of an empty one), and a state point enters
// the reduced potentials through ln(beta f) (a zero pressure gave -inf, and 0 * -inf for the empty box).
static void validateReweightingRun(std::size_t numberOfBlocks, std::span<const double> temperatures,
                                   std::span<const double> pressures, std::span<const double> reweightingTemperatures,
                                   std::pair<double, double> reweightingPressureRange)
{
  if (numberOfBlocks == 0uz)
  {
    throw std::runtime_error("[ReweightedHistogram]: 'NumberOfBlocks' must be at least 1\n");
  }
  auto requirePositive = [](std::span<const double> values, std::string_view what)
  {
    for (const double value : values)
    {
      if (!(value > 0.0) || !(value < 1e300))
      {
        throw std::runtime_error(
            std::format("[ReweightedHistogram]: {} must be positive and finite (got {})\n", what, value));
      }
    }
  };
  requirePositive(temperatures, "every temperature of the ladder");
  requirePositive(pressures, "every pressure of the ladder");
  requirePositive(reweightingTemperatures, "every reweighting temperature");
  requirePositive(std::array{reweightingPressureRange.first, reweightingPressureRange.second},
                  "the reweighting pressure range");
}

// The intramolecular energy of the fractional molecules. At lambda = 0 a fractional molecule has no
// intermolecular interactions, but its bonded terms (and its intramolecular pairs) are not scaled by lambda and
// are inside runningEnergies: the configuration is then one of N molecules plus an ideal-gas molecule, whose
// energy does not belong to the N-molecule state the sample is recorded for.
static double fractionalMoleculesIntraEnergy(const System& system)
{
  double energy = 0.0;
  std::size_t index = 0uz;
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    const std::size_t numberOfFractional = system.numberOfFractionalMoleculesPerComponent[componentId];
    if (numberOfFractional > 0uz)
    {
      // the fractional molecules occupy the first slots of their component
      const std::span<const Molecule> fractional{&system.moleculeData[index], numberOfFractional};
      energy += Interactions::computeIntraMolecularEnergy(system.forceField, system.simulationBox,
                                                          system.components[componentId].intraMolecularPotentials,
                                                          fractional, system.spanOfMoleculeAtoms())
                    .potentialEnergy();
    }
    index += system.numberOfMoleculesPerComponent[componentId];
  }
  return energy;
}

ReweightedHistogram::ReweightedHistogram(InputReader& reader)
    : random(reader.randomSeed),
      numberOfProductionCycles(reader.numberOfProductionCycles),
      numberOfPreInitializationCycles(reader.numberOfPreInitializationCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(std::max(1uz, reader.printEvery)),
      optimizeMCMovesEvery(std::max(1uz, reader.optimizeMCMovesEvery)),
      rescaleWangLandauEvery(std::max(1uz, reader.rescaleWangLandauEvery)),
      writeBinaryRestartEvery(reader.writeBinaryRestartEvery),
      numberOfBlocks(reader.numberOfBlocks),
      parallelTemperingSwapEvery(reader.parallelTemperingSwapEvery),
      sampleReweightingEvery(std::max(1uz, reader.sampleReweightingEvery)),
      reweightingTemperatures(reader.reweightingTemperatures),
      reweightingNumberOfPressures(reader.reweightingNumberOfPressures),
      whamTolerance(reader.whamTolerance),
      whamMaximumIterations(reader.whamMaximumIterations),
      betScoutMaximumCycles(reader.betScoutMaximumCycles),
      computeBET(reader.computeBET),
      autoExternalPressures(reader.autoExternalPressures),
      autoReweightingPressureRange(reader.autoReweightingPressureRange),
      autoMacroStateMaximum(reader.autoMacroStateMaximum),
      temperatures(reader.parallelTemperingTemperatures),
      pressures(reader.parallelTemperingPressures),
      numberOfTemperatures(temperatures.size()),
      numberOfPressures(pressures.size()),
      numberOfReplicas(numberOfTemperatures * numberOfPressures)
{
  System templateSystem = std::move(reader.systems.front());
  reader.systems.clear();

  // On a binary restart the restart file carries everything built below: the scouted pressure plan and
  // filling ceiling, the pressure ladder, the reweighting span and the replicas themselves. The Henry and
  // occupancy scouts are minutes of GCMC whose only product the restart file would then overwrite, and a
  // ladder placed by the scout is not there to validate before the file is read.
  if (reader.restartFromBinary)
  {
    systems.push_back(std::move(templateSystem));
    return;
  }

  if (temperatures.empty())
  {
    temperatures.push_back(templateSystem.temperature);
  }
  // an automatic pressure ladder ('ExternalPressures': 'auto') is a ladder of several pressures
  validateReweightingTemplate(templateSystem, temperatures, autoExternalPressures ? 2uz : pressures.size());

  if (autoExternalPressures || autoReweightingPressureRange || autoMacroStateMaximum)
  {
    const NitrogenBETAutoSetupResult betSetup = applyNitrogenBETAutoSetup(
        templateSystem, temperatures.front(), reader.numberOfThreads, betScoutMaximumCycles,
        autoExternalPressures || autoReweightingPressureRange, autoMacroStateMaximum);
    nitrogenBETPressurePlan = std::move(betSetup.pressurePlan);
    nitrogenBETFillingCeiling = std::move(betSetup.fillingCeiling);
    if (autoExternalPressures && nitrogenBETPressurePlan.has_value())
    {
      pressures = nitrogenBETPressurePlan->samplingPressures;
    }
    if (autoReweightingPressureRange && nitrogenBETPressurePlan.has_value())
    {
      reweightingPressureRange =
          std::make_pair(nitrogenBETPressurePlan->lowestPressure, nitrogenBETPressurePlan->highestPressure);
      if (!reader.reweightingNumberOfPressuresSpecified)
      {
        reweightingNumberOfPressures = nitrogenBETPressurePlan->whamReweightingNumberOfPressures;
      }
    }
  }

  numberOfTemperatures = temperatures.size();
  numberOfPressures = pressures.size();
  numberOfReplicas = numberOfTemperatures * numberOfPressures;

  if (temperatures.empty() || pressures.empty())
  {
    throw std::runtime_error(
        "[ReweightedHistogram]: a non-empty temperature ladder and pressure ladder are required\n");
  }
  if (numberOfReplicas < 2)
  {
    throw std::runtime_error(
        "[ReweightedHistogram]: at least two (temperature, pressure) grid points are required\n");
  }

  if (reweightingTemperatures.empty())
  {
    reweightingTemperatures = temperatures;
  }
  if (!autoReweightingPressureRange)
  {
    reweightingPressureRange =
        reader.reweightingPressureRange.value_or(std::make_pair(pressures.front(), pressures.back()));
  }

  validateReweightingRun(numberOfBlocks, temperatures, pressures, reweightingTemperatures, reweightingPressureRange);

  initializeReplicas(std::move(templateSystem));
  numberOfStepsPerCycle = std::max(1uz, systems.front().tmmc.maxMacrostate);
}

ReweightedHistogram::ReweightedHistogram(System templateSystem, std::vector<double> temperatures_,
                                         std::vector<double> pressures_, ReweightedHistogramParameters parameters)
    : random(parameters.randomSeed),
      numberOfProductionCycles(parameters.numberOfProductionCycles),
      numberOfPreInitializationCycles(parameters.numberOfPreInitializationCycles),
      numberOfInitializationCycles(parameters.numberOfInitializationCycles),
      numberOfEquilibrationCycles(parameters.numberOfEquilibrationCycles),
      printEvery(std::max(1uz, parameters.printEvery)),
      optimizeMCMovesEvery(std::max(1uz, parameters.optimizeMCMovesEvery)),
      rescaleWangLandauEvery(std::max(1uz, parameters.rescaleWangLandauEvery)),
      writeBinaryRestartEvery(parameters.writeBinaryRestartEvery),
      numberOfBlocks(parameters.numberOfBlocks),
      parallelTemperingSwapEvery(parameters.parallelTemperingSwapEvery),
      sampleReweightingEvery(std::max(1uz, parameters.sampleReweightingEvery)),
      reweightingTemperatures(std::move(parameters.reweightingTemperatures)),
      reweightingNumberOfPressures(parameters.reweightingNumberOfPressures),
      whamTolerance(parameters.whamTolerance),
      whamMaximumIterations(parameters.whamMaximumIterations),
      betScoutMaximumCycles(parameters.betScoutMaximumCycles),
      computeBET(parameters.computeBET),
      temperatures(std::move(temperatures_)),
      pressures(std::move(pressures_)),
      numberOfTemperatures(temperatures.size()),
      numberOfPressures(pressures.size()),
      numberOfReplicas(numberOfTemperatures * numberOfPressures)
{
  if (temperatures.empty() || pressures.empty())
  {
    throw std::runtime_error(
        "[ReweightedHistogram]: a non-empty temperature ladder and pressure ladder are required\n");
  }

  if (reweightingTemperatures.empty())
  {
    reweightingTemperatures = temperatures;
  }
  reweightingPressureRange =
      parameters.reweightingPressureRange.value_or(std::make_pair(pressures.front(), pressures.back()));
  const bool reweightsBeyondTheSample =
      std::ranges::any_of(reweightingTemperatures, [&](double T) { return T != temperatures.front(); }) ||
      reweightingPressureRange.first != pressures.front() || reweightingPressureRange.second != pressures.front();
  validateReweightingTemplate(templateSystem, temperatures, pressures.size(), reweightsBeyondTheSample);

  validateReweightingRun(numberOfBlocks, temperatures, pressures, reweightingTemperatures, reweightingPressureRange);

  initializeReplicas(std::move(templateSystem));
  numberOfStepsPerCycle = std::max(1uz, systems.front().tmmc.maxMacrostate);
}

void ReweightedHistogram::initializeReplicas(System templateSystem)
{
  systems.reserve(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId + 1 < numberOfReplicas; ++replicaId)
  {
    systems.push_back(templateSystem);
  }
  systems.push_back(std::move(templateSystem));

  randoms.reserve(numberOfReplicas);
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
    {
      const std::size_t replicaId = replicaIndex(temperatureIndex, pressureIndex);
      System& system = systems[replicaId];
      const double T = temperatures[temperatureIndex];
      const double P = pressures[pressureIndex];

      system.temperature = T;
      system.beta = 1.0 / (Units::KB * T);

      if (system.forceField.temperature != T)
      {
        system.forceField.temperature = T;
        system.forceField.preComputeDerivedParameters();
        system.forceField.preComputePotentialShift();
        system.forceField.preComputeTailCorrection();
      }

      pinSystemPengRobinsonPressure(system, P);
      // the log partition functions given as species, at the replica temperature
      for (Component& component : system.components) component.updateLnPartitionFunction(T);

      randoms.emplace_back(random.seed + replicaId + 1);
    }
  }

  // What belongs to the temperature of each replica, and of each temperature that is reweighted to: a flexible
  // molecule is sampled at the activity beta f / q(beta), with q its intramolecular configuration integral (see
  // 'Preparation::idealGasReferences'). The reduced potentials and the configuration swaps need it across
  // temperatures. With it, the conformations of the isolated molecules at the temperature of each replica.
  reweightingLogNormalization = Preparation::prepareReplicas(
      systems, {.ladderTemperatures = temperatures,
                .otherTemperatures = reweightingTemperatures})[systems.front().tmmcComponent()];

  reweightingSamples.resize(numberOfReplicas);
  for (std::vector<Sample>& samples : reweightingSamples)
  {
    samples.reserve(numberOfProductionCycles / sampleReweightingEvery + 1uz);
  }

  stepsPerReplica.assign(numberOfReplicas, 0uz);
  swapAttemptsPerTemperaturePair.assign(numberOfTemperatures > 0uz ? numberOfTemperatures - 1uz : 0uz, 0uz);
  swapAcceptedPerTemperaturePair.assign(numberOfTemperatures > 0uz ? numberOfTemperatures - 1uz : 0uz, 0uz);
  swapAttemptsPerPressurePair.assign(numberOfPressures > 0uz ? numberOfPressures - 1uz : 0uz, 0uz);
  swapAcceptedPerPressurePair.assign(numberOfPressures > 0uz ? numberOfPressures - 1uz : 0uz, 0uz);
}

void ReweightedHistogram::run()
{
  setup();
  runStage(SimulationStage::PreInitialization, numberOfPreInitializationCycles);
  runStage(SimulationStage::Initialization, numberOfInitializationCycles);
  runStage(SimulationStage::Equilibration, numberOfEquilibrationCycles);
  runStage(SimulationStage::Production, numberOfProductionCycles);
  output();
}

void ReweightedHistogram::setup()
{
  if (!systems.empty())
  {
    numberOfStepsPerCycle = std::max(1uz, systems.front().tmmc.maxMacrostate);
  }

  for (System& system : systems)
  {
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);
    system.refusePolarizationWithMobileFramework("The reweighted-histogram driver");
  }

  std::filesystem::create_directories("output");
  std::filesystem::create_directories("wham");
  // On a resume from a binary restart every excess pressure is recomputed before it is sampled: the restart is
  // written right after the replica exchange, which swaps the stored tensors with the configurations and marks
  // them stale, and the marks are not part of the restart. Reset to "fresh", a swapped replica sampled its
  // partner's tensor until its next scheduled recomputation.
  excessPressureIsStale.assign(systems.size(), simulationStage != SimulationStage::Uninitialized ? 1u : 0u);

  // on a binary-restart resume append to the existing output files (and skip re-printing the
  // headers) so each log continues where the interrupted run left off
  const bool resumedFromBinaryRestart = simulationStage != SimulationStage::Uninitialized;
  stream.open("output/output.reweighted_histogram.txt", resumedFromBinaryRestart ? std::ios::app : std::ios::out);
  outputJsonFileName = "output/output.reweighted_histogram.json";

  const System& front = systems.front();
  if (!resumedFromBinaryRestart)
  {
    std::print(stream, "{}", front.writeOutputHeader());
    std::print(stream, "Random seed: {}\n\n", random.seed);
    std::print(stream, "{}\n", HardwareInfo::writeInfo());
    std::print(stream, "{}", Units::printStatus());

    std::print(stream, "Reweighted histogram\n");
    std::print(stream, "========================================================================================================================\n\n");
    std::print(stream, "Number of temperatures:                      {}\n", numberOfTemperatures);
    std::print(stream, "Number of pressures:                         {}\n", numberOfPressures);
    std::print(stream, "Number of replicas / threads:                {}\n", numberOfReplicas);
    std::print(stream, "MC steps per cycle:                          {} (filling ceiling N_max)\n",
               numberOfStepsPerCycle);
    std::string temperatureLadder;
    for (double T : temperatures) temperatureLadder += std::format("{} ", T);
    std::print(stream, "{}",
               wrapText(temperatureLadder + "[K]", "Temperature ladder:                          ",
                        std::string(45, ' ')));
    if (parallelTemperingSwapEvery == 0uz)
    {
      std::print(stream, "Configuration swaps:                         disabled\n");
    }
    else
    {
      std::print(stream, "Configuration-swap sweep every:              {} cycles\n", parallelTemperingSwapEvery);
      std::print(stream, "  (sweeps alternate between the temperature and the pressure direction of the grid)\n");
    }
    std::print(stream, "Reweighting (N, U) sample every:             {} cycles\n", sampleReweightingEvery);
    const bool anyCFCMC = std::ranges::any_of(systems, [](const System& system)
    {
      return std::ranges::any_of(system.components, [](const Component& component)
                                 { return component.hasFractionalMolecule; });
    });
    if (anyCFCMC)
    {
      std::print(stream,
                 "  With CFCMC a sample is recorded whenever every fractional molecule is at λ = 0;\n"
                 "  SampleReweightingEvery applies only when no fractional molecule is present.\n");
    }
    std::print(stream, "Reweighted-isotherm temperatures:           ");
    for (double T : reweightingTemperatures)
    {
      std::print(stream, " {}", T);
    }
    std::print(stream, " [K]\n");
    std::print(stream, "Reweighted-isotherm pressure range:          {:.5e} - {:.5e} [Pa], {} log-spaced points\n\n",
               reweightingPressureRange.first, reweightingPressureRange.second, reweightingNumberOfPressures);
    if (nitrogenBETPressurePlan.has_value())
    {
      writeNitrogenBETPressurePlan(stream, *nitrogenBETPressurePlan);
    }
    if (nitrogenBETFillingCeiling.has_value())
    {
      writeNitrogenBETFillingCeiling(stream, *nitrogenBETFillingCeiling);
    }

  }

#ifdef VERSION
#define QUOTE(str) #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)
  outputJson["version"] = EXPAND_AND_QUOTE(VERSION);
#endif
  outputJson["seed"] = random.seed;
  outputJson["initialization"]["hardwareInfo"] = HardwareInfo::jsonInfo();
  outputJson["initialization"]["units"] = Units::jsonStatus();
  outputJson["initialization"]["temperatures"] = temperatures;
  outputJson["initialization"]["pressures"] = pressures;
  outputJson["initialization"]["numberOfStepsPerCycle"] = numberOfStepsPerCycle;
  outputJson["initialization"]["parallelTemperingSwapEvery"] = parallelTemperingSwapEvery;
  outputJson["initialization"]["sampleReweightingEvery"] = sampleReweightingEvery;
  outputJson["initialization"]["reweightingTemperatures"] = reweightingTemperatures;
  outputJson["initialization"]["reweightingPressureRange"] =
      std::vector<double>{reweightingPressureRange.first, reweightingPressureRange.second};
  outputJson["initialization"]["reweightingNumberOfPressures"] = reweightingNumberOfPressures;
  outputJson["initialization"]["betScoutMaximumCycles"] = betScoutMaximumCycles;
  if (nitrogenBETPressurePlan.has_value())
  {
    outputJson["initialization"]["nitrogenBETHenryCoefficient"] = nitrogenBETPressurePlan->henryCoefficientPerCell;
    outputJson["initialization"]["autoExternalPressures"] = autoExternalPressures;
    outputJson["initialization"]["autoReweightingPressureRange"] = autoReweightingPressureRange;
  }
  if (nitrogenBETFillingCeiling.has_value())
  {
    outputJson["initialization"]["autoMacroStateMaximum"] = autoMacroStateMaximum;
    outputJson["initialization"]["fillingCeiling"] = nitrogenBETFillingCeiling->maxMacrostate;
    outputJson["initialization"]["fillingCeilingMeanOccupancy"] = nitrogenBETFillingCeiling->meanOccupancy;
  }

  // interpolation grids are computed once and shared (copied) between the replicas
  // (with gradients when a gradient-based move is on: the default polynomial scheme has none, and the hybrid
  // and smart Monte Carlo moves ran without the forces of the framework)
  for (System& system : systems) system.upgradeInterpolationSchemeForGradientMoves();
  systems.front().createExternalFieldInterpolationGrid(stream, 0);
  systems.front().createFrameworkInterpolationGrids(stream);
  for (std::size_t replicaId = 1; replicaId < systems.size(); ++replicaId)
  {
    systems[replicaId].externalFieldInterpolationGrid = systems.front().externalFieldInterpolationGrid;
    systems[replicaId].interpolationGrids = systems.front().interpolationGrids;
  }

  if (!resumedFromBinaryRestart && autoExternalPressures && numberOfPressures >= 3uz &&
      nitrogenBETPressurePlan.has_value())
  {
    placeNitrogenBETPressureLadderFromHenryAndSaturation(systems, temperatures, pressures, *nitrogenBETPressurePlan,
                                                         nitrogenBETFillingCeiling, betScoutMaximumCycles, stream,
                                                         outputJson);
    outputJson["initialization"]["pressures"] = pressures;
  }

  // the pressure ladder and the replica grid as they are run: printed after the Fisher re-placement above,
  // which replaces the log-spaced ladder of an automatic setup
  if (!resumedFromBinaryRestart)
  {
    std::string pressureLadder;
    for (double P : pressures) pressureLadder += std::format("{} ", P);
    pressureLadder += "[Pa]";
    if (autoExternalPressures && numberOfPressures >= 3uz && nitrogenBETPressurePlan.has_value())
    {
      pressureLadder += " (log skeleton + Fisher extras from pre-isotherm fit)";
    }
    std::print(stream, "{}",
               wrapText(pressureLadder, "Pressure ladder:                             ", std::string(45, ' ')));
    std::print(stream, "Replica grid: replica (t, p) = t * {} + p, fugacity coefficients from Peng-Robinson\n",
               numberOfPressures);
    std::print(stream, "    replica    temperature [K]    pressure [Pa]\n");
    std::print(stream, "    ------------------------------------------------\n");
    for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
    {
      for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
      {
        const std::size_t replicaId = replicaIndex(temperatureIndex, pressureIndex);
        std::print(stream, "    {:7d}    {:15.4f}    {:13.5e}\n", replicaId, temperatures[temperatureIndex],
                   pressures[pressureIndex]);
      }
    }
    std::print(stream, "\n");
  }

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);

  // per-replica output files: each worker thread writes exclusively to its own stream
  replicaStreams.reserve(numberOfReplicas);
  replicaJsonFileNames.reserve(numberOfReplicas);
  replicaJsons.resize(numberOfReplicas);
  for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
  {
    const System& system = systems[replicaId];
    replicaStreams.emplace_back(std::format("output/output_{}_{}.reweighted_histogram.r{}.txt", system.temperature,
                                            system.input_pressure, replicaId),
                                resumedFromBinaryRestart ? std::ios::app : std::ios::out);
    replicaJsonFileNames.emplace_back(std::format("output/output_{}_{}.reweighted_histogram.r{}.json",
                                                  system.temperature, system.input_pressure, replicaId));

    if (!resumedFromBinaryRestart)
    {
      std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
      std::print(replicaStream, "{}", system.writeOutputHeader());
      std::print(replicaStream, "Reweighted histogram: replica {} of {} (temperature {} [K], pressure {} [Pa])\n",
                 replicaId, numberOfReplicas, system.temperature, system.input_pressure);
      std::print(replicaStream, "Random seed of this replica: {}\n\n", randoms[replicaId].seed);
      std::print(replicaStream, "{}\n", HardwareInfo::writeInfo());
      std::print(replicaStream, "{}", Units::printStatus());
      std::print(replicaStream, "{}", system.writeSystemStatus());
      std::print(replicaStream, "{}", system.forceField.printPseudoAtomStatus());
      std::print(replicaStream, "{}", system.forceField.printForceFieldStatus());
      std::print(replicaStream, "{}", system.writeComponentStatus());
      std::print(replicaStream, "{}", system.writeNumberOfPseudoAtoms());
    }

#ifdef VERSION
    replicaJsons[replicaId]["version"] = EXPAND_AND_QUOTE(VERSION);
#endif
    replicaJsons[replicaId]["seed"] = randoms[replicaId].seed;
    replicaJsons[replicaId]["replicaId"] = replicaId;
    replicaJsons[replicaId]["temperature"] = system.temperature;
    replicaJsons[replicaId]["pressure"] = system.input_pressure;
    // a resumed run keeps the initial conditions it restored, not the mid-run state of the restart
    if (!resumedFromBinaryRestart || !replicaJsons[replicaId]["initialization"].contains("initialConditions"))
    {
      replicaJsons[replicaId]["initialization"]["initialConditions"] = system.jsonSystemStatus();
      replicaJsons[replicaId]["initialization"]["components"] = system.jsonComponentStatus();
    }

    std::ofstream replicaJson(replicaJsonFileNames[replicaId]);
    replicaJson << replicaJsons[replicaId].dump(4);
  }
}

void ReweightedHistogram::performReplicaCycle(std::size_t replicaId, SimulationStage stage, std::size_t currentBlock)
{
  System& system = systems[replicaId];
  RandomNumber& rng = randoms[replicaId];

  // every replica is self-contained; the Gibbs-style moves that need a partner system are not
  // supported by this driver
  std::size_t fractionalMoleculeSystem = 0uz;

  // filling ceiling (MacroStateMaximumNumberOfMolecules), not the current loading: an empty replica
  // must do as many moves as a full one so cycle-based sampling and swaps compare equal work
  for (std::size_t j = 0uz; j != numberOfStepsPerCycle; ++j)
  {
    std::size_t selectedComponent = system.randomComponent(rng);

    switch (stage)
    {
      case SimulationStage::Uninitialized:
        break;
      case SimulationStage::PreInitialization:
        MC_Moves::performRandomMovePreInitialization(rng, system, system, selectedComponent, fractionalMoleculeSystem);
        break;
      case SimulationStage::Initialization:
        MC_Moves::performRandomMoveInitialization(rng, system, system, selectedComponent, fractionalMoleculeSystem);
        break;
      case SimulationStage::Equilibration:
        MC_Moves::performRandomMoveEquilibration(rng, system, system, selectedComponent, fractionalMoleculeSystem);

        // Wang-Landau biasing of the CFCMC lambda moves (all state is owned by this replica)
        system.components[selectedComponent].lambdaGC.WangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample,
            system.lambdaWangLandauIsActive(selectedComponent));
        system.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        system.reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        break;
      case SimulationStage::Production:
        MC_Moves::performRandomMoveProduction(rng, system, system, selectedComponent, fractionalMoleculeSystem,
                                              currentBlock);
        ++stepsPerReplica[replicaId];
        break;
    }

    system.components[selectedComponent].lambdaGC.sampleOccupancy(system.gcFractionalMoleculeIsActive(selectedComponent));
    system.pairSwapLambdaSampleOccupancy();
    system.reactionLambdaSampleOccupancy();
  }
}

void ReweightedHistogram::performSwapSweep(SimulationStage stage, std::size_t numberOfCycles) noexcept
{
  // the replicas keep their (temperature, pressure) state points; the configurations migrate
  // through the grid. Sweeps alternate between the temperature direction and the pressure
  // direction; within a direction the pairing offset alternates between visits, so a configuration
  // can traverse the whole grid.
  bool temperatureDirection = (swapSweeps % 2uz == 0uz);
  if (temperatureDirection && numberOfTemperatures < 2uz) temperatureDirection = false;
  if (!temperatureDirection && numberOfPressures < 2uz) temperatureDirection = true;
  const std::size_t offset = (swapSweeps / 2uz) % 2uz;

  if (temperatureDirection)
  {
    for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
    {
      for (std::size_t temperatureIndex = offset; temperatureIndex + 1 < numberOfTemperatures; temperatureIndex += 2uz)
      {
        ++swapAttempts;
        ++swapAttemptsPerTemperaturePair[temperatureIndex];
        if (MC_Moves::ParallelTemperingSwap(random, systems[replicaIndex(temperatureIndex, pressureIndex)],
                                            systems[replicaIndex(temperatureIndex + 1, pressureIndex)])
                .has_value())
        {
          excessPressureIsStale[replicaIndex(temperatureIndex, pressureIndex)] = 1u;
          excessPressureIsStale[replicaIndex(temperatureIndex + 1, pressureIndex)] = 1u;
          ++swapAccepted;
          ++swapAcceptedPerTemperaturePair[temperatureIndex];
        }
      }
    }
  }
  else
  {
    for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
    {
      for (std::size_t pressureIndex = offset; pressureIndex + 1 < numberOfPressures; pressureIndex += 2uz)
      {
        ++swapAttempts;
        ++swapAttemptsPerPressurePair[pressureIndex];
        if (MC_Moves::ParallelTemperingSwap(random, systems[replicaIndex(temperatureIndex, pressureIndex)],
                                            systems[replicaIndex(temperatureIndex, pressureIndex + 1)])
                .has_value())
        {
          excessPressureIsStale[replicaIndex(temperatureIndex, pressureIndex)] = 1u;
          excessPressureIsStale[replicaIndex(temperatureIndex, pressureIndex + 1)] = 1u;
          ++swapAccepted;
          ++swapAcceptedPerPressurePair[pressureIndex];
        }
      }
    }
  }
  ++swapSweeps;
  ++sweepsThisStage;

  // progress report: all worker threads are parked on the barrier, so this is race-free
  const std::size_t printSweeps = std::max(1uz, printEvery / std::max(1uz, parallelTemperingSwapEvery));
  if (sweepsThisStage % printSweeps == 0uz)
  {
    const std::string_view stageName = (stage == SimulationStage::PreInitialization)  ? "pre-initialization"
                                       : (stage == SimulationStage::Initialization)   ? "initialization"
                                       : (stage == SimulationStage::Equilibration)    ? "equilibration"
                                                                                      : "production";
    const std::size_t cycle = std::min(sweepsThisStage * parallelTemperingSwapEvery, numberOfCycles);
    {
      std::scoped_lock lock(outputMutex);
      std::print(stream, "Reweighted-histogram sweep {} ({}, cycle {} of {}): accepted {}/{} ({:.2f}%)\n",
                 swapSweeps, stageName, cycle, numberOfCycles, swapAccepted, swapAttempts,
                 100.0 * static_cast<double>(swapAccepted) / static_cast<double>(std::max(1uz, swapAttempts)));
      std::flush(stream);
    }

    // running snapshot of the directly-measured adsorption isotherms: all worker threads are
    // parked on the barrier, so reading the per-replica loading averages is race-free
    if (stage == SimulationStage::Production)
    {
      writeIsothermSnapshot();
    }
  }
}

void ReweightedHistogram::runStage(SimulationStage stage, std::size_t numberOfCycles)
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

  // binary restart: stages the restart file was written after are skipped entirely; the stage it
  // was written in resumes from the checkpointed cycle (its per-stage preparation already ran
  // before the checkpoint and must not be repeated)
  if (stage < simulationStage)
  {
    return;
  }
  const std::size_t startCycle = (stage == simulationStage) ? cyclesCompletedThisStage : 0uz;

  simulationStage = stage;

  // serial per-stage preparation
  if (startCycle == 0uz && stage == SimulationStage::Equilibration)
  {
    for (System& system : systems)
    {
      for (Component& component : system.components)
      {
        component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize);
        component.lambdaGC.clear();
      }
      system.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize);
      system.pairSwapLambdaClearBookkeeping();
      system.reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize);
      system.reactionLambdaClearBookkeeping();
    }
  }
  if (startCycle == 0uz && stage == SimulationStage::Production)
  {
    for (System& system : systems)
    {
      system.mc_moves_statistics.clearMoveStatistics();
      system.mc_moves_cputime.clearTimingStatistics();

      for (Component& component : system.components)
      {
        component.mc_moves_statistics.clearMoveStatistics();
        component.mc_moves_cputime.clearTimingStatistics();

        component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize);
        component.lambdaGC.clear();
      }
      system.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize);
      system.pairSwapLambdaClearBookkeeping();
      system.reactionLambdaFinalize();
      system.reactionLambdaClearBookkeeping();
      // Wang-Landau only lowers the bias: left as it ended, exp(-bias) in System::weight() overflows and
      // every weighted average of the production run becomes inf / inf
      system.normalizeLambdaBiasFactors();
    }
    std::fill(stepsPerReplica.begin(), stepsPerReplica.end(), 0uz);
  }
  if (startCycle == 0uz)
  {
    sweepsThisStage = 0uz;
  }

  {
    std::scoped_lock lock(outputMutex);
    const std::string_view stageName = (stage == SimulationStage::PreInitialization)  ? "Pre-initialization"
                                       : (stage == SimulationStage::Initialization)   ? "Initialization"
                                       : (stage == SimulationStage::Equilibration)    ? "Equilibration"
                                                                                      : "Production";
    if (startCycle == 0uz)
    {
      std::print(stream, "\n{} stage: {} cycles on {} replicas/threads\n", stageName, numberOfCycles,
                 numberOfReplicas);
    }
    else
    {
      std::print(stream, "\n{} stage: resumed from binary restart at cycle {} of {} on {} replicas/threads\n",
                 stageName, startCycle, numberOfCycles, numberOfReplicas);
    }
    std::flush(stream);
  }

  // one worker thread per replica; the barrier completion performs the swap sweeps and writes the
  // periodic binary restart file (all worker threads are parked, so the state is consistent)
  auto onAllArrived = [this, stage, numberOfCycles]() noexcept
  {
    const std::size_t completedCycles = checkpointCycle.load(std::memory_order_relaxed);
    if (parallelTemperingSwapEvery != 0uz && completedCycles % parallelTemperingSwapEvery == 0uz)
    {
      performSwapSweep(stage, numberOfCycles);
    }
    const bool stopRequested = GracefulShutdown::requested();
    const bool binaryRestartDue =
        writeBinaryRestartEvery != 0uz &&
        (completedCycles % writeBinaryRestartEvery == 0uz || completedCycles == numberOfCycles);
    if (binaryRestartDue || stopRequested)
    {
      writeBinaryRestartFile(completedCycles);
    }
    // all worker threads are parked on this barrier, so the checkpoint just written is a
    // consistent snapshot: safe to exit here on a shutdown signal
    if (stopRequested)
    {
      // std::exit skips stack unwinding: flush the text output streams explicitly
      std::flush(stream);
      for (std::ofstream& replicaStream : replicaStreams) std::flush(replicaStream);
      GracefulShutdown::exitAfterCheckpoint();
    }
  };
  std::barrier synchronizationPoint(static_cast<std::ptrdiff_t>(numberOfReplicas), onAllArrived);

  const std::size_t stageCycleOffset = absoluteCycleOffset;

  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfReplicas);
    for (std::size_t replicaId = 0; replicaId < numberOfReplicas; ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId, stage, numberOfCycles, stageCycleOffset, startCycle, &synchronizationPoint]()
          {
            System& system = systems[replicaId];

            // each thread computes the total energies of its own replica
            if (stage == SimulationStage::PreInitialization || stage == SimulationStage::Initialization)
            {
              system.precomputeTotalRigidEnergy();
            }
            system.runningEnergies = system.computeTotalEnergies();

            BlockErrorEstimation estimation(numberOfBlocks, std::max(1uz, numberOfProductionCycles));

            for (std::size_t cycle = startCycle; cycle != numberOfCycles; ++cycle)
            {
              if (stage == SimulationStage::Production)
              {
                estimation.setCurrentSample(cycle);
              }

              performReplicaCycle(replicaId, stage, estimation.currentBin);

              // time-evolution properties (number of molecules, volume): sampled over all stages,
              // indexed by the absolute cycle number; the writers gate on their own 'writeEvery'
              const std::size_t absoluteCycle = stageCycleOffset + cycle;
              system.samplePropertiesEvolution(absoluteCycle);
              if (system.propertyNumberOfMoleculesEvolution.has_value())
              {
                system.propertyNumberOfMoleculesEvolution->writeOutput(replicaId, absoluteCycle);
              }
              if (system.propertyVolumeEvolution.has_value())
              {
                system.propertyVolumeEvolution->writeOutput(replicaId, absoluteCycle);
              }

              if (stage == SimulationStage::Production)
              {
                // the virial of this configuration before it is sampled: sampleProperties adds the current excess
                // pressure to the pressure averages, and computed after it the first cycle sampled no virial and the
                // others an earlier one. It is recomputed every tenth cycle and is a lagged sample of the same chain
                // in between; after an accepted swap it is not (the replica holds a new configuration and the tensor
                // of the one it gave away), so the first cycle after the swap recomputes it as well.
                std::optional<std::pair<EnergyStatus, double3x3>> freshPressure{};
                std::chrono::steady_clock::duration freshPressureTime{};
                if (cycle % 10uz == 0uz || cycle % printEvery == 0uz || excessPressureIsStale[replicaId] != 0u)
                {
                  excessPressureIsStale[replicaId] = 0u;
                  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
                  freshPressure = system.computeMolecularPressure();
                  system.currentEnergyStatus = freshPressure->first;
                  system.currentExcessPressureTensor = freshPressure->second / system.simulationBox.volume;
                  freshPressureTime = std::chrono::steady_clock::now() - start;
                }
                // (the pressure is sampled in the cycles its tensor was computed for)
                system.sampleProperties(replicaId, estimation.currentBin, cycle, false, freshPressure.has_value());
                // the force-based RDF needs the full gradients, which the moves leave out: sampled as in MonteCarlo
                if (system.forceBasedRDFSampleDue(cycle)) system.sampleForceBasedRDFWithFullGradients(cycle, estimation.currentBin);

                // the raw (N, U) reweighting sample of this replica: the molecule count of the
                // single adsorbate component and the potential energy, tagged with the block.
                //
                // A CFCMC fractional molecule only counts once it is fully decoupled. Its partial
                // interactions are already inside runningEnergies while it is absent from
                // numberOfIntegerMoleculesPerComponent, so at any intermediate lambda the pair
                // would carry the energy of N molecules plus a fraction against a count of N, and
                // WHAM would bin those shifted energies as if they belonged to the N-molecule
                // state. At lambda = 0 the configuration is an ordinary N-molecule one.
                // Lambda is decoupled on only a small share of cycles, so when a fractional
                // molecule exists every decoupled configuration is taken instead of every nth
                // cycle. Applying both gates at once starved the analysis: on FER most replicas
                // ended a 10000-cycle production run with no samples at all.
                const bool anyFractionalMolecule =
                    std::ranges::any_of(system.components, [](const Component& component)
                                        { return component.hasFractionalMolecule; });
                // Every lambda coordinate of the component has to be at bin 0: besides the GC fractional molecule
                // a component may carry pair-swap, group-swap, or Gibbs fractional molecules with their own
                // lambda histograms (an unused histogram stays at bin 0).
                const bool fractionalMoleculesAreDecoupled = std::ranges::all_of(
                    system.components,
                    [](const Component& component)
                    {
                      return !component.hasFractionalMolecule ||
                             (component.lambdaGC.currentBin == 0uz && component.lambdaPairSwap.currentBin == 0uz &&
                              component.lambdaPairSwapCB.currentBin == 0uz &&
                              component.lambdaGroupSwap.currentBin == 0uz &&
                              component.lambdaGroupSwapCB.currentBin == 0uz &&
                              component.lambdaGibbs.currentBin == 0uz);
                    });
                if (anyFractionalMolecule ? fractionalMoleculesAreDecoupled
                                          : cycle % sampleReweightingEvery == 0uz)
                {
                  // at lambda = 0 the fractional molecule still carries its intramolecular energy
                  const double sampleEnergy =
                      system.runningEnergies.potentialEnergy() -
                      (anyFractionalMolecule ? fractionalMoleculesIntraEnergy(system) : 0.0);
                  reweightingSamples[replicaId].push_back(
                      Sample{.energy = sampleEnergy,
                             .numberOfMolecules =
                                 static_cast<std::uint32_t>(
                                     system.numberOfIntegerMoleculesPerComponent[system.tmmcComponent()]),
                             .block = static_cast<std::uint32_t>(estimation.currentBin)});
                }

                // analysis-property files (RDFs, density grid, histograms, molecule properties);
                // the writers gate on their own 'writeEvery'
                writeReplicaAnalysisOutputs(system, replicaId, cycle);

                // energy/pressure averages for the per-replica final report
                if (cycle % 10uz == 0uz || cycle % printEvery == 0uz)
                {
                  const std::pair<EnergyStatus, double3x3> molecularPressure = *freshPressure;

                  system.mc_moves_cputime.energyPressureComputation += freshPressureTime;
                  system.averageEnergies.addSample(estimation.currentBin, molecularPressure.first, system.weight());
                }
              }

              // the step sizes are tuned while the system equilibrates only: an adaptive kernel during production
              // makes the chain depend on its history and breaks detailed balance
              if (stage != SimulationStage::Production && cycle % optimizeMCMovesEvery == 0uz)
              {
                system.optimizeMCMoves();
              }

              // Wang-Landau biasing-factor adjustment (all state is owned by this replica)
              if (stage == SimulationStage::Equilibration && (cycle + 1uz) % rescaleWangLandauEvery == 0uz)
              {
                for (Component& component : system.components)
                {
                  component.lambdaGC.WangLandauIteration(
                      PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors);
                }
                system.pairSwapLambdaWangLandauIteration(
                    PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors);
                system.reactionLambdaWangLandauIteration(
                    PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors);
              }

              if (cycle % printEvery == 0uz)
              {
                // each thread writes exclusively to its own replica stream: no locking needed
                system.loadings = LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent,
                                              system.simulationBox);

                std::ostream replicaStream(replicaStreams[replicaId].rdbuf());
                switch (stage)
                {
                  case SimulationStage::PreInitialization:
                    std::print(replicaStream, "{}", system.writePreInitializationStatusReport(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Initialization:
                    std::print(replicaStream, "{}", system.writeInitializationStatusReport(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Equilibration:
                    std::print(replicaStream, "{}", system.writeEquilibrationStatusReportMC(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Production:
                  {
                    std::string status_line = std::format("Current cycle: {} out of {}\n", cycle, numberOfCycles);
                    std::print(
                        replicaStream, "{}",
                        system.writeProductionStatusReportMC(status_line, std::make_pair(cycle, numberOfCycles)));
                    break;
                  }
                  default:
                    break;
                }
                std::flush(replicaStream);
              }

              // the only synchronization point between the threads: the configuration-swap sweep
              // and the periodic binary-restart checkpoint, both performed by the barrier completion, and a
              // regular check for a shutdown signal (a run without swaps and without binary restarts never
              // reached the barrier and ignored the signal until the stage was over)
              const bool swapDue =
                  parallelTemperingSwapEvery != 0uz && (cycle + 1uz) % parallelTemperingSwapEvery == 0uz;
              const bool binaryRestartDue =
                  writeBinaryRestartEvery != 0uz &&
                  ((cycle + 1uz) % writeBinaryRestartEvery == 0uz || cycle + 1uz == numberOfCycles);
              const bool shutdownCheckDue = (cycle + 1uz) % GracefulShutdown::checkEvery == 0uz;
              if (swapDue || binaryRestartDue || shutdownCheckDue)
              {
                if (replicaId == 0uz)
                {
                  checkpointCycle.store(cycle + 1uz, std::memory_order_relaxed);
                }
                synchronizationPoint.arrive_and_wait();
              }
            }
          });
    }
    // jthreads join on scope exit
  }

  absoluteCycleOffset += numberOfCycles;

  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  switch (stage)
  {
    case SimulationStage::PreInitialization:
      totalPreInitializationSimulationTime += (t2 - t1);
      break;
    case SimulationStage::Initialization:
      totalInitializationSimulationTime += (t2 - t1);
      break;
    case SimulationStage::Equilibration:
      totalEquilibrationSimulationTime += (t2 - t1);
      break;
    case SimulationStage::Production:
      totalProductionSimulationTime += (t2 - t1);
      break;
    default:
      break;
  }
  totalSimulationTime += (t2 - t1);
}

void ReweightedHistogram::output()
{
  std::size_t numberOfSteps = std::accumulate(stepsPerReplica.begin(), stepsPerReplica.end(), 0uz);

  MCMoveCpuTime total;
  MCMoveStatistics countTotal;
  for (const System& system : systems)
  {
    total += system.mc_moves_cputime;
    countTotal += system.mc_moves_statistics;
    for (const Component& component : system.components)
    {
      countTotal += component.mc_moves_statistics;
    }
  }

  std::print(stream, "\n");
  std::print(stream, "========================================================================================================================\n");
  std::print(stream, "                             Simulation finished!\n");
  std::print(stream, "========================================================================================================================\n");
  std::print(stream, "\n");

  // energy drift check of every replica (energies recomputed in parallel, one thread per replica);
  // the final state used by the per-replica reports is refreshed in the same pass
  std::vector<RunningEnergy> recomputed(systems.size());
  {
    std::vector<std::jthread> threads;
    threads.reserve(systems.size());
    for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
    {
      threads.emplace_back(
          [this, replicaId, &recomputed]()
          {
            System& system = systems[replicaId];
            recomputed[replicaId] = system.computeTotalEnergies();

            std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
            system.currentEnergyStatus = molecularPressure.first;
            system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
            system.loadings = LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent,
                                          system.simulationBox);
          });
    }
  }

  writeReplicaFinalReports(recomputed);

  // final directly-measured adsorption isotherms (also covers runs with the swaps disabled)
  writeIsothermSnapshot();

  std::print(stream, "Energy drift per replica\n");
  std::print(stream, "========================================================================================================================\n\n");
  for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
  {
    const RunningEnergy drift = systems[replicaId].runningEnergies - recomputed[replicaId];
    std::print(stream, "    replica {:4d} (temperature {:10.4f} [K], pressure {:13.5e} [Pa]): drift {: .6e} [K]\n",
               replicaId, systems[replicaId].temperature, systems[replicaId].input_pressure,
               Units::EnergyToKelvin * drift.potentialEnergy());
  }
  std::print(stream, "\n\n");

  std::print(stream, "Production run counting of the MC moves summed over replicas and components\n");
  std::print(stream, "========================================================================================================================\n\n");
  std::print(stream, "{}", countTotal.writeMCMoveStatistics(numberOfSteps));
  std::print(stream, "\n\n");

  std::print(stream, "Replica-exchange swap statistics\n");
  std::print(stream, "========================================================================================================================\n\n");
  std::print(stream, "    sweeps:    {}\n", swapSweeps);
  std::print(stream, "    attempts:  {}\n", swapAttempts);
  std::print(stream, "    accepted:  {} ({:.4f} %)\n\n", swapAccepted,
             100.0 * static_cast<double>(swapAccepted) / static_cast<double>(std::max(1uz, swapAttempts)));

  // low acceptance for a particular pair marks a bottleneck in the grid (configurations cannot
  // migrate past it); consider a denser ladder around such a pair
  if (numberOfTemperatures > 1uz)
  {
    std::print(stream, "    temperature direction (aggregated over the pressures)\n");
    std::print(stream, "    pair (indices)     temperature [K]           attempts    accepted    acceptance\n");
    std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
    for (std::size_t temperatureIndex = 0; temperatureIndex + 1 < numberOfTemperatures; ++temperatureIndex)
    {
      std::print(stream, "    {:4d} - {:<4d}   {:10.4f} - {:<10.4f}   {:9d}   {:9d}    {:8.4f} %\n", temperatureIndex,
                 temperatureIndex + 1, temperatures[temperatureIndex], temperatures[temperatureIndex + 1],
                 swapAttemptsPerTemperaturePair[temperatureIndex], swapAcceptedPerTemperaturePair[temperatureIndex],
                 100.0 * static_cast<double>(swapAcceptedPerTemperaturePair[temperatureIndex]) /
                     static_cast<double>(std::max(1uz, swapAttemptsPerTemperaturePair[temperatureIndex])));
    }
    std::print(stream, "\n");
  }
  if (numberOfPressures > 1uz)
  {
    std::print(stream, "    pressure direction (aggregated over the temperatures)\n");
    std::print(stream, "    pair (indices)     pressure [Pa]                 attempts    accepted    acceptance\n");
    std::print(stream, "    --------------------------------------------------------------------------------------------------------------------\n");
    for (std::size_t pressureIndex = 0; pressureIndex + 1 < numberOfPressures; ++pressureIndex)
    {
      std::print(stream, "    {:4d} - {:<4d}   {:12.5e} - {:<12.5e}   {:9d}   {:9d}    {:8.4f} %\n", pressureIndex,
                 pressureIndex + 1, pressures[pressureIndex], pressures[pressureIndex + 1],
                 swapAttemptsPerPressurePair[pressureIndex], swapAcceptedPerPressurePair[pressureIndex],
                 100.0 * static_cast<double>(swapAcceptedPerPressurePair[pressureIndex]) /
                     static_cast<double>(std::max(1uz, swapAttemptsPerPressurePair[pressureIndex])));
    }
    std::print(stream, "\n");
  }
  std::print(stream, "\n");

  // the multiple-histogram reweighting analysis: combines the pooled (N, U) samples of all
  // replicas and writes the reweighted isotherms and free energies
  performReweightingAnalysis();

  std::print(stream, "Production run CPU timings of the MC moves summed over replicas and components\n");
  std::print(stream, "========================================================================================================================\n\n");
  std::print(stream, "{}", total.writeMCMoveCPUTimeStatistics(totalProductionSimulationTime));
  std::print(stream, "Pre-initialization simulation time: {:14f} [s]\n", totalPreInitializationSimulationTime.count());
  std::print(stream, "Initalization simulation time:  {:14f} [s]\n", totalInitializationSimulationTime.count());
  std::print(stream, "Equilibration simulation time:  {:14f} [s]\n", totalEquilibrationSimulationTime.count());
  std::print(stream, "Production simulation time:     {:14f} [s]\n", totalProductionSimulationTime.count());
  std::print(stream, "Reweighting analysis time:      {:14f} [s]\n", totalReweightingAnalysisTime.count());
  std::print(stream, "Total simulation time:          {:14f} [s]\n",
             (totalSimulationTime + totalReweightingAnalysisTime).count());
  std::print(stream, "\n\n");
  std::flush(stream);

  outputJson["output"]["numberOfSteps"] = numberOfSteps;
  outputJson["output"]["replicaExchange"]["sweeps"] = swapSweeps;
  outputJson["output"]["replicaExchange"]["attempts"] = swapAttempts;
  outputJson["output"]["replicaExchange"]["accepted"] = swapAccepted;
  outputJson["output"]["replicaExchange"]["attemptsPerTemperaturePair"] = swapAttemptsPerTemperaturePair;
  outputJson["output"]["replicaExchange"]["acceptedPerTemperaturePair"] = swapAcceptedPerTemperaturePair;
  outputJson["output"]["replicaExchange"]["attemptsPerPressurePair"] = swapAttemptsPerPressurePair;
  outputJson["output"]["replicaExchange"]["acceptedPerPressurePair"] = swapAcceptedPerPressurePair;
  outputJson["output"]["cpuTimings"]["preInitialization"] = totalPreInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["initialization"] = totalInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["equilibration"] = totalEquilibrationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["production"] = totalProductionSimulationTime.count();
  outputJson["output"]["cpuTimings"]["reweightingAnalysis"] = totalReweightingAnalysisTime.count();
  outputJson["output"]["cpuTimings"]["total"] = (totalSimulationTime + totalReweightingAnalysisTime).count();

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);
}

void ReweightedHistogram::writeIsothermSnapshot() const
{
  // one isotherm file per temperature: every replica along the pressure ladder contributes one
  // directly-measured point, assembled from the per-replica block averages
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    std::ofstream isotherm(std::format("wham/isotherm_{}.reweighted_histogram.txt",
                                       temperatures[temperatureIndex]),
                           std::ios::trunc);

    std::print(isotherm, "# Reweighted histogram: directly-measured adsorption isotherm at {} [K]\n",
               temperatures[temperatureIndex]);
    std::print(isotherm, "# one block per component (gnuplot 'index'), rows ordered by increasing pressure\n");
    std::print(isotherm, "# errors are the block confidence-interval errors (zero until three blocks are filled)\n");
    std::print(isotherm, "# column 1: fugacity [Pa]\n");
    std::print(isotherm, "# column 2, 3: absolute loading, error [molecules / simulation supercell]\n");
    std::print(isotherm, "# column 4, 5: absolute loading, error [molecules / unit cell]\n");
    std::print(isotherm, "# column 6, 7: absolute loading, error [mol/kg-framework]\n");
    std::print(isotherm, "# column 8, 9: absolute loading, error [mg/g-framework]\n");
    std::print(isotherm, "# column 10: pressure [Pa]\n");

    const System& front = systems[replicaIndex(temperatureIndex, 0)];
    for (std::size_t componentId = 0; componentId < front.components.size(); ++componentId)
    {
      std::print(isotherm, "\n\n# component {} ({})\n", componentId, front.components[componentId].name);

      for (std::size_t pressureIndex = 0; pressureIndex < numberOfPressures; ++pressureIndex)
      {
        const System& system = systems[replicaIndex(temperatureIndex, pressureIndex)];
        const Component& component = system.components[componentId];

        const double fugacity =
            component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.input_pressure;

        const auto [loadingAverage, loadingError] = system.averageLoadings.average();
        const double moleculesInSupercell = loadingAverage.numberOfMolecules[componentId];
        const double moleculesInSupercellError = loadingError.numberOfMolecules[componentId];

        // without a framework the per-unit-cell and per-framework-mass units are not defined
        double toMoleculesPerUnitCell = 1.0;
        double toMolePerKg = 0.0;
        double toMgPerG = 0.0;
        if (system.framework.has_value())
        {
          const int3 numberOfUnitCells = system.framework->numberOfUnitCells;
          toMoleculesPerUnitCell =
              1.0 / static_cast<double>(numberOfUnitCells.x * numberOfUnitCells.y * numberOfUnitCells.z);
          // (a framework without mass has no loading per kg, see System::frameworkMass)
          if (const std::optional<double> frameworkMass = system.frameworkMass(); frameworkMass.has_value())
          {
            toMolePerKg = 1000.0 / frameworkMass.value();
            toMgPerG = 1000.0 * component.totalMass / frameworkMass.value();
          }
        }

        std::print(isotherm,
                   "{: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e}\n",
                   fugacity, moleculesInSupercell, moleculesInSupercellError,
                   toMoleculesPerUnitCell * moleculesInSupercell, toMoleculesPerUnitCell * moleculesInSupercellError,
                   toMolePerKg * moleculesInSupercell, toMolePerKg * moleculesInSupercellError,
                   toMgPerG * moleculesInSupercell, toMgPerG * moleculesInSupercellError, system.input_pressure);
      }
    }
  }
}

void ReweightedHistogram::performReweightingAnalysis()
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

  std::print(stream, "Multiple-histogram reweighting (WHAM) analysis\n");
  std::print(stream, "========================================================================================================================\n\n");

  std::size_t totalNumberOfSamples = 0uz;
  for (const std::vector<Sample>& samples : reweightingSamples)
  {
    totalNumberOfSamples += samples.size();
  }
  if (totalNumberOfSamples == 0uz)
  {
    whamConverged = false;  // nothing was solved: an empty analysis is not a converged one
    whamResidual = 0.0;
    whamIterations = 0uz;
    whamUnconvergedBlocks = 0uz;
    whamUnconvergedBlockDetails.clear();
    std::print(stream, "    no (N, U) samples were collected (no production cycles); analysis skipped\n\n\n");
    return;
  }

  reweightedIsotherms.clear();

  // per-state (replica) thermodynamic fields, all in internal units: the reduced potential of a
  // sample (N, U) at state i is u_i = beta_i U - N ln(beta_i f_i)
  std::vector<double> stateBeta(numberOfReplicas);
  std::vector<double> stateLogBetaFugacity(numberOfReplicas);
  for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
  {
    const System& system = systems[stateId];
    const Component& component = system.components[system.tmmcComponent()];
    const double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    stateBeta[stateId] = system.beta;
    stateLogBetaFugacity[stateId] = std::log(system.beta * fugacity);
  }
  // A flexible adsorbate on a temperature ladder is sampled at the activity beta f e^{-L(beta)}, with L the ln of
  // what its insertions divide beta f by ('Component::logInsertionNormalization'): the sample (N, U) has the
  // weight -N L(beta) next to N ln(beta f), at a replica and at a state point that is reweighted to (there L is
  // the ln of the intramolecular configuration integral). At the temperature of a replica it is the normalization
  // of that replica, which holds the weight the replica was run with. Zero where no reference is needed.
  std::map<double, double> logNormalizationOfBeta;
  for (std::size_t i = 0; i < reweightingTemperatures.size() && i < reweightingLogNormalization.size(); ++i)
  {
    logNormalizationOfBeta[1.0 / (Units::KB * reweightingTemperatures[i])] = reweightingLogNormalization[i];
  }
  for (const System& system : systems)
  {
    logNormalizationOfBeta[system.beta] =
        system.components[system.tmmcComponent()].logInsertionNormalization.value_or(0.0);
  }
  auto logNormalizationAt = [&](double beta)
  {
    const auto it = logNormalizationOfBeta.find(beta);
    return it == logNormalizationOfBeta.end() ? 0.0 : it->second;
  };

  // bin the pooled samples over (N, energy bin); the mean energy of the samples inside a bin is
  // used in the reduced potentials, which removes the leading-order binning error
  double energyMinimum = 1e300;
  double energyMaximum = -1e300;
  for (const std::vector<Sample>& samples : reweightingSamples)
  {
    for (const Sample& sample : samples)
    {
      energyMinimum = std::min(energyMinimum, sample.energy);
      energyMaximum = std::max(energyMaximum, sample.energy);
    }
  }
  constexpr std::size_t numberOfEnergyBins = 512uz;
  const double deltaEnergy =
      std::max((energyMaximum - energyMinimum) / static_cast<double>(numberOfEnergyBins), 1e-10);

  std::map<std::pair<std::uint32_t, std::uint32_t>, std::size_t> binIds;
  std::vector<std::uint32_t> binMoleculeCount;      // N of the bin
  std::vector<double> binEnergySum;                 // sum of the sampled energies in the bin
  std::vector<std::vector<double>> binBlockCounts;  // per-bin sample count per block
  // per bin: whether state i has a sample there in block j (at i * numberOfBlocks + j), for the overlap test
  std::vector<std::vector<bool>> binSampledBy;
  std::vector<std::vector<double>> stateBlockCounts(numberOfReplicas,
                                                    std::vector<double>(numberOfBlocks, 0.0));

  for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
  {
    for (const Sample& sample : reweightingSamples[stateId])
    {
      const std::uint32_t energyBin = static_cast<std::uint32_t>(std::min(
          numberOfEnergyBins - 1uz,
          static_cast<std::size_t>(std::max(0.0, (sample.energy - energyMinimum) / deltaEnergy))));
      const auto [it, inserted] =
          binIds.try_emplace(std::make_pair(sample.numberOfMolecules, energyBin), binIds.size());
      if (inserted)
      {
        binMoleculeCount.push_back(sample.numberOfMolecules);
        binEnergySum.push_back(0.0);
        binBlockCounts.emplace_back(numberOfBlocks, 0.0);
        binSampledBy.emplace_back(numberOfReplicas * numberOfBlocks, false);
      }
      const std::size_t binId = it->second;
      binEnergySum[binId] += sample.energy;
      binBlockCounts[binId][sample.block] += 1.0;
      binSampledBy[binId][stateId * numberOfBlocks + sample.block] = true;
      stateBlockCounts[stateId][sample.block] += 1.0;
    }
  }
  const std::size_t numberOfBins = binIds.size();

  std::vector<double> binTotalCount(numberOfBins, 0.0);
  std::vector<double> binMeanEnergy(numberOfBins);
  for (std::size_t binId = 0; binId < numberOfBins; ++binId)
  {
    binTotalCount[binId] = std::accumulate(binBlockCounts[binId].begin(), binBlockCounts[binId].end(), 0.0);
    binMeanEnergy[binId] = binEnergySum[binId] / binTotalCount[binId];
  }

  // precomputed reduced potentials u_i(b) = beta_i U_b - N_b ln(beta_i f_i)
  std::vector<double> reducedPotential(numberOfBins * numberOfReplicas);
  for (std::size_t binId = 0; binId < numberOfBins; ++binId)
  {
    for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
    {
      reducedPotential[binId * numberOfReplicas + stateId] =
          stateBeta[stateId] * binMeanEnergy[binId] -
          static_cast<double>(binMoleculeCount[binId]) *
              (stateLogBetaFugacity[stateId] - logNormalizationAt(stateBeta[stateId]));
    }
  }

  // the code base is compiled with -ffast-math, so infinities must not occur: empty bins/states
  // carry this finite 'log of zero' sentinel instead and drop out of the sums
    constexpr double logZero = -1e300;
  constexpr double logZeroThreshold = -1e299;

  // Whether the sampled states form one group linked through (N, U) bins in which both of them have samples
  // (of one block, or of all of them). The free-energy difference of two states is measured only there. The share
  // n_k exp(g_k - u_k(b)) of a state in the mixture of a bin is no such evidence: WHAM knows the samples of a bin
  // by their sum only, and two states with disjoint histograms still have a unique fixed point. A replica that
  // stays in the empty box (0, 0) and one that stays condensed at (N_1, U_1), on either side of a nucleation
  // barrier neither crosses, converge on g_1 = -(N_1 / 2) ln(f_1 / f_0), each credited with a share of the
  // samples of the other that is far above any threshold: a condensation step at sqrt(f_0 f_1) was reported as
  // converged, whatever the free energy of the condensed phase.
  //
  // Two states at one temperature are linked by a molecule count both of them sampled, whatever the energies: their
  // reduced potentials differ by -N ln(f_i / f_j) only, the energy drops out of the WHAM equations of the pair, and
  // both draw U at that N from one conditional distribution. Asked for a shared energy bin as well (512 of them,
  // against the few samples of a count in the tail where two pressures overlap), the neighbours of an isothermal
  // pressure ladder were reported as disconnected, and with them the blocks, which have 1/B of the samples each.
  auto statesAreConnected = [&](const std::vector<double>& logStateCounts,
                                const std::vector<double>& logDensityOfStates,
                                std::optional<std::size_t> block) -> bool
  {
    std::vector<std::size_t> parent(numberOfReplicas);
    std::iota(parent.begin(), parent.end(), 0uz);
    auto find = [&](std::size_t k)
    {
      while (parent[k] != k)
      {
        parent[k] = parent[parent[k]];
        k = parent[k];
      }
      return k;
    };
    auto sampled = [&](std::size_t binId, std::size_t stateId)
    {
      const std::vector<bool>& flags = binSampledBy[binId];
      if (block.has_value()) return static_cast<bool>(flags[stateId * numberOfBlocks + block.value()]);
      for (std::size_t j = 0; j < numberOfBlocks; ++j)
      {
        if (flags[stateId * numberOfBlocks + j]) return true;
      }
      return false;
    };
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      if (logDensityOfStates[binId] <= logZeroThreshold) continue;
      std::optional<std::size_t> first;
      for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
      {
        if (logStateCounts[stateId] <= logZeroThreshold || !sampled(binId, stateId)) continue;
        if (!first.has_value())
        {
          first = stateId;
          continue;
        }
        parent[find(stateId)] = find(first.value());
      }
    }
    std::map<std::pair<std::uint32_t, double>, std::size_t> firstOfMoleculeCount;  // by (N, beta)
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      if (logDensityOfStates[binId] <= logZeroThreshold) continue;
      for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
      {
        if (logStateCounts[stateId] <= logZeroThreshold || !sampled(binId, stateId)) continue;
        const auto [it, inserted] =
            firstOfMoleculeCount.try_emplace(std::make_pair(binMoleculeCount[binId], stateBeta[stateId]), stateId);
        if (!inserted) parent[find(stateId)] = find(it->second);
      }
    }
    std::optional<std::size_t> root;
    for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
    {
      if (logStateCounts[stateId] <= logZeroThreshold) continue;
      if (!root.has_value())
      {
        root = find(stateId);
      }
      else if (find(stateId) != root.value())
      {
        return false;
      }
    }
    return true;
  };

  // Solves the WHAM self-consistent equations in log space by direct iteration:
  //     ln Omega_b = ln M_b - ln sum_i n_i exp(g_i - u_i(b))
  //     g_i        = -ln sum_b Omega_b exp(-u_i(b))
  // gauge-fixed on the first state with samples (a replica without samples in a block does not enter Omega, and
  // its g is arbitrary). Returns the number of iterations performed, the final residual, and whether the
  // sampled states overlap as one group.
  struct WhamConvergence
  {
    std::size_t iterations;
    double residual;
    bool connected;
  };
  auto solveWham = [&](const std::vector<double>& logStateCounts, const std::vector<double>& logBinCounts,
                       std::vector<double>& freeEnergies, std::vector<double>& logDensityOfStates,
                       std::size_t maximumNumberOfIterations, std::optional<std::size_t> block) -> WhamConvergence
  {
    const double tolerance = whamTolerance;
    std::vector<double> updatedFreeEnergies(numberOfReplicas);
    std::vector<bool> takesPart(numberOfReplicas, true);
    logDensityOfStates.assign(numberOfBins, logZero);
    std::size_t referenceState = 0uz;
    while (referenceState < numberOfReplicas && logStateCounts[referenceState] <= logZeroThreshold) ++referenceState;
    // no samples at all: nothing was solved, which is not a converged solution of connected states
    if (referenceState == numberOfReplicas || numberOfBins == 0uz) return {0uz, 1e30, false};

    // Omega_b from the current g_i
    auto updateDensity = [&]()
    {
      for (std::size_t binId = 0; binId < numberOfBins; ++binId)
      {
        if (logBinCounts[binId] <= logZeroThreshold)
        {
          logDensityOfStates[binId] = logZero;
          continue;
        }
        double largest = logZero;
        for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
        {
          if (logStateCounts[stateId] <= logZeroThreshold) continue;
          largest = std::max(largest, logStateCounts[stateId] + freeEnergies[stateId] -
                                          reducedPotential[binId * numberOfReplicas + stateId]);
        }
        // no sampled state contributes to the bin (it can not hold samples then; guarded all the same, the
        // logarithm of an empty sum being undefined under -ffast-math)
        if (largest <= logZeroThreshold)
        {
          logDensityOfStates[binId] = logZero;
          continue;
        }
        double sum = 0.0;
        for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
        {
          if (logStateCounts[stateId] <= logZeroThreshold) continue;
          sum += std::exp(logStateCounts[stateId] + freeEnergies[stateId] -
                          reducedPotential[binId * numberOfReplicas + stateId] - largest);
        }
        if (!(sum > 0.0))
        {
          logDensityOfStates[binId] = logZero;
          continue;
        }
        logDensityOfStates[binId] = logBinCounts[binId] - largest - std::log(sum);
      }
    };

    std::size_t iteration = 0uz;
    double residual = 1e30;
    for (; iteration < maximumNumberOfIterations && residual > tolerance; ++iteration)
    {
      updateDensity();

      for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
      {
        double largest = logZero;
        for (std::size_t binId = 0; binId < numberOfBins; ++binId)
        {
          if (logDensityOfStates[binId] <= logZeroThreshold) continue;
          largest = std::max(largest,
                             logDensityOfStates[binId] - reducedPotential[binId * numberOfReplicas + stateId]);
        }
        takesPart[stateId] = largest > logZeroThreshold;
        if (!takesPart[stateId])
        {
          updatedFreeEnergies[stateId] = 0.0;
          continue;
        }
        double sum = 0.0;
        for (std::size_t binId = 0; binId < numberOfBins; ++binId)
        {
          if (logDensityOfStates[binId] <= logZeroThreshold) continue;
          sum += std::exp(logDensityOfStates[binId] - reducedPotential[binId * numberOfReplicas + stateId] - largest);
        }
        updatedFreeEnergies[stateId] = -(largest + std::log(sum));
      }

      const double gauge = updatedFreeEnergies[referenceState];
      residual = 0.0;
      for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
      {
        if (takesPart[stateId]) updatedFreeEnergies[stateId] -= gauge;  // the others stay at zero
        if (logStateCounts[stateId] <= logZeroThreshold) continue;  // an unsampled state does not take part
        residual = std::max(residual, std::abs(updatedFreeEnergies[stateId] - freeEnergies[stateId]));
      }
      freeEnergies = updatedFreeEnergies;
    }
    // the density of states that belongs to the returned free energies (inside the loop it is one iteration and
    // one gauge shift behind them)
    updateDensity();
    return {iteration, residual, statesAreConnected(logStateCounts, logDensityOfStates, block)};
  };

  // Reweights the density of states to an arbitrary state point (beta, ln(beta f)): the weight of
  // bin b is w_b ~ Omega_b exp(-beta U_b + N_b ln(beta f)). Returns the average molecule count
  // and the effective number of independent samples supporting the estimate (an overlap
  // diagnostic: small values flag extrapolation beyond the sampled (N, U) region).
  // The isosteric heat at the same state point, from the grand-canonical fluctuation identity
  //
  //   q_st = kT - (<NU> - <N><U>) / (<N^2> - <N>^2),
  //
  // which is the energy an added molecule brings with it, at the loading the state point sits at. It needs
  // no extra sampling: the moments come from the same reweighted bins as the loading, and U is the energy
  // that was stored with every (N, U) sample. Where the pore is empty or full the number fluctuations
  // vanish and the ratio is meaningless, which the caller sees as an empty value.
  //
  // <U_IG> is the mean intramolecular energy of the molecule in the ideal gas: zero for a rigid molecule, and
  // for a flexible one a number the component has to carry ('idealGasEnergy'). Without it the heat is too low by
  // <U_IG>, as the enthalpy of adsorption of a Monte Carlo run is; that was written out without a word.
  const Component& heatComponent = systems.front().components[systems.front().tmmcComponent()];
  const double idealGasEnergyOfAdsorbate = heatComponent.idealGasEnergy.value_or(0.0);
  const bool isostericHeatLacksIdealGasEnergy =
      !heatComponent.rigid && heatComponent.atoms.size() > 1uz && !heatComponent.idealGasEnergy.has_value();
  auto reweightIsostericHeat = [&](const std::vector<double>& logDensityOfStates, double beta,
                                   double logBetaFugacity) -> std::optional<double>
  {
    const double logActivity = logBetaFugacity - logNormalizationAt(beta);
    double largest = logZero;
    std::vector<double> logWeights(numberOfBins, logZero);
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      if (logDensityOfStates[binId] <= logZeroThreshold) continue;
      logWeights[binId] = logDensityOfStates[binId] - beta * binMeanEnergy[binId] +
                          static_cast<double>(binMoleculeCount[binId]) * logActivity;
      largest = std::max(largest, logWeights[binId]);
    }
    if (largest <= logZeroThreshold) return std::nullopt;

    double partitionSum = 0.0;
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      if (logWeights[binId] <= logZeroThreshold) continue;
      partitionSum += std::exp(logWeights[binId] - largest);
    }
    const double logPartition = largest + std::log(partitionSum);

    double averageMolecules = 0.0;
    double averageEnergy = 0.0;
    double averageMoleculesSquared = 0.0;
    double averageMoleculesEnergy = 0.0;
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      if (logWeights[binId] <= logZeroThreshold) continue;
      const double weight = std::exp(logWeights[binId] - logPartition);
      const double molecules = static_cast<double>(binMoleculeCount[binId]);
      const double energy = binMeanEnergy[binId];
      averageMolecules += weight * molecules;
      averageEnergy += weight * energy;
      averageMoleculesSquared += weight * molecules * molecules;
      averageMoleculesEnergy += weight * molecules * energy;
    }

    const double varianceOfMolecules = averageMoleculesSquared - averageMolecules * averageMolecules;
    const double covariance = averageMoleculesEnergy - averageMolecules * averageEnergy;
    if (!(varianceOfMolecules > 1e-8)) return std::nullopt;

    // q_st = -(Delta H) with the same ideal-gas reference as PropertyEnthalpy: kT + <U_IG> - cov(N, U) / var(N)
    return 1.0 / beta + idealGasEnergyOfAdsorbate - covariance / varianceOfMolecules;
  };

  auto reweight = [&](const std::vector<double>& logDensityOfStates, double beta,
                      double logBetaFugacity) -> std::pair<double, double>
  {
    const double logActivity = logBetaFugacity - logNormalizationAt(beta);
    double largest = logZero;
    std::vector<double> logWeights(numberOfBins, logZero);
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      if (logDensityOfStates[binId] <= logZeroThreshold) continue;
      logWeights[binId] = logDensityOfStates[binId] - beta * binMeanEnergy[binId] +
                          static_cast<double>(binMoleculeCount[binId]) * logActivity;
      largest = std::max(largest, logWeights[binId]);
    }
    if (largest <= logZeroThreshold)
    {
      return {0.0, 0.0};
    }
    double partitionSum = 0.0;
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      if (logWeights[binId] <= logZeroThreshold) continue;
      partitionSum += std::exp(logWeights[binId] - largest);
    }
    const double logPartition = largest + std::log(partitionSum);

    double averageMolecules = 0.0;
    double inverseEffectiveSamples = 0.0;
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      if (logWeights[binId] <= logZeroThreshold) continue;
      const double weight = std::exp(logWeights[binId] - logPartition);
      averageMolecules += weight * static_cast<double>(binMoleculeCount[binId]);
      inverseEffectiveSamples += weight * weight / binTotalCount[binId];
    }
    return {averageMolecules, 1.0 / std::max(inverseEffectiveSamples, 1e-300)};
  };

  // full-data solve
  std::vector<double> logStateCounts(numberOfReplicas);
  std::vector<double> logBinCounts(numberOfBins);
  for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
  {
    const double count = std::accumulate(stateBlockCounts[stateId].begin(), stateBlockCounts[stateId].end(), 0.0);
    logStateCounts[stateId] = count > 0.0 ? std::log(count) : logZero;
  }
  for (std::size_t binId = 0; binId < numberOfBins; ++binId)
  {
    logBinCounts[binId] = std::log(binTotalCount[binId]);
  }

  std::vector<double> freeEnergies(numberOfReplicas, 0.0);
  std::vector<double> logDensityOfStates;
  const WhamConvergence convergence =
      solveWham(logStateCounts, logBinCounts, freeEnergies, logDensityOfStates, whamMaximumIterations, std::nullopt);
  const std::size_t iterations = convergence.iterations;
  const double residual = convergence.residual;
  whamIterations = iterations;
  whamResidual = residual;
  whamConverged = residual <= whamTolerance && convergence.connected;

  std::print(stream, "    pooled samples:          {}\n", totalNumberOfSamples);
  std::print(stream, "    occupied (N, U) bins:    {} ({} energy bins, deltaU = {:.6e} [K])\n", numberOfBins,
             numberOfEnergyBins, Units::EnergyToKelvin * deltaEnergy);
  std::print(stream, "    WHAM iterations:         {} (residual {:.3e}, tolerance {:.3e}, max {})\n\n", iterations,
             residual, whamTolerance, whamMaximumIterations);
  // A state point without a sample takes no part in the solution and is skipped by the overlap test: with a
  // fractional molecule that is never decoupled at the higher pressures, the whole ladder was reweighted from the
  // replicas that were left (from a single one the solve 'converges' in one iteration) without a word.
  {
    std::string unsampled;
    for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
    {
      if (logStateCounts[stateId] > logZeroThreshold) continue;
      unsampled += std::format("{}{}", unsampled.empty() ? "" : ", ", stateId);
    }
    if (!unsampled.empty())
    {
      std::print(stream, "    WARNING: no (N, U) sample was collected at state point(s) {}; the results there are\n",
                 unsampled);
      std::print(stream, "             extrapolated from the other state points (see the effective sample sizes).\n\n");
    }
  }
  if (!convergence.connected)
  {
    std::print(stream, "    WARNING: the sampled state points fall apart into groups whose (N, U) histograms share no bin;\n");
    std::print(stream, "             the free-energy offsets between the groups, and every reweighted result that spans\n");
    std::print(stream, "             them, are not measured. Densify the ladders across the gap.\n\n");
  }
  else if (!whamConverged)
  {
    std::print(stream, "    WARNING: the WHAM equations did not converge; the reweighted results are unreliable.\n");
    std::print(stream, "             Weak overlap between neighboring state points is the usual cause; use denser\n");
    std::print(stream, "             temperature/pressure ladders or longer production runs.\n\n");
  }

  // per-block solves (warm-started from the full solution) for the error bars
  std::vector<std::vector<double>> blockLogDensityOfStates(numberOfBlocks);
  std::vector<bool> blockIsValid(numberOfBlocks, false);
  whamUnconvergedBlocks = 0uz;
  whamUnconvergedBlockDetails.clear();
  for (std::size_t block = 0; block < numberOfBlocks; ++block)
  {
    std::vector<double> blockLogStateCounts(numberOfReplicas);
    std::vector<double> blockLogBinCounts(numberOfBins);
    double blockTotal = 0.0;
    for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
    {
      const double count = stateBlockCounts[stateId][block];
      blockTotal += count;
      blockLogStateCounts[stateId] = count > 0.0 ? std::log(count) : logZero;
    }
    if (blockTotal == 0.0) continue;
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      const double count = binBlockCounts[binId][block];
      blockLogBinCounts[binId] = count > 0.0 ? std::log(count) : logZero;
    }

    std::vector<double> blockFreeEnergies = freeEnergies;
    const auto [blockIterations, blockResidual, blockConnected] =
        solveWham(blockLogStateCounts, blockLogBinCounts, blockFreeEnergies, blockLogDensityOfStates[block],
                  whamMaximumIterations, block);
    blockIsValid[block] = blockResidual <= whamTolerance && blockConnected;
    if (!blockIsValid[block])
    {
      ++whamUnconvergedBlocks;
      whamUnconvergedBlockDetails.push_back(UnconvergedBlock{
          .index = block, .residual = blockResidual, .iterations = blockIterations});
    }
  }
  if (whamUnconvergedBlocks > 0uz)
  {
    std::print(stream, "    WARNING: {} of {} block WHAM solves did not converge; error bars may be underestimated.\n",
               whamUnconvergedBlocks, numberOfBlocks);
    for (const UnconvergedBlock& item : whamUnconvergedBlockDetails)
    {
      std::print(stream, "             production block {} residual {:.3e} ({} iterations)\n", item.index,
                 item.residual, item.iterations);
    }
    std::print(stream, "\n");
  }

  // mean and confidence-interval error of the reweighted molecule count at (beta, ln(beta f));
  // the error is estimated from the spread of the per-block reweighting solutions
  auto reweightWithError = [&](double beta, double logBetaFugacity) -> std::tuple<double, double, double>
  {
    const auto [mean, effectiveSamples] = reweight(logDensityOfStates, beta, logBetaFugacity);
    std::vector<double> blockValues;
    blockValues.reserve(numberOfBlocks);
    for (std::size_t block = 0; block < numberOfBlocks; ++block)
    {
      if (!blockIsValid[block]) continue;
      blockValues.push_back(reweight(blockLogDensityOfStates[block], beta, logBetaFugacity).first);
    }
    const double error = blockErrorEstimate(blockValues, mean);
    return {mean, error, effectiveSamples};
  };

  // unit conversions shared by all replicas (identical framework)
  const System& front = systems.front();
  const std::size_t adsorbateId = front.tmmcComponent();
  const Component& frontComponent = front.components[adsorbateId];
  double toMoleculesPerUnitCell = 1.0;
  double toMolePerKg = 0.0;
  double toMgPerG = 0.0;
  if (front.framework.has_value())
  {
    const int3 numberOfUnitCells = front.framework->numberOfUnitCells;
    toMoleculesPerUnitCell =
        1.0 / static_cast<double>(numberOfUnitCells.x * numberOfUnitCells.y * numberOfUnitCells.z);
    // (a framework without mass has no loading per kg, see System::frameworkMass)
    if (const std::optional<double> frameworkMass = front.frameworkMass(); frameworkMass.has_value())
    {
      toMolePerKg = 1000.0 / frameworkMass.value();
      toMgPerG = 1000.0 * frontComponent.totalMass / frameworkMass.value();
    }
  }

  // free energies per state point, with the self-consistency check: the reweighted loading at
  // every simulated state point must agree with the directly-measured average loading
  {
    std::ofstream freeEnergyFile("wham/reweighted_free_energies.reweighted_histogram.txt", std::ios::trunc);
    std::print(freeEnergyFile, "# Reweighted histogram: WHAM free energies per simulated state point\n");
    std::print(freeEnergyFile, "# g_i = -ln Xi_i (dimensionless, gauge: the first sampled state at 0)\n");
    std::print(freeEnergyFile,
               "# the direct column is the block-averaged loading measured by the replica itself;\n");
    std::print(freeEnergyFile,
               "# the reweighted column is the WHAM estimate at the same state point (self-consistency check)\n");
    std::print(freeEnergyFile, "# column 1: replica id\n");
    std::print(freeEnergyFile, "# column 2: temperature [K]\n");
    std::print(freeEnergyFile, "# column 3: pressure [Pa]\n");
    std::print(freeEnergyFile, "# column 4: fugacity [Pa]\n");
    std::print(freeEnergyFile, "# column 5: free energy g_i [-]\n");
    std::print(freeEnergyFile, "# column 6: samples [-]\n");
    std::print(freeEnergyFile, "# column 7, 8: direct loading, error [molecules / unit cell]\n");
    std::print(freeEnergyFile, "# column 9, 10: reweighted loading, error [molecules / unit cell]\n");

    std::print(stream, "    self-consistency check (loading in molecules / unit cell)\n");
    std::print(stream, "    replica    temperature [K]    pressure [Pa]        direct              reweighted\n");
    std::print(stream, "    ---------------------------------------------------------------------------------\n");
    for (std::size_t stateId = 0; stateId < numberOfReplicas; ++stateId)
    {
      const System& system = systems[stateId];
      const Component& component = system.components[adsorbateId];
      const double fugacityPa =
          component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.input_pressure;

      const auto [loadingAverage, loadingError] = system.averageLoadings.average();
      const double directLoading = toMoleculesPerUnitCell * loadingAverage.numberOfMolecules[adsorbateId];
      const double directError = toMoleculesPerUnitCell * loadingError.numberOfMolecules[adsorbateId];

      const auto [reweightedLoadingSupercell, reweightedErrorSupercell, effectiveSamples] =
          reweightWithError(stateBeta[stateId], stateLogBetaFugacity[stateId]);
      const double reweightedLoading = toMoleculesPerUnitCell * reweightedLoadingSupercell;
      const double reweightedError = toMoleculesPerUnitCell * reweightedErrorSupercell;
      const double samples = std::accumulate(stateBlockCounts[stateId].begin(), stateBlockCounts[stateId].end(), 0.0);

      std::print(freeEnergyFile, "{:6d}   {:10.4f}   {: .6e}   {: .6e}   {: .10e}   {:10.0f}   "
                                 "{: .6e} {: .6e}   {: .6e} {: .6e}\n",
                 stateId, system.temperature, system.input_pressure, fugacityPa, freeEnergies[stateId], samples,
                 directLoading, directError, reweightedLoading, reweightedError);

      std::print(stream, "    {:7d}    {:15.4f}    {:13.5e}    {:9.4f} ± {:7.4f}   {:9.4f} ± {:7.4f}\n", stateId,
                 system.temperature, system.input_pressure, directLoading, directError, reweightedLoading,
                 reweightedError);

      outputJson["output"]["reweighting"]["states"][stateId] = {
          {"temperature", system.temperature},   {"pressure", system.input_pressure},
          {"fugacity", fugacityPa},              {"freeEnergy", freeEnergies[stateId]},
          {"samples", samples},                  {"directLoading", directLoading},
          {"directLoadingError", directError},   {"reweightedLoading", reweightedLoading},
          {"reweightedLoadingError", reweightedError}, {"effectiveSamples", effectiveSamples}};
    }
    std::print(stream, "\n");
  }

  {
    std::ofstream dosFile("wham/density_of_states.reweighted_histogram.txt", std::ios::trunc);
    std::print(dosFile, "# Reweighted histogram: WHAM density of states Omega(N, U)\n");
    std::print(dosFile, "# column 1: N [molecules]\n");
    std::print(dosFile, "# column 2: mean energy of the bin [internal units]\n");
    std::print(dosFile, "# column 3: sample count in the bin\n");
    std::print(dosFile, "# column 4: ln Omega(N, U) (gauge: the first sampled state at g = 0)\n\n");
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      std::print(dosFile, "{:6d}   {: .10e}   {: .6e}   {: .10e}\n", binMoleculeCount[binId], binMeanEnergy[binId],
                 binTotalCount[binId], logDensityOfStates[binId]);
    }
  }

  // reweighted isotherms: one file per requested temperature, evaluated on a fine log-spaced
  // pressure grid (fugacity coefficients from the Peng-Robinson equation of state at every point)
  // the fluid the walkers ran with: every component at its mole fraction, so that the reweighted fugacity of the
  // adsorbate, y phi(P) P, reproduces the walkers' own y phi P at their pressure
  std::vector<EquationOfState::FluidInput> fluidInputs;
  for (const Component& component : front.components)
  {
    fluidInputs.push_back({component.criticalTemperature, component.criticalPressure, component.acentricFactor,
                           component.molFraction, component.swappable});
  }
  // A component without critical constants has no Peng-Robinson: its walkers ran with the fugacity
  // coefficient it was given (or the ideal gas's one), and the reweighting uses that same coefficient at
  // every pressure, so that the two sides of it are one fluid.
  const bool recomputesFugacityCoefficient =
      !frontComponent.fugacityCoefficientGiven &&
      EquationOfState::hasCriticalConstants(frontComponent.criticalTemperature, frontComponent.criticalPressure);
  const double keptFugacityCoefficient = frontComponent.fugacityCoefficient.value_or(1.0);
  const double logPressureMinimum = std::log(reweightingPressureRange.first);
  const double logPressureMaximum = std::log(reweightingPressureRange.second);
  const std::size_t numberOfTargetPressures =
      reweightingPressureRange.second > reweightingPressureRange.first ? reweightingNumberOfPressures : 1uz;

  outputJson["output"]["reweighting"]["isostericHeatLacksIdealGasEnergy"] = isostericHeatLacksIdealGasEnergy;
  if (isostericHeatLacksIdealGasEnergy)
  {
    std::print(stream,
               "    WARNING: '{}' is flexible and has no ideal-gas energy: the reweighted isosteric heats are\n"
               "             kT - cov(N, U) / var(N), to which the mean intramolecular energy of the ideal-gas\n"
               "             molecule, <U_IG>, still has to be added\n\n",
               heatComponent.name);
  }

  // in reduced units there is no kJ/mol (the conversion factor is zero): the heat is written in units of epsilon
  const bool reducedUnits = Units::unitSystem == Units::System::ReducedUnits;
  const double isostericHeatConversion = reducedUnits ? 1.0 : Units::EnergyToKJPerMol;
  for (const double targetTemperature : reweightingTemperatures)
  {
    std::ofstream isotherm(std::format("wham/reweighted_isotherm_{}.reweighted_histogram.txt", targetTemperature),
                           std::ios::trunc);
    std::print(isotherm, "# Reweighted histogram: reweighted adsorption isotherm at {} [K] ({})\n", targetTemperature,
               frontComponent.name);
    std::print(isotherm, "# errors from re-solving the WHAM equations per block\n");
    std::print(isotherm,
               "# the effective sample size (last column) diagnoses the overlap: small values flag\n"
               "# extrapolation beyond the sampled (N, U) region\n");
    std::print(isotherm, "# column 1: fugacity [Pa]\n");
    std::print(isotherm, "# column 2, 3: absolute loading, error [molecules / simulation supercell]\n");
    std::print(isotherm, "# column 4, 5: absolute loading, error [molecules / unit cell]\n");
    std::print(isotherm, "# column 6, 7: absolute loading, error [mol/kg-framework]\n");
    std::print(isotherm, "# column 8, 9: absolute loading, error [mg/g-framework]\n");
    std::print(isotherm, "# column 10: pressure [Pa]\n");
    std::print(isotherm, "# column 11: effective sample size [-]\n");
    std::print(isotherm,
               "# column 12: isosteric heat [{}], kT - cov(N, U) / var(N) over the reweighted\n"
               "#            distribution; 'nan' where the loading is pinned and N does not fluctuate\n{}\n",
               reducedUnits ? "reduced energy units" : "kJ/mol",
               isostericHeatLacksIdealGasEnergy
                   ? "#            WARNING: the molecule is flexible and has no ideal-gas energy: the mean\n"
                     "#            intramolecular energy of the ideal-gas molecule, <U_IG>, still has to be added\n"
                   : "");

    const double targetBeta = 1.0 / (Units::KB * targetTemperature);
    ReweightedIsotherm stored;
    stored.temperature = targetTemperature;
    stored.points.reserve(numberOfTargetPressures);
    for (std::size_t pressureIndex = 0; pressureIndex < numberOfTargetPressures; ++pressureIndex)
    {
      const double targetPressure =
          numberOfTargetPressures == 1uz
              ? reweightingPressureRange.first
              : std::exp(logPressureMinimum + static_cast<double>(pressureIndex) *
                                                  (logPressureMaximum - logPressureMinimum) /
                                                  static_cast<double>(numberOfTargetPressures - 1uz));

      const std::vector<EquationOfState::FluidResult> fluidResults = EquationOfState::computeFluidProperties(
          targetTemperature, targetPressure, fluidInputs, EquationOfState::Type::PengRobinson,
          EquationOfState::MixingRules::VanDerWaals);
      const double fugacityCoefficient = recomputesFugacityCoefficient
                                             ? fluidResults[adsorbateId].fugacityCoefficient.value_or(1.0)
                                             : keptFugacityCoefficient;
      const double fugacityPa = frontComponent.molFraction * fugacityCoefficient * targetPressure;
      const double logBetaFugacity = std::log(targetBeta * fugacityPa / Units::PressureConversionFactor);

      const auto [loading, loadingError, effectiveSamples] = reweightWithError(targetBeta, logBetaFugacity);
      const std::optional<double> isostericHeat =
          reweightIsostericHeat(logDensityOfStates, targetBeta, logBetaFugacity);
      const double loadingPerCell = toMoleculesPerUnitCell * loading;
      const double loadingErrorPerCell = toMoleculesPerUnitCell * loadingError;

      stored.points.push_back(ReweightedIsothermPoint{targetPressure, fugacityPa, loadingPerCell,
                                                      loadingErrorPerCell, effectiveSamples, isostericHeat});

      const std::string isostericHeatText =
          isostericHeat.has_value() ? std::format("{: .6e}", isostericHeatConversion * isostericHeat.value())
                                    : std::string(" nan");
      std::print(isotherm,
                 "{: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e}   "
                 "{: .6e}   {}\n",
                 fugacityPa, loading, loadingError, loadingPerCell, loadingErrorPerCell, toMolePerKg * loading,
                 toMolePerKg * loadingError, toMgPerG * loading, toMgPerG * loadingError, targetPressure,
                 effectiveSamples, isostericHeatText);
    }
    reweightedIsotherms.push_back(std::move(stored));
  }

  if (computeBET)
  {
    std::print(stream, "BET extraction\n");
    if (!front.framework.has_value())
    {
      std::print(stream, "    skipped: BET needs a framework (unit-cell mass and volume)\n\n");
    }
    else if (front.components.empty())
    {
      std::print(stream, "    skipped: BET needs an adsorbate component\n\n");
    }
    else
    {
      const BETProbeProperties probe = requireBETProbeProperties(frontComponent);
      const int3 numberOfUnitCells = front.framework->numberOfUnitCells;
      const double cells =
          static_cast<double>(numberOfUnitCells.x * numberOfUnitCells.y * numberOfUnitCells.z);
      const double unitCellVolume = front.simulationBox.volume / cells;
      // (the mass the loadings per kg are counted against: the framework with its cations of a fixed number, see
      // System::frameworkMass. The bare unit cell left the cations out, and the gravimetric area of a
      // cation-exchanged zeolite was too large by their share of the mass.)
      const double unitCellMass = front.frameworkMass().value_or(front.framework->unitCellMass * cells) / cells;
      outputJson["output"]["reweighting"]["bet"] = nlohmann::json::array();
      for (const ReweightedIsotherm& isotherm : reweightedIsotherms)
      {
        // 'SaturationPressure' is P0 at one temperature, the front of the sampled ladder (where the automatic
        // pressure placement put P0 as well): P / P0 means nothing at any other temperature
        if (isotherm.temperature != temperatures.front())
        {
          std::print(stream,
                     "    {} [K]: skipped; 'SaturationPressure' P0 = {:.5e} [Pa] is taken at {} [K], the first "
                     "temperature of the ladder\n\n",
                     isotherm.temperature, probe.saturationPressure, temperatures.front());
          continue;
        }
        std::vector<SimulatedIsothermPoint> points;
        points.reserve(isotherm.points.size());
        for (const ReweightedIsothermPoint& point : isotherm.points)
        {
          points.push_back(SimulatedIsothermPoint{point.pressure, point.moleculesPerCell});
        }
        const BETSurfaceArea bet = fitNitrogenBET(points, unitCellMass, unitCellVolume, probe);
        const FiniteLayerBETFit finiteLayer =
            fitNitrogenFiniteLayerBET(points, unitCellMass, unitCellVolume, probe);

        std::optional<double> areaError;
        std::optional<double> monolayerError;
        std::optional<double> cError;
        nlohmann::json blockAreas = nlohmann::json::array();
        nlohmann::json blockMonolayers = nlohmann::json::array();
        nlohmann::json blockCConstants = nlohmann::json::array();
        std::size_t betBlocksUsed = 0uz;
        std::size_t betBlocksFailed = 0uz;

        std::optional<double> finiteLayerAreaError;
        std::optional<double> finiteLayerMonolayerError;
        std::optional<double> finiteLayerCError;
        nlohmann::json finiteLayerBlockAreas = nlohmann::json::array();
        nlohmann::json finiteLayerBlockMonolayers = nlohmann::json::array();
        nlohmann::json finiteLayerBlockCConstants = nlohmann::json::array();
        std::size_t finiteLayerBlocksUsed = 0uz;
        std::size_t finiteLayerBlocksFailed = 0uz;

        // Jackknife classical BET (fixed Rouquerol window) and finite-layer BET (frozen n + window).
        if ((!bet.plateauReading && bet.monolayerCapacity > 0.0 && bet.windowHigh > bet.windowLow) ||
            finiteLayer.ok)
        {
          const double targetBeta = 1.0 / (Units::KB * isotherm.temperature);
          std::vector<double> blockAreaValues;
          std::vector<double> blockMonolayerValues;
          std::vector<double> blockCValues;
          std::vector<double> finiteLayerAreaValues;
          std::vector<double> finiteLayerMonolayerValues;
          std::vector<double> finiteLayerCValues;
          blockAreaValues.reserve(numberOfBlocks);
          blockMonolayerValues.reserve(numberOfBlocks);
          blockCValues.reserve(numberOfBlocks);
          finiteLayerAreaValues.reserve(numberOfBlocks);
          finiteLayerMonolayerValues.reserve(numberOfBlocks);
          finiteLayerCValues.reserve(numberOfBlocks);

          for (std::size_t block = 0; block < numberOfBlocks; ++block)
          {
            if (!blockIsValid[block]) continue;

            std::vector<SimulatedIsothermPoint> blockPoints;
            blockPoints.reserve(isotherm.points.size());
            for (const ReweightedIsothermPoint& point : isotherm.points)
            {
              const double logBetaFugacity =
                  std::log(targetBeta * point.fugacity / Units::PressureConversionFactor);
              const double loading =
                  reweight(blockLogDensityOfStates[block], targetBeta, logBetaFugacity).first;
              blockPoints.push_back(
                  SimulatedIsothermPoint{point.pressure, toMoleculesPerUnitCell * loading});
            }

            if (!bet.plateauReading && bet.monolayerCapacity > 0.0 && bet.windowHigh > bet.windowLow)
            {
              const BETSurfaceArea blockBet = fitNitrogenBETFixedWindow(
                  blockPoints, unitCellMass, unitCellVolume, bet.windowLow, bet.windowHigh, probe);
              if (!(blockBet.monolayerCapacity > 0.0) || !(blockBet.gravimetricArea > 0.0))
              {
                ++betBlocksFailed;
              }
              else
              {
                ++betBlocksUsed;
                blockAreaValues.push_back(blockBet.gravimetricArea);
                blockMonolayerValues.push_back(blockBet.monolayerCapacity);
                blockCValues.push_back(blockBet.cConstant);
                blockAreas.push_back(blockBet.gravimetricArea);
                blockMonolayers.push_back(blockBet.monolayerCapacity);
                blockCConstants.push_back(blockBet.cConstant);
              }
            }

            if (finiteLayer.ok)
            {
              const FiniteLayerBETFit blockFit = fitNitrogenFiniteLayerBETFixedLayers(
                  blockPoints, unitCellMass, unitCellVolume, finiteLayer.numberOfLayers, finiteLayer.windowLow,
                  finiteLayer.windowHigh, probe);
              if (!blockFit.ok || !(blockFit.gravimetricArea > 0.0))
              {
                ++finiteLayerBlocksFailed;
              }
              else
              {
                ++finiteLayerBlocksUsed;
                finiteLayerAreaValues.push_back(blockFit.gravimetricArea);
                finiteLayerMonolayerValues.push_back(blockFit.monolayerCapacity);
                finiteLayerCValues.push_back(blockFit.cConstant);
                finiteLayerBlockAreas.push_back(blockFit.gravimetricArea);
                finiteLayerBlockMonolayers.push_back(blockFit.monolayerCapacity);
                finiteLayerBlockCConstants.push_back(blockFit.cConstant);
              }
            }
          }

          if (betBlocksUsed >= 3uz)
          {
            areaError = blockErrorEstimate(blockAreaValues, bet.gravimetricArea);
            monolayerError = blockErrorEstimate(blockMonolayerValues, bet.monolayerCapacity);
            cError = blockErrorEstimate(blockCValues, bet.cConstant);
          }
          if (finiteLayerBlocksUsed >= 3uz)
          {
            finiteLayerAreaError = blockErrorEstimate(finiteLayerAreaValues, finiteLayer.gravimetricArea);
            finiteLayerMonolayerError =
                blockErrorEstimate(finiteLayerMonolayerValues, finiteLayer.monolayerCapacity);
            finiteLayerCError = blockErrorEstimate(finiteLayerCValues, finiteLayer.cConstant);
          }
        }

        writeNitrogenBETSummary(stream, bet, probe, "    ", areaError, monolayerError, cError);
        if (betBlocksUsed > 0uz || betBlocksFailed > 0uz)
        {
          std::print(stream,
                     "    BET block errors:       {} blocks (fixed Rouquerol window); {} fixed-window "
                     "refits failed{}\n",
                     betBlocksUsed, betBlocksFailed,
                     areaError.has_value() ? "" : " — need ≥ 3 successful blocks for a CI");
        }
        std::print(stream, "\n");
        writeNitrogenFiniteLayerBETSummary(stream, finiteLayer, probe, "    ", finiteLayerAreaError,
                                           finiteLayerMonolayerError, finiteLayerCError);
        if (finiteLayerBlocksUsed > 0uz || finiteLayerBlocksFailed > 0uz)
        {
          std::print(stream,
                     "    Finite-layer block errors: {} blocks (frozen n, window); {} fixed-layer "
                     "refits failed{}\n",
                     finiteLayerBlocksUsed, finiteLayerBlocksFailed,
                     finiteLayerAreaError.has_value() ? "" : " — need ≥ 3 successful blocks for a CI");
        }
        std::print(stream, "\n");
        const std::string betFile =
            std::format("wham/bet_{}.reweighted_histogram.txt", isotherm.temperature);
        std::ofstream table(betFile, std::ios::trunc);
        std::print(table, "# Reweighted histogram: BET plot at {} [K] ({})\n", isotherm.temperature,
                   frontComponent.name);
        writeNitrogenBETTable(table, bet, probe);
        nlohmann::json entry = nitrogenBETJson(bet, probe);
        entry["temperature"] = isotherm.temperature;
        entry["file"] = betFile;
        entry["fixedWindowBlockErrors"] = true;
        entry["betBlocksUsed"] = betBlocksUsed;
        entry["betBlocksFailed"] = betBlocksFailed;
        entry["blockGravimetricAreas"] = std::move(blockAreas);
        entry["blockMonolayerCapacities"] = std::move(blockMonolayers);
        entry["blockCConstants"] = std::move(blockCConstants);
        if (areaError.has_value())
        {
          entry["gravimetricAreaError"] = *areaError;
          entry["monolayerCapacityError"] = *monolayerError;
          entry["cConstantError"] = *cError;
        }
        nlohmann::json finiteLayerEntry = nitrogenFiniteLayerBETJson(finiteLayer, probe);
        finiteLayerEntry["fixedLayerBlockErrors"] = true;
        finiteLayerEntry["betBlocksUsed"] = finiteLayerBlocksUsed;
        finiteLayerEntry["betBlocksFailed"] = finiteLayerBlocksFailed;
        finiteLayerEntry["blockGravimetricAreas"] = std::move(finiteLayerBlockAreas);
        finiteLayerEntry["blockMonolayerCapacities"] = std::move(finiteLayerBlockMonolayers);
        finiteLayerEntry["blockCConstants"] = std::move(finiteLayerBlockCConstants);
        if (finiteLayerAreaError.has_value())
        {
          finiteLayerEntry["gravimetricAreaError"] = *finiteLayerAreaError;
          finiteLayerEntry["monolayerCapacityError"] = *finiteLayerMonolayerError;
          finiteLayerEntry["cConstantError"] = *finiteLayerCError;
        }
        entry["finiteLayer"] = std::move(finiteLayerEntry);
        outputJson["output"]["reweighting"]["bet"].push_back(std::move(entry));
        std::print(stream, "    BET plot written to {}\n\n", betFile);
      }
    }
  }

  // vapor-liquid equilibrium (bulk boxes only): at every requested subcritical temperature the
  // coexistence fugacity is located by the equal-weight criterion (Wilding): the reweighted
  // molecule-number distribution P(N) is bimodal, and the fugacity is bisected until the vapor
  // peak (N < N_cut) and the liquid peak (N >= N_cut) carry equal probability weight. The cut is
  // placed at the deepest valley of ln P(N) between the two peaks.
  if (!front.framework.has_value())
  {
    const double volume = front.simulationBox.volume;
    const std::uint32_t maximumMoleculeCount = *std::ranges::max_element(binMoleculeCount);

    std::vector<std::vector<std::size_t>> binsByMoleculeCount(maximumMoleculeCount + 1uz);
    for (std::size_t binId = 0; binId < numberOfBins; ++binId)
    {
      binsByMoleculeCount[binMoleculeCount[binId]].push_back(binId);
    }

    // ln P(N) (unnormalized) at (beta, ln(beta f)): the energy dimension is summed out
    auto moleculeDistribution = [&](const std::vector<double>& logDOS, double beta,
                                    double logBetaFugacity) -> std::vector<double>
    {
      const double logActivity = logBetaFugacity - logNormalizationAt(beta);
      std::vector<double> logProbability(maximumMoleculeCount + 1uz, logZero);
      for (std::size_t n = 0; n <= maximumMoleculeCount; ++n)
      {
        double largest = logZero;
        for (const std::size_t binId : binsByMoleculeCount[n])
        {
          if (logDOS[binId] <= logZeroThreshold) continue;
          largest = std::max(largest, logDOS[binId] - beta * binMeanEnergy[binId] +
                                          static_cast<double>(n) * logActivity);
        }
        if (largest <= logZeroThreshold) continue;
        double sum = 0.0;
        for (const std::size_t binId : binsByMoleculeCount[n])
        {
          if (logDOS[binId] <= logZeroThreshold) continue;
          sum += std::exp(logDOS[binId] - beta * binMeanEnergy[binId] +
                          static_cast<double>(n) * logActivity - largest);
        }
        logProbability[n] = largest + std::log(sum);
      }
      return logProbability;
    };

    // The phase cut at the deepest valley of ln P(N), or nothing when the distribution is one phase; see
    // findMacrostatePhaseSplit for what an unsampled N counts as.
    constexpr double bimodalDepthThreshold = 0.2;
    auto findPhaseSplit = [&](const std::vector<double>& logProbability)
    { return findMacrostatePhaseSplit(logProbability, logZeroThreshold, bimodalDepthThreshold); };

    auto logSumExpRange = [&](const std::vector<double>& logProbability, std::size_t begin,
                              std::size_t end) -> double
    {
      double largest = logZero;
      for (std::size_t n = begin; n < end; ++n)
      {
        if (logProbability[n] > logZeroThreshold) largest = std::max(largest, logProbability[n]);
      }
      if (largest <= logZeroThreshold) return logZero;
      double sum = 0.0;
      for (std::size_t n = begin; n < end; ++n)
      {
        if (logProbability[n] > logZeroThreshold) sum += std::exp(logProbability[n] - largest);
      }
      return largest + std::log(sum);
    };

    struct CoexistencePoint
    {
      double logBetaFugacity;
      double fugacityPa;
      std::optional<double> saturationPressurePa;  // nothing: the sampled N do not reach a state of known pressure
      bool saturationPressureFromEmptyBox;  // exact normalization (false: ideal-gas reference)
      double vaporDensity;                  // [kg/m^3]
      double liquidDensity;                 // [kg/m^3]
      double vaporMolecules;
      double liquidMolecules;
      double valleyDepth;
      std::size_t cut;
      std::vector<double> logProbability;
    };

    // Scans the pressure range for a bracket where the liquid-vapor weight difference changes
    // sign (both endpoints bimodal), then bisects the fugacity to the equal-weight point.
    auto solveCoexistence = [&](const std::vector<double>& logDOS,
                                double temperature) -> std::optional<CoexistencePoint>
    {
      const double beta = 1.0 / (Units::KB * temperature);

      auto scanPoint = [&](double logBetaFugacity) -> MacrostateScanPoint
      {
        const std::vector<double> logProbability = moleculeDistribution(logDOS, beta, logBetaFugacity);
        const double logTotal = logSumExpRange(logProbability, 0uz, logProbability.size());
        double meanMolecules = 0.0;
        for (std::size_t n = 0; n < logProbability.size(); ++n)
        {
          if (logProbability[n] <= logZeroThreshold) continue;
          meanMolecules += static_cast<double>(n) * std::exp(logProbability[n] - logTotal);
        }
        const std::optional<std::pair<std::size_t, double>> split = findPhaseSplit(logProbability);
        if (!split.has_value()) return {std::nullopt, meanMolecules};
        const double logVaporWeight = logSumExpRange(logProbability, 0uz, split->first);
        const double logLiquidWeight = logSumExpRange(logProbability, split->first, logProbability.size());
        return {logLiquidWeight - logVaporWeight, meanMolecules};
      };

      auto logBetaFugacityOfPressure = [&](double pressure) -> double
      {
        const std::vector<EquationOfState::FluidResult> fluidResults = EquationOfState::computeFluidProperties(
            temperature, pressure, fluidInputs, EquationOfState::Type::PengRobinson,
            EquationOfState::MixingRules::VanDerWaals);
        const double fugacityPa = frontComponent.molFraction *
                                  (recomputesFugacityCoefficient
                                       ? fluidResults[adsorbateId].fugacityCoefficient.value_or(1.0)
                                       : keptFugacityCoefficient) *
                                  pressure;
        return std::log(beta * fugacityPa / Units::PressureConversionFactor);
      };

      // The search runs over ln(beta f), the variable the density of states is reweighted in. It ran over the
      // pressure, with the fugacity of each from the stable Peng-Robinson root: past the saturation pressure of
      // that equation of state the root is the liquid's, whose fugacity hardly rises with the pressure, so the
      // scan stalled at the Peng-Robinson saturation fugacity and a model that coexists above it was never
      // bracketed, at any pressure range. It runs from the fugacity of the dilute end to the larger of the
      // ideal-gas and the Peng-Robinson fugacity of the upper end.
      constexpr std::size_t numberOfScanPoints = 200uz;
      const double scanLower = logBetaFugacityOfPressure(std::exp(logPressureMinimum));
      const double scanUpper =
          std::max(logBetaFugacityOfPressure(std::exp(logPressureMaximum)),
                   std::log(beta * frontComponent.molFraction * std::exp(logPressureMaximum) /
                            Units::PressureConversionFactor));
      const std::optional<double> located = locateEqualWeightPoint(scanLower, scanUpper, numberOfScanPoints, scanPoint);
      if (!located.has_value()) return std::nullopt;

      const double logBetaFugacity = located.value();
      const std::vector<double> logProbability = moleculeDistribution(logDOS, beta, logBetaFugacity);
      const std::optional<std::pair<std::size_t, double>> split = findPhaseSplit(logProbability);
      if (!split.has_value()) return std::nullopt;
      const std::size_t cut = split->first;

      const double logVaporWeight = logSumExpRange(logProbability, 0uz, cut);
      const double logLiquidWeight = logSumExpRange(logProbability, cut, logProbability.size());
      double vaporMolecules = 0.0;
      double liquidMolecules = 0.0;
      for (std::size_t n = 0; n < logProbability.size(); ++n)
      {
        if (logProbability[n] <= logZeroThreshold) continue;
        if (n < cut)
        {
          vaporMolecules += static_cast<double>(n) * std::exp(logProbability[n] - logVaporWeight);
        }
        else
        {
          liquidMolecules += static_cast<double>(n) * std::exp(logProbability[n] - logLiquidWeight);
        }
      }

      // Saturation pressure from beta p V = ln Xi. The absolute normalization of Xi is fixed by
      // the empty-box state (N = 0, U = 0, weight exactly one) when it was sampled; otherwise it
      // is fixed approximately at the dilute end of the scanned pressure range, where the vapor
      // is nearly ideal and beta p V ~ beta f V.
      std::optional<double> saturationPressurePa;
      // The pressure of the coexisting phases is that of either pure phase: beta p V = ln Xi_vap, the grand
      // partition function restricted to the vapor basin (Errington, J. Chem. Phys. 118, 9915 (2003)). The sum over
      // all N is Xi_vap + Xi_liq = 2 Xi_vap at the equal-weight point and would add k_B T ln 2 / V to p_sat.
      const double logPartitionAtCoexistence = logVaporWeight;
      const bool saturationPressureFromEmptyBox = logProbability[0] > logZeroThreshold;
      if (saturationPressureFromEmptyBox)
      {
        saturationPressurePa =
            ((logPartitionAtCoexistence - logProbability[0]) / (beta * volume)) * Units::PressureConversionFactor;
      }
      else
      {
        // The reference stands in for the full grand partition function only where the states below the
        // lowest sampled N carry no weight. With <N> = beta f V of the reference below that N the sum over the
        // sampled states is its truncated tail, ln Xi ~ N_min ln(beta f V) - ln N_min! in place of beta f V,
        // and the saturation pressure rose without bound as the dilute end of the range was lowered (megapascals
        // for five molecules at 100 Pa). The same holds for the vapor peak at coexistence. Where either
        // distribution still has more than a thousandth of its weight on the lowest sampled N there is no
        // state of known pressure in reach, and no saturation pressure is given.
        const std::size_t lowestSampled = static_cast<std::size_t>(
            std::ranges::find_if(logProbability, [&](double value) { return value > logZeroThreshold; }) -
            logProbability.begin());
        constexpr double logTruncatedShare = -6.907755278982137;  // ln 1e-3
        const double referenceLogBetaFugacity = logBetaFugacityOfPressure(std::exp(logPressureMinimum));
        const std::vector<double> referenceLogProbability =
            moleculeDistribution(logDOS, beta, referenceLogBetaFugacity);
        const double logPartitionAtReference =
            logSumExpRange(referenceLogProbability, 0uz, referenceLogProbability.size());
        const bool truncated =
            lowestSampled >= logProbability.size() ||
            referenceLogProbability[lowestSampled] - logPartitionAtReference > logTruncatedShare ||
            logProbability[lowestSampled] - logPartitionAtCoexistence > logTruncatedShare;
        if (!truncated)
        {
          const double referenceBetaPressureVolume = std::exp(referenceLogBetaFugacity) * volume;
          saturationPressurePa =
              ((logPartitionAtCoexistence - logPartitionAtReference + referenceBetaPressureVolume) /
               (beta * volume)) *
              Units::PressureConversionFactor;
        }
      }

      const double toKgPerCubicMeter = frontComponent.totalMass * Units::DensityConversionFactor / volume;
      return CoexistencePoint{
          .logBetaFugacity = logBetaFugacity,
          .fugacityPa = (std::exp(logBetaFugacity) / beta) * Units::PressureConversionFactor,
          .saturationPressurePa = saturationPressurePa,
          .saturationPressureFromEmptyBox = saturationPressureFromEmptyBox,
          .vaporDensity = vaporMolecules * toKgPerCubicMeter,
          .liquidDensity = liquidMolecules * toKgPerCubicMeter,
          .vaporMolecules = vaporMolecules,
          .liquidMolecules = liquidMolecules,
          .valleyDepth = split->second,
          .cut = cut,
          .logProbability = logProbability};
    };

    std::ofstream coexistence("output/vle_coexistence.reweighted_histogram.txt", std::ios::trunc);
    std::print(coexistence, "# Reweighted histogram: vapor-liquid coexistence of {} (equal-weight criterion)\n",
               frontComponent.name);
    std::print(coexistence, "# box volume {:.4f} [A^3]; errors from re-solving the WHAM equations per block\n",
               volume);
    std::print(coexistence,
               "# the saturation pressure follows from beta p V = ln Xi, normalized by the empty-box\n"
               "# state when sampled (column 13 = 0), otherwise approximately by an ideal-gas reference\n"
               "# at the dilute end of the scanned pressure range (column 13 = 1); 'nan' where the sampled\n"
               "# molecule counts reach neither (column 13 = 2)\n");
    std::print(coexistence, "# column 1: temperature [K]\n");
    std::print(coexistence, "# column 2, 3: coexistence fugacity, error [Pa]\n");
    std::print(coexistence, "# column 4, 5: saturation pressure, error [Pa]\n");
    std::print(coexistence, "# column 6, 7: vapor density, error [kg/m^3]\n");
    std::print(coexistence, "# column 8, 9: liquid density, error [kg/m^3]\n");
    std::print(coexistence, "# column 10, 11: vapor, liquid peak [molecules]\n");
    std::print(coexistence, "# column 12: ln P(N) valley depth [-] (vanishes towards the critical point)\n");
    std::print(coexistence,
               "# column 13: saturation-pressure normalization (0 exact, 1 ideal-gas reference, 2 none)\n\n");

    bool anyUnavailableSaturationPressure = false;

    std::print(stream, "    vapor-liquid coexistence (equal-weight criterion)\n");
    std::print(stream, "    temperature [K]    fugacity [Pa]        P_sat [Pa]           "
                       "rho_vap [kg/m^3]     rho_liq [kg/m^3]\n");
    std::print(stream, "    -----------------------------------------------------------------"
                       "----------------------------------\n");

    for (const double targetTemperature : reweightingTemperatures)
    {
      const std::optional<CoexistencePoint> point = solveCoexistence(logDensityOfStates, targetTemperature);
      if (!point.has_value())
      {
        std::print(coexistence,
                   "# {} [K]: no coexistence found (supercritical, pressure range does not bracket the\n"
                   "#   transition, or the (N, U) sampling does not connect the two phases)\n",
                   targetTemperature);
        std::print(stream, "    {:15.4f}    no coexistence found in the scanned pressure range\n",
                   targetTemperature);
        continue;
      }

      // per-block coexistence solves for the error bars
      std::vector<double> blockFugacities, blockPressures, blockVaporDensities, blockLiquidDensities;
      for (std::size_t block = 0; block < numberOfBlocks; ++block)
      {
        if (!blockIsValid[block]) continue;
        const std::optional<CoexistencePoint> blockPoint =
            solveCoexistence(blockLogDensityOfStates[block], targetTemperature);
        if (!blockPoint.has_value()) continue;
        blockFugacities.push_back(blockPoint->fugacityPa);
        blockVaporDensities.push_back(blockPoint->vaporDensity);
        blockLiquidDensities.push_back(blockPoint->liquidDensity);
        if (blockPoint->saturationPressurePa.has_value())
        {
          blockPressures.push_back(blockPoint->saturationPressurePa.value());
        }
      }
      const double fugacityError = blockErrorEstimate(blockFugacities, point->fugacityPa);
      const bool hasSaturationPressure = point->saturationPressurePa.has_value();
      const double pressureError =
          hasSaturationPressure ? blockErrorEstimate(blockPressures, point->saturationPressurePa.value()) : 0.0;
      const std::string saturationPressureText =
          hasSaturationPressure ? std::format("{: .6e} {: .6e}", point->saturationPressurePa.value(), pressureError)
                                : std::string(" nan nan");
      anyUnavailableSaturationPressure = anyUnavailableSaturationPressure || !hasSaturationPressure;
      const double vaporDensityError = blockErrorEstimate(blockVaporDensities, point->vaporDensity);
      const double liquidDensityError = blockErrorEstimate(blockLiquidDensities, point->liquidDensity);

      std::print(coexistence,
                 "{:10.4f}   {: .6e} {: .6e}   {}   {: .6e} {: .6e}   {: .6e} {: .6e}   "
                 "{: .6e} {: .6e}   {: .6e}   {}\n",
                 targetTemperature, point->fugacityPa, fugacityError, saturationPressureText,
                 point->vaporDensity, vaporDensityError, point->liquidDensity, liquidDensityError,
                 point->vaporMolecules, point->liquidMolecules, point->valleyDepth,
                 point->saturationPressureFromEmptyBox ? 0 : (hasSaturationPressure ? 1 : 2));

      std::print(stream, "    {:15.4f}    {: .6e} ± {:.2e}   {}{}  {:9.4f} ± {:7.4f}   "
                         "{:9.4f} ± {:7.4f}\n",
                 targetTemperature, point->fugacityPa, fugacityError,
                 hasSaturationPressure
                     ? std::format("{: .6e} ± {:.2e}", point->saturationPressurePa.value(), pressureError)
                     : std::format("{:>24}", "not available"),
                 point->saturationPressureFromEmptyBox ? ' ' : (hasSaturationPressure ? '*' : '#'),
                 point->vaporDensity, vaporDensityError, point->liquidDensity, liquidDensityError);

      // the coexistence molecule-number distribution (for inspection and finite-size scaling)
      std::ofstream distribution(std::format("output/vle_distribution_{}.reweighted_histogram.txt",
                                             targetTemperature),
                                 std::ios::trunc);
      std::print(distribution, "# Reweighted histogram: P(N) at coexistence, {} at {} [K]\n", frontComponent.name,
                 targetTemperature);
      std::print(distribution, "# coexistence fugacity {:.6e} [Pa], phase cut at N = {}\n", point->fugacityPa,
                 point->cut);
      std::print(distribution, "# column 1: N [molecules]\n");
      std::print(distribution, "# column 2: P(N) (normalized)\n");
      std::print(distribution, "# column 3: ln P(N)\n\n");
      const double logNormalization =
          logSumExpRange(point->logProbability, 0uz, point->logProbability.size());
      for (std::size_t n = 0; n < point->logProbability.size(); ++n)
      {
        if (point->logProbability[n] <= logZeroThreshold) continue;
        const double logNormalized = point->logProbability[n] - logNormalization;
        std::print(distribution, "{:6d}   {: .6e}   {: .6e}\n", n, std::exp(logNormalized), logNormalized);
      }

      nlohmann::json vleEntry;
      vleEntry["temperature"] = targetTemperature;
      vleEntry["fugacity"] = point->fugacityPa;
      vleEntry["fugacityError"] = fugacityError;
      if (hasSaturationPressure)
      {
        vleEntry["saturationPressure"] = point->saturationPressurePa.value();
        vleEntry["saturationPressureError"] = pressureError;
      }
      else
      {
        vleEntry["saturationPressure"] = nullptr;
        vleEntry["saturationPressureError"] = nullptr;
      }
      vleEntry["vaporDensity"] = point->vaporDensity;
      vleEntry["vaporDensityError"] = vaporDensityError;
      vleEntry["liquidDensity"] = point->liquidDensity;
      vleEntry["liquidDensityError"] = liquidDensityError;
      vleEntry["valleyDepth"] = point->valleyDepth;
      vleEntry["saturationPressureNormalization"] =
          point->saturationPressureFromEmptyBox ? "empty-box"
                                                : (hasSaturationPressure ? "ideal-gas reference" : "unavailable");
      outputJson["output"]["reweighting"]["vaporLiquidEquilibrium"].push_back(vleEntry);
    }
    std::print(stream, "\n");
    std::print(stream, "    (* saturation pressure normalized approximately by an ideal-gas reference at the\n");
    std::print(stream, "       dilute end of the pressure range; the empty-box state was never sampled)\n");
    if (anyUnavailableSaturationPressure)
    {
      std::print(stream, "    (# no saturation pressure: the empty box was never sampled and the sampled molecule\n");
      std::print(stream, "       counts do not reach down to the vapor, or to the dilute end of the pressure range;\n");
      std::print(stream, "       the coexistence fugacity and densities do not depend on it)\n");
    }
    std::print(stream, "    note: near coexistence the direct and reweighted loadings in the self-consistency\n");
    std::print(stream, "          table may legitimately differ; a single replica stays in one (possibly\n");
    std::print(stream, "          metastable) phase while the reweighted average covers both phases\n\n");
    std::print(stream, "    coexistence written to output/vle_coexistence.reweighted_histogram.txt\n");
    std::print(stream, "    coexistence P(N) written to output/vle_distribution_{{T}}.reweighted_histogram.txt\n\n");
  }

  outputJson["output"]["reweighting"]["pooledSamples"] = totalNumberOfSamples;
  outputJson["output"]["reweighting"]["occupiedBins"] = numberOfBins;
  outputJson["output"]["reweighting"]["iterations"] = iterations;
  outputJson["output"]["reweighting"]["residual"] = residual;
  outputJson["output"]["reweighting"]["tolerance"] = whamTolerance;
  outputJson["output"]["reweighting"]["maximumIterations"] = whamMaximumIterations;
  outputJson["output"]["reweighting"]["converged"] = whamConverged;
  outputJson["output"]["reweighting"]["unconvergedBlocks"] = whamUnconvergedBlocks;
  nlohmann::json unconvergedBlockDetails = nlohmann::json::array();
  for (const UnconvergedBlock& item : whamUnconvergedBlockDetails)
  {
    unconvergedBlockDetails.push_back(
        {{"block", item.index}, {"residual", item.residual}, {"iterations", item.iterations}});
  }
  outputJson["output"]["reweighting"]["unconvergedBlockDetails"] = unconvergedBlockDetails;
  outputJson["output"]["reweighting"]["freeEnergies"] = freeEnergies;

  std::print(stream, "    reweighted isotherms written to wham/reweighted_isotherm_{{T}}.reweighted_histogram.txt\n");
  std::print(stream, "    free energies written to wham/reweighted_free_energies.reweighted_histogram.txt\n");
  std::print(stream, "    density of states written to wham/density_of_states.reweighted_histogram.txt\n");
  if (computeBET && front.framework.has_value())
  {
    std::print(stream, "    BET plots written to wham/bet_{{T}}.reweighted_histogram.txt\n");
  }
  std::print(stream, "\n\n");
  std::flush(stream);

  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  totalReweightingAnalysisTime += (t2 - t1);
}

void ReweightedHistogram::writeReplicaFinalReports(std::vector<RunningEnergy>& recomputed)
{
  for (std::size_t replicaId = 0; replicaId < systems.size(); ++replicaId)
  {
    System& system = systems[replicaId];
    std::ostream replicaStream(replicaStreams[replicaId].rdbuf());

    std::print(replicaStream, "\n");
    std::print(replicaStream, "========================================================================================================================\n");
    std::print(replicaStream, "                             Simulation finished!\n");
    std::print(replicaStream, "========================================================================================================================\n");
    std::print(replicaStream, "\n");

    std::string status_line = std::format("Final state after {} cycles\n", numberOfProductionCycles);
    std::print(replicaStream, "{}", system.writeProductionStatusReportMC(status_line));

    const RunningEnergy drift = system.runningEnergies - recomputed[replicaId];
    replicaStream << system.runningEnergies.printMCDiff(recomputed[replicaId]);
    std::print(replicaStream, "\n\n");

    std::print(replicaStream, "Monte-Carlo moves statistics\n");
    std::print(replicaStream, "========================================================================================================================\n\n");
    std::print(replicaStream, "{}", system.writeMCMoveStatistics());

    std::print(replicaStream, "Production run CPU timings of the MC moves of this replica\n");
    std::print(replicaStream, "========================================================================================================================\n\n");
    for (std::size_t componentId{0}; const Component& component : system.components)
    {
      std::print(replicaStream, "{}",
                 component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(componentId, component.name));
      ++componentId;
    }
    std::print(replicaStream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());
    std::print(replicaStream, "\n\n");

    std::print(replicaStream, "{}",
               system.averageEnergies.writeAveragesStatistics(system.hasExternalField, system.framework,
                                                              system.components));
    if (!(system.framework.has_value() && system.framework->rigid))
    {
      std::print(replicaStream, "{}", system.averagePressure.writeAveragesStatistics());
    }
    std::print(replicaStream, "{}",
               system.averageLoadings.writeAveragesStatistics(
                   system.components, system.frameworkMass(),
                   system.framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));
    std::flush(replicaStream);

    // final flush of the analysis-property files (cycle 0 bypasses the 'writeEvery' gate)
    writeReplicaAnalysisOutputs(system, replicaId, 0uz);

    // final flush of the time-evolution files: the total cycle count is rounded up to a multiple
    // of the writer's own 'writeEvery' so its gate passes and all collected samples are written
    if (system.propertyNumberOfMoleculesEvolution.has_value() &&
        system.propertyNumberOfMoleculesEvolution->writeEvery.value_or(0uz) > 0uz)
    {
      const std::size_t writeEvery = system.propertyNumberOfMoleculesEvolution->writeEvery.value();
      system.propertyNumberOfMoleculesEvolution->writeOutput(
          replicaId, ((absoluteCycleOffset + writeEvery - 1uz) / writeEvery) * writeEvery);
    }
    if (system.propertyVolumeEvolution.has_value() && system.propertyVolumeEvolution->writeEvery.value_or(0uz) > 0uz)
    {
      const std::size_t writeEvery = system.propertyVolumeEvolution->writeEvery.value();
      system.propertyVolumeEvolution->writeOutput(
          replicaId, ((absoluteCycleOffset + writeEvery - 1uz) / writeEvery) * writeEvery);
    }

    // per-replica json statistics
    replicaJsons[replicaId]["output"]["runningEnergies"] = system.runningEnergies.jsonMC();
    replicaJsons[replicaId]["output"]["recomputedEnergies"] = recomputed[replicaId].jsonMC();
    replicaJsons[replicaId]["output"]["drift"] = drift.jsonMC();
    replicaJsons[replicaId]["output"]["MCMoveStatistics"]["system"] = system.jsonMCMoveStatistics();
    replicaJsons[replicaId]["output"]["cpuTimings"]["system"] =
        system.mc_moves_cputime.jsonSystemMCMoveCPUTimeStatistics();
    for (const Component& component : system.components)
    {
      replicaJsons[replicaId]["output"]["cpuTimings"][component.name] =
          component.mc_moves_cputime.jsonComponentMCMoveCPUTimeStatistics();
    }
    replicaJsons[replicaId]["properties"]["averageEnergies"] =
        system.averageEnergies.jsonAveragesStatistics(system.hasExternalField, system.framework, system.components);
    replicaJsons[replicaId]["properties"]["averagePressure"] = system.averagePressure.jsonAveragesStatistics();
    replicaJsons[replicaId]["properties"]["averageLoadings"] = system.averageLoadings.jsonAveragesStatistics(
        system.components, system.frameworkMass(),
        system.framework.transform([](const Framework& f) { return f.numberOfUnitCells; }));

    std::ofstream json(replicaJsonFileNames[replicaId]);
    json << replicaJsons[replicaId].dump(4);
  }
}

void ReweightedHistogram::writeBinaryRestartFile(std::size_t cyclesCompleted) noexcept
{
  cyclesCompletedThisStage = cyclesCompleted;

  ::writeBinaryRestartFile(*this);
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ReweightedHistogram::Sample& s)
{
  archive << s.energy;
  archive << s.numberOfMolecules;
  archive << s.block;
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ReweightedHistogram::Sample& s)
{
  archive >> s.energy;
  archive >> s.numberOfMolecules;
  archive >> s.block;
  return archive;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ReweightedHistogram& rh)
{
  archive << rh.versionNumber;

  archive << rh.random;

  archive << rh.numberOfProductionCycles;
  archive << rh.numberOfPreInitializationCycles;
  archive << rh.numberOfInitializationCycles;
  archive << rh.numberOfEquilibrationCycles;

  archive << rh.printEvery;
  archive << rh.optimizeMCMovesEvery;
  archive << rh.rescaleWangLandauEvery;
  archive << rh.writeBinaryRestartEvery;

  archive << rh.numberOfBlocks;
  archive << rh.parallelTemperingSwapEvery;

  archive << rh.sampleReweightingEvery;
  archive << rh.reweightingTemperatures;
  archive << rh.reweightingPressureRange;
  archive << rh.reweightingNumberOfPressures;
  archive << rh.whamTolerance;
  archive << rh.whamMaximumIterations;
  archive << rh.betScoutMaximumCycles;

  archive << rh.simulationStage;
  archive << rh.cyclesCompletedThisStage;

  archive << rh.temperatures;
  archive << rh.pressures;
  archive << rh.numberOfTemperatures;
  archive << rh.numberOfPressures;
  archive << rh.numberOfReplicas;
  archive << rh.systems;
  archive << rh.randoms;

  archive << rh.reweightingSamples;

  archive << rh.stepsPerReplica;
  archive << rh.absoluteCycleOffset;

  archive << rh.swapSweeps;
  archive << rh.sweepsThisStage;
  archive << rh.swapAttempts;
  archive << rh.swapAccepted;
  archive << rh.swapAttemptsPerTemperaturePair;
  archive << rh.swapAcceptedPerTemperaturePair;
  archive << rh.swapAttemptsPerPressurePair;
  archive << rh.swapAcceptedPerPressurePair;

  archive << rh.totalPreInitializationSimulationTime;
  archive << rh.totalInitializationSimulationTime;
  archive << rh.totalEquilibrationSimulationTime;
  archive << rh.totalProductionSimulationTime;
  archive << rh.totalReweightingAnalysisTime;
  archive << rh.totalSimulationTime;

  // the BET auto-setup state (the scouts are not rerun on a restart)
  archive << rh.computeBET;
  archive << rh.autoExternalPressures;
  archive << rh.autoReweightingPressureRange;
  archive << rh.autoMacroStateMaximum;
  archive << rh.nitrogenBETPressurePlan;
  archive << rh.nitrogenBETFillingCeiling;

  // the JSON output of the run (its header, the pre-isotherm fit of the ladder and the initial conditions of the
  // replicas), which a resumed run must not rebuild from its mid-run state
  archive << rh.outputJson.dump();
  {
    std::vector<std::string> jsons;
    for (const nlohmann::json& json : rh.replicaJsons) jsons.push_back(json.dump());
    archive << jsons;
  }
  archive << rh.reweightingLogNormalization;

  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ReweightedHistogram& rh)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > rh.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'ReweightedHistogram' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> rh.random;

  archive >> rh.numberOfProductionCycles;
  archive >> rh.numberOfPreInitializationCycles;
  archive >> rh.numberOfInitializationCycles;
  archive >> rh.numberOfEquilibrationCycles;

  archive >> rh.printEvery;
  archive >> rh.optimizeMCMovesEvery;
  archive >> rh.rescaleWangLandauEvery;
  archive >> rh.writeBinaryRestartEvery;

  archive >> rh.numberOfBlocks;
  archive >> rh.parallelTemperingSwapEvery;

  archive >> rh.sampleReweightingEvery;
  archive >> rh.reweightingTemperatures;
  archive >> rh.reweightingPressureRange;
  archive >> rh.reweightingNumberOfPressures;
  archive >> rh.whamTolerance;
  archive >> rh.whamMaximumIterations;
  archive >> rh.betScoutMaximumCycles;

  archive >> rh.simulationStage;
  archive >> rh.cyclesCompletedThisStage;

  archive >> rh.temperatures;
  archive >> rh.pressures;
  archive >> rh.numberOfTemperatures;
  archive >> rh.numberOfPressures;
  archive >> rh.numberOfReplicas;
  archive >> rh.systems;
  archive >> rh.randoms;

  archive >> rh.reweightingSamples;

  archive >> rh.stepsPerReplica;
  archive >> rh.absoluteCycleOffset;

  archive >> rh.swapSweeps;
  archive >> rh.sweepsThisStage;
  archive >> rh.swapAttempts;
  archive >> rh.swapAccepted;
  archive >> rh.swapAttemptsPerTemperaturePair;
  archive >> rh.swapAcceptedPerTemperaturePair;
  archive >> rh.swapAttemptsPerPressurePair;
  archive >> rh.swapAcceptedPerPressurePair;

  archive >> rh.totalPreInitializationSimulationTime;
  archive >> rh.totalInitializationSimulationTime;
  archive >> rh.totalEquilibrationSimulationTime;
  archive >> rh.totalProductionSimulationTime;
  archive >> rh.totalReweightingAnalysisTime;
  archive >> rh.totalSimulationTime;

  archive >> rh.computeBET;
  archive >> rh.autoExternalPressures;
  archive >> rh.autoReweightingPressureRange;
  archive >> rh.autoMacroStateMaximum;
  archive >> rh.nitrogenBETPressurePlan;
  archive >> rh.nitrogenBETFillingCeiling;
  std::string json;
  archive >> json;
  rh.outputJson = nlohmann::json::parse(json);
  std::vector<std::string> jsons;
  archive >> jsons;
  rh.replicaJsons.clear();
  for (const std::string& replicaJson : jsons) rh.replicaJsons.push_back(nlohmann::json::parse(replicaJson));
  archive >> rh.reweightingLogNormalization;

  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error("ReweightedHistogram: error in binary restart\n");
  }

  return archive;
}
