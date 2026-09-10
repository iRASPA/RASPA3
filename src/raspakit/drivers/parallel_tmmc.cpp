module;

module parallel_tmmc;

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
import equation_of_states;
import energy_status;
import running_energy;
import atom;
import int3;
import double3;
import double3x3;
import property_lambda_probability_histogram;
import transition_matrix;
import mc_moves;
import mc_moves_move_types;
import mc_moves_widom;
import mc_moves_cputime;
import mc_moves_statistics;
import cbmc;
import cbmc_chain_data;
import json;
import isotherm_bet;

// One Widom test insertion per this many production MC steps, per walker, and only while the
// walker sits in the bottom fraction of the macrostate range. The macrostates at the bottom hold
// only a small share of the steps even under a flattening bias, so they need the maximum rate to
// gather enough insertions; above the anchor region the insertions buy nothing, because the mean
// Rosenbluth weight of a filling pore is carried by ever rarer lucky insertions (measured on
// ferrierite: a relative error of 0.4 at the empty end rising past 0.9 at the full end).
static constexpr std::size_t widomSampleEvery = 1uz;
static constexpr double widomSampleFractionOfRange = 1.0;

// A counted test-insertion mean replaces the collection-matrix increment. 35 percent of the mean
// is 0.3 natural logs, while a 1D-N collection matrix is off by some ten natural logs per
// increment over the rise and an (N, λ) chain is worse (the molecule increment is a product of
// ~20 Metropolis λ-hops; a starved hop emits another 5–12 nats). Direct GCMC on ferrierite:
// reaching 35 molecules takes 0.0164 Pa, increment 15.6, against 15.5–15.6 Widom and 4.3–5.4
// from the collection matrix.
//
// The relative-error cut used to be 0.35 and it *ended the block*, not just labelled a hole.
// That is the wrong test: "slightly noisy" is not "no estimate". The (N, λ) FER canary of
// 50 k production cycles had 384–511 insertions at N = 2, 3, 4 with relative error 0.35–0.39,
// then 512 at N = 5 with relative error 0.23. Three misses in a row closed the block at N = 1
// and the isotherm sat on the λ-chain again (BET 104 m²/g, n_m = 2.30 / cell). The cut below
// is "the mean is still a mean": above it the average is one lucky insertion (relative error
// → 1 as the pore fills) and the collection-matrix tail takes over. Inside the block every
// counted increment is applied, including 0.36–0.80, so a flicker around 0.35 cannot open a
// collection-matrix hole.
static constexpr std::size_t minimumWidomInsertions = 100uz;
static constexpr double maximumWidomRelativeError = 0.85;

// The anchored block ends where the test insertions stop being a usable mean for good, which
// takes this many consecutive macrostates to establish. A single noisy or uncounted macrostate
// does not end the block: the relative error is not monotone in N. Ending at the first flicker
// truncated the block at N = 1 on ferrierite (see above).
static constexpr std::size_t widomAnchorFailureRun = 3uz;

// The analysis-property writers (RDFs, density grid, histograms, molecule properties) gate
// themselves on their own 'writeEvery'; a cycle argument of 0 forces the write (used for the
// final flush at the end of the run). The walker id keys the output filenames, so every walker
// writes its own set of files.
static void writeWalkerAnalysisOutputs(System& system, std::size_t walkerId, std::size_t cycle)
{
  if (system.propertyConventionalRadialDistributionFunction.has_value())
  {
    system.propertyConventionalRadialDistributionFunction->writeOutput(
        system.forceField, walkerId, system.simulationBox.volume, system.totalNumberOfPseudoAtoms, cycle);
  }
  if (system.propertyRadialDistributionFunction.has_value())
  {
    system.propertyRadialDistributionFunction->writeOutput(system.forceField, walkerId, system.simulationBox.volume,
                                                           system.totalNumberOfPseudoAtoms, cycle);
  }
  if (system.propertyDensityGrid.has_value())
  {
    system.propertyDensityGrid->writeOutput(walkerId, system.simulationBox, system.forceField, system.framework,
                                            system.components, cycle);
  }
  if (system.averageEnergyHistogram.has_value())
  {
    system.averageEnergyHistogram->writeOutput(walkerId, cycle);
  }
  if (system.averageNumberOfMoleculesHistogram.has_value())
  {
    system.averageNumberOfMoleculesHistogram->writeOutput(walkerId, system.components, cycle);
  }
  if (system.propertyMoleculeProperties.has_value())
  {
    system.propertyMoleculeProperties->writeOutput(walkerId, system.components, cycle);
  }
}

namespace
{
// A window whose lower bound is past the physical packing can never be grown into. The whole
// budget is spent per molecule, single-threaded, before any sampling starts, so it also sets how
// long the run takes to notice that the window is impossible.
constexpr std::size_t maxGrowAttempts = 5000uz;

bool growToMacrostate(System& system, RandomNumber& rng, std::size_t target)
{
  const std::size_t componentId = 0uz;
  while (system.numberOfIntegerMoleculesPerComponent[componentId] < target)
  {
    std::optional<ChainGrowData> growData = std::nullopt;
    bool insideBlockedPocket{false};
    std::size_t attempts = 0uz;
    do
    {
      do
      {
        growData = CBMC::growMoleculeSwapInsertion(
            rng,
            CBMC::GrowContext{system.hasExternalField, system.forceField, system.simulationBox,
                              system.interpolationGrids, system.externalFieldInterpolationGrid, system.framework,
                              system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), system.beta,
                              system.forceField.cutOffFrameworkVDW, system.forceField.cutOffMoleculeVDW,
                              system.forceField.cutOffCoulomb},
            system.components[componentId], componentId, system.numberOfMolecules(), 1.0, false, false);
        ++attempts;
        if (attempts >= maxGrowAttempts) return false;
      } while (!growData || growData->energies.potentialEnergy() > system.forceField.energyOverlapCriteria);

      std::span<const Atom> newMolecule = std::span(growData->atoms.begin(), growData->atoms.end());
      insideBlockedPocket = system.insideBlockedPockets(system.components[componentId], newMolecule);
    } while (insideBlockedPocket);

    system.insertMolecule(componentId, growData->molecule, growData->atoms);
  }
  return true;
}

// The collection matrix is built from acceptance probabilities of attempted moves, not from visit
// frequencies, so it is an unbiased estimator of the transition probabilities whatever the
// configuration it was sampled from and is kept across all stages. Dropping the initialization
// entries strands the walker: filling a micropore from empty is the only part of the run that
// records transitions at low N, and once they are gone the bias gives no reason to go back down,
// while an unbiased deletion out of a filled micropore at 77 K is of order e^-40.
void pinWalkerWindow(System& system, double temperature, std::size_t windowIndex, std::size_t minN, std::size_t maxN)
{
  system.tmmc.minMacrostate = minN;
  system.tmmc.maxMacrostate = maxN;
  system.tmmc.rezeroAfterInitialization = false;
  system.tmmc.statisticsFileName = std::format("tmmc/tmmc_statistics_{}_w{}.parallel_tmmc.txt", temperature, windowIndex);
  system.tmmc.initialize();
}

// Neighboring partition segments [boundaries[w], boundaries[w+1]] are expanded by `overlap` so
// the walkers share several interior macrostates, not only the bound state they cannot cross.
struct WindowRange
{
  std::size_t lo;
  std::size_t hi;
};

WindowRange overlappedWindow(std::size_t windowIndex, std::size_t numberOfWindows, std::size_t minMacrostate,
                             std::size_t maxMacrostate, const std::vector<std::size_t>& windowBoundaries)
{
  WindowRange range{windowBoundaries[windowIndex], windowBoundaries[windowIndex + 1uz]};
  if (numberOfWindows <= 1uz) return range;
  const std::size_t span = maxMacrostate - minMacrostate;
  const std::size_t slice = std::max(1uz, span / numberOfWindows);
  const std::size_t overlap = std::max(4uz, slice / 3uz);
  if (windowIndex > 0uz)
  {
    range.lo -= std::min(overlap, range.lo - minMacrostate);
  }
  if (windowIndex + 1uz < numberOfWindows)
  {
    range.hi = std::min(maxMacrostate, range.hi + overlap);
  }
  return range;
}
}  // namespace

ParallelTMMC::ParallelTMMC(InputReader& reader)
    : random(reader.randomSeed),
      numberOfProductionCycles(reader.numberOfProductionCycles),
      numberOfPreInitializationCycles(reader.numberOfPreInitializationCycles),
      numberOfInitializationCycles(reader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(reader.numberOfEquilibrationCycles),
      printEvery(reader.printEvery),
      optimizeMCMovesEvery(reader.optimizeMCMovesEvery),
      rescaleWangLandauEvery(reader.rescaleWangLandauEvery),
      writeBinaryRestartEvery(reader.writeBinaryRestartEvery),
      numberOfBlocks(reader.numberOfBlocks),
      reweightingNumberOfPressures(reader.reweightingNumberOfPressures),
      computeBET(reader.computeBET),
      autoReweightingPressureRange(reader.autoReweightingPressureRange),
      autoMacroStateMaximum(reader.autoMacroStateMaximum),
      temperatures(reader.parallelTemperingTemperatures),
      numberOfWindows(reader.tmmcNumberOfWindows)
{
  // a single 'ExternalTemperature' gives a one-temperature run
  if (temperatures.empty())
  {
    temperatures.push_back(reader.systems.front().temperature);
  }
  numberOfTemperatures = temperatures.size();
  numberOfWalkers = numberOfTemperatures * numberOfWindows;

  referencePressure = reader.systems.front().input_pressure;

  System templateSystem = std::move(reader.systems.front());
  reader.systems.clear();

  if (autoReweightingPressureRange)
  {
    nitrogenBETPressurePlan =
        planNitrogenBETPressures(templateSystem, temperatures.front(), reader.numberOfThreads);
    reweightingPressureRange =
        std::make_pair(nitrogenBETPressurePlan->lowestPressure, nitrogenBETPressurePlan->highestPressure);
    if (!reader.reweightingNumberOfPressuresSpecified)
    {
      reweightingNumberOfPressures = nitrogenBETPressurePlan->tmmcReweightingNumberOfPressures;
    }
  }
  else
  {
    reweightingPressureRange =
        reader.reweightingPressureRange.value_or(std::make_pair(0.01 * referencePressure, 100.0 * referencePressure));
  }

  if (autoMacroStateMaximum)
  {
    nitrogenBETFillingCeiling = scoutNitrogenBETFillingCeiling(templateSystem);
    templateSystem.tmmc.maxMacrostate =
        std::max(templateSystem.tmmc.minMacrostate + 1uz, nitrogenBETFillingCeiling->maxMacrostate);
  }
  if (templateSystem.tmmc.maxMacrostate <= templateSystem.tmmc.minMacrostate)
  {
    throw std::runtime_error(
        "[ParallelTMMC]: the template system needs a macrostate range with minimum < maximum\n");
  }
  const std::size_t macrostateSpan = templateSystem.tmmc.maxMacrostate - templateSystem.tmmc.minMacrostate;
  if (numberOfWindows > macrostateSpan)
  {
    numberOfWindows = std::max(1uz, macrostateSpan);
  }
  numberOfWalkers = numberOfTemperatures * numberOfWindows;

  initializeWalkers(std::move(templateSystem));
  numberOfStepsPerCycle = std::max(1uz, maxMacrostate);
}

ParallelTMMC::ParallelTMMC(System templateSystem, std::vector<double> temperatures_,
                           ParallelTMMCParameters parameters)
    : random(parameters.randomSeed),
      numberOfProductionCycles(parameters.numberOfProductionCycles),
      numberOfPreInitializationCycles(parameters.numberOfPreInitializationCycles),
      numberOfInitializationCycles(parameters.numberOfInitializationCycles),
      numberOfEquilibrationCycles(parameters.numberOfEquilibrationCycles),
      printEvery(parameters.printEvery),
      optimizeMCMovesEvery(parameters.optimizeMCMovesEvery),
      rescaleWangLandauEvery(parameters.rescaleWangLandauEvery),
      writeBinaryRestartEvery(parameters.writeBinaryRestartEvery),
      numberOfBlocks(parameters.numberOfBlocks),
      reweightingPressureRange(parameters.reweightingPressureRange),
      reweightingNumberOfPressures(parameters.reweightingNumberOfPressures),
      computeBET(parameters.computeBET),
      temperatures(std::move(temperatures_)),
      numberOfWindows(std::max(1uz, parameters.numberOfWindows))
{
  if (temperatures.empty())
  {
    temperatures.push_back(templateSystem.temperature);
  }
  numberOfTemperatures = temperatures.size();

  referencePressure = templateSystem.input_pressure;

  templateSystem.tmmc.doTMMC = true;
  templateSystem.tmmc.useBias = true;
  templateSystem.tmmc.useTMBias = true;
  templateSystem.tmmc.rejectOutOfBound = true;
  templateSystem.tmmc.useWangLandau = true;
  templateSystem.tmmc.updateTMEvery = std::max(1uz, parameters.tmmcUpdateEvery);

  if (templateSystem.tmmc.maxMacrostate <= templateSystem.tmmc.minMacrostate)
  {
    throw std::runtime_error(
        "[ParallelTMMC]: the template system needs a macrostate range with minimum < maximum\n");
  }
  const std::size_t macrostateSpan = templateSystem.tmmc.maxMacrostate - templateSystem.tmmc.minMacrostate;
  if (numberOfWindows > macrostateSpan)
  {
    numberOfWindows = std::max(1uz, macrostateSpan);
  }
  numberOfWalkers = numberOfTemperatures * numberOfWindows;

  initializeWalkers(std::move(templateSystem));
  numberOfStepsPerCycle = std::max(1uz, maxMacrostate);
}

void ParallelTMMC::initializeWalkers(System templateSystem)
{
  if (templateSystem.components.size() == 1uz)
  {
    const Component& component = templateSystem.components.front();
    const bool usesCFCMCSwap =
        component.mc_moves_probabilities.getProbability(Move::Types::SwapCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::SwapCBCFCMC) > 0.0;
    const std::size_t componentLambdaBins = component.lambdaGC.numberOfSamplePoints;
    if (usesCFCMCSwap && componentLambdaBins > 1uz)
    {
      if (templateSystem.tmmc.numberOfLambdaBins != 1uz &&
          templateSystem.tmmc.numberOfLambdaBins != componentLambdaBins)
      {
        throw std::runtime_error(std::format(
            "[ParallelTMMC]: transition-matrix lambda bins ({}) do not match the CFCMC component lambda bins ({})\n",
            templateSystem.tmmc.numberOfLambdaBins, componentLambdaBins));
      }
      templateSystem.tmmc.numberOfLambdaBins = componentLambdaBins;
    }
  }

  numberOfWalkers = numberOfTemperatures * numberOfWindows;

  // the macrostate windows: windowBoundaries is a partition of [min, max]; each walker is then
  // expanded so neighboring windows overlap by several interior states (needed to gauge-match ln Π)
  minMacrostate = templateSystem.tmmc.minMacrostate;
  maxMacrostate = templateSystem.tmmc.maxMacrostate;
  windowBoundaries.resize(numberOfWindows + 1uz);
  for (std::size_t windowIndex = 0; windowIndex <= numberOfWindows; ++windowIndex)
  {
    windowBoundaries[windowIndex] = minMacrostate + (windowIndex * (maxMacrostate - minMacrostate)) / numberOfWindows;
  }

  // the single declared system is replicated into one walker per (temperature, window) pair
  systems.reserve(numberOfWalkers);
  for (std::size_t walkerId = 0; walkerId + 1 < numberOfWalkers; ++walkerId)
  {
    systems.push_back(templateSystem);
  }
  systems.push_back(std::move(templateSystem));

  // walker (t, w) is pinned at temperature T_t and macrostate window w, with its own
  // random-number stream and its own transition-matrix statistics file
  randoms.reserve(numberOfWalkers);
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
    {
      const std::size_t walkerId = walkerIndex(temperatureIndex, windowIndex);
      System& system = systems[walkerId];
      const double T = temperatures[temperatureIndex];

      system.temperature = T;
      system.beta = 1.0 / (Units::KB * T);

      if (system.forceField.temperature != T)
      {
        system.forceField.temperature = T;
        system.forceField.preComputeDerivedParameters();
        system.forceField.preComputePotentialShift();
        system.forceField.preComputeTailCorrection();
      }

      // convert the reference pressure to the per-temperature reference fugacity: the fugacity
      // coefficient is recomputed with the Peng-Robinson equation of state at (T_t, P_ref)
      // (an explicitly given 'FugacityCoefficient' would only be valid at one temperature)
      for (Component& component : system.components)
      {
        component.fugacityCoefficient = std::nullopt;
      }
      system.equationOfState =
          EquationOfState(EquationOfState::Type::PengRobinson, EquationOfState::MixingRules::VanDerWaals, T,
                          referencePressure, system.simulationBox, system.heliumVoidFraction, system.components);

      // the CBMC ideal-gas conformation reservoirs are Boltzmann samples at the system temperature
      system.buildConformationReservoirs();

      // the walker is confined to its (overlapped) window; the collection matrix is kept across
      // the stages, see pinWalkerWindow
      const WindowRange range =
          overlappedWindow(windowIndex, numberOfWindows, minMacrostate, maxMacrostate, windowBoundaries);
      system.tmmc.minMacrostate = range.lo;
      system.tmmc.maxMacrostate = range.hi;
      system.tmmc.rezeroAfterInitialization = false;
      system.tmmc.statisticsFileName = std::format("tmmc/tmmc_statistics_{}_w{}.parallel_tmmc.txt", T, windowIndex);
      system.tmmc.initialize();

      randoms.emplace_back(random.seed + walkerId + 1);
    }
  }

  blockCollectionMatrices.resize(numberOfWalkers);
  stepsPerWalker.assign(numberOfWalkers, 0uz);

  // The Widom statistics live on the global macrostate grid, so the windows of one temperature
  // simply add up (like the collection matrices). Grown rather than assigned: on a resume from a
  // binary restart the deserialized statistics must survive this call.
  const std::size_t numberOfMacrostates = maxMacrostate - minMacrostate + 1uz;
  widomWeightSums.resize(numberOfWalkers);
  widomWeightSquaredSums.resize(numberOfWalkers);
  widomInsertions.resize(numberOfWalkers);
  for (std::size_t walkerId = 0; walkerId < numberOfWalkers; ++walkerId)
  {
    widomWeightSums[walkerId].resize(numberOfMacrostates, 0.0);
    widomWeightSquaredSums[walkerId].resize(numberOfMacrostates, 0.0);
    widomInsertions[walkerId].resize(numberOfMacrostates, 0uz);
  }
}

void ParallelTMMC::run()
{
  setup();
  runStage(SimulationStage::PreInitialization, numberOfPreInitializationCycles);
  runStage(SimulationStage::Initialization, numberOfInitializationCycles);
  runStage(SimulationStage::Equilibration, numberOfEquilibrationCycles);
  runStage(SimulationStage::Production, numberOfProductionCycles);
  output();
}

void ParallelTMMC::setup()
{
  numberOfStepsPerCycle = std::max(1uz, maxMacrostate);

  for (System& system : systems)
  {
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);
  }

  std::filesystem::create_directories("output");
  std::filesystem::create_directories("tmmc");

  // on a binary-restart resume append to the existing output files (and skip re-printing the
  // headers) so each log continues where the interrupted run left off
  const bool resumedFromBinaryRestart = simulationStage != SimulationStage::Uninitialized;
  stream.open("output/output.parallel_tmmc.txt", resumedFromBinaryRestart ? std::ios::app : std::ios::out);
  outputJsonFileName = "output/output.parallel_tmmc.json";

  if (!resumedFromBinaryRestart)
  {
    std::print(stream, "{}", systems.front().writeOutputHeader());
    std::print(stream, "Random seed: {}\n\n", random.seed);
    std::print(stream, "{}\n", HardwareInfo::writeInfo());
    std::print(stream, "{}", Units::printStatus());
  }

  // interpolation grids are computed once and shared (copied) between the walkers
  systems.front().createExternalFieldInterpolationGrid(stream, 0);
  systems.front().createFrameworkInterpolationGrids(stream);
  for (std::size_t walkerId = 1; walkerId < systems.size(); ++walkerId)
  {
    systems[walkerId].externalFieldInterpolationGrid = systems.front().externalFieldInterpolationGrid;
    systems[walkerId].interpolationGrids = systems.front().interpolationGrids;
  }

  // Grow windows in increasing-N order per temperature: copy the previous walker's configuration
  // and insert the extra molecules. Independent growth from empty into a high window hangs when
  // that window's lower bound exceeds the physical packing. If packing stops short of Nmax, the
  // remaining windows are dropped and the last used window is closed at the packed N.
  if (!resumedFromBinaryRestart)
  {
    std::print(stream, "Growing the initial configurations into their windows (CBMC)\n");
    std::flush(stream);

    const std::size_t layoutWindows = numberOfWindows;
    for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
    {
      for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
      {
        const std::size_t walkerId = temperatureIndex * layoutWindows + windowIndex;
        if (windowIndex > 0)
        {
          const std::size_t previousId = temperatureIndex * layoutWindows + windowIndex - 1uz;
          systems[walkerId] = systems[previousId];
          const WindowRange range =
              overlappedWindow(windowIndex, layoutWindows, minMacrostate, maxMacrostate, windowBoundaries);
          pinWalkerWindow(systems[walkerId], temperatures[temperatureIndex], windowIndex, range.lo, range.hi);
        }
        const WindowRange growRange =
            overlappedWindow(windowIndex, layoutWindows, minMacrostate, maxMacrostate, windowBoundaries);
        if (growToMacrostate(systems[walkerId], randoms[walkerId], growRange.lo)) continue;

        const std::size_t packed = systems[walkerId].numberOfIntegerMoleculesPerComponent[0];
        if (windowIndex == 0uz || packed <= minMacrostate)
        {
          throw std::runtime_error(
              "[ParallelTMMC]: CBMC could not place a molecule in the empty box; the macrostate range is empty\n");
        }

        maxMacrostate = packed;
        numberOfWindows = windowIndex;
        windowBoundaries.resize(numberOfWindows + 1uz);
        windowBoundaries[numberOfWindows] = maxMacrostate;
        const std::size_t lastId = temperatureIndex * layoutWindows + windowIndex - 1uz;
        pinWalkerWindow(systems[lastId], temperatures[temperatureIndex], windowIndex - 1uz,
                        systems[lastId].tmmc.minMacrostate, maxMacrostate);
        std::print(stream, "Physical packing is {} molecules; dropping windows above that bound\n", packed);
        break;
      }
    }

    if (numberOfWindows < layoutWindows)
    {
      std::vector<System> compacted;
      compacted.reserve(numberOfTemperatures * numberOfWindows);
      std::vector<RandomNumber> compactedRandoms;
      compactedRandoms.reserve(numberOfTemperatures * numberOfWindows);
      for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
      {
        for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
        {
          const std::size_t walkerId = temperatureIndex * layoutWindows + windowIndex;
          compacted.push_back(std::move(systems[walkerId]));
          compactedRandoms.push_back(std::move(randoms[walkerId]));
        }
      }
      systems = std::move(compacted);
      randoms = std::move(compactedRandoms);
      numberOfWalkers = systems.size();
      blockCollectionMatrices.resize(numberOfWalkers);
      stepsPerWalker.assign(numberOfWalkers, 0uz);
    }

    std::print(stream, "Initial configurations ready\n\n");
    std::flush(stream);
  }

  if (!resumedFromBinaryRestart)
  {
    std::print(stream, "Parallel transition-matrix Monte Carlo (TMMC)\n");
    std::print(stream, "===============================================================================\n\n");
    std::print(stream, "Number of temperatures:                      {}\n", numberOfTemperatures);
    std::print(stream, "Number of macrostate windows:                {}\n", numberOfWindows);
    std::print(stream, "Number of walkers / threads:                 {}\n", numberOfWalkers);
    std::print(stream, "MC steps per cycle:                          {} (global macrostate ceiling N_max)\n",
               numberOfStepsPerCycle);
    std::print(stream, "Temperature ladder:                         ");
    for (double T : temperatures)
    {
      std::print(stream, " {}", T);
    }
    std::print(stream, " [K]\n");
    std::print(stream, "Reference pressure:                          {:.5e} [Pa]\n", referencePressure);
    std::print(stream, "Macrostate range:                            [{}, {}] molecules\n", minMacrostate, maxMacrostate);
    if (systems.front().tmmc.lambdaChain())
    {
      std::print(stream,
                 "Lambda bins per molecule count:              {} (macrostate is (N, λ); isotherm uses λ = 0)\n",
                 systems.front().tmmc.lambdaBinCount());
    }
    std::print(stream, "Bias update every:                           {} steps\n", systems.front().tmmc.updateTMEvery);
    std::print(stream, "Isotherm/coexistence pressure scan:          {:.5e} - {:.5e} [Pa], {} log-spaced points\n\n",
               reweightingPressureRange.first, reweightingPressureRange.second, reweightingNumberOfPressures);
    if (nitrogenBETPressurePlan.has_value())
    {
      writeNitrogenBETPressurePlan(stream, *nitrogenBETPressurePlan);
    }
    if (nitrogenBETFillingCeiling.has_value())
    {
      writeNitrogenBETFillingCeiling(stream, *nitrogenBETFillingCeiling);
    }

    std::print(stream, "Walker grid: walker (t, w) = t * {} + w\n", numberOfWindows);
    std::print(stream, "    walker    temperature [K]    window [molecules]\n");
    std::print(stream, "    ----------------------------------------------\n");
    for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
    {
      for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
      {
        const std::size_t walkerId = walkerIndex(temperatureIndex, windowIndex);
        std::print(stream, "    {:6d}    {:15.4f}    [{}, {}]\n", walkerId, temperatures[temperatureIndex],
                   systems[walkerId].tmmc.minMacrostate, systems[walkerId].tmmc.maxMacrostate);
      }
    }
    std::print(stream, "\n");
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
  outputJson["initialization"]["referencePressure"] = referencePressure;
  outputJson["initialization"]["macrostateRange"] = std::vector<std::size_t>{minMacrostate, maxMacrostate};
  outputJson["initialization"]["numberOfStepsPerCycle"] = numberOfStepsPerCycle;
  outputJson["initialization"]["numberOfLambdaBins"] = systems.front().tmmc.lambdaBinCount();
  outputJson["initialization"]["windowBoundaries"] = windowBoundaries;
  outputJson["initialization"]["reweightingPressureRange"] =
      std::vector<double>{reweightingPressureRange.first, reweightingPressureRange.second};
  outputJson["initialization"]["reweightingNumberOfPressures"] = reweightingNumberOfPressures;
  if (nitrogenBETPressurePlan.has_value())
  {
    outputJson["initialization"]["nitrogenBETHenryCoefficient"] = nitrogenBETPressurePlan->henryCoefficientPerCell;
    outputJson["initialization"]["autoReweightingPressureRange"] = autoReweightingPressureRange;
  }
  if (nitrogenBETFillingCeiling.has_value())
  {
    outputJson["initialization"]["autoMacroStateMaximum"] = autoMacroStateMaximum;
    outputJson["initialization"]["fillingCeiling"] = nitrogenBETFillingCeiling->maxMacrostate;
    outputJson["initialization"]["fillingCeilingMeanOccupancy"] = nitrogenBETFillingCeiling->meanOccupancy;
  }

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);

  // per-walker output files: each worker thread writes exclusively to its own stream
  walkerStreams.reserve(numberOfWalkers);
  walkerJsonFileNames.reserve(numberOfWalkers);
  walkerJsons.resize(numberOfWalkers);
  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
    {
      const std::size_t walkerId = walkerIndex(temperatureIndex, windowIndex);
      const System& system = systems[walkerId];
      walkerStreams.emplace_back(
          std::format("output/output_{}_w{}.parallel_tmmc.r{}.txt", system.temperature, windowIndex, walkerId),
          resumedFromBinaryRestart ? std::ios::app : std::ios::out);
      walkerJsonFileNames.emplace_back(
          std::format("output/output_{}_w{}.parallel_tmmc.r{}.json", system.temperature, windowIndex, walkerId));

    if (!resumedFromBinaryRestart)
    {
        std::ostream walkerStream(walkerStreams[walkerId].rdbuf());
        std::print(walkerStream, "{}", system.writeOutputHeader());
        std::print(walkerStream, "Parallel TMMC: walker {} of {} (temperature {} [K], window [{}, {}] molecules)\n",
                   walkerId, numberOfWalkers, system.temperature, system.tmmc.minMacrostate, system.tmmc.maxMacrostate);
        std::print(walkerStream, "Random seed of this walker: {}\n\n", randoms[walkerId].seed);
        std::print(walkerStream, "{}\n", HardwareInfo::writeInfo());
        std::print(walkerStream, "{}", Units::printStatus());
        std::print(walkerStream, "{}", system.writeSystemStatus());
        std::print(walkerStream, "{}", system.forceField.printPseudoAtomStatus());
        std::print(walkerStream, "{}", system.forceField.printForceFieldStatus());
        std::print(walkerStream, "{}", system.writeComponentStatus());
        std::print(walkerStream, "{}", system.writeNumberOfPseudoAtoms());
    }

#ifdef VERSION
      walkerJsons[walkerId]["version"] = EXPAND_AND_QUOTE(VERSION);
#endif
      walkerJsons[walkerId]["seed"] = randoms[walkerId].seed;
      walkerJsons[walkerId]["walkerId"] = walkerId;
      walkerJsons[walkerId]["temperature"] = system.temperature;
      walkerJsons[walkerId]["window"] = std::vector<std::size_t>{system.tmmc.minMacrostate, system.tmmc.maxMacrostate};
      walkerJsons[walkerId]["initialization"]["initialConditions"] = system.jsonSystemStatus();
      walkerJsons[walkerId]["initialization"]["components"] = system.jsonComponentStatus();

      std::ofstream walkerJson(walkerJsonFileNames[walkerId]);
      walkerJson << walkerJsons[walkerId].dump(4);
    }
  }
}

void ParallelTMMC::performWalkerCycle(std::size_t walkerId, SimulationStage stage, std::size_t currentBlock)
{
  System& system = systems[walkerId];
  RandomNumber& rng = randoms[walkerId];

  // every walker is self-contained; the Gibbs-style moves that need a partner system are not
  // supported by this driver
  std::size_t fractionalMoleculeSystem = 0uz;

  // global maxMacrostate, not the current loading or the window width: every walker does the
  // same number of moves per cycle so cycle-based sampling compares equal amounts of work
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

        // Wang-Landau biasing of the CFCMC lambda moves. Off when TMMC already flattens the
        // (N, λ) chain: that Wang-Landau would double-bias λ and is replaced by tmmc.visitWangLandau.
        if (!system.tmmc.lambdaChain())
        {
          system.components[selectedComponent].lambdaGC.WangLandauIteration(
              PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample,
              system.lambdaWangLandauIsActive(selectedComponent));
          system.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
          system.reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        }
        break;
      case SimulationStage::Production:
        MC_Moves::performRandomMoveProduction(rng, system, system, selectedComponent, fractionalMoleculeSystem,
                                              currentBlock);
        ++stepsPerWalker[walkerId];

        // Widom test insertion at the current macrostate: the mean Rosenbluth weight gives the
        // exact ln Pi increment out of this macrostate (see the class documentation). The move
        // leaves the configuration untouched, so it does not disturb the walk. On an (N, λ)
        // chain the physical increment is a full-molecule insertion at λ = 0, so only sample there.
        if (stepsPerWalker[walkerId] % widomSampleEvery == 0uz)
        {
          const std::size_t currentN = system.numberOfIntegerMoleculesPerComponent[0];
          const bool atDecoupledLambda =
              !system.tmmc.lambdaChain() ||
              (!system.components.empty() && system.components.front().lambdaGC.currentBin == 0uz);
          const std::size_t sampledMacrostates = std::max(
              1uz, static_cast<std::size_t>(widomSampleFractionOfRange *
                                            static_cast<double>(widomWeightSums[walkerId].size())));
          if (atDecoupledLambda && currentN >= minMacrostate && currentN - minMacrostate < sampledMacrostates)
          {
            const double weight = MC_Moves::WidomMove(rng, system, 0uz);
            const std::size_t index = currentN - minMacrostate;
            widomWeightSums[walkerId][index] += weight;
            widomWeightSquaredSums[walkerId][index] += weight * weight;
            ++widomInsertions[walkerId][index];
          }
        }
        break;
    }

    // the TMMC state sampling: the moves already recorded the unbiased acceptance probabilities
    // into the collection matrix; here the visit histogram is updated and the flattening bias is
    // re-derived from the collection matrix every 'TMMCUpdateEvery' steps.
    //
    // The bias is already built during the initialization stage. That stage is where a walker
    // started from an empty box crosses every macrostate on its way to its equilibrium loading;
    // biasing it there is what turns that one-way filling into a random walk over the window.
    // Left unbiased, the walker arrives at the filling loading and stays: the equilibration stage
    // then only ever sees the top of the range.
    const std::size_t currentN = system.numberOfIntegerMoleculesPerComponent[0];
    const std::size_t lambdaBin =
        system.tmmc.lambdaChain() && !system.components.empty() ? system.components.front().lambdaGC.currentBin : 0uz;
    system.tmmc.updateHistogram(currentN, lambdaBin);
    system.tmmc.numberOfSteps++;
    if (stage == SimulationStage::Initialization || stage == SimulationStage::Equilibration ||
        stage == SimulationStage::Production)
    {
      // While Wang-Landau owns the bias (initialization, equilibration, first quarter of
      // production) it drives the walk and adjustBias only keeps ln Pi current for the
      // statistics files. After switchToTMBias (rest of production) visitWangLandau is inert
      // and adjustBias re-derives the flattening bias from the growing collection matrix
      // every 'tmmcUpdateEvery' steps. Either way the estimate is safe: the collection matrix
      // is built from the unbiased acceptance probabilities, so whatever bias makes the walk
      // ergodic, ln Pi stands.
      system.tmmc.visitWangLandau(currentN, lambdaBin);
      system.tmmc.adjustBias();
    }

    system.components[selectedComponent].lambdaGC.sampleOccupancy(system.containsTheFractionalMolecule);
    system.pairSwapLambdaSampleOccupancy();
    system.reactionLambdaSampleOccupancy();
  }
}

std::vector<double3> ParallelTMMC::productionCollectionMatrix(std::size_t walkerId) const
{
  if (walkerId >= systems.size() || walkerId >= productionStartCollectionMatrices.size())
  {
    throw std::runtime_error("[ParallelTMMC]: production collection-matrix snapshot is missing");
  }
  std::vector<double3> matrix = systems[walkerId].tmmc.cmatrix;
  const std::vector<double3>& start = productionStartCollectionMatrices[walkerId];
  if (start.size() != matrix.size())
  {
    throw std::runtime_error("[ParallelTMMC]: production collection-matrix snapshot has the wrong size");
  }
  for (std::size_t index = 0; index < matrix.size(); ++index)
  {
    matrix[index] -= start[index];
  }
  return matrix;
}

void ParallelTMMC::runStage(SimulationStage stage, std::size_t numberOfCycles)
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
      // drop the collection-matrix statistics of the relaxing initial configurations
      system.tmmc.clearCMatrix();

      if (!system.tmmc.lambdaChain())
      {
        for (Component& component : system.components)
        {
          component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize,
                                                 system.containsTheFractionalMolecule);
          component.lambdaGC.clear();
        }
        system.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize);
        system.pairSwapLambdaClearBookkeeping();
        system.reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize);
        system.reactionLambdaClearBookkeeping();
      }
    }
  }
  if (startCycle == 0uz && stage == SimulationStage::Production)
  {
    productionStartCollectionMatrices.resize(systems.size());
    productionStartHistograms.resize(systems.size());
    for (std::size_t walkerId = 0; walkerId < systems.size(); ++walkerId)
    {
      System& system = systems[walkerId];

      // Drop all pre-production statistics but keep the Wang-Landau bias built so far. The
      // recorded acceptance probabilities are formally unbiased under any applied bias, but
      // they are conditional averages over the configurations actually visited: during the
      // Wang-Landau exploration the molecules at low N have not yet relaxed into the deep
      // adsorption sites, their deletion acceptances are overestimated by orders of magnitude,
      // and because a row estimate is a heavy-tailed mean those early entries dominate it no
      // matter how long production runs. Production therefore restarts the statistics from
      // zero: Wang-Landau keeps driving the walk for the first quarter of production (its
      // descents sample low N with the well-bound survivor configurations), after which the
      // bias switches to the transition-matrix bias -ln Pi derived from those clean rows,
      // flattening the walk so every macrostate collects comparable statistics (Errington/Shen
      // hybrid; the switch happens inside the walker threads at wangLandauProductionCycles).
      system.tmmc.clearStatisticsKeepBias();
      productionStartCollectionMatrices[walkerId] = system.tmmc.cmatrix;
      productionStartHistograms[walkerId] = system.tmmc.histogram;

      // the Widom anchor is a production-stage measurement as well
      std::fill(widomWeightSums[walkerId].begin(), widomWeightSums[walkerId].end(), 0.0);
      std::fill(widomWeightSquaredSums[walkerId].begin(), widomWeightSquaredSums[walkerId].end(), 0.0);
      std::fill(widomInsertions[walkerId].begin(), widomInsertions[walkerId].end(), 0uz);

      system.mc_moves_statistics.clearMoveStatistics();
      system.mc_moves_cputime.clearTimingStatistics();

      for (Component& component : system.components)
      {
        component.mc_moves_statistics.clearMoveStatistics();
        component.mc_moves_cputime.clearTimingStatistics();

        if (!system.tmmc.lambdaChain())
        {
          component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize,
                                                 system.containsTheFractionalMolecule);
          component.lambdaGC.clear();
        }
      }
      if (!system.tmmc.lambdaChain())
      {
        system.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize);
        system.pairSwapLambdaClearBookkeeping();
        system.reactionLambdaFinalize();
        system.reactionLambdaClearBookkeeping();
      }
    }
    std::fill(stepsPerWalker.begin(), stepsPerWalker.end(), 0uz);
  }

  {
    std::scoped_lock lock(outputMutex);
    const std::string_view stageName = (stage == SimulationStage::PreInitialization) ? "Pre-initialization"
                                       : (stage == SimulationStage::Initialization)  ? "Initialization"
                                       : (stage == SimulationStage::Equilibration)   ? "Equilibration"
                                                                                     : "Production";
    if (startCycle == 0uz)
    {
      std::print(stream, "\n{} stage: {} cycles on {} walkers/threads\n", stageName, numberOfCycles, numberOfWalkers);
    }
    else
    {
      std::print(stream, "\n{} stage: resumed from binary restart at cycle {} of {} on {} walkers/threads\n",
                 stageName, startCycle, numberOfCycles, numberOfWalkers);
    }
    std::flush(stream);
  }

  const std::size_t stageCycleOffset = absoluteCycleOffset;

  // the walkers are fully independent within a stage; the barrier is only used for the periodic
  // binary-restart checkpoint (every walker runs the same cycle count, so they all agree on when
  // a checkpoint is due). The completion runs while all worker threads are parked, so the full
  // driver state is consistent. With 'writeBinaryRestartEvery' disabled the barrier is never used.
  auto onAllArrived = [this]() noexcept
  {
    writeBinaryRestartFile(checkpointCycle.load(std::memory_order_relaxed));
    // all worker threads are parked on this barrier, so the checkpoint just written is a
    // consistent snapshot: safe to exit here on a shutdown signal
    if (GracefulShutdown::requested())
    {
      // std::exit skips stack unwinding: flush the text output streams explicitly
      std::flush(stream);
      for (std::ofstream& walkerStream : walkerStreams) std::flush(walkerStream);
      GracefulShutdown::exitAfterCheckpoint();
    }
  };
  std::barrier synchronizationPoint(static_cast<std::ptrdiff_t>(numberOfWalkers), onAllArrived);

  // one worker thread per walker
  {
    std::vector<std::jthread> threads;
    threads.reserve(numberOfWalkers);
    for (std::size_t walkerId = 0; walkerId < numberOfWalkers; ++walkerId)
    {
      threads.emplace_back(
          [this, walkerId, stage, numberOfCycles, stageCycleOffset, startCycle, &synchronizationPoint]()
          {
            System& system = systems[walkerId];

            // each thread computes the total energies of its own walker
            if (stage == SimulationStage::PreInitialization || stage == SimulationStage::Initialization)
            {
              system.precomputeTotalRigidEnergy();
            }
            system.runningEnergies = system.computeTotalEnergies();

            BlockErrorEstimation estimation(numberOfBlocks, std::max(1uz, numberOfProductionCycles));
            // one cumulative snapshot has been pushed per block boundary crossed so far, so the
            // snapshot count restores the current block on a resume from a binary restart
            std::size_t currentBlock = blockCollectionMatrices[walkerId].size();

            // first quarter of production: Wang-Landau keeps exploring the window while the
            // freshly reset collection matrix fills with relaxed-configuration statistics;
            // afterwards the bias becomes the transition-matrix bias -ln Pi derived from those
            // rows, which flattens the walk over the whole window (on a binary-restart resume
            // past this point the switched flag is part of the checkpointed TMMC state)
            const std::size_t wangLandauProductionCycles = numberOfCycles / 4uz;

            for (std::size_t cycle = startCycle; cycle != numberOfCycles; ++cycle)
            {
              if (stage == SimulationStage::Production)
              {
                if (cycle == wangLandauProductionCycles)
                {
                  system.tmmc.switchToTMBias();
                }
                estimation.setCurrentSample(cycle);

                // cumulative production-only collection-matrix snapshot at every block boundary
                // (the per-block increments give the error bars of the analysis)
                if (estimation.currentBin != currentBlock)
                {
                  blockCollectionMatrices[walkerId].push_back(productionCollectionMatrix(walkerId));
                  currentBlock = estimation.currentBin;
                }
              }

              performWalkerCycle(walkerId, stage, estimation.currentBin);

              // time-evolution properties (number of molecules, volume): sampled over all stages,
              // indexed by the absolute cycle number; the writers gate on their own 'writeEvery'
              const std::size_t absoluteCycle = stageCycleOffset + cycle;
              system.samplePropertiesEvolution(absoluteCycle);
              if (system.propertyNumberOfMoleculesEvolution.has_value())
              {
                system.propertyNumberOfMoleculesEvolution->writeOutput(walkerId, absoluteCycle);
              }
              if (system.propertyVolumeEvolution.has_value())
              {
                system.propertyVolumeEvolution->writeOutput(walkerId, absoluteCycle);
              }

              if (stage == SimulationStage::Production)
              {
                system.sampleProperties(walkerId, estimation.currentBin, cycle);

                // analysis-property files (RDFs, density grid, histograms, molecule properties);
                // the writers gate on their own 'writeEvery'
                writeWalkerAnalysisOutputs(system, walkerId, cycle);

                // energy/pressure averages for the per-walker final report
                if (cycle % 10uz == 0uz || cycle % printEvery == 0uz)
                {
                  std::chrono::steady_clock::time_point time1 = std::chrono::steady_clock::now();
                  std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
                  system.currentEnergyStatus = molecularPressure.first;
                  system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
                  std::chrono::steady_clock::time_point time2 = std::chrono::steady_clock::now();

                  system.mc_moves_cputime.energyPressureComputation += (time2 - time1);
                  system.averageEnergies.addSample(estimation.currentBin, molecularPressure.first, system.weight());
                }
              }

              if (cycle % optimizeMCMovesEvery == 0uz)
              {
                system.optimizeMCMoves();
              }

              // Wang-Landau biasing-factor adjustment (all state is owned by this walker)
              if (stage == SimulationStage::Equilibration && cycle % rescaleWangLandauEvery == 0uz &&
                  !system.tmmc.lambdaChain())
              {
                for (Component& component : system.components)
                {
                  component.lambdaGC.WangLandauIteration(
                      PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors,
                      system.containsTheFractionalMolecule);
                }
                system.pairSwapLambdaWangLandauIteration(
                    PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors);
                system.reactionLambdaWangLandauIteration(
                    PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors);
              }

              if (cycle % printEvery == 0uz)
              {
                // each thread writes exclusively to its own walker stream: no locking needed
                system.loadings = LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent,
                                              system.simulationBox);

                std::ostream walkerStream(walkerStreams[walkerId].rdbuf());
                switch (stage)
                {
                  case SimulationStage::PreInitialization:
                    std::print(walkerStream, "{}", system.writePreInitializationStatusReport(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Initialization:
                    std::print(walkerStream, "{}", system.writeInitializationStatusReport(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Equilibration:
                    std::print(walkerStream, "{}", system.writeEquilibrationStatusReportMC(cycle, numberOfCycles));
                    break;
                  case SimulationStage::Production:
                  {
                    std::string status_line = std::format("Current cycle: {} out of {}\n", cycle, numberOfCycles);
                    std::print(walkerStream, "{}", system.writeProductionStatusReportMC(status_line));
                    break;
                  }
                  default:
                    break;
                }
                std::flush(walkerStream);

                // one combined progress line, from the first walker
                if (walkerId == 0uz)
                {
                  std::scoped_lock lock(outputMutex);
                  std::print(stream, "Parallel-TMMC cycle {} of {} (walker 0)\n", cycle, numberOfCycles);
                  std::flush(stream);
                }
              }

              // the only synchronization point between the threads: the periodic binary-restart
              // checkpoint, written by the barrier completion while all threads are parked
              if (writeBinaryRestartEvery != 0uz &&
                  ((cycle + 1uz) % writeBinaryRestartEvery == 0uz || cycle + 1uz == numberOfCycles))
              {
                if (walkerId == 0uz)
                {
                  checkpointCycle.store(cycle + 1uz, std::memory_order_relaxed);
                }
                synchronizationPoint.arrive_and_wait();
              }
            }

            // final snapshot: the cumulative production-only collection matrix at the end
            if (stage == SimulationStage::Production)
            {
              blockCollectionMatrices[walkerId].push_back(productionCollectionMatrix(walkerId));
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

void ParallelTMMC::output()
{
  std::size_t numberOfSteps = std::accumulate(stepsPerWalker.begin(), stepsPerWalker.end(), 0uz);

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
  std::print(stream, "===============================================================================\n");
  std::print(stream, "                             Simulation finished!\n");
  std::print(stream, "===============================================================================\n");
  std::print(stream, "\n");

  // energy drift check of every walker (energies recomputed in parallel, one thread per walker);
  // the final state used by the per-walker reports is refreshed in the same pass
  std::vector<RunningEnergy> recomputed(systems.size());
  {
    std::vector<std::jthread> threads;
    threads.reserve(systems.size());
    for (std::size_t walkerId = 0; walkerId < systems.size(); ++walkerId)
    {
      threads.emplace_back(
          [this, walkerId, &recomputed]()
          {
            System& system = systems[walkerId];
            recomputed[walkerId] = system.computeTotalEnergies();

            std::pair<EnergyStatus, double3x3> molecularPressure = system.computeMolecularPressure();
            system.currentEnergyStatus = molecularPressure.first;
            system.currentExcessPressureTensor = molecularPressure.second / system.simulationBox.volume;
            system.loadings = LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent,
                                          system.simulationBox);

            // final per-walker transition-matrix statistics file
            system.tmmc.writeStatistics();
          });
    }
  }

  writeWalkerFinalReports(recomputed);

  std::print(stream, "Energy drift per walker\n");
  std::print(stream, "===============================================================================\n\n");
  for (std::size_t walkerId = 0; walkerId < systems.size(); ++walkerId)
  {
    const RunningEnergy drift = systems[walkerId].runningEnergies - recomputed[walkerId];
    std::print(stream, "    walker {:4d} (temperature {:10.4f} [K], window [{}, {}]): drift {: .6e} [K]\n", walkerId,
               systems[walkerId].temperature, systems[walkerId].tmmc.minMacrostate,
               systems[walkerId].tmmc.maxMacrostate, Units::EnergyToKelvin * drift.potentialEnergy());
  }
  std::print(stream, "\n\n");

  std::print(stream, "Production run counting of the MC moves summed over walkers and components\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "{}", countTotal.writeMCMoveStatistics(numberOfSteps));
  std::print(stream, "\n\n");

  // macrostate coverage per walker: every state of the window must be visited for the stitched
  // ln Pi(N) to be reliable; the min/max visit counts diagnose the flatness of the biased walk
  std::print(stream, "Macrostate coverage per walker (production)\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "    walker    temperature [K]    window            visited    min visits    max visits\n");
  std::print(stream, "    -------------------------------------------------------------------------------\n");
  for (std::size_t walkerId = 0; walkerId < systems.size(); ++walkerId)
  {
    const System& system = systems[walkerId];

    // production-only visits: the cumulative histogram minus its production-start snapshot
    std::vector<std::size_t> histogram = system.tmmc.histogram;
    if (walkerId < productionStartHistograms.size() &&
        productionStartHistograms[walkerId].size() == histogram.size())
    {
      for (std::size_t index = 0; index < histogram.size(); ++index)
      {
        histogram[index] -= productionStartHistograms[walkerId][index];
      }
    }
    const std::size_t visited =
        static_cast<std::size_t>(std::ranges::count_if(histogram, [](std::size_t count) { return count > 0uz; }));
    const std::size_t minVisits = histogram.empty() ? 0uz : *std::ranges::min_element(histogram);
    const std::size_t maxVisits = histogram.empty() ? 0uz : *std::ranges::max_element(histogram);
    std::print(stream, "    {:6d}    {:15.4f}    [{:5d}, {:5d}]    {:4d}/{:<4d}   {:10d}    {:10d}\n", walkerId,
               system.temperature, system.tmmc.minMacrostate, system.tmmc.maxMacrostate, visited, histogram.size(),
               minVisits, maxVisits);

    outputJson["output"]["tmmc"]["coverage"][walkerId] = {
        {"temperature", system.temperature},
        {"window", std::vector<std::size_t>{system.tmmc.minMacrostate, system.tmmc.maxMacrostate}},
        {"visited", visited},
        {"states", histogram.size()},
        {"minVisits", minVisits},
        {"maxVisits", maxVisits}};
  }
  std::print(stream, "\n\n");

  // the transition-matrix analysis: combines the collection matrices of the windows per
  // temperature and writes ln Pi(N), the reweighted isotherms and the vapor-liquid coexistence
  performTransitionMatrixAnalysis();

  std::print(stream, "Production run CPU timings of the MC moves summed over walkers and components\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "{}", total.writeMCMoveCPUTimeStatistics(totalProductionSimulationTime));
  std::print(stream, "Pre-initialization simulation time: {:14f} [s]\n", totalPreInitializationSimulationTime.count());
  std::print(stream, "Initalization simulation time:  {:14f} [s]\n", totalInitializationSimulationTime.count());
  std::print(stream, "Equilibration simulation time:  {:14f} [s]\n", totalEquilibrationSimulationTime.count());
  std::print(stream, "Production simulation time:     {:14f} [s]\n", totalProductionSimulationTime.count());
  std::print(stream, "Analysis time:                  {:14f} [s]\n", totalAnalysisTime.count());
  std::print(stream, "Total simulation time:          {:14f} [s]\n", (totalSimulationTime + totalAnalysisTime).count());
  std::print(stream, "\n\n");
  std::flush(stream);

  outputJson["output"]["numberOfSteps"] = numberOfSteps;
  outputJson["output"]["cpuTimings"]["preInitialization"] = totalPreInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["initialization"] = totalInitializationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["equilibration"] = totalEquilibrationSimulationTime.count();
  outputJson["output"]["cpuTimings"]["production"] = totalProductionSimulationTime.count();
  outputJson["output"]["cpuTimings"]["analysis"] = totalAnalysisTime.count();
  outputJson["output"]["cpuTimings"]["total"] = (totalSimulationTime + totalAnalysisTime).count();

  std::ofstream json(outputJsonFileName);
  json << outputJson.dump(4);
}

void ParallelTMMC::performTransitionMatrixAnalysis()
{
  reweightedIsotherms.clear();

  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

  std::print(stream, "Transition-matrix analysis\n");
  std::print(stream, "===============================================================================\n\n");

  const std::size_t numberOfMacrostates = maxMacrostate - minMacrostate + 1uz;
  const std::size_t nLambda = std::max(1uz, systems.front().tmmc.lambdaBinCount());
  const std::size_t nChain = numberOfMacrostates * nLambda;

  // the code base is compiled with -ffast-math, so infinities must not occur: unsampled states
  // carry this finite 'log of zero' sentinel instead and drop out of the sums
  constexpr double logZero = -1e300;
  constexpr double logZeroThreshold = -1e299;

  // ln Pi from a collection matrix over a 1D chain by the detailed-balance recursion
  // ln Pi(i+1) = ln Pi(i) + ln P(i -> i+1) - ln P(i+1 -> i). For 1D TMMC the chain is N;
  // for (N, λ) it is the flattened nearest-neighbour walk. States outside the visited
  // range stay at the 'log of zero' sentinel.
  auto computeLogPi = [&](const std::vector<double3>& collectionMatrix) -> std::vector<double>
  {
    const std::size_t nStates = collectionMatrix.size();
    std::vector<double> logPi(nStates, logZero);

    auto rowTotal = [&](std::size_t index) -> double
    { return collectionMatrix[index].x + collectionMatrix[index].y + collectionMatrix[index].z; };

    std::size_t first = nStates;
    std::size_t last = 0uz;
    for (std::size_t index = 0; index < nStates; ++index)
    {
      if (rowTotal(index) > 0.0)
      {
        first = std::min(first, index);
        last = std::max(last, index);
      }
    }
    if (first >= nStates) return logPi;

    logPi[first] = 0.0;
    for (std::size_t index = first; index < last; ++index)
    {
      const double forwardTotal = rowTotal(index);
      const double reverseTotal = rowTotal(index + 1uz);
      const bool supportedLink = logPi[index] > logZeroThreshold && collectionMatrix[index].z > 0.0 &&
                                 collectionMatrix[index + 1uz].x > 0.0 && forwardTotal > 0.0 &&
                                 reverseTotal > 0.0;
      if (!supportedLink) continue;

      logPi[index + 1uz] = logPi[index] + std::log(collectionMatrix[index].z) - std::log(forwardTotal) -
                           std::log(collectionMatrix[index + 1uz].x) + std::log(reverseTotal);
    }
    return logPi;
  };

  // WHAM samples (N, U) only at λ = 0 (decoupled). Π_GC(N) ∝ Π_EE(N, k = 0).
  auto toMoleculeLogPi = [&](const std::vector<double>& chainLogPi) -> std::vector<double>
  {
    if (nLambda <= 1uz) return chainLogPi;
    std::vector<double> moleculeLogPi(numberOfMacrostates, logZero);
    for (std::size_t n = 0; n < numberOfMacrostates; ++n)
    {
      const std::size_t index = n * nLambda;
      if (index < chainLogPi.size()) moleculeLogPi[n] = chainLogPi[index];
    }
    return moleculeLogPi;
  };

  auto logSumExpRange = [&](const std::vector<double>& values, std::size_t begin, std::size_t end) -> double
  {
    double largest = logZero;
    for (std::size_t n = begin; n < end; ++n)
    {
      if (values[n] > logZeroThreshold) largest = std::max(largest, values[n]);
    }
    if (largest <= logZeroThreshold) return logZero;
    double sum = 0.0;
    for (std::size_t n = begin; n < end; ++n)
    {
      if (values[n] > logZeroThreshold) sum += std::exp(values[n] - largest);
    }
    return largest + std::log(sum);
  };

  // ln Pi(N; f) = ln Pi(N; f_ref) + N ln(f / f_ref) (exact); N is the global molecule count
  auto reweightedDistribution = [&](const std::vector<double>& logPi, double deltaLogFugacity) -> std::vector<double>
  {
    std::vector<double> logProbability(numberOfMacrostates, logZero);
    for (std::size_t index = 0; index < numberOfMacrostates; ++index)
    {
      if (logPi[index] <= logZeroThreshold) continue;
      logProbability[index] = logPi[index] + static_cast<double>(minMacrostate + index) * deltaLogFugacity;
    }
    return logProbability;
  };

  // conditional average of N over the macrostate subrange [begin, end)
  auto conditionalAverageMolecules = [&](const std::vector<double>& logProbability, std::size_t begin,
                                         std::size_t end) -> double
  {
    const double logPartition = logSumExpRange(logProbability, begin, end);
    if (logPartition <= logZeroThreshold) return 0.0;
    double average = 0.0;
    for (std::size_t index = begin; index < end; ++index)
    {
      if (logProbability[index] <= logZeroThreshold) continue;
      average += static_cast<double>(minMacrostate + index) * std::exp(logProbability[index] - logPartition);
    }
    return average;
  };

  auto averageMolecules = [&](const std::vector<double>& logProbability) -> double
  { return conditionalAverageMolecules(logProbability, 0uz, numberOfMacrostates); };

  // Places the phase cut at the deepest valley of ln P(N): the cut maximizing
  // min(peak left, peak right) - valley. Returns (cut, depth) or nothing when the distribution
  // is not bimodal (depth below threshold). Unsampled N are treated as deep-valley entries.
  constexpr double bimodalDepthThreshold = 0.2;
  auto findPhaseSplit = [&](const std::vector<double>& logProbability) -> std::optional<std::pair<std::size_t, double>>
  {
    const std::size_t size = logProbability.size();
    double lowestSampled = 1e300;
    for (const double value : logProbability)
    {
      if (value > logZeroThreshold) lowestSampled = std::min(lowestSampled, value);
    }
    std::vector<double> filled(size);
    for (std::size_t n = 0; n < size; ++n)
    {
      filled[n] = logProbability[n] > logZeroThreshold ? logProbability[n] : lowestSampled - 20.0;
    }

    std::vector<double> prefixMaximum(size), suffixMaximum(size);
    prefixMaximum[0] = filled[0];
    for (std::size_t n = 1; n < size; ++n) prefixMaximum[n] = std::max(prefixMaximum[n - 1], filled[n]);
    suffixMaximum[size - 1] = filled[size - 1];
    for (std::size_t n = size - 1; n > 0; --n) suffixMaximum[n - 1] = std::max(suffixMaximum[n], filled[n - 1]);

    std::size_t bestCut = 0uz;
    double bestDepth = -1e300;
    for (std::size_t cut = 1; cut + 1 < size; ++cut)
    {
      const double depth = std::min(prefixMaximum[cut - 1], suffixMaximum[cut + 1]) - filled[cut];
      if (depth > bestDepth)
      {
        bestDepth = depth;
        bestCut = cut;
      }
    }
    if (bestCut == 0uz || bestDepth < bimodalDepthThreshold)
    {
      return std::nullopt;
    }
    return std::make_pair(bestCut, bestDepth);
  };

  // unit conversions shared by all walkers (identical framework/box)
  const System& front = systems.front();
  const Component& frontComponent = front.components.front();
  const double volume = front.simulationBox.volume;
  double toMoleculesPerUnitCell = 1.0;
  double toMolePerKg = 0.0;
  double toMgPerG = 0.0;
  if (front.framework.has_value())
  {
    const int3 numberOfUnitCells = front.framework->numberOfUnitCells;
    toMoleculesPerUnitCell = 1.0 / static_cast<double>(numberOfUnitCells.x * numberOfUnitCells.y * numberOfUnitCells.z);
    const double frameworkMass = front.frameworkMass().value();
    toMolePerKg = 1000.0 / frameworkMass;
    toMgPerG = 1000.0 * frontComponent.totalMass / frameworkMass;
  }

  const std::vector<EquationOfState::FluidInput> fluidInputs = {
      {frontComponent.criticalTemperature, frontComponent.criticalPressure, frontComponent.acentricFactor, 1.0, true}};
  const double logPressureMinimum = std::log(reweightingPressureRange.first);
  const double logPressureMaximum = std::log(reweightingPressureRange.second);

  struct CoexistencePoint
  {
    double logBetaFugacity;
    double fugacityPa;
    double saturationPressurePa;
    bool saturationPressureFromEmptyBox;  // exact normalization (false: ideal-gas reference)
    double vaporDensity;                  // [kg/m^3]
    double liquidDensity;                 // [kg/m^3]
    double vaporMolecules;
    double liquidMolecules;
    double valleyDepth;
    std::size_t cut;
    std::vector<double> logProbability;
  };

  std::ofstream coexistence;
  const bool doVLE = !front.framework.has_value();
  if (doVLE)
  {
    coexistence.open("output/vle_coexistence.parallel_tmmc.txt", std::ios::trunc);
    std::print(coexistence, "# Parallel TMMC: vapor-liquid coexistence of {} (equal-weight criterion)\n",
               frontComponent.name);
    std::print(coexistence, "# box volume {:.4f} [A^3]; errors from the per-block collection-matrix increments\n",
               volume);
    std::print(coexistence,
               "# the saturation pressure follows from beta p V = ln Xi, normalized by the empty-box\n"
               "# state when sampled (column 13 = 0), otherwise approximately by an ideal-gas reference\n"
               "# at the dilute end of the scanned pressure range (column 13 = 1)\n");
    std::print(coexistence, "# column 1: temperature [K]\n");
    std::print(coexistence, "# column 2, 3: coexistence fugacity, error [Pa]\n");
    std::print(coexistence, "# column 4, 5: saturation pressure, error [Pa]\n");
    std::print(coexistence, "# column 6, 7: vapor density, error [kg/m^3]\n");
    std::print(coexistence, "# column 8, 9: liquid density, error [kg/m^3]\n");
    std::print(coexistence, "# column 10, 11: vapor, liquid peak [molecules]\n");
    std::print(coexistence, "# column 12: ln P(N) valley depth [-] (vanishes towards the critical point)\n");
    std::print(coexistence, "# column 13: saturation-pressure normalization (0 exact, 1 ideal-gas reference)\n\n");
  }

  bool anyApproximateNormalization = false;

  if (computeBET)
  {
    outputJson["output"]["tmmc"]["bet"] = nlohmann::json::array();
  }

  for (std::size_t temperatureIndex = 0; temperatureIndex < numberOfTemperatures; ++temperatureIndex)
  {
    const double temperature = temperatures[temperatureIndex];
    const double beta = 1.0 / (Units::KB * temperature);

    // the reference fugacity of this temperature (all windows share it)
    const System& windowFront = systems[walkerIndex(temperatureIndex, 0uz)];
    const Component& component = windowFront.components.front();
    const double referenceFugacity =
        component.molFraction * component.fugacityCoefficient.value_or(1.0) * windowFront.pressure;
    const double referenceLogFugacity = std::log(referenceFugacity);

    auto logBetaFugacityToDelta = [&](double fugacityInternal) -> double
    { return std::log(fugacityInternal) - referenceLogFugacity; };

    auto fugacityOfPressure = [&](double pressurePa) -> double
    {
      const std::vector<EquationOfState::FluidResult> fluidResults = EquationOfState::computeFluidProperties(
          temperature, pressurePa, fluidInputs, EquationOfState::Type::PengRobinson,
          EquationOfState::MixingRules::VanDerWaals);
      return fluidResults.front().fugacityCoefficient.value_or(1.0) * pressurePa / Units::PressureConversionFactor;
    };

    // Scatter one walker's collection matrix onto the global chain (zeros elsewhere).
    auto scatterWindow = [&](const System& system, const std::vector<double3>& local) -> std::vector<double3>
    {
      std::vector<double3> scattered(nChain, double3(0.0, 0.0, 0.0));
      if (system.tmmc.minMacrostate < minMacrostate) return scattered;
      const std::size_t offset = (system.tmmc.minMacrostate - minMacrostate) * nLambda;
      if (offset >= nChain) return scattered;
      const std::size_t count = std::min(local.size(), nChain - offset);
      for (std::size_t index = 0; index < count; ++index)
      {
        scattered[offset + index] = local[index];
      }
      return scattered;
    };

    // Pool the collection matrices of all windows on the global chain and reconstruct one
    // ln Π by the detailed-balance recursion. Every row entry is an accumulated sum of
    // unbiased acceptance probabilities, so summing rows across windows yields the
    // attempt-count weighted combination of their estimates of the same physical transition
    // probabilities. No gauge matching is needed and no seams appear at the window walls;
    // per-window reconstruction with overlap matching left ln Π discontinuities where the
    // poorly-sampled edge of one window met the well-sampled interior of the next.
    auto stitchLogPi = [&](auto localMatrixOfWalker) -> std::vector<double>
    {
      std::vector<double3> pooled(nChain, double3(0.0, 0.0, 0.0));
      for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
      {
        const std::size_t walkerId = walkerIndex(temperatureIndex, windowIndex);
        const System& system = systems[walkerId];
        const std::vector<double3> scattered = scatterWindow(system, localMatrixOfWalker(walkerId));
        for (std::size_t index = 0; index < nChain; ++index)
        {
          pooled[index] += scattered[index];
        }
      }
      return toMoleculeLogPi(computeLogPi(pooled));
    };

    // The Widom anchor. Pool the production test insertions of all windows per macrostate (they
    // are accumulated on the global N grid already) and turn them into the exact increment
    //
    //   ln Pi(N+1) - ln Pi(N) = ln(beta f V <W>_N / (N+1)),
    //
    // the same ratio the insertion acceptance rule uses, evaluated at the reference fugacity that
    // gauges ln Pi. A macrostate has an increment once it has enough insertions. The increment
    // extends the anchored block while the relative error of the mean is below the hopeless cut
    // (a flicker around 0.35 must not close the block; see the constant). Every counted increment
    // inside that block replaces the collection-matrix value.
    std::vector<double> widomLogIncrement(numberOfMacrostates, 0.0);
    std::vector<double> widomRelativeError(numberOfMacrostates, 0.0);
    std::vector<std::size_t> widomCount(numberOfMacrostates, 0uz);
    std::vector<bool> widomHasIncrement(numberOfMacrostates, false);
    std::vector<bool> widomPrecise(numberOfMacrostates, false);
    for (std::size_t index = 0; index + 1uz < numberOfMacrostates; ++index)
    {
      double weightSum = 0.0;
      double weightSquaredSum = 0.0;
      std::size_t count = 0uz;
      for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
      {
        const std::size_t walkerId = walkerIndex(temperatureIndex, windowIndex);
        weightSum += widomWeightSums[walkerId][index];
        weightSquaredSum += widomWeightSquaredSums[walkerId][index];
        count += widomInsertions[walkerId][index];
      }
      widomCount[index] = count;
      if (count == 0uz) continue;

      const double mean = weightSum / static_cast<double>(count);
      if (!(mean > 0.0)) continue;
      const double variance = std::max(0.0, weightSquaredSum / static_cast<double>(count) - mean * mean);
      const double relativeError = std::sqrt(variance / static_cast<double>(count)) / mean;
      widomRelativeError[index] = relativeError;
      if (count < minimumWidomInsertions) continue;

      const double moleculesAfter = static_cast<double>(minMacrostate + index + 1uz);
      widomLogIncrement[index] = std::log(beta * referenceFugacity * volume * mean / moleculesAfter);
      widomHasIncrement[index] = true;
      if (relativeError <= maximumWidomRelativeError) widomPrecise[index] = true;
    }

    std::size_t widomAnchorEnd = 0uz;  // one past the last anchored increment
    std::size_t consecutiveFailures = 0uz;
    for (std::size_t index = 0; index + 1uz < numberOfMacrostates; ++index)
    {
      if (widomPrecise[index])
      {
        consecutiveFailures = 0uz;
        widomAnchorEnd = index + 1uz;
      }
      else if (++consecutiveFailures >= widomAnchorFailureRun)
      {
        break;
      }
    }

    if (widomAnchorEnd > 0uz)
    {
      // the N = minMacrostate increment of an empty box is the Henry coefficient, reported in the
      // same units as the Widom-based value of the reweighted-histogram driver
      std::print(stream, "Widom anchor at {} K: increments out of N = {}--{} taken from test insertions\n",
                 temperature, minMacrostate, minMacrostate + widomAnchorEnd - 1uz);
      if (minMacrostate == 0uz && toMoleculesPerUnitCell > 0.0 && widomHasIncrement[0])
      {
        const double henryPerCell = std::exp(widomLogIncrement[0]) / referenceFugacity /
                                    Units::PressureConversionFactor * toMoleculesPerUnitCell;
        std::print(stream, "    Henry coefficient {:.6e} [molecules/cell/Pa] ({} test insertions)\n", henryPerCell,
                   widomCount[0]);
      }
      std::print(stream, "\n");
    }
    else
    {
      std::print(stream, "Widom anchor at {} K: no macrostate qualified, ln Pi is the collection-matrix estimate\n\n",
                 temperature);
    }

    // Replace the anchored increments and carry the accumulated correction into every macrostate
    // above them, which leaves all collection-matrix increments outside the block untouched.
    auto applyWidomAnchor = [&](std::vector<double>& logPi)
    {
      double shift = 0.0;
      // A valid insertion measurement at the first macrostate supplies its own gauge even when
      // the collection-matrix chain is disconnected there.
      if (!logPi.empty() && logPi[0] <= logZeroThreshold && widomAnchorEnd > 0uz && widomHasIncrement[0])
      {
        logPi[0] = 0.0;
      }
      for (std::size_t index = 0; index + 1uz < numberOfMacrostates; ++index)
      {
        if (logPi[index] <= logZeroThreshold) continue;

        if (index < widomAnchorEnd && widomHasIncrement[index])
        {
          const double rawNext = logPi[index + 1uz];
          const double anchoredNext = logPi[index] + widomLogIncrement[index];
          logPi[index + 1uz] = anchoredNext;
          if (rawNext > logZeroThreshold) shift = anchoredNext - rawNext;
        }
        else if (logPi[index + 1uz] > logZeroThreshold)
        {
          logPi[index + 1uz] += shift;
        }
      }
    };

    std::size_t availableBlocks = numberOfBlocks;
    for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
    {
      availableBlocks =
          std::min(availableBlocks, blockCollectionMatrices[walkerIndex(temperatureIndex, windowIndex)].size());
    }

    // Production-only statistics (like the per-block increments below). The cumulative matrix
    // still contains the Wang-Landau-era samples, taken while the walker descended through
    // low N with molecules that had not yet relaxed into the deep sites: their deletion
    // acceptances are overestimated by orders of magnitude, and because a row estimate is a
    // mean of acceptance probabilities with a heavy upper tail, those early entries dominate
    // the average no matter how long production runs. Discarding pre-production statistics
    // from the readout (the bias keeps using the cumulative matrix) removes that transient.
    std::vector<double> logPi =
        stitchLogPi([&](std::size_t walkerId) { return productionCollectionMatrix(walkerId); });
    applyWidomAnchor(logPi);

    // The blocks are anchored with the same (whole-production) test-insertion increments, so the
    // reported error over the anchored block carries the spread of the collection-matrix
    // increments above it but not the anchor's own uncertainty, which the file lists per
    // macrostate instead.
    std::vector<std::vector<double>> blockLogPi;
    blockLogPi.reserve(availableBlocks);
    for (std::size_t block = 0; block < availableBlocks; ++block)
    {
      blockLogPi.push_back(stitchLogPi([&](std::size_t walkerId)
      {
        const std::vector<std::vector<double3>>& snapshots = blockCollectionMatrices[walkerId];
        std::vector<double3> increment = snapshots[block];
        if (block > 0uz)
        {
          for (std::size_t index = 0; index < increment.size(); ++index)
          {
            increment[index] -= snapshots[block - 1uz][index];
          }
        }
        return increment;
      }));
      applyWidomAnchor(blockLogPi.back());
    }

    // total visit histogram over the windows (a pure diagnostic); (N, λ) visits are summed over λ
    std::vector<std::size_t> totalHistogram(numberOfMacrostates, 0uz);
    for (std::size_t windowIndex = 0; windowIndex < numberOfWindows; ++windowIndex)
    {
      const System& system = systems[walkerIndex(temperatureIndex, windowIndex)];
      const std::size_t offset = system.tmmc.minMacrostate - minMacrostate;
      const std::size_t windowLambda = std::max(1uz, system.tmmc.lambdaBinCount());
      for (std::size_t index = 0; index < system.tmmc.histogram.size(); ++index)
      {
        const std::size_t molIndex = offset + index / windowLambda;
        std::size_t visits = system.tmmc.histogram[index];
        const std::size_t walkerId = walkerIndex(temperatureIndex, windowIndex);
        if (walkerId < productionStartHistograms.size() &&
            index < productionStartHistograms[walkerId].size())
        {
          visits -= productionStartHistograms[walkerId][index];
        }
        if (molIndex < totalHistogram.size()) totalHistogram[molIndex] += visits;
      }
    }

    // ln Pi(N) at the reference fugacity, normalized to unit total probability, with the error
    // from the per-block solutions (normalized the same way, which removes the arbitrary gauge)
    {
      std::ofstream lnpiFile(std::format("tmmc/lnpi_{}.parallel_tmmc.txt", temperature), std::ios::trunc);
      std::print(lnpiFile, "# Parallel TMMC: macrostate probability distribution of {} at {} [K]\n",
                 frontComponent.name, temperature);
      std::print(lnpiFile, "# at the reference fugacity {:.6e} [Pa] (reference pressure {:.6e} [Pa])\n",
                 referenceFugacity * Units::PressureConversionFactor, referencePressure);
      std::print(lnpiFile, "# window collection matrices pooled on the global macrostate grid (production only)\n");
      if (nLambda > 1uz)
      {
        std::print(lnpiFile,
                   "# (N, λ) TMMC: ln Pi is the λ = 0 slice of the chain (WHAM samples the decoupled molecule)\n");
      }
      std::print(lnpiFile, "# reweight exactly with ln Pi(N; f) = ln Pi(N) + N ln(f / f_ref)\n");
      if (widomAnchorEnd > 0uz)
      {
        std::print(lnpiFile,
                   "# the increments out of N = {}--{} that have a test-insertion mean (column 7 = 1) are\n"
                   "# Widom-anchored: ln(beta f V <W>_N / (N+1)) from the test insertions of column 5 rather\n"
                   "# than from the collection matrix (the N = {} increment is the Henry coefficient); a\n"
                   "# relative error above {:.2f} ends the block after {} consecutive misses but does not\n"
                   "# leave a collection-matrix hole inside it; the increments above the block are the\n"
                   "# collection-matrix estimates\n",
                   minMacrostate, minMacrostate + widomAnchorEnd - 1uz, minMacrostate, maximumWidomRelativeError,
                   widomAnchorFailureRun);
      }
      else
      {
        std::print(lnpiFile, "# no macrostate reached {} test insertions with a usable mean: ln Pi is the\n"
                             "# collection-matrix estimate throughout\n", minimumWidomInsertions);
      }
      std::print(lnpiFile, "# column 1: N [molecules]\n");
      std::print(lnpiFile, "# column 2, 3: ln Pi(N), error (normalized to unit total probability)\n");
      std::print(lnpiFile, "# column 4: visits (summed over the windows)\n");
      std::print(lnpiFile, "# column 5, 6: test insertions at N, relative error of their mean Rosenbluth weight\n");
      std::print(lnpiFile, "# column 7: 1 when the increment out of N is Widom-anchored, 0 when it is not\n\n");

      const double logNormalization = logSumExpRange(logPi, 0uz, numberOfMacrostates);
      std::vector<double> blockLogNormalization(availableBlocks);
      for (std::size_t block = 0; block < availableBlocks; ++block)
      {
        blockLogNormalization[block] = logSumExpRange(blockLogPi[block], 0uz, numberOfMacrostates);
      }

      for (std::size_t index = 0; index < numberOfMacrostates; ++index)
      {
        if (logPi[index] <= logZeroThreshold) continue;
        const double normalized = logPi[index] - logNormalization;

        std::vector<double> blockValues;
        blockValues.reserve(availableBlocks);
        for (std::size_t block = 0; block < availableBlocks; ++block)
        {
          if (blockLogPi[block][index] > logZeroThreshold && blockLogNormalization[block] > logZeroThreshold)
          {
            blockValues.push_back(blockLogPi[block][index] - blockLogNormalization[block]);
          }
        }
        const double error = blockErrorEstimate(blockValues, normalized);

        std::print(lnpiFile, "{:6d}   {: .10e} {: .6e}   {}   {} {: .4e}   {}\n", minMacrostate + index, normalized,
                   error, totalHistogram[index], widomCount[index], widomRelativeError[index],
                   (index < widomAnchorEnd && widomHasIncrement[index]) ? 1 : 0);
      }
    }

    // the equilibrium, adsorption-branch and desorption-branch loadings at a given fugacity:
    // when the reweighted Pi(N) is bimodal the adsorption branch is the conditional average over
    // the low-density basin (N below the deepest valley; the metastable states followed on the
    // way up), the desorption branch the conditional average over the high-density basin, and
    // the equilibrium loading averages over both. Where Pi(N) is unimodal all three coincide.
    // branch order: 0 = equilibrium, 1 = adsorption, 2 = desorption
    auto branchLoadings = [&](const std::vector<double>& logPiInput,
                              double delta) -> std::pair<std::array<double, 3>, bool>
    {
      const std::vector<double> logProbability = reweightedDistribution(logPiInput, delta);
      const double equilibrium = averageMolecules(logProbability);
      const std::optional<std::pair<std::size_t, double>> split = findPhaseSplit(logProbability);
      if (!split.has_value())
      {
        return {{equilibrium, equilibrium, equilibrium}, false};
      }
      return {{equilibrium, conditionalAverageMolecules(logProbability, 0uz, split->first),
               conditionalAverageMolecules(logProbability, split->first, numberOfMacrostates)},
              true};
    };

    // reweighted isotherms <N>(f) on a log-spaced pressure grid (exact in the fugacity; the
    // fugacity coefficients from the Peng-Robinson equation of state at every point): the
    // equilibrium isotherm plus the adsorption and desorption branches (the hysteresis loop)
    {
      const std::array<std::string, 3> branchFileNames = {
          std::format("tmmc/equilibrium_isotherm_{}.parallel_tmmc.txt", temperature),
          std::format("tmmc/adsorption_isotherm_{}.parallel_tmmc.txt", temperature),
          std::format("tmmc/desorption_isotherm_{}.parallel_tmmc.txt", temperature)};
      const std::array<std::string_view, 3> branchDescriptions = {
          "equilibrium reweighted adsorption isotherm (averaged over both basins of Pi(N))",
          "adsorption branch (conditional average over the low-density basin of Pi(N);\n"
          "# the metastable states followed on adsorption - equals the equilibrium isotherm\n"
          "# where Pi(N) is unimodal)",
          "desorption branch (conditional average over the high-density basin of Pi(N);\n"
          "# the metastable states followed on desorption - equals the equilibrium isotherm\n"
          "# where Pi(N) is unimodal)"};

      std::array<std::ofstream, 3> branchFiles;
      for (std::size_t branch = 0; branch < 3uz; ++branch)
      {
        branchFiles[branch].open(branchFileNames[branch], std::ios::trunc);
        std::print(branchFiles[branch], "# Parallel TMMC: {} at {} [K] ({})\n", branchDescriptions[branch], temperature,
                   frontComponent.name);
        std::print(branchFiles[branch], "# errors from the per-block collection-matrix increments\n");
        std::print(branchFiles[branch],
                   "# the isotherm saturates artificially near the upper macrostate bound N = {}; points with\n"
                   "# <N> approaching the bound are truncated by the finite macrostate range\n",
                   maxMacrostate);
        std::print(branchFiles[branch], "# column 1: fugacity [Pa]\n");
        std::print(branchFiles[branch], "# column 2, 3: absolute loading, error [molecules/cell]\n");
        std::print(branchFiles[branch], "# column 4, 5: absolute loading, error [molecules/unit-cell]\n");
        std::print(branchFiles[branch], "# column 6, 7: absolute loading, error [mol/kg-framework]\n");
        std::print(branchFiles[branch], "# column 8, 9: absolute loading, error [mg/g-framework]\n");
        std::print(branchFiles[branch], "# column 10: pressure [Pa]\n");
        std::print(branchFiles[branch], "# column 11: bimodal (1 when Pi(N) has two basins at this pressure)\n\n");
      }

      TMMCIsotherm stored;
      stored.temperature = temperature;
      stored.equilibrium.reserve(reweightingNumberOfPressures);
      stored.adsorption.reserve(reweightingNumberOfPressures);
      stored.desorption.reserve(reweightingNumberOfPressures);

      for (std::size_t pressureIndex = 0; pressureIndex < reweightingNumberOfPressures; ++pressureIndex)
      {
        const double targetPressure =
            reweightingNumberOfPressures == 1uz
                ? reweightingPressureRange.first
                : std::exp(logPressureMinimum + static_cast<double>(pressureIndex) *
                                                    (logPressureMaximum - logPressureMinimum) /
                                                    static_cast<double>(reweightingNumberOfPressures - 1uz));
        const double fugacityInternal = fugacityOfPressure(targetPressure);
        const double delta = logBetaFugacityToDelta(fugacityInternal);

        const auto [loadings, bimodal] = branchLoadings(logPi, delta);

        std::array<std::vector<double>, 3> blockValues;
        for (std::size_t branch = 0; branch < 3uz; ++branch)
        {
          blockValues[branch].reserve(availableBlocks);
        }
        for (std::size_t block = 0; block < availableBlocks; ++block)
        {
          const std::array<double, 3> blockLoadings = branchLoadings(blockLogPi[block], delta).first;
          for (std::size_t branch = 0; branch < 3uz; ++branch)
          {
            blockValues[branch].push_back(blockLoadings[branch]);
          }
        }

        for (std::size_t branch = 0; branch < 3uz; ++branch)
        {
          const double loading = loadings[branch];
          const double loadingError = blockErrorEstimate(blockValues[branch], loading);
          std::print(branchFiles[branch],
                     "{: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e}   "
                     "{}\n",
                     fugacityInternal * Units::PressureConversionFactor, loading, loadingError,
                     toMoleculesPerUnitCell * loading, toMoleculesPerUnitCell * loadingError, toMolePerKg * loading,
                     toMolePerKg * loadingError, toMgPerG * loading, toMgPerG * loadingError, targetPressure,
                     bimodal ? 1 : 0);

          const TMMCIsothermPoint point{
              .pressure = targetPressure,
              .fugacity = fugacityInternal * Units::PressureConversionFactor,
              .moleculesPerCell = loading,
              .moleculesPerCellError = loadingError,
              .bimodal = bimodal};
          if (branch == 0uz)
          {
            stored.equilibrium.push_back(point);
          }
          else if (branch == 1uz)
          {
            stored.adsorption.push_back(point);
          }
          else
          {
            stored.desorption.push_back(point);
          }
        }
      }
      reweightedIsotherms.push_back(std::move(stored));
    }

    if (computeBET)
    {
      if (!front.framework.has_value())
      {
        if (temperatureIndex == 0uz)
        {
          std::print(stream, "BET extraction skipped: BET needs a framework (unit-cell mass and volume)\n\n");
        }
      }
      else if (front.components.empty())
      {
        if (temperatureIndex == 0uz)
        {
          std::print(stream, "BET extraction skipped: BET needs an adsorbate component\n\n");
        }
      }
      else
      {
        const BETProbeProperties probe = requireBETProbeProperties(front.components.front());
        const int3 numberOfUnitCells = front.framework->numberOfUnitCells;
        const double cells =
            static_cast<double>(numberOfUnitCells.x * numberOfUnitCells.y * numberOfUnitCells.z);
        const double unitCellVolume = front.simulationBox.volume / cells;
        const double unitCellMass = front.framework->unitCellMass;
        const TMMCIsotherm& storedIsotherm = reweightedIsotherms.back();

        auto pointsFromBranch = [&](const std::vector<TMMCIsothermPoint>& branch) -> std::vector<SimulatedIsothermPoint>
        {
          std::vector<SimulatedIsothermPoint> points;
          points.reserve(branch.size());
          for (const TMMCIsothermPoint& point : branch)
          {
            points.push_back(
                SimulatedIsothermPoint{point.pressure, toMoleculesPerUnitCell * point.moleculesPerCell});
          }
          return points;
        };

        const std::array<std::pair<std::string_view, const std::vector<TMMCIsothermPoint>*>, 3> branches = {
            {{"equilibrium", &storedIsotherm.equilibrium},
             {"adsorption", &storedIsotherm.adsorption},
             {"desorption", &storedIsotherm.desorption}}};

        std::print(stream, "BET extraction at {} [K]\n", temperature);
        nlohmann::json betJson;
        betJson["temperature"] = temperature;
        for (std::size_t branchIndex = 0; branchIndex < branches.size(); ++branchIndex)
        {
          const auto& [branchName, branchPoints] = branches[branchIndex];
          const BETSurfaceArea bet =
              fitNitrogenBET(pointsFromBranch(*branchPoints), unitCellMass, unitCellVolume, probe);

          std::optional<double> areaError;
          std::optional<double> monolayerError;
          std::optional<double> cError;
          nlohmann::json blockAreas = nlohmann::json::array();
          nlohmann::json blockMonolayers = nlohmann::json::array();
          nlohmann::json blockCConstants = nlohmann::json::array();
          std::size_t betBlocksUsed = 0uz;
          std::size_t betBlocksFailed = 0uz;

          // Jackknife the BET line inside the full-data Rouquerol window using the same
          // per-block ln Pi(N) that supply the isotherm loading errors.
          if (!bet.plateauReading && bet.monolayerCapacity > 0.0 && bet.windowHigh > bet.windowLow &&
              availableBlocks > 0uz)
          {
            std::vector<double> blockAreaValues;
            std::vector<double> blockMonolayerValues;
            std::vector<double> blockCValues;
            blockAreaValues.reserve(availableBlocks);
            blockMonolayerValues.reserve(availableBlocks);
            blockCValues.reserve(availableBlocks);

            for (std::size_t block = 0; block < availableBlocks; ++block)
            {
              std::vector<SimulatedIsothermPoint> blockPoints;
              blockPoints.reserve(branchPoints->size());
              for (const TMMCIsothermPoint& point : *branchPoints)
              {
                const double fugacityInternal = point.fugacity / Units::PressureConversionFactor;
                const double delta = logBetaFugacityToDelta(fugacityInternal);
                const double loading = branchLoadings(blockLogPi[block], delta).first[branchIndex];
                blockPoints.push_back(
                    SimulatedIsothermPoint{point.pressure, toMoleculesPerUnitCell * loading});
              }

              const BETSurfaceArea blockBet = fitNitrogenBETFixedWindow(
                  blockPoints, unitCellMass, unitCellVolume, bet.windowLow, bet.windowHigh, probe);
              if (!(blockBet.monolayerCapacity > 0.0) || !(blockBet.gravimetricArea > 0.0))
              {
                ++betBlocksFailed;
                continue;
              }
              ++betBlocksUsed;
              blockAreaValues.push_back(blockBet.gravimetricArea);
              blockMonolayerValues.push_back(blockBet.monolayerCapacity);
              blockCValues.push_back(blockBet.cConstant);
              blockAreas.push_back(blockBet.gravimetricArea);
              blockMonolayers.push_back(blockBet.monolayerCapacity);
              blockCConstants.push_back(blockBet.cConstant);
            }

            if (betBlocksUsed >= 3uz)
            {
              areaError = blockErrorEstimate(blockAreaValues, bet.gravimetricArea);
              monolayerError = blockErrorEstimate(blockMonolayerValues, bet.monolayerCapacity);
              cError = blockErrorEstimate(blockCValues, bet.cConstant);
            }
          }

          if (branchName == "equilibrium")
          {
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
          }
          const std::string betFile =
              std::format("tmmc/bet_{}_{}.parallel_tmmc.txt", branchName, temperature);
          std::ofstream table(betFile, std::ios::trunc);
          std::print(table, "# Parallel TMMC: BET plot of the {} isotherm at {} [K] ({})\n", branchName,
                     temperature, frontComponent.name);
          writeNitrogenBETTable(table, bet, probe);
          nlohmann::json entry = nitrogenBETJson(bet, probe);
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
          betJson[std::string(branchName)] = std::move(entry);
        }
        outputJson["output"]["tmmc"]["bet"].push_back(std::move(betJson));
        std::print(stream, "    BET plots written to tmmc/bet_{{branch}}_{}.parallel_tmmc.txt\n\n", temperature);
      }
    }

    outputJson["output"]["tmmc"]["temperatures"].push_back(temperature);

    // vapor-liquid coexistence (bulk boxes only) by the equal-weight criterion (Wilding): the
    // fugacity is bisected until the vapor peak (N < N_cut) and the liquid peak (N >= N_cut) of
    // the reweighted Pi(N) carry equal probability weight; the cut is placed at the deepest
    // valley of ln Pi(N) between the two peaks
    if (!doVLE) continue;

    auto solveCoexistence = [&](const std::vector<double>& logPiInput) -> std::optional<CoexistencePoint>
    {
      auto weightDifference = [&](double delta) -> std::optional<std::pair<double, std::size_t>>
      {
        const std::vector<double> logProbability = reweightedDistribution(logPiInput, delta);
        const std::optional<std::pair<std::size_t, double>> split = findPhaseSplit(logProbability);
        if (!split.has_value()) return std::nullopt;
        const double logVaporWeight = logSumExpRange(logProbability, 0uz, split->first);
        const double logLiquidWeight = logSumExpRange(logProbability, split->first, logProbability.size());
        return std::make_pair(logLiquidWeight - logVaporWeight, split->first);
      };

      // bracket search over the (log-spaced) scanned pressure range
      constexpr std::size_t numberOfScanPoints = 200uz;
      double lowerDelta = 0.0;
      double upperDelta = 0.0;
      double lowerDifference = 0.0;
      bool havePrevious = false;
      bool haveBracket = false;
      for (std::size_t scanIndex = 0; scanIndex < numberOfScanPoints; ++scanIndex)
      {
        const double pressure =
            std::exp(logPressureMinimum + static_cast<double>(scanIndex) * (logPressureMaximum - logPressureMinimum) /
                                              static_cast<double>(numberOfScanPoints - 1uz));
        const double delta = logBetaFugacityToDelta(fugacityOfPressure(pressure));
        const std::optional<std::pair<double, std::size_t>> difference = weightDifference(delta);
        if (!difference.has_value())
        {
          havePrevious = false;
          continue;
        }
        if (havePrevious && lowerDifference < 0.0 && difference->first >= 0.0)
        {
          upperDelta = delta;
          haveBracket = true;
          break;
        }
        lowerDelta = delta;
        lowerDifference = difference->first;
        havePrevious = true;
      }
      if (!haveBracket) return std::nullopt;

      // bisection on ln f to the equal-weight point
      for (std::size_t iterationIndex = 0; iterationIndex < 100uz; ++iterationIndex)
      {
        const double midpoint = 0.5 * (lowerDelta + upperDelta);
        const std::optional<std::pair<double, std::size_t>> difference = weightDifference(midpoint);
        if (!difference.has_value()) break;  // lost bimodality mid-bracket; use the current interval
        if (difference->first < 0.0)
        {
          lowerDelta = midpoint;
        }
        else
        {
          upperDelta = midpoint;
        }
        if (std::abs(difference->first) < 1e-10) break;
      }

      const double delta = 0.5 * (lowerDelta + upperDelta);
      const std::vector<double> logProbability = reweightedDistribution(logPiInput, delta);
      const std::optional<std::pair<std::size_t, double>> split = findPhaseSplit(logProbability);
      if (!split.has_value()) return std::nullopt;
      const std::size_t cut = split->first;

      const double logVaporWeight = logSumExpRange(logProbability, 0uz, cut);
      const double logLiquidWeight = logSumExpRange(logProbability, cut, logProbability.size());
      double vaporMolecules = 0.0;
      double liquidMolecules = 0.0;
      for (std::size_t index = 0; index < logProbability.size(); ++index)
      {
        if (logProbability[index] <= logZeroThreshold) continue;
        const double n = static_cast<double>(minMacrostate + index);
        if (index < cut)
        {
          vaporMolecules += n * std::exp(logProbability[index] - logVaporWeight);
        }
        else
        {
          liquidMolecules += n * std::exp(logProbability[index] - logLiquidWeight);
        }
      }

      const double coexistenceFugacityInternal = std::exp(referenceLogFugacity + delta);

      // Saturation pressure from beta p V = ln Xi. The absolute normalization of Xi is fixed by
      // the empty-box state (N = 0, weight exactly one) when the macrostate range starts at zero
      // and the state was visited; otherwise it is fixed approximately at the dilute end of the
      // scanned pressure range, where the vapor is nearly ideal and beta p V ~ beta f V.
      double saturationPressurePa;
      const double logPartitionAtCoexistence = logSumExpRange(logProbability, 0uz, logProbability.size());
      const bool saturationPressureFromEmptyBox = (minMacrostate == 0uz) && (logProbability[0] > logZeroThreshold);
      if (saturationPressureFromEmptyBox)
      {
        saturationPressurePa =
            ((logPartitionAtCoexistence - logProbability[0]) / (beta * volume)) * Units::PressureConversionFactor;
      }
      else
      {
        const double referenceFugacityDilute = fugacityOfPressure(std::exp(logPressureMinimum));
        const std::vector<double> referenceLogProbability =
            reweightedDistribution(logPiInput, logBetaFugacityToDelta(referenceFugacityDilute));
        const double logPartitionAtReference =
            logSumExpRange(referenceLogProbability, 0uz, referenceLogProbability.size());
        const double referenceBetaPressureVolume = beta * referenceFugacityDilute * volume;
        saturationPressurePa =
            ((logPartitionAtCoexistence - logPartitionAtReference + referenceBetaPressureVolume) / (beta * volume)) *
            Units::PressureConversionFactor;
      }

      const double toKgPerCubicMeter = frontComponent.totalMass * Units::DensityConversionFactor / volume;
      return CoexistencePoint{.logBetaFugacity = std::log(beta * coexistenceFugacityInternal),
                              .fugacityPa = coexistenceFugacityInternal * Units::PressureConversionFactor,
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

    if (temperatureIndex == 0uz)
    {
      std::print(stream, "    vapor-liquid coexistence (equal-weight criterion)\n");
      std::print(stream,
                 "    temperature [K]    fugacity [Pa]        P_sat [Pa]           "
                 "rho_vap [kg/m^3]     rho_liq [kg/m^3]\n");
      std::print(stream,
                 "    -----------------------------------------------------------------"
                 "----------------------------------\n");
    }

    const std::optional<CoexistencePoint> point = solveCoexistence(logPi);
    if (!point.has_value())
    {
      std::print(coexistence,
                 "# {} [K]: no coexistence found (supercritical, pressure range does not bracket the\n"
                 "#   transition, or the macrostate range does not span both phases)\n",
                 temperature);
      std::print(stream, "    {:15.4f}    no coexistence found in the scanned pressure range\n", temperature);
      continue;
    }

    // per-block coexistence solves for the error bars
    std::vector<double> blockFugacities, blockPressures, blockVaporDensities, blockLiquidDensities;
    for (std::size_t block = 0; block < availableBlocks; ++block)
    {
      const std::optional<CoexistencePoint> blockPoint = solveCoexistence(blockLogPi[block]);
      if (!blockPoint.has_value()) continue;
      blockFugacities.push_back(blockPoint->fugacityPa);
      blockVaporDensities.push_back(blockPoint->vaporDensity);
      blockLiquidDensities.push_back(blockPoint->liquidDensity);
      blockPressures.push_back(blockPoint->saturationPressurePa);
    }
    const double fugacityError = blockErrorEstimate(blockFugacities, point->fugacityPa);
    const double pressureError = blockErrorEstimate(blockPressures, point->saturationPressurePa);
    const double vaporDensityError = blockErrorEstimate(blockVaporDensities, point->vaporDensity);
    const double liquidDensityError = blockErrorEstimate(blockLiquidDensities, point->liquidDensity);

    anyApproximateNormalization = anyApproximateNormalization || !point->saturationPressureFromEmptyBox;

    std::print(coexistence,
               "{:10.4f}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   {: .6e} {: .6e}   "
               "{: .6e} {: .6e}   {: .6e}   {}\n",
               temperature, point->fugacityPa, fugacityError, point->saturationPressurePa, pressureError,
               point->vaporDensity, vaporDensityError, point->liquidDensity, liquidDensityError, point->vaporMolecules,
               point->liquidMolecules, point->valleyDepth, point->saturationPressureFromEmptyBox ? 0 : 1);

    std::print(stream,
               "    {:15.4f}    {: .6e} ± {:.2e}   {: .6e} ± {:.2e}{}  {:9.4f} ± {:7.4f}   "
               "{:9.4f} ± {:7.4f}\n",
               temperature, point->fugacityPa, fugacityError, point->saturationPressurePa, pressureError,
               point->saturationPressureFromEmptyBox ? ' ' : '*', point->vaporDensity, vaporDensityError,
               point->liquidDensity, liquidDensityError);

    // the coexistence molecule-number distribution (for inspection and finite-size scaling)
    std::ofstream distribution(std::format("output/vle_distribution_{}.parallel_tmmc.txt", temperature),
                               std::ios::trunc);
    std::print(distribution, "# Parallel TMMC: P(N) at coexistence, {} at {} [K]\n", frontComponent.name, temperature);
    std::print(distribution, "# coexistence fugacity {:.6e} [Pa], phase cut at N = {}\n", point->fugacityPa,
               minMacrostate + point->cut);
    std::print(distribution, "# column 1: N [molecules]\n");
    std::print(distribution, "# column 2: P(N) (normalized)\n");
    std::print(distribution, "# column 3: ln P(N)\n\n");
    const double logNormalization = logSumExpRange(point->logProbability, 0uz, point->logProbability.size());
    for (std::size_t index = 0; index < point->logProbability.size(); ++index)
    {
      if (point->logProbability[index] <= logZeroThreshold) continue;
      const double logNormalized = point->logProbability[index] - logNormalization;
      std::print(distribution, "{:6d}   {: .6e}   {: .6e}\n", minMacrostate + index, std::exp(logNormalized),
                 logNormalized);
    }

    nlohmann::json vleEntry;
    vleEntry["temperature"] = temperature;
    vleEntry["fugacity"] = point->fugacityPa;
    vleEntry["fugacityError"] = fugacityError;
    vleEntry["saturationPressure"] = point->saturationPressurePa;
    vleEntry["saturationPressureError"] = pressureError;
    vleEntry["vaporDensity"] = point->vaporDensity;
    vleEntry["vaporDensityError"] = vaporDensityError;
    vleEntry["liquidDensity"] = point->liquidDensity;
    vleEntry["liquidDensityError"] = liquidDensityError;
    vleEntry["valleyDepth"] = point->valleyDepth;
    vleEntry["saturationPressureNormalization"] =
        point->saturationPressureFromEmptyBox ? "empty-box" : "ideal-gas reference";
    outputJson["output"]["tmmc"]["vaporLiquidEquilibrium"].push_back(vleEntry);
  }

  if (doVLE)
  {
    std::print(stream, "\n");
    if (anyApproximateNormalization)
    {
      std::print(stream, "    (* saturation pressure normalized approximately by an ideal-gas reference at the\n");
      std::print(stream, "       dilute end of the pressure range; the empty-box state N = 0 was never sampled -\n");
      std::print(stream, "       set 'MacroStateMinimumNumberOfMolecules' to 0 for the exact normalization)\n");
    }
    std::print(stream, "    coexistence written to output/vle_coexistence.parallel_tmmc.txt\n");
    std::print(stream, "    coexistence P(N) written to output/vle_distribution_{{T}}.parallel_tmmc.txt\n\n");
  }

  std::print(stream, "    ln Pi(N) written to tmmc/lnpi_{{T}}.parallel_tmmc.txt\n");
  std::print(stream, "    equilibrium isotherms written to tmmc/equilibrium_isotherm_{{T}}.parallel_tmmc.txt\n");
  std::print(stream, "    adsorption/desorption branches (hysteresis loop) written to\n");
  std::print(stream,
             "    tmmc/adsorption_isotherm_{{T}}.parallel_tmmc.txt and "
             "tmmc/desorption_isotherm_{{T}}.parallel_tmmc.txt\n");
  if (computeBET && front.framework.has_value())
  {
    std::print(stream, "    BET plots written to tmmc/bet_{{branch}}_{{T}}.parallel_tmmc.txt\n");
  }
  std::print(stream, "\n\n");
  std::flush(stream);

  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  totalAnalysisTime += (t2 - t1);
}

void ParallelTMMC::writeWalkerFinalReports(std::vector<RunningEnergy>& recomputed)
{
  for (std::size_t walkerId = 0; walkerId < systems.size(); ++walkerId)
  {
    System& system = systems[walkerId];
    std::ostream walkerStream(walkerStreams[walkerId].rdbuf());

    std::print(walkerStream, "\n");
    std::print(walkerStream, "===============================================================================\n");
    std::print(walkerStream, "                             Simulation finished!\n");
    std::print(walkerStream, "===============================================================================\n");
    std::print(walkerStream, "\n");

    std::string status_line = std::format("Final state after {} cycles\n", numberOfProductionCycles);
    std::print(walkerStream, "{}", system.writeProductionStatusReportMC(status_line));

    const RunningEnergy drift = system.runningEnergies - recomputed[walkerId];
    walkerStream << system.runningEnergies.printMCDiff(recomputed[walkerId]);
    std::print(walkerStream, "\n\n");

    std::print(walkerStream, "Monte-Carlo moves statistics\n");
    std::print(walkerStream, "===============================================================================\n\n");
    std::print(walkerStream, "{}", system.writeMCMoveStatistics());

    std::print(walkerStream, "Production run CPU timings of the MC moves of this walker\n");
    std::print(walkerStream, "===============================================================================\n\n");
    for (std::size_t componentId{0}; const Component& component : system.components)
    {
      std::print(walkerStream, "{}",
                 component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(componentId, component.name));
      ++componentId;
    }
    std::print(walkerStream, "{}", system.mc_moves_cputime.writeMCMoveCPUTimeStatistics());
    std::print(walkerStream, "\n\n");

    // under the flattening transition-matrix bias the direct averages are flat-histogram
    // averages, not grand-canonical ones; they are reported as diagnostics only
    std::print(walkerStream, "NOTE: the walker samples under the flattening transition-matrix bias; the\n");
    std::print(walkerStream, "      averages below are biased flat-histogram averages (diagnostics only).\n");
    std::print(walkerStream, "      The physical results are in the combined analysis files.\n\n");
    std::print(
        walkerStream, "{}",
        system.averageEnergies.writeAveragesStatistics(system.hasExternalField, system.framework, system.components));
    if (!(system.framework.has_value() && system.framework->rigid))
    {
      std::print(walkerStream, "{}", system.averagePressure.writeAveragesStatistics());
    }
    std::print(walkerStream, "{}",
               system.averageLoadings.writeAveragesStatistics(
                   system.components, system.frameworkMass(),
                   system.framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));
    std::flush(walkerStream);

    // final flush of the analysis-property files (cycle 0 bypasses the 'writeEvery' gate)
    writeWalkerAnalysisOutputs(system, walkerId, 0uz);

    // final flush of the time-evolution files: the total cycle count is rounded up to a multiple
    // of the writer's own 'writeEvery' so its gate passes and all collected samples are written
    if (system.propertyNumberOfMoleculesEvolution.has_value() &&
        system.propertyNumberOfMoleculesEvolution->writeEvery.value_or(0uz) > 0uz)
    {
      const std::size_t writeEvery = system.propertyNumberOfMoleculesEvolution->writeEvery.value();
      system.propertyNumberOfMoleculesEvolution->writeOutput(
          walkerId, ((absoluteCycleOffset + writeEvery - 1uz) / writeEvery) * writeEvery);
    }
    if (system.propertyVolumeEvolution.has_value() && system.propertyVolumeEvolution->writeEvery.value_or(0uz) > 0uz)
    {
      const std::size_t writeEvery = system.propertyVolumeEvolution->writeEvery.value();
      system.propertyVolumeEvolution->writeOutput(walkerId,
                                                  ((absoluteCycleOffset + writeEvery - 1uz) / writeEvery) * writeEvery);
    }

    // per-walker json statistics
    walkerJsons[walkerId]["output"]["runningEnergies"] = system.runningEnergies.jsonMC();
    walkerJsons[walkerId]["output"]["recomputedEnergies"] = recomputed[walkerId].jsonMC();
    walkerJsons[walkerId]["output"]["drift"] = drift.jsonMC();
    walkerJsons[walkerId]["output"]["MCMoveStatistics"]["system"] = system.jsonMCMoveStatistics();
    walkerJsons[walkerId]["output"]["cpuTimings"]["system"] =
        system.mc_moves_cputime.jsonSystemMCMoveCPUTimeStatistics();
    for (const Component& component : system.components)
    {
      walkerJsons[walkerId]["output"]["cpuTimings"][component.name] =
          component.mc_moves_cputime.jsonComponentMCMoveCPUTimeStatistics();
    }
    walkerJsons[walkerId]["properties"]["averageEnergies"] =
        system.averageEnergies.jsonAveragesStatistics(system.hasExternalField, system.framework, system.components);
    walkerJsons[walkerId]["properties"]["averagePressure"] = system.averagePressure.jsonAveragesStatistics();

    std::ofstream json(walkerJsonFileNames[walkerId]);
    json << walkerJsons[walkerId].dump(4);
  }
}

void ParallelTMMC::writeBinaryRestartFile(std::size_t cyclesCompleted) noexcept
{
  cyclesCompletedThisStage = cyclesCompleted;

  ::writeBinaryRestartFile(*this);
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ParallelTMMC& ptmmc)
{
  archive << ptmmc.versionNumber;

  archive << ptmmc.random;

  archive << ptmmc.numberOfProductionCycles;
  archive << ptmmc.numberOfPreInitializationCycles;
  archive << ptmmc.numberOfInitializationCycles;
  archive << ptmmc.numberOfEquilibrationCycles;

  archive << ptmmc.printEvery;
  archive << ptmmc.optimizeMCMovesEvery;
  archive << ptmmc.rescaleWangLandauEvery;
  archive << ptmmc.writeBinaryRestartEvery;

  archive << ptmmc.numberOfBlocks;

  archive << ptmmc.reweightingPressureRange;
  archive << ptmmc.reweightingNumberOfPressures;

  archive << ptmmc.simulationStage;
  archive << ptmmc.cyclesCompletedThisStage;

  archive << ptmmc.temperatures;
  archive << ptmmc.referencePressure;
  archive << ptmmc.numberOfTemperatures;
  archive << ptmmc.numberOfWindows;
  archive << ptmmc.numberOfWalkers;

  archive << ptmmc.minMacrostate;
  archive << ptmmc.maxMacrostate;
  archive << ptmmc.windowBoundaries;

  archive << ptmmc.systems;
  archive << ptmmc.randoms;

  archive << ptmmc.blockCollectionMatrices;
  archive << ptmmc.productionStartCollectionMatrices;
  archive << ptmmc.productionStartHistograms;

  archive << ptmmc.widomWeightSums;
  archive << ptmmc.widomWeightSquaredSums;
  archive << ptmmc.widomInsertions;

  archive << ptmmc.stepsPerWalker;
  archive << ptmmc.absoluteCycleOffset;

  archive << ptmmc.totalPreInitializationSimulationTime;
  archive << ptmmc.totalInitializationSimulationTime;
  archive << ptmmc.totalEquilibrationSimulationTime;
  archive << ptmmc.totalProductionSimulationTime;
  archive << ptmmc.totalAnalysisTime;
  archive << ptmmc.totalSimulationTime;

  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ParallelTMMC& ptmmc)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > ptmmc.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'ParallelTMMC' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> ptmmc.random;

  archive >> ptmmc.numberOfProductionCycles;
  archive >> ptmmc.numberOfPreInitializationCycles;
  archive >> ptmmc.numberOfInitializationCycles;
  archive >> ptmmc.numberOfEquilibrationCycles;

  archive >> ptmmc.printEvery;
  archive >> ptmmc.optimizeMCMovesEvery;
  archive >> ptmmc.rescaleWangLandauEvery;
  archive >> ptmmc.writeBinaryRestartEvery;

  archive >> ptmmc.numberOfBlocks;

  archive >> ptmmc.reweightingPressureRange;
  archive >> ptmmc.reweightingNumberOfPressures;

  archive >> ptmmc.simulationStage;
  archive >> ptmmc.cyclesCompletedThisStage;

  archive >> ptmmc.temperatures;
  archive >> ptmmc.referencePressure;
  archive >> ptmmc.numberOfTemperatures;
  archive >> ptmmc.numberOfWindows;
  archive >> ptmmc.numberOfWalkers;

  archive >> ptmmc.minMacrostate;
  archive >> ptmmc.maxMacrostate;
  archive >> ptmmc.windowBoundaries;

  archive >> ptmmc.systems;
  archive >> ptmmc.randoms;

  archive >> ptmmc.blockCollectionMatrices;
  archive >> ptmmc.productionStartCollectionMatrices;
  archive >> ptmmc.productionStartHistograms;

  archive >> ptmmc.widomWeightSums;
  archive >> ptmmc.widomWeightSquaredSums;
  archive >> ptmmc.widomInsertions;

  archive >> ptmmc.stepsPerWalker;
  archive >> ptmmc.absoluteCycleOffset;

  archive >> ptmmc.totalPreInitializationSimulationTime;
  archive >> ptmmc.totalInitializationSimulationTime;
  archive >> ptmmc.totalEquilibrationSimulationTime;
  archive >> ptmmc.totalProductionSimulationTime;
  archive >> ptmmc.totalAnalysisTime;
  archive >> ptmmc.totalSimulationTime;

  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error("ParallelTMMC: error in binary restart\n");
  }

  return archive;
}
