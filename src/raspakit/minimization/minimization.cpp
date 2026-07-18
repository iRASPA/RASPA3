module;

module minimization;

import std;

import input_reader;
import atom;
import system;
import randomnumbers;
import mc_moves;
import property_lambda_probability_histogram;
import component;
import intra_molecular_potentials;
import int3;
import double3;
import double3x3;
import simulationbox;
import framework;
import forcefield;
import skbandpath;
import skspacegroup;
import sksymmetrycell;
import generalized_hessian;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import minimization_cell_layout;
import elastic_constants;
import normal_modes;
import phonon_kpath;
import phonon_dynamical_matrix;
import units;

namespace
{
void clipCellStep(BakerStep& step, const MinimizationDofLayout& layout, const CellMinimizationLayout& cellLayout,
                  double maximumNorm)
{
  if (cellLayout.empty()) return;
  double3x3 strain{};
  for (std::size_t a = 0; a < cellLayout.size(); ++a)
  {
    strain += step.displacement[*layout.cellDof(a)] * cellLayout.bases[a];
  }
  double normSquared = 0.0;
  for (std::size_t column = 0; column < 3; ++column)
  {
    for (std::size_t row = 0; row < 3; ++row)
    {
      normSquared += strain.mm[column][row] * strain.mm[column][row];
    }
  }
  const double norm = std::sqrt(normSquared);
  if (norm <= maximumNorm || norm == 0.0) return;
  const double scale = maximumNorm / norm;
  for (std::size_t a = 0; a < cellLayout.size(); ++a)
  {
    step.displacement[*layout.cellDof(a)] *= scale;
  }
  step.stepNorm =
      std::sqrt(std::inner_product(step.displacement.begin(), step.displacement.end(), step.displacement.begin(), 0.0));
}

bool useNegativeModeEscape(BakerStep& step, std::span<const double> gradient, const MinimizationOptions& options)
{
  if (step.negativeModes == 0 || step.lowestMode.size() != step.displacement.size())
  {
    return false;
  }
  const double collapseThreshold = std::min(1.0e-2, 0.1 * options.maximumStepLength);
  if (step.stepNorm >= collapseThreshold)
  {
    return false;
  }

  const double escapeLength = std::min(5.0e-2, options.maximumStepLength);
  const double projection = std::inner_product(gradient.begin(), gradient.end(), step.lowestMode.begin(), 0.0);
  const double sign = projection > 0.0 ? -1.0 : 1.0;
  for (std::size_t i = 0; i < step.displacement.size(); ++i)
  {
    step.displacement[i] = sign * escapeLength * step.lowestMode[i];
  }
  step.stepNorm = escapeLength;
  step.convergenceReason = "escaping lowest negative Hessian mode";
  return true;
}

// Build a symmetry-aware default phonon band path (HPKOT high-symmetry points and recommended path) for the
// framework of `system`. Symmetry is detected on the final (minimized) simulation cell together with all of
// its framework atoms: findSpaceGroup reduces the supercell to its true primitive cell internally, so the
// returned transformation maps the conventional cell directly onto the minimized simulation cell and the
// high-symmetry points can be expressed in the simulation-cell reciprocal basis without any extra supercell
// bookkeeping. Returns an empty list when no framework is present or the symmetry could not be detected, in
// which case the caller falls back to a simple reciprocal-axis star.
std::vector<PhononPathNode> symmetryDefaultPhononPath(const System& system)
{
  std::vector<PhononPathNode> nodes;
  if (!system.framework.has_value()) return nodes;

  const SimulationBox& box = system.simulationBox;

  std::vector<std::tuple<double3, std::size_t, double>> atoms;
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atoms.reserve(frameworkAtoms.size());
  for (const Atom& atom : frameworkAtoms)
  {
    double3 fractional = box.inverseCell * atom.position;
    atoms.emplace_back(fractional, static_cast<std::size_t>(atom.type), 1.0);
  }
  if (atoms.empty()) return nodes;

  std::optional<SKBandPath> bandPath =
      SKBrillouinZonePath(box.cell, atoms, /*allowPartialOccupancies=*/true, /*symmetryPrecision=*/1e-4);
  if (!bandPath) return nodes;

  auto toNode = [&](const std::string& label) -> std::optional<PhononPathNode>
  {
    std::optional<double3> coordinates = bandPath->coordinatesForLabel(label);
    if (!coordinates) return std::nullopt;
    // Map from the primitive reciprocal basis directly into the (minimized) simulation-cell reciprocal basis.
    double3 fractional = bandPath->reciprocalFractionalInAnalyzedCell(*coordinates);
    std::string display = (label == "GAMMA") ? std::string{"G"} : label;
    return PhononPathNode{.kPoint = fractional, .label = display};
  };

  std::string previousLabel;
  for (const auto& [from, to] : bandPath->segments)
  {
    std::optional<PhononPathNode> endNode = toNode(to);
    if (!endNode)
    {
      nodes.clear();
      return nodes;
    }
    if (nodes.empty())
    {
      std::optional<PhononPathNode> startNode = toNode(from);
      if (!startNode)
      {
        nodes.clear();
        return nodes;
      }
      nodes.push_back(*startNode);
    }
    else if (from != previousLabel)
    {
      // Discontinuity in the recommended path: anchor a new sub-path at `from`.
      std::optional<PhononPathNode> startNode = toNode(from);
      if (!startNode)
      {
        nodes.clear();
        return nodes;
      }
      startNode->startsNewSegment = true;
      nodes.push_back(*startNode);
    }
    nodes.push_back(*endNode);
    previousLabel = to;
  }
  return nodes;
}

// Reduce the minimized simulation cell to its primitive cell (SymmetryKit) and retarget the given System in
// place to that primitive cell. The intramolecular topology is re-derived from the retained framework
// definitions on the reduced atoms, and the configured van der Waals cutoff is kept: periodic-image replicas
// satisfy the minimum-image convention on the (small) primitive cell, while Ewald falls back to a half-box
// real-space cutoff. Returns false and leaves \p system unchanged when there is no flexible framework, no
// retained definition, or the symmetry/primitive reduction fails, in which case the caller keeps the
// simulation cell. This is called only after the elastic constants and normal modes have been computed, so
// consuming (overwriting) the minimized conventional cell is safe.
bool reduceSystemToPrimitiveFramework(System& system, double symmetryPrecision)
{
  if (!system.framework.has_value() || system.framework->rigid) return false;
  const Framework& framework = *system.framework;
  if (framework.frameworkDefinitionName.empty()) return false;

  const SimulationBox& box = system.simulationBox;
  const std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  if (frameworkAtoms.empty()) return false;

  std::vector<std::tuple<double3, std::size_t, double>> atoms;
  atoms.reserve(frameworkAtoms.size());
  std::unordered_map<std::size_t, double> chargeForType;
  for (const Atom& atom : frameworkAtoms)
  {
    atoms.emplace_back(box.inverseCell * atom.position, static_cast<std::size_t>(atom.type), 1.0);
    chargeForType.try_emplace(static_cast<std::size_t>(atom.type), atom.charge);
  }

  std::optional<SKSpaceGroup::FoundPrimitiveCellInfo> primitive =
      SKSpaceGroup::SKFindPrimitive(box.cell, atoms, /*allowPartialOccupancies=*/true, symmetryPrecision);
  if (!primitive || primitive->atoms.empty()) return false;

  const double3x3 primitiveCell = primitive->cell.unitCell();
  const SimulationBox primitiveBox(primitiveCell);

  std::vector<Atom> primitiveFractionalAtoms;
  primitiveFractionalAtoms.reserve(primitive->atoms.size());
  for (const auto& [fractional, type, occupancy] : primitive->atoms)
  {
    const auto charge = chargeForType.find(type);
    primitiveFractionalAtoms.emplace_back(fractional, charge != chargeForType.end() ? charge->second : 0.0, 1.0, 0,
                                          static_cast<std::uint16_t>(type), 0, 0, static_cast<std::uint8_t>(true));
  }

  try
  {
    Framework primitiveFramework(system.forceField, framework.name + "-primitive", primitiveBox, 1,
                                 primitiveFractionalAtoms, primitiveFractionalAtoms, int3(1, 1, 1));
    // Re-derive connectivity and all intramolecular potentials (bonds/bends/torsions/vdw/coulomb) on the
    // reduced atoms directly from the parent framework's retained type-based definitions and exclusion/
    // scaling options, i.e. without re-reading the definition file. All fallible work happens here on local
    // objects, so `system` is left untouched if the primitive framework cannot be constructed.
    primitiveFramework.rigid = framework.rigid;
    primitiveFramework.frameworkDefinitionName = framework.frameworkDefinitionName;
    primitiveFramework.bondDefinitions = framework.bondDefinitions;
    primitiveFramework.bendDefinitions = framework.bendDefinitions;
    primitiveFramework.torsionDefinitions = framework.torsionDefinitions;
    primitiveFramework.improperTorsionDefinitions = framework.improperTorsionDefinitions;
    primitiveFramework.excludeIntra12Interactions = framework.excludeIntra12Interactions;
    primitiveFramework.excludeIntra13Interactions = framework.excludeIntra13Interactions;
    primitiveFramework.excludeIntraBondInteractions = framework.excludeIntraBondInteractions;
    primitiveFramework.excludeIntraBendInteractions = framework.excludeIntraBendInteractions;
    primitiveFramework.intra14VanDerWaalsScaling = framework.intra14VanDerWaalsScaling;
    primitiveFramework.intra14ChargeChargeScaling = framework.intra14ChargeChargeScaling;
    primitiveFramework.generateIntraMolecularPotentials(system.forceField);

    // Retarget the minimized system in place to the primitive cell. Keep the (finalized) van der Waals cutoff
    // so the force constants are computed at the same range; periodic-image replicas satisfy the minimum-image
    // convention on the small primitive cell. Only the real-space Coulomb cutoff is clamped to the primitive
    // half-box (Ewald is exact for any box at a half-box real-space cutoff).
    system.forceField.cutOffFrameworkVDWAutomatic = false;
    system.forceField.cutOffMoleculeVDWAutomatic = false;
    const double3 primitiveWidths = primitiveBox.perpendicularWidths();
    const double halfSmallestWidth =
        0.5 * std::min({primitiveWidths.x, primitiveWidths.y, primitiveWidths.z}) - std::numeric_limits<double>::epsilon();
    system.forceField.cutOffCoulomb = std::min(system.forceField.cutOffCoulomb, halfSmallestWidth);

    system.rebuildForFramework(primitiveFramework, primitiveBox);
    system.framework->regenerateVanDerWaalsImageList(system.forceField, primitiveBox,
                                                     system.forceField.cutOffFrameworkVDW);
    return true;
  }
  catch (const std::exception&)
  {
    return false;
  }
}

}  // namespace

Minimization::Minimization(InputReader& inputReader)
    : options(inputReader.minimizationOptions),
      numberOfPreInitializationCycles(inputReader.numberOfPreInitializationCycles),
      numberOfInitializationCycles(inputReader.numberOfInitializationCycles),
      numberOfEquilibrationCycles(inputReader.numberOfEquilibrationCycles),
      printEvery(inputReader.printEvery),
      writeBinaryRestartEvery(inputReader.writeBinaryRestartEvery),
      rescaleWangLandauEvery(inputReader.rescaleWangLandauEvery),
      optimizeMCMovesEvery(inputReader.optimizeMCMovesEvery),
      randomSeed(inputReader.randomSeed),
      random(RandomNumber(inputReader.randomSeed)),
      systems(std::move(inputReader.systems))
{
}

Minimization::Minimization(const MinimizationOptions& options, std::vector<System> systems, bool outputToFiles)
    : options(options), outputToFiles(outputToFiles), systems(std::move(systems))
{
}

void Minimization::setup()
{
  results.assign(systems.size(), {});
  if (outputToFiles)
  {
    std::filesystem::create_directories("output");
    for (std::size_t systemIndex = 0; systemIndex < systems.size(); ++systemIndex)
    {
      streams.emplace_back(std::format("output/minimization.s{}.txt", systemIndex), std::ios::out);
    }
  }

  for (std::size_t systemIndex = 0; systemIndex < systems.size(); ++systemIndex)
  {
    System& system = systems[systemIndex];
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);
    // For cells smaller than twice the van der Waals cutoff a single minimum image is insufficient;
    // regenerate the framework van der Waals pair list with periodic-image replicas so the minimization
    // energy, gradient and Hessian (and the downstream phonon force constants) remain exact.
    if (system.framework.has_value() && !system.framework->rigid)
    {
      system.framework->regenerateVanDerWaalsImageList(system.forceField, system.simulationBox,
                                                       system.forceField.cutOffFrameworkVDW);
    }
    system.precomputeTotalRigidEnergy();
    if (outputToFiles)
    {
      std::print(streams[systemIndex], "{}", system.writeOutputHeader());
      std::print(streams[systemIndex], "{}\n", system.writeSystemStatus());
      std::print(streams[systemIndex],
                 "Baker minimization: maxSteps={} maxStep={} maxCellStep={} rmsTolerance={} maxTolerance={} "
                 "minEigenvalue={} computeElasticConstants={} elasticEigenvalueTolerance={} computeNormalModes={} "
                 "normalModeMovies={} normalModeMoviePeriods={} normalModeMoviePointsPerPeriod={} "
                 "normalModeMovieAmplitude={} computePhononDispersion={} phononDispersionPointsPerSegment={}\n\n",
                 options.maximumNumberOfSteps, options.maximumStepLength, options.maximumCellStepLength,
                 options.rmsGradientTolerance, options.maxGradientTolerance, options.minimumEigenvalue,
                 options.computeElasticConstants, options.elasticEigenvalueTolerance, options.computeNormalModes,
                 options.normalModeMovies, options.normalModeMoviePeriods, options.normalModeMoviePointsPerPeriod,
                 options.normalModeMovieAmplitude, options.computePhononDispersion,
                 options.phononDispersionPointsPerSegment);
    }
  }
}

void Minimization::performCycle()
{
  std::size_t totalNumberOfMolecules{0uz};
  std::size_t totalNumberOfComponents{0uz};
  std::size_t numberOfStepsPerCycle{0uz};

  totalNumberOfMolecules = std::transform_reduce(
      systems.begin(), systems.end(), 0uz, [](const std::size_t& acc, const std::size_t& b) { return acc + b; },
      [](const System& system) { return system.numberOfMolecules(); });
  totalNumberOfComponents = systems.front().numerOfAdsorbateComponents();

  numberOfStepsPerCycle = std::max(totalNumberOfMolecules, 20uz) * totalNumberOfComponents;

  for (std::size_t j = 0uz; j != numberOfStepsPerCycle; j++)
  {
    std::pair<std::size_t, std::size_t> selectedSystemPair = random.randomPairAdjacentIntegers(systems.size());
    System& selectedSystem = systems[selectedSystemPair.first];
    System& selectedSecondSystem = systems[selectedSystemPair.second];

    std::size_t selectedComponent = selectedSystem.randomComponent(random);

    switch (simulationStage)
    {
      case SimulationStage::Uninitialized:
        break;
      case SimulationStage::PreInitialization:
        MC_Moves::performRandomMovePreInitialization(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                     fractionalMoleculeSystem);
        break;
      case SimulationStage::Initialization:
        MC_Moves::performRandomMoveInitialization(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                  fractionalMoleculeSystem);
        break;
      case SimulationStage::Equilibration:
        MC_Moves::performRandomMoveEquilibration(random, selectedSystem, selectedSecondSystem, selectedComponent,
                                                 fractionalMoleculeSystem);

        selectedSystem.components[selectedComponent].lambdaGC.WangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample, selectedSystem.containsTheFractionalMolecule);
        selectedSecondSystem.components[selectedComponent].lambdaGC.WangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample,
            selectedSecondSystem.containsTheFractionalMolecule);

        if (selectedSystem.usesGibbsConventionalCFCMC())
        {
          selectedSystem.components[selectedComponent].lambdaGibbs.WangLandauIteration(
              PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample, true);
          selectedSecondSystem.components[selectedComponent].lambdaGibbs.WangLandauIteration(
              PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample, true);
        }

        selectedSystem.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        selectedSecondSystem.pairSwapLambdaWangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);

        selectedSystem.reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        selectedSecondSystem.reactionLambdaWangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample);
        break;
      case SimulationStage::Run:
        break;
    }

    selectedSystem.components[selectedComponent].lambdaGC.sampleOccupancy(selectedSystem.containsTheFractionalMolecule);
    selectedSecondSystem.components[selectedComponent].lambdaGC.sampleOccupancy(
        selectedSecondSystem.containsTheFractionalMolecule);
    if (selectedSystem.usesGibbsConventionalCFCMC())
    {
      selectedSystem.components[selectedComponent].lambdaGibbs.sampleOccupancy(true);
      selectedSecondSystem.components[selectedComponent].lambdaGibbs.sampleOccupancy(true);
    }
    selectedSystem.pairSwapLambdaSampleOccupancy();
    selectedSecondSystem.pairSwapLambdaSampleOccupancy();
    selectedSystem.reactionLambdaSampleOccupancy();
    selectedSecondSystem.reactionLambdaSampleOccupancy();
  }
}

void Minimization::preInitialize()
{
  if (simulationStage == SimulationStage::PreInitialization) goto continuePreInitializationStage;
  simulationStage = SimulationStage::PreInitialization;

  for (currentCycle = 0uz; currentCycle != numberOfPreInitializationCycles; ++currentCycle, ++absoluteCurrentCycle)
  {
    performCycle();

    for (System& system : systems)
    {
      system.samplePropertiesEvolution(absoluteCurrentCycle);
    }

    if (currentCycle % printEvery == 0uz && outputToFiles)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        std::print(streams[system_id], "{}",
                   system.writePreInitializationStatusReport(currentCycle, numberOfPreInitializationCycles));
        std::flush(streams[system_id]);
        ++system_id;
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // minimization driver does not currently write binary restarts
    }

  continuePreInitializationStage:;
  }
}

void Minimization::initialize()
{
  if (simulationStage == SimulationStage::Initialization) goto continueInitializationStage;
  simulationStage = SimulationStage::Initialization;

  for (currentCycle = 0uz; currentCycle != numberOfInitializationCycles; ++currentCycle, ++absoluteCurrentCycle)
  {
    performCycle();

    for (System& system : systems)
    {
      system.samplePropertiesEvolution(absoluteCurrentCycle);
    }

    if (currentCycle % printEvery == 0uz && outputToFiles)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        std::print(streams[system_id], "{}",
                   system.writeInitializationStatusReport(currentCycle, numberOfInitializationCycles));
        std::flush(streams[system_id]);
        ++system_id;
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // minimization driver does not currently write binary restarts
    }

  continueInitializationStage:;
  }
}

void Minimization::equilibrate()
{
  if (simulationStage == SimulationStage::Equilibration) goto continueEquilibrationStage;
  simulationStage = SimulationStage::Equilibration;

  for (System& system : systems)
  {
    system.runningEnergies = system.computeTotalEnergies();
    for (Component& component : system.components)
    {
      component.lambdaGC.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize,
                                             system.containsTheFractionalMolecule);
      component.lambdaGC.clear();
      if (system.usesGibbsConventionalCFCMC())
      {
        component.lambdaGibbs.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize,
                                                  true);
        component.lambdaGibbs.clear();
      }
    }
    system.pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize);
    system.pairSwapLambdaClearBookkeeping();
    system.reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize);
    system.reactionLambdaClearBookkeeping();
  }

  for (currentCycle = 0uz; currentCycle != numberOfEquilibrationCycles; ++currentCycle, ++absoluteCurrentCycle)
  {
    performCycle();

    for (System& system : systems)
    {
      system.samplePropertiesEvolution(absoluteCurrentCycle);
    }

    if (currentCycle % printEvery == 0uz && outputToFiles)
    {
      for (std::size_t system_id{0}; System& system : systems)
      {
        std::print(streams[system_id], "{}",
                   system.writeEquilibrationStatusReportMC(currentCycle, numberOfEquilibrationCycles));
        std::flush(streams[system_id]);
        ++system_id;
      }
    }

    if (currentCycle % optimizeMCMovesEvery == 0uz)
    {
      for (System& system : systems)
      {
        system.optimizeMCMoves();
      }
    }

    if (currentCycle % rescaleWangLandauEvery == 0uz)
    {
      for (System& system : systems)
      {
        for (Component& component : system.components)
        {
          component.lambdaGC.WangLandauIteration(
              PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors,
              system.containsTheFractionalMolecule);
          if (system.usesGibbsConventionalCFCMC())
          {
            component.lambdaGibbs.WangLandauIteration(
                PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors, true);
          }
        }
        system.pairSwapLambdaWangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors);
        system.reactionLambdaWangLandauIteration(
            PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors);
      }
    }

    if (currentCycle % writeBinaryRestartEvery == 0uz)
    {
      // minimization driver does not currently write binary restarts
    }

  continueEquilibrationStage:;
  }
}

void Minimization::runPhase()
{
  if (simulationStage == SimulationStage::Run) return;
  simulationStage = SimulationStage::Run;

  bool allConverged = true;

  for (std::size_t systemIndex = 0; systemIndex < systems.size(); ++systemIndex)
  {
    System& system = systems[systemIndex];
    MinimizationSystemResult& systemResult = results[systemIndex];
    const std::size_t numberOfFlexibleFrameworkAtoms =
        system.framework && !system.framework->rigid ? system.spanOfFrameworkAtoms().size() : 0;
    const CellMinimizationLayout cellLayout =
        makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
    if (!cellLayout.empty() && system.hasExternalField)
    {
      throw std::runtime_error("Variable-cell minimization does not support an external field");
    }
    const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components,
                                                                    numberOfFlexibleFrameworkAtoms, cellLayout.size());
    // Seed the authoritative rigid-group state once (single Procrustes fit); afterwards
    // applyGeneralizedDisplacement advances the stored quaternions directly.
    system.initializeGroupData();
    GeneralizedHessian hessian(layout.numDofs(), 0);
    std::vector<double> gradient(layout.numDofs(), 0.0);

    DerivativeCapabilities capabilities{};
    capabilities.energy = true;
    capabilities.gradient = true;
    capabilities.hessianPositionPosition = true;

    DerivativeResults derivatives{.gradient = gradient, .hessian = hessian};
    for (std::size_t iteration = 0; iteration < options.maximumNumberOfSteps; ++iteration)
    {
      evaluateDerivatives(system, layout, capabilities, derivatives);
      if (!std::isfinite(derivatives.energy))
      {
        throw std::runtime_error(std::format("Minimization system {} produced a non-finite energy", systemIndex));
      }
      bool negativeModeEscape = false;
      BakerStep step = computeBakerStep(hessian.positionPosition(), gradient, options);
      negativeModeEscape = !cellLayout.empty() && useNegativeModeEscape(step, gradient, options);
      clipCellStep(step, layout, cellLayout, options.maximumCellStepLength);
      if (iteration == 0)
      {
        systemResult.initialEnergy = derivatives.energy;
      }
      systemResult.iterations = iteration;
      systemResult.finalEnergy = derivatives.energy;
      systemResult.rmsGradient = step.rmsGradient;
      systemResult.maxGradient = step.maxGradient;
      systemResult.negativeModes = step.negativeModes;
      systemResult.zeroModes = step.zeroModes;
      systemResult.converged = step.converged;

      const double currentVolume = system.simulationBox.volume;
      system.updateSamplePDBMovie(systemIndex, iteration);

      double acceptedScale = 1.0;
      if (!step.converged && iteration + 1 < options.maximumNumberOfSteps)
      {
        applyGeneralizedDisplacement(system, layout, step.displacement);
        step.stepNorm *= acceptedScale;
      }

      if (iteration % std::max<std::size_t>(1, options.printEvery) == 0 || step.converged)
      {
        const std::string status = std::format(
            "Minimization s{} step {:6d}: E={: .12e} rms(g)={:.5e} max(g)={:.5e} "
            "negative={} zero={} |dx|={:.5e} alpha={:.5e} V={:.8e}{}\n",
            systemIndex, iteration, derivatives.energy, step.rmsGradient, step.maxGradient, step.negativeModes,
            step.zeroModes, step.stepNorm, acceptedScale, currentVolume,
            negativeModeEscape ? " negative-mode escape" : "");
        std::cout << status;
        if (outputToFiles)
        {
          streams[systemIndex] << status;
        }
      }

      if (step.converged)
      {
        break;
      }
    }
    if (!systemResult.converged)
    {
      allConverged = false;
    }
    else
    {
      if (options.computeElasticConstants)
      {
        systemResult.elasticConstants = computeElasticConstants(system, options.elasticEigenvalueTolerance);
      }
      if (options.computeNormalModes || options.normalModeMovies)
      {
        systemResult.normalModes = computeNormalModes(system);
        if (options.normalModeMovies && outputToFiles)
        {
          writeNormalModeMovies(system, *systemResult.normalModes, systemIndex, options.normalModeMoviePeriods,
                                options.normalModeMoviePointsPerPeriod, options.normalModeMovieAmplitude);
        }
      }
      if (options.computePhononDispersion || options.computePhononDensityOfStates)
      {
        // The phonon post-processing is the final step: the elastic constants and normal modes above are
        // already computed and stored, so the minimized conventional cell is no longer needed and can be
        // consumed here. Optionally reduce it in place to the primitive cell so both the dispersion and the
        // density of states describe the true unfolded spectrum over the primitive Brillouin zone (no zone
        // folding). The reduction is a no-op (system left unchanged) when it is unavailable (rigid/no
        // framework, no retained definition, or symmetry detection fails). An explicit dispersion path is
        // given in the simulation-cell reciprocal basis, so it is only meaningful on the simulation cell;
        // primitive reduction is used only together with the default symmetry-derived path.
        if (options.phononUsePrimitiveCell && options.phononDispersionPath.empty())
        {
          reduceSystemToPrimitiveFramework(system, options.phononPrimitiveCellSymmetryPrecision);
        }

        if (options.computePhononDispersion)
        {
          std::vector<PhononPathNode> nodes = options.phononDispersionPath;
          if (nodes.empty())
          {
            // Prefer a symmetry-aware HPKOT path derived from the framework's space group.
            nodes = symmetryDefaultPhononPath(system);
          }
          if (nodes.empty())
          {
            // Fall back to a connected G-X-G-Y-G-Z star along the reciprocal axes.
            nodes = {{double3(0.0, 0.0, 0.0), "G"}, {double3(0.5, 0.0, 0.0), "X"}, {double3(0.0, 0.0, 0.0), "G"},
                     {double3(0.0, 0.5, 0.0), "Y"}, {double3(0.0, 0.0, 0.0), "G"}, {double3(0.0, 0.0, 0.5), "Z"}};
          }
          systemResult.phononDispersion =
              computePhononDispersionAlongPath(system, nodes, options.phononDispersionPointsPerSegment);
        }

        if (options.computePhononDensityOfStates)
        {
          systemResult.phononDensityOfStates = computePhononDensityOfStates(
              system, options.phononDOSMesh, options.phononDOSNumberOfBins, options.phononDOSBroadening);
        }
      }
    }
  }

  if (!allConverged)
  {
    throw std::runtime_error("Baker minimization did not converge for all systems");
  }
}

void Minimization::run()
{
  // Semi-flexible molecules built from rigid/flexible 'Groups' are minimized with per-group rigid-body
  // degrees of freedom (six per rigid group) plus Cartesian degrees of freedom for the flexible atoms.
  // Variable-cell minimization is supported: the cell strain drives the rigid-group centers of mass
  // (the internal group geometry does not scale) and the flexible atoms individually. Electrostatics
  // is supported: the Ewald reciprocal-space and exclusion Hessians resolve rigid bodies per atom
  // (whole molecules and rigid groups alike). Higher-order cross terms (bond-bond, bond-bend,
  // bend-bend, bond/bend-torsion, inversion and out-of-plane bends) are supported as well: their
  // dense Cartesian per-term Hessian is projected onto the generalized (flexible Cartesian +
  // rigid-body) degrees of freedom.
  switch (simulationStage)
  {
    case SimulationStage::Uninitialized:
      setup();
      break;
    case SimulationStage::PreInitialization:
      goto continuePreInitializationStage;
    case SimulationStage::Initialization:
      goto continueInitializationStage;
    case SimulationStage::Equilibration:
      goto continueEquilibrationStage;
    case SimulationStage::Run:
      goto continueRunStage;
    default:
      break;
  }

continuePreInitializationStage:
  preInitialize();
continueInitializationStage:
  initialize();
continueEquilibrationStage:
  equilibrate();
continueRunStage:
  runPhase();

  tearDown();
}

void Minimization::output()
{
  const auto writeArray = [](std::ofstream& stream, std::span<const double> values, double conversion)
  {
    std::print(stream, "[");
    for (std::size_t index = 0; index < values.size(); ++index)
    {
      std::print(stream, "{:.17g}{}", conversion * values[index], index + 1 == values.size() ? "" : ", ");
    }
    std::print(stream, "]");
  };

  for (std::size_t systemIndex = 0; systemIndex < systems.size(); ++systemIndex)
  {
    const System& system = systems[systemIndex];
    const MinimizationSystemResult& result = results[systemIndex];
    std::ofstream json(std::format("output/minimization.s{}.json", systemIndex), std::ios::out);
    std::print(json,
               "{{\n  \"converged\": {},\n  \"iterations\": {},\n  \"initialEnergy\": {:.17g},\n"
               "  \"finalEnergy\": {:.17g},\n  \"rmsGradient\": {:.17g},\n  \"maxGradient\": {:.17g},\n"
               "  \"negativeModes\": {},\n  \"zeroModes\": {},\n  \"atomPositions\": [\n",
               result.converged, result.iterations, result.initialEnergy, result.finalEnergy, result.rmsGradient,
               result.maxGradient, result.negativeModes, result.zeroModes);
    const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    for (std::size_t atom = 0; atom < atoms.size(); ++atom)
    {
      std::print(json, "    [{:.17g}, {:.17g}, {:.17g}]{}", atoms[atom].position.x, atoms[atom].position.y,
                 atoms[atom].position.z, atom + 1 == atoms.size() ? "" : ",");
      std::print(json, "\n");
    }
    std::print(json, "  ]");
    if (result.elasticConstants)
    {
      const ElasticConstantsResult& elastic = *result.elasticConstants;
      const bool reduced = Units::unitSystem == Units::System::ReducedUnits;
      const double pressureConversion = reduced ? 1.0 : 1.0e-9 * Units::PressureConversionFactor;
      const double complianceConversion = reduced ? 1.0 : 1.0 / pressureConversion;
      std::print(json, ",\n  \"elasticConstants\": {{\n");
      std::print(json, "    \"voigtOrder\": [\"xx\", \"yy\", \"zz\", \"yz\", \"xz\", \"xy\"],\n");
      std::print(json, "    \"pressureUnit\": \"{}\",\n", reduced ? "reduced" : "GPa");
      std::print(json, "    \"born\": ");
      writeArray(json, elastic.born, pressureConversion);
      std::print(json, ",\n    \"relaxation\": ");
      writeArray(json, elastic.relaxation, pressureConversion);
      std::print(json, ",\n    \"pressureCorrection\": ");
      writeArray(json, elastic.pressureCorrection, pressureConversion);
      std::print(json, ",\n    \"stiffness\": ");
      writeArray(json, elastic.stiffness, pressureConversion);
      std::print(json, ",\n    \"stabilityEigenvalues\": ");
      writeArray(json, elastic.stabilityEigenvalues, pressureConversion);
      std::print(json, ",\n    \"discardedInternalModes\": {},\n    \"complianceAvailable\": {}",
                 elastic.discardedInternalModes, elastic.complianceAvailable);
      if (elastic.complianceAvailable)
      {
        std::print(json, ",\n    \"complianceUnit\": \"{}\",\n    \"compliance\": ",
                   reduced ? "inverse reduced pressure" : "GPa^-1");
        writeArray(json, elastic.compliance, complianceConversion);
        std::print(json, ",\n    \"youngModuli\": ");
        writeArray(json, elastic.youngModuli, pressureConversion);
        std::print(json, ",\n    \"poissonRatios\": ");
        writeArray(json, elastic.poissonRatios, 1.0);
        std::print(json,
                   ",\n    \"bulkModulus\": {{\"voigt\": {:.17g}, \"reuss\": {:.17g}, \"hill\": {:.17g}}},"
                   "\n    \"shearModulus\": {{\"voigt\": {:.17g}, \"reuss\": {:.17g}, \"hill\": {:.17g}}}",
                   pressureConversion * elastic.bulkModulusVoigt, pressureConversion * elastic.bulkModulusReuss,
                   pressureConversion * elastic.bulkModulusHill, pressureConversion * elastic.shearModulusVoigt,
                   pressureConversion * elastic.shearModulusReuss, pressureConversion * elastic.shearModulusHill);
      }
      std::print(json, "\n  }}");
    }
    if (result.normalModes)
    {
      const NormalModesResult& modes = *result.normalModes;
      const bool reduced = Units::unitSystem == Units::System::ReducedUnits;
      std::print(json, ",\n  \"normalModes\": {{\n");
      std::print(json,
                 "    \"numberOfModes\": {},\n    \"negativeModes\": {},\n    \"zeroModes\": {},\n"
                 "    \"discardedRotationalDofs\": {},\n",
                 modes.numberOfModes, modes.negativeModes, modes.zeroModes, modes.discardedRotationalDofs);
      std::print(json, "    \"eigenvalueUnit\": \"{}\",\n    \"eigenvalues\": ",
                 reduced ? "reduced omega^2" : "ps^-2");
      writeArray(json, modes.eigenvalues, 1.0);
      std::print(json, ",\n    \"frequencyUnit\": \"{}\",\n    \"frequencies\": ", reduced ? "reduced" : "cm^-1");
      writeArray(json, normalModeFrequencies(modes), 1.0);
      std::print(json, "\n  }}");
    }
    if (result.phononDispersion)
    {
      const PhononDispersionResult& dispersion = *result.phononDispersion;
      const bool reduced = Units::unitSystem == Units::System::ReducedUnits;
      const std::size_t numberOfBands =
          dispersion.modes.empty() ? 0 : dispersion.modes.front().eigenvalues.size();
      std::print(json, ",\n  \"phononDispersion\": {{\n");
      std::print(json,
                 "    \"numberOfKPoints\": {},\n    \"numberOfBands\": {},\n"
                 "    \"eigenvalueUnit\": \"{}\",\n    \"frequencyUnit\": \"{}\",\n    \"kPoints\": [\n",
                 dispersion.path.size(), numberOfBands, reduced ? "reduced omega^2" : "ps^-2",
                 reduced ? "reduced" : "cm^-1");
      for (std::size_t index = 0; index < dispersion.path.size(); ++index)
      {
        const PhononKPoint& point = dispersion.path[index];
        std::print(json,
                   "      {{\"label\": \"{}\", \"kFractional\": [{:.17g}, {:.17g}, {:.17g}], "
                   "\"pathCoordinate\": {:.17g},\n        \"eigenvalues\": ",
                   point.label, point.kFractional.x, point.kFractional.y, point.kFractional.z, point.pathCoordinate);
        writeArray(json, dispersion.modes[index].eigenvalues, 1.0);
        std::print(json, ",\n        \"frequencies\": ");
        writeArray(json, phononFrequenciesWavenumber(dispersion.modes[index].eigenvalues), 1.0);
        std::print(json, "}}{}\n", index + 1 == dispersion.path.size() ? "" : ",");
      }
      std::print(json, "    ]\n  }}");

      // Gnuplot-friendly band file: one row per k-point, "pathCoordinate band0 band1 ...".
      std::ofstream bands(std::format("output/phonon_dispersion.s{}.txt", systemIndex), std::ios::out);
      std::print(bands, "# phonon dispersion, frequencies in {}\n", reduced ? "reduced units" : "cm^-1");
      std::print(bands, "# column 1: path coordinate, columns 2..{}: bands; labeled k-points annotated with #\n",
                 numberOfBands + 1);
      for (std::size_t index = 0; index < dispersion.path.size(); ++index)
      {
        const PhononKPoint& point = dispersion.path[index];
        std::print(bands, "{:.17g}", point.pathCoordinate);
        for (const double frequency : phononFrequenciesWavenumber(dispersion.modes[index].eigenvalues))
        {
          std::print(bands, " {:.17g}", frequency);
        }
        if (!point.label.empty()) std::print(bands, "  # {}", point.label);
        std::print(bands, "\n");
      }
    }
    if (result.phononDensityOfStates)
    {
      const PhononDensityOfStates& dos = *result.phononDensityOfStates;
      const bool reduced = Units::unitSystem == Units::System::ReducedUnits;
      std::print(json, ",\n  \"phononDensityOfStates\": {{\n");
      std::print(json,
                 "    \"mesh\": [{}, {}, {}],\n    \"numberOfQPoints\": {},\n"
                 "    \"frequencyUnit\": \"{}\",\n    \"frequency\": ",
                 dos.mesh.x, dos.mesh.y, dos.mesh.z, dos.numberOfQPoints, reduced ? "reduced" : "cm^-1");
      writeArray(json, dos.frequency, 1.0);
      std::print(json, ",\n    \"dos\": ");
      writeArray(json, dos.dos, 1.0);
      std::print(json, "\n  }}");

      // Gnuplot-friendly two-column file: "frequency dos".
      std::ofstream dosStream(std::format("output/phonon_dos.s{}.txt", systemIndex), std::ios::out);
      std::print(dosStream, "# phonon density of states, frequencies in {}\n", reduced ? "reduced units" : "cm^-1");
      std::print(dosStream, "# q-mesh {} x {} x {} ({} q-points); integral of column 2 = number of branches (3N)\n",
                 dos.mesh.x, dos.mesh.y, dos.mesh.z, dos.numberOfQPoints);
      std::print(dosStream, "# column 1: frequency, column 2: DOS\n");
      for (std::size_t index = 0; index < dos.frequency.size(); ++index)
      {
        std::print(dosStream, "{:.17g} {:.17g}\n", dos.frequency[index], dos.dos[index]);
      }
    }
    std::print(json, "\n}}\n");
  }
}

void Minimization::tearDown()
{
  if (outputToFiles)
  {
    output();
    for (std::size_t systemIndex = 0; systemIndex < results.size(); ++systemIndex)
    {
      const MinimizationSystemResult& result = results[systemIndex];
      std::print(streams[systemIndex],
                 "\nFinal minimization status: converged={} iterations={} initialEnergy={:.12e} "
                 "finalEnergy={:.12e} rmsGradient={:.5e} maxGradient={:.5e}\n",
                 result.converged, result.iterations, result.initialEnergy, result.finalEnergy, result.rmsGradient,
                 result.maxGradient);
      if (result.elasticConstants)
      {
        std::print(streams[systemIndex], "{}", writeElasticConstants(*result.elasticConstants));
      }
      if (result.normalModes)
      {
        std::print(streams[systemIndex], "{}", writeNormalModes(*result.normalModes));
      }
      if (result.phononDispersion)
      {
        std::print(streams[systemIndex], "{}", writePhononDispersion(*result.phononDispersion));
      }
      if (result.phononDensityOfStates)
      {
        std::print(streams[systemIndex], "{}", writePhononDensityOfStates(*result.phononDensityOfStates));
      }
    }
  }
}
