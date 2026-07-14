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
import generalized_hessian;
import minimization_dof_layout;
import minimization_evaluate_derivatives;

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
    system.precomputeTotalRigidEnergy();
    if (outputToFiles)
    {
      std::print(streams[systemIndex], "{}", system.writeOutputHeader());
      std::print(streams[systemIndex], "{}\n", system.writeSystemStatus());
      std::print(streams[systemIndex],
                 "Baker minimization: maxSteps={} maxStep={} rmsTolerance={} maxTolerance={} minEigenvalue={}\n\n",
                 options.maximumNumberOfSteps, options.maximumStepLength, options.rmsGradientTolerance,
                 options.maxGradientTolerance, options.minimumEigenvalue);
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
    const MinimizationDofLayout layout =
        buildMinimizationDofLayout(system.moleculeData, system.components, numberOfFlexibleFrameworkAtoms);
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
      const BakerStep step = computeBakerStep(hessian.positionPosition(), gradient, options);
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

      if (iteration % std::max<std::size_t>(1, options.printEvery) == 0 || step.converged)
      {
        const std::string status = std::format(
            "Minimization s{} step {:6d}: E={: .12e} rms(g)={:.5e} max(g)={:.5e} "
            "negative={} zero={} |dx|={:.5e}\n",
            systemIndex, iteration, derivatives.energy, step.rmsGradient, step.maxGradient, step.negativeModes,
            step.zeroModes, step.stepNorm);
        std::cout << status;
        if (outputToFiles)
        {
          streams[systemIndex] << status;
        }
      }

      system.updateSamplePDBMovie(systemIndex, iteration);

      if (step.converged)
      {
        break;
      }
      if (iteration + 1 < options.maximumNumberOfSteps)
      {
        applyGeneralizedDisplacement(system, layout, step.displacement);
      }
    }
    if (!systemResult.converged)
    {
      allConverged = false;
    }
  }

  if (!allConverged)
  {
    throw std::runtime_error("Baker minimization did not converge for all systems");
  }
}

void Minimization::run()
{
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
    std::print(json, "  ]\n}}\n");
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
    }
  }
}
