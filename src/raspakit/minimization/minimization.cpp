module;

module minimization;

import std;

import input_reader;
import atom;
import system;
import generalized_hessian;
import minimization_dof_layout;
import minimization_evaluate_derivatives;

Minimization::Minimization(InputReader &inputReader)
    : options(inputReader.minimizationOptions), systems(std::move(inputReader.systems))
{
}

Minimization::Minimization(const MinimizationOptions &options, std::vector<System> systems, bool outputToFiles)
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
    System &system = systems[systemIndex];
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

void Minimization::run()
{
  setup();
  bool allConverged = true;

  for (std::size_t systemIndex = 0; systemIndex < systems.size(); ++systemIndex)
  {
    System &system = systems[systemIndex];
    MinimizationSystemResult &systemResult = results[systemIndex];
    const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
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
        const std::string status =
            std::format("Minimization s{} step {:6d}: E={: .12e} rms(g)={:.5e} max(g)={:.5e} "
                        "negative={} zero={} |dx|={:.5e}\n",
                        systemIndex, iteration, derivatives.energy, step.rmsGradient, step.maxGradient,
                        step.negativeModes, step.zeroModes, step.stepNorm);
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

  tearDown();
  if (!allConverged)
  {
    throw std::runtime_error("Baker minimization did not converge for all systems");
  }
}

void Minimization::output()
{
  for (std::size_t systemIndex = 0; systemIndex < systems.size(); ++systemIndex)
  {
    const System &system = systems[systemIndex];
    const MinimizationSystemResult &result = results[systemIndex];
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
      const MinimizationSystemResult &result = results[systemIndex];
      std::print(streams[systemIndex],
                 "\nFinal minimization status: converged={} iterations={} initialEnergy={:.12e} "
                 "finalEnergy={:.12e} rmsGradient={:.5e} maxGradient={:.5e}\n",
                 result.converged, result.iterations, result.initialEnergy, result.finalEnergy, result.rmsGradient,
                 result.maxGradient);
    }
  }
}
