module;

module system;

import std;

import double3;
import double3x3;
import atom;
import framework;
import component;
import simulationbox;
import forcefield;
import units;
import loading_data;
import averages;
import enthalpy_of_adsorption_data;
import pressure_data;
import energy_status;
import energy_status_inter;
import property_simulationbox;
import average_energy_type;
import property_energy;
import property_pressure;
import property_loading;
import property_enthalpy;
import property_lambda_probability_histogram;
import property_widom;
import property_temperature;
import reaction;
import reactions;
import running_energy;
import move_statistics;
import mc_moves_move_types;
import mc_moves_statistics;
import mc_moves_cputime;
import json;
import integrators;
import integrators_compute;

// System output: status reports and JSON/text formatting for simulation results.

std::string System::writeOutputHeader() const
{
  std::ostringstream stream;

  std::print(stream, "Compiler and run-time data\n");
  std::print(stream, "===============================================================================\n");

#ifdef VERSION
#define QUOTE(str) #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)
  std::print(stream, "RASPA {}\n\n", EXPAND_AND_QUOTE(VERSION));
#endif

  // ThreadPool &pool = ThreadPool::instance();
  // const std::size_t numberOfHelperThreads = pool.getThreadCount();

  // switch(pool.threadingType)
  //{
  //   case ThreadPool::ThreadingType::Serial:
  //     std::print(stream, "Parallelization: Serial, 1 thread\n");
  //     break;
  //   case ThreadPool::ThreadingType::OpenMP:
  //     std::print(stream, "Parallelization: OpenMP, {} threads\n", numberOfHelperThreads + 1);
  //     break;
  //   case ThreadPool::ThreadingType::ThreadPool:
  //     std::print(stream, "Parallelization: ThreadPool, {} threads\n", numberOfHelperThreads + 1);
  //     break;
  //   case ThreadPool::ThreadingType::GPU_Offload:
  //     std::print(stream, "Parallelization: GPU-Offload\n");
  //     break;
  // }
  // std::print(stream, "\n");

  return stream.str();
}

std::string System::writeNumberOfPseudoAtoms() const
{
  std::ostringstream stream;

  std::print(stream, "Number of pseudo-atoms\n");
  std::print(stream, "===============================================================================\n\n");

  for (std::size_t componentId = 0; const Component& c : components)
  {
    std::print(stream, "Component {:3d} ({})\n", componentId, c.name);
    std::print(stream, "-------------------------------------------------------------------------------\n");
    for (std::size_t index = 0; const std::size_t number_of_pseudo_atoms : numberOfPseudoAtoms[componentId])
    {
      std::print(stream, "    index {:3d} ({}): {} atoms\n", index, forceField.pseudoAtoms[index].name,
                 number_of_pseudo_atoms);
      ++index;
    }
    std::print(stream, "\n");
    ++componentId;
  }

  std::print(stream, "Total number of pseudo-atoms:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (std::size_t index = 0; const std::size_t number_of_pseudo_atoms : totalNumberOfPseudoAtoms)
  {
    std::print(stream, "    index {:3d} ({}): {} atoms\n", index, forceField.pseudoAtoms[index].name,
               number_of_pseudo_atoms);
    ++index;
  }

  std::print(stream, "\n\n\n\n");

  return stream.str();
}

std::string System::writeInitializationStatusReport(std::size_t currentCycle, std::size_t numberOfCycles) const
{
  std::ostringstream stream;

  std::print(stream, "Initialization: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::print(stream, "{}\n", simulationBox.printStatus());
  std::print(stream, "Net charge: {:12.8f}\n", netCharge);
  std::print(stream, "{}", forceField.printCutOffAutoStatus());
  std::print(stream, "\n");

  for (std::size_t componentId{0}; const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {:3d} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n",
                 componentId, c.name, lambda, runningEnergies.dudlambda(lambda, c.lambdaGC.dUdlambdaGroupId), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {:3d} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", componentId, c.name,
                 lambda, occupancy, averageOccupancy);
    }
    std::print(stream, "    net charge: {:12.8f} [e]\n", netChargePerComponent[componentId]);
    ++componentId;
  }

  if (usesReactionConventionalCFCMC())
  {
    for (const Reaction& reaction : reactions.list)
    {
      const double lambda = reaction.currentLambda;
      if (reaction.isSerialRxCFC())
      {
        std::print(stream,
                   "reaction {:3d} lambda: {: g} dUdlambda: {: g} fractional side: {} "
                   "occupancy reactants: ({:3f}) products: ({:3f})\n",
                   reaction.id, lambda, reactionDUdlambda(reaction),
                   reaction.fractionalSideIsReactants ? "reactants" : "products", reaction.lambda.occupancy(),
                   reaction.lambdaProductSide.occupancy());
        continue;
      }
      const PropertyLambdaProbabilityHistogram& histogram = activeReactionLambdaHistogram(reaction);
      const double averageOccupancy = histogram.occupancy();
      std::print(stream, "reaction {:3d} lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", reaction.id,
                 lambda, reactionDUdlambda(reaction),
                 static_cast<double>(hasReactionFractionalMolecules()),
                 averageOccupancy);
    }
    std::print(stream, "\n");
  }

  std::print(stream, "Amount of molecules per component:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (std::size_t componentId{0}; const Component& c : components)
  {
    std::print(stream, "{}",
               loadings.printStatus(componentId, c.name, c.totalMass, c.amountOfExcessMolecules, frameworkMass(),
                                    framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));
    ++componentId;
  }
  std::print(stream, "\n");

  stream << runningEnergies.printMC();

  std::print(stream, "\n");

  return stream.str();
}

std::string System::writeEquilibrationStatusReportMC(std::size_t currentCycle, std::size_t numberOfCycles) const
{
  std::ostringstream stream;

  std::print(stream, "Equilibration: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::print(stream, "{}\n", simulationBox.printStatus());
  std::print(stream, "Net charge: {:12.8f}\n", netCharge);
  std::print(stream, "{}", forceField.printCutOffAutoStatus());
  std::print(stream, "\n");

  for (std::size_t componentId = 0; const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", componentId,
                 c.name, lambda, runningEnergies.dudlambda(lambda, c.lambdaGC.dUdlambdaGroupId), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", componentId, c.name,
                 c.lambdaGC.lambdaValue(), occupancy, averageOccupancy);
    }
    std::print(stream, "    net charge: {:12.8f} [e]\n", netChargePerComponent[componentId]);
    ++componentId;
  }

  if (usesReactionConventionalCFCMC())
  {
    for (const Reaction& reaction : reactions.list)
    {
      const double lambda = reaction.currentLambda;
      if (reaction.isSerialRxCFC())
      {
        std::print(stream,
                   "reaction {:3d} lambda: {: g} dUdlambda: {: g} fractional side: {} "
                   "occupancy reactants: ({:3f}) products: ({:3f})\n",
                   reaction.id, lambda, reactionDUdlambda(reaction),
                   reaction.fractionalSideIsReactants ? "reactants" : "products", reaction.lambda.occupancy(),
                   reaction.lambdaProductSide.occupancy());
        continue;
      }
      const PropertyLambdaProbabilityHistogram& histogram = activeReactionLambdaHistogram(reaction);
      const double averageOccupancy = histogram.occupancy();
      std::print(stream, "reaction {:3d} lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", reaction.id,
                 lambda, reactionDUdlambda(reaction),
                 static_cast<double>(hasReactionFractionalMolecules()),
                 averageOccupancy);
    }
    std::print(stream, "\n");
  }

  std::print(stream, "Amount of molecules per component:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (std::size_t componentId{0}; const Component& c : components)
  {
    std::print(stream, "{}",
               loadings.printStatus(componentId, c.name, c.totalMass, c.amountOfExcessMolecules, frameworkMass(),
                                    framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));
    ++componentId;
  }
  std::print(stream, "\n");

  stream << runningEnergies.printMC();

  std::print(stream, "\n");

  return stream.str();
}

std::string System::writeEquilibrationStatusReportMD(std::size_t currentCycle, std::size_t numberOfCycles) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;

  std::print(stream, "Equilibration: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::print(stream, "{}\n", simulationBox.printStatus());
  double3 linear_momentum = Integrators::computeLinearMomentum(moleculeData);
  std::print(stream, "Linear momentum: {:12.8f} {:12.8f} {:12.8f}\n", linear_momentum.x, linear_momentum.y,
             linear_momentum.z);
  double3 com_velocity = Integrators::computeCenterOfMassVelocity(moleculeData);
  std::print(stream, "Center of mass velocity: {:12.8f} {:12.8f} {:12.8f}\n", com_velocity.x, com_velocity.y,
             com_velocity.z);
  double3 com = Integrators::computeCenterOfMass(moleculeData);
  std::print(stream, "Center of mass: {:12.8f} {:12.8f} {:12.8f}\n", com.x, com.y, com.z);
  std::print(stream, "Net charge: {:12.8f}\n", netCharge);
  std::print(stream, "\n");

  double translationalKineticEnergy =
      Integrators::computeTranslationalKineticEnergy(moleculeData, spanOfMoleculeAtoms(), spanOfMoleculeDynamics(),
                                                     components);
  double translationalTemperature =
      2.0 * translationalKineticEnergy /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint));
  double rotationalKineticEnergy = Integrators::computeRotationalKineticEnergy(moleculeData, components);
  double rotationalTemperature =
      rotationalDegreesOfFreedom > 0
          ? 2.0 * rotationalKineticEnergy / (Units::KB * static_cast<double>(rotationalDegreesOfFreedom))
          : 0.0;
  double overallTemperature =
      2.0 * (translationalKineticEnergy + rotationalKineticEnergy) /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint +
                                       rotationalDegreesOfFreedom));
  std::print(stream, "Temperature: {: .6e}\n", overallTemperature);
  std::print(stream, "Translational temperature: {: .6e}\n", translationalTemperature);
  std::print(stream, "Rotational temperature: {: .6e}\n\n", rotationalTemperature);

  std::print(stream, "Translational constraint degrees of freedom center of mass: {}\n",
             translationalCenterOfMassConstraint);
  std::print(stream, "Translational degrees of freedom molecules: {}\n", translationalDegreesOfFreedom);
  std::print(stream, "Total translational degrees of freedom molecules: {}\n",
             translationalDegreesOfFreedom - translationalCenterOfMassConstraint);
  std::print(stream, "Rotational degrees of freedom molecules: {}\n\n", rotationalDegreesOfFreedom);

  std::print(stream, "Potential energy:   {: .6e}\n", conv * runningEnergies.potentialEnergy());
  std::print(stream, "Kinetic energy:     {: .6e}\n",
             conv * (runningEnergies.translationalKineticEnergy + runningEnergies.rotationalKineticEnergy));
  std::print(stream, "Nose-Hoover energy: {: .6e}\n", conv * runningEnergies.NoseHooverEnergy);
  std::print(stream, "Conserved energy:   {: .6e}\n", conv * runningEnergies.conservedEnergy());
  double drift = std::abs(conv * (conservedEnergy - referenceEnergy) / referenceEnergy);
  std::print(stream, "Drift: {:.6e} Average drift: {:.6e}\n\n", drift,
             accumulatedDrift / static_cast<double>(std::max(currentCycle, 1uz)));

  std::print(stream, "\n");

  for (std::size_t componentId{0}; const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", componentId,
                 c.name, lambda, runningEnergies.dudlambda(lambda, c.lambdaGC.dUdlambdaGroupId), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", componentId, c.name,
                 c.lambdaGC.lambdaValue(), occupancy, averageOccupancy);
    }
    std::print(stream, "    net charge: {:12.8f} [e]\n", netChargePerComponent[componentId]);
    ++componentId;
  }

  if (usesReactionConventionalCFCMC())
  {
    for (const Reaction& reaction : reactions.list)
    {
      const double lambda = reaction.currentLambda;
      if (reaction.isSerialRxCFC())
      {
        std::print(stream,
                   "reaction {:3d} lambda: {: g} dUdlambda: {: g} fractional side: {} "
                   "occupancy reactants: ({:3f}) products: ({:3f})\n",
                   reaction.id, lambda, reactionDUdlambda(reaction),
                   reaction.fractionalSideIsReactants ? "reactants" : "products", reaction.lambda.occupancy(),
                   reaction.lambdaProductSide.occupancy());
        continue;
      }
      const PropertyLambdaProbabilityHistogram& histogram = activeReactionLambdaHistogram(reaction);
      const double averageOccupancy = histogram.occupancy();
      std::print(stream, "reaction {:3d} lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", reaction.id,
                 lambda, reactionDUdlambda(reaction),
                 static_cast<double>(hasReactionFractionalMolecules()),
                 averageOccupancy);
    }
    std::print(stream, "\n");
  }

  std::print(stream, "Amount of molecules per component:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (std::size_t componentId{0}; const Component& c : components)
  {
    std::print(stream, "{}",
               loadings.printStatus(componentId, c.name, c.totalMass, c.amountOfExcessMolecules, frameworkMass(),
                                    framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));
    ++componentId;
  }
  std::print(stream, "\n");

  return stream.str();
}

std::string System::writeProductionStatusReportMC(const std::string& statusLine) const
{
  std::ostringstream stream;

  std::print(stream, "{}", statusLine);
  std::print(stream, "===============================================================================\n\n");

  auto [simulation_box, average_simulation_box] = averageSimulationBox.averageSimulationBox();
  std::print(stream, "{}\n", simulationBox.printStatus(simulation_box, average_simulation_box));

  std::print(stream, "Net charge: {:12.8f}\n", netCharge);
  std::print(stream, "{}", forceField.printCutOffAutoStatus());
  std::print(stream, "\n");

  for (std::size_t componentId{0}; const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", componentId,
                 c.name, lambda, runningEnergies.dudlambda(lambda, c.lambdaGC.dUdlambdaGroupId), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", componentId, c.name,
                 c.lambdaGC.lambdaValue(), occupancy, averageOccupancy);
    }
    std::print(stream, "    net charge: {:12.8f} [e]\n", netChargePerComponent[componentId]);
    ++componentId;
  }

  if (usesReactionConventionalCFCMC())
  {
    for (const Reaction& reaction : reactions.list)
    {
      const double lambda = reaction.currentLambda;
      if (reaction.isSerialRxCFC())
      {
        std::print(stream,
                   "reaction {:3d} lambda: {: g} dUdlambda: {: g} fractional side: {} "
                   "occupancy reactants: ({:3f}) products: ({:3f})\n",
                   reaction.id, lambda, reactionDUdlambda(reaction),
                   reaction.fractionalSideIsReactants ? "reactants" : "products", reaction.lambda.occupancy(),
                   reaction.lambdaProductSide.occupancy());
        continue;
      }
      const PropertyLambdaProbabilityHistogram& histogram = activeReactionLambdaHistogram(reaction);
      const double averageOccupancy = histogram.occupancy();
      std::print(stream, "reaction {:3d} lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", reaction.id,
                 lambda, reactionDUdlambda(reaction),
                 static_cast<double>(hasReactionFractionalMolecules()),
                 averageOccupancy);
    }
    std::print(stream, "\n");
  }

  std::print(stream, "Amount of molecules per component:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::pair<LoadingData, LoadingData> loadingData = averageLoadings.result();
  for (std::size_t componentId{0}; const Component& c : components)
  {
    std::print(stream, "{}",
               loadings.printStatus(componentId, c.name, c.totalMass, c.amountOfExcessMolecules,
                                    loadingData.first, loadingData.second, frameworkMass(),
                                    framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));
    ++componentId;
  }
  std::print(stream, "\n");
  double conv = Units::EnergyToKelvin;

  if (!(framework.has_value() && framework->rigid))
  {
    std::pair<PressureData, PressureData> average_pressure = averagePressure.result();

    switch (Units::unitSystem)
    {
      case Units::System::RASPA:
      {
        double3x3 pressureTensor = 1e-5 * Units::PressureConversionFactor * average_pressure.first.totalPressureTensor;
        double3x3 pressureTensorError = 1e-5 * Units::PressureConversionFactor * average_pressure.second.totalPressureTensor;
        std::print(stream, "Average pressure tensor: \n");
        std::print(stream, "-------------------------------------------------------------------------------\n");
        std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.ax,
                   pressureTensor.bx, pressureTensor.cx, pressureTensorError.ax, pressureTensorError.bx,
                   pressureTensorError.cx);
        std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.ay,
                   pressureTensor.by, pressureTensor.cy, pressureTensorError.ay, pressureTensorError.by,
                   pressureTensorError.cy);
        std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.az,
                   pressureTensor.bz, pressureTensor.cz, pressureTensorError.az, pressureTensorError.bz,
                   pressureTensorError.cz);
        std::print(stream, "Ideal-gas pressure:  {: .6e} +/ {:.6e} [bar]\n",
                   1e-5 * Units::PressureConversionFactor * average_pressure.first.idealGasPressure,
                   1e-5 * Units::PressureConversionFactor * average_pressure.second.idealGasPressure);
        std::print(stream, "Excess pressure:     {: .6e} +/ {:.6e} [bar]\n",
                   1e-5 * Units::PressureConversionFactor * average_pressure.first.excessPressure,
                   1e-5 * Units::PressureConversionFactor * average_pressure.second.excessPressure);
        std::print(stream, "Pressure:            {: .6e} +/ {:.6e} [bar]\n\n",
                   1e-5 * Units::PressureConversionFactor * average_pressure.first.totalPressure, 
                   1e-5 * Units::PressureConversionFactor * average_pressure.second.totalPressure);
      }
      break;
      case Units::System::ReducedUnits:
      {
        double3x3 pressureTensor = average_pressure.first.totalPressureTensor;
        double3x3 pressureTensorError = average_pressure.second.totalPressureTensor;
        std::print(stream, "Average pressure tensor: \n");
        std::print(stream, "-------------------------------------------------------------------------------\n");
        std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [{}]\n", pressureTensor.ax,
                   pressureTensor.bx, pressureTensor.cx, pressureTensorError.ax, pressureTensorError.bx,
                   pressureTensorError.cx, Units::unitOfPressureString);
        std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [{}]\n", pressureTensor.ay,
                   pressureTensor.by, pressureTensor.cy, pressureTensorError.ay, pressureTensorError.by,
                   pressureTensorError.cy, Units::unitOfPressureString);
        std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [{}]\n", pressureTensor.az,
                   pressureTensor.bz, pressureTensor.cz, pressureTensorError.az, pressureTensorError.bz,
                   pressureTensorError.cz, Units::unitOfPressureString);
        std::print(stream, "Ideal-gas pressure:  {: .6e} +/ {:.6e} [{}]\n", average_pressure.first.idealGasPressure,
                   average_pressure.second.idealGasPressure, Units::unitOfPressureString);
        std::print(stream, "Excess pressure:     {: .6e} +/ {:.6e} [{}]\n", average_pressure.first.excessPressure, 
                   average_pressure.second.excessPressure,
                   Units::unitOfPressureString);
        std::print(stream, "Pressure:            {: .6e} +/ {:.6e} [{}]\n\n", average_pressure.first.totalPressure, 
                   average_pressure.second.totalPressure, Units::unitOfPressureString);
      }
      break;
    }
  }

  std::pair<EnergyStatus, EnergyStatus> energyData = averageEnergies.result();
  std::print(stream, "Total potential energy{}  {: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString, conv * currentEnergyStatus.totalEnergy.energy,
             conv * energyData.first.totalEnergy.energy, conv * energyData.second.totalEnergy.energy,
             Units::displayedUnitOfEnergyString);
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::print(stream, "ExternalField-molecule\n");
  std::print(stream, "    Van der Waals{}       {: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString,
             conv * currentEnergyStatus.externalFieldMoleculeEnergy.VanDerWaals.energy,
             conv * energyData.first.externalFieldMoleculeEnergy.VanDerWaals.energy,
             conv * energyData.second.externalFieldMoleculeEnergy.VanDerWaals.energy,
             Units::displayedUnitOfEnergyString);
  std::print(stream, "Framework-molecule\n");
  std::print(stream, "    Van der Waals{}       {: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString,
             conv * currentEnergyStatus.frameworkMoleculeEnergy.VanDerWaals.energy,
             conv * energyData.first.frameworkMoleculeEnergy.VanDerWaals.energy,
             conv * energyData.second.frameworkMoleculeEnergy.VanDerWaals.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Van der Waals (Tail){}{: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString,
             conv * currentEnergyStatus.frameworkMoleculeEnergy.VanDerWaalsTailCorrection.energy,
             conv * energyData.first.frameworkMoleculeEnergy.VanDerWaalsTailCorrection.energy,
             conv * energyData.second.frameworkMoleculeEnergy.VanDerWaalsTailCorrection.energy,
             Units::displayedUnitOfEnergyString);
  std::print(stream, "    Coulombic Real{}      {: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString,
             conv * currentEnergyStatus.frameworkMoleculeEnergy.CoulombicReal.energy,
             conv * energyData.first.frameworkMoleculeEnergy.CoulombicReal.energy,
             conv * energyData.second.frameworkMoleculeEnergy.CoulombicReal.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Coulombic Fourier{}   {: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString,
             conv * currentEnergyStatus.frameworkMoleculeEnergy.CoulombicFourier.energy,
             conv * energyData.first.frameworkMoleculeEnergy.CoulombicFourier.energy,
             conv * energyData.second.frameworkMoleculeEnergy.CoulombicFourier.energy,
             Units::displayedUnitOfEnergyString);
  std::print(stream, "Molecule-molecule\n");
  std::print(stream, "    Van der Waals{}       {: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString, conv * currentEnergyStatus.interEnergy.VanDerWaals.energy,
             conv * energyData.first.interEnergy.VanDerWaals.energy,
             conv * energyData.second.interEnergy.VanDerWaals.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Van der Waals (Tail){}{: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString,
             conv * currentEnergyStatus.interEnergy.VanDerWaalsTailCorrection.energy,
             conv * energyData.first.interEnergy.VanDerWaalsTailCorrection.energy,
             conv * energyData.second.interEnergy.VanDerWaalsTailCorrection.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Coulombic Real{}      {: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString, conv * currentEnergyStatus.interEnergy.CoulombicReal.energy,
             conv * energyData.first.interEnergy.CoulombicReal.energy,
             conv * energyData.second.interEnergy.CoulombicReal.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Coulombic Fourier{}   {: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString,
             conv * currentEnergyStatus.interEnergy.CoulombicFourier.energy,
             conv * energyData.first.interEnergy.CoulombicFourier.energy,
             conv * energyData.second.interEnergy.CoulombicFourier.energy, Units::displayedUnitOfEnergyString);
  std::print(stream, "    Molecule Intra{}      {: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString, conv * currentEnergyStatus.intraEnergy.total().energy,
             conv * energyData.first.intraEnergy.total().energy, conv * energyData.second.intraEnergy.total().energy,
             Units::displayedUnitOfEnergyString);
  std::print(stream, "Polarization energy{}     {: .6e} ({: .6e} +/- {:.6e}) [{}]\n",
             Units::displayedUnitOfEnergyConversionString, conv * currentEnergyStatus.polarizationEnergy.energy,
             conv * energyData.first.polarizationEnergy.energy, conv * energyData.second.polarizationEnergy.energy,
             Units::displayedUnitOfEnergyString);

  std::print(stream, "\n");

  return stream.str();
}

std::string System::writeProductionStatusReportMD(std::size_t currentCycle, std::size_t numberOfCycles) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;

  std::print(stream, "Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::pair<SimulationBox, SimulationBox> simulationBoxData = averageSimulationBox.averageSimulationBox();
  std::print(stream, "{}", simulationBox.printStatus(simulationBoxData.first, simulationBoxData.second));
  std::print(stream, "\n");

  double3 linear_momentum = Integrators::computeLinearMomentum(moleculeData);
  std::print(stream, "Linear momentum: {:12.8f} {:12.8f} {:12.8f}\n", linear_momentum.x, linear_momentum.y,
             linear_momentum.z);
  double3 com_velocity = Integrators::computeCenterOfMassVelocity(moleculeData);
  std::print(stream, "Center of mass velocity: {:12.8f} {:12.8f} {:12.8f}\n", com_velocity.x, com_velocity.y,
             com_velocity.z);
  double3 com = Integrators::computeCenterOfMass(moleculeData);
  std::print(stream, "Center of mass: {:12.8f} {:12.8f} {:12.8f}\n", com.x, com.y, com.z);
  std::print(stream, "Net charge: {:12.8f}\n", netCharge);
  std::print(stream, "Time run: {:g} [ps]  {:g} [ns]\n\n", static_cast<double>(currentCycle) * timeStep,
             static_cast<double>(currentCycle) * timeStep / 1000.0);

  double translational_kinetic_energy = Integrators::computeTranslationalKineticEnergy(
      moleculeData, spanOfMoleculeAtoms(), spanOfMoleculeDynamics(), components);
  double translational_temperature =
      2.0 * translational_kinetic_energy /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint));
  double rotational_kinetic_energy = Integrators::computeRotationalKineticEnergy(moleculeData, components);
  double rotational_temperature =
      rotationalDegreesOfFreedom > 0
          ? 2.0 * rotational_kinetic_energy / (Units::KB * static_cast<double>(rotationalDegreesOfFreedom))
          : 0.0;
  double overall_temperature =
      2.0 * (translational_kinetic_energy + rotational_kinetic_energy) /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint +
                                       rotationalDegreesOfFreedom));
  std::pair<double, double> average_temperature = averageTemperature.averageTemperature();
  std::pair<double, double> average_translational_temperature = averageTranslationalTemperature.averageTemperature();
  std::pair<double, double> average_rotational_temperature = averageRotationalTemperature.averageTemperature();

  std::print(stream, "Temperature: {: .6e} ({: .6e} +/- {:.6e})\n", overall_temperature, average_temperature.first,
             average_temperature.second);
  std::print(stream, "Translational temperature: {: .6e} ({: .6e} +/- {:.6e})\n", translational_temperature,
             average_translational_temperature.first, average_translational_temperature.second);
  std::print(stream, "Rotational temperature: {: .6e} ({: .6e} +/- {:.6e})\n\n", rotational_temperature,
             average_rotational_temperature.first, average_rotational_temperature.second);

  std::print(stream, "Translational constraint degrees of freedom center of mass: {}\n",
             translationalCenterOfMassConstraint);
  std::print(stream, "Translational degrees of freedom molecules: {}\n", translationalDegreesOfFreedom);
  std::print(stream, "Total translational degrees of freedom molecules: {}\n",
             translationalDegreesOfFreedom - translationalCenterOfMassConstraint);
  std::print(stream, "Rotational degrees of freedom molecules: {}\n\n", rotationalDegreesOfFreedom);

  std::print(stream, "Potential energy:   {: .6e}\n", conv * runningEnergies.potentialEnergy());
  std::print(stream, "Kinetic energy:     {: .6e}\n",
             conv * (runningEnergies.translationalKineticEnergy + runningEnergies.rotationalKineticEnergy));
  std::print(stream, "Nose-Hoover energy: {: .6e}\n", conv * runningEnergies.NoseHooverEnergy);
  std::print(stream, "Conserved energy:   {: .6e}\n", conv * runningEnergies.conservedEnergy());
  double drift = std::abs(conv * (conservedEnergy - referenceEnergy) / referenceEnergy);
  std::print(stream, "Drift: {:.6e} Average drift: {:.6e}\n\n", drift,
             accumulatedDrift / static_cast<double>(std::max(currentCycle, 1uz)));

  std::pair<EnergyStatus, EnergyStatus> energyData = averageEnergies.result();
  std::print(stream, "Total potential energy:   {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
             conv * currentEnergyStatus.totalEnergy.energy, conv * energyData.first.totalEnergy.energy,
             conv * energyData.second.totalEnergy.energy);
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::print(stream, "ExternalField-molecule\n");
  std::print(stream, "    Van der Waals:        {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
             conv * currentEnergyStatus.externalFieldMoleculeEnergy.VanDerWaals.energy,
             conv * energyData.first.externalFieldMoleculeEnergy.VanDerWaals.energy,
             conv * energyData.second.externalFieldMoleculeEnergy.VanDerWaals.energy);
  std::print(stream, "Framework-molecule\n");
  std::print(stream, "    Van der Waals:        {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
             conv * currentEnergyStatus.frameworkMoleculeEnergy.VanDerWaals.energy,
             conv * energyData.first.frameworkMoleculeEnergy.VanDerWaals.energy,
             conv * energyData.second.frameworkMoleculeEnergy.VanDerWaals.energy);
  std::print(stream, "    Van der Waals (Tail): {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
             conv * currentEnergyStatus.frameworkMoleculeEnergy.VanDerWaalsTailCorrection.energy,
             conv * energyData.first.frameworkMoleculeEnergy.VanDerWaalsTailCorrection.energy,
             conv * energyData.second.frameworkMoleculeEnergy.VanDerWaalsTailCorrection.energy);
  std::print(stream, "    Coulombic Real:       {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
             conv * currentEnergyStatus.frameworkMoleculeEnergy.CoulombicReal.energy,
             conv * energyData.first.frameworkMoleculeEnergy.CoulombicReal.energy,
             conv * energyData.second.frameworkMoleculeEnergy.CoulombicReal.energy);
  std::print(stream, "    Coulombic Fourier:    {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
             conv * currentEnergyStatus.frameworkMoleculeEnergy.CoulombicFourier.energy,
             conv * energyData.first.frameworkMoleculeEnergy.CoulombicFourier.energy,
             conv * energyData.second.frameworkMoleculeEnergy.CoulombicFourier.energy);
  std::print(stream, "Molecule-molecule\n");
  std::print(stream, "    Van der Waals:        {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
             conv * currentEnergyStatus.interEnergy.VanDerWaals.energy,
             conv * energyData.first.interEnergy.VanDerWaals.energy,
             conv * energyData.second.interEnergy.VanDerWaals.energy);
  std::print(stream, "    Van der Waals (Tail): {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
             conv * currentEnergyStatus.interEnergy.VanDerWaalsTailCorrection.energy,
             conv * energyData.first.interEnergy.VanDerWaalsTailCorrection.energy,
             conv * energyData.second.interEnergy.VanDerWaalsTailCorrection.energy);
  std::print(stream, "    Coulombic Real:       {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
             conv * currentEnergyStatus.interEnergy.CoulombicReal.energy,
             conv * energyData.first.interEnergy.CoulombicReal.energy,
             conv * energyData.second.interEnergy.CoulombicReal.energy);
  std::print(stream, "    Coulombic Fourier:    {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
             conv * currentEnergyStatus.interEnergy.CoulombicFourier.energy,
             conv * energyData.first.interEnergy.CoulombicFourier.energy,
             conv * energyData.second.interEnergy.CoulombicFourier.energy);
  std::print(stream, "    Molecule Intra:       {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
             conv * currentEnergyStatus.intraEnergy.total().energy, conv * energyData.first.intraEnergy.total().energy,
             conv * energyData.second.intraEnergy.total().energy);

  std::print(stream, "\n");
  for (std::size_t componentId = 0; const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", componentId,
                 c.name, lambda, runningEnergies.dudlambda(lambda, c.lambdaGC.dUdlambdaGroupId), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", componentId, c.name,
                 lambda, occupancy, averageOccupancy);
    }
    std::print(stream, "    net charge: {:12.8f} [e]\n", netChargePerComponent[componentId]);
    ++componentId;
  }
  std::print(stream, "\n");

  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::pair<LoadingData, LoadingData> loadingData = averageLoadings.result();
  for (std::size_t componentId{0}; const Component& c : components)
  {
    std::print(stream, "{}",
               loadings.printStatus(componentId, c.name, c.totalMass, c.amountOfExcessMolecules,
                                    loadingData.first, loadingData.second, frameworkMass(),
                                    framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));
    ++componentId;
  }
  std::print(stream, "\n");

  std::pair<PressureData, PressureData> average_pressure = averagePressure.result();

  double3x3 pressureTensor = 1e-5 * Units::PressureConversionFactor * average_pressure.first.totalPressureTensor;
  double3x3 pressureTensorError = 1e-5 * Units::PressureConversionFactor * average_pressure.second.totalPressureTensor;

  std::print(stream, "Average pressure tensor: \n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.ax, pressureTensor.bx,
             pressureTensor.cx, pressureTensorError.ax, pressureTensorError.bx, pressureTensorError.cx);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.ay, pressureTensor.by,
             pressureTensor.cy, pressureTensorError.ay, pressureTensorError.by, pressureTensorError.cy);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.az, pressureTensor.bz,
             pressureTensor.cz, pressureTensorError.az, pressureTensorError.bz, pressureTensorError.cz);

  std::print(stream, "Ideal-gas pressure:  {: .6e} +/ {:.6e} [bar]\n",
             1e-5 * Units::PressureConversionFactor * average_pressure.first.idealGasPressure,
             1e-5 * Units::PressureConversionFactor * average_pressure.second.idealGasPressure);
  std::print(stream, "Excess pressure:     {: .6e} +/ {:.6e} [bar]\n",
             1e-5 * Units::PressureConversionFactor * average_pressure.first.excessPressure,
             1e-5 * Units::PressureConversionFactor * average_pressure.second.excessPressure);
  std::print(stream, "Pressure:            {: .6e} +/ {:.6e} [bar]\n\n",
             1e-5 * Units::PressureConversionFactor * average_pressure.first.totalPressure, 
             1e-5 * Units::PressureConversionFactor * average_pressure.second.totalPressure);

  return stream.str();
}

std::string System::writeSystemStatus() const
{
  std::ostringstream stream;

  std::print(stream, "System definitions\n");
  std::print(stream, "===============================================================================\n\n");

  std::print(stream, "Temperature:          {} [{}]\n", temperature, Units::unitOfTemperatureString);
  std::print(stream, "Beta:                 {} [-]\n", beta);
  std::print(stream, "Pressure:             {} [{}]\n", pressure * Units::PressureConversionFactor,
             Units::unitOfPressureString);
  std::print(stream, "Helium void fraction: {} [-]\n\n", heliumVoidFraction);

  stream << simulationBox.printStatus();
  std::print(stream, "\n\n\n");

  std::print(stream, "Property measurement settings\n");
  std::print(stream, "===============================================================================\n\n");
  if (averageEnergyHistogram.has_value())
  {
    stream << averageEnergyHistogram->printSettings();
  }
  if (propertyMoleculeProperties.has_value())
  {
    stream << propertyMoleculeProperties->printSettings();
  }
  std::print(stream, "\n\n\n");

  return stream.str();
}

nlohmann::json System::jsonSystemStatus() const
{
  nlohmann::json system;
  system["temperature"] = temperature;
  system["beta"] = beta;
  system["pressure"] = pressure * Units::PressureConversionFactor;

  system.merge_patch(simulationBox.jsonStatus());
  return system;
}

std::string System::writeComponentStatus() const
{
  std::ostringstream stream;

  std::print(stream, "Component definitions\n");
  std::print(stream, "===============================================================================\n\n");
  if (framework.has_value())
  {
    std::print(stream, "{}", framework->printStatus(forceField));
  }
  for (std::size_t componentId{0}; const Component& component : components)
  {
    std::print(stream, "{}", component.printStatus(componentId, forceField, input_pressure));

    ++componentId;
  }
  std::print(stream, "\n\n\n\n");

  return stream.str();
}

nlohmann::json System::jsonComponentStatus() const
{
  nlohmann::json status;
  if (framework.has_value())
  {
    status[framework->name] = framework->jsonStatus();
  }
  for (const Component& component : components)
  {
    status[component.name] = component.jsonStatus();
  }

  return status;
}

void System::writeComponentFittingStatus(std::ostream& stream,
                                         const std::vector<std::pair<double, double>>& rawData) const
{
  std::print(stream, "Found {} data points\n", rawData.size());
  for (const std::pair<double, double>& data : rawData)
  {
    std::print(stream, "pressure: {:.8e}  loading: {}\n", data.first, data.second);
  }
  std::print(stream, "\n");

  if (!rawData.empty())
  {
    std::pair<double, double> pressureRange = std::make_pair(rawData.front().first, rawData.back().first);
    std::print(stream, "Lowest pressure:     {:.8e}\n", pressureRange.first);
    std::print(stream, "Highest pressure:    {:.8e}\n", pressureRange.second);
  }
  std::print(stream, "\n\n");
}
void System::writeCPUTimeStatistics(std::ostream& stream) const
{
  std::print(stream, "Sampling properties:        {:14f} [s]\n", mc_moves_cputime.propertySampling.count());
  std::print(stream, "Pressure computation:       {:14f} [s]\n\n", mc_moves_cputime.energyPressureComputation.count());

  for (std::size_t componentId = 0; const Component& component : components)
  {
    std::print(stream, "{}", component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(componentId, component.name));
  }
}
std::string System::writeMCMoveStatistics() const
{
  std::ostringstream stream;

  std::print(stream, "{}", mc_moves_statistics.writeMCMoveStatistics());
  for (std::size_t componentId = 0; const Component& component : components)
  {
    std::print(stream, "Component {} [{}]\n", componentId, component.name);

    std::print(stream, "{}", component.mc_moves_statistics.writeMCMoveStatistics());

    if (component.hasFractionalMolecule)
    {
      double imposedChemicalPotential = std::log(beta * component.fugacityCoefficient.value_or(1.0) * 
                                                 component.molFraction * pressure) / beta;
      double imposedFugacity = component.fugacityCoefficient.value_or(1.0) * component.molFraction * pressure;

      std::print(stream, "{}",
                 component.lambdaGC.writeAveragesStatistics(beta, imposedChemicalPotential, imposedFugacity));
      std::print(stream, "{}",
                 component.lambdaGC.writeDUdLambdaStatistics(beta, imposedChemicalPotential, imposedFugacity));
    }

    if (component.mc_moves_probabilities.getProbability(Move::Types::Widom) > 0.0)
    {
      double imposedChemicalPotential =
          std::log(beta * component.molFraction * component.fugacityCoefficient.value_or(1.0) * pressure) / beta;
      double imposedFugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * pressure;
      std::print(stream, "{}",
                 component.averageRosenbluthWeights.writeAveragesRosenbluthWeightStatistics(
                     temperature, simulationBox.volume, frameworkMass(),
                     framework.transform([](const Framework& f) { return f.numberOfUnitCells; })));
      std::print(stream, "{}",
                 component.averageRosenbluthWeights.writeAveragesChemicalPotentialStatistics(
                     beta, imposedChemicalPotential, imposedFugacity));
    }

    for (std::size_t i = 0; i != component.atoms.size(); ++i)
    {
      std::print(stream, "{}", component.cbmc_moves_statistics[i].writeMCMoveStatistics());
    }

    ++componentId;
  }

  if (usesReactionConventionalCFCMC())
  {
    for (const Reaction& reaction : reactions.list)
    {
      if (reaction.isSerialRxCFC())
      {
        std::print(stream, "reaction {} lambda statistics (reactant side, occupancy {:.6f}):\n", reaction.id,
                   reaction.lambda.occupancy());
        std::print(stream, "{}", reaction.lambda.writeAveragesStatistics(beta, std::nullopt, std::nullopt));
        std::print(stream, "{}", reaction.lambda.writeDUdLambdaStatistics(beta, std::nullopt, std::nullopt));

        std::print(stream, "reaction {} lambda statistics (product side, occupancy {:.6f}):\n", reaction.id,
                   reaction.lambdaProductSide.occupancy());
        std::print(stream, "{}",
                   reaction.lambdaProductSide.writeAveragesStatistics(beta, std::nullopt, std::nullopt));
        std::print(stream, "{}",
                   reaction.lambdaProductSide.writeDUdLambdaStatistics(beta, std::nullopt, std::nullopt));
        continue;
      }
      const PropertyLambdaProbabilityHistogram& histogram = activeReactionLambdaHistogram(reaction);
      std::print(stream, "reaction {} lambda statistics:\n", reaction.id);
      std::print(stream, "{}", histogram.writeAveragesStatistics(beta, std::nullopt, std::nullopt));
      std::print(stream, "{}", histogram.writeDUdLambdaStatistics(beta, std::nullopt, std::nullopt));
    }
  }

  std::print(stream, "\n\n");

  return stream.str();
}
nlohmann::json System::jsonMCMoveStatistics() const
{
  nlohmann::json status;

  status["system"] = mc_moves_statistics.jsonMCMoveStatistics();
  for (const Component& component : components)
  {
    status[component.name] = component.mc_moves_statistics.jsonMCMoveStatistics();

    if (component.hasFractionalMolecule)
    {
      double imposedChemicalPotential = std::log(beta * component.molFraction * pressure) / beta;
      double imposedFugacity = component.molFraction * pressure;

      status["lambdaStatistics"]["CFCMC"] =
          component.lambdaGC.jsonAveragesStatistics(beta, imposedChemicalPotential, imposedFugacity);
      status["lambdaStatistics"]["thermodynamicIntegration"] =
          component.lambdaGC.jsonDUdLambdaStatistics(beta, imposedChemicalPotential, imposedFugacity);
    }
  }

  return status;
}
