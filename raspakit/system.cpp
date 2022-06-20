module;

module system;

import randomnumbers;
import double3;
import cubic;
import atom;
import component;
import simulationbox;
import forcefield;
import double3x3;
import units;
import loadings;
import averages;
import skparser;
import skposcarparser;
import skstructure;
import skasymmetricatom;
import skcell;
import sample_movies;
import print;
import enthalpy_of_adsorption;
import energy_status;
import energy_status_inter;
import energy_status_intra;
import cbmc;
import property_simulationbox;
import property_energy;
import property_loading;
import property_enthalpy;
import lambda;
import property_widom;
import property_dudlambda;

import <numbers>;
import <complex>;
import <vector>;
import <random>;
import <span>;
import <tuple>;
import <iostream>;
import <fstream>;
import <streambuf>;
import <filesystem>;
import <optional>;
import <cmath>;

System::System(size_t s, ForceField forcefield, std::vector<Component> c, std::vector<size_t> initialNumberOfMolecules, size_t numberOfBlocks) :
               systemId(s),
               components(c),
               loadings(c.size()),
               averageLoadings(numberOfBlocks, c.size()),
               averageEnthalpiesOfAdsorption(numberOfBlocks, c.size()),
               swapableComponents(),
               numberOfMoleculesPerComponent(initialNumberOfMolecules),
               numberOfIntegerMoleculesPerComponent(c.size()),
               idealGasEnergiesPerComponent(c.size()),
               forceField(forcefield),
               numberOfPseudoAtoms(c.size(), std::vector<size_t>(forceField.pseudoAtoms.size())),
               totalNumberOfPseudoAtoms(forceField.pseudoAtoms.size()),
               simulationBox(0.0, 0.0, 0.0, 90.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0),
               averageSimulationBox(numberOfBlocks),
               atomPositions({}),
               atomVelocities({}),
               atomForces({}),
               runningEnergies(c.size()),
               averageEnergies(numberOfBlocks, c.size()),
               sampleMovie(systemId, forceField, simulationBox, atomPositions),
               netCharge(c.size())
{
    
}

System::System(System&& s) noexcept :
    systemId(s.systemId),
    components(s.components),
    loadings(s.loadings),
    averageLoadings(s.averageLoadings),
    averageEnthalpiesOfAdsorption(s.averageEnthalpiesOfAdsorption),
    swapableComponents(s.swapableComponents),
    numberOfMoleculesPerComponent(s.numberOfMoleculesPerComponent),
    numberOfIntegerMoleculesPerComponent(s.numberOfIntegerMoleculesPerComponent),
    numberOfFractionalMoleculesPerComponent(s.numberOfFractionalMoleculesPerComponent),
    numberOfReactionMoleculesPerComponent(s.numberOfReactionMoleculesPerComponent),
    numberOfReactionFractionalMoleculesPerComponent(s.numberOfReactionFractionalMoleculesPerComponent),
    idealGasEnergiesPerComponent(s.idealGasEnergiesPerComponent),
    forceField(s.forceField),
    numberOfPseudoAtoms(s.numberOfPseudoAtoms),
    totalNumberOfPseudoAtoms(s.totalNumberOfPseudoAtoms),
    simulationBox(std::move(s.simulationBox)),
    averageSimulationBox(std::move(s.averageSimulationBox)),
    atomPositions(std::move(s.atomPositions)),
    atomVelocities(std::move(s.atomVelocities)),
    atomForces(std::move(s.atomForces)),
    runningEnergies(s.runningEnergies),
    averageEnergies(s.averageEnergies),
    sampleMovie(std::move(s.sampleMovie)),
    netCharge(std::move(s.netCharge)),
    outputFile(std::move(s.outputFile))
{
}

void System::addComponent(const Component&& component) noexcept(false)
{
  numberOfMoleculesPerComponent.resize(components.size() + 1);
  numberOfIntegerMoleculesPerComponent.resize(components.size() + 1);
  numberOfFractionalMoleculesPerComponent.resize(components.size() + 1);
  numberOfReactionMoleculesPerComponent.resize(components.size() + 1);
  numberOfReactionFractionalMoleculesPerComponent.resize(components.size() + 1);
  idealGasEnergiesPerComponent.resize(components.size() + 1);

  numberOfPseudoAtoms.resize(components.size() + 1, std::vector<size_t>(forceField.pseudoAtoms.size()));
  totalNumberOfPseudoAtoms.resize(forceField.pseudoAtoms.size());

  runningEnergies.resize(components.size() + 1);
  averageEnergies.resize(components.size() + 1);

  loadings.resize(components.size() + 1);
  averageLoadings.resize(components.size() + 1);

  netCharge.resize(components.size() + 1);

  switch (component.type)
  {
    case Component::Type::Framework:
    {
      for(const Atom& atom: component.atoms)
      {
          atomPositions.push_back(atom);
      }
      numberOfFrameworkAtoms += component.atoms.size();
      numberOfMoleculesPerComponent[component.componentId] += 1;
      numberOfFractionalMoleculesPerComponent[component.componentId] = 0;
      numberOfIntegerMoleculesPerComponent[component.componentId] += 1;
      ++numberOfFrameworks;
      if(component.simulationBox.has_value())
      {
        // For multiple framework, the simulation box is the union of the boxes
        simulationBox = max(simulationBox, component.simulationBox.value());
      }
    }
    default:
      break;
  }

  // Move the component to the system
  components.push_back(std::move(component));
}

void System::insertFractionalMolecule(size_t selectedComponent, std::vector<Atom> atoms)
{
    std::vector<Atom>::const_iterator iterator = iteratorForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
    atomPositions.insert(iterator, atoms.begin(), atoms.end());
    numberOfMoleculesPerComponent[selectedComponent] += 1;
}

void System::insertMolecule(size_t selectedComponent, std::vector<Atom> atoms)
{
    std::vector<Atom>::const_iterator iterator = iteratorForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
    atomPositions.insert(iterator, atoms.begin(), atoms.end());
    numberOfMoleculesPerComponent[selectedComponent] += 1;
    numberOfIntegerMoleculesPerComponent[selectedComponent] += 1;

    // Update the number of pseudo atoms per type (used for tail-corrections)
    for(const Atom& atom: atoms)
    {
      numberOfPseudoAtoms[selectedComponent][static_cast<size_t>(atom.type)] += 1;
      totalNumberOfPseudoAtoms[static_cast<size_t>(atom.type)] += 1;
    }
   
    size_t index = 0; // indexOfFirstMolecule(selectedComponent);
    for (size_t componentId = 0; componentId < components.size(); componentId++)
    {
        for (size_t i = 0; i < numberOfMoleculesPerComponent[componentId]; ++i)
        {
            for (size_t j = 0; j < components[componentId].atoms.size(); ++j)
            {
                atomPositions[index].moleculeId = static_cast<short>(i);
                atomPositions[index].componentId = static_cast<short>(componentId);
                ++index;
            }
        }
    }
}

void System::deleteMolecule(size_t selectedComponent, size_t selectedMolecule, const std::span<Atom>& molecule)
{
    // Update the number of pseudo atoms per type (used for tail-corrections)
    for(const Atom& atom: molecule)
    {
      numberOfPseudoAtoms[selectedComponent][static_cast<size_t>(atom.type)] -= 1;
      totalNumberOfPseudoAtoms[static_cast<size_t>(atom.type)] -= 1;
    }

    std::vector<Atom>::const_iterator iterator = iteratorForMolecule(selectedComponent, selectedMolecule);
    atomPositions.erase(iterator, iterator + static_cast<std::vector<Atom>::difference_type>(molecule.size()));

    numberOfMoleculesPerComponent[selectedComponent] -= 1;
    numberOfIntegerMoleculesPerComponent[selectedComponent] -= 1;


    size_t index = 0; // indexOfFirstMolecule(selectedComponent);
    for (size_t componentId = 0; componentId < components.size(); componentId++)
    {
        for (size_t i = 0; i < numberOfMoleculesPerComponent[componentId]; ++i)
        {
            for (size_t j = 0; j < components[componentId].atoms.size(); ++j)
            {
                atomPositions[index].moleculeId = static_cast<short>(i);
                atomPositions[index].componentId = static_cast<short>(componentId);
                ++index;
            }
        }
    }
}

bool System::checkMoleculeIds()
{
    size_t index = 0; // indexOfFirstMolecule(selectedComponent);
    for (size_t componentId = 0; componentId < components.size(); componentId++)
    {
        for (size_t i = 0; i < numberOfMoleculesPerComponent[componentId]; ++i)
        {
            for (size_t j = 0; j < components[componentId].atoms.size(); ++j)
            {
                if (atomPositions[index].moleculeId != static_cast<int>(i)) return false;
                if (atomPositions[index].componentId != static_cast<int>(componentId)) return false;
                ++index;
            }
        }
    }
    return true;
}

void System::createInitialMolecules()
{
    // Make vectors same size and assume to number of atoms are zero if too short.
    numberOfMoleculesPerComponent.resize(components.size(), 0);

    for (size_t componentId = 0; const Component & component : components)
    {
        size_t moleculeIdx = 0;
        if (component.swapable)
        {
            numberOfMoleculesPerComponent[componentId] = 0;
            for (size_t i = 0; i < numberOfFractionalMoleculesPerComponent[componentId]; ++i)
            {
                std::optional<ChainData> growData = std::nullopt;
                do
                {
                    growData = growMoleculeSwapInsertion(componentId, numberOfMoleculesPerComponent[componentId], 0.0);

                } while (!growData || growData->RosenbluthWeight < 1.0);

                insertFractionalMolecule(componentId, growData->atom);
                moleculeIdx++;
            }
        }

        for (size_t i = 0; i < component.initialNumberOfMolecules; ++i)
        {
            std::optional<ChainData> growData = std::nullopt;
            do
            {
                growData = growMoleculeSwapInsertion(componentId, numberOfMoleculesPerComponent[componentId], 1.0);

            } while(!growData || growData->RosenbluthWeight < 1.0);

            insertMolecule(componentId, growData->atom);

            moleculeIdx++;
        }
        componentId++;
    }
}

size_t System::randomMoleculeOfComponent(size_t selectedComponent)
{
    return size_t(RandomNumber::Uniform() * static_cast<double>(numberOfMoleculesPerComponent[selectedComponent]));
}

size_t System::randomIntegerMoleculeOfComponent(size_t selectedComponent)
{
    return numberOfFractionalMoleculesPerComponent[selectedComponent] + 
           size_t(RandomNumber::Uniform() * static_cast<double>(numberOfIntegerMoleculesPerComponent[selectedComponent]));
}


std::vector<Atom>::const_iterator System::iteratorForMolecule(size_t selectedComponent, size_t selectedMolecule)
{
    size_t index{ 0 };
    for (size_t i = 0; i < selectedComponent; ++i)
    {
        size_t size = components[i].atoms.size();
        index += size * numberOfMoleculesPerComponent[i];
    }
    size_t size = components[selectedComponent].atoms.size();
    index += size * selectedMolecule;
    return atomPositions.cbegin() + static_cast<std::vector<Atom>::difference_type>(index);
}

std::span<const Atom> System::spanOfFrameworkAtoms() const
{
  return std::span(atomPositions.begin(), atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfFrameworkAtoms));
}

std::span<const Atom> System::spanOfMoleculeAtoms() const
{
  return std::span(atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfFrameworkAtoms), atomPositions.end());
}

std::span<Atom> System::spanOfMolecule(size_t selectedComponent, size_t selectedMolecule)
{
    size_t index{ 0 };
    for (size_t i = 0; i < selectedComponent; ++i)
    {
        size_t size = components[i].atoms.size();
        index += size * numberOfMoleculesPerComponent[i];
    }
    size_t size = components[selectedComponent].atoms.size();
    index += size * selectedMolecule;
    return std::span(&atomPositions[index], size);
}

size_t System::indexOfFirstMolecule(size_t selectedComponent)
{
    size_t index{ 0 };
    for (size_t i = 0; i < selectedComponent; ++i)
    {
        size_t size = components[i].atoms.size();
        index += size * numberOfMoleculesPerComponent[i];
    }
    return index;
}

void System::determineSwapableComponents()
{
    for (Component& component : components)
    {
        if (component.probabilitySwapMove_CBMC > 0.0)
        {
            component.swapable = true;
        }
        
        if (component.probabilitySwapMove_CFCMC > 0.0)
        {
            component.swapable = true;
        }

        if (component.probabilitySwapMove_CFCMC_CBMC > 0.0)
        {
            component.swapable = true;
        }

        if(component.swapable)
        {
            swapableComponents.push_back(component.componentId);
        }
    }
}

void System::determineFractionalComponents()
{
    for (size_t i = 0; i < components.size(); ++i)
    {
        if (components[i].probabilitySwapMove_CFCMC> 0.0)
        {
            numberOfFractionalMoleculesPerComponent[i] = 1;
            components[i].hasFractionalMolecule = true;
        }
        if (components[i].probabilitySwapMove_CFCMC_CBMC > 0.0)
        {
            numberOfFractionalMoleculesPerComponent[i] = 1;
            components[i].hasFractionalMolecule = true;
        }

        if (components[i].probabilityWidomMove_CFCMC > 0.0)
        {
            numberOfFractionalMoleculesPerComponent[i] = 1;
            components[i].hasFractionalMolecule = true;
        }
        if (components[i].probabilityWidomMove_CFCMC_CBMC > 0.0)
        {
            numberOfFractionalMoleculesPerComponent[i] = 1;
            components[i].hasFractionalMolecule = true;
        }
    }
}

void System::rescaleMoveProbabilities()
{
    for (Component& component : components)
    {
        component.normalizeMoveProbabilties();
    }
}

void System::removeRedundantMoves()
{
  for (Component& component : components)
  {
      // WidomMove_CFCMC already done when using SwapMove_CFCMC
      if(component.probabilityWidomMove_CFCMC > 0.0 && component.probabilitySwapMove_CFCMC > 0.0)
      {
          component.probabilityWidomMove_CFCMC = 0.0;
      }

      // WidomMove_CFCMC_CBMC already done when using SwapMove_CFCMC_CBMC
      if(component.probabilityWidomMove_CFCMC_CBMC > 0.0 && component.probabilitySwapMove_CFCMC_CBMC > 0.0)
      {
          component.probabilityWidomMove_CFCMC_CBMC = 0.0;
      }
  }
}

void System::rescaleMolarFractions()
{
    double totalMolfraction = 0.0;
    double numberOfSwapableComponents = 0.0;
    for (const Component& component : components)
    {
        if (component.swapable)
        {
            totalMolfraction += component.molFraction;
            numberOfSwapableComponents += 1.0;
        }
    }

    if (totalMolfraction > 0.0)
    {
        for (Component& component : components)
        {
            if (component.swapable)
            {
                component.molFraction /= totalMolfraction;
            }
        }
    }
    else
    {
        for (Component& component : components)
        {
            if (component.swapable)
            {
                component.molFraction /= numberOfSwapableComponents;
            }
        }
    }
}

void System::computeFrameworkDensity()
{
    for (Component& component : components)
    {
        if(component.type == Component::Type::Framework)
        {
            frameworkMass = frameworkMass.value_or(0.0) + component.mass;
        }
    }
}

void System::computeNumberOfPseudoAtoms()
{
  for(size_t i = 0; i != components.size(); ++i)
  {
    std::fill(numberOfPseudoAtoms[i].begin(), numberOfPseudoAtoms[i].end(), 0);
  }
  std::fill(totalNumberOfPseudoAtoms.begin(), totalNumberOfPseudoAtoms.end(), 0);

  for(const Atom& atom: atomPositions)
  {
    size_t componentId = static_cast<size_t>(atom.componentId);
    size_t type = static_cast<size_t>(atom.type);
    numberOfPseudoAtoms[componentId][type] += 1;
    totalNumberOfPseudoAtoms[type] += 1;
  } 
}


void System::createOutputFile()
{
    std::filesystem::path cwd = std::filesystem::current_path();

    std::filesystem::path directoryName = cwd / std::print("Output/System_{}/", systemId);
    std::filesystem::path fileName = cwd / std::print("Output/System_{}/output_{}_{}.data",
        systemId, simulationBox.temperature, simulationBox.pressure * Units::PressureConversionFactor);
    std::filesystem::create_directories(directoryName);

    outputFile = std::ofstream(fileName, std::ios::out);
}

void System::closeOutputFile()
{
     if (outputFile.is_open()) outputFile.close();
}

void System::writeToOutputFile(const std::string &s) 
{
    std::print(outputFile, s);
}

void System::writeOutputHeader()
{
    std::print(outputFile, "Compiler and run-time data\n");
    std::print(outputFile, "===============================================================================\n");

    std::print(outputFile, "RASPA 3.0.0\n\n");
}

void System::writeInitializationStatusReport(size_t currentCycle, size_t numberOfCycles)
{
    std::print(outputFile, "Initialization: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
    std::print(outputFile, "===============================================================================\n\n");

    std::print(outputFile, simulationBox.printStatus());
    std::print(outputFile, "\n");

    std::print(outputFile, "Amount of molecules per component :\n");
    std::print(outputFile, "-------------------------------------------------------------------------------\n");
    for (int compA = 0; Component & c : components)
    {
        std::print(outputFile, loadings.printStatus(c, frameworkMass));
        ++compA;
    }
    outputFile << "\n";
    double conv = Units::EnergyToKelvin;
    std::print(outputFile, "Total potential energy:   {: .6e}\n", conv * runningEnergies.totalEnergy);
    std::print(outputFile, "    Van der Waals:        {: .6e}\n", conv * runningEnergies.interEnergy.VanDerWaals);
    std::print(outputFile, "    Van der Waals (Tail): {: .6e}\n", conv * runningEnergies.interEnergy.VanDerWaalsTailCorrection);
    std::print(outputFile, "    Coulombic Real:       {: .6e}\n", conv * runningEnergies.interEnergy.CoulombicReal);
    std::print(outputFile, "    Coulombic Fourier:    {: .6e}\n", conv * runningEnergies.interEnergy.CoulombicFourier);
    std::print(outputFile, "    Intra:                {: .6e}\n", conv * runningEnergies.intraEnergy.total());
    std::print(outputFile, "    dU/dlambda:           {: .6e}\n", conv * runningEnergies.dUdlambda);

    outputFile << "\n\n\n\n";
}

void System::writeEquilibrationStatusReport(size_t currentCycle, size_t numberOfCycles)
{
    std::print(outputFile, "Equilibration: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
    std::print(outputFile, "===============================================================================\n\n");

    std::print(outputFile, simulationBox.printStatus());
    std::print(outputFile, "\n");

    std::print(outputFile, "Amount of molecules per component :\n");
    std::print(outputFile, "-------------------------------------------------------------------------------\n");
    for (int compA = 0; Component & c : components)
    {
        std::print(outputFile, loadings.printStatus(c, frameworkMass));
        ++compA;
    }
    outputFile << "\n";
    double conv = Units::EnergyToKelvin;
    std::print(outputFile, "Total potential energy:   {: .6e} [K]\n", conv * runningEnergies.totalEnergy);
    std::print(outputFile, "    Van der Waals:        {: .6e} [K]\n", conv * runningEnergies.interEnergy.VanDerWaals);
    std::print(outputFile, "    Van der Waals (Tail): {: .6e} [K]\n", conv * runningEnergies.interEnergy.VanDerWaalsTailCorrection);
    std::print(outputFile, "    Coulombic Real:       {: .6e} [K]\n", conv * runningEnergies.interEnergy.CoulombicReal);
    std::print(outputFile, "    Coulombic Fourier:    {: .6e} [K]\n", conv * runningEnergies.interEnergy.CoulombicFourier);
    std::print(outputFile, "    Intra:                {: .6e} [K]\n", conv * runningEnergies.intraEnergy.total());
    std::print(outputFile, "    dU/dlambda:           {: .6e}\n", conv * runningEnergies.dUdlambda);

    outputFile << "\n\n\n\n";
}

void System::writeProductionStatusReport(size_t currentCycle, size_t numberOfCycles)
{
    std::print(outputFile, "Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
    std::print(outputFile, "===============================================================================\n\n");

    std::pair<SimulationBox, SimulationBox> simulationBoxData = averageSimulationBox.averageSimulationBox();
    std::print(outputFile, simulationBox.printStatus(simulationBoxData.first, simulationBoxData.second));
    std::print(outputFile, "\n");
    
    std::print(outputFile, "Amount of molecules per component :\n");
    std::print(outputFile, "-------------------------------------------------------------------------------\n");
    std::pair<Loadings, Loadings> loadingData = averageLoadings.averageLoading();
    for (int compA = 0; Component & c : components)
    {
        std::print(outputFile, loadings.printStatus(c, loadingData.first, loadingData.second, frameworkMass));
        ++compA;
    }
    std::print(outputFile, "\n");
    double conv = Units::EnergyToKelvin;
    std::pair<EnergyStatus, EnergyStatus> energyData = averageEnergies.averageEnergy();
    std::print(outputFile, "Total potential energy :  {: .6e} ({: .6e} +/- {: .6e})\n",
        conv * runningEnergies.totalEnergy, conv * energyData.first.totalEnergy, conv * energyData.second.totalEnergy);
    std::print(outputFile, "    Van der Waals:        {: .6e} ({: .6e} +/- {: .6e})\n", 
        conv * runningEnergies.interEnergy.VanDerWaals,
        conv * energyData.first.interEnergy.VanDerWaals,
        conv * energyData.second.interEnergy.VanDerWaals);
    std::print(outputFile, "    Van der Waals (Tail): {: .6e} ({: .6e} +/- {: .6e})\n", 
        conv * runningEnergies.interEnergy.VanDerWaalsTailCorrection,
        conv * energyData.first.interEnergy.VanDerWaalsTailCorrection,
        conv * energyData.second.interEnergy.VanDerWaalsTailCorrection);
    std::print(outputFile, "    Coulombic Real:       {: .6e} ({: .6e} +/- {: .6e})\n",
        conv * runningEnergies.interEnergy.CoulombicReal,
        conv * energyData.first.interEnergy.CoulombicReal,
        conv * energyData.second.interEnergy.CoulombicReal);
    std::print(outputFile, "    Coulombic Fourier:    {: .6e} ({: .6e} +/- {: .6e})\n",
        conv * runningEnergies.interEnergy.CoulombicFourier, 
        conv * energyData.first.interEnergy.CoulombicFourier, 
        conv * energyData.second.interEnergy.CoulombicFourier);
    std::print(outputFile, "    Molecule Intra:       {: .6e} ({: .6e} +/- {: .6e})\n",
        conv * runningEnergies.intraEnergy.total(),
        conv * energyData.first.intraEnergy.total(),
        conv * energyData.second.intraEnergy.total());
    std::print(outputFile, "    dU/dlambda:           {: .6e}\n", conv * runningEnergies.dUdlambda);
    
    std::print(outputFile, "\n\n\n\n");
}

void System::writeComponentStatus()
{
    std::print(outputFile, "Component definitions\n");
    std::print(outputFile, "===============================================================================\n\n");
    for (const Component& component : components)
    {
        std::print(outputFile, component.printStatus(forceField));
    }
    std::print(outputFile, "\n\n\n\n");
}

void System::sampleProperties(size_t currentBlock)
{
   double w = weight();

   averageSimulationBox.addSample(currentBlock, simulationBox, w);

   averageEnergies.addSample(currentBlock, runningEnergies, w);

   loadings = Loadings(components.size(), numberOfIntegerMoleculesPerComponent, simulationBox);
   averageLoadings.addSample(currentBlock, loadings, w);

   EnthalpyOfAdsorptionTerms enthalpyTerms = EnthalpyOfAdsorptionTerms(swapableComponents, numberOfIntegerMoleculesPerComponent, runningEnergies.totalEnergy, simulationBox.temperature);
   averageEnthalpiesOfAdsorption.addSample(currentBlock, enthalpyTerms, w);

   for(Component &component : components)
   {
     double density  = static_cast<double>(numberOfIntegerMoleculesPerComponent[component.componentId]) / simulationBox.volume;

     component.lambda.sampleHistogram(currentBlock, density);
     component.lambda.sampleHistogram();
     component.averageRosenbluthWeights.addDensitySample(currentBlock, density, w);
     component.lambda.sampledUdLambdaHistogram(currentBlock, runningEnergies.dUdlambda);
     component.lambda.dUdlambdaBookKeeping.addDensitySample(currentBlock, density, w);
   }

}
