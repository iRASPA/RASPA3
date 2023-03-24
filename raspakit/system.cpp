module;

module system;

import randomnumbers;
import int3;
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
import property_pressure;
import property_loading;
import property_enthalpy;
import lambda;
import property_widom;
import property_dudlambda;
import energy_factor;
import running_energy;
import threadpool;
import isotherm;
import multi_site_isotherm;
import pressure_range;
import bond_potential;
import move_statistics;
import mc_moves_probabilities_system;
import mc_moves_probabilities_particles;

import <cstddef>;
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
import <chrono>;

System::System(size_t id, double T, double P, ForceField forcefield, std::vector<Component> c, std::vector<size_t> initialNumberOfMolecules, size_t numberOfBlocks) :
    systemId(id), 
    temperature(T),
    pressure(P / Units::PressureConversionFactor),
    input_pressure(P),
    beta(1.0 / (Units::KB * T)),
    components(c),
    loadings(c.size()),
    averageLoadings(numberOfBlocks, c.size()),
    averageEnthalpiesOfAdsorption(numberOfBlocks, c.size()),
    swapableComponents(),
    initialNumberOfMolecules(initialNumberOfMolecules),
    numberOfMoleculesPerComponent(c.size()),
    numberOfIntegerMoleculesPerComponent(c.size()),
    numberOfFractionalMoleculesPerComponent(c.size()),
    numberOfReactionMoleculesPerComponent(c.size()),
    numberOfReactionFractionalMoleculesPerComponent(c.size()),
    idealGasEnergiesPerComponent(c.size()),
    forceField(forcefield),
    numberOfPseudoAtoms(c.size(), std::vector<size_t>(forceField.pseudoAtoms.size())),
    totalNumberOfPseudoAtoms(forceField.pseudoAtoms.size()),
    averageSimulationBox(numberOfBlocks),
    atomPositions({}),
    runningEnergies(),
    averageEnergies(numberOfBlocks, c.size()),
    currentEnergyStatus(c.size()),
    averagePressure(numberOfBlocks),
    netCharge(c.size()),
    lambda(numberOfBlocks, 41)
{

  for (Component& component : components)
  {
     if(component.type == Component::Type::Framework)
     {  
        numberOfMoleculesPerComponent[component.componentId] += 1;
        numberOfFractionalMoleculesPerComponent[component.componentId] = 0;
        numberOfIntegerMoleculesPerComponent[component.componentId] += 1;
        ++numberOfFrameworks;

        // add Framework-atoms to the atomPositions-list
        for (const Atom& atom : component.atoms)
        {
            atomPositions.push_back(atom);
        }
        numberOfFrameworkAtoms += component.atoms.size();

        if (component.rigid)
        {
            numberOfRigidFrameworkAtoms += component.atoms.size();
        }
     }

     // For multiple framework, the simulation box is the union of the boxes
     simulationBox = max(simulationBox, component.simulationBox.scaled(component.numberOfUnitCells));
  }

  registerEwaldFourierEnergySingleIon(double3(0.0, 0.0, 0.0), 0.0);
  removeRedundantMoves();
  determineSwapableComponents();
  determineFractionalComponents();
  rescaleMoveProbabilities();
  rescaleMolarFractions();
  computeComponentFluidProperties();
  computeFrameworkDensity();
  computeNumberOfPseudoAtoms();

  double3 perpendicularWidths = simulationBox.perpendicularWidths();
  forceField.initializeEwaldParameters(perpendicularWidths);

  createInitialMolecules();

  averageEnthalpiesOfAdsorption.resize(swapableComponents.size());
}


System::System(size_t s, ForceField forcefield, std::vector<Component> c, [[maybe_unused]] std::vector<size_t> initialNumberOfMolecules, size_t numberOfBlocks) :
               systemId(s),
               components(c),
               loadings(c.size()),
               averageLoadings(numberOfBlocks, c.size()),
               averageEnthalpiesOfAdsorption(numberOfBlocks, c.size()),
               swapableComponents(),
               initialNumberOfMolecules(c.size()),
               numberOfMoleculesPerComponent(c.size()),
               numberOfIntegerMoleculesPerComponent(c.size()),
               idealGasEnergiesPerComponent(c.size()),
               forceField(forcefield),
               numberOfPseudoAtoms(c.size(), std::vector<size_t>(forceField.pseudoAtoms.size())),
               totalNumberOfPseudoAtoms(forceField.pseudoAtoms.size()),
               simulationBox(0.0, 0.0, 0.0, 90.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0),
               averageSimulationBox(numberOfBlocks),
               atomPositions({}),
               runningEnergies(),
               averageEnergies(numberOfBlocks, c.size()),
               currentEnergyStatus(c.size()),
               averagePressure(numberOfBlocks),
               //sampleMovie(systemId, forceField, simulationBox, atomPositions),
               netCharge(c.size()),
               mc_moves_probabilities(),
               lambda(numberOfBlocks, 41)
{
    
}

System::System(System&& s) noexcept :
    systemId(s.systemId),
    temperature(s.temperature),
    pressure(s.pressure),
    input_pressure(s.input_pressure),
    beta(s.beta),
    numberOfFrameworks(s.numberOfFrameworks),
    components(std::move(s.components)),
    loadings(s.loadings),
    averageLoadings(s.averageLoadings),
    averageEnthalpiesOfAdsorption(s.averageEnthalpiesOfAdsorption),
    swapableComponents(s.swapableComponents),
    initialNumberOfMolecules(s.initialNumberOfMolecules),
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
    runningEnergies(s.runningEnergies),
    averageEnergies(s.averageEnergies),
    currentEnergyStatus(s.currentEnergyStatus),
    averagePressure(s.averagePressure),
    //sampleMovie(std::move(s.sampleMovie)),
    netCharge(std::move(s.netCharge)),
    mc_moves_probabilities(s.mc_moves_probabilities),
    lambda(s.lambda),
    columnNumberOfGridPoints(s.columnNumberOfGridPoints),
    columnTotalPressure(s.columnTotalPressure),
    columnPressureGradient(s.columnPressureGradient),
    columnVoidFraction(s.columnVoidFraction),
    columnParticleDensity(s.columnParticleDensity),
    columnEntranceVelocity(s.columnEntranceVelocity),
    columnLength(s.columnLength),
    columnTimeStep(s.columnTimeStep),
    columnNumberOfTimeSteps(s.columnNumberOfTimeSteps),
    columnAutoNumberOfTimeSteps(s.columnAutoNumberOfTimeSteps),
    mixturePredictionMethod(s.mixturePredictionMethod),
    pressure_range(s.pressure_range)
{
}

System::System(const System&& s) noexcept :
    systemId(s.systemId),
    temperature(s.temperature),
    pressure(s.pressure),
    input_pressure(s.input_pressure),
    beta(s.beta),
    numberOfFrameworks(s.numberOfFrameworks),
    components(std::move(s.components)),
    loadings(s.loadings),
    averageLoadings(s.averageLoadings),
    averageEnthalpiesOfAdsorption(s.averageEnthalpiesOfAdsorption),
    swapableComponents(s.swapableComponents),
    initialNumberOfMolecules(s.initialNumberOfMolecules),
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
    //atomVelocities(std::move(s.atomVelocities)),
    //atomForces(std::move(s.atomForces)),
    runningEnergies(s.runningEnergies),
    averageEnergies(s.averageEnergies),
    currentEnergyStatus(s.currentEnergyStatus),
    averagePressure(s.averagePressure),
    //sampleMovie(std::move(s.sampleMovie)),
    netCharge(std::move(s.netCharge)),
    mc_moves_probabilities(s.mc_moves_probabilities),
    lambda(s.lambda),
    columnNumberOfGridPoints(s.columnNumberOfGridPoints),
    columnTotalPressure(s.columnTotalPressure),
    columnPressureGradient(s.columnPressureGradient),
    columnVoidFraction(s.columnVoidFraction),
    columnParticleDensity(s.columnParticleDensity),
    columnEntranceVelocity(s.columnEntranceVelocity),
    columnLength(s.columnLength),
    columnTimeStep(s.columnTimeStep),
    columnNumberOfTimeSteps(s.columnNumberOfTimeSteps),
    columnAutoNumberOfTimeSteps(s.columnAutoNumberOfTimeSteps),
    mixturePredictionMethod(s.mixturePredictionMethod),
    pressure_range(s.pressure_range)
{
}

void System::addComponent(const Component&& component) noexcept(false)
{
  initialNumberOfMolecules.resize(components.size() + 1);
  numberOfMoleculesPerComponent.resize(components.size() + 1);
  numberOfIntegerMoleculesPerComponent.resize(components.size() + 1);
  numberOfFractionalMoleculesPerComponent.resize(components.size() + 1);
  numberOfReactionMoleculesPerComponent.resize(components.size() + 1);
  numberOfReactionFractionalMoleculesPerComponent.resize(components.size() + 1);
  idealGasEnergiesPerComponent.resize(components.size() + 1);

  numberOfPseudoAtoms.resize(components.size() + 1, std::vector<size_t>(forceField.pseudoAtoms.size()));
  totalNumberOfPseudoAtoms.resize(forceField.pseudoAtoms.size());

  averageEnergies.resize(components.size() + 1);

  loadings.resize(components.size() + 1);
  averageLoadings.resize(components.size() + 1);

  netCharge.resize(components.size() + 1);

  // Move the component to the system
  components.push_back(std::move(component));
}

// read the component from file if needed and initialize
void System::initializeComponents()
{
  // sort components, order: rigid frameworks, flexible frameworks, adsorbates, cations
  std::sort(components.begin(), components.end(), [](const Component& a, const Component& b)
        {
          if (a.type == b.type)
          {
            if (a.rigid) return true; else return false;
          }
          else
          {
            return a.type < b.type;
          }
        });

  // read components from file
  for (Component& component : components)
  {
      if (component.filenameData.has_value())
      {
         switch (component.type)
         {
         case Component::Type::Adsorbate:
         case Component::Type::Cation:
           component.readComponent(forceField, component.filenameData.value());
           break;
         case Component::Type::Framework:
           component.readFramework(forceField, component.filenameData.value());

           numberOfMoleculesPerComponent[component.componentId] += 1;
           numberOfFractionalMoleculesPerComponent[component.componentId] = 0;
           numberOfIntegerMoleculesPerComponent[component.componentId] += 1;
           ++numberOfFrameworks;

           // add Framework-atoms to the atomPositions-list
           for (const Atom& atom : component.atoms)
           {
             atomPositions.push_back(atom);
           }
           numberOfFrameworkAtoms += component.atoms.size();

           if (component.rigid)
           {
               numberOfRigidFrameworkAtoms += component.atoms.size();
           }
           
           // For multiple framework, the simulation box is the union of the boxes
           simulationBox = max(simulationBox, component.simulationBox.scaled(component.numberOfUnitCells));
           break;
         }
      }
  }
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
    //numberOfMoleculesPerComponent.resize(components.size(), 0);

    for (size_t componentId = 0; const Component & component : components)
    {
        //size_t moleculeIdx = 0;
        if (component.swapable)
        {
            numberOfMoleculesPerComponent[componentId] = 0;
            for (size_t i = 0; i < numberOfFractionalMoleculesPerComponent[componentId]; ++i)
            {
                std::optional<ChainData> growData = std::nullopt;
                do
                {
                    growData = growMoleculeSwapInsertion(forceField.cutOffVDW, forceField.cutOffCoulomb, 
                                        componentId, numberOfMoleculesPerComponent[componentId], 0.0);

                } while (!growData || growData->RosenbluthWeight < 1.0);

                insertFractionalMolecule(componentId, growData->atom);
                //moleculeIdx++;
            }
        }

        for (size_t i = 0; i < initialNumberOfMolecules[componentId]; ++i)
        {
            std::optional<ChainData> growData = std::nullopt;
            do
            {
                growData = growMoleculeSwapInsertion(forceField.cutOffVDW, forceField.cutOffCoulomb,
                                    componentId, numberOfMoleculesPerComponent[componentId], 1.0);

            } while(!growData || growData->RosenbluthWeight < 1.0);

            insertMolecule(componentId, growData->atom);

            //moleculeIdx++;
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

std::span<const Atom> System::spanOfRigidFrameworkAtoms() const
{
    return std::span(atomPositions.begin(), atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfRigidFrameworkAtoms));
}


std::span<const Atom> System::spanOfFlexibleAtoms() const
{
    return std::span(atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfRigidFrameworkAtoms), atomPositions.end());
}


std::span<const Atom> System::spanOfMoleculeAtoms() const
{
  return std::span(atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfFrameworkAtoms), atomPositions.end());
}

std::span<Atom> System::spanOfMoleculeAtoms()
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

const std::span<const Atom> System::spanOfMolecule(size_t selectedComponent, size_t selectedMolecule) const
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
        if (component.mc_moves_probabilities.probabilitySwapMove_CBMC > 0.0)
        {
            component.swapable = true;
        }
        
        if (component.mc_moves_probabilities.probabilitySwapMove_CFCMC > 0.0)
        {
            component.swapable = true;
        }

        if (component.mc_moves_probabilities.probabilitySwapMove_CFCMC_CBMC > 0.0)
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
        if (components[i].mc_moves_probabilities.probabilitySwapMove_CFCMC> 0.0)
        {
            numberOfFractionalMoleculesPerComponent[i] = 1;
            components[i].hasFractionalMolecule = true;
        }
        if (components[i].mc_moves_probabilities.probabilitySwapMove_CFCMC_CBMC > 0.0)
        {
            numberOfFractionalMoleculesPerComponent[i] = 1;
            components[i].hasFractionalMolecule = true;
        }

        if (components[i].mc_moves_probabilities.probabilityWidomMove_CFCMC > 0.0)
        {
            numberOfFractionalMoleculesPerComponent[i] = 1;
            components[i].hasFractionalMolecule = true;
        }
        if (components[i].mc_moves_probabilities.probabilityWidomMove_CFCMC_CBMC > 0.0)
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
    component.mc_moves_probabilities.probabilityVolumeMove = mc_moves_probabilities.probabilityVolumeMove;
    component.mc_moves_probabilities.probabilityGibbsVolumeMove = mc_moves_probabilities.probabilityGibbsVolumeMove;
  
    component.mc_moves_probabilities.normalizeMoveProbabilties();
  }
}

void System::removeRedundantMoves()
{
  for (Component& component : components)
  {
      // WidomMove_CFCMC already done when using SwapMove_CFCMC
      if(component.mc_moves_probabilities.probabilityWidomMove_CFCMC > 0.0 && component.mc_moves_probabilities.probabilitySwapMove_CFCMC > 0.0)
      {
          component.mc_moves_probabilities.probabilityWidomMove_CFCMC = 0.0;
      }

      // WidomMove_CFCMC_CBMC already done when using SwapMove_CFCMC_CBMC
      if(component.mc_moves_probabilities.probabilityWidomMove_CFCMC_CBMC > 0.0 && component.mc_moves_probabilities.probabilitySwapMove_CFCMC_CBMC > 0.0)
      {
          component.mc_moves_probabilities.probabilityWidomMove_CFCMC_CBMC = 0.0;
      }
  }
}

void System::optimizeMCMoves()
{
  mc_moves_probabilities.optimizeAcceptance();
  for (Component& component : components)
  {
    component.mc_moves_probabilities.optimizeMCMoves();
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

void System::writeOutputHeader(std::ostream &outputFile) const
{
  std::print(outputFile, "Compiler and run-time data\n");
  std::print(outputFile, "===============================================================================\n");

  std::print(outputFile, "RASPA 3.0.0\n\n");

  ThreadPool &pool = ThreadPool::instance();
  const size_t numberOfHelperThreads = pool.getThreadCount();

  switch(pool.threadingType)
  {
    case ThreadPool::ThreadingType::Serial:
      std::print(outputFile, "Parallization: Serial, 1 thread\n");
      break;
    case ThreadPool::ThreadingType::OpenMP:
      std::print(outputFile, "Parallization: OpenMP, {} threads\n", numberOfHelperThreads + 1);
      break;
    case ThreadPool::ThreadingType::ThreadPool:
      std::print(outputFile, "Parallization: ThreadPool, {} threads\n", numberOfHelperThreads + 1);
      break;
    case ThreadPool::ThreadingType::GPU_Offload:
      std::print(outputFile, "Parallization: GPU-Offload\n");
      break;
  } 
  std::print(outputFile, "\n");
}

void System::writeInitializationStatusReport(std::ostream &stream, [[maybe_unused]] size_t currentCycle, [[maybe_unused]] size_t numberOfCycles) const
{
  std::print(stream, "Initialization: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  simulationBox.printStatus(stream);
  std::print(stream, "\n");

  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (const Component & c : components)
  {
    loadings.printStatus(stream, c, frameworkMass);
  }
  std::print(stream, "\n");
  double conv = Units::EnergyToKelvin;

  //size_t numberOfMolecules = std::reduce(numberOfIntegerMoleculesPerComponent.begin(),  numberOfIntegerMoleculesPerComponent.end());
  //double currentExcessPressure = runningEnergies.totalEnergy.forceFactor / (3.0 * simulationBox.volume);
  //double currentIdealPressure =  static_cast<double>(numberOfMolecules)/(simulationBox.Beta * simulationBox.volume);
  //double currentPressure = currentIdealPressure + currentExcessPressure;
  //std::print(outputFile, "Pressure:             {: .6e} [bar]\n", 1e-5 * Units::PressureConversionFactor * currentPressure);
  std::print(stream, "dU/dlambda:           {: .6e} [K]\n\n", conv * runningEnergies.dUdlambda);

  std::print(stream, "Total potential energy:      {: .6e} [K]\n", conv * runningEnergies.total());
  std::print(stream, "    framework-molecule VDW:  {: .6e} [K]\n", conv * runningEnergies.frameworkMoleculeVDW);
  std::print(stream, "    framework-molecule Real: {: .6e} [K]\n", conv * runningEnergies.frameworkMoleculeCharge);
  std::print(stream, "    molecule-molecule VDW:   {: .6e} [K]\n", conv * runningEnergies.moleculeMoleculeVDW);
  std::print(stream, "    molecule-molecule Real:  {: .6e} [K]\n", conv * runningEnergies.moleculeMoleculeCharge);
  std::print(stream, "    Van der Waals (Tail):    {: .6e} [K]\n", conv * runningEnergies.tail);
  std::print(stream, "    Coulombic Fourier:       {: .6e} [K]\n", conv * runningEnergies.ewald);
  std::print(stream, "    Intra VDW:               {: .6e} [K]\n", conv * runningEnergies.intraVDW);
  std::print(stream, "    Intra Charge:            {: .6e} [K]\n", conv * runningEnergies.intraCoul);
  std::print(stream, "    Polarization:            {: .6e} [K]\n", conv * runningEnergies.polarization);

  std::print(stream, "\n\n\n\n");
}

void System::writeEquilibrationStatusReport(std::ostream &stream, [[maybe_unused]] size_t currentCycle, [[maybe_unused]] size_t numberOfCycles) const
{
  std::print(stream, "Equilibration: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  simulationBox.printStatus(stream);
  std::print(stream, "\n");

  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (const Component & c : components)
  {
    loadings.printStatus(stream, c, frameworkMass);
  }
  std::print(stream, "\n");
  double conv = Units::EnergyToKelvin;

  std::print(stream, "dU/dlambda:           {: .6e} [K]\n\n", conv * runningEnergies.dUdlambda);

  std::print(stream, "Total potential energy:      {: .6e} [K]\n", conv * runningEnergies.total());
  std::print(stream, "    framework-molecule VDW:  {: .6e} [K]\n", conv * runningEnergies.frameworkMoleculeVDW);
  std::print(stream, "    framework-molecule Real: {: .6e} [K]\n", conv * runningEnergies.frameworkMoleculeCharge);
  std::print(stream, "    molecule-molecule VDW:   {: .6e} [K]\n", conv * runningEnergies.moleculeMoleculeVDW);
  std::print(stream, "    molecule-molecule Real:  {: .6e} [K]\n", conv * runningEnergies.moleculeMoleculeCharge);
  std::print(stream, "    Van der Waals (Tail):    {: .6e} [K]\n", conv * runningEnergies.tail);
  std::print(stream, "    Coulombic Fourier:       {: .6e} [K]\n", conv * runningEnergies.ewald);
  std::print(stream, "    Intra VDW:               {: .6e} [K]\n", conv * runningEnergies.intraVDW);
  std::print(stream, "    Intra Charge:            {: .6e} [K]\n", conv * runningEnergies.intraCoul);
  std::print(stream, "    Polarization:            {: .6e} [K]\n", conv * runningEnergies.polarization);

  std::print(stream, "\n\n\n\n");
}

void System::writeProductionStatusReport(std::ostream &stream, [[maybe_unused]] size_t currentCycle, [[maybe_unused]] size_t numberOfCycles) const
{
  std::print(stream, "Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::pair<SimulationBox, SimulationBox> simulationBoxData = averageSimulationBox.averageSimulationBox();
  std::print(stream, simulationBox.printStatus(simulationBoxData.first, simulationBoxData.second));
  std::print(stream, "\n");
  
  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::pair<Loadings, Loadings> loadingData = averageLoadings.averageLoading();
  for (const Component & c : components)
  {
    loadings.printStatus(stream, c, loadingData.first, loadingData.second, frameworkMass);
  }
  std::print(stream, "\n");
  double conv = Units::EnergyToKelvin;

  std::pair<double3x3, double3x3> currentPressureTensor = averagePressure.averagePressureTensor();
  double3x3 pressureTensor = 1e-5 * Units::PressureConversionFactor * currentPressureTensor.first;
  double3x3 pressureTensorError = 1e-5 * Units::PressureConversionFactor * currentPressureTensor.second;
  std::print(stream, "Average pressure tensor: \n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", 
          pressureTensor.ax, pressureTensor.bx, pressureTensor.cx, pressureTensorError.ax, pressureTensorError.bx, pressureTensorError.cx);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", 
          pressureTensor.ay, pressureTensor.by, pressureTensor.cy, pressureTensorError.ay, pressureTensorError.by, pressureTensorError.cy);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", 
          pressureTensor.az, pressureTensor.bz, pressureTensor.cz, pressureTensorError.az, pressureTensorError.bz, pressureTensorError.cz);
  std::pair<double, double> idealGasPressure = averagePressure.averageIdealGasPressure();
  std::pair<double, double> excessPressure = averagePressure.averageExcessPressure();
  std::pair<double, double> p = averagePressure.averagePressure();
  std::print(stream, "Ideal-gas pressure:  {: .6e} +/ {:.6e} [bar]\n", 
          1e-5 * Units::PressureConversionFactor * idealGasPressure.first, 1e-5 * Units::PressureConversionFactor * idealGasPressure.second);
  std::print(stream, "Excess pressure:     {: .6e} +/ {:.6e} [bar]\n", 
          1e-5 * Units::PressureConversionFactor * excessPressure.first, 1e-5 * Units::PressureConversionFactor * excessPressure.second);
  std::print(stream, "Pressure:            {: .6e} +/ {:.6e} [bar]\n\n", 
          1e-5 * Units::PressureConversionFactor * p.first, 1e-5 * Units::PressureConversionFactor * p.second);

  std::print(stream, "dU/dlambda:           {: .6e} [K]\n\n", conv * runningEnergies.dUdlambda);

  std::pair<EnergyStatus, EnergyStatus> energyData = averageEnergies.averageEnergy();
  std::print(stream, "Total potential energy :  {: .6e} ({: .6e} +/- {:.6e}) [K]\n",
      conv * currentEnergyStatus.totalEnergy.energy, 
      conv * energyData.first.totalEnergy.energy, 
      conv * energyData.second.totalEnergy.energy);
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
      conv * currentEnergyStatus.intraEnergy.total().energy,
      conv * energyData.first.intraEnergy.total().energy,
      conv * energyData.second.intraEnergy.total().energy);
  
  std::print(stream, "\n\n\n");
}

void System::writeComponentStatus(std::ostream &stream) const
{
  std::print(stream, "Component definitions\n");
  std::print(stream, "===============================================================================\n\n");
  for (const Component& component : components)
  {
    component.printStatus(stream, forceField);
  }
  std::print(stream, "\n\n\n\n");
}

void System::writeComponentFittingStatus(std::ostream &stream, [[maybe_unused]] const std::vector<std::pair<double, double>> &rawData) const
{
  std::print(stream, "Found {} data points\n", rawData.size());
  for(const std::pair<double, double> &data : rawData)
  {
    std::print(stream, "pressure: {:.8e}  loading: {}\n", data.first, data.second);
  }
  std::print(stream, "\n");

  if(!rawData.empty())
  {
    std::pair<double, double> pressureRange = std::make_pair(rawData.front().first, rawData.back().first);
    std::print(stream,"Lowest pressure:     {:.8e}\n", pressureRange.first);
    std::print(stream,"Highest pressure:    {:.8e}\n", pressureRange.second);
  }
  std::print(stream, "\n\n");
}

void System::sampleProperties(size_t currentBlock)
{
   std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
   double w = weight();

   averageSimulationBox.addSample(currentBlock, simulationBox, w);

   loadings = Loadings(components.size(), numberOfIntegerMoleculesPerComponent, simulationBox);
   averageLoadings.addSample(currentBlock, loadings, w);

   EnthalpyOfAdsorptionTerms enthalpyTerms = EnthalpyOfAdsorptionTerms(swapableComponents, numberOfIntegerMoleculesPerComponent, runningEnergies.total(), temperature);
   averageEnthalpiesOfAdsorption.addSample(currentBlock, enthalpyTerms, w);

   size_t numberOfMolecules = std::reduce(numberOfIntegerMoleculesPerComponent.begin(),  numberOfIntegerMoleculesPerComponent.end());
   double currentIdealPressure =  static_cast<double>(numberOfMolecules)/(beta * simulationBox.volume);

   averagePressure.addSample(currentBlock, currentIdealPressure, currentExcessPressureTensor, w);

   for(Component &component : components)
   {
     double density  = static_cast<double>(numberOfIntegerMoleculesPerComponent[component.componentId]) / simulationBox.volume;

     component.lambda.sampleHistogram(currentBlock, density);
     component.averageRosenbluthWeights.addDensitySample(currentBlock, density, w);
   }

   lambda.sampledUdLambdaHistogram(currentBlock, runningEnergies.dUdlambda);
   //lambda.dUdlambdaBookKeeping.addDensitySample(currentBlock, density, w);

   std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

  cpuTime_Sampling += (t2 - t1);
}

void System::writeCPUTimeStatistics(std::ostream &stream) const
{
  std::print(stream, "Sampling properties:        {:14f} [s]\n", cpuTime_Sampling.count());
  std::print(stream, "Pressure computation:       {:14f} [s]\n\n", cpuTime_Pressure.count());
  for (int componentId = 0; const Component& component : components)
  {
    std::print(stream, "Component {} [{}]\n", componentId, component.name);
    std::print(stream, component.mc_moves_probabilities.writeMCMoveCPUTimeStatistics());
    componentId++;
  }

  mc_moves_probabilities.writeMCMoveCPUTimeStatistics();

 
}

std::vector<Atom> System::scaledCenterOfMassPositions(double scale) const
{
  std::vector<Atom> scaledAtoms;
  scaledAtoms.reserve(atomPositions.size());

  for(size_t componentId = 0; componentId < components.size(); ++componentId)
  {
    for(size_t i = 0; i < numberOfMoleculesPerComponent[componentId]; ++i)
    {
      const std::span<const Atom> span = spanOfMolecule(componentId, i);

      double totalMass = 0.0;
      double3 com(0.0, 0.0, 0.0);
      for(const Atom& atom: span)
      {
        double mass = forceField.pseudoAtoms[static_cast<size_t>(atom.type)].mass;
        com += mass * atom.position;
        totalMass += mass;
      }
      com = com / totalMass;

      double3 d = com * (scale - 1.0);

      // create copy
      for(Atom atom: span)
      {
        atom.position += d; 
        scaledAtoms.push_back(atom);
      }
    }
  }
  return scaledAtoms;
}

void System::clearMoveStatistics()
{
  mc_moves_probabilities.clear();
}

void System::clearTimingStatistics()
{
  mc_moves_probabilities.clearTimingStatistics();
}

