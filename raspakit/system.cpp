module;

module system;

import <cstddef>;
import <numbers>;
import <complex>;
import <vector>;
import <random>;
import <span>;
import <tuple>;
import <iostream>;
import <ostream>;
import <fstream>;
import <streambuf>;
import <filesystem>;
import <optional>;
import <cmath>;
import <chrono>;
import <algorithm>;
import <numeric>;
import <format>;
import <print>;
import <exception>;
import <source_location>;

import archive;
import randomnumbers;
import stringutils;
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
import enthalpy_of_adsorption;
import energy_status;
import energy_status_inter;
import energy_status_intra;
import property_simulationbox;
import property_energy;
import property_pressure;
import property_loading;
import property_enthalpy;
import property_lambda_probability_histogram;
import property_widom;
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
import mc_moves_cputime;
import reaction;
import reactions;
import cbmc;
import cbmc_chain_data;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;


// construct System programmatically
System::System(size_t id, double T, double P, ForceField forcefield, std::vector<Component> c, 
               std::vector<size_t> initialNumberOfMolecules, size_t numberOfBlocks) :
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
    numberOfGCFractionalMoleculesPerComponent_CFCMC(c.size()),
    numberOfPairGCFractionalMoleculesPerComponent_CFCMC(c.size()),
    numberOfGibbsFractionalMoleculesPerComponent_CFCMC(c.size()),
    numberOfReactionFractionalMoleculesPerComponent_CFCMC(),
    idealGasEnergiesPerComponent(c.size()),
    forceField(forcefield),
    hasExternalField(false),
    numberOfPseudoAtoms(c.size(), std::vector<size_t>(forceField.pseudoAtoms.size())),
    totalNumberOfPseudoAtoms(forceField.pseudoAtoms.size()),
    averageSimulationBox(numberOfBlocks),
    atomPositions({}),
    runningEnergies(),
    averageEnergies(numberOfBlocks, c.size()),
    currentEnergyStatus(c.size()),
    averagePressure(numberOfBlocks),
    netCharge(c.size()),
    mc_moves_probabilities(),
    mc_moves_statistics(),
    reactions(),
    tmmc(),
    computeConventionalRadialDistributionFunction(false),
    writeConventionalRadialDistributionFunctionEvery(1000),
    conventionalRadialDistributionFunctionHistogramSize(128),
    conventionalRadialDistributionFunctionRange(15.0),
    conventionalRadialDistributionFunction(numberOfBlocks, forceField.pseudoAtoms.size(), 
                                           conventionalRadialDistributionFunctionHistogramSize,
                                           conventionalRadialDistributionFunctionRange),
    computeRadialDistributionFunction(false),
    writeRadialDistributionFunctionEvery(1000),
    radialDistributionFunctionHistogramSize(128),
    radialDistributionFunctionRange(15.0),
    radialDistributionFunction(numberOfBlocks, forceField.pseudoAtoms.size(), radialDistributionFunctionHistogramSize,
                               radialDistributionFunctionRange)
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

  equationOfState = EquationOfState(EquationOfState::Type::PengRobinson, 
                                    EquationOfState::MultiComponentMixingRules::VanDerWaals, 
                                    T, P, simulationBox, HeliumVoidFraction, components);

  Interactions::computeEwaldFourierEnergySingleIon(eik_x, eik_y, eik_z, eik_xy,
                                                   forceField, simulationBox,
                                                   double3(0.0, 0.0, 0.0), 1.0);
  removeRedundantMoves();
  determineSwapableComponents();
  determineFractionalComponents();
  rescaleMoveProbabilities();
  rescaleMolarFractions();
  computeFrameworkDensity();
  computeNumberOfPseudoAtoms();

  double3 perpendicularWidths = simulationBox.perpendicularWidths();
  forceField.initializeEwaldParameters(perpendicularWidths);

  RandomNumber random(1400);
  createInitialMolecules(random);

  averageEnthalpiesOfAdsorption.resize(swapableComponents.size());
}


// used in 'input_reader.cpp' for Box and Framework
System::System(size_t s, ForceField forcefield, std::vector<Component> c, 
               [[maybe_unused]] std::vector<size_t> initialNumberOfMolecules, size_t numberOfBlocks) :
    systemId(s),
    components(c),
    loadings(c.size()),
    averageLoadings(numberOfBlocks, c.size()),
    averageEnthalpiesOfAdsorption(numberOfBlocks, c.size()),
    swapableComponents(),
    initialNumberOfMolecules(c.size()),
    numberOfMoleculesPerComponent(c.size()),
    numberOfIntegerMoleculesPerComponent(c.size()),
    numberOfFractionalMoleculesPerComponent(c.size()),
    numberOfGCFractionalMoleculesPerComponent_CFCMC(c.size()),
    numberOfPairGCFractionalMoleculesPerComponent_CFCMC(c.size()),
    numberOfGibbsFractionalMoleculesPerComponent_CFCMC(c.size()),
    numberOfReactionFractionalMoleculesPerComponent_CFCMC(),
    idealGasEnergiesPerComponent(c.size()),
    forceField(forcefield),
    hasExternalField(false),
    numberOfPseudoAtoms(c.size(), std::vector<size_t>(forceField.pseudoAtoms.size())),
    totalNumberOfPseudoAtoms(forceField.pseudoAtoms.size()),
    simulationBox(0.0, 0.0, 0.0, 90.0 * std::numbers::pi / 180.0, 
                  90.0 * std::numbers::pi / 180.0, 90.0 * std::numbers::pi / 180.0),
    averageSimulationBox(numberOfBlocks),
    atomPositions({}),
    runningEnergies(),
    averageEnergies(numberOfBlocks, c.size()),
    currentEnergyStatus(c.size()),
    averagePressure(numberOfBlocks),
    //sampleMovie(systemId, forceField, simulationBox, atomPositions),
    netCharge(c.size()),
    mc_moves_probabilities(),
    mc_moves_statistics(),
    reactions(),
    tmmc(),
    computeConventionalRadialDistributionFunction(false),
    writeConventionalRadialDistributionFunctionEvery(1000),
    conventionalRadialDistributionFunctionHistogramSize(128),
    conventionalRadialDistributionFunctionRange(15.0),
    conventionalRadialDistributionFunction(numberOfBlocks, forceField.pseudoAtoms.size(), 
                                           conventionalRadialDistributionFunctionHistogramSize,
                                           conventionalRadialDistributionFunctionRange),
    computeRadialDistributionFunction(false),
    writeRadialDistributionFunctionEvery(1000),
    radialDistributionFunctionHistogramSize(128),
    radialDistributionFunctionRange(15.0),
    radialDistributionFunction(numberOfBlocks, forceField.pseudoAtoms.size(), radialDistributionFunctionHistogramSize,
                               radialDistributionFunctionRange)
{
    
}

void System::addComponent(const Component&& component) noexcept(false)
{
  // Move the component to the system
  components.push_back(std::move(component));

  initialNumberOfMolecules.resize(components.size());
  numberOfMoleculesPerComponent.resize(components.size());
  numberOfIntegerMoleculesPerComponent.resize(components.size());
  numberOfFractionalMoleculesPerComponent.resize(components.size());

  numberOfGCFractionalMoleculesPerComponent_CFCMC.resize(components.size());
  numberOfPairGCFractionalMoleculesPerComponent_CFCMC.resize(components.size());
  numberOfGibbsFractionalMoleculesPerComponent_CFCMC.resize(components.size());
  
  // FIX: add reactions

  idealGasEnergiesPerComponent.resize(components.size());

  numberOfPseudoAtoms.resize(components.size(), std::vector<size_t>(forceField.pseudoAtoms.size()));
  totalNumberOfPseudoAtoms.resize(forceField.pseudoAtoms.size());

  averageEnergies.resize(components.size());

  loadings.resize(components.size());
  averageLoadings.resize(components.size());

  netCharge.resize(components.size());

}

// read the component from file if needed and initialize
void System::initializeComponents()
{
  // sort rigid frameworks first, before flexible frameworks
  std::sort(components.begin(), components.begin() + std::make_signed_t<std::size_t>(numberOfFrameworks), 
     [](const Component& a, const Component& b)
        {
            return a.rigid > b.rigid;
        });

  // read components from file
  for (Component& component : components)
  {
    switch (component.type)
    {
    case Component::Type::Adsorbate:
    case Component::Type::Cation:
      if (component.filenameData.has_value())
      {
       component.readComponent(forceField, component.filenameData.value());
      }
      break;
    case Component::Type::Framework:
      if (component.filenameData.has_value())
      {
        component.readFramework(forceField, component.filenameData.value());
      }

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

void System::insertFractionalMolecule(size_t selectedComponent, std::vector<Atom> atoms, size_t moleculeId)
{
  double l = 0.0;
  for (Atom& atom : atoms)
  {
    if(components[selectedComponent].lambdaGC.computeDUdlambda)
    {
      atom.moleculeId = static_cast<uint16_t>(moleculeId);
      atom.groupId = uint8_t{ 1 };
    }
    atom.setScaling(l);
  }
  std::vector<Atom>::const_iterator iterator = 
    iteratorForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  atomPositions.insert(iterator, atoms.begin(), atoms.end());
  numberOfMoleculesPerComponent[selectedComponent] += 1;

}

void System::insertMolecule(size_t selectedComponent, std::vector<Atom> atoms)
{
  std::vector<Atom>::const_iterator iterator = 
    iteratorForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  atomPositions.insert(iterator, atoms.begin(), atoms.end());
  numberOfMoleculesPerComponent[selectedComponent] += 1;
  numberOfIntegerMoleculesPerComponent[selectedComponent] += 1;

  // Update the number of pseudo atoms per type (used for tail-corrections)
  for(Atom& atom: atoms)
  {
    atom.moleculeId = static_cast<uint16_t>(numberOfMoleculesPerComponent[selectedComponent]);
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
        atomPositions[index].moleculeId = static_cast<uint16_t>(i);
        atomPositions[index].componentId = static_cast<uint8_t>(componentId);
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
        atomPositions[index].moleculeId = static_cast<uint16_t>(i);
        atomPositions[index].componentId = static_cast<uint8_t>(componentId);
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
        if (atomPositions[index].moleculeId != static_cast<uint32_t>(i)) return false;
        if (atomPositions[index].componentId != static_cast<uint8_t>(componentId)) return false;
        ++index;
      }
    }
  }
  return true;
}

void System::createInitialMolecules([[maybe_unused]] RandomNumber &random)
{
  for (size_t componentId = 0; const Component & component : components)
  {
    if (component.swapable)
    {
      numberOfMoleculesPerComponent[componentId] = 0;
      for (size_t i = 0; i < numberOfFractionalMoleculesPerComponent[componentId]; ++i)
      {
        std::optional<ChainData> growData = std::nullopt;
        do
        {
          std::vector<Atom> atoms = 
            components[componentId].recenteredCopy(0.0, numberOfMoleculesPerComponent[componentId]);
          Component::GrowType growType  = components[componentId].growType;
          growData = CBMC::growMoleculeSwapInsertion(random, this->hasExternalField, this->components, this->forceField, this->simulationBox, 
                                   this->spanOfFrameworkAtoms(), this->spanOfMoleculeAtoms(), this->beta,
                                   growType, forceField.cutOffVDW, forceField.cutOffCoulomb, componentId,
                                   numberOfMoleculesPerComponent[componentId], 0.0, atoms, numberOfTrialDirections);

        } while (!growData || growData->energies.total() > forceField.overlapCriteria);

        insertFractionalMolecule(componentId, growData->atom, i);
      }
    }

    for (size_t i = 0; i < initialNumberOfMolecules[componentId]; ++i)
    {
      std::optional<ChainData> growData = std::nullopt;
      do
      {
        std::vector<Atom> atoms = 
          components[componentId].recenteredCopy(1.0, numberOfMoleculesPerComponent[componentId]);
        Component::GrowType growType  = components[componentId].growType;
        growData = CBMC::growMoleculeSwapInsertion(random, this->hasExternalField, this->components, this->forceField, this->simulationBox, 
                                 this->spanOfFrameworkAtoms(), this->spanOfMoleculeAtoms(), this->beta, 
                                 growType, forceField.cutOffVDW, forceField.cutOffCoulomb, componentId,
                                 numberOfMoleculesPerComponent[componentId], 1.0, atoms, numberOfTrialDirections);

      } while(!growData || growData->energies.total() > forceField.overlapCriteria);

      insertMolecule(componentId, growData->atom);
    }
    componentId++;
  }
}

size_t System::randomMoleculeOfComponent(RandomNumber &random, size_t selectedComponent)
{
  return size_t(random.uniform() * static_cast<double>(numberOfMoleculesPerComponent[selectedComponent]));
}

size_t System::randomIntegerMoleculeOfComponent(RandomNumber &random, size_t selectedComponent)
{
  return numberOfFractionalMoleculesPerComponent[selectedComponent] + 
         size_t(random.uniform() * static_cast<double>(numberOfIntegerMoleculesPerComponent[selectedComponent]));
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
  return std::span(atomPositions.begin(), 
                   atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfFrameworkAtoms)); 
}

std::span<Atom> System::spanOfFrameworkAtoms()
{
  return std::span(atomPositions.begin(), 
                   atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfFrameworkAtoms));
}

std::span<const Atom> System::spanOfRigidFrameworkAtoms() const
{
  return std::span(atomPositions.begin(), 
                   atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfRigidFrameworkAtoms));
}


std::span<const Atom> System::spanOfFlexibleAtoms() const
{
  return std::span(atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfRigidFrameworkAtoms), 
                   atomPositions.end());
}


std::span<const Atom> System::spanOfMoleculeAtoms() const
{
  return std::span(atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfFrameworkAtoms), 
                   atomPositions.end());
}

std::span<Atom> System::spanOfMoleculeAtoms()
{
  return std::span(atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfFrameworkAtoms), 
                   atomPositions.end());
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
    if (component.mc_moves_probabilities.probabilitySwapMove_CBMC > 0.0 || 
        component.mc_moves_probabilities.probabilitySwapMove_CFCMC > 0.0 ||
        component.mc_moves_probabilities.probabilitySwapMove_CFCMC_CBMC > 0.0)
    {
      component.swapable = true;
    }

    if (component.mc_moves_probabilities.probabilityGibbsSwapMove_CBMC > 0.0 ||
        component.mc_moves_probabilities.probabilityGibbsSwapMove_CFCMC > 0.0 ||
        component.mc_moves_probabilities.probabilityGibbsSwapMove_CFCMC_CBMC > 0.0)
    {
      component.swapable = true;
    }
    
    if(component.swapable)
    {
      swapableComponents.push_back(component.componentId);
    }
  }
}


// determine the required number of fractional molecules 
void System::determineFractionalComponents()
{
  for (size_t i = 0; i < components.size(); ++i)
  {
    numberOfFractionalMoleculesPerComponent[i] = 0;
    numberOfGCFractionalMoleculesPerComponent_CFCMC[i] = 0;
    numberOfGibbsFractionalMoleculesPerComponent_CFCMC[i] = 0;
   

    if (components[i].mc_moves_probabilities.probabilitySwapMove_CFCMC> 0.0 || 
        components[i].mc_moves_probabilities.probabilityWidomMove_CFCMC > 0.0 ||
        components[i].mc_moves_probabilities.probabilitySwapMove_CFCMC_CBMC > 0.0 || 
        components[i].mc_moves_probabilities.probabilityWidomMove_CFCMC_CBMC > 0.0)
    {
      numberOfFractionalMoleculesPerComponent[i] += 1;
      numberOfGCFractionalMoleculesPerComponent_CFCMC[i] = 1;
      components[i].hasFractionalMolecule = true;
    }

    // Gibbs
    if (components[i].mc_moves_probabilities.probabilityGibbsSwapMove_CFCMC > 0.0 || 
        components[i].mc_moves_probabilities.probabilityGibbsSwapMove_CFCMC_CBMC > 0.0)
    {
      numberOfFractionalMoleculesPerComponent[i] += 1;
      numberOfGibbsFractionalMoleculesPerComponent_CFCMC[i] = 1;
      components[i].hasFractionalMolecule = true;
    }
  }

  for (size_t reactionId{ 0 }; [[maybe_unused]] const Reaction& reaction : reactions.list)
  {
    for (size_t j = 0; j < components.size(); ++j)
    {
      numberOfReactionFractionalMoleculesPerComponent_CFCMC[reactionId][j] = 0;
    }

    //for ([[maybe_unused]]  size_t componentId{ 0 };  const size_t [[maybe_unused]]  stoichiometry : reaction.reactantStoichiometry)
    //{
    //  //numberOfReactionFractionalMoleculesPerComponent[reactionId][componentId] += stoichiometry;
    //  ++componentId;
    //}
    //for ([[maybe_unused]] size_t componentId{ 0 };  const [[maybe_unused]] size_t stoichiometry : reaction.productStoichiometry)
    //{
    //  //numberOfReactionFractionalMoleculesPerComponent[reactionId][componentId] += stoichiometry;
    //  ++componentId;
    //}
    ++reactionId;
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
    if(component.mc_moves_probabilities.probabilityWidomMove_CFCMC > 0.0 && 
       component.mc_moves_probabilities.probabilitySwapMove_CFCMC > 0.0)
    {
      component.mc_moves_probabilities.probabilityWidomMove_CFCMC = 0.0;
    }

    // WidomMove_CFCMC_CBMC already done when using SwapMove_CFCMC_CBMC
    if(component.mc_moves_probabilities.probabilityWidomMove_CFCMC_CBMC > 0.0 && 
       component.mc_moves_probabilities.probabilitySwapMove_CFCMC_CBMC > 0.0)
    {
      component.mc_moves_probabilities.probabilityWidomMove_CFCMC_CBMC = 0.0;
    }
  }
}

void System::optimizeMCMoves()
{
  mc_moves_statistics.optimizeAcceptance();
  for (Component& component : components)
  {
    component.mc_moves_statistics.optimizeMCMoves();
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

std::vector<Atom> System::randomConfiguration(RandomNumber &random, size_t selectedComponent, 
                                              const std::span<const Atom> molecule)
{
  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());
  double3 position = simulationBox.randomPosition(random);
  size_t startingBead = components[selectedComponent].startingBead;
  for (size_t i = 0; i != molecule.size(); ++i)
  {
    copied_atoms[i].position = position + 
                               randomRotationMatrix * (molecule[i].position - molecule[startingBead].position);
  }
  return copied_atoms;
}

std::string System::writeOutputHeader() const
{
  std::ostringstream stream;

  std::print(stream, "Compiler and run-time data\n");
  std::print(stream, "===============================================================================\n");

  std::print(stream, "RASPA 3.0.0\n\n");

  ThreadPool &pool = ThreadPool::instance();
  const size_t numberOfHelperThreads = pool.getThreadCount();

  switch(pool.threadingType)
  {
    case ThreadPool::ThreadingType::Serial:
      std::print(stream, "Parallization: Serial, 1 thread\n");
      break;
    case ThreadPool::ThreadingType::OpenMP:
      std::print(stream, "Parallization: OpenMP, {} threads\n", numberOfHelperThreads + 1);
      break;
    case ThreadPool::ThreadingType::ThreadPool:
      std::print(stream, "Parallization: ThreadPool, {} threads\n", numberOfHelperThreads + 1);
      break;
    case ThreadPool::ThreadingType::GPU_Offload:
      std::print(stream, "Parallization: GPU-Offload\n");
      break;
  } 
  std::print(stream, "\n");

  return stream.str();
}

std::string System::writeInitializationStatusReport(size_t currentCycle, size_t numberOfCycles) const
{
  std::ostringstream stream;

  std::print(stream, "Initialization: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::print(stream, "{}", simulationBox.printStatus());
  std::print(stream, "\n");

  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (const Component & c : components)
  {
    std::print(stream, "{}", loadings.printStatus(c, frameworkMass));
  }
  std::print(stream, "\n");
  double conv = Units::EnergyToKelvin;

  //size_t numberOfMolecules = std::reduce(numberOfIntegerMoleculesPerComponent.begin(),  numberOfIntegerMoleculesPerComponent.end());
  //double currentExcessPressure = runningEnergies.totalEnergy.forceFactor / (3.0 * simulationBox.volume);
  //double currentIdealPressure =  static_cast<double>(numberOfMolecules)/(simulationBox.Beta * simulationBox.volume);
  //double currentPressure = currentIdealPressure + currentExcessPressure;
  //std::print(outputFile, "Pressure:             {: .6e} [bar]\n", 1e-5 * Units::PressureConversionFactor * currentPressure);

  for (const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();
    
    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", 
                         c.componentId, c.name, lambda, runningEnergies.dudlambda(lambda), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", 
                 c.componentId, c.name, c.lambdaGC.lambdaValue(), occupancy, averageOccupancy);
    }
  }
  std::print(stream, "\n");

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

  return stream.str();
}

std::string System::writeEquilibrationStatusReport(size_t currentCycle, size_t numberOfCycles) const
{
  std::ostringstream stream;

  std::print(stream, "Equilibration: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::print(stream, "{}", simulationBox.printStatus());
  std::print(stream, "\n");

  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (const Component & c : components)
  {
    std::print(stream, "{}", loadings.printStatus(c, frameworkMass));
  }
  std::print(stream, "\n");
  double conv = Units::EnergyToKelvin;

  for (const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    { 
      std::print(stream, "component {} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", 
                 c.componentId, c.name, lambda, runningEnergies.dudlambda(lambda), occupancy, averageOccupancy);
    }
    else 
    {
      std::print(stream, "component {} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", 
                 c.componentId, c.name, c.lambdaGC.lambdaValue(), occupancy, averageOccupancy);
    }
  }
  std::print(stream, "\n");

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

  return stream.str();
}

std::string System::writeProductionStatusReport(size_t currentCycle, size_t numberOfCycles) const
{
  std::ostringstream stream;

  std::print(stream, "Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::pair<SimulationBox, SimulationBox> simulationBoxData = averageSimulationBox.averageSimulationBox();
  std::print(stream, "{}", simulationBox.printStatus(simulationBoxData.first, simulationBoxData.second));
  std::print(stream, "\n");
  
  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::pair<Loadings, Loadings> loadingData = averageLoadings.averageLoading();
  for (const Component & c : components)
  {
    std::print(stream, "{}", loadings.printStatus(c, loadingData.first, loadingData.second, frameworkMass));
  }
  std::print(stream, "\n");
  double conv = Units::EnergyToKelvin;

  std::pair<double3x3, double3x3> currentPressureTensor = averagePressure.averagePressureTensor();
  double3x3 pressureTensor = 1e-5 * Units::PressureConversionFactor * currentPressureTensor.first;
  double3x3 pressureTensorError = 1e-5 * Units::PressureConversionFactor * currentPressureTensor.second;
  std::print(stream, "Average pressure tensor: \n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", 
          pressureTensor.ax, pressureTensor.bx, pressureTensor.cx, 
          pressureTensorError.ax, pressureTensorError.bx, pressureTensorError.cx);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", 
          pressureTensor.ay, pressureTensor.by, pressureTensor.cy, 
          pressureTensorError.ay, pressureTensorError.by, pressureTensorError.cy);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", 
          pressureTensor.az, pressureTensor.bz, pressureTensor.cz, 
          pressureTensorError.az, pressureTensorError.bz, pressureTensorError.cz);
  std::pair<double, double> idealGasPressure = averagePressure.averageIdealGasPressure();
  std::pair<double, double> excessPressure = averagePressure.averageExcessPressure();
  std::pair<double, double> p = averagePressure.averagePressure();
  std::print(stream, "Ideal-gas pressure:  {: .6e} +/ {:.6e} [bar]\n", 
          1e-5 * Units::PressureConversionFactor * idealGasPressure.first, 
          1e-5 * Units::PressureConversionFactor * idealGasPressure.second);
  std::print(stream, "Excess pressure:     {: .6e} +/ {:.6e} [bar]\n", 
          1e-5 * Units::PressureConversionFactor * excessPressure.first, 
          1e-5 * Units::PressureConversionFactor * excessPressure.second);
  std::print(stream, "Pressure:            {: .6e} +/ {:.6e} [bar]\n\n", 
          1e-5 * Units::PressureConversionFactor * p.first, 
          1e-5 * Units::PressureConversionFactor * p.second);

  for (const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", 
                 c.componentId, c.name, lambda, runningEnergies.dudlambda(lambda), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", 
                 c.componentId, c.name, c.lambdaGC.lambdaValue(), occupancy, averageOccupancy);
    }
  }
  std::print(stream, "\n");

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

  return stream.str();
}

std::string System::writeComponentStatus() const
{
  std::ostringstream stream;

  std::print(stream, "Component definitions\n");
  std::print(stream, "===============================================================================\n\n");
  for (const Component& component : components)
  {
    std::print(stream, "{}", component.printStatus(forceField));
  }
  std::print(stream, "\n\n\n\n");

  return stream.str();
}

void System::writeComponentFittingStatus(std::ostream &stream, 
                                         const std::vector<std::pair<double, double>> &rawData) const
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

void System::sampleProperties(size_t currentBlock, size_t currentCycle)
{
  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  double w = weight();

  averageSimulationBox.addSample(currentBlock, simulationBox, w);

  loadings = Loadings(components.size(), numberOfIntegerMoleculesPerComponent, simulationBox);
  averageLoadings.addSample(currentBlock, loadings, w);

  EnthalpyOfAdsorptionTerms enthalpyTerms = 
    EnthalpyOfAdsorptionTerms(swapableComponents, numberOfIntegerMoleculesPerComponent, 
                              runningEnergies.total(), temperature);
  averageEnthalpiesOfAdsorption.addSample(currentBlock, enthalpyTerms, w);

  size_t numberOfMolecules = 
    std::reduce(numberOfIntegerMoleculesPerComponent.begin(),  numberOfIntegerMoleculesPerComponent.end());
  double currentIdealPressure =  static_cast<double>(numberOfMolecules)/(beta * simulationBox.volume);

  averagePressure.addSample(currentBlock, currentIdealPressure, currentExcessPressureTensor, w);

  for(Component &component : components)
  {
    double componentDensity  = 
      static_cast<double>(numberOfIntegerMoleculesPerComponent[component.componentId]) / simulationBox.volume;

    double lambda = component.lambdaGC.lambdaValue();
    double dudlambda = runningEnergies.dudlambda(lambda);
    component.lambdaGC.sampleHistogram(currentBlock, componentDensity, dudlambda, containsTheFractionalMolecule, w);

    component.averageRosenbluthWeights.addDensitySample(currentBlock, componentDensity, w);
  }

  if (currentCycle % 10uz == 0uz)
  {
    if(computeConventionalRadialDistributionFunction)
    {
      conventionalRadialDistributionFunction.sample(simulationBox, spanOfFrameworkAtoms(),
                                                    spanOfMoleculeAtoms(), currentBlock);
    }

    if(computeRadialDistributionFunction)
    {
      radialDistributionFunction.sample(simulationBox, spanOfFrameworkAtoms(),
                                        spanOfMoleculeAtoms(), currentBlock);
    }
  }


  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

  mc_moves_cputime.propertySampling += (t2 - t1);
}

void System::writeCPUTimeStatistics(std::ostream &stream) const
{
  std::print(stream, "Sampling properties:        {:14f} [s]\n", mc_moves_cputime.propertySampling.count());
  std::print(stream, "Pressure computation:       {:14f} [s]\n\n", mc_moves_cputime.energyPressureComputation.count());


 for (size_t componentId = 0; const Component& component : components)
 {
   std::print(stream, "{}", component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(componentId, component.name));
 }
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
  mc_moves_statistics.clear();
}

inline std::pair<EnergyStatus, double3x3> pair_acc(const std::pair<EnergyStatus, double3x3> &lhs, 
                                                   const std::pair<EnergyStatus, double3x3> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

void System::precomputeTotalRigidEnergy() noexcept
{
  rigidEnergies.zero();
  Interactions::computeEwaldFourierRigidEnergy(eik_x, eik_y, eik_z, eik_xy, 
                                               fixedFrameworkStoredEik, forceField, simulationBox,
                                               spanOfRigidFrameworkAtoms(), rigidEnergies);
}

void System::recomputeTotalEnergies() noexcept
{
  runningEnergies.zero();

  if(fixedFrameworkStoredEik.empty())
  {
    precomputeTotalRigidEnergy();
  }

  std::span<const Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<const Atom> moleculeAtomPositions = spanOfMoleculeAtoms();
  [[maybe_unused]]std::span<const Atom> flexibleAtomPositions = spanOfFlexibleAtoms();

  Interactions::computeFrameworkMoleculeEnergy(forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions, runningEnergies);
  Interactions::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtomPositions, runningEnergies);

  Interactions::computeFrameworkMoleculeTailEnergy(forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions, runningEnergies);
  Interactions::computeInterMolecularTailEnergy(forceField, simulationBox,moleculeAtomPositions, runningEnergies);

  Interactions::computeEwaldFourierEnergy(eik_x, eik_y, eik_z, eik_xy,
                                          fixedFrameworkStoredEik, storedEik,
                                          forceField, simulationBox,
                                          components, numberOfMoleculesPerComponent,
                                          flexibleAtomPositions, atomPositions, runningEnergies);

  // correct for the energy of rigid parts
  runningEnergies -= rigidEnergies;
}


RunningEnergy System::computeTotalEnergies() noexcept
{
  RunningEnergy runningEnergy{};

  if(fixedFrameworkStoredEik.empty())
  {
    precomputeTotalRigidEnergy();
  }

  std::span<const Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<const Atom> moleculeAtomPositions = spanOfMoleculeAtoms();
  [[maybe_unused]]std::span<const Atom> flexibleAtomPositions = spanOfFlexibleAtoms();

  Interactions::computeFrameworkMoleculeEnergy(forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions, runningEnergy);
  Interactions::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtomPositions, runningEnergy);

  Interactions::computeFrameworkMoleculeTailEnergy(forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions, runningEnergy);
  Interactions::computeInterMolecularTailEnergy(forceField, simulationBox, moleculeAtomPositions, runningEnergy);

  Interactions::computeEwaldFourierEnergy(eik_x, eik_y, eik_z, eik_xy, 
                                          fixedFrameworkStoredEik, storedEik,
                                          forceField, simulationBox,
                                          components, numberOfMoleculesPerComponent,
                                          flexibleAtomPositions, atomPositions, runningEnergy);

  // correct for the energy of rigid parts
  runningEnergy -= rigidEnergies;

  return runningEnergy;
}

void System::computeTotalGradients() noexcept
{
}

std::pair<EnergyStatus, double3x3> System::computeMolecularPressure() noexcept
{
  for(Atom& atom: atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeFrameworkMoleculeEnergyStrainDerivative(forceField, components,
                                                      simulationBox, spanOfFrameworkAtoms(), spanOfMoleculeAtoms());


  pressureInfo = pair_acc(pressureInfo, Interactions::computeInterMolecularEnergyStrainDerivative(forceField, components, simulationBox,
                                                      spanOfMoleculeAtoms()));
  pressureInfo = pair_acc(pressureInfo, Interactions::computeEwaldFourierEnergyStrainDerivative(eik_x, eik_y, eik_z, eik_xy,
                                                                    fixedFrameworkStoredEik, storedEik, forceField, simulationBox,
                                                                    components, numberOfMoleculesPerComponent, atomPositions));

  pressureInfo.first.sumTotal();

  // Correct rigid molecule contribution using the constraints forces
  double3x3 correctionTerm;
  for(size_t componentId = 0; componentId < components.size(); ++componentId)
  {
    if(components[componentId].rigid)
    {
      for(size_t i = 0; i < numberOfMoleculesPerComponent[componentId]; ++i)
      {
        std::span<Atom> span = spanOfMolecule(componentId, i);

        double totalMass = 0.0;
        double3 com(0.0, 0.0, 0.0);
        for(const Atom& atom: span)
        {
          double mass = forceField.pseudoAtoms[static_cast<size_t>(atom.type)].mass;
          com += mass * atom.position;
          totalMass += mass;
        }
        com = com / totalMass;

        for(const Atom& atom: span)
        {
          correctionTerm.ax += (atom.position.x - com.x) * atom.gradient.x;
          correctionTerm.ay += (atom.position.x - com.x) * atom.gradient.y;
          correctionTerm.az += (atom.position.x - com.x) * atom.gradient.z;

          correctionTerm.bx += (atom.position.y - com.y) * atom.gradient.x;
          correctionTerm.by += (atom.position.y - com.y) * atom.gradient.y;
          correctionTerm.bz += (atom.position.y - com.y) * atom.gradient.z;

          correctionTerm.cx += (atom.position.z - com.z) * atom.gradient.x;
          correctionTerm.cy += (atom.position.z - com.z) * atom.gradient.y;
          correctionTerm.cz += (atom.position.z - com.z) * atom.gradient.z;
        }
      }
    }
  }
  pressureInfo.second = -(pressureInfo.second - correctionTerm);

  return pressureInfo;
}

void System::MD_Loop()
{
}

std::vector<Atom> 
System::equilibratedMoleculeRandomInBox(RandomNumber &random, size_t selectedComponent, std::span<Atom> molecule, 
                                        double scaling, size_t moleculeId) const
{
  size_t startingBead = components[selectedComponent].startingBead;
  double3 center = molecule[startingBead].position;
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());

  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  double3 position = simulationBox.randomPosition(random);

  for (size_t i = 0; i != copied_atoms.size(); ++i)
  {
    copied_atoms[i].setScaling(scaling);
    copied_atoms[i].position = position + randomRotationMatrix * (molecule[i].position - center);
    copied_atoms[i].moleculeId = static_cast<uint32_t>(moleculeId);
  }
  return copied_atoms;
}

inline std::string formatMoveStatistics(const std::string name, const MoveStatistics<double>& move)
{
  std::ostringstream stream;

  std::print(stream, "{} total:        {:10}\n", name, move.counts);
  std::print(stream, "{} constructed:  {:10}\n", name, move.constructed);
  std::print(stream, "{} accepted:     {:10}\n", name, move.accepted);
  std::print(stream, "{} fraction:     {:10f}\n", name, move.accepted / std::max(1.0, double(move.counts)));
  std::print(stream, "{} max-change:   {:10f}\n\n", name, move.maxChange);

  return stream.str();
}

inline std::string formatMoveStatistics(const std::string name, const MoveStatistics<double3> &move)
{
  std::ostringstream stream;

  std::print(stream, "{} total:        {:10} {:10} {:10}\n", 
                     name, move.counts.x, move.counts.y, move.counts.z);
  std::print(stream, "{} constructed:  {:10} {:10} {:10}\n", 
                     name, move.constructed.x, move.constructed.y, move.constructed.z);
  std::print(stream, "{} accepted:     {:10} {:10} {:10}\n", 
                     name, move.accepted.x, move.accepted.y, move.accepted.z);
  std::print(stream, "{} fraction:     {:10f} {:10f} {:10f}\n", 
                     name, 
                     move.accepted.x / std::max(1.0, double(move.counts.x)),
                     move.accepted.y / std::max(1.0, double(move.counts.y)), 
                     move.accepted.z / std::max(1.0, double(move.counts.z)));
  std::print(stream, "{} max-change:   {:10f} {:10f} {:10f}\n\n", 
                     name, move.maxChange.x, move.maxChange.y, move.maxChange.z);

  return stream.str();
}

std::string System::writeMCMoveStatistics() const
{
  std::ostringstream stream;

  if (mc_moves_statistics.volumeMove.totalCounts > 0.0) 
    std::print(stream, "{}", formatMoveStatistics( "Volume", mc_moves_statistics.volumeMove));
  if (mc_moves_statistics.GibbsVolumeMove.totalCounts > 0.0) 
    std::print(stream, "{}", formatMoveStatistics("Gibbs Volume", mc_moves_statistics.GibbsVolumeMove));

  for (size_t componentId = 0; const Component& component: components)
  {
    std::print(stream,"Component {} [{}]\n", componentId, component.name);

    std::print(stream, "{}", component.mc_moves_statistics.writeMCMoveStatistics());

    if(component.hasFractionalMolecule)
    {
      double imposedChemicalPotential = std::log(beta * component.molFraction * pressure) / beta;
      double imposedFugacity = component.molFraction * pressure;

      std::print(stream, "{}", 
                 component.lambdaGC.writeAveragesStatistics(beta, imposedChemicalPotential, imposedFugacity));
      std::print(stream, "{}", 
                 component.lambdaGC.writeDUdLambdaStatistics(beta, imposedChemicalPotential, imposedFugacity));
    }

    if(component.mc_moves_probabilities.probabilityWidomMove > 0.0)
    {
      double imposedChemicalPotential = std::log(beta * component.molFraction * pressure) / beta;
      double imposedFugacity = component.molFraction * pressure;
      std::print(stream, "{}", 
        component.averageRosenbluthWeights.writeAveragesStatistics(beta, imposedChemicalPotential, imposedFugacity));
    }

    ++componentId;
  }

 
  std::print(stream, "\n\n");

  return stream.str();
}



Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const System &s)
{
  archive << s.versionNumber;

  archive << s.systemId;
  archive << s.temperature;
  archive << s.pressure;
  archive << s.input_pressure;
  archive << s.beta;
  archive << s.HeliumVoidFraction;
  archive << s.numberOfFrameworks;
  archive << s.numberOfFrameworkAtoms;
  archive << s.numberOfRigidFrameworkAtoms;
  archive << s.components;
  archive << s.loadings;
  archive << s.averageLoadings;
  archive << s.averageEnthalpiesOfAdsorption;
  archive << s.swapableComponents;
  archive << s.initialNumberOfMolecules;
  archive << s.numberOfMoleculesPerComponent;
  archive << s.numberOfIntegerMoleculesPerComponent;
  archive << s.numberOfFractionalMoleculesPerComponent;
  archive << s.numberOfGCFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfPairGCFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfGibbsFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfReactionFractionalMoleculesPerComponent_CFCMC;
  archive << s.idealGasEnergiesPerComponent;
  archive << s.forceField;
  archive << s.hasExternalField;
  archive << s.numberOfPseudoAtoms;
  archive << s.totalNumberOfPseudoAtoms;
  archive << s.frameworkMass;
  archive << s.timeStep;
  archive << s.simulationBox;
  archive << s.averageSimulationBox;
  archive << s.atomPositions;
  archive << s.runningEnergies;
  archive << s.rigidEnergies;
  archive << s.averageEnergies;
  archive << s.currentExcessPressureTensor;
  archive << s.currentEnergyStatus;
  archive << s.averagePressure;
  archive << s.numberOfTrialDirections;
  archive << s.eik_xy;
  archive << s.eik_x;
  archive << s.eik_y;
  archive << s.eik_z;
  archive << s.storedEik;
  archive << s.fixedFrameworkStoredEik;
  archive << s.totalEik;
  archive << s.CoulombicFourierEnergySingleIon;
  archive << s.netCharge;
  archive << s.mc_moves_probabilities;
  archive << s.mc_moves_statistics;
  archive << s.mc_moves_cputime;
  archive << s.mc_moves_count;
  archive << s.reactions;
  archive << s.tmmc;
  archive << s.columnNumberOfGridPoints;
  archive << s.columnTotalPressure;
  archive << s.columnPressureGradient;
  archive << s.columnVoidFraction;
  archive << s.columnParticleDensity;
  archive << s.columnEntranceVelocity;
  archive << s.columnLength;
  archive << s.columnTimeStep;
  archive << s.columnNumberOfTimeSteps;
  archive << s.columnAutoNumberOfTimeSteps;
  archive << s.mixturePredictionMethod;
  archive << s.pressure_range;
  archive << s.numberOfCarrierGases;
  archive << s.carrierGasComponent;
  archive << s.maxIsothermTerms;
  archive << s.containsTheFractionalMolecule;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, System &s)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > s.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'System' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> s.systemId;
  archive >> s.temperature;
  archive >> s.pressure;
  archive >> s.input_pressure;
  archive >> s.beta;
  archive >> s.HeliumVoidFraction;
  archive >> s.numberOfFrameworks;
  archive >> s.numberOfFrameworkAtoms;
  archive >> s.numberOfRigidFrameworkAtoms;
  archive >> s.components;
  archive >> s.loadings;
  archive >> s.averageLoadings;
  archive >> s.averageEnthalpiesOfAdsorption;
  archive >> s.swapableComponents;
  archive >> s.initialNumberOfMolecules;
  archive >> s.numberOfMoleculesPerComponent;
  archive >> s.numberOfIntegerMoleculesPerComponent;
  archive >> s.numberOfFractionalMoleculesPerComponent;
  archive >> s.numberOfGCFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfPairGCFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfGibbsFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfReactionFractionalMoleculesPerComponent_CFCMC;
  archive >> s.idealGasEnergiesPerComponent;
  archive >> s.forceField;
  archive >> s.hasExternalField;
  archive >> s.numberOfPseudoAtoms;
  archive >> s.totalNumberOfPseudoAtoms;
  archive >> s.frameworkMass;
  archive >> s.timeStep;
  archive >> s.simulationBox;
  archive >> s.averageSimulationBox;
  archive >> s.atomPositions;
  archive >> s.runningEnergies;
  archive >> s.rigidEnergies;
  archive >> s.averageEnergies;
  archive >> s.currentExcessPressureTensor;
  archive >> s.currentEnergyStatus;
  archive >> s.averagePressure;
  archive >> s.numberOfTrialDirections;
  archive >> s.eik_xy;
  archive >> s.eik_x;
  archive >> s.eik_y;
  archive >> s.eik_z;
  archive >> s.storedEik;
  archive >> s.fixedFrameworkStoredEik;
  archive >> s.totalEik;
  archive >> s.CoulombicFourierEnergySingleIon;
  archive >> s.netCharge;
  archive >> s.mc_moves_probabilities;
  archive >> s.mc_moves_statistics;
  archive >> s.mc_moves_cputime;
  archive >> s.mc_moves_count;
  archive >> s.reactions;
  archive >> s.tmmc;
  archive >> s.columnNumberOfGridPoints;
  archive >> s.columnTotalPressure;
  archive >> s.columnPressureGradient;
  archive >> s.columnVoidFraction;
  archive >> s.columnParticleDensity;
  archive >> s.columnEntranceVelocity;
  archive >> s.columnLength;
  archive >> s.columnTimeStep;
  archive >> s.columnNumberOfTimeSteps;
  archive >> s.columnAutoNumberOfTimeSteps;
  archive >> s.mixturePredictionMethod;
  archive >> s.pressure_range;
  archive >> s.numberOfCarrierGases;
  archive >> s.carrierGasComponent;
  archive >> s.maxIsothermTerms;
  archive >> s.containsTheFractionalMolecule;

  return archive;
}

