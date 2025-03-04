module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <mdspan>
#include <numbers>
#include <numeric>
#include <optional>
#include <ostream>
#include <print>
#include <random>
#include <source_location>
#include <span>
#include <streambuf>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>
#endif

module system;

#ifndef USE_LEGACY_HEADERS
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
import <exception>;
import <source_location>;
import <map>;
import <vector>;
import <array>;
import <string>;
import <string_view>;
import <print>;
import <mdspan>;
#endif

import archive;
import randomnumbers;
import stringutils;
import int3;
import double3;
import double3x3;
import double3x3x3;
import simd_quatd;
import cubic;
import atom;
import framework;
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
import skatom;
import skcell;
import sample_movies;
import enthalpy_of_adsorption;
import energy_factor;
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
import property_temperature;
import property_msd;
import energy_factor;
import running_energy;
import threadpool;
import isotherm;
import multi_site_isotherm;
import pressure_range;
import bond_potential;
import move_statistics;
import mc_moves_probabilities;
import mc_moves_move_types;
import mc_moves_cputime;
import reaction;
import reactions;
import cbmc;
import cbmc_chain_data;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import equation_of_states;
import thermostat;
import json;
import integrators;
import integrators_compute;
import integrators_update;

// construct System programmatically
/*! \brief Brief description.
 *         Brief description continued.
 *
 *  Detailed description starts here.
 */
System::System(size_t id, ForceField forcefield, std::optional<SimulationBox> box, double T, std::optional<double> P,
               double heliumVoidFraction, std::vector<Framework> f, std::vector<Component> c,
               std::vector<size_t> initialNumberOfMolecules, size_t numberOfBlocks,
               const MCMoveProbabilities& systemProbabilities, std::optional<size_t> sampleMoviesEvery)
    : systemId(id),
      temperature(T),
      pressure(P.value_or(0.0) / Units::PressureConversionFactor),
      input_pressure(P.value_or(0.0)),
      beta(1.0 / (Units::KB * T)),
      heliumVoidFraction(heliumVoidFraction),
      frameworkComponents(f),
      components(c),
      loadings(c.size()),
      swappableComponents(),
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
      atomPositions({}),
      moleculePositions({}),
      runningEnergies(),
      currentEnergyStatus(1, f.size(), c.size()),
      netChargePerComponent(c.size()),
      mc_moves_probabilities(systemProbabilities),
      mc_moves_statistics(),
      reactions(),
      tmmc(),
      averageEnergies(numberOfBlocks, 1, f.size(), c.size()),
      averageLoadings(numberOfBlocks, c.size()),
      averageEnthalpiesOfAdsorption(numberOfBlocks, c.size()),
      averageTemperature(numberOfBlocks),
      averageTranslationalTemperature(numberOfBlocks),
      averageRotationalTemperature(numberOfBlocks),
      averagePressure(numberOfBlocks),
      averageSimulationBox(numberOfBlocks)
{
  if (box.has_value())
  {
    simulationBox = box.value();
  }

  removeRedundantMoves();
  determineSwappableComponents();
  determineFractionalComponents();
  rescaleMoveProbabilities();
  rescaleMolarFractions();
  computeNumberOfPseudoAtoms();

  createFrameworks();
  determineSimulationBox();

  forceField.initializeEwaldParameters(simulationBox);

  CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      eik_x, eik_y, eik_z, eik_xy, forceField, simulationBox, double3(0.0, 0.0, 0.0), 1.0);

  precomputeTotalRigidEnergy();

  RandomNumber random(1400);

  translationalCenterOfMassConstraint = 0;
  translationalDegreesOfFreedom = 0;
  rotationalDegreesOfFreedom = 0;

  createInitialMolecules(random);

  equationOfState =
      EquationOfState(EquationOfState::Type::PengRobinson, EquationOfState::MultiComponentMixingRules::VanDerWaals, T,
                      P.value_or(0.0), simulationBox, heliumVoidFraction, components);

  averageEnthalpiesOfAdsorption.resize(swappableComponents.size());

  if (sampleMoviesEvery.has_value())
  {
    samplePDBMovie = SampleMovie(id, sampleMoviesEvery.value());
  }
}

System::System(size_t id, double T, std::optional<double> P, double heliumVoidFraction, std::vector<Framework> f,
               std::vector<Component> c)
    : systemId(id),
      temperature(T),
      pressure(P.value_or(0.0) / Units::PressureConversionFactor),
      input_pressure(P.value_or(0.0)),
      beta(1.0 / (Units::KB * T)),
      heliumVoidFraction(heliumVoidFraction),
      frameworkComponents(f),
      components(c)
{
}

void System::createFrameworks()
{
  netChargeFramework = 0.0;
  for (Framework& framework : frameworkComponents)
  {
    const std::vector<Atom>& atoms = framework.atoms;
    for (const Atom& atom : atoms)
    {
      atomPositions.push_back(atom);
      electricPotential.push_back(0.0);
      electricField.push_back(double3(0.0, 0.0, 0.0));
      electricFieldNew.push_back(double3(0.0, 0.0, 0.0));
    }
    numberOfFrameworkAtoms += atoms.size();
    numberOfRigidFrameworkAtoms += atoms.size();
    netChargeFramework += framework.netCharge;
    netCharge += framework.netCharge;
  }
}

std::optional<double> System::frameworkMass() const
{
  if (frameworkComponents.empty()) return std::nullopt;

  double mass = 0.0;
  for (const Framework& framework : frameworkComponents)
  {
    mass += framework.mass;
  }
  for (size_t i = 0; i < components.size(); ++i)
  {
    if (components[i].type == Component::Type::Cation)
    {
      mass += components[i].totalMass * static_cast<double>(numberOfIntegerMoleculesPerComponent[i]);
    }
  }
  return mass;
}

void System::determineSimulationBox()
{
  for (Framework& framework : frameworkComponents)
  {
    // For multiple framework, the simulation box is the union of the boxes
    simulationBox = max(simulationBox, framework.simulationBox.scaled(framework.numberOfUnitCells));
  }
}

void System::insertFractionalMolecule(size_t selectedComponent, [[maybe_unused]] const Molecule& molecule,
                                      std::vector<Atom> atoms, size_t moleculeId)
{
  double l = 0.0;
  for (Atom& atom : atoms)
  {
    atom.moleculeId = static_cast<uint16_t>(moleculeId);
    atom.groupId = uint8_t{components[selectedComponent].lambdaGC.computeDUdlambda};
    atom.setScaling(l);
  }
  std::vector<Atom>::const_iterator iterator = iteratorForMolecule(selectedComponent, 0);
  atomPositions.insert(iterator, atoms.begin(), atoms.end());

  std::vector<Molecule>::iterator moleculeIterator = indexForMolecule(selectedComponent, 0);
  moleculePositions.insert(moleculeIterator, molecule);

  electricPotential.resize(electricPotential.size() + atoms.size());
  electricField.resize(electricField.size() + atoms.size());
  electricFieldNew.resize(electricFieldNew.size() + atoms.size());

  numberOfMoleculesPerComponent[selectedComponent] += 1;

  netCharge += components[selectedComponent].netCharge;
  netChargeAdsorbates += components[selectedComponent].netCharge;
  netChargePerComponent[selectedComponent] += components[selectedComponent].netCharge;

  translationalDegreesOfFreedom += components[selectedComponent].translationalDegreesOfFreedom;
  rotationalDegreesOfFreedom += components[selectedComponent].rotationalDegreesOfFreedom;

  // set moleculesIds
  size_t index = numberOfFrameworkAtoms;  // indexOfFirstMolecule(selectedComponent);
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

/// Inserts a molecule into the vector of atoms.
///
/// Note: updates the numberOfMoleculesPerComponent, numberOfIntegerMoleculesPerComponent,
///       numberOfPseudoAtoms, totalNumberOfPseudoAtoms.
/// - Parameters:
///   - selectedComponent: the index of the component
///   - atoms: vector of atoms to be inserted
/// - returns:
void System::insertMolecule(size_t selectedComponent, [[maybe_unused]] const Molecule& molecule,
                            std::vector<Atom> atoms)
{
  std::vector<Atom>::const_iterator iterator =
      iteratorForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  atomPositions.insert(iterator, atoms.begin(), atoms.end());

  std::vector<Molecule>::iterator moleculeIterator =
      indexForMolecule(selectedComponent, numberOfMoleculesPerComponent[selectedComponent]);
  moleculePositions.insert(moleculeIterator, molecule);

  electricPotential.resize(electricPotential.size() + atoms.size());
  electricField.resize(electricField.size() + atoms.size());
  electricFieldNew.resize(electricFieldNew.size() + atoms.size());

  numberOfMoleculesPerComponent[selectedComponent] += 1;
  numberOfIntegerMoleculesPerComponent[selectedComponent] += 1;

  netCharge += components[selectedComponent].netCharge;
  netChargeAdsorbates += components[selectedComponent].netCharge;
  netChargePerComponent[selectedComponent] += components[selectedComponent].netCharge;

  translationalDegreesOfFreedom += components[selectedComponent].translationalDegreesOfFreedom;
  rotationalDegreesOfFreedom += components[selectedComponent].rotationalDegreesOfFreedom;

  // Update the number of pseudo atoms per type (used for tail-corrections)
  for (Atom& atom : atoms)
  {
    atom.moleculeId = static_cast<uint16_t>(numberOfMoleculesPerComponent[selectedComponent]);
    numberOfPseudoAtoms[selectedComponent][static_cast<size_t>(atom.type)] += 1;
    totalNumberOfPseudoAtoms[static_cast<size_t>(atom.type)] += 1;
  }

  size_t index = numberOfFrameworkAtoms;
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

void System::deleteMolecule(size_t selectedComponent, size_t selectedMolecule, const std::span<Atom> molecule)
{
  // Update the number of pseudo atoms per type (used for tail-corrections)
  for (const Atom& atom : molecule)
  {
    numberOfPseudoAtoms[selectedComponent][static_cast<size_t>(atom.type)] -= 1;
    totalNumberOfPseudoAtoms[static_cast<size_t>(atom.type)] -= 1;
  }

  std::vector<Atom>::const_iterator iterator = iteratorForMolecule(selectedComponent, selectedMolecule);
  atomPositions.erase(iterator, iterator + static_cast<std::vector<Atom>::difference_type>(molecule.size()));

  std::vector<Molecule>::iterator moleculeIterator = indexForMolecule(selectedComponent, selectedMolecule);
  moleculePositions.erase(moleculeIterator, moleculeIterator + 1);

  electricPotential.resize(electricPotential.size() - molecule.size());
  electricFieldNew.resize(electricFieldNew.size() - molecule.size());

  numberOfMoleculesPerComponent[selectedComponent] -= 1;
  numberOfIntegerMoleculesPerComponent[selectedComponent] -= 1;

  netCharge -= components[selectedComponent].netCharge;
  netChargeAdsorbates -= components[selectedComponent].netCharge;
  netChargePerComponent[selectedComponent] -= components[selectedComponent].netCharge;

  translationalDegreesOfFreedom -= components[selectedComponent].translationalDegreesOfFreedom;
  rotationalDegreesOfFreedom -= components[selectedComponent].rotationalDegreesOfFreedom;

  size_t index = numberOfFrameworkAtoms;
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

void System::checkMoleculeIds()
{
  std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();

  size_t index = 0;  // indexOfFirstMolecule(selectedComponent);
  for (size_t componentId = 0; componentId < components.size(); componentId++)
  {
    for (size_t i = 0; i < numberOfMoleculesPerComponent[componentId]; ++i)
    {
      for (size_t j = 0; j < components[componentId].atoms.size(); ++j)
      {
        if (moleculeAtoms[index].moleculeId != static_cast<uint32_t>(i))
        {
          throw std::runtime_error(std::format("Wrong molecule-id detected {} for component {} molecule {}\n",
                                               moleculeAtoms[index].moleculeId, componentId, i));
        }
        if (moleculeAtoms[index].componentId != static_cast<uint8_t>(componentId))
        {
          throw std::runtime_error(std::format("Wrong component-id detected {} for component {} molecule {}\n",
                                               moleculeAtoms[index].componentId, componentId, i));
        }
        ++index;
      }
    }
  }

  for (size_t componentId = 0; componentId < components.size(); componentId++)
  {
    if (numberOfGCFractionalMoleculesPerComponent_CFCMC[componentId] > 0)
    {
      size_t indexFractionalMolecule = indexOfGCFractionalMoleculesPerComponent_CFCMC(componentId);
      std::span<Atom> fractionalMolecule = spanOfMolecule(componentId, indexFractionalMolecule);

      for (const Atom& atom : fractionalMolecule)
      {
        if (components[componentId].lambdaGC.computeDUdlambda)
        {
          if (static_cast<size_t>(atom.groupId) == 0)
          {
            throw std::runtime_error(std::format("Wrong group-id detected! (0 where it should be 1)\n"));
          }
        }
        else
        {
          if (static_cast<size_t>(atom.groupId) == 1)
          {
            throw std::runtime_error(std::format("Wrong group-id detected! (1 where it should be 0)\n"));
          }
        }
      }
    }
  }
}

void System::createInitialMolecules([[maybe_unused]] RandomNumber& random)
{
  for (size_t componentId = 0; const Component& component : components)
  {
    if (component.swappable)
    {
      numberOfMoleculesPerComponent[componentId] = 0;
      for (size_t i = 0; i < numberOfFractionalMoleculesPerComponent[componentId]; ++i)
      {
        std::optional<ChainData> growData = std::nullopt;
        do
        {
          Component::GrowType growType = components[componentId].growType;
          growData = CBMC::growMoleculeSwapInsertion(
              random, frameworkComponents, components[componentId], hasExternalField, components, forceField,
              simulationBox, spanOfFrameworkAtoms(), spanOfMoleculeAtoms(), beta, growType,
              forceField.cutOffFrameworkVDW, forceField.cutOffMoleculeVDW, forceField.cutOffCoulomb, componentId,
              numberOfMoleculesPerComponent[componentId], 0.0, 1uz, numberOfTrialDirections);
        } while (!growData || growData->energies.potentialEnergy() > forceField.overlapCriteria);

        insertFractionalMolecule(componentId, growData->molecule, growData->atom, i);
      }
    }

    for (size_t i = 0; i < initialNumberOfMolecules[componentId]; ++i)
    {
      std::optional<ChainData> growData = std::nullopt;
      bool inside_blocked_pocket{false};
      do
      {
        do
        {
          Component::GrowType growType = components[componentId].growType;
          growData = CBMC::growMoleculeSwapInsertion(
              random, frameworkComponents, components[componentId], hasExternalField, components, forceField,
              simulationBox, spanOfFrameworkAtoms(), spanOfMoleculeAtoms(), beta, growType,
              forceField.cutOffFrameworkVDW, forceField.cutOffMoleculeVDW, forceField.cutOffCoulomb, componentId,
              numberOfMoleculesPerComponent[componentId], 1.0, 0uz, numberOfTrialDirections);

        } while (!growData || growData->energies.potentialEnergy() > forceField.overlapCriteria);

        std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());
        inside_blocked_pocket = insideBlockedPockets(components[componentId], newMolecule);

      } while (inside_blocked_pocket);

      insertMolecule(componentId, growData->molecule, growData->atom);
    }
    componentId++;
  }
}

bool System::insideBlockedPockets(const Component& component, std::span<const Atom> molecule_atoms) const
{
  for (const Framework& framework : frameworkComponents)
  {
    for (size_t i = 0; i != component.blockingPockets.size(); ++i)
    {
      double radius_squared = component.blockingPockets[i].w * component.blockingPockets[i].w;
      double3 pos =
          framework.simulationBox.cell *
          double3(component.blockingPockets[i].x, component.blockingPockets[i].y, component.blockingPockets[i].z);
      for (const Atom& atom : molecule_atoms)
      {
        double lambda = atom.scalingVDW;
        double3 dr = atom.position - pos;
        dr = framework.simulationBox.applyPeriodicBoundaryConditions(dr);
        if (dr.length_squared() < lambda * radius_squared)
        {
          return true;
        }
      }
    }
  }
  return false;
}

size_t System::randomMoleculeOfComponent(RandomNumber& random, size_t selectedComponent)
{
  return size_t(random.uniform() * static_cast<double>(numberOfMoleculesPerComponent[selectedComponent]));
}

size_t System::randomIntegerMoleculeOfComponent(RandomNumber& random, size_t selectedComponent)
{
  return numberOfFractionalMoleculesPerComponent[selectedComponent] +
         size_t(random.uniform() * static_cast<double>(numberOfIntegerMoleculesPerComponent[selectedComponent]));
}

std::vector<Atom>::iterator System::iteratorForMolecule(size_t selectedComponent, size_t selectedMolecule)
{
  size_t index{0};
  for (size_t i = 0; i < selectedComponent; ++i)
  {
    size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule + numberOfFrameworkAtoms;
  return atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(index);
}

std::vector<Molecule>::iterator System::indexForMolecule(size_t selectedComponent, size_t selectedMolecule)
{
  size_t index{0};
  for (size_t i = 0; i < selectedComponent; ++i)
  {
    index += numberOfMoleculesPerComponent[i];
  }
  index += selectedMolecule;
  return moleculePositions.begin() + static_cast<std::vector<Atom>::difference_type>(index);
}

size_t System::moleculeIndexOfComponent(size_t selectedComponent, size_t selectedMolecule)
{
  size_t index{0};
  for (size_t i = 0; i < selectedComponent; ++i)
  {
    index += numberOfMoleculesPerComponent[i];
  }
  index += selectedMolecule;
  return index;
}

std::span<const Atom> System::spanOfFrameworkAtoms() const
{
  return std::span(atomPositions.begin(), numberOfFrameworkAtoms);
}

std::span<Atom> System::spanOfFrameworkAtoms() { return std::span(atomPositions.begin(), numberOfFrameworkAtoms); }

std::span<const Atom> System::spanOfRigidFrameworkAtoms() const
{
  return std::span(atomPositions.begin(), numberOfFrameworkAtoms);
}

std::span<const Atom> System::spanOfFlexibleAtoms() const
{
  return std::span(atomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(numberOfFrameworkAtoms),
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

std::span<double> System::spanOfMoleculeElectricPotential()
{
  return std::span(
      electricPotential.begin() + static_cast<std::vector<double3>::difference_type>(numberOfFrameworkAtoms),
      electricPotential.end());
}

std::span<double3> System::spanOfMoleculeElectricField()
{
  return std::span(electricField.begin() + static_cast<std::vector<double3>::difference_type>(numberOfFrameworkAtoms),
                   electricField.end());
}

std::span<double3> System::spanOfMoleculeElectricFieldNew()
{
  return std::span(
      electricFieldNew.begin() + static_cast<std::vector<double3>::difference_type>(numberOfFrameworkAtoms),
      electricFieldNew.end());
}

std::span<Atom> System::spanOfMolecule(size_t selectedComponent, size_t selectedMolecule)
{
  size_t index{0};
  for (size_t i = 0; i < selectedComponent; ++i)
  {
    size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule;
  return std::span(&atomPositions[index + numberOfFrameworkAtoms], size);
}

const std::span<const Atom> System::spanOfMolecule(size_t selectedComponent, size_t selectedMolecule) const
{
  size_t index{0};
  for (size_t i = 0; i < selectedComponent; ++i)
  {
    size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule;
  return std::span(&atomPositions[index + numberOfFrameworkAtoms], size);
}

std::span<double3> System::spanElectricFieldNew(size_t selectedComponent, size_t selectedMolecule)
{
  size_t index{0};
  for (size_t i = 0; i < selectedComponent; ++i)
  {
    size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule;
  return std::span(&electricFieldNew[index + numberOfFrameworkAtoms], size);
}

const std::span<const double3> System::spanElectricFieldNew(size_t selectedComponent, size_t selectedMolecule) const
{
  size_t index{0};
  for (size_t i = 0; i < selectedComponent; ++i)
  {
    size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  size_t size = components[selectedComponent].atoms.size();
  index += size * selectedMolecule;
  return std::span(&electricFieldNew[index + numberOfFrameworkAtoms], size);
}

size_t System::indexOfFirstMolecule(size_t selectedComponent)
{
  size_t index{0};
  for (size_t i = 0; i < selectedComponent; ++i)
  {
    size_t size = components[i].atoms.size();
    index += size * numberOfMoleculesPerComponent[i];
  }
  return index + numberOfFrameworkAtoms;
}

void System::determineSwappableComponents()
{
  for (Component& component : components)
  {
    if (component.mc_moves_probabilities.getProbability(MoveTypes::Swap) > 0.0 ||
        component.mc_moves_probabilities.getProbability(MoveTypes::SwapCBMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(MoveTypes::SwapCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(MoveTypes::SwapCBCFCMC) > 0.0)
    {
      component.swappable = true;
    }

    if (component.mc_moves_probabilities.getProbability(MoveTypes::GibbsSwapCBMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(MoveTypes::GibbsSwapCFCMC) > 0.0)
    {
      component.swappable = true;
    }

    if (component.swappable)
    {
      swappableComponents.push_back(component.componentId);
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

    if (components[i].mc_moves_probabilities.getProbability(MoveTypes::SwapCFCMC) > 0.0 ||
        components[i].mc_moves_probabilities.getProbability(MoveTypes::WidomCFCMC) > 0.0 ||
        components[i].mc_moves_probabilities.getProbability(MoveTypes::SwapCBCFCMC) > 0.0 ||
        components[i].mc_moves_probabilities.getProbability(MoveTypes::WidomCBCFCMC) > 0.0)
    {
      numberOfFractionalMoleculesPerComponent[i] += 1;
      numberOfGCFractionalMoleculesPerComponent_CFCMC[i] = 1;
      components[i].hasFractionalMolecule = true;
    }

    // Gibbs
    if (components[i].mc_moves_probabilities.getProbability(MoveTypes::GibbsSwapCFCMC) > 0.0)
    {
      numberOfFractionalMoleculesPerComponent[i] += 1;
      numberOfGibbsFractionalMoleculesPerComponent_CFCMC[i] = 1;
      components[i].hasFractionalMolecule = true;
    }
  }

  for (size_t reactionId{0}; [[maybe_unused]] const Reaction& reaction : reactions.list)
  {
    for (size_t j = 0; j < components.size(); ++j)
    {
      numberOfReactionFractionalMoleculesPerComponent_CFCMC[reactionId][j] = 0;
    }

    // for ([[maybe_unused]]  size_t componentId{ 0 };  const size_t [[maybe_unused]]  stoichiometry :
    // reaction.reactantStoichiometry)
    //{
    //   //numberOfReactionFractionalMoleculesPerComponent[reactionId][componentId] += stoichiometry;
    //   ++componentId;
    // }
    // for ([[maybe_unused]] size_t componentId{ 0 };  const [[maybe_unused]] size_t stoichiometry :
    // reaction.productStoichiometry)
    //{
    //   //numberOfReactionFractionalMoleculesPerComponent[reactionId][componentId] += stoichiometry;
    //   ++componentId;
    // }
    ++reactionId;
  }
}

void System::rescaleMoveProbabilities()
{
  for (Component& component : components)
  {
    component.mc_moves_probabilities.join(mc_moves_probabilities);
  }
}

void System::removeRedundantMoves()
{
  for (Component& component : components)
  {
    component.mc_moves_probabilities.removeRedundantMoves();
  }
}

void System::optimizeMCMoves()
{
  mc_moves_statistics.optimizeMCMoves();
  for (Component& component : components)
  {
    component.mc_moves_statistics.optimizeMCMoves();
  }
}

void System::rescaleMolarFractions()
{
  double totalMolfraction = 0.0;
  double numberOfSwappableComponents = 0.0;
  for (const Component& component : components)
  {
    if (component.swappable)
    {
      totalMolfraction += component.molFraction;
      numberOfSwappableComponents += 1.0;
    }
  }

  if (totalMolfraction > 0.0)
  {
    for (Component& component : components)
    {
      if (component.swappable)
      {
        component.molFraction /= totalMolfraction;
      }
    }
  }
  else
  {
    for (Component& component : components)
    {
      if (component.swappable)
      {
        component.molFraction /= numberOfSwappableComponents;
      }
    }
  }
}

void System::computeNumberOfPseudoAtoms()
{
  for (size_t i = 0; i != components.size(); ++i)
  {
    std::fill(numberOfPseudoAtoms[i].begin(), numberOfPseudoAtoms[i].end(), 0);
  }
  std::fill(totalNumberOfPseudoAtoms.begin(), totalNumberOfPseudoAtoms.end(), 0);

  for (const Atom& atom : atomPositions)
  {
    size_t componentId = static_cast<size_t>(atom.componentId);
    size_t type = static_cast<size_t>(atom.type);
    numberOfPseudoAtoms[componentId][type] += 1;
    totalNumberOfPseudoAtoms[type] += 1;
  }
}

std::vector<Atom> System::randomConfiguration(RandomNumber& random, size_t selectedComponent,
                                              const std::span<const Atom> molecule)
{
  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());
  double3 position = simulationBox.randomPosition(random);
  size_t startingBead = components[selectedComponent].startingBead;
  for (size_t i = 0; i != molecule.size(); ++i)
  {
    copied_atoms[i].position =
        position + randomRotationMatrix * (molecule[i].position - molecule[startingBead].position);
  }
  return copied_atoms;
}

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
  // const size_t numberOfHelperThreads = pool.getThreadCount();

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

  for (size_t i = 0; const Component& c : components)
  {
    std::print(stream, "Component {:3d} ({})\n", c.componentId, c.name);
    std::print(stream, "-------------------------------------------------------------------------------\n");
    for (size_t index = 0; const size_t number_of_pseudo_atoms : numberOfPseudoAtoms[i])
    {
      std::print(stream, "    index {:3d} ({}): {} atoms\n", index, forceField.pseudoAtoms[index].name,
                 number_of_pseudo_atoms);
      ++index;
    }
    std::print(stream, "\n");
    ++i;
  }

  std::print(stream, "Total number of pseudo-atoms:\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (size_t index = 0; const size_t number_of_pseudo_atoms : totalNumberOfPseudoAtoms)
  {
    std::print(stream, "    index {:3d} ({}): {} atoms\n", index, forceField.pseudoAtoms[index].name,
               number_of_pseudo_atoms);
    ++index;
  }

  std::print(stream, "\n\n\n\n");

  return stream.str();
}

std::string System::writeInitializationStatusReport(size_t currentCycle, size_t numberOfCycles) const
{
  std::ostringstream stream;

  std::print(stream, "Initialization: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::print(stream, "{}\n", simulationBox.printStatus());
  std::print(stream, "net charge: {:12.8f}\n", netCharge);
  std::print(stream, "\n");

  for (size_t i = 0; const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {:3d} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n",
                 c.componentId, c.name, lambda, runningEnergies.dudlambda(lambda), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {:3d} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", c.componentId, c.name,
                 c.lambdaGC.lambdaValue(), occupancy, averageOccupancy);
    }
    std::print(stream, "    net charge: {:12.8f} [e]\n", netChargePerComponent[i]);
    ++i;
  }
  std::print(stream, "\n");

  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (const Component& c : components)
  {
    std::print(stream, "{}", loadings.printStatus(c, frameworkMass()));
  }
  std::print(stream, "\n");

  stream << runningEnergies.printMC();

  std::print(stream, "\n");

  return stream.str();
}

std::string System::writeEquilibrationStatusReportMC(size_t currentCycle, size_t numberOfCycles) const
{
  std::ostringstream stream;

  std::print(stream, "Equilibration: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::print(stream, "{}\n", simulationBox.printStatus());
  std::print(stream, "net charge: {:12.8f}\n", netCharge);
  std::print(stream, "\n");

  for (size_t i = 0; const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", c.componentId,
                 c.name, lambda, runningEnergies.dudlambda(lambda), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", c.componentId, c.name,
                 c.lambdaGC.lambdaValue(), occupancy, averageOccupancy);
    }
    std::print(stream, "    net charge: {:12.8f} [e]\n", netChargePerComponent[i]);
    ++i;
  }
  std::print(stream, "\n");

  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (const Component& c : components)
  {
    std::print(stream, "{}", loadings.printStatus(c, frameworkMass()));
  }
  std::print(stream, "\n");

  stream << runningEnergies.printMC();

  std::print(stream, "\n");

  return stream.str();
}

std::string System::writeEquilibrationStatusReportMD(size_t currentCycle, size_t numberOfCycles) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;

  std::print(stream, "Equilibration: Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::print(stream, "{}\n", simulationBox.printStatus());
  double3 linear_momentum = Integrators::computeLinearMomentum(moleculePositions);
  std::print(stream, "Linear momentum: {:12.8f} {:12.8f} {:12.8f}\n", linear_momentum.x, linear_momentum.y,
             linear_momentum.z);
  double3 com_velocity = Integrators::computeCenterOfMassVelocity(moleculePositions);
  std::print(stream, "Center of mass velocity: {:12.8f} {:12.8f} {:12.8f}\n", com_velocity.x, com_velocity.y,
             com_velocity.z);
  double3 com = Integrators::computeCenterOfMass(moleculePositions);
  std::print(stream, "Center of mass: {:12.8f} {:12.8f} {:12.8f}\n", com.x, com.y, com.z);
  std::print(stream, "Net charge: {:12.8f}\n", netCharge);
  std::print(stream, "\n");

  double translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculePositions);
  double translationalTemperature =
      2.0 * translationalKineticEnergy /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint));
  double rotationalKineticEnergy = Integrators::computeRotationalKineticEnergy(moleculePositions, components);
  double rotationalTemperature =
      2.0 * rotationalKineticEnergy / (Units::KB * static_cast<double>(rotationalDegreesOfFreedom));
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
  std::print(stream, "Conserved energy:   {: .6e}\n", conservedEnergy);
  double drift = std::abs(conv * (conservedEnergy - referenceEnergy) / referenceEnergy);
  std::print(stream, "Drift: {:.6e} Average drift: {:.6e}\n\n", drift,
             accumulatedDrift / static_cast<double>(std::max(currentCycle, 1uz)));

  stream << runningEnergies.printMD();

  std::print(stream, "\n");

  for (size_t i = 0; const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", c.componentId,
                 c.name, lambda, runningEnergies.dudlambda(lambda), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", c.componentId, c.name,
                 c.lambdaGC.lambdaValue(), occupancy, averageOccupancy);
    }
    std::print(stream, "    net charge: {:12.8f} [e]\n", netChargePerComponent[i]);
    ++i;
  }
  std::print(stream, "\n");

  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  for (const Component& c : components)
  {
    std::print(stream, "{}", loadings.printStatus(c, frameworkMass()));
  }
  std::print(stream, "\n");

  return stream.str();
}

std::string System::writeProductionStatusReportMC(size_t currentCycle, size_t numberOfCycles) const
{
  std::ostringstream stream;

  std::print(stream, "Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::pair<SimulationBox, SimulationBox> simulationBoxData = averageSimulationBox.averageSimulationBox();
  std::print(stream, "{}\n", simulationBox.printStatus(simulationBoxData.first, simulationBoxData.second));
  std::print(stream, "net charge: {:12.8f}\n", netCharge);
  std::print(stream, "\n");

  for (size_t i = 0; const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", c.componentId,
                 c.name, lambda, runningEnergies.dudlambda(lambda), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", c.componentId, c.name,
                 c.lambdaGC.lambdaValue(), occupancy, averageOccupancy);
    }
    std::print(stream, "    net charge: {:12.8f} [e]\n", netChargePerComponent[i]);
    ++i;
  }
  std::print(stream, "\n");

  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::pair<Loadings, Loadings> loadingData = averageLoadings.averageLoading();
  for (const Component& c : components)
  {
    std::print(stream, "{}", loadings.printStatus(c, loadingData.first, loadingData.second, frameworkMass()));
  }
  std::print(stream, "\n");
  double conv = Units::EnergyToKelvin;

  std::pair<double3x3, double3x3> currentPressureTensor = averagePressure.averagePressureTensor();

  switch (Units::unitSystem)
  {
    case Units::System::RASPA:
    {
      double3x3 pressureTensor = 1e-5 * Units::PressureConversionFactor * currentPressureTensor.first;
      double3x3 pressureTensorError = 1e-5 * Units::PressureConversionFactor * currentPressureTensor.second;
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
                 1e-5 * Units::PressureConversionFactor * p.first, 1e-5 * Units::PressureConversionFactor * p.second);
    }
    break;
    case Units::System::ReducedUnits:
    {
      double3x3 pressureTensor = currentPressureTensor.first;
      double3x3 pressureTensorError = currentPressureTensor.second;
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
      std::pair<double, double> idealGasPressure = averagePressure.averageIdealGasPressure();
      std::pair<double, double> excessPressure = averagePressure.averageExcessPressure();
      std::pair<double, double> p = averagePressure.averagePressure();
      std::print(stream, "Ideal-gas pressure:  {: .6e} +/ {:.6e} [{}]\n", idealGasPressure.first,
                 idealGasPressure.second, Units::unitOfPressureString);
      std::print(stream, "Excess pressure:     {: .6e} +/ {:.6e} [{}]\n", excessPressure.first, excessPressure.second,
                 Units::unitOfPressureString);
      std::print(stream, "Pressure:            {: .6e} +/ {:.6e} [{}]\n\n", p.first, p.second,
                 Units::unitOfPressureString);
    }
    break;
  }

  std::pair<EnergyStatus, EnergyStatus> energyData = averageEnergies.averageEnergy();
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

  std::print(stream, "\n");

  return stream.str();
}

std::string System::writeProductionStatusReportMD(size_t currentCycle, size_t numberOfCycles) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;

  std::print(stream, "Current cycle: {} out of {}\n", currentCycle, numberOfCycles);
  std::print(stream, "===============================================================================\n\n");

  std::pair<SimulationBox, SimulationBox> simulationBoxData = averageSimulationBox.averageSimulationBox();
  std::print(stream, "{}", simulationBox.printStatus(simulationBoxData.first, simulationBoxData.second));
  std::print(stream, "\n");

  double3 linear_momentum = Integrators::computeLinearMomentum(moleculePositions);
  std::print(stream, "Linear momentum: {:12.8f} {:12.8f} {:12.8f}\n", linear_momentum.x, linear_momentum.y,
             linear_momentum.z);
  double3 com_velocity = Integrators::computeCenterOfMassVelocity(moleculePositions);
  std::print(stream, "Center of mass velocity: {:12.8f} {:12.8f} {:12.8f}\n", com_velocity.x, com_velocity.y,
             com_velocity.z);
  double3 com = Integrators::computeCenterOfMass(moleculePositions);
  std::print(stream, "Center of mass: {:12.8f} {:12.8f} {:12.8f}\n", com.x, com.y, com.z);
  std::print(stream, "Net charge: {:12.8f}\n", netCharge);
  std::print(stream, "Time run: {:g} [ps]  {:g} [ns]\n\n", static_cast<double>(currentCycle) * timeStep,
             static_cast<double>(currentCycle) * timeStep / 1000.0);

  double translational_kinetic_energy = Integrators::computeTranslationalKineticEnergy(moleculePositions);
  double translational_temperature =
      2.0 * translational_kinetic_energy /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint));
  double rotational_kinetic_energy = Integrators::computeRotationalKineticEnergy(moleculePositions, components);
  double rotational_temperature =
      2.0 * rotational_kinetic_energy / (Units::KB * static_cast<double>(rotationalDegreesOfFreedom));
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
  std::print(stream, "Conserved energy:   {: .6e}\n", conservedEnergy);
  double drift = std::abs(conv * (conservedEnergy - referenceEnergy) / referenceEnergy);
  std::print(stream, "Drift: {:.6e} Average drift: {:.6e}\n\n", drift,
             accumulatedDrift / static_cast<double>(std::max(currentCycle, 1uz)));

  std::pair<EnergyStatus, EnergyStatus> energyData = averageEnergies.averageEnergy();
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
  for (size_t i = 0; const Component& c : components)
  {
    double occupancy = static_cast<double>(containsTheFractionalMolecule);
    double averageOccupancy = c.lambdaGC.occupancy();
    double lambda = c.lambdaGC.lambdaValue();

    if (c.lambdaGC.computeDUdlambda)
    {
      std::print(stream, "component {} ({}) lambda: {: g} dUdlambda: {: g} occupancy: {: g} ({:3f})\n", c.componentId,
                 c.name, lambda, runningEnergies.dudlambda(lambda), occupancy, averageOccupancy);
    }
    else
    {
      std::print(stream, "component {} ({}) lambda: {: g} occupancy: {: g} ({:3f})\n", c.componentId, c.name,
                 c.lambdaGC.lambdaValue(), occupancy, averageOccupancy);
    }
    std::print(stream, "    net charge: {:12.8f} [e]\n", netChargePerComponent[i]);
    ++i;
  }
  std::print(stream, "\n");

  std::print(stream, "Amount of molecules per component :\n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::pair<Loadings, Loadings> loadingData = averageLoadings.averageLoading();
  for (const Component& c : components)
  {
    std::print(stream, "{}", loadings.printStatus(c, loadingData.first, loadingData.second, frameworkMass()));
  }
  std::print(stream, "\n");

  std::pair<double3x3, double3x3> currentPressureTensor = averagePressure.averagePressureTensor();
  double3x3 pressureTensor = 1e-5 * Units::PressureConversionFactor * currentPressureTensor.first;
  double3x3 pressureTensorError = 1e-5 * Units::PressureConversionFactor * currentPressureTensor.second;
  std::print(stream, "Average pressure tensor: \n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.ax, pressureTensor.bx,
             pressureTensor.cx, pressureTensorError.ax, pressureTensorError.bx, pressureTensorError.cx);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.ay, pressureTensor.by,
             pressureTensor.cy, pressureTensorError.ay, pressureTensorError.by, pressureTensorError.cy);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.az, pressureTensor.bz,
             pressureTensor.cz, pressureTensorError.az, pressureTensorError.bz, pressureTensorError.cz);
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
             1e-5 * Units::PressureConversionFactor * p.first, 1e-5 * Units::PressureConversionFactor * p.second);

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
  for (const Framework& component : frameworkComponents)
  {
    std::print(stream, "{}", component.printStatus(forceField));
  }
  for (const Component& component : components)
  {
    std::print(stream, "{}", component.printStatus(forceField));
  }
  std::print(stream, "\n\n\n\n");

  return stream.str();
}

nlohmann::json System::jsonComponentStatus() const
{
  nlohmann::json status;
  for (const Framework& component : frameworkComponents)
  {
    status[component.name] = component.jsonStatus();
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

void System::sampleProperties(size_t currentBlock, size_t currentCycle)
{
  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  double w = weight();

  averageSimulationBox.addSample(currentBlock, simulationBox, w);

  double translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculePositions);
  double translationalTemperature =
      2.0 * translationalKineticEnergy /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint));
  averageTranslationalTemperature.addSample(currentBlock, translationalTemperature, w);

  double rotationalKineticEnergy = Integrators::computeRotationalKineticEnergy(moleculePositions, components);
  double rotationalTemperature =
      2.0 * rotationalKineticEnergy / (Units::KB * static_cast<double>(rotationalDegreesOfFreedom));
  averageRotationalTemperature.addSample(currentBlock, rotationalTemperature, w);

  double overallTemperature =
      2.0 * (translationalKineticEnergy + rotationalKineticEnergy) /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint +
                                       rotationalDegreesOfFreedom));
  averageTemperature.addSample(currentBlock, overallTemperature, w);

  loadings = Loadings(components.size(), numberOfIntegerMoleculesPerComponent, simulationBox);
  averageLoadings.addSample(currentBlock, loadings, w);

  EnthalpyOfAdsorptionTerms enthalpyTerms = EnthalpyOfAdsorptionTerms(
      swappableComponents, numberOfIntegerMoleculesPerComponent, runningEnergies.potentialEnergy(), temperature);
  averageEnthalpiesOfAdsorption.addSample(currentBlock, enthalpyTerms, w);

  size_t numberOfMolecules =
      std::reduce(numberOfIntegerMoleculesPerComponent.begin(), numberOfIntegerMoleculesPerComponent.end());
  double currentIdealPressure = static_cast<double>(numberOfMolecules) / (beta * simulationBox.volume);

  averagePressure.addSample(currentBlock, currentIdealPressure, currentExcessPressureTensor, w);

  for (Component& component : components)
  {
    double componentDensity =
        static_cast<double>(numberOfIntegerMoleculesPerComponent[component.componentId]) / simulationBox.volume;

    double lambda = component.lambdaGC.lambdaValue();
    double dudlambda = runningEnergies.dudlambda(lambda);
    component.lambdaGC.sampleHistogram(currentBlock, componentDensity, dudlambda, containsTheFractionalMolecule, w);

    component.averageRosenbluthWeights.addDensitySample(currentBlock, componentDensity, w);
  }

  if (samplePDBMovie.has_value())
  {
    samplePDBMovie->update(forceField, systemId, simulationBox, spanOfMoleculeAtoms(), currentCycle);
  }

  if (propertyConventionalRadialDistributionFunction.has_value())
  {
    propertyConventionalRadialDistributionFunction->sample(simulationBox, spanOfFrameworkAtoms(), spanOfMoleculeAtoms(),
                                                           currentCycle, currentBlock);
  }

  if (propertyRadialDistributionFunction.has_value())
  {
    precomputeTotalGradients();
    Integrators::updateCenterOfMassAndQuaternionGradients(moleculePositions, spanOfMoleculeAtoms(), components);
    propertyRadialDistributionFunction->sample(simulationBox, spanOfFrameworkAtoms(), moleculePositions,
                                               spanOfMoleculeAtoms(), currentCycle, currentBlock);
  }

  if (averageEnergyHistogram.has_value())
  {
    averageEnergyHistogram->addSample(
        currentBlock, currentCycle,
        {runningEnergies.potentialEnergy(), runningEnergies.frameworkMoleculeVDW + runningEnergies.moleculeMoleculeVDW,
         runningEnergies.frameworkMoleculeCharge + runningEnergies.moleculeMoleculeCharge +
             runningEnergies.ewald_fourier + runningEnergies.ewald_self + runningEnergies.ewald_exclusion,
         runningEnergies.polarization},
        w);
  }

  if (averageNumberOfMoleculesHistogram.has_value())
  {
    averageNumberOfMoleculesHistogram->addSample(currentBlock, currentCycle, numberOfIntegerMoleculesPerComponent, w);
  }

  if (propertyMSD.has_value())
  {
    propertyMSD->addSample(currentCycle, components, numberOfMoleculesPerComponent, moleculePositions);
  }

  if (propertyVACF.has_value())
  {
    propertyVACF->addSample(currentCycle, components, numberOfMoleculesPerComponent, moleculePositions);
  }

  if (propertyDensityGrid.has_value())
  {
    propertyDensityGrid->sample(frameworkComponents, simulationBox, spanOfMoleculeAtoms(), currentCycle);
  }

  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();

  mc_moves_cputime.propertySampling += (t2 - t1);
}

void System::writeCPUTimeStatistics(std::ostream& stream) const
{
  std::print(stream, "Sampling properties:        {:14f} [s]\n", mc_moves_cputime.propertySampling.count());
  std::print(stream, "Pressure computation:       {:14f} [s]\n\n", mc_moves_cputime.energyPressureComputation.count());

  for (size_t componentId = 0; const Component& component : components)
  {
    std::print(stream, "{}", component.mc_moves_cputime.writeMCMoveCPUTimeStatistics(componentId, component.name));
  }
}

std::pair<std::vector<Molecule>, std::vector<Atom>> System::scaledCenterOfMassPositions(double scale) const
{
  std::vector<Molecule> scaledMolecules;
  std::vector<Atom> scaledAtoms;
  scaledMolecules.reserve(scaledMolecules.size());
  scaledAtoms.reserve(atomPositions.size());

  for (size_t componentId = 0; componentId < components.size(); ++componentId)
  {
    for (size_t i = 0; i < numberOfMoleculesPerComponent[componentId]; ++i)
    {
      const std::span<const Atom> span = spanOfMolecule(componentId, i);
      const Molecule& molecule = moleculePositions[i];

      double totalMass = 0.0;
      double3 com(0.0, 0.0, 0.0);
      for (const Atom& atom : span)
      {
        double mass = forceField.pseudoAtoms[static_cast<size_t>(atom.type)].mass;
        com += mass * atom.position;
        totalMass += mass;
      }
      com = com / totalMass;

      double3 d = com * (scale - 1.0);
      scaledMolecules.push_back({com * scale, molecule.orientation, molecule.mass, componentId, span.size()});

      // create copy
      for (Atom atom : span)
      {
        atom.position += d;
        scaledAtoms.push_back(atom);
      }
    }
  }
  return {scaledMolecules, scaledAtoms};
}

inline std::pair<EnergyStatus, double3x3> pair_acc(const std::pair<EnergyStatus, double3x3>& lhs,
                                                   const std::pair<EnergyStatus, double3x3>& rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

void System::precomputeTotalRigidEnergy() noexcept
{
  Interactions::precomputeEwaldFourierRigid(eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, forceField,
                                            simulationBox, spanOfRigidFrameworkAtoms());
}

void System::precomputeTotalGradients() noexcept
{
  runningEnergies = Integrators::updateGradients(spanOfMoleculeAtoms(), spanOfFrameworkAtoms(), forceField,
                                                 simulationBox, components, eik_x, eik_y, eik_z, eik_xy, totalEik,
                                                 fixedFrameworkStoredEik, numberOfMoleculesPerComponent);
}

RunningEnergy System::computeTotalEnergies() noexcept
{
  RunningEnergy runningEnergy{};

  if (fixedFrameworkStoredEik.empty())
  {
    precomputeTotalRigidEnergy();
  }

  std::span<const Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();

  if (forceField.computePolarization)
  {
    std::span<double3> moleculeElectricField = spanOfMoleculeElectricField();

    std::fill(moleculeElectricField.begin(), moleculeElectricField.end(), double3(0.0, 0.0, 0.0));

    RunningEnergy frameworkMoleculeEnergy = Interactions::computeFrameworkMoleculeElectricField(
        forceField, simulationBox, moleculeElectricField, frameworkAtomPositions, moleculeAtomPositions);

    RunningEnergy intermolecularEnergy = Interactions::computeInterMolecularElectricField(
        forceField, simulationBox, moleculeElectricField, moleculeAtomPositions);

    RunningEnergy frameworkMoleculeTailEnergy = Interactions::computeFrameworkMoleculeTailEnergy(
        forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions);
    RunningEnergy intermolecularTailEnergy =
        Interactions::computeInterMolecularTailEnergy(forceField, simulationBox, moleculeAtomPositions);

    RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierElectricField(
        eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, forceField, simulationBox,
        moleculeElectricField, components, numberOfMoleculesPerComponent, moleculeAtomPositions);

    RunningEnergy polarizationEnergy = computePolarizationEnergy();

    return frameworkMoleculeEnergy + intermolecularEnergy + frameworkMoleculeTailEnergy + intermolecularTailEnergy +
           ewaldEnergy + polarizationEnergy;
  }
  else
  {
    RunningEnergy frameworkMoleculeEnergy = Interactions::computeFrameworkMoleculeEnergy(
        forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions);
    RunningEnergy intermolecularEnergy =
        Interactions::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtomPositions);

    RunningEnergy frameworkMoleculeTailEnergy = Interactions::computeFrameworkMoleculeTailEnergy(
        forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions);
    RunningEnergy intermolecularTailEnergy =
        Interactions::computeInterMolecularTailEnergy(forceField, simulationBox, moleculeAtomPositions);

    RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierEnergy(
        eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, forceField, simulationBox, components,
        numberOfMoleculesPerComponent, moleculeAtomPositions);

    return frameworkMoleculeEnergy + intermolecularEnergy + frameworkMoleculeTailEnergy + intermolecularTailEnergy +
           ewaldEnergy;
  }
}

RunningEnergy System::computePolarizationEnergy() noexcept
{
  RunningEnergy energy{};

  std::span<const Atom> moleculeAtomPositions = spanOfMoleculeAtoms();
  std::span<double3> moleculeElectricField = spanOfMoleculeElectricField();

  for (size_t i = 0; i < moleculeAtomPositions.size(); ++i)
  {
    size_t type = moleculeAtomPositions[i].type;
    double polarizability = forceField.pseudoAtoms[type].polarizability / Units::CoulombicConversionFactor;
    energy.polarization -= 0.5 * polarizability * double3::dot(moleculeElectricField[i], moleculeElectricField[i]);
  }

  return energy;
}

void System::computeTotalElectricPotential() noexcept
{
  if (fixedFrameworkStoredEik.empty())
  {
    precomputeTotalRigidEnergy();
  }

  std::span<Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();
  std::span<double> moleculeElectricPotential = spanOfMoleculeElectricPotential();

  std::fill(moleculeElectricPotential.begin(), moleculeElectricPotential.end(), 0.0);

  Interactions::computeInterMolecularElectricPotential(forceField, simulationBox, moleculeElectricPotential,
                                                       moleculeAtomPositions);

  Interactions::computeFrameworkMoleculeElectricPotential(forceField, simulationBox, moleculeElectricPotential,
                                                          frameworkAtomPositions, moleculeAtomPositions);

  Interactions::computeEwaldFourierElectricPotential(eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik,
                                                     moleculeElectricPotential, forceField, simulationBox, components,
                                                     numberOfMoleculesPerComponent, moleculeAtomPositions);
}

void System::computeTotalElectricField() noexcept
{
  if (fixedFrameworkStoredEik.empty())
  {
    precomputeTotalRigidEnergy();
  }

  std::span<Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();
  std::span<double3> moleculeElectricField = spanOfMoleculeElectricField();

  std::fill(moleculeElectricField.begin(), moleculeElectricField.end(), double3(0.0, 0.0, 0.0));

  Interactions::computeInterMolecularElectricField(forceField, simulationBox, moleculeElectricField,
                                                   moleculeAtomPositions);

  Interactions::computeFrameworkMoleculeElectricField(forceField, simulationBox, moleculeElectricField,
                                                      frameworkAtomPositions, moleculeAtomPositions);

  Interactions::computeEwaldFourierElectricField(eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik,
                                                 forceField, simulationBox, moleculeElectricField, components,
                                                 numberOfMoleculesPerComponent, moleculeAtomPositions);
}

std::pair<EnergyStatus, double3x3> System::computeMolecularPressure() noexcept
{
  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeFrameworkMoleculeEnergyStrainDerivative(
      forceField, frameworkComponents, components, simulationBox, spanOfFrameworkAtoms(), spanOfMoleculeAtoms());

  pressureInfo = pair_acc(pressureInfo, Interactions::computeInterMolecularEnergyStrainDerivative(
                                            forceField, components, simulationBox, spanOfMoleculeAtoms()));

  pressureInfo = pair_acc(
      pressureInfo, Interactions::computeEwaldFourierEnergyStrainDerivative(
                        eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, forceField, simulationBox,
                        frameworkComponents, components, numberOfMoleculesPerComponent, spanOfMoleculeAtoms(),
                        CoulombicFourierEnergySingleIon, netChargeFramework, netChargePerComponent));

  pressureInfo.first.sumTotal();

  double pressureTailCorrection = 0.0;
  double preFactor = 2.0 * std::numbers::pi / simulationBox.volume;
  for (std::vector<Atom>::iterator it1 = atomPositions.begin(); it1 != atomPositions.end(); ++it1)
  {
    size_t typeA = static_cast<size_t>(it1->type);
    double scalingVDWA = it1->scalingVDW;

    pressureTailCorrection += scalingVDWA * scalingVDWA * preFactor * forceField(typeA, typeA).tailCorrectionPressure;

    for (std::vector<Atom>::iterator it2 = it1 + 1; it2 != atomPositions.end(); ++it2)
    {
      size_t typeB = static_cast<size_t>(it2->type);
      double scalingVDWB = it2->scalingVDW;

      pressureTailCorrection +=
          scalingVDWA * scalingVDWB * 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionPressure;
    }
  }

  pressureInfo.second.ax -= pressureTailCorrection;
  pressureInfo.second.by -= pressureTailCorrection;
  pressureInfo.second.cz -= pressureTailCorrection;

  // Correct rigid molecule contribution using the constraints forces
  double3x3 correctionTerm{};
  for (size_t componentId = 0; componentId < components.size(); ++componentId)
  {
    if (components[componentId].rigid)
    {
      for (size_t i = 0; i < numberOfMoleculesPerComponent[componentId]; ++i)
      {
        std::span<Atom> span = spanOfMolecule(componentId, i);

        double totalMass = 0.0;
        double3 com(0.0, 0.0, 0.0);
        for (const Atom& atom : span)
        {
          double mass = forceField.pseudoAtoms[static_cast<size_t>(atom.type)].mass;
          com += mass * atom.position;
          totalMass += mass;
        }
        com = com / totalMass;

        for (const Atom& atom : span)
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

  double temp = 0.5 * (pressureInfo.second.ay + pressureInfo.second.bx);
  pressureInfo.second.ay = pressureInfo.second.bx = temp;
  temp = 0.5 * (pressureInfo.second.az + pressureInfo.second.cx);
  pressureInfo.second.az = pressureInfo.second.cx = temp;
  temp = 0.5 * (pressureInfo.second.bz + pressureInfo.second.cy);
  pressureInfo.second.bz = pressureInfo.second.cy = temp;

  return pressureInfo;
}

void System::checkCartesianPositions()
{
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();

  size_t index{};
  for (Molecule& molecule : moleculePositions)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
    if (components[molecule.componentId].rigid)
    {
      simd_quatd q = molecule.orientation;
      double3x3 M = double3x3::buildRotationMatrixInverse(q);

      for (size_t i = 0; i != span.size(); i++)
      {
        double3 expandedPosition =
            molecule.centerOfMassPosition + M * components[molecule.componentId].atoms[i].position;
        if ((std::abs(span[i].position.x - expandedPosition.x) > 1e-5) ||
            (std::abs(span[i].position.y - expandedPosition.y) > 1e-5) ||
            (std::abs(span[i].position.z - expandedPosition.z) > 1e-5))
        {
          throw std::runtime_error(
              std::format("Difference detected between atom position ({} {} {}) and position generated from "
                          "quaternion ({} {} {})\n",
                          span[i].position.x, span[i].position.y, span[i].position.z, expandedPosition.x,
                          expandedPosition.y, expandedPosition.z));
        }
      }
    }
    index += molecule.numberOfAtoms;
  }
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

inline std::string formatMoveStatistics(const std::string name, const MoveStatistics<double3>& move)
{
  std::ostringstream stream;

  std::print(stream, "{} total:        {:10} {:10} {:10}\n", name, move.counts.x, move.counts.y, move.counts.z);
  std::print(stream, "{} constructed:  {:10} {:10} {:10}\n", name, move.constructed.x, move.constructed.y,
             move.constructed.z);
  std::print(stream, "{} accepted:     {:10} {:10} {:10}\n", name, move.accepted.x, move.accepted.y, move.accepted.z);
  std::print(
      stream, "{} fraction:     {:10f} {:10f} {:10f}\n", name, move.accepted.x / std::max(1.0, double(move.counts.x)),
      move.accepted.y / std::max(1.0, double(move.counts.y)), move.accepted.z / std::max(1.0, double(move.counts.z)));
  std::print(stream, "{} max-change:   {:10f} {:10f} {:10f}\n\n", name, move.maxChange.x, move.maxChange.y,
             move.maxChange.z);

  return stream.str();
}

std::string System::writeMCMoveStatistics() const
{
  std::ostringstream stream;

  std::print(stream, "{}", mc_moves_statistics.writeMCMoveStatistics());
  for (size_t componentId = 0; const Component& component : components)
  {
    std::print(stream, "Component {} [{}]\n", componentId, component.name);

    std::print(stream, "{}", component.mc_moves_statistics.writeMCMoveStatistics());

    if (component.hasFractionalMolecule)
    {
      double imposedChemicalPotential = std::log(beta * component.molFraction * pressure) / beta;
      double imposedFugacity = component.molFraction * pressure;

      std::print(stream, "{}",
                 component.lambdaGC.writeAveragesStatistics(beta, imposedChemicalPotential, imposedFugacity));
      std::print(stream, "{}",
                 component.lambdaGC.writeDUdLambdaStatistics(beta, imposedChemicalPotential, imposedFugacity));
    }

    if (component.mc_moves_probabilities.getProbability(MoveTypes::Widom) > 0.0)
    {
      double imposedChemicalPotential = std::log(beta * component.molFraction * pressure) / beta;
      double imposedFugacity = component.molFraction * pressure;
      std::print(
          stream, "{}",
          component.averageRosenbluthWeights.writeAveragesStatistics(beta, imposedChemicalPotential, imposedFugacity));
    }

    ++componentId;
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

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const System& s)
{
  archive << s.versionNumber;

  archive << s.systemId;
  archive << s.temperature;
  archive << s.pressure;
  archive << s.input_pressure;
  archive << s.beta;
  archive << s.heliumVoidFraction;
  archive << s.numberOfFrameworks;
  archive << s.numberOfFrameworkAtoms;
  archive << s.numberOfRigidFrameworkAtoms;
  archive << s.frameworkComponents;
  archive << s.components;
  archive << s.equationOfState;
  archive << s.thermostat;
  archive << s.loadings;
  archive << s.swappableComponents;
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
  archive << s.translationalCenterOfMassConstraint;
  archive << s.translationalDegreesOfFreedom;
  archive << s.rotationalDegreesOfFreedom;
  archive << s.timeStep;
  archive << s.simulationBox;
  archive << s.atomPositions;
  archive << s.moleculePositions;
  archive << s.electricPotential;
  archive << s.electricField;
  archive << s.electricFieldNew;
  archive << s.conservedEnergy;
  archive << s.referenceEnergy;
  archive << s.rigidEnergies;
  archive << s.runningEnergies;
  archive << s.currentExcessPressureTensor;
  archive << s.currentEnergyStatus;
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
  archive << s.netChargeFramework;
  archive << s.netChargeAdsorbates;
  archive << s.netChargePerComponent;
  archive << s.mc_moves_probabilities;
  archive << s.mc_moves_statistics;
  archive << s.mc_moves_cputime;
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
  archive << s.averageEnergies;
  archive << s.averageLoadings;
  archive << s.averageEnthalpiesOfAdsorption;
  archive << s.averageTemperature;
  archive << s.averageTranslationalTemperature;
  archive << s.averageRotationalTemperature;
  archive << s.averagePressure;
  archive << s.averageSimulationBox;
  archive << s.averageEnergyHistogram;
  archive << s.propertyConventionalRadialDistributionFunction;
  // archive << s.propertyRadialDistributionFunction;
  // archive << s.propertyDensityGrid;

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, System& s)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > s.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'System' at line {} in file {}\n", location.line(), location.file_name()));
  }

  archive >> s.systemId;
  archive >> s.temperature;
  archive >> s.pressure;
  archive >> s.input_pressure;
  archive >> s.beta;
  archive >> s.heliumVoidFraction;
  archive >> s.numberOfFrameworks;
  archive >> s.numberOfFrameworkAtoms;
  archive >> s.numberOfRigidFrameworkAtoms;
  archive >> s.frameworkComponents;
  archive >> s.components;
  archive >> s.equationOfState;
  archive >> s.thermostat;
  archive >> s.loadings;
  archive >> s.swappableComponents;
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
  archive >> s.translationalCenterOfMassConstraint;
  archive >> s.translationalDegreesOfFreedom;
  archive >> s.rotationalDegreesOfFreedom;
  archive >> s.timeStep;
  archive >> s.simulationBox;
  archive >> s.atomPositions;
  archive >> s.moleculePositions;
  archive >> s.electricPotential;
  archive >> s.electricField;
  archive >> s.electricFieldNew;
  archive >> s.conservedEnergy;
  archive >> s.referenceEnergy;
  archive >> s.rigidEnergies;
  archive >> s.runningEnergies;
  archive >> s.currentExcessPressureTensor;
  archive >> s.currentEnergyStatus;
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
  archive >> s.netChargeFramework;
  archive >> s.netChargeAdsorbates;
  archive >> s.netChargePerComponent;
  archive >> s.mc_moves_probabilities;
  archive >> s.mc_moves_statistics;
  archive >> s.mc_moves_cputime;
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
  archive >> s.averageEnergies;
  archive >> s.averageLoadings;
  archive >> s.averageEnthalpiesOfAdsorption;
  archive >> s.averageTemperature;
  archive >> s.averageTranslationalTemperature;
  archive >> s.averageRotationalTemperature;
  archive >> s.averagePressure;
  archive >> s.averageSimulationBox;
  archive >> s.averageEnergyHistogram;
  archive >> s.propertyConventionalRadialDistributionFunction;
  // archive >> s.propertyRadialDistributionFunction;
  // archive >> s.propertyDensityGrid;

  return archive;
}

void System::writeRestartFile()
{
  nlohmann::json j;

  j["simulationBox"] = simulationBox;
  j["atomPositions"] = atomPositions;
  j["moleculePositions"] = moleculePositions;

  // use pretty print and indent of 2
  // std::cout << std::setw(2) << j << std::endl;
}

void System::readRestartFile()
{
  nlohmann::json j;

  simulationBox = j["simulationBox"];
  atomPositions = j["atomPositions"];
  moleculePositions = j["moleculePositions"];
}

std::string System::repr() const { return std::string("system test"); }

// std::vector<double> createVDWGrid(const ForceField &forcefield, const std::vector<Framework> & frameworkComponents,
// size_t typeA);

void System::createInterpolationGrids()
{
  double percent = 100.0 / static_cast<double>(128 * gridPseudoAtomIndices.size());

  for (size_t pseudo_atom_index : gridPseudoAtomIndices)
  {
    InterpolationEnergyGrid grid(int3(128, 128, 128));

    std::mdspan<double, std::dextents<size_t, 4>, std::layout_left> data(grid.data.data(), 128, 128, 128, 8);

    double teller{};
    for (size_t i = 0; i <= 127; i++)
    {
      teller += 1.0;
      for (size_t j = 0; j <= 127; j++)
      {
        for (size_t k = 0; k <= 127; k++)
        {
          double3 s = double3(static_cast<double>(i) / static_cast<double>(127),
                              static_cast<double>(j) / static_cast<double>(127),
                              static_cast<double>(k) / static_cast<double>(127));

          double3 pos = simulationBox.cell * s;

          auto [energy, gradient, hessian, third_derivative] = Interactions::calculateThirdDerivativeAtPositionVDW(
              forceField, simulationBox, pos, pseudo_atom_index, spanOfFrameworkAtoms());

          double energy_fractional = energy;
          double3 gradient_fractional = simulationBox.cell.transpose() * gradient;
          double3x3 hessian_fractional = simulationBox.cell.transpose() * hessian * simulationBox.cell;
          double third_derivative_fractional_xyz =
              simulationBox.cell.ax * simulationBox.cell.bx * simulationBox.cell.cx * third_derivative.m111 +
              simulationBox.cell.ax * simulationBox.cell.bx * simulationBox.cell.cy * third_derivative.m112 +
              simulationBox.cell.ax * simulationBox.cell.bx * simulationBox.cell.cz * third_derivative.m113 +
              simulationBox.cell.ax * simulationBox.cell.by * simulationBox.cell.cx * third_derivative.m121 +
              simulationBox.cell.ax * simulationBox.cell.by * simulationBox.cell.cy * third_derivative.m122 +
              simulationBox.cell.ax * simulationBox.cell.by * simulationBox.cell.cz * third_derivative.m132;

          data[i, j, k, 0] = energy_fractional;
          data[i, j, k, 1] = gradient_fractional.x;
          data[i, j, k, 2] = gradient_fractional.y;
          data[i, j, k, 3] = gradient_fractional.z;
          data[i, j, k, 4] = hessian_fractional.ay;
          data[i, j, k, 5] = hessian_fractional.az;
          data[i, j, k, 6] = hessian_fractional.bz;
          data[i, j, k, 7] = third_derivative_fractional_xyz;
        }
      }
      if (static_cast<size_t>(teller * percent) > static_cast<size_t>((teller - 1.0) * percent))
      {
        std::cout << "Percentage finished: " << static_cast<size_t>(teller * percent) << std::endl;
        ;
      }
    }

    interpolationGrids[pseudo_atom_index] = grid;
  }

  // for(const Framework &framework: frameworks)
  //{
  // }
}
