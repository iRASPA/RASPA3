module;

#ifdef USE_LEGACY_HEADERS
#include <iostream>
#include <ostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <array>
#include <map>
#include <string>
#include <span>
#include <optional>
#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <exception>
#include <iterator>
#include <chrono>
#include <cstddef>
#include <type_traits>
#if defined(_WIN32)
#include <cassert.h>
#endif
#include <exception>
#include <source_location>
#include <complex>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

#if !defined(_WIN32)
#include <assert.h>
#endif

module component;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <ostream>;
import <sstream>;
import <algorithm>;
import <vector>;
import <array>;
import <map>;
import <string>;
import <span>;
import <optional>;
import <filesystem>;
import <fstream>;
import <cstdlib>;
import <exception>;
import <iterator>;
import <chrono>;
import <cstddef>;
import <type_traits>;
#if defined(_WIN32)
import <cassert>;
#endif
import <exception>;
import <source_location>;
import <complex>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import int3;
import double3;
import double3x3;
import randomnumbers;
import archive;
import json;
import skposcarparser;
import characterset;
import stringutils;
import skparser;
import skposcarparser;
import skstructure;
import skasymmetricatom;
import skatomcopy;
import skcell;
import skspacegroup;
import forcefield;
import atom;
import property_lambda_probability_histogram;
import property_widom;
import multi_site_isotherm;
import simulationbox;
import cif_reader;
import move_statistics;
import bond_potential;
import mc_moves_probabilities_particles;
import mc_moves_statistics_particles;
import mc_moves_cputime;
import mc_moves_count;


// default constructor, needed for binary restart-file
Component::Component()
{
}

// create Component in 'inputreader.cpp'
Component::Component(Component::Type type, size_t currentComponent, const ForceField &forceField, const std::string &componentName,
                     std::optional<const std::string> fileName,
                     size_t numberOfBlocks, size_t numberOfLambdaBins,
                     const MCMoveProbabilitiesParticles &particleProbalities) noexcept(false) :
                     type(type), 
                     componentId(currentComponent), 
                     name(componentName),
                     filenameData(fileName),
                     lambdaGC(numberOfBlocks, numberOfLambdaBins),
                     lambdaGibbs(numberOfBlocks, numberOfLambdaBins),
                     mc_moves_probabilities(particleProbalities),
                     averageRosenbluthWeights(numberOfBlocks)
{
  if (filenameData.has_value())
  {
    readComponent(forceField, filenameData.value());
  }
}

// create programmatically an 'adsorbate' component
Component::Component(size_t componentId, const ForceField &forceField, std::string componentName,
                     double T_c, double P_c, double w, std::vector<Atom> definedAtoms,
                     size_t numberOfBlocks, size_t numberOfLambdaBins,
                     const MCMoveProbabilitiesParticles &particleProbalities) noexcept(false) :
    type(Type::Adsorbate),
    componentId(componentId),
    name(componentName),
    criticalTemperature(T_c),
    criticalPressure(P_c),
    acentricFactor(w),
    definedAtoms(definedAtoms),
    atoms(definedAtoms),
    lambdaGC(numberOfBlocks, numberOfLambdaBins),
    lambdaGibbs(numberOfBlocks, numberOfLambdaBins),
    mc_moves_probabilities(particleProbalities),
    averageRosenbluthWeights(numberOfBlocks)
{
  mass = 0.0;
  for (const Atom& atom : atoms)
  {
    size_t atomType = static_cast<size_t>(atom.type);
    mass += forceField.pseudoAtoms[atomType].mass;
  }
}



// read the component from the molecule-file
void Component::readComponent(const ForceField& forceField, const std::string& fileName)
{
  const std::string defaultMoleculeFileName = fileName + ".json";

  std::string moleculeFileName = defaultMoleculeFileName;
  if(!std::filesystem::exists(std::filesystem::path{moleculeFileName}))
  {
    const char* env_p = std::getenv("RASPA_DIR");
    if (env_p)
    {
      moleculeFileName = env_p + std::string("/") + defaultMoleculeFileName;
    }
  }

  if (!std::filesystem::exists(moleculeFileName)) 
  {
    throw std::runtime_error(std::format("[Component reader]: File '{}' not found\n", moleculeFileName));
  }

  std::filesystem::path moleculePathfile = std::filesystem::path(moleculeFileName);
  std::ifstream moleculeStream{ moleculePathfile };
  if (!moleculeStream) 
  {
    throw std::runtime_error(std::format("[Component reader] File '{}' exists, but error opening file\n", moleculeFileName));
  }

  nlohmann::basic_json<nlohmann::raspa_map> parsed_data{};

  try
  {
    parsed_data = nlohmann::json::parse(moleculeStream);
  }
  catch (nlohmann::json::parse_error& ex)
  {
    throw std::runtime_error(std::format("[Component reader]: Parse error of file {} at byte {}\n{}\n",
                                         std::format("{}.json", fileName), ex.byte, ex.what()));
  }

  try
  {
    criticalTemperature = parsed_data.value("CriticalTemperature", 0.0);
  }
  catch (nlohmann::json::exception& ex)
  {
    throw std::runtime_error(std::format("[Component reader]: item 'CriticalTemperature' listed as {} must be floating point number\n{}\n",
           parsed_data["CriticalTemperature"].dump(), ex.what()));
  }

  try
  {
    criticalPressure = parsed_data.value("CriticalPressure", 0.0);
  }
  catch (nlohmann::json::exception& ex)
  {
    throw std::runtime_error(std::format("[Component reader]: item 'CriticalPressure' listed as {} must be floating point number\n{}\n",
           parsed_data["CriticalPressure"].dump(), ex.what()));
  }

  try
  {
    acentricFactor = parsed_data.value("AcentricFactor", 0.0);
  }
  catch (nlohmann::json::exception& ex)
  {
    throw std::runtime_error(std::format("[Component reader]: item 'AcentricFactor' listed as {} must be floating point number\n{}\n",
           parsed_data["AcentricFactor"].dump(), ex.what()));
  }

  size_t numberOfPseudoAtoms = parsed_data["PseudoAtoms"].size();

  definedAtoms.clear();
  definedAtoms.reserve(numberOfPseudoAtoms);

  if(numberOfPseudoAtoms == 0)
  {
    throw std::runtime_error(std::format("[Component reader]: key 'PseudoAtoms' empty]\n"));
  }

  if(!parsed_data.contains("PseudoAtoms"))
  {
    throw std::runtime_error(std::format("[Component reader]: No pseudo-atoms found [keyword 'PseudoAtoms' missing]\n"));
  }

  mass = 0.0;
  for (auto& [_, item ] : parsed_data["PseudoAtoms"].items())
  {
    if(!item.is_array())
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
    }


    if(item.size() != 2)
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be an array with two elements, "
                "the pseudo-atom-name and an array with the x,y,z positions\n", item.dump()));
    }

    if(!item[0].is_string())
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be an string (the name of the pseudo-atom)\n", item[0].dump()));
    }
    std::string pseudoAtomName = item[0].get<std::string>();

    // find atom-type based on read 'atomTypeString'
    std::optional<size_t> index = forceField.findPseudoAtom(pseudoAtomName);
    if(!index.has_value())
    {
      throw std::runtime_error(std::format("[Component reader]: unknown pseudo-atom '{}', please lookup type in in 'pseudo_atoms.json'\n", pseudoAtomName));
    }
    size_t pseudoAtomType = index.value();

    if(!item[1].is_array())
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be an array (with the positions)\n", item[1].dump()));
    }

    if(item[1].size() != 3)
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be an array with three elements, "
                "the x,y,z positions\n", item[1].dump()));
    }

    std::vector<double> position{};
    try
    {
      position = item[1].is_array() ? item[1].get<std::vector<double>>() : std::vector<double>{};
    }
    catch (nlohmann::json::exception& ex)
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be array of three floating point numbers \n{}\n",
             item[1].dump(), ex.what()));
    }

    mass += forceField.pseudoAtoms[pseudoAtomType].mass;
    double charge = forceField.pseudoAtoms[pseudoAtomType].charge;
    double scaling = 1.0;

    definedAtoms.emplace_back(double3(position[0], position[1], position[2]) , charge, scaling, 0, static_cast<uint16_t>(pseudoAtomType),
                           static_cast<uint8_t>(componentId), 0);
  }

  atoms = definedAtoms;
}


std::string Component::printStatus(const ForceField& forceField) const
{
  std::ostringstream stream;

  std::print(stream, "Component {} [{}]\n\n", componentId, name);

  std::print(stream, "    Critical temperature:  {} [K]\n", criticalTemperature);
  std::print(stream, "    Critical pressure:     {} [Pa]\n", criticalPressure);
  std::print(stream, "    Acentric factor:       {} [-]\n\n", acentricFactor);

  std::print(stream, "    Mol-fraction:                 {} [-]\n", molFraction);
  std::print(stream << std::boolalpha, "    Swapable:                     {}\n\n", swapable);
  std::print(stream, "    Mass:                         {} [-]\n", mass);
  if(fugacityCoefficient.has_value())
  {
    std::print(stream, "    Fugacity coefficient:         {} [-]\n", fugacityCoefficient.value());
  }
  std::print(stream, "    Bulk fluid density:           {} [-]\n", bulkFluidDensity);
  std::print(stream, "    Compressibility:              {} [-]\n", compressibility);
  std::print(stream, "    Excess molecules:             {} [-]\n\n", amountOfExcessMolecules);

  std::print(stream, "    number Of Atoms:  {}\n", atoms.size());
  
  for (size_t i = 0; i < definedAtoms.size(); ++i)
  {
    size_t atomType = static_cast<size_t>(definedAtoms[i].type);
    std::string atomTypeString = forceField.pseudoAtoms[atomType].name;
    std::print(stream, "    {:3d}: {:6} position {:8.5f} {:8.5f} {:8.5f}, charge {:8.5f}\n", 
                   i, atomTypeString, definedAtoms[i].position.x, definedAtoms[i].position.y, 
                   definedAtoms[i].position.z, definedAtoms[i].charge);
  }
  std::print(stream, "    net-charge:  {}\n", netCharge);
  std::print(stream, "\n");
 
  const MCMoveProbabilitiesParticles &mc = mc_moves_probabilities; 
  std::print(stream, "    Translation-move probability:             {} [-]\n", mc.probabilityTranslationMove);
  std::print(stream, "    Random translation-move probability:      {} [-]\n", mc.probabilityRandomTranslationMove);
  std::print(stream, "    Rotation-move probability:                {} [-]\n", mc.probabilityRotationMove);
  std::print(stream, "    Random rotation-move probability:         {} [-]\n", mc.probabilityRandomRotationMove);
  std::print(stream, "    Volume-move probability:                  {} [-]\n", mc.probabilityVolumeMove);
  std::print(stream, "    Reinsertion (CBMC) probability:           {} [-]\n", mc.probabilityReinsertionMove_CBMC);
  std::print(stream, "    Identity-change (CBMC) probability:       {} [-]\n", mc.probabilityIdentityChangeMove_CBMC);
  std::print(stream, "    Swap-move (CBMC) probability:             {} [-]\n", mc.probabilitySwapMove_CBMC);
  std::print(stream, "    Swap-move (CFCMC) probability:            {} [-]\n", mc.probabilitySwapMove_CFCMC);
  std::print(stream, "    Swap-move (CFCMC/CBMC) probability:       {} [-]\n", mc.probabilitySwapMove_CFCMC_CBMC);
  std::print(stream, "    Gibbs Volume-move probability:            {} [-]\n", mc.probabilityGibbsVolumeMove);
  std::print(stream, "    Gibbs Swap-move (CBMC) probability:       {} [-]\n", mc.probabilityGibbsSwapMove_CBMC);
  std::print(stream, "    Gibbs Swap-move (CFCMC) probability:      {} [-]\n", mc.probabilityGibbsSwapMove_CFCMC);
  std::print(stream, "    Gibbs Swap-move (CFCMC/CBMC) probability: {} [-]\n", mc.probabilityGibbsSwapMove_CFCMC_CBMC);
  std::print(stream, "    Widom probability:                        {} [-]\n", mc.probabilityWidomMove);
  std::print(stream, "    Widom (CFCMC) probability:                {} [-]\n", mc.probabilityWidomMove_CFCMC);
  std::print(stream, "    Widom (CFCMC/CBMC) probability:           {} [-]\n", mc.probabilityWidomMove_CFCMC_CBMC);
  std::print(stream, "\n");

  std::print(stream, "    number of bonds: {}\n", bonds.size());
  for (size_t i = 0; i < bonds.size(); ++i)
  {
    std::print(stream, "        {}", bonds[i].print());
  }
  std::print(stream, "\n");

  return stream.str();
}

std::vector<double3> Component::randomlyRotatedPositionsAroundStartingBead(RandomNumber &random) const
{
  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  std::vector<double3> randomPositions{};
  std::transform(std::begin(atoms), std::end(atoms),
          std::back_inserter(randomPositions), [&](const Atom& atom) 
          {return randomRotationMatrix * (atom.position - atoms[startingBead].position); });
  return randomPositions;
}

std::vector<Atom> Component::recenteredCopy(double scaling, size_t moleculeId) const
{
  std::vector<Atom> new_atoms(atoms);

  for (size_t i = 0; i < atoms.size(); ++i)
  {
    new_atoms[i] = Atom(atoms[i].position - atoms[startingBead].position, 
                 atoms[i].charge, scaling, static_cast<uint32_t>(moleculeId), static_cast<uint16_t>(atoms[i].type), 
                 static_cast<uint8_t>(componentId), static_cast<uint8_t>(atoms[i].groupId));
  }

  return new_atoms;
}

std::vector<Atom> Component::copyAtoms(std::span<Atom> molecule, double scaling, size_t moleculeId) const
{
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());
  for (size_t i = 0; i != atoms.size(); ++i)
  {
    copied_atoms[i].setScaling(scaling);
    copied_atoms[i].position = molecule[i].position - molecule[startingBead].position;
    copied_atoms[i].moleculeId = static_cast<uint32_t>(moleculeId);
  }
  return copied_atoms;
}

std::vector<Atom> Component::copyAtomsRandomlyRotatedAt(RandomNumber &random, double3 position, 
                                   std::span<Atom> molecule, double scaling, size_t moleculeId) const
{
  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());
  for (size_t i = 0; i != atoms.size(); ++i)
  {
    copied_atoms[i].setScaling(scaling);
    copied_atoms[i].position = position + randomRotationMatrix * 
                                          (molecule[i].position - molecule[startingBead].position);
    copied_atoms[i].moleculeId = static_cast<uint32_t>(moleculeId);
  }
  return copied_atoms;
}

std::vector<Atom> Component::copiedAtoms(std::span<Atom> molecule) const
{
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());
  for (size_t i = 0; i != atoms.size(); ++i)
  {
    copied_atoms[i].position = molecule[i].position - molecule[startingBead].position;
  }
  return copied_atoms;
}


std::string Component::printBreakthroughStatus() const
{
  std::ostringstream stream;

  std::print(stream, "Component {} [{}]\n", componentId, name);
  if(isCarrierGas)
  {
    std::print(stream, "    carrier-gas\n");

    std::print(stream, "{}", isotherm.print());
  }
  std::print(stream, "    mol-fraction in the gas:   {} [-]\n", molFraction);
  if(!isCarrierGas)
  {
    std::print(stream, "    mass-transfer coefficient: {} [1/s]\n", massTransferCoefficient);
    std::print(stream, "    diffusion coefficient:     {} [m^2/s]\n", axialDispersionCoefficient);

    std::print(stream, "{}", isotherm.print());
  }

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Component &c)
{
  archive << c.versionNumber;

  archive << c.type;
  archive << c.growType;

  archive << c.componentId;
  archive << c.name;
  archive << c.filenameData;
  archive << c.filename;

  archive << c.rigid;

  archive << c.criticalTemperature;
  archive << c.criticalPressure;
  archive << c.acentricFactor;
  archive << c.molFraction;
  archive << c.swapable;
  archive << c.partialPressure;

  archive << c.mass;
  archive << c.fugacityCoefficient;
  archive << c.amountOfExcessMolecules;
  archive << c.bulkFluidDensity;
  archive << c.compressibility;

  archive << c.idealGasRosenbluthWeight;
  archive << c.idealGasEnergy;

  archive << c.netCharge;
  archive << c.startingBead;
  archive << c.definedAtoms;
  archive << c.atoms;

  archive << c.initialNumberOfMolecules;

  archive << c.lambdaGC;
  archive << c.lambdaGibbs;
  archive << c.hasFractionalMolecule;

  archive << c.chiralCenters;
  archive << c.bonds;
  archive << c.connectivityTable;
  //std::vector<std::pair<size_t, size_t>> bondDipoles{};
  //std::vector<std::tuple<size_t, size_t, size_t>> bends{};
  //std::vector<std::pair<size_t, size_t>>  UreyBradley{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> inversionBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> Torsion{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> ImproperTorsions{};
  //std::vector<std::tuple<size_t, size_t, size_t>> bondBonds{};
  //std::vector<std::tuple<size_t, size_t, size_t>> stretchBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> stretchTorsions{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendTorsions{};
  //std::vector<std::pair<size_t, size_t>> intraVDW{};
  //std::vector<std::pair<size_t, size_t>> intraCoulomb{};
  //std::vector<std::pair<size_t, size_t>> excludedIntraCoulomb{};
  //std::vector<std::pair<size_t, std::vector<size_t>>> configMoves{};

  archive << c.mc_moves_probabilities;
  archive << c.mc_moves_statistics;
  archive << c.mc_moves_cputime;
  archive << c.mc_moves_count;

  archive << c.averageRosenbluthWeights;

  //MultiSiteIsotherm isotherm{};      // isotherm information
  archive << c.massTransferCoefficient;
  archive << c.axialDispersionCoefficient;
  archive << c.isCarrierGas;

  archive << c.columnPressure;
  archive << c.columnLoading;
  archive << c.columnError;

  archive << c.lnPartitionFunction;

  archive << c.pressureScale;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Component &c)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > c.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Component' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> c.type;
  archive >> c.growType;

  archive >> c.componentId;
  archive >> c.name;
  archive >> c.filenameData;
  archive >> c.filename;

  archive >> c.rigid;

  archive >> c.criticalTemperature;
  archive >> c.criticalPressure;
  archive >> c.acentricFactor;
  archive >> c.molFraction;
  archive >> c.swapable;
  archive >> c.partialPressure;

  archive >> c.mass;
  archive >> c.fugacityCoefficient;
  archive >> c.amountOfExcessMolecules;
  archive >> c.bulkFluidDensity;
  archive >> c.compressibility;

  archive >> c.idealGasRosenbluthWeight;
  archive >> c.idealGasEnergy;

  archive >> c.netCharge;
  archive >> c.startingBead;
  archive >> c.definedAtoms;
  archive >> c.atoms;

  archive >> c.initialNumberOfMolecules;

  archive >> c.lambdaGC;
  archive >> c.lambdaGibbs;
  archive >> c.hasFractionalMolecule;

  archive >> c.chiralCenters;
  archive >> c.bonds;
  archive >> c.connectivityTable;
  //std::vector<std::pair<size_t, size_t>> bondDipoles{};
  //std::vector<std::tuple<size_t, size_t, size_t>> bends{};
  //std::vector<std::pair<size_t, size_t>>  UreyBradley{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> inversionBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> Torsion{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> ImproperTorsions{};
  //std::vector<std::tuple<size_t, size_t, size_t>> bondBonds{};
  //std::vector<std::tuple<size_t, size_t, size_t>> stretchBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendBends{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> stretchTorsions{};
  //std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendTorsions{};
  //std::vector<std::pair<size_t, size_t>> intraVDW{};
  //std::vector<std::pair<size_t, size_t>> intraCoulomb{};
  //std::vector<std::pair<size_t, size_t>> excludedIntraCoulomb{};
  //std::vector<std::pair<size_t, std::vector<size_t>>> configMoves{};

  archive >> c.mc_moves_probabilities;
  archive >> c.mc_moves_statistics;
  archive >> c.mc_moves_cputime;
  archive >> c.mc_moves_count;

  archive >> c.averageRosenbluthWeights;

  //MultiSiteIsotherm isotherm{};      // isotherm information
  archive >> c.massTransferCoefficient;
  archive >> c.axialDispersionCoefficient;
  archive >> c.isCarrierGas;

  archive >> c.columnPressure;
  archive >> c.columnLoading;
  archive >> c.columnError;

  archive >> c.lnPartitionFunction;

  archive >> c.pressureScale;

  return archive;
}

std::string Component::repr() const
{
  return std::string("Component test");
}
