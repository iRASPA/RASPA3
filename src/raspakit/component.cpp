module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <optional>
#include <ostream>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#if defined(_WIN32)
#include <cassert.h>
#endif
#include <complex>
#include <exception>
#include <source_location>
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
import <utility>;
import <functional>;
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
import simd_quatd;
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
import hdf5;

// default constructor, needed for binary restart-file
Component::Component() {}

// create Component in 'inputreader.cpp'
Component::Component(Component::Type type, size_t currentComponent, const ForceField &forceField,
                     const std::string &componentName, std::optional<const std::string> fileName, size_t numberOfBlocks,
                     size_t numberOfLambdaBins, const MCMoveProbabilitiesParticles &particleProbalities) noexcept(false)
    : type(type),
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
Component::Component(size_t componentId, const ForceField &forceField, std::string componentName, double T_c,
                     double P_c, double w, std::vector<Atom> atomList, size_t numberOfBlocks, size_t numberOfLambdaBins,
                     const MCMoveProbabilitiesParticles &particleProbalities) noexcept(false)
    : type(Type::Adsorbate),
      componentId(componentId),
      name(componentName),
      criticalTemperature(T_c),
      criticalPressure(P_c),
      acentricFactor(w),
      lambdaGC(numberOfBlocks, numberOfLambdaBins),
      lambdaGibbs(numberOfBlocks, numberOfLambdaBins),
      mc_moves_probabilities(particleProbalities),
      averageRosenbluthWeights(numberOfBlocks)
{
  totalMass = 0.0;
  for (const Atom &atom : atomList)
  {
    size_t atomType = static_cast<size_t>(atom.type);
    double mass = forceField.pseudoAtoms[atomType].mass;
    totalMass += mass;
    definedAtoms.push_back({atom, mass});
  }

  computeRigidProperties();
}

// read the component from the molecule-file
void Component::readComponent(const ForceField &forceField, const std::string &fileName)
{
  const std::string defaultMoleculeFileName = fileName + ".json";

  std::string moleculeFileName = defaultMoleculeFileName;
  if (!std::filesystem::exists(std::filesystem::path{moleculeFileName}))
  {
    const char *env_p = std::getenv("RASPA_DIR");
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
  std::ifstream moleculeStream{moleculePathfile};
  if (!moleculeStream)
  {
    throw std::runtime_error(
        std::format("[Component reader] File '{}' exists, but error opening file\n", moleculeFileName));
  }

  nlohmann::basic_json<nlohmann::raspa_map> parsed_data{};

  try
  {
    parsed_data = nlohmann::json::parse(moleculeStream);
  }
  catch (nlohmann::json::parse_error &ex)
  {
    throw std::runtime_error(std::format("[Component reader]: Parse error of file {} at byte {}\n{}\n",
                                         std::format("{}.json", fileName), ex.byte, ex.what()));
  }

  try
  {
    criticalTemperature = parsed_data.value("CriticalTemperature", 0.0);
  }
  catch (nlohmann::json::exception &ex)
  {
    throw std::runtime_error(
        std::format("[Component reader]: item 'CriticalTemperature' listed as {} must be floating point number\n{}\n",
                    parsed_data["CriticalTemperature"].dump(), ex.what()));
  }

  try
  {
    criticalPressure = parsed_data.value("CriticalPressure", 0.0);
  }
  catch (nlohmann::json::exception &ex)
  {
    throw std::runtime_error(
        std::format("[Component reader]: item 'CriticalPressure' listed as {} must be floating point number\n{}\n",
                    parsed_data["CriticalPressure"].dump(), ex.what()));
  }

  try
  {
    acentricFactor = parsed_data.value("AcentricFactor", 0.0);
  }
  catch (nlohmann::json::exception &ex)
  {
    throw std::runtime_error(
        std::format("[Component reader]: item 'AcentricFactor' listed as {} must be floating point number\n{}\n",
                    parsed_data["AcentricFactor"].dump(), ex.what()));
  }

  size_t jsonNumberOfPseudoAtoms = parsed_data["PseudoAtoms"].size();

  definedAtoms.clear();
  definedAtoms.reserve(jsonNumberOfPseudoAtoms);

  if (jsonNumberOfPseudoAtoms == 0)
  {
    throw std::runtime_error(std::format("[Component reader]: key 'PseudoAtoms' empty]\n"));
  }

  if (!parsed_data.contains("PseudoAtoms"))
  {
    throw std::runtime_error(
        std::format("[Component reader]: No pseudo-atoms found [keyword 'PseudoAtoms' missing]\n"));
  }

  for (auto &[_, item] : parsed_data["PseudoAtoms"].items())
  {
    if (!item.is_array())
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
    }

    if (item.size() != 2)
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an array with two elements, "
                      "the pseudo-atom-name and an array with the x,y,z positions\n",
                      item.dump()));
    }

    if (!item[0].is_string())
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an string (the name of the pseudo-atom)\n", item[0].dump()));
    }
    std::string pseudoAtomName = item[0].get<std::string>();

    // find atom-type based on read 'atomTypeString'
    std::optional<size_t> index = forceField.findPseudoAtom(pseudoAtomName);
    if (!index.has_value())
    {
      throw std::runtime_error(
          std::format("[Component reader]: unknown pseudo-atom '{}', please lookup type in in 'pseudo_atoms.json'\n",
                      pseudoAtomName));
    }
    size_t pseudoAtomType = index.value();

    if (!item[1].is_array())
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an array (with the positions)\n", item[1].dump()));
    }

    if (item[1].size() != 3)
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an array with three elements, "
                      "the x,y,z positions\n",
                      item[1].dump()));
    }

    std::vector<double> position{};
    try
    {
      position = item[1].is_array() ? item[1].get<std::vector<double>>() : std::vector<double>{};
    }
    catch (nlohmann::json::exception &ex)
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be array of three floating point numbers \n{}\n",
                      item[1].dump(), ex.what()));
    }

    double mass = forceField.pseudoAtoms[pseudoAtomType].mass;
    double charge = forceField.pseudoAtoms[pseudoAtomType].charge;
    double scaling = 1.0;

    definedAtoms.push_back({Atom(double3(position[0], position[1], position[2]), charge, scaling, 0,
                                 static_cast<uint16_t>(pseudoAtomType), static_cast<uint8_t>(componentId), 0),
                            mass});
  }

  totalMass = 0.0;
  for (auto [_, mass] : definedAtoms)
  {
    totalMass += mass;
  }

  computeRigidProperties();
}

void Component::computeRigidProperties()
{
  double3 com{};
  double total_mass{};
  for (auto &[atom, mass] : definedAtoms)
  {
    com += mass * atom.position;
    total_mass += mass;
  }
  com = com / total_mass;

  double3x3 inertiaTensor{};
  for (auto &[atom, mass] : definedAtoms)
  {
    double3 dr = atom.position - com;
    inertiaTensor.ax += mass * dr.x * dr.x;
    inertiaTensor.bx += mass * dr.y * dr.x;
    inertiaTensor.cx += mass * dr.z * dr.x;
    inertiaTensor.ay += mass * dr.x * dr.y;
    inertiaTensor.by += mass * dr.y * dr.y;
    inertiaTensor.cy += mass * dr.z * dr.y;
    inertiaTensor.az += mass * dr.x * dr.z;
    inertiaTensor.bz += mass * dr.y * dr.z;
    inertiaTensor.cz += mass * dr.z * dr.z;
  }

  // the local body frame is taken to be that in which the rotational inertia tensor is diagonal
  double3 eigenvalues{};
  double3x3 eigenvectors{};
  inertiaTensor.EigenSystemSymmetric(eigenvalues, eigenvectors);

  inertiaVector = double3{0.0, 0.0, 0.0};
  atoms.clear();
  for (auto [atom, mass] : definedAtoms)
  {
    double3 dr = atom.position - com;
    double3 pos = eigenvectors.transpose() * dr;
    if (std::abs(pos.x) < 1e-8) pos.x = 0.0;
    if (std::abs(pos.y) < 1e-8) pos.y = 0.0;
    if (std::abs(pos.z) < 1e-8) pos.z = 0.0;

    inertiaVector.x += mass * (pos.y * pos.y + pos.z * pos.z);
    inertiaVector.y += mass * (pos.x * pos.x + pos.z * pos.z);
    inertiaVector.z += mass * (pos.x * pos.x + pos.y * pos.y);

    // correct the position
    atom.position = pos;
    atoms.push_back(atom);
  }

  // set axis system: Ixx >= Iyy >= Izz
  double rot_xyz = std::max({inertiaVector.x, inertiaVector.y, inertiaVector.z});
  if (rot_xyz >= inertiaVector.x)
  {
    if (inertiaVector.y >= rot_xyz)
    {
      for (auto &atom : atoms)
      {
        double temp = atom.position.x;
        atom.position.x = atom.position.y;
        atom.position.y = -temp;
      }
      inertiaVector.y = inertiaVector.x;
      inertiaVector.x = rot_xyz;
    }
    else if (inertiaVector.z >= rot_xyz)
    {
      for (auto &atom : atoms)
      {
        double temp = atom.position.x;
        atom.position.x = atom.position.z;
        atom.position.z = -temp;
      }
      inertiaVector.z = inertiaVector.x;
      inertiaVector.x = rot_xyz;
    }
  }
  if (inertiaVector.z > inertiaVector.y)
  {
    for (auto &atom : atoms)
    {
      double temp = atom.position.y;
      atom.position.y = atom.position.z;
      atom.position.z = -temp;
    }
    double temp = inertiaVector.z;
    inertiaVector.z = inertiaVector.y;
    inertiaVector.y = temp;
  }

  double rotlim = std::max(1.0e-2, inertiaVector.x + inertiaVector.y + inertiaVector.z) * 1.0e-5;

  inverseInertiaVector.x = (inertiaVector.x < rotlim) ? 0.0 : 1.0 / inertiaVector.x;
  inverseInertiaVector.y = (inertiaVector.y < rotlim) ? 0.0 : 1.0 / inertiaVector.y;
  inverseInertiaVector.z = (inertiaVector.z < rotlim) ? 0.0 : 1.0 / inertiaVector.z;

  double rotall = inertiaVector.x + inertiaVector.y + inertiaVector.z > 1.0e-5
                      ? inertiaVector.x + inertiaVector.y + inertiaVector.z
                      : 1.0;

  size_t index = 0;
  if (inertiaVector.x / rotall < 1.0e-5) ++index;
  if (inertiaVector.y / rotall < 1.0e-5) ++index;
  if (inertiaVector.z / rotall < 1.0e-5) ++index;

  translationalDegreesOfFreedom = 3;
  rotationalDegreesOfFreedom = 3;
  if (inertiaVector.x / rotall < 1.0e-5) --rotationalDegreesOfFreedom;
  if (inertiaVector.y / rotall < 1.0e-5) --rotationalDegreesOfFreedom;
  if (inertiaVector.z / rotall < 1.0e-5) --rotationalDegreesOfFreedom;

  shapeType = Component::Shape{index};
}

std::vector<Atom> Component::rotatePositions(const simd_quatd &q) const
{
  double3x3 rotationMatrix = double3x3::buildRotationMatrixInverse(q);
  std::vector<Atom> rotatedAtoms{};
  for (size_t i = 0; i < atoms.size(); ++i)
  {
    Atom a = atoms[i];
    a.position = rotationMatrix * atoms[i].position;
    rotatedAtoms.push_back(a);
  }
  return rotatedAtoms;
}

double3 Component::computeCenterOfMass(std::vector<Atom> atom_list) const
{
  std::vector<std::pair<Atom, double>> a{definedAtoms};
  for (size_t i = 0; i != a.size(); ++i)
  {
    a[i].first.position = atom_list[i].position;
  }
  double3 com{};
  double total_mass{};
  for (const auto &[atom, mass] : a)
  {
    com += mass * atom.position;
    total_mass += mass;
  }
  return com / total_mass;
}

std::string Component::printStatus(const ForceField &forceField) const
{
  std::ostringstream stream;

  std::print(stream, "Component {} [{}]\n\n", componentId, name);

  std::print(stream, "    Critical temperature:  {} [K]\n", criticalTemperature);
  std::print(stream, "    Critical pressure:     {} [Pa]\n", criticalPressure);
  std::print(stream, "    Acentric factor:       {} [-]\n\n", acentricFactor);

  std::print(stream, "    Mol-fraction:                 {} [-]\n", molFraction);
  std::print(stream << std::boolalpha, "    Swapable:                     {}\n\n", swapable);
  std::print(stream, "    Mass:                         {} [-]\n", totalMass);
  if (fugacityCoefficient.has_value())
  {
    std::print(stream, "    Fugacity coefficient:         {} [-]\n", fugacityCoefficient.value());
  }
  std::print(stream, "    Bulk fluid density:           {} [-]\n", bulkFluidDensity);
  std::print(stream, "    Compressibility:              {} [-]\n", compressibility);
  std::print(stream, "    Excess molecules:             {} [-]\n\n", amountOfExcessMolecules);

  std::print(stream, "    Number Of Atoms:    {}\n", atoms.size());
  std::print(stream, "    CBMC starting bead: {}\n", startingBead);
  for (size_t i = 0; i != atoms.size(); ++i)
  {
    size_t atomType = static_cast<size_t>(atoms[i].type);
    std::string atomTypeString = forceField.pseudoAtoms[atomType].name;
    std::print(stream, "    {:3d}: {:6} position {:8.5f} {:8.5f} {:8.5f}, charge {:8.5f}\n", i, atomTypeString,
               atoms[i].position.x, atoms[i].position.y, atoms[i].position.z, atoms[i].charge);
  }
  std::print(stream, "    Diagonalized inertia-vector:      {:10.8f} {:10.8f} {:10.8f}\n", inertiaVector.x,
             inertiaVector.y, inertiaVector.z);
  std::print(stream, "    Translational degrees of freedom: {}\n", translationalDegreesOfFreedom);
  std::print(stream, "    Rotational degrees of freedom:    {}\n", rotationalDegreesOfFreedom);
  switch (shapeType)
  {
    case Shape::NonLinear:
      std::print(stream, "    Shape of the molecule is: 'Non-linear'\n");
      break;
    case Shape::Linear:
      std::print(stream, "    Shape of the molecule is: 'Linear'\n");
      break;
    case Shape::Point:
      std::print(stream, "    Shape of the molecule is: 'Point-particle'\n");
      break;
  }
  std::print(stream, "    Net-charge:      {:10.8f}\n", netCharge);
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
  std::print(stream, "    Parallel Tempering Swap:                  {} [-]\n", mc.probabilityParallelTemperingSwap);
  std::print(stream, "\n");

  std::print(stream, "    number of bonds: {}\n", bonds.size());
  for (size_t i = 0; i < bonds.size(); ++i)
  {
    std::print(stream, "        {}", bonds[i].print());
  }
  std::print(stream, "\n");
  return stream.str();
}

void Component::logStatus(HDF5Writer &hdf5, const ForceField &forceField) const
{
  std::string group = std::format("{}_{}", componentId, name);
  hdf5.createGroup(group);
  hdf5.writeMetaInfo(group, "Critical temperature [K]", criticalTemperature);
  hdf5.writeMetaInfo(group, "Critical pressure [Pa]", criticalPressure);
  hdf5.writeMetaInfo(group, "Acentric factor [-]", acentricFactor);
  hdf5.writeMetaInfo(group, "Mol fraction [-]", molFraction);
  hdf5.writeMetaInfo(group, "Swappable", std::format("{}", swapable));
  hdf5.writeMetaInfo(group, "Mass", totalMass);
  if (fugacityCoefficient.has_value())
  {
    hdf5.writeMetaInfo(group, "Fugacity coefficient [-]", fugacityCoefficient.value());
  }
  hdf5.writeMetaInfo(group, "Bulk fluid density [-]", bulkFluidDensity);
  hdf5.writeMetaInfo(group, "Compressibility [-]", compressibility);
  hdf5.writeMetaInfo(group, "Excess molecules [-]", amountOfExcessMolecules);
  hdf5.writeMetaInfo(group, "Number of Atoms", atoms.size());
  hdf5.writeMetaInfo(group, "CBMC starting bead", startingBead);

  std::vector<double> position(3 * atoms.size());
  std::vector<double> charge(atoms.size());
  std::vector<std::string> typenames(atoms.size());
  for (size_t i = 0; i < atoms.size(); i++)
  {
    size_t atomType = static_cast<size_t>(atoms[i].type);
    typenames[i] = forceField.pseudoAtoms[atomType].name;
    position[3 * i] = atoms[i].position.x;
    position[3 * i + 1] = atoms[i].position.y;
    position[3 * i + 2] = atoms[i].position.z;
    charge[i] = atoms[i].charge;
  }

  hdf5.createDataset<double>(group, "positions", {atoms.size(), 3}, {{"dimensions", "(numberOfAtoms, 3)"}});
  hdf5.writeVector(group, "positions", position);

  hdf5.createDataset<double>(group, "charges", {atoms.size()}, {{"dimensions", "(numberOfAtoms, )"}});
  hdf5.writeVector(group, "charges", charge);

  hdf5.createStringDataset(group, "types", {atoms.size()}, 8);
  hdf5.writeVector<std::string>(group, "types", typenames);

  hdf5.createDataset<double>(group, "Diagonalized inertia-vector", {3}, {{"dimensions", "(3, )"}});
  hdf5.writeVector(group, "Diagonalized inertia-vector",
                   std::vector<double>{inertiaVector.x, inertiaVector.y, inertiaVector.z});
  hdf5.writeMetaInfo(group, "Translational degrees of freedom", translationalDegreesOfFreedom);
  hdf5.writeMetaInfo(group, "Rotational degrees of freedom", rotationalDegreesOfFreedom);
  switch (shapeType)
  {
    case Shape::NonLinear:

      hdf5.writeMetaInfo(group, "Shape", "Non-linear");
      break;
    case Shape::Linear:
      hdf5.writeMetaInfo(group, "Shape", "Linear");
      break;
    case Shape::Point:
      hdf5.writeMetaInfo(group, "Shape", "Point");
      break;
  }

  const MCMoveProbabilitiesParticles &mc = mc_moves_probabilities;
  hdf5.writeMetaInfo(group, "Translation-move probability [-]", mc.probabilityTranslationMove);
  hdf5.writeMetaInfo(group, "Random translation-move probability [-]", mc.probabilityRandomTranslationMove);
  hdf5.writeMetaInfo(group, "Rotation-move probability [-]", mc.probabilityRotationMove);
  hdf5.writeMetaInfo(group, "Random rotation-move probability [-]", mc.probabilityRandomRotationMove);
  hdf5.writeMetaInfo(group, "Volume-move probability [-]", mc.probabilityVolumeMove);
  hdf5.writeMetaInfo(group, "Reinsertion (CBMC) probability [-]", mc.probabilityReinsertionMove_CBMC);
  hdf5.writeMetaInfo(group, "Identity-change (CBMC) probability [-]", mc.probabilityIdentityChangeMove_CBMC);
  hdf5.writeMetaInfo(group, "Swap-move (CBMC) probability [-]", mc.probabilitySwapMove_CBMC);
  hdf5.writeMetaInfo(group, "Swap-move (CFCMC) probability [-]", mc.probabilitySwapMove_CFCMC);
  hdf5.writeMetaInfo(group, "Swap-move (CFCMC/CBMC) probability [-]", mc.probabilitySwapMove_CFCMC_CBMC);
  hdf5.writeMetaInfo(group, "Gibbs Volume-move probability [-]", mc.probabilityGibbsVolumeMove);
  hdf5.writeMetaInfo(group, "Gibbs Swap-move (CBMC) probability [-]", mc.probabilityGibbsSwapMove_CBMC);
  hdf5.writeMetaInfo(group, "Gibbs Swap-move (CFCMC) probability [-]", mc.probabilityGibbsSwapMove_CFCMC);
  hdf5.writeMetaInfo(group, "Gibbs Swap-move (CFCMC/CBMC) probability [-]", mc.probabilityGibbsSwapMove_CFCMC_CBMC);
  hdf5.writeMetaInfo(group, "Widom probability [-]", mc.probabilityWidomMove);
  hdf5.writeMetaInfo(group, "Widom (CFCMC) probability [-]", mc.probabilityWidomMove_CFCMC);
  hdf5.writeMetaInfo(group, "Widom (CFCMC/CBMC) probability [-]", mc.probabilityWidomMove_CFCMC_CBMC);
  hdf5.writeMetaInfo(group, "Parallel Tempering Swap [-]", mc.probabilityParallelTemperingSwap);

  hdf5.writeMetaInfo(group, "Number of Bonds", bonds.size());
  std::vector<std::string> bondTypes(bonds.size());
  for (size_t i = 0; i < bonds.size(); ++i)
  {
    bondTypes[i] = bonds[i].print();
  }
  hdf5.createStringDataset(group, "bondtypes", {bonds.size()}, 128);
  hdf5.writeVector<std::string>(group, "bondtypes", bondTypes);
}

std::vector<Atom> Component::rotateRandomlyAroundCenterOfMass([[maybe_unused]] const simd_quatd &q)
{
  double3 com{};

  double total_mass{};
  for (auto &[atom, mass] : definedAtoms)
  {
    com += mass * atom.position;
    totalMass += mass;
  }
  com = com / total_mass;

  // TODO
  return {};
}

std::vector<double3> Component::randomlyRotatedPositionsAroundStartingBead(RandomNumber &random) const
{
  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  std::vector<double3> randomPositions{};
  std::transform(std::begin(atoms), std::end(atoms), std::back_inserter(randomPositions), [&](const Atom &atom)
                 { return randomRotationMatrix * (atom.position - atoms[startingBead].position); });
  return randomPositions;
}

std::vector<Atom> Component::recenteredCopy(double scaling, size_t moleculeId) const
{
  std::vector<Atom> new_atoms(atoms.size());

  for (size_t i = 0; i != atoms.size(); ++i)
  {
    new_atoms[i] = Atom(atoms[i].position - atoms[startingBead].position, atoms[i].charge, scaling,
                        static_cast<uint32_t>(moleculeId), static_cast<uint16_t>(atoms[i].type),
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
                                                        std::span<Atom> molecule, double scaling,
                                                        size_t moleculeId) const
{
  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());
  for (size_t i = 0; i != atoms.size(); ++i)
  {
    copied_atoms[i].setScaling(scaling);
    copied_atoms[i].position =
        position + randomRotationMatrix * (molecule[i].position - molecule[startingBead].position);
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

std::pair<Molecule, std::vector<Atom>> Component::equilibratedMoleculeRandomInBox(
    RandomNumber &random, const SimulationBox &simulationBox) const
{
  simd_quatd q = random.randomSimdQuatd();
  double3x3 M = double3x3::buildRotationMatrixInverse(q);
  double3 com = simulationBox.randomPosition(random);

  std::vector<Atom> trial_atoms(atoms);
  for (size_t i = 0; i != atoms.size(); i++)
  {
    trial_atoms[i].position = com + M * atoms[i].position;
  }
  /*
  size_t startingBead = components[selectedComponent].startingBead;
  double3 center = components[selectedComponent].atoms[startingBead].position;
  std::vector<Atom> copied_atoms(components[selectedComponent].atoms.begin(),
  components[selectedComponent].atoms.end());

  double3x3 randomRotationMatrix = random.randomRotationMatrix();
  double3 position = simulationBox.randomPosition(random);

  for (size_t i = 0; i != copied_atoms.size(); ++i)
  {
    copied_atoms[i].setScaling(scaling);
    copied_atoms[i].position = position + randomRotationMatrix * (components[selectedComponent].atoms[i].position -
  center); copied_atoms[i].moleculeId = static_cast<uint32_t>(moleculeId);
  }
  */
  return {{com, q}, trial_atoms};
}

std::pair<Molecule, std::vector<Atom>> Component::translate(const Molecule &molecule, std::span<Atom> molecule_atoms,
                                                            double3 displacement) const
{
  std::vector<Atom> trialAtoms(molecule_atoms.begin(), molecule_atoms.end());
  Molecule trialMolecule = molecule;

  if (rigid)
  {
    simd_quatd q = trialMolecule.orientation;
    double3x3 M = double3x3::buildRotationMatrixInverse(q);

    trialMolecule.centerOfMassPosition += displacement;
    double3 com = trialMolecule.centerOfMassPosition;
    for (size_t i = 0; i != trialAtoms.size(); ++i)
    {
      trialAtoms[i].position = com + M * atoms[i].position;
    }
  }
  else
  {
    std::transform(molecule_atoms.begin(), molecule_atoms.end(), trialAtoms.begin(),
                   [&](Atom a)
                   {
                     a.position += displacement;
                     return a;
                   });
  }

  return {trialMolecule, trialAtoms};
}

std::pair<Molecule, std::vector<Atom>> Component::rotate(const Molecule &molecule, std::span<Atom> molecule_atoms,
                                                         simd_quatd rotation) const
{
  std::vector<Atom> trialAtoms(molecule_atoms.begin(), molecule_atoms.end());
  Molecule trialMolecule = molecule;

  if (rigid)
  {
    simd_quatd q = rotation * trialMolecule.orientation;
    double3x3 M = double3x3::buildRotationMatrixInverse(q);

    trialMolecule.orientation = q;
    double3 com = trialMolecule.centerOfMassPosition;
    for (size_t i = 0; i != trialAtoms.size(); ++i)
    {
      trialAtoms[i].position = com + M * atoms[i].position;
    }
  }
  else
  {
    double3x3 rotationMatrix = double3x3(rotation);
    std::transform(molecule_atoms.begin(), molecule_atoms.end(), trialAtoms.begin(),
                   [&](Atom a)
                   {
                     a.position = rotationMatrix * (a.position - molecule_atoms[startingBead].position) +
                                  molecule_atoms[startingBead].position;
                     return a;
                   });
  }

  return {trialMolecule, trialAtoms};
}

std::string Component::printBreakthroughStatus() const
{
  std::ostringstream stream;

  std::print(stream, "Component {} [{}]\n", componentId, name);
  if (isCarrierGas)
  {
    std::print(stream, "    carrier-gas\n");

    std::print(stream, "{}", isotherm.print());
  }
  std::print(stream, "    mol-fraction in the gas:   {} [-]\n", molFraction);
  if (!isCarrierGas)
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

  archive << c.totalMass;
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
  // std::vector<std::pair<size_t, size_t>> bondDipoles{};
  // std::vector<std::tuple<size_t, size_t, size_t>> bends{};
  // std::vector<std::pair<size_t, size_t>>  UreyBradley{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> inversionBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> Torsion{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> ImproperTorsions{};
  // std::vector<std::tuple<size_t, size_t, size_t>> bondBonds{};
  // std::vector<std::tuple<size_t, size_t, size_t>> stretchBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> stretchTorsions{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendTorsions{};
  // std::vector<std::pair<size_t, size_t>> intraVDW{};
  // std::vector<std::pair<size_t, size_t>> intraCoulomb{};
  // std::vector<std::pair<size_t, size_t>> excludedIntraCoulomb{};
  // std::vector<std::pair<size_t, std::vector<size_t>>> configMoves{};

  archive << c.mc_moves_probabilities;
  archive << c.mc_moves_statistics;
  archive << c.mc_moves_cputime;
  archive << c.mc_moves_count;

  archive << c.averageRosenbluthWeights;

  // MultiSiteIsotherm isotherm{};      // isotherm information
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
  if (versionNumber > c.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Component' at line {} in file {}\n", location.line(),
                                         location.file_name()));
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

  archive >> c.totalMass;
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
  // std::vector<std::pair<size_t, size_t>> bondDipoles{};
  // std::vector<std::tuple<size_t, size_t, size_t>> bends{};
  // std::vector<std::pair<size_t, size_t>>  UreyBradley{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> inversionBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> Torsion{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> ImproperTorsions{};
  // std::vector<std::tuple<size_t, size_t, size_t>> bondBonds{};
  // std::vector<std::tuple<size_t, size_t, size_t>> stretchBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendBends{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> stretchTorsions{};
  // std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendTorsions{};
  // std::vector<std::pair<size_t, size_t>> intraVDW{};
  // std::vector<std::pair<size_t, size_t>> intraCoulomb{};
  // std::vector<std::pair<size_t, size_t>> excludedIntraCoulomb{};
  // std::vector<std::pair<size_t, std::vector<size_t>>> configMoves{};

  archive >> c.mc_moves_probabilities;
  archive >> c.mc_moves_statistics;
  archive >> c.mc_moves_cputime;
  archive >> c.mc_moves_count;

  archive >> c.averageRosenbluthWeights;

  // MultiSiteIsotherm isotherm{};      // isotherm information
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

std::string Component::repr() const { return std::string("Component test"); }
