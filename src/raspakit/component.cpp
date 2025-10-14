module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <unordered_map>
#include <numeric>
#include <optional>
#include <ostream>
#include <print>
#include <source_location>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#endif

#if !defined(_WIN32)
#include <assert.h>
#endif

module component;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import simd_quatd;
import double3;
import double3x3;
import randomnumbers;
import archive;
import json;
import units;
import skposcarparser;
import characterset;
import stringutils;
import skparser;
import skposcarparser;
import skstructure;
import skatom;
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
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_statistics;
import mc_moves_cputime;
import bond_potential;
import urey_bradley_potential;
import bend_potential;
import inversion_bend_potential;
import out_of_plane_bend_potential;
import torsion_potential;
import bond_bond_potential;
import bond_bend_potential;
import bond_torsion_potential;
import bend_bend_potential;
import bend_torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import intra_molecular_potentials;
import vdwparameters;
import json;

// default constructor, needed for binary restart-file
Component::Component() {}

// create Component in 'inputreader.cpp'
Component::Component(Component::Type type, std::size_t currentComponent, const ForceField &forceField,
                     const std::string &componentName, std::optional<const std::string> fileName,
                     std::size_t numberOfBlocks, std::size_t numberOfLambdaBins,
                     const MCMoveProbabilities &particleProbabilities, std::optional<double> fugacityCoefficient,
                     bool thermodynamicIntegration) noexcept(false)
    : type(type),
      componentId(currentComponent),
      name(componentName),
      filenameData(fileName),
      fugacityCoefficient(fugacityCoefficient),
      lambdaGC(numberOfBlocks, numberOfLambdaBins),
      lambdaGibbs(numberOfBlocks, numberOfLambdaBins),
      mc_moves_probabilities(particleProbabilities),
      averageRosenbluthWeights(numberOfBlocks)
{
  if (filenameData.has_value())
  {
    readComponent(forceField, filenameData.value());
  }
  lambdaGC.computeDUdlambda = thermodynamicIntegration;

  cbmc_moves_statistics = std::vector<CBMCMoveStatistics>(definedAtoms.size());
}

// create programmatically an 'adsorbate' component
Component::Component(std::size_t componentId, const ForceField &forceField, std::string componentName, double T_c,
                     double P_c, double w, std::vector<Atom> atomList, const ConnectivityTable &connectivityTable,
                     const Potentials::IntraMolecularPotentials &intraMolecularPotentials, std::size_t numberOfBlocks,
                     std::size_t numberOfLambdaBins, const MCMoveProbabilities &particleProbabilities,
                     std::optional<double> fugacityCoefficient, bool thermodynamicIntegration) noexcept(false)
    : type(Type::Adsorbate),
      componentId(componentId),
      name(componentName),
      criticalTemperature(T_c),
      criticalPressure(P_c),
      acentricFactor(w),
      fugacityCoefficient(fugacityCoefficient),
      connectivityTable(connectivityTable),
      intraMolecularPotentials(intraMolecularPotentials),
      lambdaGC(numberOfBlocks, numberOfLambdaBins),
      lambdaGibbs(numberOfBlocks, numberOfLambdaBins),
      mc_moves_probabilities(particleProbabilities),
      cbmc_moves_statistics(atomList.size()),
      averageRosenbluthWeights(numberOfBlocks)
{
  totalMass = 0.0;
  netCharge = 0.0;
  for (const Atom &atom : atomList)
  {
    std::size_t atomType = static_cast<std::size_t>(atom.type);
    double mass = forceField.pseudoAtoms[atomType].mass;
    totalMass += mass;
    netCharge += atom.charge;
    definedAtoms.push_back({atom, mass});
  }

  computeRigidProperties();
  lambdaGC.computeDUdlambda = thermodynamicIntegration;

  if (!intraMolecularPotentials.bonds.empty())
  {
    rigid = false;
    growType = Component::GrowType::Flexible;
  }
}

// read the component from the molecule-file
void Component::readComponent(const ForceField &forceField, const std::string &fileName)
{
  const std::string defaultMoleculeFileName = addExtension(fileName, ".json");

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

  for (auto &[_, item] : parsed_data["BlockingPockets"].items())
  {
    if (!item.is_array())
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
    }

    if (item.size() != 4)
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an array with four elements, "
                      "an array with the x,y,z positions, and a radius\n",
                      item.dump()));
    }

    std::vector<double> data = item.is_array() ? item.get<std::vector<double>>() : std::vector<double>{};
    blockingPockets.push_back(double4(data[0], data[1], data[2], data[3]));
  }

  std::size_t jsonNumberOfPseudoAtoms = parsed_data["PseudoAtoms"].size();

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
    std::optional<std::size_t> index = forceField.findPseudoAtom(pseudoAtomName);
    if (!index.has_value())
    {
      throw std::runtime_error(
          std::format("[Component reader]: unknown pseudo-atom '{}', please lookup type in in 'force_field.json'\n",
                      pseudoAtomName));
    }
    std::size_t pseudoAtomType = index.value();

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

    definedAtoms.push_back(
        {Atom(double3(position[0], position[1], position[2]), charge, scaling, 0,
              static_cast<std::uint16_t>(pseudoAtomType), static_cast<std::uint8_t>(componentId), 0, 0),
         mass});
  }

  if (parsed_data.contains("StartingBead"))
  {
    if (!parsed_data["StartingBead"].is_number_integer())
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an integer\n", parsed_data["StartingBead"].dump()));
    }

    std::int64_t starting_bead = parsed_data["StartingBead"].get<std::int64_t>();
    if (starting_bead >= 0)
    {
      startingBead = static_cast<std::size_t>(starting_bead);
    }
  }

  totalMass = 0.0;
  netCharge = 0.0;
  for (auto [atom, mass] : definedAtoms)
  {
    totalMass += mass;
    netCharge += atom.charge;
  }

  if (parsed_data.contains("Type"))
  {
    if (!parsed_data["Type"].is_string())
    {
      throw std::runtime_error(std::format(
          "[Component reader]: item {} must be an string (the name of the pseudo-atom)\n", parsed_data["Type"].dump()));
    }
    std::string typeString = parsed_data["Type"].get<std::string>();
    if (caseInSensStringCompare(typeString, "Flexible"))
    {
      rigid = false;
      growType = GrowType::Flexible;
    }
    else if (caseInSensStringCompare(typeString, "Rigid"))
    {
      rigid = true;
      growType = GrowType::Rigid;
    }
  }

  computeRigidProperties();

  connectivityTable = readConnectivityTable(definedAtoms.size(), parsed_data);
  if (growType == GrowType::Flexible)
  {
    intraMolecularPotentials.bonds = readBondPotentials(forceField, parsed_data);
    intraMolecularPotentials.bends = readBendPotentials(forceField, parsed_data);
    intraMolecularPotentials.torsions = readTorsionPotentials(forceField, parsed_data);
    intraMolecularPotentials.vanDerWaals = readVanDerWaalsPotentials(forceField, parsed_data);
    intraMolecularPotentials.coulombs = readCoulombPotentials(forceField, parsed_data);

    partialReinsertionFixedAtoms = readPartialReinsertionFixedAtoms(parsed_data);
  }
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

  std::size_t index = 0;
  if (inertiaVector.x / rotall < 1.0e-5) ++index;
  if (inertiaVector.y / rotall < 1.0e-5) ++index;
  if (inertiaVector.z / rotall < 1.0e-5) ++index;

  if(rigid)
  {
    translationalDegreesOfFreedom = 3;
    rotationalDegreesOfFreedom = 3;
    if (inertiaVector.x / rotall < 1.0e-5) --rotationalDegreesOfFreedom;
    if (inertiaVector.y / rotall < 1.0e-5) --rotationalDegreesOfFreedom;
    if (inertiaVector.z / rotall < 1.0e-5) --rotationalDegreesOfFreedom;

    shapeType = Component::Shape{index};
  }
  else
  {
    translationalDegreesOfFreedom = 3 * definedAtoms.size();
    rotationalDegreesOfFreedom = 0;
    shapeType = Component::Shape::NonLinear;
  }
}

std::vector<Atom> Component::rotatePositions(const simd_quatd &q) const
{
  double3x3 rotationMatrix = double3x3::buildRotationMatrixInverse(q);
  std::vector<Atom> rotatedAtoms{};
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    Atom a = atoms[i];
    a.position = atoms[startingBead].position + rotationMatrix * (atoms[i].position - atoms[startingBead].position);
    rotatedAtoms.push_back(a);
  }
  return rotatedAtoms;
}

double3 Component::computeCenterOfMass(std::span<Atom> atom_list) const
{
  std::vector<std::pair<Atom, double>> a{definedAtoms};
  for (std::size_t i = 0; i != a.size(); ++i)
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

std::string Component::printStatus(const ForceField &forceField, double inputPressure) const
{
  std::ostringstream stream;

  std::print(stream, "Component {} [{}]\n\n", componentId, name);

  if (type == Component::Type::Adsorbate)
  {
    std::print(stream, "    Type: Adsorbate\n\n");
  }
  else if (type == Component::Type::Cation)
  {
    std::print(stream, "    Type: Cation\n\n");
  }

  std::print(stream, "    Critical temperature:  {} [{}]\n", criticalTemperature, Units::unitOfTemperatureString);
  std::print(stream, "    Critical pressure:     {} [{}]\n", criticalPressure, Units::unitOfPressureString);
  std::print(stream, "    Acentric factor:       {} [-]\n\n", acentricFactor);

  std::print(stream, "    Mol-fraction:                 {} [-]\n", molFraction);
  std::print(stream << std::boolalpha, "    Swappable:                   {}\n\n", swappable);
  std::print(stream, "    Mass:                         {} [-]\n", totalMass);
  std::print(stream, "    Fugacity:                     {} [Pa]\n",
             molFraction * fugacityCoefficient.value_or(1.0) * inputPressure);
  if (fugacityCoefficient.has_value())
  {
    std::print(stream, "    Fugacity coefficient:         {} [-]\n", fugacityCoefficient.value());
  }
  std::print(stream, "    Bulk fluid density:           {} [-]\n", bulkFluidDensity);
  std::print(stream, "    Compressibility:              {} [-]\n", compressibility);
  std::print(stream, "    Excess molecules:             {} [-]\n\n", amountOfExcessMolecules);

  std::print(stream, "    Number Of Atoms:              {}\n", atoms.size());
  std::print(stream, "    CBMC starting bead:           {}\n", startingBead);
  if (idealGasRosenbluthWeight.has_value())
  {
    std::print(stream, "    Ideal gas Rosenbluth weight:  {:10.8f}\n", idealGasRosenbluthWeight.value());
  }
  for (std::size_t i = 0; i != atoms.size(); ++i)
  {
    std::size_t atomType = static_cast<std::size_t>(atoms[i].type);
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
  std::print(stream, "    Net-charge:      {:12.8f} [e]\n", netCharge);
  std::print(stream, "\n");

  const std::unordered_map<MoveTypes, double> normalizedProbabilities = mc_moves_probabilities.normalizedMap();
  std::print(stream, "    Move probabilities:\n");
  for (auto &[moveType, probability] : normalizedProbabilities)
  {
    std::print(stream, "    {:<30} {:8.6f} [-]\n", moveNames[moveType] + ":", probability);
  }
  std::print(stream, "\n");

  std::print(stream, "    number of blocking-pockets: {}\n", blockingPockets.size());
  for (std::size_t i = 0; i < blockingPockets.size(); ++i)
  {
    std::print(stream, "        fractional s_x,s_y,s_z: {},{},{} radius: {}\n", blockingPockets[i].x,
               blockingPockets[i].y, blockingPockets[i].z, blockingPockets[i].w);
  }
  std::print(stream, "\n");

  switch (growType)
  {
    case GrowType::Rigid:
      std::print(stream, "    molecule is modelled as 'rigid'\n\n");
      break;
    case GrowType::Flexible:
      std::print(stream, "    molecule is modelled as 'flexible'\n\n");
      break;
    default:
      std::unreachable();
  }

  if (growType == GrowType::Flexible)
  {
    std::print(stream, "    connectivity:\n");
    std::print(stream, "{}\n", connectivityTable.print("        "));

    std::size_t number_of_total_bonds = connectivityTable.findAllBonds().size();
    std::print(stream, "    number of bond potentials: {} (out of {})\n", intraMolecularPotentials.bonds.size(),
               number_of_total_bonds);
    for (std::size_t i = 0; i < intraMolecularPotentials.bonds.size(); ++i)
    {
      std::print(stream, "        {}", intraMolecularPotentials.bonds[i].print());
    }
    std::print(stream, "\n");

    if (!intraMolecularPotentials.ureyBradleys.empty())
    {
      std::print(stream, "    number of Urey-Bradley potentials: {}\n", intraMolecularPotentials.ureyBradleys.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.ureyBradleys.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.ureyBradleys[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bends.empty())
    {
      std::size_t number_of_total_bends = connectivityTable.findAllBends().size();
      std::print(stream, "    number of bend potentials: {} (out of {})\n", intraMolecularPotentials.bends.size(),
                 number_of_total_bends);
      for (std::size_t i = 0; i < intraMolecularPotentials.bends.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bends[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.inversionBends.empty())
    {
      std::print(stream, "    number of inversion-bend potentials: {}\n",
                 intraMolecularPotentials.inversionBends.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.inversionBends.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.inversionBends[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.outOfPlaneBends.empty())
    {
      std::print(stream, "    number of out-of-plane bend potentials: {}\n",
                 intraMolecularPotentials.outOfPlaneBends.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.outOfPlaneBends.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.outOfPlaneBends[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.torsions.empty())
    {
      std::size_t number_of_total_torsions = connectivityTable.findAllTorsions().size();
      std::print(stream, "    number of torsion potentials: {} (out of {})\n", intraMolecularPotentials.torsions.size(),
                 number_of_total_torsions);
      for (std::size_t i = 0; i < intraMolecularPotentials.torsions.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.torsions[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.improperTorsions.empty())
    {
      std::print(stream, "    number of improper-torsion potentials: {}\n",
                 intraMolecularPotentials.improperTorsions.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.improperTorsions.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.improperTorsions[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bondBonds.empty())
    {
      std::print(stream, "    number of bond-bond potentials: {}\n", intraMolecularPotentials.bondBonds.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.bondBonds.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bondBonds[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bondBends.empty())
    {
      std::print(stream, "    number of bond-bend potentials: {}\n", intraMolecularPotentials.bondBends.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.bondBends.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bondBends[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bondTorsions.empty())
    {
      std::print(stream, "    number of bond-torsion potentials: {}\n", intraMolecularPotentials.bondTorsions.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.bondTorsions.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bondTorsions[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bendBends.empty())
    {
      std::print(stream, "    number of bend-bend potentials: {}\n", intraMolecularPotentials.bendBends.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.bendBends.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bendBends[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.bendTorsions.empty())
    {
      std::print(stream, "    number of bend-torsion potentials: {}\n", intraMolecularPotentials.bendTorsions.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.bendTorsions.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.bendTorsions[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.vanDerWaals.empty())
    {
      std::print(stream, "    number of Van der Waals potentials: {}\n", intraMolecularPotentials.vanDerWaals.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.vanDerWaals.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.vanDerWaals[i].print());
      }
      std::print(stream, "\n");
    }

    if (!intraMolecularPotentials.coulombs.empty())
    {
      std::print(stream, "    number of coulomb potentials: {}\n", intraMolecularPotentials.coulombs.size());
      for (std::size_t i = 0; i < intraMolecularPotentials.coulombs.size(); ++i)
      {
        std::print(stream, "        {}", intraMolecularPotentials.coulombs[i].print());
      }
      std::print(stream, "\n");
    }

    if (!partialReinsertionFixedAtoms.empty())
    {
      std::print(stream, "    number of partial-reinsertion moves: {}\n", partialReinsertionFixedAtoms.size());
      for (std::size_t i = 0; i < partialReinsertionFixedAtoms.size(); ++i)
      {
        std::print(stream, "        fixed atoms: {}\n", partialReinsertionFixedAtoms[i]);
      }
      std::print(stream, "\n");
    }
  }

  return stream.str();
}

nlohmann::json Component::jsonStatus() const
{
  nlohmann::json status;

  status["name"] = name;
  status["id"] = componentId;
  status["criticalTemperature"] = criticalTemperature;
  status["criticalPressure"] = criticalPressure;
  status["acentricFactor"] = acentricFactor;
  status["molFraction"] = molFraction;
  status["swappable"] = swappable;
  status["mass"] = totalMass;
  status["n_atoms"] = atoms.size();
  status["cbmcStartingBead"] = startingBead;
  status["diagonalizedInertiaVector"] = inertiaVector;
  status["translationalDOF"] = translationalDegreesOfFreedom;
  status["rotationalDOF"] = rotationalDegreesOfFreedom;
  status["netCharge"] = netCharge;

  if (fugacityCoefficient.has_value())
  {
    status["fugacityCoefficient"] = fugacityCoefficient.value();
  }

  switch (shapeType)
  {
    case Shape::NonLinear:
      status["shape"] = "nonlinear";
      break;
    case Shape::Linear:
      status["shape"] = "linear";
      break;
    case Shape::Point:
      status["shape"] = "point-particle";
      break;
  }

  nlohmann::json moves;
  std::unordered_map<MoveTypes, double> normalizedProbabilities = mc_moves_probabilities.normalizedMap();
  for (auto &[moveType, probability] : normalizedProbabilities)
  {
    moves[moveNames[moveType]] = probability;
  }
  status["moveProbabilities"] = moves;

  status["n_bonds"] = intraMolecularPotentials.bonds.size();
  std::vector<std::string> bondTypes(intraMolecularPotentials.bonds.size());
  for (std::size_t i = 0; i < intraMolecularPotentials.bonds.size(); ++i)
  {
    bondTypes[i] = intraMolecularPotentials.bonds[i].print();
  }
  status["bondTypes"] = bondTypes;
  return status;
}

std::vector<Atom> Component::copiedAtoms(std::span<Atom> molecule) const
{
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());
  for (std::size_t i = 0; i != atoms.size(); ++i)
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
  for (std::size_t i = 0; i != atoms.size(); i++)
  {
    trial_atoms[i].position = com + M * atoms[i].position;
  }

  return {{com, q, totalMass, componentId, atoms.size()}, trial_atoms};
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
    for (std::size_t i = 0; i != trialAtoms.size(); ++i)
    {
      trialAtoms[i].position = com + M * atoms[i].position;
    }
  }
  else
  {
    trialMolecule.centerOfMassPosition += displacement;
    for (Atom &atom : trialAtoms)
    {
      atom.position += displacement;
    }
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
    for (std::size_t i = 0; i != trialAtoms.size(); ++i)
    {
      trialAtoms[i].position = com + M * atoms[i].position;
    }
  }
  else
  {
    double3x3 rotationMatrix = double3x3(rotation);
    for (std::size_t i = 0; i < trialAtoms.size(); ++i)
    {
      trialAtoms[i].position = rotationMatrix * (molecule_atoms[i].position - molecule_atoms[startingBead].position) +
                               molecule_atoms[startingBead].position;
    }
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

ConnectivityTable Component::readConnectivityTable(std::size_t size,
                                                   const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  ConnectivityTable table(size);

  if (parsed_data.contains("Connectivity"))
  {
    for (auto &[_, item] : parsed_data["Connectivity"].items())
    {
      if (!item.is_array())
      {
        throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
      }

      if (item.size() != 2)
      {
        throw std::runtime_error(
            std::format("[Component reader]: item {} must be an array with two unsigned integers.\n", item.dump()));
      }

      try
      {
        std::vector<std::size_t> identifiers = item.get<std::vector<std::size_t>>();

        table[identifiers[0], identifiers[1]] = true;
        table[identifiers[1], identifiers[0]] = true;
      }
      catch (std::exception const &e)
      {
        throw std::runtime_error(std::format("Error in connectivities ({}): {}\n", item.dump(), e.what()));
      }
    }
  }
  return table;
}

std::vector<BondPotential> Component::readBondPotentials(const ForceField &forceField,
                                                         const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<BondPotential> bond_potentials{};

  std::vector<std::array<std::size_t, 2>> found_bonds = connectivityTable.findAllBonds();

  if (parsed_data.contains("Bonds"))
  {
    intraMolecularPotentials.bonds.reserve(found_bonds.size());

    // fill in each of the bonds
    for (const std::array<std::size_t, 2> found_bond : found_bonds)
    {
      // check which bond-definition matches on this definition
      // only one is allowed (no ambiguities)
      bool bond_identified{ false };

      for (auto &[_, item] : parsed_data["Bonds"].items())
      {
        if (!item.is_array())
        {
          throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
        }

        if (item.size() != 3)
        {
          throw std::runtime_error(
              std::format("[Component reader]: item {} must be an array with at least three elements, "
                          "(1) an array with the two identifiers of the bond, (2) the potential type, "
                          "and (3) an array containg the potential parameters\n",
                          item.dump()));
        }

        try
        {
          // parse [1]: the potential name
          std::string potential_name = item[1].get<std::string>();

          // parse [2]: the parameters
          std::vector<double> potential_parameters =
              item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};

          if (item[0][0].is_string() && item[0][1].is_string())
          {
            std::vector<std::string> identifiers = item[0].get<std::vector<std::string>>();

            std::optional<std::size_t> type_A = forceField.findPseudoAtom(identifiers[0]);
            std::optional<std::size_t> type_B = forceField.findPseudoAtom(identifiers[1]);

            if (!(type_A.has_value() && type_B.has_value()))
            {
              throw std::runtime_error(
                  std::format("[Component reader]: unknown pseudo-atom {} {}\n", identifiers[0], identifiers[1]));
            }

            if ((atoms[found_bond[0]].type == type_A.value() && atoms[found_bond[1]].type == type_B.value()) ||
                (atoms[found_bond[0]].type == type_B.value() && atoms[found_bond[1]].type == type_A.value()))
            {
              if (bond_identified)
              {
                throw std::runtime_error(std::format("Error in Bond-potential: ambiguity in define bond-potential (bond {}-{} matches on multiple definitions)\n",
                      found_bond[0], found_bond[1]));
              }

              bond_identified = true;

              BondPotential bond_potential =
                  BondPotential({found_bond[0], found_bond[1]}, BondPotential::definitionForString.at(potential_name),
                                potential_parameters);
              bond_potentials.push_back(bond_potential);
            }
          }
          else if (item[0][0].is_number_unsigned() && item[0][1].is_number_unsigned())
          {
            std::size_t A = item[0][0].get<std::size_t>();
            std::size_t B = item[0][1].get<std::size_t>();
            if ((found_bond[0] == A && found_bond[1] == B) ||
                (found_bond[0] == B && found_bond[1] == A))
            {
              if (bond_identified)
              {
                throw std::runtime_error(std::format("Error in Bond-potential: ambiguity in define bond-potential (bond {}-{} matches on multiple definitions)\n",
                      found_bond[0], found_bond[1]));
              }

              bond_identified = true;

              // parse using unsigned integer identifiers
              std::vector<std::size_t> bond =
                  item[0].is_array() ? item[0].get<std::vector<std::size_t>>() : std::vector<std::size_t>{};

              BondPotential bond_potential = BondPotential(
                  {bond[0], bond[1]}, BondPotential::definitionForString.at(potential_name), potential_parameters);
              bond_potentials.push_back(bond_potential);
            }
          }

        }
        catch (std::exception const &e)
        {
          throw std::runtime_error(std::format("Error in Bond-potential ({}): {}\n", item.dump(), e.what()));
        }


      }
    }
  }

  if (found_bonds.size() != bond_potentials.size())
  {
    throw std::runtime_error(std::format("Error in Bond-potential: not all bond-potentials defined ({} missing)\n",
                                         found_bonds.size() - bond_potentials.size()));
  }

  return bond_potentials;
}

std::vector<BendPotential> Component::readBendPotentials(const ForceField &forceField,
                                                         const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<BendPotential> bend_potentials{};

  std::vector<std::array<std::size_t, 3>> found_bends = connectivityTable.findAllBends();

  if (parsed_data.contains("Bends"))
  {
    intraMolecularPotentials.bends.reserve(found_bends.size());

    // fill in each of the bends
    for (const std::array<std::size_t, 3> found_bend : found_bends)
    {
      // check which bend-definition matches on this definition
      // only one is allowed (no ambiguities)
      bool bend_identified{ false };

      for (auto &[_, item] : parsed_data["Bends"].items())
      {
        if (!item.is_array())
        {
          throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
        }

        if (item.size() != 3)
        {
          throw std::runtime_error(
              std::format("[Component reader]: item {} must be an array with at least three elements, "
                          "(1) an array with the three identifiers of the bend, (2) the potential type, "
                          "and (3) an array containg the potential parameters\n",
                          item.dump()));
        }

        try
        {
          // parse [1]: the potential name
          std::string potential_name = item[1].get<std::string>();

          // parse [2]: the parameters
          std::vector<double> potential_parameters =
              item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};

          if (item[0][0].is_string() && item[0][1].is_string() && item[0][2].is_string())
          {
            std::vector<std::string> identifiers = item[0].get<std::vector<std::string>>();

            std::optional<std::size_t> type_A = forceField.findPseudoAtom(identifiers[0]);
            std::optional<std::size_t> type_B = forceField.findPseudoAtom(identifiers[1]);
            std::optional<std::size_t> type_C = forceField.findPseudoAtom(identifiers[2]);

            if (!(type_A.has_value() && type_B.has_value() && type_C.has_value()))
            {
              throw std::runtime_error(
                  std::format("[Component reader]: unknown pseudo-atom {} {} {}\n", identifiers[0], identifiers[1], identifiers[2]));
            }

            if ((atoms[found_bend[0]].type == type_A.value() && atoms[found_bend[1]].type == type_B.value() && atoms[found_bend[2]].type == type_C.value()) ||
                (atoms[found_bend[0]].type == type_C.value() && atoms[found_bend[1]].type == type_B.value() && atoms[found_bend[2]].type == type_A.value()))
            {
              if (bend_identified)
              {
                throw std::runtime_error(std::format("Error in Bend-potential: ambiguity in define bend-potential (bend {}-{}-{} matches on multiple definitions)\n",
                      found_bend[0], found_bend[1], found_bend[2]));
              }

              bend_identified = true;

              BendPotential bend_potential =
                  BendPotential({found_bend[0], found_bend[1], found_bend[2]}, BendPotential::definitionForString.at(potential_name),
                                potential_parameters);
              bend_potentials.push_back(bend_potential);
            }
          }
          else if (item[0][0].is_number_unsigned() && item[0][1].is_number_unsigned() && item[0][2].is_number_unsigned())
          {
            std::size_t A = item[0][0].get<std::size_t>();
            std::size_t B = item[0][1].get<std::size_t>();
            std::size_t C = item[0][2].get<std::size_t>();
            if ((found_bend[0] == A && found_bend[1] == B && found_bend[2] == C) ||
                (found_bend[0] == C && found_bend[1] == B && found_bend[2] == A))
            {
              if (bend_identified)
              {
                throw std::runtime_error(std::format("Error in Bend-potential: ambiguity in define bend-potential (bend {}-{}-{} matches on multiple definitions)\n",
                      found_bend[0], found_bend[1], found_bend[2]));
              }

              bend_identified = true;

              // parse using unsigned integer identifiers
              std::vector<std::size_t> bend =
                  item[0].is_array() ? item[0].get<std::vector<std::size_t>>() : std::vector<std::size_t>{};

              BendPotential bend_potential = BendPotential(
                  {bend[0], bend[1], bend[2]}, BendPotential::definitionForString.at(potential_name), potential_parameters);
              bend_potentials.push_back(bend_potential);
            }
          }

        }
        catch (std::exception const &e)
        {
          throw std::runtime_error(std::format("Error in Bend-potential ({}): {}\n", item.dump(), e.what()));
        }


      }
    }
  }

  if (found_bends.size() != bend_potentials.size())
  {
    throw std::runtime_error(std::format("Error in Bend-potential: not all bend-potentials defined ({} missing)\n",
                                         found_bends.size() - bend_potentials.size()));
  }

  return bend_potentials;
}

std::vector<TorsionPotential> Component::readTorsionPotentials(
    const ForceField &forceField, const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<TorsionPotential> torsion_potentials{};

  std::vector<std::array<std::size_t, 4>> found_torsions = connectivityTable.findAllTorsions();

  if (parsed_data.contains("Torsions"))
  {
    intraMolecularPotentials.torsions.reserve(found_torsions.size());

    // fill in each of the torsions
    for (const std::array<std::size_t, 4> found_torsion : found_torsions)
    {
      // check which torsion-definition matches on this definition
      // only one is allowed (no ambiguities)
      bool torsion_identified{ false };

      for (auto &[_, item] : parsed_data["Torsions"].items())
      {
        if (!item.is_array())
        {
          throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
        }

        if (item.size() != 3)
        {
          throw std::runtime_error(
              std::format("[Component reader]: item {} must be an array with at least three elements, "
                          "(1) an array with the four identifiers of the torsion, (2) the potential type, "
                          "and (3) an array containg the potential parameters\n",
                          item.dump()));
        }

        try
        {
          // parse [1]: the potential name
          std::string potential_name = item[1].get<std::string>();

          // parse [2]: the parameters
          std::vector<double> potential_parameters =
              item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};

          if (item[0][0].is_string() && item[0][1].is_string() && item[0][2].is_string())
          {
            std::vector<std::string> identifiers = item[0].get<std::vector<std::string>>();

            std::optional<std::size_t> type_A = forceField.findPseudoAtom(identifiers[0]);
            std::optional<std::size_t> type_B = forceField.findPseudoAtom(identifiers[1]);
            std::optional<std::size_t> type_C = forceField.findPseudoAtom(identifiers[2]);
            std::optional<std::size_t> type_D = forceField.findPseudoAtom(identifiers[3]);

            if (!(type_A.has_value() && type_B.has_value() && type_C.has_value() && type_D.has_value()))
            {
              throw std::runtime_error(
                  std::format("[Component reader]: unknown pseudo-atom {} {} {} {}\n", identifiers[0], identifiers[1], identifiers[2], identifiers[3]));
            }

            if ((atoms[found_torsion[0]].type == type_A.value() && atoms[found_torsion[1]].type == type_B.value() && 
                 atoms[found_torsion[2]].type == type_C.value() && atoms[found_torsion[3]].type == type_D.value()) ||
                (atoms[found_torsion[0]].type == type_D.value() && atoms[found_torsion[1]].type == type_C.value() && 
                 atoms[found_torsion[2]].type == type_B.value() && atoms[found_torsion[3]].type == type_A.value()))
            {
              if (torsion_identified)
              {
                throw std::runtime_error(std::format("Error in Torsion-potential: ambiguity in define torsion-potential (torsion {}-{}-{}-{} matches on multiple definitions)\n",
                      found_torsion[0], found_torsion[1], found_torsion[2], found_torsion[3]));
              }

              torsion_identified = true;

              TorsionPotential torsion_potential =
                  TorsionPotential({found_torsion[0], found_torsion[1], found_torsion[2], found_torsion[3]}, TorsionPotential::definitionForString.at(potential_name),
                                potential_parameters);
              torsion_potentials.push_back(torsion_potential);
            }
          }
          else if (item[0][0].is_number_unsigned() && item[0][1].is_number_unsigned() && item[0][2].is_number_unsigned() && item[0][3].is_number_unsigned())
          {
            std::size_t A = item[0][0].get<std::size_t>();
            std::size_t B = item[0][1].get<std::size_t>();
            std::size_t C = item[0][2].get<std::size_t>();
            std::size_t D = item[0][3].get<std::size_t>();
            if ((found_torsion[0] == A && found_torsion[1] == B && found_torsion[2] == C && found_torsion[3] == D) ||
                (found_torsion[0] == D && found_torsion[1] == C && found_torsion[2] == B && found_torsion[3] == A))
            {
              if (torsion_identified)
              {
                throw std::runtime_error(std::format("Error in Torsion-potential: ambiguity in define torsion-potential (torsion {}-{}-{}-{} matches on multiple definitions)\n",
                      found_torsion[0], found_torsion[1], found_torsion[2], found_torsion[3]));
              }

              torsion_identified = true;

              // parse using unsigned integer identifiers
              std::vector<std::size_t> torsion =
                  item[0].is_array() ? item[0].get<std::vector<std::size_t>>() : std::vector<std::size_t>{};

              TorsionPotential torsion_potential = TorsionPotential(
                  {torsion[0], torsion[1], torsion[2], torsion[3]}, TorsionPotential::definitionForString.at(potential_name), potential_parameters);
              torsion_potentials.push_back(torsion_potential);
            }
          }

        }
        catch (std::exception const &e)
        {
          throw std::runtime_error(std::format("Error in Torsion-potential ({}): {}\n", item.dump(), e.what()));
        }


      }
    }
  }

  if (found_torsions.size() != torsion_potentials.size())
  {
    throw std::runtime_error(std::format("Error in Torsion-potential: not all torsion-potentials defined ({} missing)\n",
                                         found_torsions.size() - torsion_potentials.size()));
  }

  return torsion_potentials;
}

std::vector<VanDerWaalsPotential> Component::readVanDerWaalsPotentials(
    const ForceField &forceField, [[maybe_unused]] const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<VanDerWaalsPotential> van_der_waals_potentials{};

  if (parsed_data.contains("Intra14VanDerWaalsScalingValue"))
  {
    if(parsed_data["Intra14VanDerWaalsScalingValue"].is_number_float())
    {
      double scaling = parsed_data["Intra14VanDerWaalsScalingValue"].get<double>();

      if(scaling > 0.0)
      {
        std::vector<std::array<std::size_t, 4>> found_14_van_der_waals = connectivityTable.findAllTorsions();
        for (std::array<std::size_t, 4> &found_14_van_der_waal : found_14_van_der_waals)
        {
          std::size_t A = found_14_van_der_waal[0];
          std::size_t B = found_14_van_der_waal[3];
          std::size_t typeA = static_cast<std::size_t>(atoms[A].type);
          std::size_t typeB = static_cast<std::size_t>(atoms[B].type);
          
          [[maybe_unused]] VDWParameters::Type potentialType = forceField(typeA, typeB).type;
          double4 parameters = forceField(typeA, typeB).parameters;
          double shift = forceField(typeA, typeB).shift;

          // FIX: unit conversion
          VanDerWaalsPotential potential = VanDerWaalsPotential(
              {A, B}, VanDerWaalsType::LennardJones,
              {parameters.x * Units::EnergyToKelvin, parameters.y, parameters.z, parameters.w}, scaling);

          van_der_waals_potentials.push_back(potential);
        }
      }
    }
  }

  std::vector<std::array<std::size_t, 2>> found_van_der_waals = connectivityTable.findAllVanDerWaals();

  for (std::array<std::size_t, 2> &found_van_der_waal : found_van_der_waals)
  {
    std::size_t A = found_van_der_waal[0];
    std::size_t B = found_van_der_waal[1];
    std::size_t typeA = static_cast<std::size_t>(atoms[A].type);
    std::size_t typeB = static_cast<std::size_t>(atoms[B].type);

    [[maybe_unused]] VDWParameters::Type potentialType = forceField(typeA, typeB).type;
    double4 parameters = forceField(typeA, typeB).parameters;

    // FIX: unit conversion
    VanDerWaalsPotential potential = VanDerWaalsPotential(
        {A, B}, VanDerWaalsType::LennardJones,
        {parameters.x * Units::EnergyToKelvin, parameters.y, parameters.z, parameters.w}, 1.0);

    van_der_waals_potentials.push_back(potential);
  }

  return van_der_waals_potentials;
}

std::vector<CoulombPotential> Component::readCoulombPotentials(
    const ForceField &forceField, [[maybe_unused]] const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<CoulombPotential> coulomb_potentials{};

  if(!forceField.useCharge) return coulomb_potentials;

  if (parsed_data.contains("Intra14ChargeChargeScalingValue"))
  {
    if(parsed_data["Intra14ChargeChargeScalingValue"].is_number_float())
    {
      double scaling = parsed_data["Intra14ChargeChargeScalingValue"].get<double>();

      if(scaling > 0.0)
      {
        std::vector<std::array<std::size_t, 4>> found_14_coulombs = connectivityTable.findAllTorsions();
        for (std::array<std::size_t, 4> &found_14_coulomb : found_14_coulombs)
        {
          std::size_t A = found_14_coulomb[0];
          std::size_t B = found_14_coulomb[3];
          double chargeA = atoms[A].charge;
          double chargeB = atoms[B].charge;
          
          CoulombPotential potential = CoulombPotential({A, B}, CoulombType::Coulomb, chargeA, chargeB, scaling);

          coulomb_potentials.push_back(potential);
        }
      }
    }
  }

  std::vector<std::array<std::size_t, 2>> found_coulombs = connectivityTable.findAllVanDerWaals();

  for (std::array<std::size_t, 2> &found_coulomb : found_coulombs)
  {
    std::size_t A = found_coulomb[0];
    std::size_t B = found_coulomb[1];
    double chargeA = atoms[A].charge;
    double chargeB = atoms[B].charge;

    CoulombPotential potential = CoulombPotential({A, B}, CoulombType::Coulomb, chargeA, chargeB, 1.0);

    coulomb_potentials.push_back(potential);
  }

  return coulomb_potentials;
}

std::vector<std::vector<std::size_t>> Component::readPartialReinsertionFixedAtoms(
    const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data)
{
  std::vector<std::vector<std::size_t>> config_moves{};

  if (parsed_data.contains("Partial-reinsertion"))
  {
    auto item = parsed_data["Partial-reinsertion"];
    if (!item.is_array())
    {
      throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
    }

    try
    {
      config_moves = parsed_data["Partial-reinsertion"].get<std::vector<std::vector<std::size_t>>>();
    }
    catch (std::exception const &e)
    {
      throw std::runtime_error(std::format("Error in defined partial reinsertion ({}): {}\n", item.dump(), e.what()));
    }
  }

  for (const std::vector<std::size_t> &config_move : config_moves)
  {
    if (!connectivityTable.checkIsConnectedSubgraph(config_move))
    {
      throw std::runtime_error(
          std::format("Error in defined partial reinsertion ({} is not connected)\n", config_move));
    }
  }

  for (const std::vector<std::size_t> &config_move : config_moves)
  {
    std::vector<std::size_t> beads_already_placed = config_move;

    do
    {
      std::optional<std::vector<std::size_t>> nextBeads =
          connectivityTable.checkValidityNextBeads(beads_already_placed);
      if (nextBeads.has_value())
      {
        beads_already_placed.insert(beads_already_placed.end(), nextBeads->begin(), nextBeads->end());
      }
      else
      {
        throw std::runtime_error(
            std::format("Error in defined partial reinsertion\n"
                        "{} does not satisfy requirement that all branches need to be grown at the same time\n",
                        config_move));
      }
    } while (beads_already_placed.size() < connectivityTable.numberOfBeads);
  }

  return config_moves;
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

  archive << c.blockingPockets;

  archive << c.rigid;
  archive << c.translationalDegreesOfFreedom;
  archive << c.rotationalDegreesOfFreedom;

  archive << c.criticalTemperature;
  archive << c.criticalPressure;
  archive << c.acentricFactor;
  archive << c.molFraction;
  archive << c.swappable;
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

  archive << c.inertiaVector;
  archive << c.inverseInertiaVector;
  archive << c.shapeType;
  archive << c.atoms;

  archive << c.connectivityTable;
  archive << c.intraMolecularPotentials;
  archive << c.grownAtoms;
  archive << c.partialReinsertionFixedAtoms;

  archive << c.initialNumberOfMolecules;

  archive << c.lambdaGC;
  archive << c.lambdaGibbs;
  archive << c.hasFractionalMolecule;

  archive << c.intraMolecularPotentials;

  archive << c.mc_moves_probabilities;
  archive << c.mc_moves_statistics;
  archive << c.mc_moves_cputime;

  archive << c.cbmc_moves_statistics;

  archive << c.averageRosenbluthWeights;

  archive << c.lnPartitionFunction;

  archive << c.isotherm;
  archive << c.massTransferCoefficient;
  archive << c.axialDispersionCoefficient;
  archive << c.isCarrierGas;

  archive << c.columnPressure;
  archive << c.columnLoading;
  archive << c.columnError;

  archive << c.pressureScale;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Component &c)
{
  std::uint64_t versionNumber;
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

  archive >> c.blockingPockets;

  archive >> c.rigid;
  archive >> c.translationalDegreesOfFreedom;
  archive >> c.rotationalDegreesOfFreedom;

  archive >> c.criticalTemperature;
  archive >> c.criticalPressure;
  archive >> c.acentricFactor;
  archive >> c.molFraction;
  archive >> c.swappable;
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

  archive >> c.inertiaVector;
  archive >> c.inverseInertiaVector;
  archive >> c.shapeType;
  archive >> c.atoms;

  archive >> c.connectivityTable;
  archive >> c.intraMolecularPotentials;
  archive >> c.grownAtoms;
  archive >> c.partialReinsertionFixedAtoms;

  archive >> c.initialNumberOfMolecules;

  archive >> c.lambdaGC;
  archive >> c.lambdaGibbs;
  archive >> c.hasFractionalMolecule;

  archive >> c.intraMolecularPotentials;

  archive >> c.mc_moves_probabilities;
  archive >> c.mc_moves_statistics;
  archive >> c.mc_moves_cputime;

  archive >> c.cbmc_moves_statistics;

  archive >> c.averageRosenbluthWeights;

  archive >> c.lnPartitionFunction;

  archive >> c.isotherm;
  archive >> c.massTransferCoefficient;
  archive >> c.axialDispersionCoefficient;
  archive >> c.isCarrierGas;

  archive >> c.columnPressure;
  archive >> c.columnLoading;
  archive >> c.columnError;

  archive >> c.pressureScale;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Component: Error in binary restart\n"));
  }
#endif

  return archive;
}

std::string Component::repr() const { return std::string("Component test"); }
