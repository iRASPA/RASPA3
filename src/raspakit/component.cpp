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
#include <map>
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
import internal_potentials;
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

  connectivityTable = ConnectivityTable(definedAtoms.size());

  for (const BondPotential &bond_potential : internalPotentials.bonds)
  {
    std::size_t A = bond_potential.identifiers[0];
    std::size_t B = bond_potential.identifiers[1];
    connectivityTable[A, B] = true;
    connectivityTable[B, A] = true;
  }
}

// create programmatically an 'adsorbate' component
Component::Component(std::size_t componentId, const ForceField &forceField, std::string componentName, double T_c,
                     double P_c, double w, std::vector<Atom> atomList, std::size_t numberOfBlocks,
                     std::size_t numberOfLambdaBins, const MCMoveProbabilities &particleProbabilities,
                     std::optional<double> fugacityCoefficient, bool thermodynamicIntegration) noexcept(false)
    : type(Type::Adsorbate),
      componentId(componentId),
      name(componentName),
      criticalTemperature(T_c),
      criticalPressure(P_c),
      acentricFactor(w),
      fugacityCoefficient(fugacityCoefficient),
      lambdaGC(numberOfBlocks, numberOfLambdaBins),
      lambdaGibbs(numberOfBlocks, numberOfLambdaBins),
      mc_moves_probabilities(particleProbabilities),
      averageRosenbluthWeights(numberOfBlocks),
      connectivityTable(atomList.size())
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

  for (const BondPotential &bond_potential : internalPotentials.bonds)
  {
    std::size_t A = bond_potential.identifiers[0];
    std::size_t B = bond_potential.identifiers[1];
    connectivityTable[A, B] = true;
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
          std::format("[Component reader]: unknown pseudo-atom '{}', please lookup type in in 'pseudo_atoms.json'\n",
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

    definedAtoms.push_back({Atom(double3(position[0], position[1], position[2]), charge, scaling, 0,
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
    if(starting_bead >= 0)
    {
      startingBead = static_cast<std::size_t>(starting_bead);
    }

  }

  if (parsed_data.contains("Type"))
  {
    if (!parsed_data["Type"].is_string())
    {
      throw std::runtime_error(
          std::format("[Component reader]: item {} must be an string (the name of the pseudo-atom)\n", parsed_data["Type"].dump()));
    }
    std::string typeString = parsed_data["Type"].get<std::string>();
    if (caseInSensStringCompare(typeString, "Flexible"))
    {
      growType = GrowType::Flexible;
    }
    else if (caseInSensStringCompare(typeString, "Rigid"))
    {
      growType = GrowType::Rigid;
    }
  }

  // Read bonds
  if (parsed_data.contains("Bonds"))
  {
    internalPotentials.bonds.reserve(parsed_data["Bonds"].size());

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
        std::vector<std::size_t> identifiers =
            item[0].is_array() ? item[0].get<std::vector<std::size_t>>() : std::vector<std::size_t>{};
        std::string potential_name = item[1].get<std::string>();
        std::vector<double> potential_parameters =
            item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};
        BondPotential bond = BondPotential({identifiers[0], identifiers[1]},
                                           BondPotential::definitionForString.at(potential_name), potential_parameters);

        internalPotentials.bonds.push_back(bond);
      }
      catch (std::exception const &e)
      {
        throw std::runtime_error(std::format("Error in Bond-potential ({}): {}\n", item.dump(), e.what()));
      }
    }
  }

  // Read bends
  if (parsed_data.contains("Bends"))
  {
    internalPotentials.bends.reserve(parsed_data["Bends"].size());

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
        std::vector<std::size_t> identifiers =
            item[0].is_array() ? item[0].get<std::vector<std::size_t>>() : std::vector<std::size_t>{};
        std::string potential_name = item[1].get<std::string>();
        std::vector<double> potential_parameters =
            item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};

        BendPotential bend = BendPotential({identifiers[0], identifiers[1], identifiers[2]},
                                           BendPotential::definitionForString.at(potential_name), potential_parameters);

        internalPotentials.bends.push_back(bend);
      }
      catch (std::exception const &e)
      {
        throw std::runtime_error(std::format("Error in Bend-potential ({}): {}\n", item.dump(), e.what()));
      }
    }
  }

  // Read torsions
  if (parsed_data.contains("Torsions"))
  {
    internalPotentials.torsions.reserve(parsed_data["Torsions"].size());

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
                        "(1) an array with the four identifiers of the dihedrals, (2) the potential type, "
                        "and (3) an array containg the potential parameters\n",
                        item.dump()));
      }

      try
      {
        std::vector<std::size_t> identifiers =
            item[0].is_array() ? item[0].get<std::vector<std::size_t>>() : std::vector<std::size_t>{};
        std::string potential_name = item[1].get<std::string>();
        std::vector<double> potential_parameters =
            item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};

        TorsionPotential torsion =
            TorsionPotential({identifiers[0], identifiers[1], identifiers[2], identifiers[3]},
                             TorsionPotential::definitionForString.at(potential_name), potential_parameters);

        internalPotentials.torsions.push_back(torsion);
      }
      catch (std::exception const &e)
      {
        throw std::runtime_error(std::format("Error in Torsion-potential ({}): {}\n", item.dump(), e.what()));
      }
    }
  }

  // Read Intra Van der Waals
  if (parsed_data.contains("VanDerWaals"))
  {
    internalPotentials.torsions.reserve(parsed_data["VanDerWaals"].size());

    for (auto &[_, item] : parsed_data["VanDerWaals"].items())
    {
      if (!item.is_array())
      {
        throw std::runtime_error(std::format("[Component reader]: item {} must be an array\n", item.dump()));
      }

      if (item.size() == 2)
      {
        // parse 'identifier, identifier'
        std::size_t A = item[0].get<std::size_t>();
        std::size_t typeA = atoms[A].type;
        std::size_t B = item[1].get<std::size_t>();
        std::size_t typeB = atoms[B].type;
        VDWParameters parameters = forceField(A, B);
      }
      else if (item.size() == 3)
      {
        // parse '[identifier, identifier], 1.0'
      }
      else if (item.size() == 5)
      {
        // parse '[identifier, identifier], 1.0, "LENNARD_JONES", [epsilon, sigma]'
      }

      try
      {
        std::vector<std::size_t> identifiers =
            item[0].is_array() ? item[0].get<std::vector<std::size_t>>() : std::vector<std::size_t>{};
        std::string potential_name = item[1].get<std::string>();
        std::vector<double> potential_parameters =
            item[2].is_array() ? item[2].get<std::vector<double>>() : std::vector<double>{};

        TorsionPotential torsion =
            TorsionPotential({identifiers[0], identifiers[1], identifiers[2], identifiers[3]},
                             TorsionPotential::definitionForString.at(potential_name), potential_parameters);

        internalPotentials.torsions.push_back(torsion);
      }
      catch (std::exception const &e)
      {
        throw std::runtime_error(std::format("Error in Torsion-potential ({}): {}\n", item.dump(), e.what()));
      }
    }
  }

  totalMass = 0.0;
  netCharge = 0.0;
  for (auto [atom, mass] : definedAtoms)
  {
    totalMass += mass;
    netCharge += atom.charge;
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

  std::size_t index = 0;
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
  for (std::size_t i = 0; i < atoms.size(); ++i)
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

std::string Component::printStatus(const ForceField &forceField) const
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
  if (fugacityCoefficient.has_value())
  {
    std::print(stream, "    Fugacity coefficient:         {} [-]\n", fugacityCoefficient.value());
  }
  std::print(stream, "    Bulk fluid density:           {} [-]\n", bulkFluidDensity);
  std::print(stream, "    Compressibility:              {} [-]\n", compressibility);
  std::print(stream, "    Excess molecules:             {} [-]\n\n", amountOfExcessMolecules);

  std::print(stream, "    Number Of Atoms:    {}\n", atoms.size());
  std::print(stream, "    CBMC starting bead: {}\n", startingBead);
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

  const std::map<MoveTypes, double> normalizedProbabilities = mc_moves_probabilities.normalizedMap();
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

  switch(growType)
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

  if(growType == GrowType::Flexible)
  {
    std::print(stream, "    connectivity:\n");
    std::print(stream, "{}\n", connectivityTable.print("        "));

    std::print(stream, "    number of bond potentials: {}\n", internalPotentials.bonds.size());
    for (std::size_t i = 0; i < internalPotentials.bonds.size(); ++i)
    {
      std::print(stream, "        {}", internalPotentials.bonds[i].print());
    }
    std::print(stream, "\n");

    if(!internalPotentials.ureyBradleys.empty())
    {
      std::print(stream, "    number of Urey-Bradley potentials: {}\n", internalPotentials.ureyBradleys.size());
      for (std::size_t i = 0; i < internalPotentials.ureyBradleys.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.ureyBradleys[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.bends.empty())
    {
      std::print(stream, "    number of bend potentials: {}\n", internalPotentials.bends.size());
      for (std::size_t i = 0; i < internalPotentials.bends.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.bends[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.inversionBends.empty())
    {
      std::print(stream, "    number of inversion-bend potentials: {}\n", internalPotentials.inversionBends.size());
      for (std::size_t i = 0; i < internalPotentials.inversionBends.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.inversionBends[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.outOfPlaneBends.empty())
    {
      std::print(stream, "    number of out-of-plane bend potentials: {}\n", internalPotentials.outOfPlaneBends.size());
      for (std::size_t i = 0; i < internalPotentials.outOfPlaneBends.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.outOfPlaneBends[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.torsions.empty())
    {
      std::print(stream, "    number of torsion potentials: {}\n", internalPotentials.torsions.size());
      for (std::size_t i = 0; i < internalPotentials.torsions.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.torsions[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.improperTorsions.empty())
    {
      std::print(stream, "    number of improper-torsion potentials: {}\n", internalPotentials.improperTorsions.size());
      for (std::size_t i = 0; i < internalPotentials.improperTorsions.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.improperTorsions[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.bondBonds.empty())
    {
      std::print(stream, "    number of bond-bond potentials: {}\n", internalPotentials.bondBonds.size());
      for (std::size_t i = 0; i < internalPotentials.bondBonds.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.bondBonds[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.bondBends.empty())
    {
      std::print(stream, "    number of bond-bend potentials: {}\n", internalPotentials.bondBends.size());
      for (std::size_t i = 0; i < internalPotentials.bondBends.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.bondBends[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.bondTorsions.empty())
    {
      std::print(stream, "    number of bond-torsion potentials: {}\n", internalPotentials.bondTorsions.size());
      for (std::size_t i = 0; i < internalPotentials.bondTorsions.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.bondTorsions[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.bendBends.empty())
    {
      std::print(stream, "    number of bend-bend potentials: {}\n", internalPotentials.bendBends.size());
      for (std::size_t i = 0; i < internalPotentials.bendBends.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.bendBends[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.bendTorsions.empty())
    {
      std::print(stream, "    number of bend-torsion potentials: {}\n", internalPotentials.bendTorsions.size());
      for (std::size_t i = 0; i < internalPotentials.bendTorsions.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.bendTorsions[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.vanDerWaals.empty())
    {
      std::print(stream, "    number of Van der Waals potentials: {}\n", internalPotentials.vanDerWaals.size());
      for (std::size_t i = 0; i < internalPotentials.vanDerWaals.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.vanDerWaals[i].print());
      }
      std::print(stream, "\n");
    }

    if(!internalPotentials.coulombs.empty())
    {
      std::print(stream, "    number of coulomb potentials: {}\n", internalPotentials.coulombs.size());
      for (std::size_t i = 0; i < internalPotentials.coulombs.size(); ++i)
      {
        std::print(stream, "        {}", internalPotentials.coulombs[i].print());
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
  std::map<MoveTypes, double> normalizedProbabilities = mc_moves_probabilities.normalizedMap();
  for (auto &[moveType, probability] : normalizedProbabilities)
  {
    moves[moveNames[moveType]] = probability;
  }
  status["moveProbabilities"] = moves;

  status["n_bonds"] = internalPotentials.bonds.size();
  std::vector<std::string> bondTypes(internalPotentials.bonds.size());
  for (std::size_t i = 0; i < internalPotentials.bonds.size(); ++i)
  {
    bondTypes[i] = internalPotentials.bonds[i].print();
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
    for (std::size_t i = 0; i != trialAtoms.size(); ++i)
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

  archive << c.blockingPockets;

  archive << c.rigid;

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
  archive << c.atoms;

  archive << c.initialNumberOfMolecules;

  archive << c.lambdaGC;
  archive << c.lambdaGibbs;
  archive << c.hasFractionalMolecule;

  archive << c.internalPotentials;

  archive << c.connectivityTable;

  archive << c.mc_moves_probabilities;
  archive << c.mc_moves_statistics;
  archive << c.mc_moves_cputime;

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
  archive >> c.atoms;

  archive >> c.initialNumberOfMolecules;

  archive >> c.lambdaGC;
  archive >> c.lambdaGibbs;
  archive >> c.hasFractionalMolecule;

  archive >> c.internalPotentials;

  archive >> c.connectivityTable;

  archive >> c.mc_moves_probabilities;
  archive >> c.mc_moves_statistics;
  archive >> c.mc_moves_cputime;

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
