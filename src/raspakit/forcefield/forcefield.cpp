module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <print>
#include <set>
#include <source_location>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>
#endif

module forcefield;

#ifdef USE_STD_IMPORT
import std;
#endif

import archive;
import json;
import skelement;
import units;
import int3;
import double3;
import double4;
import double3x3;
import stringutils;
import pseudo_atom;
import vdwparameters;
import potential_correction_vdw;
import potential_correction_pressure;
import simulationbox;
import json;

int3 parseInt3(const std::string& item, auto json)
{
  if (json.is_array())
  {
    if (json.size() != 3)
    {
      throw std::runtime_error(
          std::format("[Input reader]: key '{}', value {} should be array of 3 integer numbers\n", item, json.dump()));
    }
    int3 value{};
    try
    {
      value.x = json[0].template get<std::int32_t>();
      value.y = json[1].template get<std::int32_t>();
      value.z = json[2].template get<std::int32_t>();
      return value;
    }
    catch (nlohmann::json::exception& ex)
    {
      throw std::runtime_error(
          std::format("[Input reader]: key '{}', value {} should be array of 3 integer numbers\n", item, json.dump()));
    }
  }
  throw std::runtime_error(
      std::format("[Input reader]: key '{}', value {} should be array of 3 integer  numbers\n", item, json.dump()));
}

uint3 parseUint3(const std::string& item, auto json)
{
  if (json.is_array())
  {
    if (json.size() != 3)
    {
      throw std::runtime_error(
          std::format("[Input reader]: key '{}', value {} should be array of 3 integer numbers\n", item, json.dump()));
    }
    uint3 value{};
    try
    {
      value.x = json[0].template get<std::size_t>();
      value.y = json[1].template get<std::size_t>();
      value.z = json[2].template get<std::size_t>();
      return value;
    }
    catch (nlohmann::json::exception& ex)
    {
      throw std::runtime_error(
          std::format("[Input reader]: key '{}', value {} should be array of 3 integer numbers\n", item, json.dump()));
    }
  }
  throw std::runtime_error(
      std::format("[Input reader]: key '{}', value {} should be array of 3 integer  numbers\n", item, json.dump()));
}

double3 parseDouble3(const std::string& item, auto json)
{
  if (json.is_array())
  {
    if (json.size() != 3)
    {
      throw std::runtime_error(std::format(
          "[Input reader]: key '{}', value {} should be array of 3 floatng point numbers\n", item, json.dump()));
    }
    double3 value{};
    try
    {
      value.x = json[0].template get<double>();
      value.y = json[1].template get<double>();
      value.z = json[2].template get<double>();
      return value;
    }
    catch (nlohmann::json::exception& ex)
    {
      throw std::runtime_error(std::format(
          "[Input reader]: key '{}', value {} should be array of 3 floating point numbers\n", item, json.dump()));
    }
  }
  throw std::runtime_error(std::format(
      "[Input reader]: key '{}', value {} should be array of 3 floating point numbers\n", item, json.dump()));
}


ForceField::ForceField(std::vector<PseudoAtom> pseudoAtoms, std::vector<VDWParameters> selfInteractions,
                       MixingRule mixingRule, double cutOffFrameworkVDW, double cutOffMoleculeVDW, double cutOffCoulomb,
                       bool shifted, bool applyTailCorrections, bool useCharge) noexcept(false)
    : data((pseudoAtoms.size() + 1) * (pseudoAtoms.size() + 1), VDWParameters(0.0, 0.0)),
      shiftPotentials((pseudoAtoms.size() + 1) * (pseudoAtoms.size() + 1), shifted),
      tailCorrections((pseudoAtoms.size() + 1) * (pseudoAtoms.size() + 1), applyTailCorrections),
      mixingRule(mixingRule),
      cutOffFrameworkVDW(cutOffFrameworkVDW),
      cutOffMoleculeVDW(cutOffMoleculeVDW),
      cutOffCoulomb(cutOffCoulomb),
      numberOfPseudoAtoms(pseudoAtoms.size() + 1),
      pseudoAtoms(pseudoAtoms),
      useCharge(useCharge)
{
  PseudoAtom unitPseudoAtom("unit", false, 1.0, 1.0, 0.0, 0, false);
  this->pseudoAtoms.push_back(unitPseudoAtom);

  VDWParameters unitVDWParameters(double4(1.0, 1.0, 0.0, 0.0), 0.0, 0.0, 0.0, VDWParameters::Type::LennardJones);
  selfInteractions.push_back(unitVDWParameters);

  for (std::size_t i = 0; i < selfInteractions.size(); ++i)
  {
    data[i + i * numberOfPseudoAtoms] = selfInteractions[i];
  }

  applyMixingRule();
  preComputePotentialShift();
  preComputeTailCorrection();
}

ForceField::ForceField(std::string filePath)
{
  // Construct the path to the force field file
  std::filesystem::path forceFieldPathfile = std::filesystem::path(filePath);

  // Check if the file exists
  if (!std::filesystem::exists(forceFieldPathfile))
  {
    throw std::runtime_error(std::format("[Forcefield reader]: File {} does not exist.\n", filePath));
  }

  // Open the file
  std::ifstream forceFieldStream{forceFieldPathfile};
  if (!forceFieldStream)
  {
    throw std::runtime_error(std::format("[Forcefield reader]: File {} is empty.\n", filePath));
  }

  // Parse the JSON data
  nlohmann::basic_json<nlohmann::raspa_map> parsed_data{};
  try
  {
    parsed_data = nlohmann::json::parse(forceFieldStream);
  }
  catch (nlohmann::json::parse_error& ex)
  {
    throw std::runtime_error(
        std::format("[Forcefield reader]: Parse error of file {} at byte {}\n{}\n", filePath, ex.byte, ex.what()));
  }

  // validate input
  validateInput(parsed_data);

  // read pseudo atoms
  if (!parsed_data.contains("PseudoAtoms"))
  {
    throw std::runtime_error("[Forcefield reader]: No pseudo-atoms found [keyword 'PseudoAtoms' missing]\n");
  }

  numberOfPseudoAtoms = parsed_data["PseudoAtoms"].size();
  if (numberOfPseudoAtoms == 0)
  {
    throw std::runtime_error("[ReadPseudoAtoms]: key 'PseudoAtoms' empty]\n");
  }

  pseudoAtoms.reserve(numberOfPseudoAtoms);
  for (const auto& [_, item] : parsed_data["PseudoAtoms"].items())
  {
    std::string jsonName = item.value("name", "");
    bool jsonFrameworkType = item.value("framework", false);
    std::int64_t oxidationState = item.value("oxidationState", 0);
    double jsonMass = item.value("mass", 0.0);
    std::string jsonElement = item.value("element", "C");
    double jsonCharge = item.value("charge", 0.0);
    double jsonPolarizibility = item.value("polarizibility", 0.0);
    bool jsonPrintToOutput = item.value("print_to_output", true);
    std::string jsonSource = item.value("source", "");

    std::size_t atomicNumber = PredefinedElements::atomicNumberData.contains(jsonElement)
                                   ? static_cast<std::size_t>(PredefinedElements::atomicNumberData.at(jsonElement))
                                   : 1;

    PseudoAtom pseudo_atom = PseudoAtom(jsonName, jsonFrameworkType, jsonMass, jsonCharge, jsonPolarizibility,
                                        atomicNumber, jsonPrintToOutput, jsonSource);
    pseudo_atom.oxidationState = oxidationState;

    pseudoAtoms.push_back(pseudo_atom);
  }

  // Read and set truncation methods
  bool shiftPotential = false;
  bool tailCorrection = false;
  if (parsed_data.value("TruncationMethod", "") == "shifted")
  {
    shiftPotential = true;
  }
  tailCorrection = parsed_data.value("TailCorrections", false);

  data.resize(numberOfPseudoAtoms * numberOfPseudoAtoms, VDWParameters());
  shiftPotentials.resize(numberOfPseudoAtoms * numberOfPseudoAtoms, shiftPotential);
  tailCorrections.resize(numberOfPseudoAtoms * numberOfPseudoAtoms, tailCorrection);

  // Read self-interactions
  for (const auto& [_, item] : parsed_data["SelfInteractions"].items())
  {
    std::string jsonName = item.value("name", "");
    std::string jsonType = item.value("type", "lennard-jones");
    std::string jsonSource = item.value("source", "");

    std::vector<double> scannedJsonParameters;
    try
    {
      scannedJsonParameters = item.value("parameters", std::vector<double>{});
    }
    catch (const nlohmann::json::exception& ex)
    {
      throw std::runtime_error(
          std::format("[ReadForceFieldSelfInteractions]: parameters {} must be array of numbers \n{}\n",
                      item["parameters"].dump(), ex.what()));
    }

    std::optional<std::size_t> index = ForceField::findPseudoAtom(jsonName);
    if (!index.has_value())
    {
      throw std::runtime_error(
          std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define\n", jsonName));
    }

    if (scannedJsonParameters.size() < 2)
    {
      throw std::runtime_error(
          std::format("[ReadForceFieldSelfInteractions]: incorrect vdw parameters {}\n", item["parameters"].dump()));
    }

    double param0 = scannedJsonParameters[0];
    double param1 = scannedJsonParameters[1];
    VDWParameters::Type type = VDWParameters::stringToEnum(jsonType);

    data[index.value() * numberOfPseudoAtoms + index.value()] = VDWParameters(param0, param1, type);
  }

  // Set mixing rule and cut-off values
  if (parsed_data.value("MixingRule", "") == "Lorentz-Berthelot")
  {
    mixingRule = MixingRule::Lorentz_Berthelot;
  }
  if (parsed_data.value("MixingRule", "") == "Jorgensen")
  {
    mixingRule = MixingRule::Jorgensen;
  }

  // Apply mixing rule and precompute potentials
  applyMixingRule();

  // Read binary interactions
  for (const auto& [_, item] : parsed_data["BinaryInteractions"].items())
  {
    std::string jsonType = item.value("type", "lennard-jones");
    std::string jsonSource = item.value("source", "");

    std::vector<std::string> scannedJsonParameterNames;
    try
    {
      scannedJsonParameterNames = item.value("names", std::vector<std::string>{});
    }
    catch (const nlohmann::json::exception& ex)
    {
      throw std::runtime_error(
          std::format("[ReadForceFieldBinaryInteractions]: names {} must be array of two strings \n{}\n",
                      item["names"].dump(), ex.what()));
    }

    if (scannedJsonParameterNames.size() != 2)
    {
      throw std::runtime_error(
          std::format("[ReadForceFieldBinaryInteractions]: incorrect number of strings {}\n", item["names"].dump()));
    }

    std::optional<std::size_t> indexA = ForceField::findPseudoAtom(scannedJsonParameterNames[0]);
    if (!indexA.has_value())
    {
      throw std::runtime_error(
          std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define'\n",
                      scannedJsonParameterNames[0]));
    }
    std::optional<std::size_t> indexB = ForceField::findPseudoAtom(scannedJsonParameterNames[1]);
    if (!indexB.has_value())
    {
      throw std::runtime_error(
          std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define'\n",
                      scannedJsonParameterNames[1]));
    }

    std::vector<double> scannedJsonParameters;
    try
    {
      scannedJsonParameters = item.value("parameters", std::vector<double>{});
    }
    catch (const nlohmann::json::exception& ex)
    {
      throw std::runtime_error(
          std::format("[ReadForceFieldBinaryInteractions]: parameters {} must be array of numbers \n{}\n",
                      item["parameters"].dump(), ex.what()));
    }

    if (scannedJsonParameters.size() < 2)
    {
      throw std::runtime_error(
          std::format("[ReadForceFieldSelfInteractions]: incorrect vdw parameters {}\n", item["parameters"].dump()));
    }

    double param0 = scannedJsonParameters[0];
    double param1 = scannedJsonParameters[1];
    VDWParameters::Type type = VDWParameters::stringToEnum(jsonType);

    data[indexA.value() * numberOfPseudoAtoms + indexB.value()] = VDWParameters(param0, param1, type);
    data[indexB.value() * numberOfPseudoAtoms + indexA.value()] = VDWParameters(param0, param1, type);
  }

  useCharge = parsed_data.value("ChargeMethod", "Ewald") != "None";

  if (parsed_data.contains("EwaldPrecision"))
  {
    if (parsed_data["EwaldPrecision"].is_string())
    {
      std::string ewaldPrecisionString = parsed_data["EwaldPrecision"].get<std::string>();

      if (caseInSensStringCompare(ewaldPrecisionString, "High"))
      {
        EwaldPrecision = 1e-8;
        automaticEwald = true;
      }
      else if (caseInSensStringCompare(ewaldPrecisionString, "Medium"))
      {
        EwaldPrecision = 1e-6;
        automaticEwald = true;
      }
      else if (caseInSensStringCompare(ewaldPrecisionString, "Low"))
      {
        EwaldPrecision = 1e-4;
        automaticEwald = true;
      }
      else
      {
      }
    }
    else if (parsed_data["EwaldPrecision"].is_number_float())
    {
      EwaldPrecision = parsed_data["EwaldPrecision"].get<double>();
      automaticEwald = true;
    }
  }

  if (parsed_data.contains("EwaldParameters"))
  {
    if (parsed_data["EwaldParameters"].is_array())
    {
      if (parsed_data["EwaldParameters"].size() != 4)
      {
        throw std::runtime_error(
            std::format("[ForceField reader]: key '{}', value {} should be array of one floaying point number and 3 "
                        "integer numbers\n",
                        "EwaldParameters", parsed_data["EwaldParameters"].dump()));
      }
      double alpha;
      int3 value{};
      try
      {
        alpha = parsed_data["EwaldParameters"][0].template get<double>();
        value.x = parsed_data["EwaldParameters"][1].template get<std::int32_t>();
        value.y = parsed_data["EwaldParameters"][2].template get<std::int32_t>();
        value.z = parsed_data["EwaldParameters"][3].template get<std::int32_t>();
      }
      catch (nlohmann::json::exception& ex)
      {
        throw std::runtime_error(
            std::format("[ForceField reader]: key '{}', value {} should be array of one floaying point number and 3 "
                        "integer numbers\n",
                        "EwaldParameters", parsed_data["EwaldParameters"].dump()));
      }

      EwaldAlpha = alpha;
      numberOfWaveVectors = value;
      automaticEwald = false;
    }
  }

  if (parsed_data.contains("CutOff"))
  {
    if (parsed_data["CutOff"].is_string())
    {
      std::string cutOffCoulombString = parsed_data["CutOff"].get<std::string>();

      if (caseInSensStringCompare(cutOffCoulombString, "auto"))
      {
        cutOffFrameworkVDWAutomatic = true;
        cutOffMoleculeVDWAutomatic = true;
      }
    }

    if (parsed_data["CutOff"].is_number_float())
    {
      cutOffFrameworkVDWAutomatic = false;
      cutOffMoleculeVDWAutomatic = false;
      cutOffFrameworkVDW = parsed_data["CutOff"].get<double>();
      cutOffMoleculeVDW = parsed_data["CutOff"].get<double>();
    }
  }

  if (parsed_data.contains("CutOffVDW"))
  {
    if (parsed_data["CutOffVDW"].is_string())
    {
      std::string cutOffString = parsed_data["CutOffVDW"].get<std::string>();

      if (caseInSensStringCompare(cutOffString, "auto"))
      {
        cutOffFrameworkVDWAutomatic = true;
        cutOffMoleculeVDWAutomatic = true;
      }
    }

    if (parsed_data["CutOffVDW"].is_number_float())
    {
      cutOffFrameworkVDWAutomatic = false;
      cutOffMoleculeVDWAutomatic = false;
      cutOffFrameworkVDW = parsed_data["CutOffVDW"].get<double>();
      cutOffMoleculeVDW = parsed_data["CutOffVDW"].get<double>();
    }
  }

  if (parsed_data.contains("CutOffFrameworkVDW"))
  {
    if (parsed_data["CutOffFrameworkVDW"].is_string())
    {
      std::string cutOffCoulombString = parsed_data["CutOffFrameworkVDW"].get<std::string>();

      if (caseInSensStringCompare(cutOffCoulombString, "auto"))
      {
        cutOffFrameworkVDWAutomatic = true;
      }
    }

    if (parsed_data["CutOffFrameworkVDW"].is_number_float())
    {
      cutOffFrameworkVDWAutomatic = false;
      cutOffFrameworkVDW = parsed_data["CutOffFrameworkVDW"].get<double>();
    }
  }

  if (parsed_data.contains("CutOffMoleculeVDW"))
  {
    if (parsed_data["CutOffMoleculeVDW"].is_string())
    {
      std::string cutOffMoleculeVDWString = parsed_data["CutOffMoleculeVDW"].get<std::string>();

      if (caseInSensStringCompare(cutOffMoleculeVDWString, "auto"))
      {
        cutOffMoleculeVDWAutomatic = true;
      }
    }

    if (parsed_data["CutOffMoleculeVDW"].is_number_float())
    {
      cutOffMoleculeVDWAutomatic = false;
      cutOffMoleculeVDW = parsed_data["CutOffMoleculeVDW"].get<double>();
    }
  }

  if (parsed_data.contains("CutOffCoulomb"))
  {
    if (parsed_data["CutOffCoulomb"].is_string())
    {
      std::string cutOffCoulombString = parsed_data["CutOffCoulomb"].get<std::string>();

      if (caseInSensStringCompare(cutOffCoulombString, "auto"))
      {
        cutOffCoulombAutomatic = true;
      }
    }

    if (parsed_data["CutOffCoulomb"].is_number_float())
    {
      cutOffCoulombAutomatic = false;
      cutOffCoulomb = parsed_data["CutOffCoulomb"].get<double>();
    }
  }

  if (parsed_data.contains("ShiftedPotentialPairs"))
  {
    if (parsed_data["ShiftedPotentialPairs"].is_array())
    {
      for (auto& pair : parsed_data["ShiftedPotentialPairs"])
      {
        if (pair.size() != 2)
        {
          throw std::runtime_error(std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n",
                                               "ShiftedPotentialPairs", pair.dump()));
        }

        std::string stringA, stringB;
        try
        {
          stringA = pair[0].template get<std::string>();
          stringB = pair[1].template get<std::string>();
        }
        catch (nlohmann::json::exception& ex)
        {
          throw std::runtime_error(std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n",
                                               "ShiftedPotentialPairs", pair.dump()));
        }

        std::optional<std::size_t> indexA = findPseudoAtom(stringA);
        if (!indexA.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringA, pair.dump()));
        }

        std::optional<std::size_t> indexB = findPseudoAtom(stringB);
        if (!indexB.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringB, pair.dump()));
        }

        shiftPotentials[indexA.value() * numberOfPseudoAtoms + indexB.value()] = true;
        shiftPotentials[indexB.value() * numberOfPseudoAtoms + indexA.value()] = true;
      }
    }
  }

  if (parsed_data.contains("TruncatedPotentialPairs"))
  {
    if (parsed_data["TruncatedPotentialPairs"].is_array())
    {
      for (auto& pair : parsed_data["TruncatedPotentialPairs"])
      {
        if (pair.size() != 2)
        {
          throw std::runtime_error(std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n",
                                               "TruncatedPotentialPairs", pair.dump()));
        }

        std::string stringA, stringB;
        try
        {
          stringA = pair[0].template get<std::string>();
          stringB = pair[1].template get<std::string>();
        }
        catch (nlohmann::json::exception& ex)
        {
          throw std::runtime_error(std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n",
                                               "TruncatedPotentialPairs", pair.dump()));
        }

        std::optional<std::size_t> indexA = findPseudoAtom(stringA);
        if (!indexA.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringA, pair.dump()));
        }

        std::optional<std::size_t> indexB = findPseudoAtom(stringB);
        if (!indexB.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringB, pair.dump()));
        }

        shiftPotentials[indexA.value() * numberOfPseudoAtoms + indexB.value()] = false;
        shiftPotentials[indexB.value() * numberOfPseudoAtoms + indexA.value()] = false;
      }
    }
  }

  if (parsed_data.contains("TailCorrectionPairs"))
  {
    if (parsed_data["TailCorrectionPairs"].is_array())
    {
      for (auto& pair : parsed_data["TailCorrectionPairs"])
      {
        if (pair.size() != 2)
        {
          throw std::runtime_error(std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n",
                                               "TailCorrectionPairs", pair.dump()));
        }

        std::string stringA, stringB;
        try
        {
          stringA = pair[0].template get<std::string>();
          stringB = pair[1].template get<std::string>();
        }
        catch (nlohmann::json::exception& ex)
        {
          throw std::runtime_error(std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n",
                                               "TailCorrectionPairs", pair.dump()));
        }

        std::optional<std::size_t> indexA = findPseudoAtom(stringA);
        if (!indexA.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringA, pair.dump()));
        }

        std::optional<std::size_t> indexB = findPseudoAtom(stringB);
        if (!indexB.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringB, pair.dump()));
        }

        tailCorrections[indexA.value() * numberOfPseudoAtoms + indexB.value()] = true;
        tailCorrections[indexB.value() * numberOfPseudoAtoms + indexA.value()] = true;
      }
    }
  }

  if (parsed_data.contains("NoTailCorrectionPairs"))
  {
    if (parsed_data["NoTailCorrectionPairs"].is_array())
    {
      for (auto& pair : parsed_data["NoTailCorrectionPairs"])
      {
        if (pair.size() != 2)
        {
          throw std::runtime_error(std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n",
                                               "NoTailCorrectionPairs", pair.dump()));
        }

        std::string stringA, stringB;
        try
        {
          stringA = pair[0].template get<std::string>();
          stringB = pair[1].template get<std::string>();
        }
        catch (nlohmann::json::exception& ex)
        {
          throw std::runtime_error(std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n",
                                               "NoTailCorrectionPairs", pair.dump()));
        }

        std::optional<std::size_t> indexA = findPseudoAtom(stringA);
        if (!indexA.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringA, pair.dump()));
        }

        std::optional<std::size_t> indexB = findPseudoAtom(stringB);
        if (!indexB.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringB, pair.dump()));
        }

        tailCorrections[indexA.value() * numberOfPseudoAtoms + indexB.value()] = false;
        tailCorrections[indexB.value() * numberOfPseudoAtoms + indexA.value()] = false;
      }
    }
  }

  // after knowing the tail-correction settings and the cutoff, compute the shifts and tail-corrections for each pair
  preComputePotentialShift();
  preComputeTailCorrection();

  std::vector<std::string> pseudoAtomStringGrids =
      parsed_data.value("UseInterpolationGrids", std::vector<std::string>{});

  for (const std::string& pseudo_atom_name : pseudoAtomStringGrids)
  {
    std::optional<std::size_t> index = findPseudoAtom(pseudo_atom_name);

    if (!index.has_value())
    {
      throw std::runtime_error(std::format("[ReadForceFieldSelfInteractions]: unknown pseudo atom {} in {}\n",
                                           pseudo_atom_name, parsed_data["UseInterpolationGrids"].dump()));
    }

    if (index.has_value())
    {
      gridPseudoAtomIndices.push_back(index.value());
    }
  }

  if (parsed_data.contains("SpacingVDWGrid"))
  {
    spacingVDWGrid = parsed_data.value("SpacingVDWGrid", parsed_data["SpacingVDWGrid"]);
  }

  if (parsed_data.contains("SpacingCoulombGrid"))
  {
    spacingCoulombGrid = parsed_data.value("SpacingCoulombGrid", parsed_data["SpacingCoulombGrid"]);
  }

  if (parsed_data.contains("NumberOfVDWGridPoints"))
  {
    numberOfVDWGridPoints = parseUint3("NumberOfVDWGridPoints", parsed_data["NumberOfVDWGridPoints"]);
  }
  if (parsed_data.contains("NumberOfCoulombGridPoints"))
  {
    numberOfCoulombGridPoints = parseUint3("NumberOfCoulombGridPoints", parsed_data["NumberOfCoulombGridPoints"]);
  }
  if (parsed_data.contains("NumberOfExternalFieldGridPoints"))
  {
    numberOfExternalFieldGridPoints = parseUint3("NumberOfExternalFieldGridPoints", parsed_data["NumberOfExternalFieldGridPoints"]);
  }

  if (parsed_data.contains("NumberOfGridTestPoints"))
  {
    numberOfGridTestPoints = parsed_data.value("NumberOfGridTestPoints", parsed_data["NumberOfGridTestPoints"]);
  }

  if (parsed_data.contains("ExternalFieldGridFileName"))
  {
    externalFieldGridFileName = parsed_data["ExternalFieldGridFileName"].get<std::string>();
  }

  if (parsed_data.contains("WriteExternalFieldInterpolationGrid"))
  {
    writeExternalFieldInterpolationGrid = parsed_data["WriteExternalFieldInterpolationGrid"].get<bool>();
  }

  if (parsed_data.contains("ExternalFieldPotentialEnergySurface"))
  {
    std::string potentialEnergySurfaceString = parsed_data["ExternalFieldPotentialEnergySurface"].get<std::string>();

    if (caseInSensStringCompare(potentialEnergySurfaceString, "GridFile"))
    {
      potentialEnergySurfaceType = PotentialEnergySurfaceType::GridFile;
    }
    if (caseInSensStringCompare(potentialEnergySurfaceString, "SecondOrderPolynomialTestFunction"))
    {
      potentialEnergySurfaceType = PotentialEnergySurfaceType::SecondOrderPolynomialTestFunction;
    }
    if (caseInSensStringCompare(potentialEnergySurfaceString, "ThirdOrderPolynomialTestFunction"))
    {
      potentialEnergySurfaceType = PotentialEnergySurfaceType::ThirdOrderPolynomialTestFunction;
    }
    if (caseInSensStringCompare(potentialEnergySurfaceString, "FourthOrderPolynomialTestFunction"))
    {
      potentialEnergySurfaceType = PotentialEnergySurfaceType::FourthOrderPolynomialTestFunction;
    }
    if (caseInSensStringCompare(potentialEnergySurfaceString, "FifthOrderPolynomialTestFunction"))
    {
      potentialEnergySurfaceType = PotentialEnergySurfaceType::FifthOrderPolynomialTestFunction;
    }
    if (caseInSensStringCompare(potentialEnergySurfaceString, "SixthOrderPolynomialTestFunction"))
    {
      potentialEnergySurfaceType = PotentialEnergySurfaceType::SixthOrderPolynomialTestFunction;
    }
    if (caseInSensStringCompare(potentialEnergySurfaceString, "ExponentialNonPolynomialTestFunction"))
    {
      potentialEnergySurfaceType = PotentialEnergySurfaceType::ExponentialNonPolynomialTestFunction;
    }
    if (caseInSensStringCompare(potentialEnergySurfaceString, "MullerBrown"))
    {
      potentialEnergySurfaceType = PotentialEnergySurfaceType::MullerBrown;
    }
    if (caseInSensStringCompare(potentialEnergySurfaceString, "Eckhardt"))
    {
      potentialEnergySurfaceType = PotentialEnergySurfaceType::Eckhardt;
    }
    if (caseInSensStringCompare(potentialEnergySurfaceString, "GonzalezSchlegel"))
    {
      potentialEnergySurfaceType = PotentialEnergySurfaceType::GonzalezSchlegel;
    }
  }

  if (parsed_data.contains("ExternalPotentialEnergySurfaceOrigin"))
  {
    potentialEnergySurfaceOrigin = parseDouble3("ExternalPotentialEnergySurfaceOrigin", parsed_data["ExternalPotentialEnergySurfaceOrigin"]);
  }

  if (parsed_data.contains("InterpolationScheme"))
  {
    std::size_t scheme = parsed_data.value("InterpolationScheme", parsed_data["InterpolationScheme"]);
    switch (scheme)
    {
      case 1:
        interpolationSchemeAuto = false;
        interpolationScheme = InterpolationScheme::Polynomial;
        break;
      case 3:
        interpolationSchemeAuto = false;
        interpolationScheme = InterpolationScheme::Tricubic;
        break;
      case 5:
        interpolationSchemeAuto = false;
        interpolationScheme = InterpolationScheme::Triquintic;
        break;
      default:
        throw std::runtime_error(std::format(
            "[ReadForceFieldSelfInteractions]: unknown grid interpolation scheme {} in {} (options: 3 or 5)\n", scheme,
            parsed_data["InterpolationScheme"].dump()));
        break;
    }
  }
}

void ForceField::applyMixingRule()
{
  switch (mixingRule)
  {
    case MixingRule::Lorentz_Berthelot:
      for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
      {
        for (std::size_t j = i + 1; j < numberOfPseudoAtoms; ++j)
        {
          if (data[i * numberOfPseudoAtoms + i].type == VDWParameters::Type::LennardJones &&
              data[j * numberOfPseudoAtoms + j].type == VDWParameters::Type::LennardJones)
          {
            double mix0 = std::sqrt(data[i * numberOfPseudoAtoms + i].parameters.x *
                                    data[j * numberOfPseudoAtoms + j].parameters.x);
            double mix1 =
                0.5 * (data[i * numberOfPseudoAtoms + i].parameters.y + data[j * numberOfPseudoAtoms + j].parameters.y);

            data[i * numberOfPseudoAtoms + j].parameters.x = mix0;
            data[i * numberOfPseudoAtoms + j].parameters.y = mix1;
            data[i * numberOfPseudoAtoms + j].type = VDWParameters::Type::LennardJones;
            data[j * numberOfPseudoAtoms + i].parameters.x = mix0;
            data[j * numberOfPseudoAtoms + i].parameters.y = mix1;
            data[j * numberOfPseudoAtoms + i].type = VDWParameters::Type::LennardJones;
          }
        }
      }
      break;
    case MixingRule::Jorgensen:
      for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
      {
        for (std::size_t j = i + 1; j < numberOfPseudoAtoms; ++j)
        {
          if (data[i * numberOfPseudoAtoms + i].type == VDWParameters::Type::LennardJones &&
              data[j * numberOfPseudoAtoms + j].type == VDWParameters::Type::LennardJones)
          {
            double mix0 = std::sqrt(data[i * numberOfPseudoAtoms + i].parameters.x *
                                    data[j * numberOfPseudoAtoms + j].parameters.x);
            double mix1 = std::sqrt(data[i * numberOfPseudoAtoms + i].parameters.y *
                                    data[j * numberOfPseudoAtoms + j].parameters.y);

            data[i * numberOfPseudoAtoms + j].parameters.x = mix0;
            data[i * numberOfPseudoAtoms + j].parameters.y = mix1;
            data[i * numberOfPseudoAtoms + j].type = VDWParameters::Type::LennardJones;

            data[j * numberOfPseudoAtoms + i].parameters.x = mix0;
            data[j * numberOfPseudoAtoms + i].parameters.y = mix1;
            data[j * numberOfPseudoAtoms + i].type = VDWParameters::Type::LennardJones;
          }
        }
      }
      break;
  }

  // set all interactions without interaction energy or length to none interactions
  for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (std::size_t j = i; j < numberOfPseudoAtoms; ++j)
    {
      if ((data[i * numberOfPseudoAtoms + i].type == VDWParameters::Type::LennardJones &&
           (data[i * numberOfPseudoAtoms + i].parameters.x == 0.0 ||
            data[i * numberOfPseudoAtoms + i].parameters.y == 0.0)) ||
          (data[j * numberOfPseudoAtoms + j].type == VDWParameters::Type::LennardJones &&
           (data[j * numberOfPseudoAtoms + j].parameters.x == 0.0 ||
            data[j * numberOfPseudoAtoms + j].parameters.y == 0.0)))
      {
        data[i * numberOfPseudoAtoms + j].type = VDWParameters::Type::None;
        data[j * numberOfPseudoAtoms + i].type = VDWParameters::Type::None;
      }
    }
  }
}

double ForceField::cutOffVDW(std::size_t i, std::size_t j) const
{
  if (pseudoAtoms[i].framework || pseudoAtoms[j].framework)
  {
    return cutOffFrameworkVDW;
  }

  return cutOffMoleculeVDW;
}

void ForceField::preComputePotentialShift()
{
  for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (std::size_t j = 0; j < numberOfPseudoAtoms; ++j)
    {
      if (shiftPotentials[i * numberOfPseudoAtoms + j])
      {
        double cut_off_vdw = cutOffVDW(i, j);
        data[i * numberOfPseudoAtoms + j].computeShiftAtCutOff(cut_off_vdw);
      }
    }
  }
}

void ForceField::preComputeTailCorrection()
{
  for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (std::size_t j = 0; j < numberOfPseudoAtoms; ++j)
    {
      data[i * numberOfPseudoAtoms + j].tailCorrectionEnergy = 0.0;
      data[i * numberOfPseudoAtoms + j].tailCorrectionPressure = 0.0;

      if (tailCorrections[i * numberOfPseudoAtoms + j])
      {
        double cut_off_vdw = cutOffVDW(i, j);
        VDWParameters::Type potentialType = data[i * numberOfPseudoAtoms + j].type;
        double4 parameters = data[i * numberOfPseudoAtoms + j].parameters;
        data[i * numberOfPseudoAtoms + j].tailCorrectionEnergy = Potentials::potentialCorrectionVDW(potentialType, parameters, cut_off_vdw);
        data[i * numberOfPseudoAtoms + j].tailCorrectionPressure = Potentials::potentialCorrectionPressure(potentialType, parameters, cut_off_vdw);
      }
    }
  }
}

std::optional<ForceField> ForceField::readForceField(std::optional<std::string> directoryName,
                                                     std::string forceFieldFileName) noexcept(false)
{
  // try to look in directory 'directoryName' if set, otherwise the local directory
  std::string filePath = directoryName.value_or(".") + "/" + forceFieldFileName;
  std::filesystem::path forceFieldPathfile = std::filesystem::path(filePath);
  if (!std::filesystem::exists(forceFieldPathfile))
  {
    // if not found, try the install directory and directory 'directoryName' in 'share/raspa3/forcefields'
    const char* env_p = std::getenv("RASPA_DIR");
    forceFieldPathfile = std::filesystem::path(std::string(env_p) + "/share/raspa3/forcefields/" + filePath);
    if (!std::filesystem::exists(forceFieldPathfile))
    {
      return std::nullopt;
    }
  }

  std::ifstream forceFieldStream{forceFieldPathfile};
  if (!forceFieldStream)
  {
    return std::nullopt;
  }
  return ForceField(forceFieldPathfile.string());
}

std::string ForceField::printPseudoAtomStatus() const
{
  std::ostringstream stream;

  std::print(stream, "Pseudo-atoms\n");
  std::print(stream, "===============================================================================\n\n");

  for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    std::print(stream, "{:3d} - {:8} framework-atom: {}\n", i, pseudoAtoms[i].name, pseudoAtoms[i].framework);
  }
  std::print(stream, "\n");
  for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    std::print(stream, "{:3d} - {:8} mass: {:8.5f}, charge: {:8.5f}\n", i, pseudoAtoms[i].name, pseudoAtoms[i].mass,
               pseudoAtoms[i].charge);
  }
  std::print(stream, "\n");
  for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    std::print(stream, "{:3d} - {:8} polarizability: {:8.5f}\n", i, pseudoAtoms[i].name, pseudoAtoms[i].polarizability);
  }
  std::print(stream, "\n");

  for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    if (!pseudoAtoms[i].source.empty())
    {
      std::print(stream, "{:3d} - {}\n", i, pseudoAtoms[i].source);
    }
  }
  std::print(stream, "\n");

  return stream.str();
}

std::vector<nlohmann::json> ForceField::jsonPseudoAtomStatus() const
{
  std::vector<nlohmann::json> jsonPseudoAtoms(numberOfPseudoAtoms);

  for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    jsonPseudoAtoms[i]["name"] = pseudoAtoms[i].name;
    jsonPseudoAtoms[i]["mass"] = pseudoAtoms[i].mass;
    jsonPseudoAtoms[i]["charge"] = pseudoAtoms[i].charge;

    if (!pseudoAtoms[i].source.empty())
    {
      jsonPseudoAtoms[i]["source"] = pseudoAtoms[i].source;
    }
  }
  return jsonPseudoAtoms;
}

std::string ForceField::printCutOffAutoStatus() const
{
  std::ostringstream stream;

  if (cutOffFrameworkVDWAutomatic)
  {
    std::print(stream, "Cutoff Framework-Molecule VDW: {:9.5f} [{}]\n", cutOffFrameworkVDW,
               Units::displayedUnitOfLengthString);
  }
  if (cutOffMoleculeVDWAutomatic)
  {
    std::print(stream, "Cutoff Molecule-Molecule VDW: {:9.5f} [{}]\n", cutOffFrameworkVDW,
               Units::displayedUnitOfLengthString);
  }
  if (cutOffCoulombAutomatic)
  {
    std::print(stream, "Cutoff Coulomb: {:9.5f} [{}]\n", cutOffCoulomb, Units::displayedUnitOfLengthString);

    std::print(stream, "Ewald alpha: {}\n", EwaldAlpha);
    std::print(stream, "Ewald k-vectors: {} {} {}\n", numberOfWaveVectors.x, numberOfWaveVectors.y,
               numberOfWaveVectors.z);
  }
  if (cutOffFrameworkVDWAutomatic || cutOffFrameworkVDWAutomatic || cutOffCoulombAutomatic)
  {
    std::print(stream, "\n");
  }

  return stream.str();
}

std::string ForceField::printForceFieldStatus() const
{
  std::ostringstream stream;

  std::print(stream, "Force field status\n");
  std::print(stream, "===============================================================================\n\n");

  if (cutOffFrameworkVDWAutomatic)
  {
    std::print(stream, "Cutoff Framework-Molecule VDW:  auto\n");
  }
  else
  {
    std::print(stream, "Cutoff Framework-Molecule VDW: {:9.5f} [{}]\n", cutOffFrameworkVDW,
               Units::displayedUnitOfLengthString);
  }
  if (cutOffMoleculeVDWAutomatic)
  {
    std::print(stream, "Cutoff Molecule-Molecule VDW:   auto\n");
  }
  else
  {
    std::print(stream, "Cutoff Molecule-Molecule VDW:  {:9.5f} [{}]\n", cutOffFrameworkVDW,
               Units::displayedUnitOfLengthString);
  }
  if (cutOffCoulombAutomatic)
  {
    std::print(stream, "Cutoff Coulomb:                 auto\n");
  }
  else
  {
    std::print(stream, "Cutoff Coulomb:                {:9.5f} [{}]\n\n", cutOffCoulomb,
               Units::displayedUnitOfLengthString);
  }

  std::print(stream, "Overlap-criteria VDW:          {: .6e} [{}]\n\n", energyOverlapCriteria,
             Units::displayedUnitOfEnergyString);

  for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (std::size_t j = i; j < numberOfPseudoAtoms; ++j)
    {
      switch (data[i * numberOfPseudoAtoms + j].type)
      {
        case VDWParameters::Type::LennardJones:
          std::print(stream, "{:8} - {:8} {} p₀{}: {:9.5f} [{}], p₁: {:8.5f} [{}]\n", pseudoAtoms[i].name,
                     pseudoAtoms[j].name, "Lennard-Jones", Units::displayedUnitOfEnergyConversionString,
                     Units::EnergyToKelvin * data[i * numberOfPseudoAtoms + j].parameters.x,
                     Units::displayedUnitOfEnergyString, data[i * numberOfPseudoAtoms + j].parameters.y,
                     Units::displayedUnitOfLengthString);
          std::print(stream, "{:33} shift: {:9.5f} [{}], tailcorrections: {}\n", std::string(""),
                     Units::EnergyToKelvin * data[i * numberOfPseudoAtoms + j].shift,
                     Units::displayedUnitOfEnergyString,
                     tailCorrections[i * numberOfPseudoAtoms + j] ? "true" : "false");
          break;
        case VDWParameters::Type::None:
          std::print(stream, "{:8} - {:8} None\n", pseudoAtoms[i].name, pseudoAtoms[j].name);
          break;
        default:
          break;
      }
    }
  }
  std::print(stream, "\n");

  if (!gridPseudoAtomIndices.empty())
  {
    if (numberOfVDWGridPoints.has_value())
    {
      std::print(stream, "Number of Van Der Waals grid points: {}x{}x{}\n", numberOfVDWGridPoints->x,
                 numberOfVDWGridPoints->y, numberOfVDWGridPoints->z);
    }
    else
    {
      std::print(stream, "Spacing of the Van Der Waals grid: {}\n", spacingVDWGrid);
    }
    if (numberOfCoulombGridPoints.has_value())
    {
      std::print(stream, "Number of Coulomb grid points: {}x{}x{}\n", numberOfCoulombGridPoints->x,
                 numberOfCoulombGridPoints->y, numberOfCoulombGridPoints->z);
    }
    else
    {
      std::print(stream, "Spacing of the Coulomb grid: {}\n", spacingCoulombGrid);
    }
    switch (interpolationScheme)
    {
      case InterpolationScheme::Polynomial:
        std::print(stream, "Interpolation-scheme: quintic\n");
        break;
      case InterpolationScheme::Tricubic:
        std::print(stream, "Interpolation-scheme: tricubic\n");
        break;
      case InterpolationScheme::Triquintic:
        std::print(stream, "Interpolation-scheme: triquintic\n");
        break;
    }
    std::print(stream, "\n");
  }

  if (automaticEwald)
  {
    std::print(stream, "Ewald precision: {}\n", EwaldPrecision);
    std::print(stream, "Ewald alpha: {}\n", EwaldAlpha);
    std::print(stream, "Ewald k-vectors: {} {} {}\n", numberOfWaveVectors.x, numberOfWaveVectors.y,
               numberOfWaveVectors.z);
  }
  else
  {
    std::print(stream, "Ewald alpha: {}\n", EwaldAlpha);
    std::print(stream, "Ewald k-vectors: {} {} {}\n", numberOfWaveVectors.x, numberOfWaveVectors.y,
               numberOfWaveVectors.z);
  }
  std::print(stream, "\n\n");

  return stream.str();
}

nlohmann::json ForceField::jsonForceFieldStatus() const
{
  nlohmann::json status;
  std::size_t n_interactions = static_cast<std::size_t>(static_cast<double>(numberOfPseudoAtoms) *
                                                        (static_cast<double>(numberOfPseudoAtoms) + 1.0) / 2.0);
  std::vector<nlohmann::json> interactions(n_interactions);

  std::size_t count = 0;
  for (std::size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (std::size_t j = i; j < numberOfPseudoAtoms; ++j)
    {
      switch (data[i * numberOfPseudoAtoms + j].type)
      {
        case VDWParameters::Type::LennardJones:
          interactions[count]["typeA"] = pseudoAtoms[i].name;
          interactions[count]["typeB"] = pseudoAtoms[j].name;
          interactions[count]["potential"] = "Lennard-Jones";
          interactions[count]["ε/kʙ [K]"] = data[i * numberOfPseudoAtoms + j].parameters.x * Units::EnergyToKelvin;
          interactions[count]["σ/kʙ [Å]"] = data[i * numberOfPseudoAtoms + j].parameters.y;
          interactions[count]["shift [K]"] = data[i * numberOfPseudoAtoms + j].shift * Units::EnergyToKelvin;
          interactions[count]["tailCorrections"] = tailCorrections[i * numberOfPseudoAtoms + j];
          break;
        default:
          break;
      }
      count++;
    }
  }
  status["interactions"] = interactions;
  if (automaticEwald)
  {
    status["Ewald"]["precision"] = EwaldPrecision;
  }
  status["Ewald"]["alpha"] = EwaldAlpha;
  status["Ewald"]["kVectors"] = {numberOfWaveVectors.x, numberOfWaveVectors.y, numberOfWaveVectors.z};

  return status;
}

std::optional<std::size_t> ForceField::findPseudoAtom(const std::string& name) const
{
  std::vector<PseudoAtom>::const_iterator match =
      std::find_if(pseudoAtoms.begin(), pseudoAtoms.end(), [&name](const PseudoAtom& x) { return x.name == name; });
  if (match != std::end(pseudoAtoms))
  {
    return static_cast<std::size_t>(std::distance(pseudoAtoms.begin(), match));
  }

  return std::nullopt;
}

std::optional<std::size_t> ForceField::findPseudoAtom(const std::vector<PseudoAtom> pseudoAtoms,
                                                      const std::string& name)
{
  std::vector<PseudoAtom>::const_iterator match =
      std::find_if(pseudoAtoms.begin(), pseudoAtoms.end(), [&name](const PseudoAtom& x) { return x.name == name; });
  if (match != std::end(pseudoAtoms))
  {
    return static_cast<std::size_t>(std::distance(pseudoAtoms.begin(), match));
  }

  return std::nullopt;
}

void ForceField::initializeEwaldParameters(const SimulationBox& simulationBox)
{
  double3 perpendicularWidths = simulationBox.perpendicularWidths();

  if (automaticEwald)
  {
    // compute the alpha-parameter and max k-vectors from the relative precision
    double eps = std::min(std::fabs(EwaldPrecision), 0.5);

    double tol = std::sqrt(std::abs(std::log(eps * cutOffCoulomb)));

    EwaldAlpha = std::sqrt(std::abs(std::log(eps * cutOffCoulomb * tol))) / cutOffCoulomb;
    double tol1 = std::sqrt(-std::log(eps * cutOffCoulomb * (2.0 * tol * EwaldAlpha) * (2.0 * tol * EwaldAlpha)));

    numberOfWaveVectors =
        int3(static_cast<std::int32_t>(std::rint(0.25 + perpendicularWidths.x * EwaldAlpha * tol1 / std::numbers::pi)),
             static_cast<std::int32_t>(std::rint(0.25 + perpendicularWidths.y * EwaldAlpha * tol1 / std::numbers::pi)),
             static_cast<std::int32_t>(std::rint(0.25 + perpendicularWidths.z * EwaldAlpha * tol1 / std::numbers::pi)));

    std::size_t maxNumberOfWaveVector =
        static_cast<std::size_t>(std::max({numberOfWaveVectors.x, numberOfWaveVectors.y, numberOfWaveVectors.z}));
    reciprocalIntegerCutOffSquared = maxNumberOfWaveVector * maxNumberOfWaveVector;
  }
}

void ForceField::initializeAutomaticCutOff(const SimulationBox& simulationBox)
{
  double3 perpendicularWidths = simulationBox.perpendicularWidths();

  double smallest_perpendicular_width = std::min({perpendicularWidths.x, perpendicularWidths.y, perpendicularWidths.z});

  if (cutOffFrameworkVDWAutomatic)
  {
    cutOffFrameworkVDW = 0.5 * smallest_perpendicular_width - std::numeric_limits<double>::epsilon();
  }

  if (cutOffMoleculeVDWAutomatic)
  {
    cutOffMoleculeVDW = 0.5 * smallest_perpendicular_width - std::numeric_limits<double>::epsilon();
  }

  if (cutOffCoulombAutomatic && automaticEwald)
  {
    cutOffCoulomb = 0.5 * smallest_perpendicular_width - std::numeric_limits<double>::epsilon();

    // compute the alpha-parameter and max k-vectors from the relative precision
    double eps = std::min(std::fabs(EwaldPrecision), 0.5);

    double tol = std::sqrt(std::abs(std::log(eps * cutOffCoulomb)));

    EwaldAlpha = std::sqrt(std::abs(std::log(eps * cutOffCoulomb * tol))) / cutOffCoulomb;

    double tol1 = std::sqrt(-std::log(eps * cutOffCoulomb * (2.0 * tol * EwaldAlpha) * (2.0 * tol * EwaldAlpha)));

    numberOfWaveVectors =
        int3(static_cast<std::int32_t>(std::rint(0.25 + perpendicularWidths.x * EwaldAlpha * tol1 / std::numbers::pi)),
             static_cast<std::int32_t>(std::rint(0.25 + perpendicularWidths.y * EwaldAlpha * tol1 / std::numbers::pi)),
             static_cast<std::int32_t>(std::rint(0.25 + perpendicularWidths.z * EwaldAlpha * tol1 / std::numbers::pi)));

    std::size_t maxNumberOfWaveVector =
        static_cast<std::size_t>(std::max({numberOfWaveVectors.x, numberOfWaveVectors.y, numberOfWaveVectors.z}));
    reciprocalIntegerCutOffSquared = maxNumberOfWaveVector * maxNumberOfWaveVector;
  }
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ForceField& f)
{
  archive << f.versionNumber;

  archive << f.data;
  archive << f.shiftPotentials;
  archive << f.tailCorrections;
  archive << f.cutOffFrameworkVDWAutomatic;
  archive << f.cutOffFrameworkVDW;
  archive << f.cutOffMoleculeVDWAutomatic;
  archive << f.cutOffMoleculeVDW;
  archive << f.cutOffCoulombAutomatic;
  archive << f.cutOffCoulomb;
  archive << f.dualCutOff;

  archive << f.numberOfPseudoAtoms;
  archive << f.pseudoAtoms;

  archive << f.chargeMethod;

  archive << f.EwaldPrecision;
  archive << f.EwaldAlpha;
  archive << f.numberOfWaveVectors;
  archive << f.reciprocalIntegerCutOffSquared;
  archive << f.reciprocalCutOffSquared;
  archive << f.automaticEwald;
  archive << f.useCharge;
  archive << f.omitEwaldFourier;

  archive << f.energyOverlapCriteria;

  archive << f.numberOfTrialDirections;
  archive << f.numberOfTorsionTrialDirections;
  archive << f.numberOfFirstBeadPositions;
  archive << f.numberOfTrialMovesPerOpenBead;
  archive << f.minimumRosenbluthFactor;

  archive << f.useDualCutOff;
  archive << f.omitInterInteractions;

  archive << f.computePolarization;
  archive << f.omitInterPolarization;

  archive << f.potentialEnergySurfaceType;
  archive << f.potentialEnergySurfaceOrigin;

  archive << f.gridPseudoAtomIndices;
  archive << f.spacingVDWGrid;
  archive << f.spacingCoulombGrid;
  archive << f.numberOfVDWGridPoints;
  archive << f.numberOfCoulombGridPoints;
  archive << f.numberOfGridTestPoints;
  archive << f.interpolationScheme;
  archive << f.writeFrameworkInterpolationGrids;

  archive << f.useExternalFieldGrid;
  archive << f.externalFieldGridFileName;
  archive << f.numberOfExternalFieldGridPoints;
  archive << f.writeExternalFieldInterpolationGrid;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ForceField& f)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > f.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'ForceField' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> f.data;
  archive >> f.shiftPotentials;
  archive >> f.tailCorrections;
  archive >> f.cutOffFrameworkVDWAutomatic;
  archive >> f.cutOffFrameworkVDW;
  archive >> f.cutOffMoleculeVDWAutomatic;
  archive >> f.cutOffMoleculeVDW;
  archive >> f.cutOffCoulombAutomatic;
  archive >> f.cutOffCoulomb;
  archive >> f.dualCutOff;

  archive >> f.numberOfPseudoAtoms;
  archive >> f.pseudoAtoms;

  archive >> f.chargeMethod;

  archive >> f.EwaldPrecision;
  archive >> f.EwaldAlpha;
  archive >> f.numberOfWaveVectors;
  archive >> f.reciprocalIntegerCutOffSquared;
  archive >> f.reciprocalCutOffSquared;
  archive >> f.automaticEwald;
  archive >> f.useCharge;
  archive >> f.omitEwaldFourier;

  archive >> f.energyOverlapCriteria;

  archive >> f.numberOfTrialDirections;
  archive >> f.numberOfTorsionTrialDirections;
  archive >> f.numberOfFirstBeadPositions;
  archive >> f.numberOfTrialMovesPerOpenBead;
  archive >> f.minimumRosenbluthFactor;

  archive >> f.useDualCutOff;
  archive >> f.omitInterInteractions;

  archive >> f.computePolarization;
  archive >> f.omitInterPolarization;

  archive >> f.potentialEnergySurfaceType;
  archive >> f.potentialEnergySurfaceOrigin;

  archive >> f.gridPseudoAtomIndices;
  archive >> f.spacingVDWGrid;
  archive >> f.spacingCoulombGrid;
  archive >> f.numberOfVDWGridPoints;
  archive >> f.numberOfCoulombGridPoints;
  archive >> f.numberOfGridTestPoints;
  archive >> f.interpolationScheme;
  archive >> f.writeFrameworkInterpolationGrids;

  archive >> f.useExternalFieldGrid;
  archive >> f.externalFieldGridFileName;
  archive >> f.numberOfExternalFieldGridPoints;
  archive >> f.writeExternalFieldInterpolationGrid;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("ForceField: Error in binary restart\n"));
  }
#endif

  return archive;
}

bool ForceField::operator==(const ForceField& other) const
{
  // first the cheap ones
  if (cutOffFrameworkVDW != other.cutOffFrameworkVDW || cutOffMoleculeVDW != other.cutOffMoleculeVDW ||
      cutOffCoulomb != other.cutOffCoulomb || dualCutOff != other.dualCutOff ||
      numberOfPseudoAtoms != other.numberOfPseudoAtoms || energyOverlapCriteria != other.energyOverlapCriteria ||
      EwaldPrecision != other.EwaldPrecision || EwaldAlpha != other.EwaldAlpha ||
      numberOfWaveVectors != other.numberOfWaveVectors || automaticEwald != other.automaticEwald ||
      useCharge != other.useCharge || omitEwaldFourier != other.omitEwaldFourier ||
      minimumRosenbluthFactor != other.minimumRosenbluthFactor ||
      energyOverlapCriteria != other.energyOverlapCriteria || useDualCutOff != other.useDualCutOff ||
      chargeMethod != other.chargeMethod)
  {
    return false;
  }
  for (std::size_t idx = 0; idx < shiftPotentials.size(); idx++)
  {
    if (shiftPotentials[idx] != other.shiftPotentials[idx] || tailCorrections[idx] != other.tailCorrections[idx] ||
        data[idx] != other.data[idx])
    {
      return false;
    }
  }
  for (auto it1 = pseudoAtoms.begin(), it2 = pseudoAtoms.begin(); it1 != pseudoAtoms.end(); ++it1, ++it2)
  {
    if (it1 != it2)
    {
      return false;
    }
  }
  return true;
}

const std::set<std::string, ForceField::InsensitiveCompare> ForceField::options = {"MixingRule",
                                                                                   "TruncationMethod",
                                                                                   "TailCorrections",
                                                                                   "TruncatedPotentialPairs",
                                                                                   "ShiftedPotentialPairs",
                                                                                   "TailCorrectionPairs",
                                                                                   "NoTailCorrectionPairs",
                                                                                   "CutOff",
                                                                                   "CutOffVDW",
                                                                                   "CutOffFrameworkVDW",
                                                                                   "CutOffMoleculeVDW",
                                                                                   "CutOffCoulomb",
                                                                                   "ChargeMethod",
                                                                                   "PseudoAtoms",
                                                                                   "SelfInteractions",
                                                                                   "BinaryInteractions",
                                                                                   "EwaldPrecision",
                                                                                   "EwaldParameters",
                                                                                   "ReciprocalCutOff",
                                                                                   "ReciprocalIntegerCutOff",
                                                                                   "UseInterpolationGrids",
                                                                                   "SpacingVDWGrid",
                                                                                   "SpacingCoulombGrid",
                                                                                   "NumberOfVDWGridPoints",
                                                                                   "NumberOfGridTestPoints",
                                                                                   "ExternalFieldGridFileName",
                                                                                   "NumberOfExternalFieldGridPoints",
                                                                                   "UseExternalFieldGrid",
                                                                                   "ExternalFieldPotentialEnergySurface",
                                                                                   "ExternalPotentialEnergySurfaceOrigin",
                                                                                   "WriteExternalFieldInterpolationGrid",
                                                                                   "InterpolationScheme"};

void ForceField::validateInput(const nlohmann::basic_json<nlohmann::raspa_map>& parsed_data)
{
  for (auto& [key, _] : parsed_data.items())
  {
    if (!options.contains(key))
    {
      throw std::runtime_error(std::format("[ForceField Error] : Unknown input '{}'\n", key));
    }
  }
}

ForceField ForceField::makeZeoliteForceField(double rc, bool shifted, bool tailCorrections, bool useEwald)
{
  return ForceField(
      {{"-", false, 0.0, 0.0, 0.0, 0, false},
       {"Si", true, 28.0855, 2.05, 0.0, 14, false},
       {"Al", true, 26.982, 2.05, 0.0, 13, false},
       {"O", true, 15.999, -1.025, 0.0, 8, false},
       {"Na+", false, 12.0, 1.0, 0.0, 6, false},
       {"Cl-", false, 15.9994, -1.0, 0.0, 8, false},
       {"CH4", false, 16.04246, 0.0, 0.0, 6, false},
       {"C_co2", false, 12.0, 0.6512, 0.2, 6, false},
       {"O_co2", false, 15.9994, -0.3256, 0.1, 8, false},
       {"Ow", false, 15.9996, 0.0, 0.0, 8, false},
       {"Hw", false, 1.0008, 0.241, 0.0, 1, false},
       {"Lw", false, 0.0, -0.241, 0.0, 0, false},
       {"probe-He", false, 4.002602, 0.0, 0.0, 2, false},
       {"probe-Ar", false, 39.948, 0.0, 0.0, 18, false},
       {"probe-CH4", false, 16.04246, 0.0, 0.0, 6, false},
       {"probe-N2", false, 14.00674, 0.0, 0.0, 6, false}},
      {{1.0, 1.0},
       {22.0, 2.30},
       {22.0, 2.30},
       {53.0, 3.30},
       {15.0966, 2.65755},
       {142.562, 3.51932},
       {158.5, 3.72},
       {29.933, 2.745},
       {85.671, 3.017},
       {89.633, 3.097},
       {0.0, 1.0},
       {0.0, 1.0},
       {10.9, 2.64},
       {124.070, 3.38},
       {158.5, 3.72},
       {91.5, 3.681}},
       ForceField::MixingRule::Lorentz_Berthelot, rc, rc, rc, shifted, tailCorrections, useEwald);
}

ForceField ForceField::makeMetalOrganicFrameworkForceField(double rc, bool shifted, bool tailCorrections, bool useEwald)
{
  return ForceField(
      {{"-",  false,  0.0,        0.0, 0.0,   0, false},
       {"H",  true,   1.00794,    0.0, 0.0,   1, false}, //   1
       {"He", true,   4.002602,   0.0, 0.0,   2, false}, //   2
       {"Li", true,   6.9421,     0.0, 0.0,   3, false}, //   3
       {"Be", true,   9.012182,   0.0, 0.0,   4, false}, //   4
       {"B",  true,  10.881,      0.0, 0.0,   5, false}, //   5
       {"C",  true,  12.0107,     0.0, 0.0,   6, false}, //   6
       {"N",  true,  14.0067,     0.0, 0.0,   7, false}, //   7
       {"O",  true,  15.9994,     0.0, 0.0,   8, false}, //   8
       {"F",  true,  18.9984032,  0.0, 0.0,   9, false}, //   9
       {"Ne", true,  20.1797,     0.0, 0.0,  10, false}, //  10
       {"Na", true,  22.98976928, 0.0, 0.0,  11, false}, //  11
       {"Mg", true,  24.305,      0.0, 0.0,  12, false}, //  12
       {"Al", true,  26.9815386,  0.0, 0.0,  13, false}, //  13
       {"Si", true,  28.0855,     0.0, 0.0,  14, false}, //  14
       {"P",  true,  30.973762,   0.0, 0.0,  15, false}, //  15
       {"S",  true,  32.065,      0.0, 0.0,  16, false}, //  16
       {"Cl", true,  35.453,      0.0, 0.0,  17, false}, //  17
       {"Ar", true,  39.948,      0.0, 0.0,  18, false}, //  18
       {"K",  true,  39.0983,     0.0, 0.0,  19, false}, //  19
       {"Ca", true,  40.078,      0.0, 0.0,  20, false}, //  20
       {"Sc", true,  44.955912,   0.0, 0.0,  21, false}, //  21
       {"Ti", true,  47.867,      0.0, 0.0,  22, false}, //  22
       {"V",  true,  50.9415,     0.0, 0.0,  23, false}, //  23
       {"Cr", true,  51.9961,     0.0, 0.0,  24, false}, //  24
       {"Mn", true,  54.939045,   0.0, 0.0,  25, false}, //  25
       {"Fe", true,  55.845,      0.0, 0.0,  26, false}, //  26
       {"Co", true,  58.933195,   0.0, 0.0,  27, false}, //  27
       {"Ni", true,  58.6934,     0.0, 0.0,  28, false}, //  28
       {"Cu", true,  63.546,      0.0, 0.0,  29, false}, //  29
       {"Zn", true,  65.38,       0.0, 0.0,  30, false}, //  30
       {"Ga", true,  69.723,      0.0, 0.0,  31, false}, //  31
       {"Ge", true,  72.64,       0.0, 0.0,  32, false}, //  32
       {"As", true,  74.9216,     0.0, 0.0,  33, false}, //  33
       {"Se", true,  78.96,       0.0, 0.0,  34, false}, //  34
       {"Br", true,  79.904,      0.0, 0.0,  35, false}, //  35
       {"Kr", true,  83.798,      0.0, 0.0,  36, false}, //  36
       {"Rb", true,  85.4678,     0.0, 0.0,  37, false}, //  37
       {"Sr", true,  87.62,       0.0, 0.0,  38, false}, //  38
       {"Y",  true,  88.90585,    0.0, 0.0,  39, false}, //  39
       {"Zr", true,  91.224,      0.0, 0.0,  40, false}, //  40
       {"Nb", true,  92.90638,    0.0, 0.0,  41, false}, //  41
       {"Mo", true,  95.96,       0.0, 0.0,  42, false}, //  42
       {"Tc", true,  98.0,        0.0, 0.0,  43, false}, //  43
       {"Ru", true, 101.07,       0.0, 0.0,  44, false}, //  44
       {"Rh", true, 102.59055,    0.0, 0.0,  45, false}, //  45
       {"Pd", true, 106.42,       0.0, 0.0,  46, false}, //  46
       {"Ag", true, 107.8682,     0.0, 0.0,  47, false}, //  47
       {"Cd", true, 112.411,      0.0, 0.0,  48, false}, //  48
       {"In", true, 114.818,      0.0, 0.0,  49, false}, //  49
       {"Sn", true, 118.71,       0.0, 0.0,  50, false}, //  50
       {"Sb", true, 121.76,       0.0, 0.0,  51, false}, //  51
       {"Te", true, 127.6,        0.0, 0.0,  52, false}, //  52
       {"I",  true, 126.90447,    0.0, 0.0,  53, false}, //  53
       {"Xe", true, 131.293,      0.0, 0.0,  54, false}, //  54
       {"Cs", true, 132.9054519,  0.0, 0.0,  55, false}, //  55
       {"Ba", true, 137.327,      0.0, 0.0,  56, false}, //  56
       {"La", true, 138.90547,    0.0, 0.0,  57, false}, //  57
       {"Ce", true, 140.116,      0.0, 0.0,  58, false}, //  58
       {"Pr", true, 140.90765,    0.0, 0.0,  59, false}, //  59
       {"Nd", true, 144.242,      0.0, 0.0,  60, false}, //  60
       {"Pm", true, 145.0,        0.0, 0.0,  61, false}, //  61
       {"Sm", true, 150.36,       0.0, 0.0,  62, false}, //  62
       {"Eu", true, 151.964,      0.0, 0.0,  63, false}, //  63
       {"Gd", true, 157.25,       0.0, 0.0,  64, false}, //  64
       {"Tb", true, 158.92535,    0.0, 0.0,  65, false}, //  65
       {"Dy", true, 162.5,        0.0, 0.0,  66, false}, //  66
       {"Ho", true, 164.93032,    0.0, 0.0,  67, false}, //  67
       {"Er", true, 167.259,      0.0, 0.0,  68, false}, //  68
       {"Tm", true, 168.93421,    0.0, 0.0,  69, false}, //  69
       {"Yb", true, 173.054,      0.0, 0.0,  70, false}, //  70
       {"Lu", true, 174.9668,     0.0, 0.0,  71, false}, //  71
       {"Hf", true, 178.49,       0.0, 0.0,  72, false}, //  72
       {"Ta", true, 180.94788,    0.0, 0.0,  73, false}, //  73
       {"W",  true, 183.84,       0.0, 0.0,  74, false}, //  74
       {"Re", true, 186.207,      0.0, 0.0,  75, false}, //  75
       {"Os", true, 190.23,       0.0, 0.0,  76, false}, //  76
       {"Ir", true, 192.217,      0.0, 0.0,  77, false}, //  77
       {"Pt", true, 195.084,      0.0, 0.0,  78, false}, //  78
       {"Au", true, 196.966569,   0.0, 0.0,  79, false}, //  79
       {"Hg", true, 200.59,       0.0, 0.0,  80, false}, //  80
       {"Tl", true, 204.3833,     0.0, 0.0,  81, false}, //  81
       {"Pb", true, 207.2,        0.0, 0.0,  82, false}, //  82
       {"Bi", true, 208.9804,     0.0, 0.0,  83, false}, //  83
       {"Po", true, 210.0,        0.0, 0.0,  84, false}, //  84
       {"At", true, 210.8,        0.0, 0.0,  85, false}, //  85
       {"Rn", true, 222.0,        0.0, 0.0,  86, false}, //  86
       {"Fr", true, 223.0,        0.0, 0.0,  87, false}, //  87
       {"Ra", true, 226.0,        0.0, 0.0,  88, false}, //  88
       {"Ac", true, 227.0,        0.0, 0.0,  89, false}, //  89
       {"Th", true, 232.03806,    0.0, 0.0,  90, false}, //  90
       {"Pa", true, 231.03588,    0.0, 0.0,  91, false}, //  91
       {"U",  true, 238.02891,    0.0, 0.0,  92, false}, //  92
       {"Np", true, 237.0,        0.0, 0.0,  93, false}, //  93
       {"Pu", true, 244.0,        0.0, 0.0,  94, false}, //  94
       {"Am", true, 243.0,        0.0, 0.0,  95, false}, //  95
       {"Cm", true, 247.0,        0.0, 0.0,  96, false}, //  96
       {"Bk", true, 247.0,        0.0, 0.0,  97, false}, //  97
       {"Cf", true, 251.0,        0.0, 0.0,  98, false}, //  98
       {"Es", true, 252.0,        0.0, 0.0,  99, false}, //  99
       {"Fm", true, 257.0,        0.0, 0.0, 100, false}, // 100
       {"Md", true, 258.0,        0.0, 0.0, 101, false}, // 101
       {"No", true, 259.0,        0.0, 0.0, 102, false}, // 102
       {"Lr", true, 262.0,        0.0, 0.0, 103, false}, // 103
       {"Rf", true, 261.0,        0.0, 0.0, 104, false}, // 104
       {"Db", true, 268.0,        0.0, 0.0, 105, false}, // 105
       {"Sg", true, 269.0,        0.0, 0.0, 106, false}, // 106
       {"Bh", true, 270.0,        0.0, 0.0, 107, false}, // 107
       {"Hs", true, 269.0,        0.0, 0.0, 108, false}, // 108
       {"Mt", true, 278.0,        0.0, 0.0, 109, false}, // 109
       {"Ds", true, 281.0,        0.0, 0.0, 110, false}, // 110
       {"Rg", true, 281.0,        0.0, 0.0, 111, false}, // 111
       {"Cn", true, 285.0,        0.0, 0.0, 112, false}, // 112
       {"Nh", true, 286.0,        0.0, 0.0, 113, false}, // 113
       {"Fl", true, 289.0,        0.0, 0.0, 114, false}, // 114
       {"Mc", true, 288.0,        0.0, 0.0, 115, false}, // 115
       {"Lv", true, 293.0,        0.0, 0.0, 116, false}, // 116
       {"Ts", true, 294.0,        0.0, 0.0, 117, false}, // 117
       {"Og", true, 294.0,        0.0, 0.0, 118, false}, // 118
       {"probe-He",  false,  4.002602, 0.0, 0.0,  2, false},
       {"probe-Ar",  false, 39.948,    0.0, 0.0, 18, false},
       {"probe-CH4", false, 16.04246,  0.0, 0.0,  6, false},
       {"probe-N2",  false, 14.00674,  0.0, 0.0,  6, false}},
      {{1.0, 1.0},
       {7.64893, 2.84642}, //   1 "H"   DREIDING S.L. Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909.
       {10.9, 2.64},       //   2 "He"
       {1.0, 1.0}, //   3 "Li"
       {42.7736, 2.44552}, //   4 "Be"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {47.8058, 3.58141}, //   5 "B"   DREIDING S.L. Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909.
       {47.8562, 3.47299}, //   6 "C"   DREIDING S.L. Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909.
       {38.9492, 3.26256}, //   7 "N"   DREIDING S.L. Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909.
       {48.1581, 3.03315}, //   8 "O"   DREIDING S.L. Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909.
       {36.4834, 3.0932},  //   9 "F"   DREIDING S.L. Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909.
       {21.1352, 2.88918}, //  10 "Ne"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {1.0, 1.0}, //  11 "Na"
       {55.8574, 2.69141}, //  12 "Mg"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {155.998, 3.91105}, //  13 "Al"  DREIDING S.L. Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909.
       {155.998, 3.80414}, //  14 "Si"  DREIDING S.L. Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909.
       {161.03, 3.69723},  //  15 "P"   DREIDING S.L. Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909.
       {173.107, 3.59032}, //  16 "S"   DREIDING S.L. Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909.
       {142.562, 3.51932}, //  17 "Cl"  DREIDING S.L. Mayo et al., J. Phys. Chem. 1990, 94, 8897-8909.
       {1.0, 1.0}, //  18 "Ar"
       {1.0, 1.0}, //  19 "K"
       {1.0, 1.0}, //  20 "Ca"
       {9.56117, 2.93551}, //  21 "Sc"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {8.55473, 2.8286},  //  22 "Ti"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {8.05151, 2.80099}, //  23 "V"   UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {7.54829, 2.69319}, //  24 "Cr"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {6.54185, 2.63795}, //  25 "Mn"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {6.54185, 2.5943},  //  26 "Fe"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {7.04507, 2.55866}, //  27 "Co"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {7.54829, 2.52481}, //  28 "Ni"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {2.5161, 3.11369},  //  29 "Cu"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {27.677074, 4.04},  //  30 "Zn"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {208.836, 3.90481}, //  31 "Ga"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {1.0, 1.0}, //  32 "Ge"
       {1.0, 1.0}, //  33 "As"
       {1.0, 1.0}, //  34 "Se"
       {186.191, 3.51905}, //  35 "Br"
       {1.0, 1.0}, //  36 "Kr"
       {1.0, 1.0}, //  37 "Rb"
       {1.0, 1.0}, //  38 "Sr"
       {1.0, 1.0}, //  39 "Y"
       {34.7221, 2.78317}, //  40 "Zr"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {1.0, 1.0}, //  41 "Nb"
       {1.0, 1.0}, //  42 "Mo"
       {1.0, 1.0}, //  43 "Tc"
       {1.0, 1.0}, //  44 "Ru"
       {1.0, 1.0}, //  45 "Rh"
       {1.0, 1.0}, //  46 "Pd"
       {18.1159, 2.80455}, //  47 "Ag"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {114.734, 2.53728}, //  48 "Cd"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {301.428, 3.97608}, //  49 "In"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {1.0, 1.0}, //  50 "Sn"
       {225.946, 3.93777}, //  51 "Sb"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {200.281, 3.98232}, //  52 "Te"  UFF A.K. Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024-10035.
       {1.0, 1.0}, //  53 "I"
       {1.0, 1.0}, //  54 "Xe"
       {1.0, 1.0}, //  55 "Cs"
       {1.0, 1.0}, //  56 "Ba"
       {1.0, 1.0}, //  57 "La"
       {1.0, 1.0}, //  58 "Ce"
       {1.0, 1.0}, //  59 "Pr"
       {1.0, 1.0}, //  60 "Nd"
       {1.0, 1.0}, //  61 "Pm"
       {1.0, 1.0}, //  62 "Sm"
       {1.0, 1.0}, //  63 "Eu"
       {1.0, 1.0}, //  64 "Gd"
       {1.0, 1.0}, //  65 "Tb"
       {1.0, 1.0}, //  66 "Dy"
       {1.0, 1.0}, //  67 "Ho"
       {1.0, 1.0}, //  68 "Er"
       {1.0, 1.0}, //  69 "Tm"
       {1.0, 1.0}, //  70 "Yb"
       {1.0, 1.0}, //  71 "Lu"
       {1.0, 1.0}, //  72 "Hf"
       {1.0, 1.0}, //  73 "Ta"
       {1.0, 1.0}, //  74 "W"
       {1.0, 1.0}, //  75 "Re"
       {1.0, 1.0}, //  76 "Os"
       {1.0, 1.0}, //  77 "Ir"
       {1.0, 1.0}, //  78 "Pt"
       {1.0, 1.0}, //  79 "Au"
       {1.0, 1.0}, //  80 "Hg"
       {1.0, 1.0}, //  81 "Tl"
       {1.0, 1.0}, //  82 "Pb"
       {1.0, 1.0}, //  83 "Bi"
       {1.0, 1.0}, //  84 "Po"
       {1.0, 1.0}, //  85 "At"
       {1.0, 1.0}, //  86 "Rn"
       {1.0, 1.0}, //  87 "Fr"
       {1.0, 1.0}, //  88 "Ra"
       {1.0, 1.0}, //  89 "Ac"
       {1.0, 1.0}, //  90 "Th"
       {1.0, 1.0}, //  91 "Pa"
       {1.0, 1.0}, //  92 "U"
       {1.0, 1.0}, //  93 "Np"
       {1.0, 1.0}, //  94 "Pu"
       {1.0, 1.0}, //  95 "Am"
       {1.0, 1.0}, //  96 "Cm"
       {1.0, 1.0}, //  97 "Bk"
       {1.0, 1.0}, //  98 "Cf"
       {1.0, 1.0}, //  99 "Es"
       {1.0, 1.0}, // 100 "Fm"
       {1.0, 1.0}, // 101 "Md"
       {1.0, 1.0}, // 102 "No"
       {1.0, 1.0}, // 103 "Lr"
       {1.0, 1.0}, // 104 "Rf"
       {1.0, 1.0}, // 105 "Db"
       {1.0, 1.0}, // 106 "Sg"
       {1.0, 1.0}, // 107 "Bh"
       {1.0, 1.0}, // 108 "Hs"
       {1.0, 1.0}, // 109 "Mt"
       {1.0, 1.0}, // 110 "Ds"
       {1.0, 1.0}, // 111 "Rg"
       {1.0, 1.0}, // 112 "Cn"
       {1.0, 1.0}, // 113 "Nh"
       {1.0, 1.0}, // 114 "Fl"
       {1.0, 1.0}, // 115 "Mc"
       {1.0, 1.0}, // 116 "Lv"
       {1.0, 1.0}, // 117 "Ts"
       {1.0, 1.0}, // 118 "Og"
       {10.9, 2.64},    // probe-He  J.O. Hirschfelder et al., Molecular Theory of Gases and Liquids, Wiley, New York, 1954, p. 1114.
       {124.070, 3.38}, // probe-Ar
       {158.5, 3.72},   // probe-CH4
       {91.5, 3.681}},  // probe-N2
       ForceField::MixingRule::Lorentz_Berthelot, rc, rc, rc, shifted, tailCorrections, useEwald);
}

ForceField ForceField::makeZeoPlusPlusForceField(double rc, bool shifted, bool tailCorrections, bool useEwald)
{
  return ForceField({{"-",  false, 0.0, 0.0, 0.0, 0, false},
                     {"H",  false, 1.0, 0.0, 0.0, 1, false},
                     {"D",  false, 1.0, 0.0, 0.0, 1, false},
                     {"He", false, 1.0, 0.0, 0.0, 1, false},
                     {"Li", false, 1.0, 0.0, 0.0, 1, false},
                     {"Be", false, 1.0, 0.0, 0.0, 1, false},
                     {"B",  false, 1.0, 0.0, 0.0, 1, false},
                     {"C",  false, 1.0, 0.0, 0.0, 1, false},
                     {"N",  false, 1.0, 0.0, 0.0, 1, false},
                     {"O",  false, 15.999, 0.0, 0.0, 1, false},
                     {"F",  false, 1.0, 0.0, 0.0, 1, false},
                     {"Ne", false, 1.0, 0.0, 0.0, 1, false},
                     {"Na", false, 1.0, 0.0, 0.0, 1, false},
                     {"Mg", false, 1.0, 0.0, 0.0, 1, false},
                     {"Al", false, 1.0, 0.0, 0.0, 1, false},
                     {"Si", false, 28.0855, 0.0, 0.0, 1, false},
                     {"P",  false, 1.0, 0.0, 0.0, 1, false},
                     {"S",  false, 1.0, 0.0, 0.0, 1, false},
                     {"Cl", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ar", false, 1.0, 0.0, 0.0, 1, false},
                     {"K",  false, 1.0, 0.0, 0.0, 1, false},
                     {"Ca", false, 1.0, 0.0, 0.0, 1, false},
                     {"Sc", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ti", false, 1.0, 0.0, 0.0, 1, false},
                     {"V",  false, 1.0, 0.0, 0.0, 1, false},
                     {"Cr", false, 1.0, 0.0, 0.0, 1, false},
                     {"Mn", false, 1.0, 0.0, 0.0, 1, false},
                     {"Fe", false, 1.0, 0.0, 0.0, 1, false},
                     {"Co", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ni", false, 1.0, 0.0, 0.0, 1, false},
                     {"Cu", false, 1.0, 0.0, 0.0, 1, false},
                     {"Zn", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ga", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ge", false, 1.0, 0.0, 0.0, 1, false},
                     {"As", false, 1.0, 0.0, 0.0, 1, false},
                     {"Se", false, 1.0, 0.0, 0.0, 1, false},
                     {"Br", false, 1.0, 0.0, 0.0, 1, false},
                     {"Kr", false, 1.0, 0.0, 0.0, 1, false},
                     {"Rb", false, 1.0, 0.0, 0.0, 1, false},
                     {"Sr", false, 1.0, 0.0, 0.0, 1, false},
                     {"Y",  false, 1.0, 0.0, 0.0, 1, false},
                     {"Zr", false, 1.0, 0.0, 0.0, 1, false},
                     {"Nb", false, 1.0, 0.0, 0.0, 1, false},
                     {"Mo", false, 1.0, 0.0, 0.0, 1, false},
                     {"Tc", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ru", false, 1.0, 0.0, 0.0, 1, false},
                     {"Rh", false, 1.0, 0.0, 0.0, 1, false},
                     {"Pd", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ag", false, 1.0, 0.0, 0.0, 1, false},
                     {"Cd", false, 1.0, 0.0, 0.0, 1, false},
                     {"In", false, 1.0, 0.0, 0.0, 1, false},
                     {"Sn", false, 1.0, 0.0, 0.0, 1, false},
                     {"Sb", false, 1.0, 0.0, 0.0, 1, false},
                     {"Te", false, 1.0, 0.0, 0.0, 1, false},
                     {"I",  false, 1.0, 0.0, 0.0, 1, false},
                     {"Xe", false, 1.0, 0.0, 0.0, 1, false},
                     {"Cs", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ba", false, 1.0, 0.0, 0.0, 1, false},
                     {"La", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ce", false, 1.0, 0.0, 0.0, 1, false},
                     {"Pr", false, 1.0, 0.0, 0.0, 1, false},
                     {"Nd", false, 1.0, 0.0, 0.0, 1, false},
                     {"Pm", false, 1.0, 0.0, 0.0, 1, false},
                     {"Sm", false, 1.0, 0.0, 0.0, 1, false},
                     {"Eu", false, 1.0, 0.0, 0.0, 1, false},
                     {"Gd", false, 1.0, 0.0, 0.0, 1, false},
                     {"Tb", false, 1.0, 0.0, 0.0, 1, false},
                     {"Dy", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ho", false, 1.0, 0.0, 0.0, 1, false},
                     {"Er", false, 1.0, 0.0, 0.0, 1, false},
                     {"Tm", false, 1.0, 0.0, 0.0, 1, false},
                     {"Yb", false, 1.0, 0.0, 0.0, 1, false},
                     {"Lu", false, 1.0, 0.0, 0.0, 1, false},
                     {"Hf", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ta", false, 1.0, 0.0, 0.0, 1, false},
                     {"W",  false, 1.0, 0.0, 0.0, 1, false},
                     {"Re", false, 1.0, 0.0, 0.0, 1, false},
                     {"Os", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ir", false, 1.0, 0.0, 0.0, 1, false},
                     {"Pt", false, 1.0, 0.0, 0.0, 1, false},
                     {"Au", false, 1.0, 0.0, 0.0, 1, false},
                     {"Hg", false, 1.0, 0.0, 0.0, 1, false},
                     {"Tl", false, 1.0, 0.0, 0.0, 1, false},
                     {"Pb", false, 1.0, 0.0, 0.0, 1, false},
                     {"Bi", false, 1.0, 0.0, 0.0, 1, false},
                     {"Po", false, 1.0, 0.0, 0.0, 1, false},
                     {"At", false, 1.0, 0.0, 0.0, 1, false},
                     {"Rn", false, 1.0, 0.0, 0.0, 1, false},
                     {"Fr", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ra", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ac", false, 1.0, 0.0, 0.0, 1, false},
                     {"Th", false, 1.0, 0.0, 0.0, 1, false},
                     {"Pa", false, 1.0, 0.0, 0.0, 1, false},
                     {"U",  false, 1.0, 0.0, 0.0, 1, false},
                     {"Np", false, 1.0, 0.0, 0.0, 1, false},
                     {"Pu", false, 1.0, 0.0, 0.0, 1, false},
                     {"Am", false, 1.0, 0.0, 0.0, 1, false},
                     {"Cm", false, 1.0, 0.0, 0.0, 1, false},
                     {"Bk", false, 1.0, 0.0, 0.0, 1, false},
                     {"Cf", false, 1.0, 0.0, 0.0, 1, false},
                     {"Es", false, 1.0, 0.0, 0.0, 1, false},
                     {"Fm", false, 1.0, 0.0, 0.0, 1, false},
                     {"Md", false, 1.0, 0.0, 0.0, 1, false},
                     {"No", false, 1.0, 0.0, 0.0, 1, false},
                     {"Lr", false, 1.0, 0.0, 0.0, 1, false},
                     {"Rf", false, 1.0, 0.0, 0.0, 1, false},
                     {"Db", false, 1.0, 0.0, 0.0, 1, false},
                     {"Sg", false, 1.0, 0.0, 0.0, 1, false},
                     {"Bh", false, 1.0, 0.0, 0.0, 1, false},
                     {"Hs", false, 1.0, 0.0, 0.0, 1, false},
                     {"Mt", false, 1.0, 0.0, 0.0, 1, false},
                     {"Ds", false, 1.0, 0.0, 0.0, 1, false},
                     {"probe-N2", false, 14.00674, 0.0, 0.0, 6, false}},
                    {{1.0, 2.0 * 1.00},                 // custom
                     {1.0, 2.0 * 1.09},
                     {1.0, 2.0 * 1.09},
                     {1.0, 2.0 * 1.40},
                     {1.0, 2.0 * 1.82},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 1.70},
                     {1.0, 2.0 * 1.55},
                     {1.0, 2.0 * 1.52},
                     {1.0, 2.0 * 1.47},
                     {1.0, 2.0 * 1.54},
                     {1.0, 2.0 * 2.27},
                     {1.0, 2.0 * 1.73},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.10},
                     {1.0, 2.0 * 1.80},
                     {1.0, 2.0 * 1.80},
                     {1.0, 2.0 * 1.75},
                     {1.0, 2.0 * 1.88},
                     {1.0, 2.0 * 2.75},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 1.63},
                     {1.0, 2.0 * 1.40},
                     {1.0, 2.0 * 1.39},
                     {1.0, 2.0 * 1.87},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 1.85},
                     {1.0, 2.0 * 1.90},
                     {1.0, 2.0 * 1.85},
                     {1.0, 2.0 * 2.02},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 1.63},
                     {1.0, 2.0 * 1.72},
                     {1.0, 2.0 * 1.58},
                     {1.0, 2.0 * 1.93},
                     {1.0, 2.0 * 2.17},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.06},
                     {1.0, 2.0 * 1.98},
                     {1.0, 2.0 * 2.16},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 1.72},
                     {1.0, 2.0 * 1.66},
                     {1.0, 2.0 * 1.55},
                     {1.0, 2.0 * 1.96},
                     {1.0, 2.0 * 2.02},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 1.86},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {1.0, 2.0 * 2.00},
                     {91.5, 3.681}},             // N2
                    ForceField::MixingRule::Lorentz_Berthelot, rc, rc, rc, shifted, tailCorrections, useEwald);
}
