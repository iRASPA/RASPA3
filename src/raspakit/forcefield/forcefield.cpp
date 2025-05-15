module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <numbers>
#include <optional>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <cstring>
#include <string_view>
#include <type_traits>
#include <vector>
#include <set>
#endif

module forcefield;

#ifndef USE_LEGACY_HEADERS
import <cstdint>;
import <filesystem>;
import <fstream>;
import <cstdlib>;
import <iostream>;
import <sstream>;
import <vector>;
import <array>;
import <map>;
import <cmath>;
import <string>;
import <string_view>;
import <optional>;
import <numbers>;
import <algorithm>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <type_traits>;
import <iterator>;
import <functional>;
import <print>;
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
      value.x = json[0].template get<int32_t>();
      value.y = json[1].template get<int32_t>();
      value.z = json[2].template get<int32_t>();
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

  for (size_t i = 0; i < selfInteractions.size(); ++i)
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
    int64_t oxidationState = item.value("oxidationState", 0);
    double jsonMass = item.value("mass", 0.0);
    std::string jsonElement = item.value("element", "C");
    double jsonCharge = item.value("charge", 0.0);
    double jsonPolarizibility = item.value("polarizibility", 0.0);
    bool jsonPrintToOutput = item.value("print_to_output", true);
    std::string jsonSource = item.value("source", "");

    size_t atomicNumber = PredefinedElements::atomicNumberData.contains(jsonElement)
                              ? static_cast<size_t>(PredefinedElements::atomicNumberData.at(jsonElement))
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

    std::optional<size_t> index = ForceField::findPseudoAtom(jsonName);
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

    std::optional<size_t> indexA = ForceField::findPseudoAtom(scannedJsonParameterNames[0]);
    if (!indexA.has_value())
    {
      throw std::runtime_error(
          std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define'\n",
                      scannedJsonParameterNames[0]));
    }
    std::optional<size_t> indexB = ForceField::findPseudoAtom(scannedJsonParameterNames[1]);
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
    if(parsed_data["EwaldPrecision"].is_string())
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
    if(parsed_data["EwaldParameters"].is_array())
    {
      if (parsed_data["EwaldParameters"].size() != 4)
      {
        throw std::runtime_error(
            std::format("[ForceField reader]: key '{}', value {} should be array of one floaying point number and 3 integer numbers\n", 
                        "EwaldParameters", parsed_data["EwaldParameters"].dump()));
      }
      double alpha;
      int3 value{};
      try
      {
        alpha = parsed_data["EwaldParameters"][0].template get<double>();
        value.x = parsed_data["EwaldParameters"][1].template get<int32_t>();
        value.y = parsed_data["EwaldParameters"][2].template get<int32_t>();
        value.z = parsed_data["EwaldParameters"][3].template get<int32_t>();
      }
      catch (nlohmann::json::exception& ex)
      {
        throw std::runtime_error(
            std::format("[ForceField reader]: key '{}', value {} should be array of one floaying point number and 3 integer numbers\n", 
                        "EwaldParameters", parsed_data["EwaldParameters"].dump()));
      }

      EwaldAlpha = alpha;
      numberOfWaveVectors = value;
      automaticEwald = false;
    }
  }

  if (parsed_data.contains("CutOff"))
  {
    if(parsed_data["CutOff"].is_string())
    {
      std::string cutOffCoulombString = parsed_data["CutOff"].get<std::string>();

      if (caseInSensStringCompare(cutOffCoulombString, "auto"))
      {
        cutOffFrameworkVDWAutomatic = true;
        cutOffMoleculeVDWAutomatic = true;
      }
    }

    if(parsed_data["CutOff"].is_number_float())
    {
      cutOffFrameworkVDWAutomatic = false;
      cutOffMoleculeVDWAutomatic = false;
      cutOffFrameworkVDW = parsed_data["CutOff"].get<double>();
      cutOffMoleculeVDW = parsed_data["CutOff"].get<double>();
    }
  }

  if (parsed_data.contains("CutOffFrameworkVDW"))
  {
    if(parsed_data["CutOffFrameworkVDW"].is_string())
    {
      std::string cutOffCoulombString = parsed_data["CutOffFrameworkVDW"].get<std::string>();

      if (caseInSensStringCompare(cutOffCoulombString, "auto"))
      {
        cutOffFrameworkVDWAutomatic = true;
      }
    }

    if(parsed_data["CutOffFrameworkVDW"].is_number_float())
    {
      cutOffFrameworkVDWAutomatic = false;
      cutOffFrameworkVDW = parsed_data["CutOffFrameworkVDW"].get<double>();
    }
  }

  if (parsed_data.contains("CutOffMoleculeVDW"))
  {
    if(parsed_data["CutOffMoleculeVDW"].is_string())
    {
      std::string cutOffCoulombString = parsed_data["CutOffMoleculeVDW"].get<std::string>();

      if (caseInSensStringCompare(cutOffCoulombString, "auto"))
      {
        cutOffMoleculeVDWAutomatic = true;
      }
    }

    if(parsed_data["CutOffMoleculeVDW"].is_number_float())
    {
      cutOffMoleculeVDWAutomatic = false;
      cutOffMoleculeVDW = parsed_data["CutOffMoleculeVDW"].get<double>();
    }
  }

  if (parsed_data.contains("CutOffCoulomb"))
  {
    if(parsed_data["CutOffCoulomb"].is_string())
    {
      std::string cutOffCoulombString = parsed_data["CutOffCoulomb"].get<std::string>();

      if (caseInSensStringCompare(cutOffCoulombString, "auto"))
      {
        cutOffCoulombAutomatic = true;
      }
    }

    if(parsed_data["CutOffCoulomb"].is_number_float())
    {
      cutOffCoulombAutomatic = false;
      cutOffCoulomb = parsed_data["CutOffCoulomb"].get<double>();
    }
  }

  if (parsed_data.contains("ShiftedPotentialPairs"))
  {
    if(parsed_data["ShiftedPotentialPairs"].is_array())
    {
      for(auto &pair : parsed_data["ShiftedPotentialPairs"])
      {
        if (pair.size() != 2)
        {
          throw std::runtime_error(
              std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n", "ShiftedPotentialPairs", pair.dump()));
        }

        std::string stringA, stringB;
        try
        {
          stringA = pair[0].template get<std::string>();
          stringB = pair[1].template get<std::string>();
        }
        catch (nlohmann::json::exception& ex)
        {
          throw std::runtime_error(
              std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n", "ShiftedPotentialPairs", pair.dump()));
        }

        std::optional<size_t> indexA = findPseudoAtom(stringA);
        if (!indexA.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringA, pair.dump()));
        }

        std::optional<size_t> indexB = findPseudoAtom(stringB);
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
    if(parsed_data["TruncatedPotentialPairs"].is_array())
    {
      for(auto &pair : parsed_data["TruncatedPotentialPairs"])
      {
        if (pair.size() != 2)
        {
          throw std::runtime_error(
              std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n", "TruncatedPotentialPairs", pair.dump()));
        }

        std::string stringA, stringB;
        try
        {
          stringA = pair[0].template get<std::string>();
          stringB = pair[1].template get<std::string>();
        }
        catch (nlohmann::json::exception& ex)
        {
          throw std::runtime_error(
              std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n", "TruncatedPotentialPairs", pair.dump()));
        }


        std::optional<size_t> indexA = findPseudoAtom(stringA);
        if (!indexA.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringA, pair.dump()));
        }

        std::optional<size_t> indexB = findPseudoAtom(stringB);
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
    if(parsed_data["TailCorrectionPairs"].is_array())
    {
      for(auto &pair : parsed_data["TailCorrectionPairs"])
      {
        if (pair.size() != 2)
        {
          throw std::runtime_error(
              std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n", "TailCorrectionPairs", pair.dump()));
        }

        std::string stringA, stringB;
        try
        {
          stringA = pair[0].template get<std::string>();
          stringB = pair[1].template get<std::string>();
        }
        catch (nlohmann::json::exception& ex)
        {
          throw std::runtime_error(
              std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n", "TailCorrectionPairs", pair.dump()));
        }

        std::optional<size_t> indexA = findPseudoAtom(stringA);
        if (!indexA.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringA, pair.dump()));
        }

        std::optional<size_t> indexB = findPseudoAtom(stringB);
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
    if(parsed_data["NoTailCorrectionPairs"].is_array())
    {
      for(auto &pair : parsed_data["NoTailCorrectionPairs"])
      {
        if (pair.size() != 2)
        {
          throw std::runtime_error(
              std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n", "NoTailCorrectionPairs", pair.dump()));
        }

        std::string stringA, stringB;
        try
        {
          stringA = pair[0].template get<std::string>();
          stringB = pair[1].template get<std::string>();
        }
        catch (nlohmann::json::exception& ex)
        {
          throw std::runtime_error(
              std::format("[ForceField reader]: key '{}', value {} should be array of 2 strings\n", "NoTailCorrectionPairs", pair.dump()));
        }

        std::optional<size_t> indexA = findPseudoAtom(stringA);
        if (!indexA.has_value())
        {
          throw std::runtime_error(std::format("[ForceField]: unknown pseudo atom {} in {}\n", stringA, pair.dump()));
        }

        std::optional<size_t> indexB = findPseudoAtom(stringB);
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
    std::optional<size_t> index = findPseudoAtom(pseudo_atom_name);

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
    numberOfVDWGridPoints = parseInt3("NumberOfVDWGridPoints", parsed_data["NumberOfVDWGridPoints"]);
  }
  if (parsed_data.contains("NumberOfCoulombGridPoints"))
  {
    numberOfCoulombGridPoints = parseInt3("NumberOfCoulombGridPoints", parsed_data["NumberOfCoulombGridPoints"]);
  }

  if (parsed_data.contains("NumberOfGridTestPoints"))
  {
    numberOfGridTestPoints = parsed_data.value("NumberOfGridTestPoints", parsed_data["NumberOfGridTestPoints"]);
  }

  if (parsed_data.contains("InterpolationScheme"))
  {
    size_t scheme = parsed_data.value("InterpolationScheme", parsed_data["InterpolationScheme"]);
    switch (scheme)
    {
      case 3:
        interpolationScheme = InterpolationScheme::Tricubic;
        break;
      case 5:
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
      for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
      {
        for (size_t j = i + 1; j < numberOfPseudoAtoms; ++j)
        {
          if (data[i * numberOfPseudoAtoms + i].type == VDWParameters::Type::LennardJones &&
              data[i * numberOfPseudoAtoms + i].type == VDWParameters::Type::LennardJones)
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
      for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
      {
        for (size_t j = i + 1; j < numberOfPseudoAtoms; ++j)
        {
          if (data[i * numberOfPseudoAtoms + i].type == VDWParameters::Type::LennardJones &&
              data[i * numberOfPseudoAtoms + i].type == VDWParameters::Type::LennardJones)
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
}

double ForceField::cutOffVDW(size_t i, size_t j) const
{
  if (pseudoAtoms[i].framework || pseudoAtoms[j].framework)
  {
    return cutOffFrameworkVDW;
  }

  return cutOffMoleculeVDW;
}

void ForceField::preComputePotentialShift()
{
  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (size_t j = 0; j < numberOfPseudoAtoms; ++j)
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
  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (size_t j = 0; j < numberOfPseudoAtoms; ++j)
    {
      data[i * numberOfPseudoAtoms + j].tailCorrectionEnergy = 0.0;
      data[i * numberOfPseudoAtoms + j].tailCorrectionPressure = 0.0;

      if (tailCorrections[i * numberOfPseudoAtoms + j])
      {
        data[i * numberOfPseudoAtoms + j].tailCorrectionEnergy = Potentials::potentialCorrectionVDW(*this, i, j);
        data[i * numberOfPseudoAtoms + j].tailCorrectionPressure = Potentials::potentialCorrectionPressure(*this, i, j);
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

  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    std::print(stream, "{:3d} - {:8} framework-atom: {}\n", i, pseudoAtoms[i].name, pseudoAtoms[i].framework);
  }
  std::print(stream, "\n");
  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    std::print(stream, "{:3d} - {:8} mass: {:8.5f}, charge: {:8.5f}\n", i, pseudoAtoms[i].name, pseudoAtoms[i].mass,
               pseudoAtoms[i].charge);
  }
  std::print(stream, "\n");
  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    std::print(stream, "{:3d} - {:8} polarizability: {:8.5f}\n", i, pseudoAtoms[i].name, pseudoAtoms[i].polarizability);
  }
  std::print(stream, "\n");

  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
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

  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
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

std::string ForceField::printForceFieldStatus() const
{
  std::ostringstream stream;

  std::print(stream, "Force field status\n");
  std::print(stream, "===============================================================================\n\n");

  if(cutOffFrameworkVDWAutomatic)
  {
    std::print(stream, "Cutoff Framework-Molecule VDW:  auto\n");
  }
  else
  {
    std::print(stream, "Cutoff Framework-Molecule VDW: {:9.5f} [{}]\n", cutOffFrameworkVDW,
               Units::displayedUnitOfLengthString);
  }
  if(cutOffMoleculeVDWAutomatic)
  {
    std::print(stream, "Cutoff Molecule-Molecule VDW:   auto\n");
  }
  else
  {
    std::print(stream, "Cutoff Molecule-Molecule VDW:  {:9.5f} [{}]\n", cutOffFrameworkVDW,
               Units::displayedUnitOfLengthString);
  }
  if(cutOffCoulombAutomatic)
  {
    std::print(stream, "Cutoff Coulomb:                 auto\n");
  }
  else
  {
    std::print(stream, "Cutoff Coulomb:                {:9.5f} [{}]\n\n", cutOffCoulomb,
               Units::displayedUnitOfLengthString);
  }

  std::print(stream, "Overlap-criteria VDW:          {: .6e} [{}]\n\n", overlapCriteria,
             Units::displayedUnitOfEnergyString);

  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (size_t j = i; j < numberOfPseudoAtoms; ++j)
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
  size_t n_interactions = static_cast<size_t>(static_cast<double>(numberOfPseudoAtoms) *
                                              (static_cast<double>(numberOfPseudoAtoms) + 1.0) / 2.0);
  std::vector<nlohmann::json> interactions(n_interactions);

  size_t count = 0;
  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (size_t j = i; j < numberOfPseudoAtoms; ++j)
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

std::optional<size_t> ForceField::findPseudoAtom(const std::string& name) const
{
  std::vector<PseudoAtom>::const_iterator match =
      std::find_if(pseudoAtoms.begin(), pseudoAtoms.end(), [&name](const PseudoAtom& x) { return x.name == name; });
  if (match != std::end(pseudoAtoms))
  {
    return static_cast<size_t>(std::distance(pseudoAtoms.begin(), match));
  }

  return std::nullopt;
}

std::optional<size_t> ForceField::findPseudoAtom(const std::vector<PseudoAtom> pseudoAtoms, const std::string& name)
{
  std::vector<PseudoAtom>::const_iterator match =
      std::find_if(pseudoAtoms.begin(), pseudoAtoms.end(), [&name](const PseudoAtom& x) { return x.name == name; });
  if (match != std::end(pseudoAtoms))
  {
    return static_cast<size_t>(std::distance(pseudoAtoms.begin(), match));
  }

  return std::nullopt;
}

void ForceField::initializeEwaldParameters(const SimulationBox& simulationBox)
{
  double3 perpendicularWidths = simulationBox.perpendicularWidths();

  if (automaticEwald)
  {
    // compute the alpha-parameter and max k-vectors from the relative precision
    double eps = std::min(fabs(EwaldPrecision), 0.5);

    double tol = std::sqrt(std::abs(std::log(eps * cutOffCoulomb)));

    EwaldAlpha = std::sqrt(std::abs(std::log(eps * cutOffCoulomb * tol))) / cutOffCoulomb;
    double tol1 = std::sqrt(-std::log(eps * cutOffCoulomb * (2.0 * tol * EwaldAlpha) * (2.0 * tol * EwaldAlpha)));

    numberOfWaveVectors =
        int3(static_cast<int32_t>(rint(0.25 + perpendicularWidths.x * EwaldAlpha * tol1 / std::numbers::pi)),
             static_cast<int32_t>(rint(0.25 + perpendicularWidths.y * EwaldAlpha * tol1 / std::numbers::pi)),
             static_cast<int32_t>(rint(0.25 + perpendicularWidths.z * EwaldAlpha * tol1 / std::numbers::pi)));

    size_t maxNumberOfWaveVector =
        static_cast<size_t>(std::max({numberOfWaveVectors.x, numberOfWaveVectors.y, numberOfWaveVectors.z}));
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

  archive << f.overlapCriteria;

  archive << f.EwaldPrecision;
  archive << f.EwaldAlpha;
  archive << f.numberOfWaveVectors;
  archive << f.reciprocalIntegerCutOffSquared;
  archive << f.reciprocalCutOffSquared;
  archive << f.automaticEwald;
  archive << f.useCharge;
  archive << f.omitEwaldFourier;
  archive << f.minimumRosenbluthFactor;
  archive << f.energyOverlapCriteria;
  archive << f.useDualCutOff;

  archive << f.omitInterInteractions;

  archive << f.computePolarization;
  archive << f.omitInterPolarization;

  archive << f.hasExternalField;
  archive << f.potentialEnergySurfaceType;
  archive << f.potentialEnergySurfaceOrigin;

  archive << f.gridPseudoAtomIndices;
  archive << f.spacingVDWGrid;
  archive << f.spacingCoulombGrid;
  archive << f.numberOfVDWGridPoints;
  archive << f.numberOfCoulombGridPoints;
  archive << f.numberOfGridTestPoints;
  archive << f.interpolationScheme;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ForceField& f)
{
  uint64_t versionNumber;
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

  archive >> f.overlapCriteria;

  archive >> f.EwaldPrecision;
  archive >> f.EwaldAlpha;
  archive >> f.numberOfWaveVectors;
  archive >> f.reciprocalIntegerCutOffSquared;
  archive >> f.reciprocalCutOffSquared;
  archive >> f.automaticEwald;
  archive >> f.useCharge;
  archive >> f.omitEwaldFourier;
  archive >> f.minimumRosenbluthFactor;
  archive >> f.energyOverlapCriteria;
  archive >> f.useDualCutOff;

  archive >> f.omitInterInteractions;

  archive >> f.computePolarization;
  archive >> f.omitInterPolarization;

  archive >> f.hasExternalField;
  archive >> f.potentialEnergySurfaceType;
  archive >> f.potentialEnergySurfaceOrigin;

  archive >> f.gridPseudoAtomIndices;
  archive >> f.spacingVDWGrid;
  archive >> f.spacingCoulombGrid;
  archive >> f.numberOfVDWGridPoints;
  archive >> f.numberOfCoulombGridPoints;
  archive >> f.numberOfGridTestPoints;
  archive >> f.interpolationScheme;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
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
      numberOfPseudoAtoms != other.numberOfPseudoAtoms || overlapCriteria != other.overlapCriteria ||
      EwaldPrecision != other.EwaldPrecision || EwaldAlpha != other.EwaldAlpha ||
      numberOfWaveVectors != other.numberOfWaveVectors || automaticEwald != other.automaticEwald ||
      useCharge != other.useCharge || omitEwaldFourier != other.omitEwaldFourier ||
      minimumRosenbluthFactor != other.minimumRosenbluthFactor ||
      energyOverlapCriteria != other.energyOverlapCriteria || useDualCutOff != other.useDualCutOff ||
      chargeMethod != other.chargeMethod)
  {
    return false;
  }
  for (size_t idx = 0; idx < shiftPotentials.size(); idx++)
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

const std::set<std::string, ForceField::InsensitiveCompare> ForceField::options = {
    "MixingRule",
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
    "NumberOfGridTestPoints",
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
