module;

#ifdef USE_LEGACY_HEADERS
#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <string>
#include <string_view>
#include <optional>
#include <numbers>
#include <algorithm>
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <exception>
#include <source_location>
#include <complex>
#include <type_traits>
#include <iterator>
#include <functional>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

module forcefield;

#ifndef USE_LEGACY_HEADERS
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
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import archive;
import json;
import skelement;
import units;
import int3;
import double3;
import double4;
import stringutils;
import pseudo_atom;
import vdwparameters;


ForceField::ForceField(std::vector<PseudoAtom> pseudoAtoms, std::vector<VDWParameters> selfInteractions, 
                       [[maybe_unused]] MixingRule mixingRule, double cutOff, bool shifted, 
                       bool applyTailCorrections) noexcept(false) :
    data(pseudoAtoms.size()* pseudoAtoms.size(), VDWParameters(0.0, 0.0)),
    shiftPotentials(pseudoAtoms.size()* pseudoAtoms.size(), shifted),
    tailCorrections(pseudoAtoms.size()* pseudoAtoms.size(), applyTailCorrections),
    cutOffVDW(cutOff), 
    cutOffCoulomb(cutOff),
    numberOfPseudoAtoms(pseudoAtoms.size()), 
    pseudoAtoms(pseudoAtoms)
{
  for (size_t i = 0; i < selfInteractions.size(); ++i)
  {
    data[i + i * numberOfPseudoAtoms] = selfInteractions[i];
  }

  applyMixingRule();
  preComputePotentialShift();
  preComputeTailCorrection();
}

void ForceField::applyMixingRule()
{
  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (size_t j = i + 1; j < numberOfPseudoAtoms; ++j)
    {
      double mix0 = 
        std::sqrt(data[i * numberOfPseudoAtoms + i].parameters.x * data[j * numberOfPseudoAtoms + j].parameters.x);
      double mix1 = 
        0.5 * (data[i * numberOfPseudoAtoms + i].parameters.y + data[j * numberOfPseudoAtoms + j].parameters.y);
      
      data[i * numberOfPseudoAtoms + j].parameters.x = mix0;
      data[i * numberOfPseudoAtoms + j].parameters.y = mix1;
      data[j * numberOfPseudoAtoms + i].parameters.x = mix0;
      data[j * numberOfPseudoAtoms + i].parameters.y = mix1;
    }
  }
}

void ForceField::preComputePotentialShift()
{
  for (size_t i = 0; i < data.size(); ++i)
  {
    if (shiftPotentials[i])
    {
      data[i].computeShiftAtCutOff(cutOffVDW);
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

      if(tailCorrections[i * numberOfPseudoAtoms + j])
      {
        switch (data[i * numberOfPseudoAtoms + j].type)
        {
        case VDWParameters::Type::LennardJones:
        {
          double arg1 = data[i * numberOfPseudoAtoms + j].parameters.x;
          double arg2 = data[i * numberOfPseudoAtoms + j].parameters.y;
          double term3 = (arg2/cutOffVDW) * (arg2/cutOffVDW) * (arg2/cutOffVDW);
          double term6 = term3 * term3;
          data[i * numberOfPseudoAtoms + j].tailCorrectionEnergy = 
            (4.0 / 3.0) * arg1 * arg2 * arg2 * arg2 * ((1.0 / 3.0) * term6 * term3 - term3);
          break;
        }
        default:
          data[i * numberOfPseudoAtoms + j].tailCorrectionEnergy = 0.0;
          break;
        }
      }
    }
  }
}


std::optional<ForceField> ForceField::readForceField(std::optional<std::string> directoryName, std::string forceFieldFileName) noexcept(false)
{
  // try to look in directory 'directoryName' if set, otherwise the local directory
  std::filesystem::path forceFieldPathfile = std::filesystem::path(directoryName.value_or(".") + "/" + forceFieldFileName);
  if (!std::filesystem::exists(forceFieldPathfile)) 
  {
    // if not found, try the install directory and directory 'directoryName' in 'share/raspa3/forcefields'
    const char* env_p = std::getenv("RASPA_DIR");
    forceFieldPathfile = std::filesystem::path(std::string(env_p) + "/share/raspa3/forcefields/" 
                                              + directoryName.value_or(".") + "/" + forceFieldFileName);
    if (!std::filesystem::exists(forceFieldPathfile))
    {
      return std::nullopt;
    }
  }

  std::ifstream forceFieldStream{ forceFieldPathfile };
  if (!forceFieldStream) 
  {
    return std::nullopt;
  }

  nlohmann::basic_json<nlohmann::raspa_map> parsed_data{};

  try
  {
    parsed_data = nlohmann::json::parse(forceFieldStream);
  }
  catch (nlohmann::json::parse_error& ex)
  {
    throw std::runtime_error(std::format("[Forcefield reader]: Parse error of file {} at byte {}\n{}\n", forceFieldFileName, ex.byte, ex.what()));
  }

  if(!parsed_data.contains("PseudoAtoms"))
  {
    throw std::runtime_error(std::format("[Forcefield reader]: No pseudo-atoms found [keyword 'PseudoAtoms' missing]\n"));
  }
  size_t numberOfPseudoAtoms = parsed_data["PseudoAtoms"].size();

  if(numberOfPseudoAtoms == 0)
  {
    throw std::runtime_error(std::format("[ReadPseudoAtoms]: key 'PseudoAtoms' empty]\n"));
  }

  std::vector<PseudoAtom> jsonPseudoAtoms{};
  jsonPseudoAtoms.reserve(numberOfPseudoAtoms);

  for (auto& [_, item] : parsed_data["PseudoAtoms"].items())
  {
    std::string jsonName = item["name"].is_string() ? item["name"].get<std::string>() : std::string{};
    double jsonMass = item["mass"].is_number() ? item["mass"].get<double>() : 0.0;
    std::string jsonElement = item["element"].is_string() ? item["element"].get<std::string>() : "C";
    double jsonCharge = item["charge"].is_number() ? item["charge"].get<double>() : 0.0;
    size_t jsonPrintToOutput = item["print_to_output"].is_boolean() ? item["print_to_output"].get<bool>() : true;
    std::string jsonSource = item["source"].is_string() ? item["source"].get<std::string>() : std::string{};

    size_t atomicNumber{ 1 };
    auto it = PredefinedElements::atomicNumberData.find(jsonElement);
    if (it != PredefinedElements::atomicNumberData.end())
    {
      atomicNumber = static_cast<size_t>(it->second);
    }

    jsonPseudoAtoms.emplace_back(jsonName, jsonMass, jsonCharge, atomicNumber, jsonPrintToOutput, jsonSource);
  }

  std::vector<VDWParameters> jsonSelfInteractions(numberOfPseudoAtoms);

  if(!parsed_data.contains("SelfInteractions"))
  {
    throw std::runtime_error(std::format("[ReadForceFieldSelfInteractions]: No pseudo-atoms found [keyword 'SelfInteractions' missing]\n"));
  }
  size_t jsonNumberOfPseudoAtoms = parsed_data["SelfInteractions"].size();

  if(jsonNumberOfPseudoAtoms == 0)
  {
    throw std::runtime_error(std::format("[ReadForceFieldSelfInteractions]: key 'SelfInteractions' empty]\n"));
  }

  for (auto& [_, item] : parsed_data["SelfInteractions"].items())
  {
    std::string jsonName = item["name"].is_string() ? item["name"].get<std::string>() : std::string{};
    std::string jsonType = item["type"].is_string() ? item["type"].get<std::string>() : "lennard-jones";
    std::string jsonSource = item["source"].is_string() ? item["source"].get<std::string>() : std::string{};
    std::vector<double> scannedJsonParameters{};
    try
    {
      scannedJsonParameters = item["parameters"].is_array() ? item["parameters"].get<std::vector<double>>() : std::vector<double>{};
    }
    catch (nlohmann::json::exception& ex)
    {
      throw std::runtime_error(std::format("[ReadForceFieldSelfInteractions]: parameters {} must be array of numbers \n{}\n", 
             item["parameters"].dump(), ex.what()));
    }

    std::optional<size_t> index = ForceField::findPseudoAtom(jsonPseudoAtoms, jsonName);

    if(!index.has_value())
    {
      throw std::runtime_error(std::format("[ReadForceFieldSelfInteractions]: unknown pseudo-atom '{}', please define in 'pseudo_atoms.json'\n", jsonName));
    }

    if(scannedJsonParameters.size() < 2)
    {
      throw std::runtime_error(std::format("[ReadForceFieldSelfInteractions]: incorrect vdw parameters {}\n", item["parameters"].dump()));
    }

    double param0 = scannedJsonParameters[0];
    double param1 = scannedJsonParameters[1];

    jsonSelfInteractions[index.value()] = VDWParameters(param0, param1);
  }

  MixingRule jsonMixingRule{MixingRule::Lorentz_Berthelot};
  if(parsed_data["MixingRule"].is_string())
  {
    if(caseInSensStringCompare(parsed_data["MixingRule"].get<std::string>(), "Lorentz-Berthelot"))
    {
    }
  }

  bool jsonShiftPotentials{true};
  if(parsed_data["TruncationMethod"].is_string())
  {
    if(caseInSensStringCompare(parsed_data["TruncationMethod"].get<std::string>(), "shifted"))
    {
      jsonShiftPotentials = true;
    }
    else if(caseInSensStringCompare(parsed_data["TruncationMethod"].get<std::string>(), "truncated"))
    {
      jsonShiftPotentials = false;
    }
  }

  bool jsonTailCorrections{false};
  if(parsed_data["TailCorrections"].is_boolean())
  {
    if(parsed_data["TailCorrections"].get<bool>())
    {
      jsonTailCorrections = true;
    }
  }

  double cutOff = 12.0;
  return ForceField{jsonPseudoAtoms, jsonSelfInteractions, jsonMixingRule, cutOff, jsonShiftPotentials, jsonTailCorrections};
}


std::string ForceField::printPseudoAtomStatus() const
{
  std::ostringstream stream;
 
  std::print(stream, "Pseudo-atoms\n");
  std::print(stream, "===============================================================================\n\n");

  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
      std::print(stream, "{:3d} - {:8} mass: {:8.5f}, charge: {:8.5f}\n", 
                         i, pseudoAtoms[i].name, pseudoAtoms[i].mass, pseudoAtoms[i].charge);
  }
  std::print(stream, "\n");

  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    if(!pseudoAtoms[i].source.empty())
    {
      std::print(stream, "{:3d} - {:8} {}\n", i, pseudoAtoms[i].name, pseudoAtoms[i].source);
    }
  }
  std::print(stream, "\n");

  return stream.str();
}

std::string ForceField::printForceFieldStatus() const
{
  std::ostringstream stream;

  std::print(stream, "Force field status\n");
  std::print(stream, "===============================================================================\n\n");

  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (size_t j = i; j < numberOfPseudoAtoms; ++j)
    {
      switch (data[i * numberOfPseudoAtoms + j].type)
      {
      case VDWParameters::Type::LennardJones:
        std::print(stream, "{:8} - {:8} {} p₀/kʙ: {:9.5f} [K], p₁: {:8.5f} [Å]\n",
            pseudoAtoms[i].name, pseudoAtoms[j].name, "Lennard-Jones",
            Units::EnergyToKelvin * data[i * numberOfPseudoAtoms + j].parameters.x,
            data[i * numberOfPseudoAtoms + j].parameters.y);
        std::print(stream, "{:33} shift: {:9.5f} [K], tailcorrections: {}\n",
            std::string(""),
            Units::EnergyToKelvin * data[i * numberOfPseudoAtoms + j].shift,
            tailCorrections[i * numberOfPseudoAtoms + j] ? "true" : "false");
        break;
      default:
        break;
      }
    }
  }
  std::print(stream, "\n");
  if(automaticEwald)
  {
    std::print(stream, "Ewald precision: {}\n", EwaldPrecision);
    std::print(stream, "Ewald alpha: {}\n", EwaldAlpha);
    std::print(stream, "Ewald k-vectors: {} {} {}\n", 
                       numberOfWaveVectors.x, numberOfWaveVectors.y, numberOfWaveVectors.z);
  }
  else
  {
    std::print(stream, "Ewald alpha: {}\n", EwaldAlpha);
    std::print(stream, "Ewald k-vectors: {} {} {}\n", 
                       numberOfWaveVectors.x, numberOfWaveVectors.y, numberOfWaveVectors.z);
  }
  std::print(stream, "\n\n");

  return stream.str();
}


std::optional<size_t> ForceField::findPseudoAtom(const std::string& name) const
{
  std::vector<PseudoAtom>::const_iterator match = std::find_if(
        pseudoAtoms.begin(), pseudoAtoms.end(),
        [&name](const PseudoAtom& x) { return x.name == name; });
  if (match != std::end(pseudoAtoms))
  {
    return static_cast<size_t>(std::distance(pseudoAtoms.begin(), match));
  }
 
  return std::nullopt;
}

std::optional<size_t> ForceField::findPseudoAtom(const std::vector<PseudoAtom> pseudoAtoms, const std::string& name)
{
  std::vector<PseudoAtom>::const_iterator match = std::find_if(
        pseudoAtoms.begin(), pseudoAtoms.end(),
        [&name](const PseudoAtom& x) { return x.name == name; });
  if (match != std::end(pseudoAtoms))
  {
    return static_cast<size_t>(std::distance(pseudoAtoms.begin(), match));
  }
 
  return std::nullopt;
}

void ForceField::initializeEwaldParameters(double3 perpendicularWidths)
{
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

    //numberOfWavevectors = ((kx_max_unsigned + 1) * (2 * ky_max_unsigned + 1) * (2 * kz_max_unsigned + 1));
    //if (ReciprocalCutOffSquared[i] < 0.0)
    //    ReciprocalCutOffSquared[i] = SQR(1.05 * MAX3(kvec[i].x, kvec[i].y, kvec[i].z));
  }
}


Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ForceField &f)
{
  archive << f.versionNumber;

  archive << f.data;
  archive << f.shiftPotentials;
  archive << f.tailCorrections;
  archive << f.cutOffVDW;
  archive << f.cutOffCoulomb;
  archive << f.dualCutOff;

  archive << f.numberOfPseudoAtoms;
  archive << f.pseudoAtoms;

  archive << f.chargeMethod;

  archive << f.overlapCriteria;

  archive << f.EwaldPrecision;
  archive << f.EwaldAlpha;
  archive << f.numberOfWaveVectors;
  archive << f.automaticEwald;
  archive << f.noCharges;
  archive << f.omitEwaldFourier;
  archive << f.minimumRosenbluthFactor;
  archive << f.energyOverlapCriteria;
  archive << f.useDualCutOff;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ForceField &f)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > f.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'ForceField' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> f.data;
  archive >> f.shiftPotentials;
  archive >> f.tailCorrections;
  archive >> f.cutOffVDW;
  archive >> f.cutOffCoulomb;
  archive >> f.dualCutOff;

  archive >> f.numberOfPseudoAtoms;
  archive >> f.pseudoAtoms;

  archive >> f.chargeMethod;

  archive >> f.overlapCriteria;

  archive >> f.EwaldPrecision;
  archive >> f.EwaldAlpha;
  archive >> f.numberOfWaveVectors;
  archive >> f.automaticEwald;
  archive >> f.noCharges;
  archive >> f.omitEwaldFourier;
  archive >> f.minimumRosenbluthFactor;
  archive >> f.energyOverlapCriteria;
  archive >> f.useDualCutOff;

  return archive;
}

bool ForceField::operator==(const ForceField& other) const
{
  // first the cheap ones
  if (cutOffVDW != other.cutOffVDW || cutOffCoulomb != other.cutOffCoulomb || dualCutOff != other.dualCutOff ||
      numberOfPseudoAtoms != other.numberOfPseudoAtoms || overlapCriteria != other.overlapCriteria ||
      EwaldPrecision != other.EwaldPrecision || EwaldAlpha != other.EwaldAlpha ||
      numberOfWaveVectors != other.numberOfWaveVectors || automaticEwald != other.automaticEwald ||
      noCharges != other.noCharges || omitEwaldFourier != other.omitEwaldFourier ||
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
