module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
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
#include <string_view>
#include <type_traits>
#include <vector>
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

ForceField::ForceField(std::vector<PseudoAtom> pseudoAtoms, std::vector<VDWParameters> selfInteractions,
                       [[maybe_unused]] MixingRule mixingRule, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
                       double cutOffCoulomb, bool shifted, bool applyTailCorrections, bool useCharge) noexcept(false)
    : data(pseudoAtoms.size() * pseudoAtoms.size(), VDWParameters(0.0, 0.0)),
      shiftPotentials(pseudoAtoms.size() * pseudoAtoms.size(), shifted),
      tailCorrections(pseudoAtoms.size() * pseudoAtoms.size(), applyTailCorrections),
      cutOffFrameworkVDW(cutOffFrameworkVDW),
      cutOffMoleculeVDW(cutOffMoleculeVDW),
      cutOffCoulomb(cutOffCoulomb),
      numberOfPseudoAtoms(pseudoAtoms.size()),
      pseudoAtoms(pseudoAtoms),
      useCharge(useCharge)
{
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

  // Validate and read pseudo atoms
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

  // Get charge method
  useCharge = parsed_data.value("ChargeMethod", "Ewald") != "Ewald";

  cutOffFrameworkVDW = parsed_data.value("CutOffFrameworkVDW", 12.0);
  cutOffMoleculeVDW = parsed_data.value("CutOffMoleculeVDW", 12.0);
  cutOffCoulomb = parsed_data.value("CutOffCoulomb", 12.0);

  preComputePotentialShift();
  preComputeTailCorrection();
}

void ForceField::applyMixingRule()
{
  if (mixingRule == MixingRule::Lorentz_Berthelot)
  {
    for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
    {
      for (size_t j = i + 1; j < numberOfPseudoAtoms; ++j)
      {
        if (data[i * numberOfPseudoAtoms + i].type == VDWParameters::Type::LennardJones && data[i * numberOfPseudoAtoms + i].type == VDWParameters::Type::LennardJones)
        {
          double mix0 =
              std::sqrt(data[i * numberOfPseudoAtoms + i].parameters.x * data[j * numberOfPseudoAtoms + j].parameters.x);
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
        data[i * numberOfPseudoAtoms + j].tailCorrectionEnergy = potentialCorrectionVDW(*this, i, j);
        data[i * numberOfPseudoAtoms + j].tailCorrectionPressure = potentialCorrectionPressure(*this, i, j);
        /*
        switch (data[i * numberOfPseudoAtoms + j].type)
        {
          case VDWParameters::Type::LennardJones:
          {
            double arg1 = data[i * numberOfPseudoAtoms + j].parameters.x;
            double arg2 = data[i * numberOfPseudoAtoms + j].parameters.y;
            double cut_off_vdw = cutOffVDW(i, j);
            double term3 = (arg2 / cut_off_vdw) * (arg2 / cut_off_vdw) * (arg2 / cut_off_vdw);
            double term6 = term3 * term3;
            data[i * numberOfPseudoAtoms + j].tailCorrectionEnergy =
                (4.0 / 3.0) * arg1 * arg2 * arg2 * arg2 * ((1.0 / 3.0) * term6 * term3 - term3);
            data[i * numberOfPseudoAtoms + j].tailCorrectionPressure =
                (8.0 / 3.0) * arg1 * arg2 * arg2 * arg2 * ((2.0 / 3.0) * term6 * term3 - term3);
            break;
          }
          default:
            data[i * numberOfPseudoAtoms + j].tailCorrectionEnergy = 0.0;
            data[i * numberOfPseudoAtoms + j].tailCorrectionPressure = 0.0;
            break;
        }
      */
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

  std::print(stream, "Cutoff Framework-Molecule VDW: {:9.5f} [{}]\n", cutOffFrameworkVDW, Units::displayedUnitOfLengthString);
  std::print(stream, "Cutoff Molecule-Molecule VDW:  {:9.5f} [{}]\n", cutOffFrameworkVDW, Units::displayedUnitOfLengthString);
  std::print(stream, "Cutoff Coulomb:                {:9.5f} [{}]\n\n", cutOffCoulomb, Units::displayedUnitOfLengthString);

  std::print(stream, "Overlap-criteria VDW:          {: .6e} [{}]\n\n", overlapCriteria, Units::displayedUnitOfEnergyString);

  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (size_t j = i; j < numberOfPseudoAtoms; ++j)
    {
      switch (data[i * numberOfPseudoAtoms + j].type)
      {
        case VDWParameters::Type::LennardJones:
          std::print(stream, "{:8} - {:8} {} p₀{}: {:9.5f} [{}], p₁: {:8.5f} [{}]\n",
                     pseudoAtoms[i].name,
                     pseudoAtoms[j].name, 
                     "Lennard-Jones",
                     Units::displayedUnitOfEnergyConversionString,
                     Units::EnergyToKelvin * data[i * numberOfPseudoAtoms + j].parameters.x,
                     Units::displayedUnitOfEnergyString,
                     data[i * numberOfPseudoAtoms + j].parameters.y,
                     Units::displayedUnitOfLengthString);
          std::print(stream, "{:33} shift: {:9.5f} [{}], tailcorrections: {}\n", 
                     std::string(""),
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
          interactions[count]["typeB"] = pseudoAtoms[i].name;
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
  archive << f.cutOffFrameworkVDW;
  archive << f.cutOffMoleculeVDW;
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
  archive >> f.cutOffFrameworkVDW;
  archive >> f.cutOffMoleculeVDW;
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
