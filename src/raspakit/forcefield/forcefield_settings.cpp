module;

module forcefield_settings;

import std;

import archive;
import json;
import units;
import uint3;
import stringutils;
import pseudo_atom;

namespace
{
uint3 parseUint3Setting(const std::string& item, const nlohmann::basic_json<nlohmann::raspa_map>& json)
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
}  // namespace

const std::set<std::string> ForceFieldSettings::options = {"UseInterpolationGrids",
                                                           "SpacingVDWGrid",
                                                           "SpacingCoulombGrid",
                                                           "NumberOfVDWGridPoints",
                                                           "NumberOfCoulombGridPoints",
                                                           "NumberOfGridTestPoints",
                                                           "NumberOfExternalFieldGridPoints",
                                                           "UseExternalFieldGrid",
                                                           "WriteExternalFieldInterpolationGrid",
                                                           "InterpolationScheme",
                                                           "UseDualCutOff",
                                                           "DualCutOff",
                                                           "UseRecoilGrowth",
                                                           "RecoilGrowthMaximumRecoilLength",
                                                           "RecoilGrowthNumberOfTrialDirections",
                                                           "NumberOfTrialDirections",
                                                           "NumberOfTorsionTrialDirections",
                                                           "NumberOfFirstBeadPositions",
                                                           "NumberOfTrialMovesPerOpenBead",
                                                           "CBMCRingCrankshaftProbability",
                                                           "CBMCRingTiltProbability"};

bool ForceFieldSettings::isOption(const std::string& key)
{
  return std::ranges::any_of(options, [&key](const std::string& option) { return caseInSensStringCompare(option, key); });
}

void ForceFieldSettings::readFromJSON(const nlohmann::basic_json<nlohmann::raspa_map>& parsed_data,
                                      const std::vector<PseudoAtom>& pseudoAtoms, double smallestExplicitCutOff)
{
  std::vector<std::string> pseudoAtomStringGrids =
      parsed_data.value("UseInterpolationGrids", std::vector<std::string>{});

  for (const std::string& pseudo_atom_name : pseudoAtomStringGrids)
  {
    auto match = std::ranges::find_if(pseudoAtoms, [&](const PseudoAtom& x) { return x.name == pseudo_atom_name; });
    if (match == pseudoAtoms.end())
    {
      throw std::runtime_error(std::format("[ReadForceFieldSelfInteractions]: unknown pseudo atom {} in {}\n",
                                           pseudo_atom_name, parsed_data["UseInterpolationGrids"].dump()));
    }
    gridPseudoAtomIndices.push_back(static_cast<std::size_t>(std::distance(pseudoAtoms.begin(), match)));
  }

  if (parsed_data.contains("UseRecoilGrowth"))
  {
    useRecoilGrowth = parsed_data["UseRecoilGrowth"].get<bool>();
  }

  if (parsed_data.contains("RecoilGrowthMaximumRecoilLength"))
  {
    recoilGrowthMaximumRecoilLength =
        parsed_data.value("RecoilGrowthMaximumRecoilLength", recoilGrowthMaximumRecoilLength);
  }

  if (parsed_data.contains("RecoilGrowthNumberOfTrialDirections"))
  {
    recoilGrowthNumberOfTrialDirections =
        parsed_data.value("RecoilGrowthNumberOfTrialDirections", recoilGrowthNumberOfTrialDirections);
  }

  if (recoilGrowthMaximumRecoilLength < 1 || recoilGrowthNumberOfTrialDirections < 1)
  {
    throw std::runtime_error(
        std::format("[ForceField reader]: 'RecoilGrowthMaximumRecoilLength' ({}) and "
                    "'RecoilGrowthNumberOfTrialDirections' ({}) must both be at least 1\n",
                    recoilGrowthMaximumRecoilLength, recoilGrowthNumberOfTrialDirections));
  }

  // The recoil feelers are exhaustive depth-first searches: every direction that tests open at a step
  // is probed by up to k^(l-1) trial placements, each of which runs the full base sampler and the
  // torsion selection. The cost per growth step therefore scales as k^l; beyond l = 2 it grows fast and
  // silently, so say so at parse time.
  if (useRecoilGrowth && recoilGrowthMaximumRecoilLength >= 3)
  {
    double feelerCost = 1.0;
    for (std::size_t i = 0; i != recoilGrowthMaximumRecoilLength; ++i)
    {
      feelerCost *= static_cast<double>(recoilGrowthNumberOfTrialDirections);
    }
    std::print(std::cerr,
               "[ForceField reader]: warning: recoil growth with 'RecoilGrowthMaximumRecoilLength' {} and "
               "'RecoilGrowthNumberOfTrialDirections' {} probes up to k^l = {:g} trial placements per growth step "
               "(each a full base-conformation sample plus a torsion selection); the cost grows exponentially "
               "in the recoil length. l = 2 is usually sufficient.\n",
               recoilGrowthMaximumRecoilLength, recoilGrowthNumberOfTrialDirections, feelerCost);
  }

  if (parsed_data.contains("NumberOfTrialDirections"))
  {
    numberOfTrialDirections = parsed_data.value("NumberOfTrialDirections", numberOfTrialDirections);
  }

  if (parsed_data.contains("NumberOfTorsionTrialDirections"))
  {
    numberOfTorsionTrialDirections =
        parsed_data.value("NumberOfTorsionTrialDirections", numberOfTorsionTrialDirections);
  }

  if (parsed_data.contains("NumberOfFirstBeadPositions"))
  {
    numberOfFirstBeadPositions = parsed_data.value("NumberOfFirstBeadPositions", numberOfFirstBeadPositions);
  }

  if (numberOfTrialDirections < 1 || numberOfTorsionTrialDirections < 1 || numberOfFirstBeadPositions < 1)
  {
    throw std::runtime_error(
        std::format("[ForceField reader]: 'NumberOfTrialDirections' ({}), 'NumberOfTorsionTrialDirections' ({}) "
                    "and 'NumberOfFirstBeadPositions' ({}) must all be at least 1\n",
                    numberOfTrialDirections, numberOfTorsionTrialDirections, numberOfFirstBeadPositions));
  }

  if (parsed_data.contains("NumberOfTrialMovesPerOpenBead"))
  {
    numberOfTrialMovesPerOpenBead = parsed_data.value("NumberOfTrialMovesPerOpenBead", numberOfTrialMovesPerOpenBead);
  }

  if (parsed_data.contains("UseDualCutOff"))
  {
    useDualCutOff = parsed_data["UseDualCutOff"].get<bool>();
  }

  if (parsed_data.contains("DualCutOff"))
  {
    dualCutOff = parsed_data.value("DualCutOff", dualCutOff);
  }

  if (useDualCutOff)
  {
    // The inner cut-off must lie strictly inside every explicitly set full cut-off, otherwise the
    // 'correction' from the inner to the full cut-offs is not a correction at all. Automatic cut-offs are
    // only known once the system is built and are not checked here.
    if (dualCutOff <= 0.0 || dualCutOff >= smallestExplicitCutOff)
    {
      throw std::runtime_error(std::format(
          "[ForceField reader]: 'DualCutOff' {} must be positive and smaller than every full cut-off (smallest: {})\n",
          dualCutOff, smallestExplicitCutOff));
    }
  }

  if (parsed_data.contains("CBMCRingCrankshaftProbability"))
  {
    cbmcRingCrankshaftProbability = parsed_data.value("CBMCRingCrankshaftProbability", cbmcRingCrankshaftProbability);
    if (cbmcRingCrankshaftProbability < 0.0 || cbmcRingCrankshaftProbability > 1.0)
    {
      throw std::runtime_error(std::format("[ForceField reader]: 'CBMCRingCrankshaftProbability' {} not in [0, 1]\n",
                                           cbmcRingCrankshaftProbability));
    }
  }

  if (parsed_data.contains("CBMCRingTiltProbability"))
  {
    cbmcRingTiltProbability = parsed_data.value("CBMCRingTiltProbability", cbmcRingTiltProbability);
    if (cbmcRingTiltProbability < 0.0 || cbmcRingTiltProbability > 1.0)
    {
      throw std::runtime_error(
          std::format("[ForceField reader]: 'CBMCRingTiltProbability' {} not in [0, 1]\n", cbmcRingTiltProbability));
    }
  }

  if (parsed_data.contains("UseExternalFieldGrid"))
  {
    useExternalFieldGrid = parsed_data.value("UseExternalFieldGrid", parsed_data["UseExternalFieldGrid"]);
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
    numberOfVDWGridPoints = parseUint3Setting("NumberOfVDWGridPoints", parsed_data["NumberOfVDWGridPoints"]);
  }
  if (parsed_data.contains("NumberOfCoulombGridPoints"))
  {
    numberOfCoulombGridPoints =
        parseUint3Setting("NumberOfCoulombGridPoints", parsed_data["NumberOfCoulombGridPoints"]);
  }
  if (parsed_data.contains("NumberOfExternalFieldGridPoints"))
  {
    numberOfExternalFieldGridPoints =
        parseUint3Setting("NumberOfExternalFieldGridPoints", parsed_data["NumberOfExternalFieldGridPoints"]);
  }

  if (parsed_data.contains("NumberOfGridTestPoints"))
  {
    numberOfGridTestPoints = parsed_data.value("NumberOfGridTestPoints", parsed_data["NumberOfGridTestPoints"]);
  }

  if (parsed_data.contains("WriteExternalFieldInterpolationGrid"))
  {
    writeExternalFieldInterpolationGrid = parsed_data["WriteExternalFieldInterpolationGrid"].get<bool>();
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

std::string ForceFieldSettings::printSamplingStatus() const
{
  std::ostringstream stream;

  std::print(stream, "CBMC first-bead trial positions:     {}\n", numberOfFirstBeadPositions);
  std::print(stream, "CBMC trial directions:               {}\n", numberOfTrialDirections);
  std::print(stream, "CBMC torsion trial directions:       {}\n", numberOfTorsionTrialDirections);
  std::print(stream, "CBMC trial moves per open bead:      {}\n", numberOfTrialMovesPerOpenBead);
  std::print(stream, "CBMC ring crankshaft probability:    {:g}\n", cbmcRingCrankshaftProbability);
  std::print(stream, "CBMC ring tilt probability:          {:g}\n", cbmcRingTiltProbability);
  std::print(stream, "CBMC dual cut-off:                   {}\n", useDualCutOff ? "yes" : "no");
  if (useDualCutOff)
  {
    std::print(stream, "CBMC inner cut-off:                 {:9.5f} [{}]\n", dualCutOff,
               Units::displayedUnitOfLengthString);
  }
  std::print(stream, "Chain growth scheme:                 {}\n",
             useRecoilGrowth ? "recoil growth" : "configurational bias");
  if (useRecoilGrowth)
  {
    std::print(stream, "Recoil-growth trial directions (k):  {}\n", recoilGrowthNumberOfTrialDirections);
    std::print(stream, "Recoil-growth recoil length (l):     {}\n", recoilGrowthMaximumRecoilLength);
  }
  std::print(stream, "\n");

  return stream.str();
}

std::string ForceFieldSettings::printInterpolationGridStatus() const
{
  if (gridPseudoAtomIndices.empty()) return {};

  std::ostringstream stream;

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

  return stream.str();
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ForceFieldSettings& s)
{
  archive << s.versionNumber;

  archive << s.numberOfTrialDirections;
  archive << s.numberOfTorsionTrialDirections;
  archive << s.numberOfFirstBeadPositions;
  archive << s.numberOfTrialMovesPerOpenBead;
  archive << s.minimumRosenbluthFactor;
  archive << s.cbmcRingCrankshaftProbability;
  archive << s.cbmcRingTiltProbability;

  archive << s.useRecoilGrowth;
  archive << s.recoilGrowthMaximumRecoilLength;
  archive << s.recoilGrowthNumberOfTrialDirections;

  archive << s.useDualCutOff;
  archive << s.dualCutOff;

  archive << s.gridPseudoAtomIndices;
  archive << s.spacingVDWGrid;
  archive << s.spacingCoulombGrid;
  archive << s.numberOfVDWGridPoints;
  archive << s.numberOfCoulombGridPoints;
  archive << s.numberOfGridTestPoints;
  archive << s.interpolationSchemeAuto;
  archive << s.interpolationScheme;
  archive << s.writeFrameworkInterpolationGrids;

  archive << s.useExternalFieldGrid;
  archive << s.numberOfExternalFieldGridPoints;
  archive << s.writeExternalFieldInterpolationGrid;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ForceFieldSettings& s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > s.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'ForceFieldSettings' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> s.numberOfTrialDirections;
  archive >> s.numberOfTorsionTrialDirections;
  archive >> s.numberOfFirstBeadPositions;
  archive >> s.numberOfTrialMovesPerOpenBead;
  archive >> s.minimumRosenbluthFactor;
  archive >> s.cbmcRingCrankshaftProbability;
  archive >> s.cbmcRingTiltProbability;

  archive >> s.useRecoilGrowth;
  archive >> s.recoilGrowthMaximumRecoilLength;
  archive >> s.recoilGrowthNumberOfTrialDirections;

  archive >> s.useDualCutOff;
  archive >> s.dualCutOff;

  archive >> s.gridPseudoAtomIndices;
  archive >> s.spacingVDWGrid;
  archive >> s.spacingCoulombGrid;
  archive >> s.numberOfVDWGridPoints;
  archive >> s.numberOfCoulombGridPoints;
  archive >> s.numberOfGridTestPoints;
  archive >> s.interpolationSchemeAuto;
  archive >> s.interpolationScheme;
  archive >> s.writeFrameworkInterpolationGrids;

  archive >> s.useExternalFieldGrid;
  archive >> s.numberOfExternalFieldGridPoints;
  archive >> s.writeExternalFieldInterpolationGrid;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("ForceFieldSettings: Error in binary restart\n"));
  }
#endif

  return archive;
}
