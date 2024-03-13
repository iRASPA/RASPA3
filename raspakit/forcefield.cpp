
module;

module forcefield;

import <filesystem>;
import <fstream>;
import <cstdlib>;
import <iostream>;
import <sstream>;
import <vector>;
import <map>;
import <cmath>;
import <string>;
import <string_view>;
import <optional>;
import <numbers>;
import <algorithm>;
import <print>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;

import skelement;
import units;
import int3;
import double3;
import double4;
import stringutils;


ForceField::ForceField(std::vector<PseudoAtom> pseudoAtoms, std::vector<VDWParameters> selfInteractions, 
                       [[maybe_unused]] MixingRule mixingRule, double cutOff, bool shifted, 
                       bool applyTailCorrections) noexcept(false) :
    data(selfInteractions.size()* selfInteractions.size(), VDWParameters(0.0, 0.0)),
    shiftPotentials(selfInteractions.size()* selfInteractions.size(), shifted),
    tailCorrections(selfInteractions.size()* selfInteractions.size(), applyTailCorrections),
    cutOffVDW(cutOff), 
    numberOfPseudoAtoms(pseudoAtoms.size()), 
    pseudoAtoms(pseudoAtoms)
{
  for (size_t i = 0; i < selfInteractions.size(); ++i)
  {
    data[i + i * numberOfPseudoAtoms] = selfInteractions[i];
  }

  for (size_t i = 0; i < selfInteractions.size(); ++i)
  {
    for (size_t j = i + 1; j < selfInteractions.size(); ++j)
    {
      double mix0 = std::sqrt(selfInteractions[i].parameters.x * selfInteractions[j].parameters.x);
      double mix1 = 0.5 * (selfInteractions[i].parameters.y + selfInteractions[j].parameters.y);

      data[i + numberOfPseudoAtoms * j] = VDWParameters(mix0, mix1);
      data[j + numberOfPseudoAtoms * i] = VDWParameters(mix0, mix1);
    }
  }

  for (size_t i = 0; i < selfInteractions.size(); ++i)
  {
    for (size_t j = 0; j < selfInteractions.size(); ++j)
    {
      if(shiftPotentials[i * numberOfPseudoAtoms + j])
      {
        data[i * numberOfPseudoAtoms + j].computeShiftAtCutOff(cutOffVDW);
      }
    }
  }

  preComputeTailCorrection();
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


ForceField::ForceField(size_t systemId) noexcept(false)
{
  const char* env_p = std::getenv("RASPA_DIR");
  const std::string ext{".def"};

  const std::string defaultPseudoAtomsFileName{"pseudo_atoms" + ext};
  const std::string defaultForceFieldMixingRulesFileName{"force_field_mixing_rules" + ext};
  const std::string defaultForceFieldOverwriteFileName{"force_field" + ext};

  std::string pseudoAtomsFileName{std::format("{}_{}{}", "pseudo_atoms", systemId, ext)};
  std::string forceFieldMixingRulesFileName{std::format("{}_{}{}", "force_field_mixing_rules", systemId, ext)};
  std::string forceFieldOverwriteFileName{std::format("{}_{}{}", "force_field", systemId, ext)};

  if(!std::filesystem::exists(std::filesystem::path{pseudoAtomsFileName}))
  {
    pseudoAtomsFileName = defaultPseudoAtomsFileName;
  }
  if(!std::filesystem::exists(std::filesystem::path{forceFieldMixingRulesFileName}))
  {
    forceFieldMixingRulesFileName = defaultForceFieldMixingRulesFileName;
  }
  if(!std::filesystem::exists(std::filesystem::path{forceFieldOverwriteFileName}))
  {
    forceFieldOverwriteFileName = defaultForceFieldOverwriteFileName;
  }

  if (env_p) 
  {
    if(!std::filesystem::exists(std::filesystem::path{pseudoAtomsFileName}))
    {
      pseudoAtomsFileName = env_p + defaultPseudoAtomsFileName;
    }
    if(!std::filesystem::exists(std::filesystem::path{forceFieldMixingRulesFileName}))
    {
      forceFieldMixingRulesFileName = env_p + defaultForceFieldMixingRulesFileName;
    }
    if(!std::filesystem::exists(std::filesystem::path{forceFieldOverwriteFileName}))
    {
      forceFieldOverwriteFileName = env_p + defaultForceFieldOverwriteFileName;
    }
  }

  std::cout << "Reading force field: " << pseudoAtomsFileName << " " << forceFieldMixingRulesFileName << std::endl;
  ReadPseudoAtoms(pseudoAtomsFileName);
  ReadForceFieldMixing(forceFieldMixingRulesFileName);
  preComputeTailCorrection();
}

void ForceField::ReadPseudoAtoms(std::string pseudoAtomsFileName) noexcept(false)
{
  std::filesystem::path pseudoAtomsPathfile = std::filesystem::path(pseudoAtomsFileName);
  if (!std::filesystem::exists(pseudoAtomsPathfile)) 
  {
    throw std::runtime_error(std::format("File '{}' not found", pseudoAtomsFileName));
  }

  std::ifstream pseudoAtomsFile{ pseudoAtomsPathfile };
  if (!pseudoAtomsFile) 
  {
    throw std::runtime_error(std::format("File '{}' exists, but error opening file", pseudoAtomsFileName));
  }

  std::string str{};

  //skip comment line
  std::getline(pseudoAtomsFile, str);

  // read number of pseudo-atoms
  size_t n;
  std::getline(pseudoAtomsFile, str);
  std::istringstream nuberOfAtomsstream(str);
  nuberOfAtomsstream >> n;
  if (n < 0 || n>10000) throw std::runtime_error("Incorrect amount of pseudo=atoms");
  numberOfPseudoAtoms = n;
  pseudoAtoms = std::vector< PseudoAtom>{};

  //skip comment line
  std::getline(pseudoAtomsFile, str);

  for (size_t i = 0; i < numberOfPseudoAtoms; i++)
  {
    std::string name, print, as, chem, oxidation, mass, charge;

    std::getline(pseudoAtomsFile, str);
    std::istringstream atomStream(str);

    atomStream >> name;       // read name
    atomStream >> print;      // read print
    atomStream >> as;         // read as
    atomStream >> chem;       // read name
    atomStream >> oxidation;  // read name
    atomStream >> mass;       // read name
    atomStream >> charge;     // read charge

    size_t atomicNumber{ 1 };
    auto it = PredefinedElements::atomicNumberData.find(chem);
    // FIX
    if (it != PredefinedElements::atomicNumberData.end()) atomicNumber = static_cast<size_t>(it->second);

    pseudoAtoms.emplace_back(PseudoAtom(name, std::stod(mass), std::stod(charge), atomicNumber, true));
  }

  data = std::vector<VDWParameters>(numberOfPseudoAtoms * numberOfPseudoAtoms);
  tailCorrections = std::vector<bool>(numberOfPseudoAtoms * numberOfPseudoAtoms, false);
  shiftPotentials = std::vector<bool>(numberOfPseudoAtoms * numberOfPseudoAtoms, true);
}

void ForceField::ReadForceFieldMixing(std::string forceFieldMixingFileName) noexcept(false)
{
  std::filesystem::path forceFieldPathfile = std::filesystem::path(forceFieldMixingFileName);
  if (!std::filesystem::exists(forceFieldPathfile)) 
  {
    throw std::runtime_error(std::format("File '{}' not found", forceFieldMixingFileName));
  }

  std::ifstream forceFieldFile{ forceFieldPathfile };
  if (!forceFieldFile) 
  {
    throw std::runtime_error(std::format("File '{}' exists, but error opening file", forceFieldMixingFileName));
  }

  std::string str{};

  //skip comment line
  std::getline(forceFieldFile, str);

  // read shifted or trunacted
  std::getline(forceFieldFile, str);
  
  if (caseInSensStringCompare(str, "truncated"))
  {
    std::fill(shiftPotentials.begin(), shiftPotentials.end(), false);
  }
  if (caseInSensStringCompare(str, "shifted"))
  {
    std::fill(shiftPotentials.begin(), shiftPotentials.end(), true);
  }

  //skip comment line
  std::getline(forceFieldFile, str);

  // read tail-corrections yes/no
  std::getline(forceFieldFile, str);
  if (caseInSensStringCompare(str, "yes"))
  {
    std::fill(tailCorrections.begin(), tailCorrections.end(), true);
  }
  if (caseInSensStringCompare(str, "no"))
  {
    std::fill(tailCorrections.begin(), tailCorrections.end(), false);
  }

  //skip comment line
  std::getline(forceFieldFile, str);

  // read number of self-interactions
  size_t numberOfSelfInteractions{ 0 };
  std::getline(forceFieldFile, str);
  std::istringstream numberOfSelfInteractionStream(str);
  numberOfSelfInteractionStream >> numberOfSelfInteractions;
  if (numberOfSelfInteractions < 0 || numberOfSelfInteractions>10000) 
    throw std::runtime_error("Incorrect amount of self-interactions");


  //skip comment line
  std::getline(forceFieldFile, str);

  for (size_t i = 0; i < numberOfSelfInteractions; i++)
  {
    std::string name, ffType, param[10];

    std::getline(forceFieldFile, str);
    std::istringstream my_stream(str);

    my_stream >> name;        // read name
    my_stream >> ffType;      // read VDW-type
    my_stream >> param[0];    // read epsilon
    my_stream >> param[1];    // read sigma

    data[i * numberOfPseudoAtoms + i] = VDWParameters(std::stod(param[0]) * Units::KelvinToEnergy, std::stod(param[1]));
  }

  for (size_t i = 0; i < numberOfSelfInteractions; ++i)
  {
    for (size_t j = i + 1; j < numberOfSelfInteractions; ++j)
    {
      double mix0 = 
        std::sqrt(data[i * numberOfPseudoAtoms + i].parameters.x * data[j * numberOfPseudoAtoms + j].parameters.x);
      double mix1 = 
        0.5 * (data[i * numberOfPseudoAtoms + i].parameters.y + data[j * numberOfPseudoAtoms + j].parameters.y);
      
      data[i * numberOfPseudoAtoms + j] = VDWParameters(mix0, mix1);
      data[j * numberOfPseudoAtoms + i] = VDWParameters(mix0, mix1);
    }
  }

  for (size_t i = 0; i < data.size(); ++i)
  {
    if (shiftPotentials[i])
    {
      data[i].computeShiftAtCutOff(cutOffVDW);
    }
  }
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
  std::vector<PseudoAtom>::const_iterator it = std::find_if(
        pseudoAtoms.begin(), pseudoAtoms.end(),
        [&name](const PseudoAtom& x) { return x.name == name; });
  if (it != std::end(pseudoAtoms))
  {
    return (it - pseudoAtoms.begin());
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

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VDWParameters &p)
{
  archive << p.parameters;
  archive << p.shift;
  archive << p.tailCorrectionEnergy;
  archive << p.type;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VDWParameters &p)
{
  archive >> p.parameters;
  archive >> p.shift;
  archive >> p.tailCorrectionEnergy;
  archive >> p.type;

  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PseudoAtom &a)
{
  archive << a.versionNumber;
  archive << a.name;
  archive << a.mass;
  archive << a.charge;
  archive << a.atomicNumber;
  archive << a.printToPDB;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PseudoAtom &a)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > a.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PseudoAtom' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> a.name;
  archive >> a.mass;
  archive >> a.charge;
  archive >> a.atomicNumber;
  archive >> a.printToPDB;

  return archive;
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
