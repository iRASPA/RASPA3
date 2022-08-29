module;

module forcefield;

import skelement;
import units;
import double4;
import print;

import <filesystem>;
import <fstream>;
import <cstdlib>;
import <iostream>;
import <sstream>;
import <format>;
import <exception>;
import <vector>;
import <map>;
import <cmath>;
import <string>;
import <string_view>;


ForceField::ForceField(std::string pseudoAtomsFileName,
    std::string forceFieldMixingFileName, [[maybe_unused]] std::string forceFieldOverwriteFileName) noexcept(false)
{
  ReadPseudoAtoms(pseudoAtomsFileName);
  ReadForceFieldMixing(forceFieldMixingFileName);
}

void ForceField::ReadPseudoAtoms(std::string pseudoAtomsFileName) noexcept(false)
{
  const char* env_p = std::getenv("RASPA_DIR");
  if (!env_p) throw std::runtime_error("Environment variable 'RASPA_DIR' not set");

  std::filesystem::path pseudoAtomsPathfile = std::filesystem::path(pseudoAtomsFileName);
  if (!std::filesystem::exists(pseudoAtomsPathfile)) pseudoAtomsPathfile = std::filesystem::path(env_p) / pseudoAtomsFileName;

  if (!std::filesystem::exists(pseudoAtomsPathfile)) throw std::runtime_error("File 'pseudo_atoms.def' not found");

  std::ifstream pseudoAtomsFile{ pseudoAtomsPathfile };
  if (!pseudoAtomsFile) throw std::runtime_error("File 'pseudo_atoms.def' exists, but error opening file");

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

    pseudoAtoms.emplace_back(PseudoAtom{ .name = name, .mass = std::stod(mass), .charge = std::stod(charge), .atomicNumber = atomicNumber });
  }

  data = std::vector<VDWParameters>(numberOfPseudoAtoms * numberOfPseudoAtoms);
}

void ForceField::ReadForceFieldMixing(std::string forceFieldMixingFileName) noexcept(false)
{
  const char* env_p = std::getenv("RASPA_DIR");

  std::filesystem::path forceFieldPathfile = std::filesystem::path(forceFieldMixingFileName);
  if (!std::filesystem::exists(forceFieldPathfile)) forceFieldPathfile = std::filesystem::path(env_p) / forceFieldMixingFileName;

  if (!std::filesystem::exists(forceFieldPathfile)) throw std::runtime_error("File 'force_field_mixing_rules.def' not found");

  std::ifstream forceFieldFile{ forceFieldPathfile };
  if (!forceFieldFile) throw std::runtime_error("File 'force_field_mixing_rules.def' exists, but error opening file");

  std::string str{};

  //skip comment line
  std::getline(forceFieldFile, str);

  // read shifted or trunacted
  std::getline(forceFieldFile, str);
  //std::istringstream my_stream(str);

  //skip comment line
  std::getline(forceFieldFile, str);

  // read tail-corrections yes/no
  std::getline(forceFieldFile, str);
  //std::istringstream my_stream(str);

  //skip comment line
  std::getline(forceFieldFile, str);

  // read number of self-interactions
  size_t numberOfSelfInteractions{ 0 };
  std::getline(forceFieldFile, str);
  std::istringstream numberOfSelfInteractionStream(str);
  numberOfSelfInteractionStream >> numberOfSelfInteractions;
  if (numberOfSelfInteractions < 0 || numberOfSelfInteractions>10000) throw std::runtime_error("Incorrect amount of self-interactions");


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
      double mix0 = std::sqrt (data[i * numberOfPseudoAtoms + i].parameters.x * data[j * numberOfPseudoAtoms + j].parameters.x);
      double mix1 = 0.5 * (data[i * numberOfPseudoAtoms + i].parameters.y + data[j * numberOfPseudoAtoms + j].parameters.y);
      
      data[i * numberOfPseudoAtoms + j] = VDWParameters(mix0, mix1);
      data[j * numberOfPseudoAtoms + i] = VDWParameters(mix0, mix1);
    }
  }
}


void ForceField::printPseudoAtomStatus(std::ostream &stream) const
{
  std::print(stream, "Pseudo-atoms\n");
  std::print(stream, "===============================================================================\n\n");

  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
      std::print(stream, "{:3d} - {:8} mass: {:8.5f}, charge: {:8.5f}\n", i, pseudoAtoms[i].name, pseudoAtoms[i].mass, pseudoAtoms[i].charge);
  }
  std::print(stream, "\n");
}

void ForceField::printForceFieldStatus(std::ostream &stream) const
{
  std::print(stream, "Force field status\n");
  std::print(stream, "===============================================================================\n\n");

  double conv = Units::EnergyToKelvin;
  for (size_t i = 0; i < numberOfPseudoAtoms; ++i)
  {
    for (size_t j = i; j < numberOfPseudoAtoms; ++j)
    {
      switch (data[i * numberOfPseudoAtoms + j].type)
      {
      case VDWParameters::Type::LennardJones:
        std::print(stream, "{:8} - {:8} {} p_0/k_B: {:8.5f} [K], p_1 {:8.5f} [A]\n",
            pseudoAtoms[i].name, pseudoAtoms[j].name, "Lennard-Jones",
            conv * data[i * numberOfPseudoAtoms + j].parameters.x,
            data[i * numberOfPseudoAtoms + j].parameters.y);
        std::print(stream, "{:33} shift: {:8.5f} [K], tailcorrections: {}\n",
            "",
            conv * data[i * numberOfPseudoAtoms + j].shift,
            data[i * numberOfPseudoAtoms + j].tailCorrection);
        break;
      default:
        break;
      }
    }
  }
  std::print(stream, "\n");
}
