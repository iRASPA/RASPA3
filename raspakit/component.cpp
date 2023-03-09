module;

module component;

import int3;
import double3;
import double3x3;
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
import lambda;
import print;
import property_widom;
import multi_site_isotherm;
import simulationbox;
import cif_reader;
import move_statistics;
import bond_potential;

import <iostream>;
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
#if defined(_WIN32)
  import <cassert>;
#else 
  #include <assert.h>
#endif

Component::Component(size_t componentId, std::string componentName, double mass, SimulationBox simulationBox, double T_c, double P_c, double w,
    std::vector<Atom> definedAtoms, size_t numberOfBlocks) noexcept(false) :
    type(Type::Adsorbate),
    simulationBox(simulationBox),
    componentId(componentId),
    name(componentName),
    criticalTemperature(T_c),
    criticalPressure(P_c),
    acentricFactor(w),
    mass(mass),
    definedAtoms(definedAtoms),
    lambda(numberOfBlocks, 41),
    averageRosenbluthWeights(numberOfBlocks)
{
  atoms = definedAtoms;
}

Component::Component(size_t componentId, std::string fileName, double mass, SimulationBox simulationBox, size_t spaceGroupHallNumber, std::vector<Atom> definedAtoms, int3 numberOfUnitCells, size_t numberOfBlocks) noexcept(false) :
    type(Type::Framework),
    simulationBox(simulationBox),
    spaceGroupHallNumber(spaceGroupHallNumber),
    numberOfUnitCells(numberOfUnitCells),
    componentId(componentId),
    name(fileName),
    filenameData(fileName),
    mass(mass),
    definedAtoms(definedAtoms),
    lambda(numberOfBlocks, 41),
    averageRosenbluthWeights(numberOfBlocks)
{
  // expand the fractional atoms based on the space-group
  SKSpaceGroup spaceGroup = SKSpaceGroup(spaceGroupHallNumber);
  std::vector<Atom> expandedAtoms;
  expandedAtoms.reserve(definedAtoms.size() * 256);
  for (const Atom& atom : definedAtoms)
  {
    Atom atomCopy = atom;

    std::vector<double3> listOfPositions = spaceGroup.listOfSymmetricPositions(atom.position);
    for (const double3& pos : listOfPositions)
    {
      atomCopy.position = simulationBox.unitCell * pos.fract();
      expandedAtoms.push_back(atomCopy);
    }
  }

  // eliminate duplicates
  std::vector<Atom> unitCellAtoms;
  for (size_t i = 0; i < expandedAtoms.size(); ++i)
  {
    bool overLap = false;
    for (size_t j = i + 1; j < expandedAtoms.size(); ++j)
    {
      double3 dr = expandedAtoms[i].position - expandedAtoms[j].position;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double rr = double3::dot(dr, dr);
      if (rr < 0.1)
      {
        overLap = true;
        break;
      }
    }
    if (!overLap)
    {
      unitCellAtoms.push_back(expandedAtoms[i]);
    }
  }



  for (int32_t i = 0; i < numberOfUnitCells.x; ++i)
  {
    for (int32_t j = 0; j < numberOfUnitCells.y; ++j)
    {
      for (int32_t k = 0; k < numberOfUnitCells.z; ++k)
      {
        for (const Atom& atom : unitCellAtoms)
        {
          Atom atomCopy = atom;
          atomCopy.position += simulationBox.unitCell * double3(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
          atoms.push_back(atomCopy);
        }
      }
    }
  }


  for (size_t i = 0; i < atoms.size(); ++i)
  {
    atoms[i].componentId = static_cast<short>(componentId);
    atoms[i].moleculeId = 0;
  }
}

Component::Component(Component::Type type, size_t currentComponent, const std::string &componentName,
                     std::optional<const std::string> fileName, size_t numberOfBlocks) noexcept(false) :
                     type(type), 
                     componentId(currentComponent), 
                     name(componentName),
                     filenameData(fileName),
                     lambda(numberOfBlocks, 41),
                     averageRosenbluthWeights(numberOfBlocks)
{
}

template<typename T>
std::vector<T> parseListOfParameters(const std::string& arguments, size_t lineNumber)
{
    std::vector<T> list{};

    std::string str;
    std::istringstream ss(arguments);

    while (ss >> str)
    {
        if (trim(str).rfind("//", 0) == 0)
        {
            if (list.empty())
            {
                throw std::runtime_error(std::print("No values could be read at line: {}", lineNumber));
            }
            return list;
        }
        T value;
        std::istringstream s(str);
        if (s >> value)
        {
            list.push_back(value);
        }
        else
        {
            if (list.empty())
            {
                throw std::runtime_error(std::print("No values could be read at line: {}", lineNumber));
            }
            return list;
        }

    };

    return list;
}

// TODO: extend to flexible molecules
void Component::readComponent(const ForceField& forceField, const std::string& fileName)
{
  const char* env_p = std::getenv("RASPA_DIR");
  const std::string& moleculeFileName = fileName + ".def";

  std::filesystem::path moleculePathfile = std::filesystem::path(moleculeFileName);
  if (!std::filesystem::exists(moleculePathfile)) moleculePathfile = std::filesystem::path(env_p) / moleculeFileName;

  if (!std::filesystem::exists(moleculePathfile)) throw std::runtime_error(std::print("File '{}' not found", moleculeFileName));

  std::ifstream moleculeFile{ moleculePathfile };
  if (!moleculeFile) throw std::runtime_error(std::print("[Component] File '{}' exists, but error opening file", moleculeFileName));

  std::string str{};
  //size_t lineNumber{ 0 };

  // skip comment line
  std::getline(moleculeFile, str);
  
  // read critical temperature
  std::getline(moleculeFile, str);
  std::istringstream critical_temperature_stream(str);
  critical_temperature_stream >> criticalTemperature;
  if (criticalTemperature < 0.0) throw std::runtime_error("Incorrect critical temperature");

  // read critical pressure
  std::getline(moleculeFile, str);
  std::istringstream critical_pressure_stream(str);
  critical_pressure_stream >> criticalPressure;
  if (criticalTemperature < 0.0) throw std::runtime_error("Incorrect critical pressure");

  // read acentric factor
  std::getline(moleculeFile, str);
  std::istringstream acentric_factor_stream(str);
  acentric_factor_stream >> acentricFactor;
  

  // skip comment line
  std::getline(moleculeFile, str);
  
  // read number of pseudo-atoms
  size_t n;
  std::getline(moleculeFile, str);
  std::istringstream my_stream(str);
  my_stream >> n;
  if (n < 0 || n>10000) throw std::runtime_error("Incorrect amount of pseudo=atoms");

  definedAtoms.resize(n);

  // skip comment line
  std::getline(moleculeFile, str);
  
  // read Rigid or Flexible
  std::getline(moleculeFile, str);
  std::istringstream rigidStream(str);
  std::string rigidString;
  rigidStream >> rigidString;
  if (caseInSensStringCompare(str, "flexible"))
  {
    rigid = false;
  }
  
  // skip comment line
  std::getline(moleculeFile, str);

  // read atomic positions
  this->mass = 0.0;
  for (size_t i = 0; i < definedAtoms.size(); ++i)
  {
      int id;
      std::string atomTypeString;
      double3 pos;
      std::getline(moleculeFile, str);
      std::istringstream atomStream(str);
      atomStream >> id >> atomTypeString >> pos.x >> pos.y >> pos.z;

      // find atom-type based on read 'atomTypeString'
      auto it = std::find_if(forceField.pseudoAtoms.begin(), forceField.pseudoAtoms.end(), 
                             [&](const PseudoAtom &atom) 
                             {
                               return atomTypeString == atom.name; 
                             });
      
      if (it == forceField.pseudoAtoms.end())
      {
        throw std::runtime_error(std::print("readComponent: Atom-string '{}' not found (define them first in 'the pseudo_atoms.def' file)", atomTypeString));
      }
      
      size_t pseudoAtomType = static_cast<size_t>(std::distance(forceField.pseudoAtoms.begin(), it));
      this->mass += forceField.pseudoAtoms[pseudoAtomType].mass;
      double charge = forceField.pseudoAtoms[pseudoAtomType].charge;
      double scaling = 1.0;

      definedAtoms[i] = Atom(pos, charge, scaling, static_cast<short>(pseudoAtomType), static_cast<short>(componentId), 0);
  }

  atoms = definedAtoms;

  // skip comment line
  std::getline(moleculeFile, str);

  // read number of intra-molecular interactions
  std::getline(moleculeFile, str);
  std::istringstream numberOfInteractionsString(str);
  size_t numberOfChiralCenters, numberOfBonds, numberOfBondDipoles, numberOfBends, numberOfUreyBradleys;
  numberOfInteractionsString >> numberOfChiralCenters >> numberOfBonds >> numberOfBondDipoles >> numberOfBends >> numberOfUreyBradleys;

  if(numberOfBonds > 0)
  {
    // skip comment line
    std::getline(moleculeFile, str);

    bonds.resize(numberOfBonds);
    for (size_t i = 0; i < numberOfBonds; ++i)
    {
      size_t idA, idB;
      std::string bondTypeString, parameterString;
      std::getline(moleculeFile, str);
      std::istringstream bondTypeStream(str);

      bondTypeStream >> idA >> idB >> bondTypeString;
      std::getline(bondTypeStream, parameterString);

      std::vector<double> parameters = parseListOfParameters<double>(parameterString, 0);

      bonds[i] = BondPotential(BondPotential::bondDefinitionForString[bondTypeString], std::make_pair(idA, idB));
      std::copy(parameters.begin(), parameters.end(), bonds[i].parameters.begin());
    }
  }
}

void Component::readFramework([[maybe_unused]] const ForceField& forceField, [[maybe_unused]] const std::string& fileName)
{
  const char* env_p = std::getenv("RASPA_DIR");

  const std::string frameworkFileName = fileName + ".cif";

  std::filesystem::path frameworkPathfile = std::filesystem::path(frameworkFileName);
  if (!std::filesystem::exists(frameworkPathfile)) frameworkPathfile = std::filesystem::path(env_p) / frameworkFileName;

  if (!std::filesystem::exists(frameworkPathfile)) throw std::runtime_error(std::print("File '{}' not found", frameworkFileName));

  std::ifstream t(frameworkPathfile);
  std::string fileContent((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());
 
  CIFReader parser = CIFReader(fileContent, forceField);
  simulationBox = parser.simulationBox;
  definedAtoms = parser.fractionalAtoms;

  // expand the fractional atoms based on the space-group
  size_t number = parser._spaceGroupHallNumber.value_or(1);
  SKSpaceGroup spaceGroup = SKSpaceGroup(number);
  std::vector<Atom> expandedAtoms;
  expandedAtoms.reserve(parser.fractionalAtoms.size() * 256);
  for (const Atom &atom : definedAtoms)
  {
    Atom atomCopy = atom;
    std::vector<double3> listOfPositions = spaceGroup.listOfSymmetricPositions(atom.position);
    for (const double3& pos : listOfPositions)
    {
      atomCopy.position = parser.simulationBox.unitCell * pos.fract();
      expandedAtoms.push_back(atomCopy);
    }
  }

  // eliminate duplicates
  std::vector<Atom> unitCellAtoms;
  for (size_t i = 0; i < expandedAtoms.size(); ++i)
  {
    bool overLap = false;
    for (size_t j = i + 1; j < expandedAtoms.size(); ++j)
    {
        double3 dr = expandedAtoms[i].position - expandedAtoms[j].position;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        double rr = double3::dot(dr, dr);
        if (rr < 0.1)
        {
            overLap = true;
            break;
        }
      
    }
    if (!overLap)
    {
      unitCellAtoms.push_back(expandedAtoms[i]);
    }
  }

  

  for (int32_t i = 0; i < numberOfUnitCells.x; ++i)
  {
      for (int32_t j = 0; j < numberOfUnitCells.y; ++j)
      {
          for (int32_t k = 0; k < numberOfUnitCells.z; ++k)
          {
             for (const Atom& atom : unitCellAtoms)
             {
               Atom atomCopy = atom;
               atomCopy.position += simulationBox.unitCell * double3(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
               atoms.push_back(atomCopy);
             }
          }
      }
  }

  mass = 0.0;
  for (const Atom& atom : atoms)
  {
      size_t atomType = static_cast<size_t>(atom.type);
      mass += forceField.pseudoAtoms[atomType].mass;
  }

  for (size_t i = 0; i < atoms.size(); ++i)
  {
    atoms[i].componentId = static_cast<short>(componentId);
    atoms[i].moleculeId = 0;
  }  
}

void Component::printStatus(std::ostream &stream, const ForceField& forceField) const
{
  std::print(stream, "Component {} [{}]\n\n", componentId, name);

  std::print(stream, "    Critical temperature:  {} [K]\n", criticalTemperature);
  std::print(stream, "    Critical pressure:     {} [Pa]\n", criticalPressure);
  std::print(stream, "    Acentric factor:       {} [-]\n\n", acentricFactor);

  std::print(stream, "    Mol-fraction:                 {} [-]\n", molFraction);
  std::print(stream << std::boolalpha, "    Swapable:                     {}\n\n", swapable);
  std::print(stream, "    Mass:                         {} [-]\n", mass);
  std::print(stream << std::boolalpha, "    Compute fugacity-coefficient: {}\n", computeFugacityCoefficient);
  std::print(stream, "    Fugacity coefficient:         {} [-]\n", fugacityCoefficient);
  std::print(stream, "    Bulk fluid density:           {} [-]\n", bulkFluidDensity);
  std::print(stream, "    Compressibility:              {} [-]\n", compressibility);
  std::print(stream, "    Excess molecules:             {} [-]\n\n", amountOfExcessMolecules);

  std::print(stream, "    number Of Atoms:  {}\n", atoms.size());
  
  for (size_t i = 0; i < definedAtoms.size(); ++i)
  {
    size_t atomType = static_cast<size_t>(definedAtoms[i].type);
    std::string atomTypeString = forceField.pseudoAtoms[atomType].name;
    std::print(stream, "    {:3d}: {:6} position {:8.5f} {:8.5f} {:8.5f}, charge {:8.5f}\n", 
                   i, atomTypeString, definedAtoms[i].position.x, definedAtoms[i].position.y, definedAtoms[i].position.z, definedAtoms[i].charge);
  }
  std::print(stream, "    net-charge:  {}\n", netCharge);
  std::print(stream, "\n");
  
  std::print(stream, "    Translation-move probability:             {} [-]\n", probabilityTranslationMove);
  std::print(stream, "    Random translation-move probability:      {} [-]\n", probabilityRandomTranslationMove);
  std::print(stream, "    Rotation-move probability:                {} [-]\n", probabilityRotationMove);
  std::print(stream, "    Random rotation-move probability:         {} [-]\n", probabilityRandomRotationMove);
  std::print(stream, "    Volume-move probability:                  {} [-]\n", probabilityVolumeMove);
  std::print(stream, "    Reinsertion (CBMC) probability:           {} [-]\n", probabilityReinsertionMove_CBMC);
  std::print(stream, "    Identity-change (CBMC) probability:       {} [-]\n", probabilityIdentityChangeMove_CBMC);
  std::print(stream, "    Swap-move (CBMC) probability:             {} [-]\n", probabilitySwapMove_CBMC);
  std::print(stream, "    Swap-move (CFCMC) probability:            {} [-]\n", probabilitySwapMove_CFCMC);
  std::print(stream, "    Swap-move (CFCMC/CBMC) probability:       {} [-]\n", probabilitySwapMove_CFCMC_CBMC);
  std::print(stream, "    Gibbs Volume-move probability:            {} [-]\n", probabilityGibbsVolumeMove);
  std::print(stream, "    Gibbs Swap-move (CBMC) probability:       {} [-]\n", probabilityGibbsSwapMove_CBMC);
  std::print(stream, "    Gibbs Swap-move (CFCMC) probability:      {} [-]\n", probabilityGibbsSwapMove_CFCMC);
  std::print(stream, "    Gibbs Swap-move (CFCMC/CBMC) probability: {} [-]\n", probabilityGibbsSwapMove_CFCMC_CBMC);
  std::print(stream, "    Widom probability:                        {} [-]\n", probabilityWidomMove);
  std::print(stream, "    Widom (CFCMC) probability:                {} [-]\n", probabilityWidomMove_CFCMC);
  std::print(stream, "    Widom (CFCMC/CBMC) probability:           {} [-]\n", probabilityWidomMove_CFCMC_CBMC);
  std::print(stream, "\n");

  std::print(stream, "    number of bonds: {}\n", bonds.size());
  for (size_t i = 0; i < bonds.size(); ++i)
  {
      std::print(stream, "        {}", bonds[i].print());
  }
  std::print(stream, "\n");
}

void Component::normalizeMoveProbabilties()
{
    double totalProbability = 
        probabilityTranslationMove +
        probabilityRandomTranslationMove +
        probabilityRotationMove +
        probabilityRandomRotationMove +
        probabilityVolumeMove +
        probabilityReinsertionMove_CBMC +
        probabilityIdentityChangeMove_CBMC +
        probabilitySwapMove_CBMC +
        probabilitySwapMove_CFCMC +
        probabilitySwapMove_CFCMC_CBMC +
        probabilityGibbsVolumeMove +
        probabilityGibbsSwapMove_CBMC +
        probabilityGibbsSwapMove_CFCMC +
        probabilityGibbsSwapMove_CFCMC_CBMC +
        probabilityWidomMove +
        probabilityWidomMove_CFCMC +
        probabilityWidomMove_CFCMC_CBMC;


    if (totalProbability > 1e-5)
    {
        probabilityTranslationMove /= totalProbability;
        probabilityRandomTranslationMove /= totalProbability;
        probabilityRotationMove /= totalProbability;
        probabilityRandomRotationMove /= totalProbability;
        probabilityVolumeMove /= totalProbability;
        probabilityReinsertionMove_CBMC /= totalProbability;
        probabilityIdentityChangeMove_CBMC /= totalProbability;
        probabilitySwapMove_CBMC /= totalProbability;
        probabilitySwapMove_CFCMC /= totalProbability;
        probabilitySwapMove_CFCMC_CBMC /= totalProbability;
        probabilityGibbsVolumeMove /=  totalProbability;
        probabilityGibbsSwapMove_CBMC /= totalProbability;
        probabilityGibbsSwapMove_CFCMC /= totalProbability;
        probabilityGibbsSwapMove_CFCMC_CBMC /= totalProbability;
        probabilityWidomMove /= totalProbability;
        probabilityWidomMove_CFCMC /= totalProbability;
        probabilityWidomMove_CFCMC_CBMC /= totalProbability;
    }

    accumulatedProbabilityTranslationMove = probabilityTranslationMove;
    accumulatedProbabilityRandomTranslationMove = probabilityRandomTranslationMove;
    accumulatedProbabilityRotationMove = probabilityRotationMove;
    accumulatedProbabilityRandomRotationMove = probabilityRandomRotationMove;
    accumulatedProbabilityVolumeMove = probabilityVolumeMove;
    accumulatedProbabilityReinsertionMove_CBMC = probabilityReinsertionMove_CBMC;
    accumulatedProbabilityIdentityChangeMove_CBMC = probabilityIdentityChangeMove_CBMC;
    accumulatedProbabilitySwapMove_CBMC = probabilitySwapMove_CBMC;
    accumulatedProbabilitySwapMove_CFCMC = probabilitySwapMove_CFCMC;
    accumulatedProbabilitySwapMove_CFCMC_CBMC = probabilitySwapMove_CFCMC_CBMC;
    accumulatedProbabilityGibbsVolumeMove = probabilityGibbsVolumeMove;
    accumulatedProbabilityGibbsSwapMove_CBMC = probabilityGibbsSwapMove_CBMC;
    accumulatedProbabilityGibbsSwapMove_CFCMC = probabilityGibbsSwapMove_CFCMC;
    accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC = probabilityGibbsSwapMove_CFCMC_CBMC;
    accumulatedProbabilityWidomMove = probabilityWidomMove;
    accumulatedProbabilityWidomMove_CFCMC = probabilityWidomMove_CFCMC;
    accumulatedProbabilityWidomMove_CFCMC_CBMC = probabilityWidomMove_CFCMC_CBMC;

    accumulatedProbabilityRandomTranslationMove += accumulatedProbabilityTranslationMove;
    accumulatedProbabilityRotationMove += accumulatedProbabilityRandomTranslationMove;
    accumulatedProbabilityRandomRotationMove += accumulatedProbabilityRotationMove;
    accumulatedProbabilityVolumeMove += accumulatedProbabilityRandomRotationMove;
    accumulatedProbabilityReinsertionMove_CBMC += accumulatedProbabilityVolumeMove;
    accumulatedProbabilityIdentityChangeMove_CBMC += accumulatedProbabilityReinsertionMove_CBMC;
    accumulatedProbabilitySwapMove_CBMC += accumulatedProbabilityIdentityChangeMove_CBMC;
    accumulatedProbabilitySwapMove_CFCMC += accumulatedProbabilitySwapMove_CBMC;
    accumulatedProbabilitySwapMove_CFCMC_CBMC += accumulatedProbabilitySwapMove_CFCMC;
    accumulatedProbabilityGibbsVolumeMove += accumulatedProbabilitySwapMove_CFCMC_CBMC;
    accumulatedProbabilityGibbsSwapMove_CBMC += accumulatedProbabilityGibbsVolumeMove;
    accumulatedProbabilityGibbsSwapMove_CFCMC += accumulatedProbabilityGibbsSwapMove_CBMC;
    accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC += accumulatedProbabilityGibbsSwapMove_CFCMC;
    accumulatedProbabilityWidomMove += accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC;
    accumulatedProbabilityWidomMove_CFCMC += accumulatedProbabilityWidomMove;
    accumulatedProbabilityWidomMove_CFCMC_CBMC += accumulatedProbabilityWidomMove_CFCMC;

}

std::string formatStatistics(const std::string name, const MoveStatistics<double>& move)
{
    std::ostringstream stream;
    std::print(stream, "    {} total:        {:10}\n", name, move.counts);
    std::print(stream, "    {} constructed:  {:10}\n", name, move.constructed);
    std::print(stream, "    {} accepted:     {:10}\n", name, move.accepted);
    std::print(stream, "    {} fraction:     {:10f}\n", name, move.accepted / std::max(1.0, double(move.counts)));
    std::print(stream, "    {} max-change:   {:10f}\n\n", name, move.maxChange);
    return stream.str();
}

std::string formatStatistics(const std::string name, const MoveStatistics<double3> &move)
{
    std::ostringstream stream;
    std::print(stream, "    {} total:        {:10} {:10} {:10}\n", name, move.counts.x, move.counts.y, move.counts.z);
    std::print(stream, "    {} constructed:  {:10} {:10} {:10}\n", name, move.constructed.x, move.constructed.y, move.constructed.z);
    std::print(stream, "    {} accepted:     {:10} {:10} {:10}\n", name, move.accepted.x, move.accepted.y, move.accepted.z);
    std::print(stream, "    {} fraction:     {:10f} {:10f} {:10f}\n", name, move.accepted.x / std::max(1.0, double(move.counts.x)),
                   move.accepted.y / std::max(1.0, double(move.counts.y)), move.accepted.z / std::max(1.0, double(move.counts.z)));
    std::print(stream, "    {} max-change:   {:10f} {:10f} {:10f}\n\n", name, move.maxChange.x, move.maxChange.y, move.maxChange.z);
    return stream.str();
}

const std::string Component::writeMCMoveStatistics() const
{
    std::ostringstream stream;
    if(probabilityTranslationMove > 0.0) std::print(stream, formatStatistics("Translation", statistics_TranslationMove));
    if(probabilityRandomTranslationMove > 0.0) std::print(stream, formatStatistics("Random translation", statistics_RandomTranslationMove));
    if(probabilityRotationMove > 0.0) std::print(stream, formatStatistics("Rotation", statistics_RotationMove));
    if(probabilityRandomRotationMove > 0.0) std::print(stream, formatStatistics("Random rotation", statistics_RandomRotationMove));
    if(probabilityReinsertionMove_CBMC > 0.0) std::print(stream, formatStatistics("Reinsertion(CBMC)", statistics_ReinsertionMove_CBMC));
    if(probabilityIdentityChangeMove_CBMC > 0.0) std::print(stream, formatStatistics("Identity Swap (CBMC)", statistics_IdentityChangeMove_CBMC));
    if(probabilitySwapMove_CBMC > 0.0) std::print(stream, formatStatistics("Swap Insertion (CBMC)", statistics_SwapInsertionMove_CBMC));
    if(probabilitySwapMove_CBMC > 0.0) std::print(stream, formatStatistics("Swap Deletion (CBMC)", statistics_SwapDeletionMove_CBMC));
    if(probabilitySwapMove_CFCMC > 0.0) std::print(stream, formatStatistics("Swap (CFCMC)", statistics_SwapMove_CFCMC));
    if(probabilitySwapMove_CFCMC_CBMC > 0.0) std::print(stream, formatStatistics("Swap (CB/CFCMC)", statistics_SwapMove_CFCMC_CBMC)); 
    if(probabilityWidomMove > 0.0) std::print(stream, formatStatistics("Widom (CBMC)", statistics_WidomMove_CBMC));
    if(probabilityWidomMove_CFCMC > 0.0) std::print(stream, formatStatistics("Widom (CFCMC)", statistics_WidomMove_CFCMC));
    if(probabilityWidomMove_CFCMC_CBMC > 0.0) std::print(stream, formatStatistics("Widom (CB/CFCMC)", statistics_WidomMove_CFCMC_CBMC));

    return stream.str();
}


const std::string Component::writeMCMoveCPUTimeStatistics() const
{
    std::ostringstream stream;
    if(cpuTime_TranslationMove.count() > 0.0)
    {
      std::print(stream, "    Translation move:       {:14f} [s]\n", cpuTime_TranslationMove.count());
      std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_TranslationMove_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_TranslationMove_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_TranslationMove.count() - cpuTime_TranslationMove_NonEwald.count() - cpuTime_TranslationMove_Ewald.count());
    }
    if(cpuTime_RandomTranslationMove.count() > 0.0)
    {
      std::print(stream, "    Random translation move:    {:14f} [s]\n", cpuTime_RandomTranslationMove.count());
      std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_RandomTranslationMove_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_RandomTranslationMove_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_RandomTranslationMove.count() - 
                         cpuTime_RandomTranslationMove_NonEwald.count() - cpuTime_RandomTranslationMove_Ewald.count());
    }
    if(cpuTime_RotationMove.count() > 0.0)
    {
      std::print(stream, "    Rotation move:          {:14f} [s]\n", cpuTime_RotationMove.count());
      std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_RotationMove_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_RotationMove_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_RotationMove.count() - cpuTime_RotationMove_NonEwald.count() - cpuTime_RotationMove_Ewald.count());
    }
    if(cpuTime_RandomRotationMove.count() > 0.0)
    {
      std::print(stream, "    Random rotation move:       {:14f} [s]\n", cpuTime_RandomRotationMove.count());
      std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_RandomRotationMove_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_RandomRotationMove_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_RandomRotationMove.count() - cpuTime_RandomRotationMove_NonEwald.count() - cpuTime_RandomRotationMove_Ewald.count());
    }
    if(cpuTime_ReinsertionMove_CBMC.count() > 0.0)
    {
      std::print(stream, "    Reinsertion (CBMC):     {:14f} [s]\n", cpuTime_ReinsertionMove_CBMC.count());
      std::print(stream, "        Grow Non-Ewald:         {:14f} [s]\n", cpuTime_ReinsertionGrowMove_CBMC_NonEwald.count());
      std::print(stream, "        Retrace Non-Ewald:      {:14f} [s]\n", cpuTime_ReinsertionRetraceMove_CBMC_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_ReinsertionMove_CBMC_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_ReinsertionMove_CBMC.count() - cpuTime_ReinsertionGrowMove_CBMC_NonEwald.count() - 
                                                                         cpuTime_ReinsertionRetraceMove_CBMC_NonEwald.count() - cpuTime_ReinsertionMove_CBMC_Ewald.count());
    }
    if(cpuTime_IdentityChangeMove_CBMC.count() > 0.0)
    {
      std::print(stream, "    Identity-change (CBMC): {:14f} [s]\n", cpuTime_IdentityChangeMove_CBMC.count());
      std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_IdentityChangeMove_CBMC_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_IdentityChangeMove_CBMC_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_IdentityChangeMove_CBMC.count() - cpuTime_IdentityChangeMove_CBMC_NonEwald.count() - cpuTime_IdentityChangeMove_CBMC_Ewald.count());
    }
    if(cpuTime_SwapInsertionMove_CBMC.count() > 0.0)
    {
      std::print(stream, "    Swap insert (CBMC):     {:14f} [s]\n", cpuTime_SwapInsertionMove_CBMC.count());
      std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_SwapInsertionMove_CBMC_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_SwapInsertionMove_CBMC_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_SwapInsertionMove_CBMC.count() - cpuTime_SwapInsertionMove_CBMC_NonEwald.count() - cpuTime_SwapInsertionMove_CBMC_Ewald.count());
    }
    if(cpuTime_SwapDeletionMove_CBMC.count() > 0.0)
    {
      std::print(stream, "    Swap delete (CBMC):     {:14f} [s]\n", cpuTime_SwapDeletionMove_CBMC.count());
      std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_SwapDeletionMove_CBMC_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_SwapDeletionMove_CBMC_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_SwapDeletionMove_CBMC.count() - cpuTime_SwapDeletionMove_CBMC_NonEwald.count() - cpuTime_SwapDeletionMove_CBMC_Ewald.count());
    }
    if(cpuTime_SwapMove_CFCMC.count() > 0.0)
    {
      std::print(stream, "    Swap (CFCMC):           {:14f} [s]\n", cpuTime_SwapMove_CFCMC.count());
      std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_SwapMove_CFCMC_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_SwapMove_CFCMC_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_SwapMove_CFCMC.count() - cpuTime_SwapMove_CFCMC_NonEwald.count() - cpuTime_SwapMove_CFCMC_Ewald.count());
    }
    if(cpuTime_SwapMove_CFCMC_CBMC.count() > 0.0)
    {
      std::print(stream, "    Swap (CB/CFCMC):        {:14f} [s]\n", cpuTime_SwapMove_CFCMC_CBMC.count());
      std::print(stream, "        Ins. Non-Ewald:         {:14f} [s]\n", cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald.count());
      std::print(stream, "        Ins. Ewald:             {:14f} [s]\n", cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald.count());
      std::print(stream, "        Ins. Grow Non-Ewald:    {:14f} [s]\n", cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_NonEwald.count());
      std::print(stream, "        Ins. Grow Ewald:        {:14f} [s]\n", cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_Ewald.count());
      std::print(stream, "        Del. Retrace Non-Ewald  {:14f} [s]\n", cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_NonEwald.count());
      std::print(stream, "        Del. Retrace Ewald:     {:14f} [s]\n", cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_Ewald.count());
      std::print(stream, "        Del. Non-Ewald:         {:14f} [s]\n", cpuTime_SwapDeletionMove_CFCMC_CBMC_NonEwald.count());
      std::print(stream, "        Del. Ewald:             {:14f} [s]\n", cpuTime_SwapDeletionMove_CFCMC_CBMC_Ewald.count());
      std::print(stream, "        Lambda Non-Ewald:       {:14f} [s]\n", cpuTime_SwapLambdaMove_CFCMC_CBMC_NonEwald.count());
      std::print(stream, "        Lambda Ewald:           {:14f} [s]\n", cpuTime_SwapLambdaMove_CFCMC_CBMC_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_SwapMove_CFCMC_CBMC.count() -
                                                                         cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald.count() - 
                                                                         cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald.count() -
                                                                         cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_NonEwald.count() -
                                                                         cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_Ewald.count() -
                                                                         cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_NonEwald.count() -
                                                                         cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_Ewald.count() -
                                                                         cpuTime_SwapDeletionMove_CFCMC_CBMC_NonEwald.count() -
                                                                         cpuTime_SwapDeletionMove_CFCMC_CBMC_Ewald.count() -
                                                                         cpuTime_SwapLambdaMove_CFCMC_CBMC_NonEwald.count() -
                                                                         cpuTime_SwapLambdaMove_CFCMC_CBMC_Ewald.count());
    }

    if(cpuTime_WidomMove_CBMC.count() > 0.0)
    {
      std::print(stream, "    Widom:                  {:14f} [s]\n", cpuTime_WidomMove_CBMC.count());
      std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_WidomMove_CBMC_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_WidomMove_CBMC_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_WidomMove_CBMC.count() - cpuTime_WidomMove_CBMC_NonEwald.count() - cpuTime_WidomMove_CBMC_Ewald.count());
    }
    if(cpuTime_WidomMove_CFCMC.count() > 0.0)
    {
      std::print(stream, "    Widom (CFCMC):          {:14f} [s]\n", cpuTime_WidomMove_CFCMC.count());
      std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_WidomMove_CFCMC_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_WidomMove_CFCMC_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_WidomMove_CFCMC.count() - cpuTime_WidomMove_CFCMC_NonEwald.count() - cpuTime_WidomMove_CFCMC_Ewald.count());
    }
    if(cpuTime_WidomMove_CFCMC_CBMC.count() > 0.0)
    {
      std::print(stream, "    Widom (CB/CFCMC):       {:14f} [s]\n", cpuTime_WidomMove_CFCMC_CBMC.count());
      std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_WidomMove_CFCMC_CBMC_NonEwald.count());
      std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_WidomMove_CFCMC_CBMC_Ewald.count());
      std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_WidomMove_CFCMC_CBMC.count() - cpuTime_WidomMove_CFCMC_CBMC_NonEwald.count() - cpuTime_WidomMove_CFCMC_CBMC_Ewald.count());
    }
    std::print(stream, "\n");
    
    return stream.str(); 
}


std::vector<double3> Component::randomlyRotatedPositionsAroundStartingBead() const
{
    double3x3 randomRotationMatrix = double3x3::randomRotationMatrix();
    std::vector<double3> randomPositions{};
    std::transform(std::begin(atoms), std::end(atoms),
            std::back_inserter(randomPositions), [&](const Atom& atom) {return randomRotationMatrix * (atom.position - atoms[startingBead].position); });
    return randomPositions;
}

std::vector<Atom> Component::newAtoms(double scaling, size_t moleculeId) const
{
    std::vector<Atom> new_atoms(atoms);

#if defined(_WIN32)
    //assert(atomPositions.size() == numberOfAtoms);
    //assert(atomCharges.size() == numberOfAtoms);
#else
    //_LIBCPP_ASSERT(atomPositions.size() == numberOfAtoms, "wrong number of atoms");
    //_LIBCPP_ASSERT(atomCharges.size() == numberOfAtoms, "wrong number of atoms");
#endif

    for (size_t i = 0; i < atoms.size(); ++i)
    {
        new_atoms[i] = Atom(atoms[i].position - atoms[startingBead].position, 
                   atoms[i].charge, scaling, static_cast<short>(atoms[i].type), 
                   static_cast<short>(componentId), static_cast<int>(moleculeId));
    }

    return new_atoms;
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

void Component::clearMoveStatistics()
{
  statistics_TranslationMove.clear();
  statistics_RandomTranslationMove.clear();
  statistics_RotationMove.clear();
  statistics_RandomRotationMove.clear();
  statistics_ReinsertionMove_CBMC.clear();
  statistics_IdentityChangeMove_CBMC.clear();
  statistics_SwapInsertionMove_CBMC.clear();
  statistics_SwapDeletionMove_CBMC.clear();
  statistics_SwapMove_CFCMC.clear();
  statistics_SwapMove_CFCMC_CBMC.clear();
  statistics_WidomMove_CBMC.clear();
  statistics_WidomMove_CFCMC.clear();
  statistics_WidomMove_CFCMC_CBMC.clear();
}

void Component::optimizeMCMoves()
{
    statistics_TranslationMove.optimizeAcceptance();
    statistics_RotationMove.optimizeAcceptance();
}

void Component::clearTimingStatistics()
{
    cpuTime_TranslationMove = std::chrono::duration<double>(0.0);
    cpuTime_RandomTranslationMove = std::chrono::duration<double>(0.0);
    cpuTime_RotationMove = std::chrono::duration<double>(0.0);
    cpuTime_RandomRotationMove = std::chrono::duration<double>(0.0);
    cpuTime_ReinsertionMove_CBMC = std::chrono::duration<double>(0.0);
    cpuTime_IdentityChangeMove_CBMC = std::chrono::duration<double>(0.0);
    cpuTime_SwapInsertionMove_CBMC = std::chrono::duration<double>(0.0);
    cpuTime_SwapDeletionMove_CBMC = std::chrono::duration<double>(0.0);
    cpuTime_SwapMove_CFCMC = std::chrono::duration<double>(0.0);
    cpuTime_SwapMove_CFCMC_CBMC = std::chrono::duration<double>(0.0);
    cpuTime_WidomMove_CBMC = std::chrono::duration<double>(0.0);
    cpuTime_WidomMove_CFCMC = std::chrono::duration<double>(0.0);
    cpuTime_WidomMove_CFCMC_CBMC = std::chrono::duration<double>(0.0);

    cpuTime_TranslationMove_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_RandomTranslationMove_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_RotationMove_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_RandomRotationMove_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_ReinsertionGrowMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_ReinsertionRetraceMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_IdentityChangeMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_SwapInsertionMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_SwapDeletionMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_SwapMove_CFCMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_SwapDeletionMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_SwapLambdaMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_WidomMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_WidomMove_CFCMC_NonEwald = std::chrono::duration<double>(0.0);
    cpuTime_WidomMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);

    cpuTime_TranslationMove_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_RandomTranslationMove_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_RotationMove_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_RandomRotationMove_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_ReinsertionMove_CBMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_IdentityChangeMove_CBMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_SwapInsertionMove_CBMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_SwapDeletionMove_CBMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_SwapMove_CFCMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_SwapDeletionMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_SwapLambdaMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_WidomMove_CBMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_WidomMove_CFCMC_Ewald = std::chrono::duration<double>(0.0);
    cpuTime_WidomMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
}

void Component::printBreakthroughStatus(std::ostream &stream) const
{
  std::print(stream, "Component {} [{}]\n", componentId, name);
  if(isCarrierGas)
  {
    std::print(stream, "    carrier-gas\n");

    std::print(stream, isotherm.print());
  }
  std::print(stream, "    mol-fraction in the gas:   {} [-]\n", molFraction);
  if(!isCarrierGas)
  {
    std::print(stream, "    mass-transfer coefficient: {} [1/s]\n", massTransferCoefficient);
    std::print(stream, "    diffusion coefficient:     {} [m^2/s]\n", axialDispersionCoefficient);

    std::print(stream, isotherm.print());
  }
}
