module;

#if !defined(_WIN32)
#include <assert.h>
#endif

module component;

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
import <cstddef>;
#if defined(_WIN32)
import <cassert>;
#endif

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
import property_lambda_probability_histogram;
import print;
import property_widom;
import multi_site_isotherm;
import simulationbox;
import cif_reader;
import move_statistics;
import bond_potential;
import mc_moves_probabilities_particles;



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
    lambdaGC(numberOfBlocks, 41),
    lambdaGibbs(numberOfBlocks, 41),
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
    lambdaGC(numberOfBlocks, 41),
    lambdaGibbs(numberOfBlocks, 41),
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
    atoms[i].componentId = static_cast<std::byte>(componentId);
    atoms[i].moleculeId = 0;
  }
}

Component::Component(Component::Type type, size_t currentComponent, const std::string &componentName,
                     std::optional<const std::string> fileName, size_t numberOfBlocks) noexcept(false) :
                     type(type), 
                     componentId(currentComponent), 
                     name(componentName),
                     filenameData(fileName),
                     lambdaGC(numberOfBlocks, 41),
                     lambdaGibbs(numberOfBlocks, 41),
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

      definedAtoms[i] = Atom(pos, charge, scaling, static_cast<short>(pseudoAtomType), static_cast<std::byte>(componentId), 0);
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
    atoms[i].componentId = static_cast<std::byte>(componentId);
    atoms[i].moleculeId = 0;
  }  
}

std::string Component::printStatus(const ForceField& forceField) const
{
  std::ostringstream stream;

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
  
  std::print(stream, "    Translation-move probability:             {} [-]\n", mc_moves_probabilities.probabilityTranslationMove);
  std::print(stream, "    Random translation-move probability:      {} [-]\n", mc_moves_probabilities.probabilityRandomTranslationMove);
  std::print(stream, "    Rotation-move probability:                {} [-]\n", mc_moves_probabilities.probabilityRotationMove);
  std::print(stream, "    Random rotation-move probability:         {} [-]\n", mc_moves_probabilities.probabilityRandomRotationMove);
  std::print(stream, "    Volume-move probability:                  {} [-]\n", mc_moves_probabilities.probabilityVolumeMove);
  std::print(stream, "    Reinsertion (CBMC) probability:           {} [-]\n", mc_moves_probabilities.probabilityReinsertionMove_CBMC);
  std::print(stream, "    Identity-change (CBMC) probability:       {} [-]\n", mc_moves_probabilities.probabilityIdentityChangeMove_CBMC);
  std::print(stream, "    Swap-move (CBMC) probability:             {} [-]\n", mc_moves_probabilities.probabilitySwapMove_CBMC);
  std::print(stream, "    Swap-move (CFCMC) probability:            {} [-]\n", mc_moves_probabilities.probabilitySwapMove_CFCMC);
  std::print(stream, "    Swap-move (CFCMC/CBMC) probability:       {} [-]\n", mc_moves_probabilities.probabilitySwapMove_CFCMC_CBMC);
  std::print(stream, "    Gibbs Volume-move probability:            {} [-]\n", mc_moves_probabilities.probabilityGibbsVolumeMove);
  std::print(stream, "    Gibbs Swap-move (CBMC) probability:       {} [-]\n", mc_moves_probabilities.probabilityGibbsSwapMove_CBMC);
  std::print(stream, "    Gibbs Swap-move (CFCMC) probability:      {} [-]\n", mc_moves_probabilities.probabilityGibbsSwapMove_CFCMC);
  std::print(stream, "    Gibbs Swap-move (CFCMC/CBMC) probability: {} [-]\n", mc_moves_probabilities.probabilityGibbsSwapMove_CFCMC_CBMC);
  std::print(stream, "    Widom probability:                        {} [-]\n", mc_moves_probabilities.probabilityWidomMove);
  std::print(stream, "    Widom (CFCMC) probability:                {} [-]\n", mc_moves_probabilities.probabilityWidomMove_CFCMC);
  std::print(stream, "    Widom (CFCMC/CBMC) probability:           {} [-]\n", mc_moves_probabilities.probabilityWidomMove_CFCMC_CBMC);
  std::print(stream, "\n");

  std::print(stream, "    number of bonds: {}\n", bonds.size());
  for (size_t i = 0; i < bonds.size(); ++i)
  {
      std::print(stream, "        {}", bonds[i].print());
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
                   static_cast<std::byte>(componentId), static_cast<int>(moleculeId));
    }

    return new_atoms;
}

std::vector<Atom> Component::copyAtoms(std::span<Atom> molecule, double scaling, size_t moleculeId) const
{
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());
  for (size_t i = 0; i != atoms.size(); ++i)
  {
    copied_atoms[i].setScaling(scaling);
    copied_atoms[i].position = molecule[i].position - molecule[startingBead].position;
    copied_atoms[i].moleculeId = static_cast<int>(moleculeId);
  }
  return copied_atoms;
}

std::vector<Atom> Component::copyAtomsRandomlyRotatedAt(double3 position, std::span<Atom> molecule, double scaling, size_t moleculeId) const
{
  double3x3 randomRotationMatrix = double3x3::randomRotationMatrix();
  std::vector<Atom> copied_atoms(molecule.begin(), molecule.end());
  for (size_t i = 0; i != atoms.size(); ++i)
  {
    copied_atoms[i].setScaling(scaling);
    copied_atoms[i].position = position + randomRotationMatrix * (molecule[i].position - molecule[startingBead].position);
    copied_atoms[i].moleculeId = static_cast<int>(moleculeId);
  }
  return copied_atoms;
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
