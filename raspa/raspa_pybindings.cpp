#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#endif

#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <exception>

#pragma clang diagnostic pop

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#endif

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import double3;
import double3x3;
import atom;
import simulationbox;
import pseudo_atom;
import vdwparameters;
import forcefield;
import running_energy;
import mc_moves_probabilities;
import mc_moves_move_types;
import move_statistics;
import framework;
import component;
import system;
import randomnumbers;
import monte_carlo;
import input_reader;
import property_loading;
import loadings;
import property_lambda_probability_histogram;
import connectivity_table;
import intra_molecular_potentials;

PYBIND11_MODULE(raspalib, m)
{
  pybind11::class_<int3>(m, "int3")
      .def(pybind11::init<int32_t, int32_t, int32_t>(), pybind11::arg("x") = 0, pybind11::arg("y") = 0,
           pybind11::arg("z") = 0)
      .def_readwrite("x", &int3::x)
      .def_readwrite("y", &int3::y)
      .def_readwrite("z", &int3::z);

  pybind11::class_<double3>(m, "double3")
      .def(pybind11::init<double, double, double>(), pybind11::arg("x") = 0.0, pybind11::arg("y") = 0.0,
           pybind11::arg("z") = 0.0)
      .def_readwrite("x", &double3::x)
      .def_readwrite("y", &double3::y)
      .def_readwrite("z", &double3::z);

  pybind11::class_<RandomNumber>(m, "RandomNumber").def(pybind11::init<std::size_t>(), pybind11::arg("seed") = 12);

  pybind11::class_<RunningEnergy>(m, "RunningEnergy")
      .def(pybind11::init<>())
      .def_readwrite("moleculeMoleculeVDW", &RunningEnergy::moleculeMoleculeVDW)
      .def_readwrite("frameworkMoleculeVDW", &RunningEnergy::frameworkMoleculeVDW);

  pybind11::class_<Atom>(m, "Atom")
      .def(pybind11::init<>())
      .def(pybind11::init<double3, double, double, uint32_t, uint16_t, uint8_t, bool, bool>(),
           pybind11::arg("position"), pybind11::arg("charge") = 0.0, pybind11::arg("lambda") = 0.0,
           pybind11::arg("moleculeId") = 0, pybind11::arg("type") = 0, pybind11::arg("componentId") = 0,
           pybind11::arg("groupId") = 0, pybind11::arg("isFractional") = 0)
      .def_readwrite("position", &Atom::position)
      .def("__repr__", &Atom::repr);

  pybind11::class_<SimulationBox> simulationBox(m, "SimulationBox");
  pybind11::enum_<SimulationBox::Type>(simulationBox, "Type")
      .value("Rectangular", SimulationBox::Type::Rectangular)
      .value("Triclinic", SimulationBox::Type::Triclinic)
      .export_values();
  simulationBox
      .def(pybind11::init<double, double, double>(), pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"))
      .def(pybind11::init<double, double, double, double, double, double>(), pybind11::arg("a"), pybind11::arg("b"),
           pybind11::arg("c"), pybind11::arg("alpha"), pybind11::arg("beta"), pybind11::arg("gamma"))
      .def(pybind11::init<double, double, double, double, double, double, SimulationBox::Type>(), pybind11::arg("a"),
           pybind11::arg("b"), pybind11::arg("c"), pybind11::arg("alpha"), pybind11::arg("beta"),
           pybind11::arg("gamma"), pybind11::arg("type"))
      .def(pybind11::init<double3x3, SimulationBox::Type>())
      .def_readwrite("type", &SimulationBox::type)
      .def_readonly("lengthA", &SimulationBox::lengthA)
      .def_readonly("lengthB", &SimulationBox::lengthB)
      .def_readonly("lengthC", &SimulationBox::lengthC);

  pybind11::class_<VDWParameters>(m, "VDWParameters")
      .def(pybind11::init<double, double>(), pybind11::arg("epsilon"), pybind11::arg("sigma"));

  pybind11::class_<PseudoAtom>(m, "PseudoAtom")
      .def(pybind11::init<std::string, bool, double, double, double, std::size_t, bool, std::string>(),
           pybind11::arg("name"), pybind11::arg("frameworkType"), pybind11::arg("mass"), pybind11::arg("charge"),
           pybind11::arg("polarizability") = 0.0, pybind11::arg("atomicNumber"), pybind11::arg("printToPDB") = true,
           pybind11::arg("source") = "");

  pybind11::class_<ForceField> forceField(m, "ForceField");
  forceField
      .def(pybind11::init<std::vector<PseudoAtom>, std::vector<VDWParameters>, ForceField::MixingRule, double, double,
                          double, bool, bool, bool>(),
           pybind11::arg("pseudoAtoms"), pybind11::arg("parameters"), pybind11::arg("mixingRule"),
           pybind11::arg("cutOffFrameworkVDW") = 12.0, pybind11::arg("cutOffMoleculeVDW") = 12.0,
           pybind11::arg("cutOffCoulomb") = 12.0, pybind11::arg("shifted") = true,
           pybind11::arg("tailCorrections") = false, pybind11::arg("useCharge") = true)
      .def(pybind11::init<std::string>(), pybind11::arg("fileName") = "force_field.json")
      .def("__repr__", &ForceField::repr)
      .def_readonly("pseudoAtoms", &ForceField::pseudoAtoms)
      .def_readonly("vdwParameters", &ForceField::data)
      .def_readwrite("useCharge", &ForceField::useCharge);

  pybind11::enum_<ForceField::MixingRule>(forceField, "MixingRule")
      .value("Lorentz_Berthelot", ForceField::MixingRule::Lorentz_Berthelot)
      .export_values();

  pybind11::class_<Framework> framework(m, "Framework");

  pybind11::enum_<Framework::UseChargesFrom>(framework, "UseChargesFrom")
      .value("PseudoAtoms", Framework::UseChargesFrom::PseudoAtoms)
      .value("CIF_File", Framework::UseChargesFrom::CIF_File)
      .value("ChargeEquilibration", Framework::UseChargesFrom::ChargeEquilibration);

  framework
      .def(pybind11::init<std::size_t, const ForceField &, std::string, SimulationBox, std::size_t, std::vector<Atom>,
                          int3>(),
           pybind11::arg("frameworkId"), pybind11::arg("forceField"), pybind11::arg("componentName"),
           pybind11::arg("simulationBox"), pybind11::arg("spaceGroupHallNumber"), pybind11::arg("definedAtoms"),
           pybind11::arg("numberOfUnitCells"))
      .def(pybind11::init<std::size_t, const ForceField &, const std::string &, std::optional<const std::string>, int3,
                          Framework::UseChargesFrom>(),
           pybind11::arg("frameworkId"), pybind11::arg("forceField"), pybind11::arg("componentName"),
           pybind11::arg("fileName"), pybind11::arg("numberOfUnitCells"),
           pybind11::arg("useChargesFrom") = Framework::UseChargesFrom::PseudoAtoms)
      .def_readonly("name", &Framework::name)
      .def("__repr__", &Framework::repr);

  pybind11::class_<MCMoveProbabilities>(m, "MCMoveProbabilities")
      .def(pybind11::init<double, double, double, double, double, double, double, double, double, double, double,
                          double, double, double, double, double, double, double, double>(),
           pybind11::arg("translationProbability") = 0.0, pybind11::arg("randomTranslationProbability") = 0.0,
           pybind11::arg("rotationProbability") = 0.0, pybind11::arg("randomRotationProbability") = 0.0,
           pybind11::arg("volumeChangeProbability") = 0.0, pybind11::arg("reinsertionCBMCProbability") = 0.0,
           pybind11::arg("identityChangeProbability") = 0.0, pybind11::arg("swapProbability") = 0.0,
           pybind11::arg("swapCBMCProbability") = 0.0, pybind11::arg("swapCFCMCProbability") = 0.0,
           pybind11::arg("swapCBCFCMCProbability") = 0.0, pybind11::arg("gibbsVolumeChangeProbability") = 0.0,
           pybind11::arg("gibbsSwapCBMCProbability") = 0.0, pybind11::arg("gibbsSwapCFCMCProbability") = 0.0,
           pybind11::arg("widomProbability") = 0.0, pybind11::arg("widomCFCMCProbability") = 0.0,
           pybind11::arg("widomCBCFCMCProbability") = 0.0, pybind11::arg("parallelTemperingProbability") = 0.0,
           pybind11::arg("hybridMCProbability") = 0.0)
      .def("join", &MCMoveProbabilities::join);

  pybind11::class_<PropertyLambdaProbabilityHistogram>(m, "PropertyLambdaProbabilityHistogram")
      .def(pybind11::init<>())
      .def_readonly("biasFactor", &PropertyLambdaProbabilityHistogram::biasFactor)
      .def_readonly("histogram", &PropertyLambdaProbabilityHistogram::histogram)
      .def("normalizedAverageProbabilityHistogram",
           &PropertyLambdaProbabilityHistogram::normalizedAverageProbabilityHistogram);

  pybind11::class_<ConnectivityTable>(m, "ConnectivityTable").def(pybind11::init<>());

  pybind11::class_<Potentials::IntraMolecularPotentials>(m, "IntraMolecularPotentials").def(pybind11::init<>());

  // define before component init to prevent failing default argument
  pybind11::class_<Component> component(m, "Component");

  pybind11::enum_<Component::Type>(component, "Type")
      .value("Adsorbate", Component::Type::Adsorbate)
      .value("Cation", Component::Type::Cation)
      .export_values();

  component
      .def(pybind11::init<std::size_t, const ForceField &, std::string, double, double, double, std::vector<Atom>,
                          const ConnectivityTable &, const Potentials::IntraMolecularPotentials &, std::size_t,
                          std::size_t, const MCMoveProbabilities &, std::optional<double>, bool>(),
           pybind11::arg("componentId"), pybind11::arg("forceField"), pybind11::arg("componentName"),
           pybind11::arg("criticalTemperature"), pybind11::arg("criticalPressure"), pybind11::arg("acentricFactor"),
           pybind11::arg("definedAtoms"), pybind11::arg("connectivitytable") = ConnectivityTable(),
           pybind11::arg("intraMolecularPotentials") = Potentials::IntraMolecularPotentials(),
           pybind11::arg("numberOfBlocks") = 5, pybind11::arg("numberOfLambdaBins") = 41,
           pybind11::arg("particleProbabilities"), pybind11::arg("fugacityCoefficient") = std::nullopt,
           pybind11::arg("thermodynamicIntegration") = false)
      .def(pybind11::init<Component::Type, std::size_t, const ForceField &, std::string &, std::string, std::size_t,
                          std::size_t, const MCMoveProbabilities &, std::optional<double>, bool>(),
           pybind11::arg("type") = Component::Type::Adsorbate, pybind11::arg("componentId"),
           pybind11::arg("forceField"), pybind11::arg("componentName"), pybind11::arg("fileName"),
           pybind11::arg("numberOfBlocks") = 5, pybind11::arg("numberOfLambdaBins") = 41,
           pybind11::arg("particleProbabilities") = MCMoveProbabilities(),
           pybind11::arg("fugacityCoefficient") = std::nullopt, pybind11::arg("thermodynamicIntegration") = false)
      .def_readonly("name", &Component::name)
      .def_readonly("lambdaGC", &Component::lambdaGC)
      .def_readonly("mc_moves_statistics", &Component::mc_moves_statistics)
      .def("__repr__", &Component::repr);

  pybind11::class_<Loadings>(m, "Loadings")
      .def(pybind11::init<std::size_t>())
      .def_readonly("numberOfMolecules", &Loadings::numberOfMolecules)
      .def("printStatus",
           static_cast<std::string (Loadings::*)(const Component &, std::optional<double>, std::optional<int3>) const>(
               &Loadings::printStatus))
      .def("printStatus", static_cast<std::string (Loadings::*)(const Component &, const Loadings &, const Loadings &,
                                                                std::optional<double>, std::optional<int3>) const>(
                              &Loadings::printStatus));

  pybind11::class_<PropertyLoading>(m, "PropertyLoading")
      .def(pybind11::init<std::size_t, std::size_t>())
      .def("averageLoading", &PropertyLoading::averageLoading)
      .def("averageLoadingNumberOfMolecules", &PropertyLoading::averageLoadingNumberOfMolecules)
      .def("writeAveragesStatistics", &PropertyLoading::writeAveragesStatistics)
      .def("__repr__", &PropertyLoading::repr);

  pybind11::class_<System>(m, "System")
      .def(pybind11::init<std::size_t, ForceField, std::optional<SimulationBox>, double, std::optional<double>, double,
                          std::optional<Framework>, std::vector<Component>, std::vector<std::vector<double3>>,
                          std::vector<std::size_t>, std::size_t, MCMoveProbabilities, std::optional<std::size_t>>(),
           pybind11::arg("systemId"), pybind11::arg("forceField"), pybind11::arg("simulationBox") = std::nullopt,
           pybind11::arg("externalTemperature"), pybind11::arg("externalPressure") = std::nullopt,
           pybind11::arg("heliumVoidFraction") = 0.0, pybind11::arg("frameworkComponents") = std::vector<Framework>(),
           pybind11::arg("components"), pybind11::arg("initialPositions") = std::vector<std::vector<double3>>(),
           pybind11::arg("initialNumberOfMolecules") = std::vector<std::size_t>(), pybind11::arg("numberOfBlocks") = 5,
           pybind11::arg("systemProbabilities") = MCMoveProbabilities(),
           pybind11::arg("sampleMoviesEvery") = std::nullopt)
      .def("computeTotalEnergies", &System::computeTotalEnergies)
      .def("frameworkMass", &System::frameworkMass)
      .def_readonly("inputPressure", &System::input_pressure)
      .def_readonly("components", &System::components)
      .def_readonly("loadings", &System::loadings)
      .def_readonly("averageLoadings", &System::averageLoadings)
      .def_readwrite("atomData", &System::atomData)
      .def("writeMCMoveStatistics", &System::writeMCMoveStatistics)
      .def("__repr__", &System::repr);

  pybind11::class_<InputReader> inputReader(m, "InputReader");
  inputReader.def(pybind11::init<const std::string>(), pybind11::arg("fileName") = "simulation.json")
      .def_readonly("numberOfBlocks", &InputReader::numberOfBlocks)
      .def_readonly("numberOfCycles", &InputReader::numberOfCycles)
      .def_readonly("numberOfInitializationCycles", &InputReader::numberOfInitializationCycles)
      .def_readonly("numberOfEquilibrationCycles", &InputReader::numberOfEquilibrationCycles)
      .def_readonly("printEvery", &InputReader::printEvery)
      .def_readonly("writeBinaryRestartEvery", &InputReader::writeBinaryRestartEvery)
      .def_readonly("rescaleWangLandauEvery", &InputReader::rescaleWangLandauEvery)
      .def_readonly("optimizeMCMovesEvery", &InputReader::optimizeMCMovesEvery)
      .def_readonly("writeEvery", &InputReader::writeEvery)
      .def_readonly("forceField", &InputReader::forceField)
      .def_readonly("systems", &InputReader::systems);

  pybind11::enum_<InputReader::SimulationType>(inputReader, "SimulationType")
      .value("MonteCarlo", InputReader::SimulationType::MonteCarlo)
      .value("MonteCarloTransitionMatrix", InputReader::SimulationType::MonteCarloTransitionMatrix)
      .value("MolecularDynamics", InputReader::SimulationType::MolecularDynamics)
      .value("Minimization", InputReader::SimulationType::Minimization)
      .value("Test", InputReader::SimulationType::Test)
      .value("Breakthrough", InputReader::SimulationType::Breakthrough)
      .value("MixturePrediction", InputReader::SimulationType::MixturePrediction)
      .value("Fitting", InputReader::SimulationType::Fitting)
      .value("ParallelTempering", InputReader::SimulationType::ParallelTempering)
      .export_values();

  pybind11::class_<MonteCarlo> mc(m, "MonteCarlo");
  mc.def(pybind11::init<std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t,
                        std::vector<System> &, RandomNumber &, std::size_t, bool>(),
         pybind11::arg("numberOfCycles"), pybind11::arg("numberOfInitializationCycles"),
         pybind11::arg("numberOfEquilibrationCycles") = 0, pybind11::arg("printEvery") = 5000,
         pybind11::arg("writeBinaryRestartEvery") = 5000, pybind11::arg("rescaleWangLandauEvery") = 1000,
         pybind11::arg("optimizeMCMovesEvery") = 100, pybind11::arg("systems"),
         pybind11::arg("randomSeed") = RandomNumber(12), pybind11::arg("numberOfBlocks") = 5,
         pybind11::arg("outputToFiles") = false)
      .def(pybind11::init<InputReader &>(), pybind11::arg("inputReader"))
      .def("run", &MonteCarlo::run)
      .def("initialize", &MonteCarlo::initialize)
      .def("equilibrate", &MonteCarlo::equilibrate)
      .def("production", &MonteCarlo::production)
      .def("cycle", &MonteCarlo::performCycle)
      .def_readonly("systems", &MonteCarlo::systems)
      .def_readwrite("simulationStage", &MonteCarlo::simulationStage);

  pybind11::enum_<MonteCarlo::SimulationStage>(mc, "SimulationStage")
      .value("Uninitialized", MonteCarlo::SimulationStage::Uninitialized)
      .value("Initialization", MonteCarlo::SimulationStage::Initialization)
      .value("Equilibration", MonteCarlo::SimulationStage::Equilibration)
      .value("Production", MonteCarlo::SimulationStage::Production)
      .export_values();
}
