#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <print>
#include <string>
#include <tuple>
#include <vector>
#include <functional>
#endif

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/native_enum.h>

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
import double4;
import double3x3;
import units;
import atom;
import cif_reader;
import simulationbox;
import pseudo_atom;
import vdwparameters;
import forcefield;
import running_energy;
import mc_moves_probabilities;
import mc_moves_move_types;
import move_statistics;
import mc_moves_statistics;
import cif_reader;
import framework;
import component;
import system;
import randomnumbers;
import monte_carlo;
import molecular_dynamics;
import input_reader;
import property_loading;
import loadings;
import energy_factor;
import energy_status;
import sample_movies;
import property_energy;
import average_energy_type;
import property_energy_histogram;
import property_lambda_probability_histogram;
import property_density_grid;
import property_conventional_rdf;
import property_rdf;
import property_number_of_molecules_evolution;
import property_volume_evolution;
import property_widom;
import connectivity_table;
import intra_molecular_potentials;


template<typename T>
std::pair<T, T> operator*(const double& a, const std::pair<T, T>& b)
{
  return std::make_pair(a * b.first, a * b.second);
}

PYBIND11_MODULE(raspalib, m)
{
  pybind11::class_<int3>(m, "int3")
      .def(pybind11::init<int32_t, int32_t, int32_t>(), pybind11::arg("x") = 0, pybind11::arg("y") = 0,
           pybind11::arg("z") = 0)
      .def_readwrite("x", &int3::x)
      .def_readwrite("y", &int3::y)
      .def_readwrite("z", &int3::z)
      .def("__repr__", &int3::repr);

  pybind11::class_<double3>(m, "double3")
      .def(pybind11::init<double, double, double>(), pybind11::arg("x") = 0.0, pybind11::arg("y") = 0.0,
           pybind11::arg("z") = 0.0)
      .def_readwrite("x", &double3::x)
      .def_readwrite("y", &double3::y)
      .def_readwrite("z", &double3::z)
      .def("__repr__", &double3::repr);

  pybind11::class_<double4>(m, "double4")
      .def(pybind11::init<double, double, double, double>(), pybind11::arg("x") = 0.0, pybind11::arg("y") = 0.0,
           pybind11::arg("z") = 0.0, pybind11::arg("w") = 0.0)
      .def_readwrite("x", &double4::x)
      .def_readwrite("y", &double4::y)
      .def_readwrite("z", &double4::z)
      .def_readwrite("w", &double4::w);

  pybind11::class_<RandomNumber>(m, "RandomNumber").def(pybind11::init<std::size_t>(), pybind11::arg("seed") = 12);

  pybind11::class_<RunningEnergy>(m, "RunningEnergy")
      .def(pybind11::init<>())
      .def_readwrite("moleculeMoleculeVDW", &RunningEnergy::moleculeMoleculeVDW)
      .def_readwrite("frameworkMoleculeVDW", &RunningEnergy::frameworkMoleculeVDW)
      .def("__repr__", &RunningEnergy::repr);


  pybind11::class_<Atom>(m, "Atom")
      .def(pybind11::init<>())
      .def(pybind11::init<double3, double, double, std::uint32_t, std::uint16_t, std::uint8_t, std::uint8_t, std::uint8_t>(),
           pybind11::arg("position"), pybind11::arg("charge") = 0.0, pybind11::arg("scaling") = 1.0,
           pybind11::arg("moleculeId") = 0, pybind11::arg("type") = 0, pybind11::arg("componentId") = 0,
           pybind11::arg("groupId") = false, pybind11::arg("isFractional") = false)
      .def_readwrite("position", &Atom::position)
      .def("__repr__", &Atom::repr);

  pybind11::class_<SimulationBox> simulationBox(m, "SimulationBox");
  pybind11::native_enum<SimulationBox::Type>(simulationBox, "SimulationBoxType", "enum.IntEnum")
      .value("Rectangular", SimulationBox::Type::Rectangular)
      .value("Triclinic", SimulationBox::Type::Triclinic)
      .finalize();
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
           pybind11::arg("polarizability") = 0.0, pybind11::arg("atomicNumber") = 1, pybind11::arg("printToPDB") = true,
           pybind11::arg("source") = "");

  pybind11::class_<ForceField> forceField(m, "ForceField");
  pybind11::native_enum<ForceField::MixingRule>(forceField, "MixingRule", "enum.IntEnum")
      .value("Lorentz_Berthelot", ForceField::MixingRule::Lorentz_Berthelot)
      .value("Jorgensen", ForceField::MixingRule::Jorgensen)
      .finalize();
  forceField.def(pybind11::init<std::vector<PseudoAtom>, std::vector<VDWParameters>, ForceField::MixingRule, 
                 double, double, double, bool, bool, bool>(),
           pybind11::arg("pseudoAtoms"), pybind11::arg("parameters"), 
           pybind11::arg("mixingRule") = ForceField::MixingRule::Lorentz_Berthelot,
           pybind11::arg("cutOffFrameworkVDW") = 12.0, pybind11::arg("cutOffMoleculeVDW") = 12.0,
           pybind11::arg("cutOffCoulomb") = 12.0, pybind11::arg("shifted") = true,
           pybind11::arg("tailCorrections") = false, pybind11::arg("useCharge") = true)
      .def("findPseudoAtom", static_cast<std::optional<std::size_t> (ForceField::*)(const std::string &) const>(&ForceField::findPseudoAtom))
      .def("__repr__", &ForceField::repr)
      .def_readonly("pseudoAtoms", &ForceField::pseudoAtoms)
      .def_readonly("vdwParameters", &ForceField::data)
      .def_readwrite("useCharge", &ForceField::useCharge);

  pybind11::class_<CIFReader> cifReader(m, "CIFReader");
  cifReader
      .def("expandDefinedAtomsToUnitCell", &CIFReader::expandDefinedAtomsToUnitCell);

  pybind11::class_<Framework> framework(m, "Framework");

  pybind11::native_enum<CIFReader::UseChargesFrom>(framework, "UseChargesFrom", "enum.IntEnum")
      .value("PseudoAtoms", CIFReader::UseChargesFrom::PseudoAtoms)
      .value("CIF_File", CIFReader::UseChargesFrom::CIF_File)
      .value("ChargeEquilibration", CIFReader::UseChargesFrom::ChargeEquilibration)
      .finalize();

  framework
      .def(pybind11::init<std::size_t, const ForceField &, std::string, SimulationBox, std::size_t, 
                          const std::vector<Atom> &, int3>(),
           pybind11::arg("frameworkId"), pybind11::arg("forceField"), pybind11::arg("componentName"),
           pybind11::arg("simulationBox"), pybind11::arg("spaceGroupHallNumber"), 
           pybind11::arg("definedAtoms"), 
           pybind11::arg("numberOfUnitCells"))
      .def_readonly("name", &Framework::name)
      .def("print", &Framework::printStatus)
      .def("__repr__", &Framework::repr);

  pybind11::class_<MCMoveProbabilities>(m, "MCMoveProbabilities")
      .def(pybind11::init<double, double, double, double, double, double, double, double, double, double, double,
                          double, double, double, double, double, double, double, double, double>(),
           pybind11::arg("translationProbability") = 0.0, pybind11::arg("randomTranslationProbability") = 0.0,
           pybind11::arg("rotationProbability") = 0.0, pybind11::arg("randomRotationProbability") = 0.0,
           pybind11::arg("volumeChangeProbability") = 0.0, pybind11::arg("reinsertionCBMCProbability") = 0.0,
           pybind11::arg("partialReinsertionCBMCProbability") = 0.0,
           pybind11::arg("identityChangeProbability") = 0.0, pybind11::arg("swapProbability") = 0.0,
           pybind11::arg("swapCBMCProbability") = 0.0, pybind11::arg("swapCFCMCProbability") = 0.0,
           pybind11::arg("swapCBCFCMCProbability") = 0.0, pybind11::arg("gibbsVolumeChangeProbability") = 0.0,
           pybind11::arg("gibbsSwapCBMCProbability") = 0.0, pybind11::arg("gibbsSwapCFCMCProbability") = 0.0,
           pybind11::arg("widomProbability") = 0.0, pybind11::arg("widomCFCMCProbability") = 0.0,
           pybind11::arg("widomCBCFCMCProbability") = 0.0, pybind11::arg("parallelTemperingProbability") = 0.0,
           pybind11::arg("hybridMCProbability") = 0.0)
      .def("join", &MCMoveProbabilities::join)
      .def("__repr__", &MCMoveProbabilities::repr);


  pybind11::class_<PropertyLambdaProbabilityHistogram>(m, "PropertyLambdaProbabilityHistogram")
      .def(pybind11::init<>())
      .def_readwrite("biasFactor", &PropertyLambdaProbabilityHistogram::biasFactor)
      .def_readwrite("histogram", &PropertyLambdaProbabilityHistogram::histogram)
      .def("result", &PropertyLambdaProbabilityHistogram::result);


  pybind11::class_<ConnectivityTable>(m, "ConnectivityTable")
    .def(pybind11::init<>());

  pybind11::class_<Potentials::IntraMolecularPotentials>(m, "IntraMolecularPotentials")
    .def(pybind11::init<>());

  pybind11::class_<MoveStatistics<double>>(m, "MoveStatisticsDouble")
    .def_readonly("counts", &MoveStatistics<double>::counts)
    .def_readonly("constructed", &MoveStatistics<double>::constructed)
    .def_readonly("accepted", &MoveStatistics<double>::accepted)
    .def_readonly("allCounts", &MoveStatistics<double>::allCounts)
    .def_readonly("totalCounts", &MoveStatistics<double>::totalCounts)
    .def_readonly("totalConstructed", &MoveStatistics<double>::totalConstructed)
    .def_readonly("totalAccepted", &MoveStatistics<double>::totalAccepted)
    .def_readwrite("maxChange", &MoveStatistics<double>::maxChange)
    .def_readwrite("targetAcceptance", &MoveStatistics<double>::targetAcceptance)
    .def_readwrite("lowerLimit", &MoveStatistics<double>::lowerLimit)
    .def_readwrite("upperLimit", &MoveStatistics<double>::upperLimit)
    .def_readwrite("optimize", &MoveStatistics<double>::optimize);

  pybind11::class_<MoveStatistics<double3>>(m, "MoveStatisticsDouble3")
    //.def_property_readonly("counts", [](MoveStatistics<double3> const &s) { return std::vector<double>{s.counts.x, s.counts.y, s.counts.z};})
    .def_readonly("counts", &MoveStatistics<double3>::counts)
    .def_readonly("constructed", &MoveStatistics<double3>::constructed)
    .def_readonly("accepted", &MoveStatistics<double3>::accepted)
    .def_readonly("allCounts", &MoveStatistics<double3>::allCounts)
    .def_readonly("totalCounts", &MoveStatistics<double3>::totalCounts)
    .def_readonly("totalConstructed", &MoveStatistics<double3>::totalConstructed)
    .def_readonly("totalAccepted", &MoveStatistics<double3>::totalAccepted)
    .def_readwrite("maxChange", &MoveStatistics<double3>::maxChange)
    .def_readwrite("targetAcceptance", &MoveStatistics<double3>::targetAcceptance)
    .def_readwrite("lowerLimit", &MoveStatistics<double3>::lowerLimit)
    .def_readwrite("upperLimit", &MoveStatistics<double3>::upperLimit)
    .def_readwrite("optimize", &MoveStatistics<double3>::optimize);

  pybind11::class_<Move> move(m, "Move");
    move.def(pybind11::init());
  pybind11::native_enum<Move::Types>(move, "Types", "enum.IntEnum")
      .value("Translation", Move::Types::Translation)
      .value("RandomTranslation", Move::Types::RandomTranslation)
      .value("Rotation", Move::Types::Rotation)
      .value("RandomRotation", Move::Types::RandomRotation)
      .value("VolumeChange", Move::Types::VolumeChange)
      .value("ReinsertionCBMC", Move::Types::ReinsertionCBMC)
      .value("PartialReinsertionCBMC", Move::Types::PartialReinsertionCBMC)
      .value("IdentityChangeCBMC", Move::Types::IdentityChangeCBMC)
      .value("Swap", Move::Types::Swap)
      .value("SwapCBMC", Move::Types::SwapCBMC)
      .value("SwapCFCMC", Move::Types::SwapCFCMC)
      .value("SwapCBCFCMC", Move::Types::SwapCBCFCMC)
      .value("GibbsVolume", Move::Types::GibbsVolume)
      .value("GibbsSwapCBMC", Move::Types::GibbsSwapCBMC)
      .value("GibbsSwapCFCMC", Move::Types::GibbsSwapCFCMC)
      .value("Widom", Move::Types::Widom)
      .value("WidomCFCMC", Move::Types::WidomCFCMC)
      .value("WidomCBCFCMC", Move::Types::WidomCBCFCMC)
      .value("ParallelTempering", Move::Types::ParallelTempering)
      .value("HybridMC", Move::Types::HybridMC)
      .finalize();

  pybind11::class_<MCMoveStatistics>(m, "MCMoveStatistics")
    //.def("__getitem__", [](MCMoveStatistics &self, std::size_t index) { return self[index]; });
    .def("__getitem__", [](MCMoveStatistics &self, Move::Types i) { return self[i]; });

  // define before component init to prevent failing default argument
  pybind11::class_<Component> component(m, "Component");

  pybind11::native_enum<Component::Type>(component, "Type", "enum.IntEnum")
      .value("Adsorbate", Component::Type::Adsorbate)
      .value("Cation", Component::Type::Cation)
      .finalize();

  component
      .def(pybind11::init<std::size_t, const ForceField &, std::string, double, double, double, std::vector<Atom>,
                          const ConnectivityTable &, const Potentials::IntraMolecularPotentials &, std::size_t,
                          std::size_t, const MCMoveProbabilities &, std::optional<double>, bool, std::vector<double4>>(),
           pybind11::arg("componentId"), pybind11::arg("forceField"), pybind11::arg("componentName"),
           pybind11::arg("criticalTemperature"), pybind11::arg("criticalPressure"), pybind11::arg("acentricFactor"),
           pybind11::arg("definedAtoms") = std::vector<Atom>(), pybind11::arg("connectivityTable") = ConnectivityTable(),
           pybind11::arg("intraMolecularPotentials") = Potentials::IntraMolecularPotentials(),
           pybind11::arg("numberOfBlocks") = 5, pybind11::arg("numberOfLambdaBins") = 41,
           pybind11::arg("particleProbabilities") = MCMoveProbabilities(), pybind11::arg("fugacityCoefficient") = std::nullopt,
           pybind11::arg("thermodynamicIntegration") = false,
           pybind11::arg("blockingPockets") = std::vector<double4>())
      .def_readonly("name", &Component::name)
      .def_readwrite("lambdaHistogram", &Component::lambdaGC)
      .def_readonly("mc_moves_probabilities", &Component::mc_moves_probabilities)
      .def_readonly("mc_moves_statistics", &Component::mc_moves_statistics)
      .def_readwrite("blockingPockets", &Component::blockingPockets)
      .def_readwrite("averageRosenbluthWeights", &Component::averageRosenbluthWeights)
      .def("printStatus", &Component::printStatus)
      .def("__repr__", &Component::repr);


  pybind11::class_<Loadings>(m, "Loadings")
      .def(pybind11::init<std::size_t>())
      .def_readonly("numberOfMolecules", &Loadings::numberOfMolecules)
      .def_readonly("numberDensities", &Loadings::numberDensities)
      .def_readonly("inverseNumberDensities", &Loadings::inverseNumberDensities)
      .def("printStatus",
           static_cast<std::string (Loadings::*)(const Component &, std::optional<double>, std::optional<int3>) const>(
               &Loadings::printStatus))
      .def("printStatus", static_cast<std::string (Loadings::*)(const Component &, const Loadings &, const Loadings &,
                                                                std::optional<double>, std::optional<int3>) const>(
                              &Loadings::printStatus));

  pybind11::class_<SampleMovie>(m, "SampleMovie")
      .def(pybind11::init<std::size_t, std::size_t, bool>(),
           pybind11::arg("systemId"), pybind11::arg("sampleEvery"), pybind11::arg("restrictToBox") = true);

  pybind11::class_<PropertyLoading>(m, "PropertyLoading")
      .def(pybind11::init<std::size_t, std::size_t>())
      .def("result", &PropertyLoading::result)
      .def("averageLoadingNumberOfMolecules", &PropertyLoading::averageLoadingNumberOfMolecules)
      .def("writeAveragesStatistics", &PropertyLoading::writeAveragesStatistics)
      .def("__repr__", &PropertyLoading::repr);

  pybind11::class_<Potentials::EnergyFactor>(m, "EnergyFactor")
      .def_readonly("energy", &Potentials::EnergyFactor::energy)
      .def_readonly("dUdlambda", &Potentials::EnergyFactor::dUdlambda);

  pybind11::class_<EnergyStatus>(m, "EnergyStatus")
      .def_readwrite("totalEnergy", &EnergyStatus::totalEnergy)
      .def("__repr__", &EnergyStatus::repr);

  // convert result to units of Kelvin
  pybind11::class_<PropertyEnergy>(m, "PropertyEnergy")
      .def("result", [](PropertyEnergy& p) { return Units::EnergyToKelvin * p.result();})
      .def("__repr__", &PropertyEnergy::repr);

  pybind11::class_<AverageEnergyType>(m, "AverageEnergyType")
      .def(pybind11::init<double, double, double, double>(), 
           pybind11::arg("totalEnergy") = 0.0, pybind11::arg("VanDerWaalsEnergy") = 0.0,
           pybind11::arg("CoulombEnergy") = 0.0, pybind11::arg("polarizationEnergy") = 0.0)
      .def_readwrite("totalEnergy", &AverageEnergyType::totalEnergy)
      .def_readwrite("VanDerWaalsEnergy", &AverageEnergyType::VanDerWaalsEnergy)
      .def_readwrite("CoulombEnergy", &AverageEnergyType::CoulombEnergy)
      .def_readwrite("polarizationEnergy", &AverageEnergyType::polarizationEnergy);


  pybind11::class_<PropertyDensityGrid> property_energy_grid(m, "PropertyDensityGrid");

  pybind11::native_enum<PropertyDensityGrid::Normalization>(property_energy_grid, "Normalization", "enum.IntEnum")
      .value("Max", PropertyDensityGrid::Normalization::Max)
      .value("NumberDensity", PropertyDensityGrid::Normalization::NumberDensity)
      .finalize();

  pybind11::native_enum<PropertyDensityGrid::Binning>(property_energy_grid, "Binning", "enum.IntEnum")
      .value("Standard", PropertyDensityGrid::Binning::Standard)
      .value("Equitable", PropertyDensityGrid::Binning::Equitable)
      .finalize();

  property_energy_grid.def(pybind11::init<std::size_t, std::size_t, int3, std::size_t, 
                                          std::size_t, std::vector<std::size_t>,
                                          PropertyDensityGrid::Normalization, PropertyDensityGrid::Binning>(),
           pybind11::arg("numberOfFrameworks"), pybind11::arg("numberOfComponents"),
           pybind11::arg("numberOfGridPoints"), pybind11::arg("sampleEvery"), 
           pybind11::arg("writeEvery"), pybind11::arg("densityGridPseudoAtomsList"),
           pybind11::arg("normalizationType") = PropertyDensityGrid::Normalization::Max,
           pybind11::arg("binningMode") = PropertyDensityGrid::Binning::Standard);


  // results in units of Kelvin
  pybind11::class_<PropertyEnergyHistogram> energy_histogram(m, "PropertyEnergyHistogram");
    energy_histogram
      .def(pybind11::init<std::size_t, std::size_t, std::pair<double, double>, std::size_t, std::size_t>(),
                 pybind11::arg("numberOfBlocks"), pybind11::arg("numberOfBins"), pybind11::arg("valueRange"),
                 pybind11::arg("sampleEvery"), pybind11::arg("writeEvery"))
      .def("result", &PropertyEnergyHistogram::result);

  pybind11::class_<PropertyConventionalRadialDistributionFunction> conv_rdf(m, "PropertyConventionalRadialDistributionFunction");
    conv_rdf.def(pybind11::init<std::size_t, std::size_t, std::size_t, double, std::size_t, std::size_t>(),
                 pybind11::arg("numberOfBlocks"), pybind11::arg("numberOfPseudoAtoms"),
                 pybind11::arg("numberOfBins"), pybind11::arg("range"),
                 pybind11::arg("sampleEvery"), pybind11::arg("writeEvery"))
       .def("result", &PropertyConventionalRadialDistributionFunction::result);

  pybind11::class_<PropertyRadialDistributionFunction> rdf(m, "PropertyRadialDistributionFunction");
    rdf.def(pybind11::init<std::size_t, std::size_t, std::size_t, double, std::size_t, std::size_t>(),
            pybind11::arg("numberOfBlocks"), pybind11::arg("numberOfPseudoAtoms"),
            pybind11::arg("numberOfBins"), pybind11::arg("range"),
            pybind11::arg("sampleEvery"), pybind11::arg("writeEvery"))
       .def("result", &PropertyRadialDistributionFunction::result);

  pybind11::class_<PropertyNumberOfMoleculesEvolution>(m, "PropertyNumberOfMoleculesEvolution")
    .def(pybind11::init<std::size_t, std::size_t, std::size_t, std::optional<std::size_t>>(),
            pybind11::arg("numberOfCycles"), pybind11::arg("numberOfComponents"),
            pybind11::arg("sampleEvery"), pybind11::arg("writeEvery") = std::nullopt)
      .def_readonly("result", &PropertyNumberOfMoleculesEvolution::result);

  pybind11::class_<PropertyVolumeEvolution>(m, "PropertyVolumeEvolution")
    .def(pybind11::init<std::size_t, std::size_t, std::optional<std::size_t>>(),
            pybind11::arg("numberOfCycles"), pybind11::arg("sampleEvery"), 
            pybind11::arg("writeEvery") = std::nullopt)
      .def_readonly("result", &PropertyVolumeEvolution::result);

  pybind11::class_<PropertyWidom>(m, "PropertyWidom");

  pybind11::class_<System>(m, "System")
      .def(pybind11::init<std::size_t, ForceField, std::optional<SimulationBox>, bool, double, std::optional<double>, double,
                          std::optional<Framework>, std::vector<Component>, std::vector<std::vector<double3>>,
                          std::vector<std::size_t>, std::size_t, MCMoveProbabilities>(),
           pybind11::arg("systemId"), pybind11::arg("forceField"), pybind11::arg("simulationBox") = std::nullopt,
           pybind11::arg("hasExternalField") = false, pybind11::arg("externalTemperature") = 298.0, 
           pybind11::arg("externalPressure") = std::nullopt, pybind11::arg("heliumVoidFraction") = 0.0,
           pybind11::arg("frameworkComponents") = std::nullopt,
           pybind11::arg("components") = std::vector<Component>(), pybind11::arg("initialPositions") = std::vector<std::vector<double3>>(),
           pybind11::arg("initialNumberOfMolecules") = std::vector<std::size_t>(), pybind11::arg("numberOfBlocks") = 5,
           pybind11::arg("systemProbabilities") = MCMoveProbabilities())
      .def("computeTotalEnergies", &System::computeTotalEnergies)
      .def("frameworkMass", &System::frameworkMass)
      .def_readonly("inputPressure", &System::input_pressure)
      .def_readonly("components", &System::components)
      .def_readonly("loadings", &System::loadings)
      .def_readwrite("averageLoadings", &System::averageLoadings)
      .def_readwrite("samplePDBMovie", &System::samplePDBMovie)
      .def_readwrite("averageEnergyHistogram", &System::averageEnergyHistogram)
      .def_readwrite("averageEnergies", &System::averageEnergies)
      .def_readwrite("densityGrid", &System::propertyDensityGrid)
      .def_readwrite("conventionalRadialDistributionFunction", &System::propertyConventionalRadialDistributionFunction)
      .def_readwrite("radialDistributionFunction", &System::propertyRadialDistributionFunction)
      .def_readwrite("propertyNumberOfMoleculesEvolution", &System::propertyNumberOfMoleculesEvolution)
      .def_readwrite("propertyVolumeEvolution", &System::propertyVolumeEvolution)
      .def_readwrite("atomData", &System::atomData)
      .def("writeMCMoveStatistics", &System::writeMCMoveStatistics)
      .def("__repr__", &System::repr);


  pybind11::class_<MonteCarlo>(m, "MonteCarlo")
      .def(pybind11::init<std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t,
                        std::vector<System> &, std::optional<std::size_t>, std::size_t, bool>(),
         pybind11::arg("numberOfCycles"), pybind11::arg("numberOfInitializationCycles"),
         pybind11::arg("numberOfEquilibrationCycles") = 0, pybind11::arg("printEvery") = 5000,
         pybind11::arg("writeBinaryRestartEvery") = 5000, pybind11::arg("rescaleWangLandauEvery") = 1000,
         pybind11::arg("optimizeMCMovesEvery") = 100, pybind11::arg("systems") = std::vector<System>(),
         pybind11::arg("randomSeed") = std::nullopt, pybind11::arg("numberOfBlocks") = 5,
         pybind11::arg("outputToFiles") = false)
      .def(pybind11::init<InputReader &>(), pybind11::arg("inputReader"))
      .def("run", &MonteCarlo::run)
      .def("initialize", &MonteCarlo::initialize, pybind11::arg("call_back_function") = pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def("equilibrate", &MonteCarlo::equilibrate, pybind11::arg("call_back_function") = pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def("production", &MonteCarlo::production, pybind11::arg("call_back_function") = pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def("cycle", &MonteCarlo::performCycle)
      .def_readonly("systems", &MonteCarlo::systems);

  pybind11::class_<MolecularDynamics>(m, "MolecularDynamics")
    .def(pybind11::init<std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t,
                        std::vector<System> &, std::optional<std::size_t>, std::size_t, bool>(),
         pybind11::arg("numberOfCycles"), pybind11::arg("numberOfInitializationCycles"),
         pybind11::arg("numberOfEquilibrationCycles") = 0, pybind11::arg("printEvery") = 5000,
         pybind11::arg("writeBinaryRestartEvery") = 5000, pybind11::arg("rescaleWangLandauEvery") = 1000,
         pybind11::arg("optimizeMCMovesEvery") = 100, pybind11::arg("systems") = std::vector<System>(),
         pybind11::arg("randomSeed") = std::nullopt, pybind11::arg("numberOfBlocks") = 5,
         pybind11::arg("outputToFiles") = false)
      .def("run", &MolecularDynamics::run)
      .def("initialize", &MolecularDynamics::initialize, pybind11::arg("call_back_function") = pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def("equilibrate", &MolecularDynamics::equilibrate, pybind11::arg("call_back_function") = pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def("production", &MolecularDynamics::production, pybind11::arg("call_back_function") = pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def_readonly("systems", &MolecularDynamics::systems);

}

namespace pybind11 {
namespace detail {

template <>
struct type_caster<int3> 
{
  // This macro inserts a lot of boilerplate code and sets the type hint.
  // `io_name` is used to specify different type hints for arguments and return values.
  // The signature of our negate function would then look like:
  // `negate(Sequence[float]) -> tuple[float, float]`
  PYBIND11_TYPE_CASTER(int3, io_name("Sequence[int]", "tuple[int, int, int]"));

  // C++ -> Python: convert `Point2D` to `tuple[float, float]`. The second and third arguments
  // are used to indicate the return value policy and parent object (for
  // return_value_policy::reference_internal) and are often ignored by custom casters.
  // The return value should reflect the type hint specified by the second argument of `io_name`.
  static handle
  cast(const int3 &number, return_value_policy /*policy*/, handle /*parent*/) 
  {
      return pybind11::make_tuple(number.x, number.y, number.z).release();
  }

  // Python -> C++: convert a `PyObject` into a `int3` and return false upon failure. The
  // second argument indicates whether implicit conversions should be allowed.
  // The accepted types should reflect the type hint specified by the first argument of
  // `io_name`.
  bool load(handle src, bool /*convert*/)
  {
    // Check if handle is a Sequence
    if (!pybind11::isinstance<pybind11::sequence>(src)) 
    {
      return false;
    }
    auto seq = pybind11::reinterpret_borrow<pybind11::sequence>(src);
    // Check if exactly two values are in the Sequence
    if (seq.size() != 3) 
    {
      return false;
    }
    // Check if each element is either an int
    for (auto item : seq)
    {
      if (!pybind11::isinstance<pybind11::int_>(item))
      {
        return false;
      }
    }
    value.x = seq[0].cast<std::int32_t>();
    value.y = seq[1].cast<std::int32_t>();
    value.z = seq[2].cast<std::int32_t>();
    return true;
  }
};

template <>
struct type_caster<double3> 
{
  // This macro inserts a lot of boilerplate code and sets the type hint.
  // `io_name` is used to specify different type hints for arguments and return values.
  // The signature of our negate function would then look like:
  // `negate(Sequence[float]) -> tuple[float, float]`
  PYBIND11_TYPE_CASTER(double3, io_name("Sequence[float]", "tuple[float, float, float]"));

  // C++ -> Python: convert `Point2D` to `tuple[float, float]`. The second and third arguments
  // are used to indicate the return value policy and parent object (for
  // return_value_policy::reference_internal) and are often ignored by custom casters.
  // The return value should reflect the type hint specified by the second argument of `io_name`.
  static handle
  cast(const double3 &number, return_value_policy /*policy*/, handle /*parent*/) {
      return pybind11::make_tuple(number.x, number.y, number.z).release();
  }

  // Python -> C++: convert a `PyObject` into a `double3` and return false upon failure. The
  // second argument indicates whether implicit conversions should be allowed.
  // The accepted types should reflect the type hint specified by the first argument of
  // `io_name`.
  bool load(handle src, bool /*convert*/) 
  {
    if (pybind11::isinstance<pybind11::float_>(src) || pybind11::isinstance<pybind11::int_>(src))
    {
      value.x = src.cast<double>();
      value.y = src.cast<double>();
      value.z = src.cast<double>();
      return true;
    }

    // Check if handle is a Sequence
    if (pybind11::isinstance<pybind11::sequence>(src)) 
    {
      auto seq = pybind11::reinterpret_borrow<pybind11::sequence>(src);

      // Check if exactly one value are in the Sequence
      if (seq.size() == 1)
      {
        if (!pybind11::isinstance<pybind11::float_>(seq[0]) && !pybind11::isinstance<pybind11::int_>(seq[0]))
        {
          return false;
        }

        value.x = seq[0].cast<double>();
        value.y = seq[0].cast<double>();
        value.z = seq[0].cast<double>();
        return true;
      }

      // Check if exactly three values are in the Sequence
      if (seq.size() == 3)
      {
        // Check if each element is either a float or an int
        for (auto item : seq) 
        {
          if (!pybind11::isinstance<pybind11::float_>(item) && !pybind11::isinstance<pybind11::int_>(item))
          {
            return false;
          }
        }
        value.x = seq[0].cast<double>();
        value.y = seq[1].cast<double>();
        value.z = seq[2].cast<double>();
        return true;
      }
    }

    return false;
  }
};

template <>
struct type_caster<double4> 
{
  // This macro inserts a lot of boilerplate code and sets the type hint.
  // `io_name` is used to specify different type hints for arguments and return values.
  // The signature of our negate function would then look like:
  // `negate(Sequence[float]) -> tuple[float, float]`
  PYBIND11_TYPE_CASTER(double4, io_name("Sequence[float]", "tuple[float, float, float,float]"));

  // C++ -> Python: convert `Point2D` to `tuple[float, float]`. The second and third arguments
  // are used to indicate the return value policy and parent object (for
  // return_value_policy::reference_internal) and are often ignored by custom casters.
  // The return value should reflect the type hint specified by the second argument of `io_name`.
  static handle
  cast(const double4 &number, return_value_policy /*policy*/, handle /*parent*/) {
      return pybind11::make_tuple(number.x, number.y, number.z, number.w).release();
  }

  // Python -> C++: convert a `PyObject` into a `double4` and return false upon failure. The
  // second argument indicates whether implicit conversions should be allowed.
  // The accepted types should reflect the type hint specified by the first argument of
  // `io_name`.
  bool load(handle src, bool /*convert*/) 
  {
    if (pybind11::isinstance<pybind11::float_>(src) || pybind11::isinstance<pybind11::int_>(src))
    {
      value.x = src.cast<double>();
      value.y = src.cast<double>();
      value.z = src.cast<double>();
      value.w = src.cast<double>();
      return true;
    }

    // Check if handle is a Sequence
    if (pybind11::isinstance<pybind11::sequence>(src)) 
    {
      auto seq = pybind11::reinterpret_borrow<pybind11::sequence>(src);

      // Check if exactly one value are in the Sequence
      if (seq.size() == 1)
      {
        if (!pybind11::isinstance<pybind11::float_>(seq[0]) && !pybind11::isinstance<pybind11::int_>(seq[0]))
        {
          return false;
        }

        value.x = seq[0].cast<double>();
        value.y = seq[0].cast<double>();
        value.z = seq[0].cast<double>();
        value.w = seq[0].cast<double>();
        return true;
      }

      // Check if exactly three values are in the Sequence
      if (seq.size() == 4)
      {
        // Check if each element is either a float or an int
        for (auto item : seq) 
        {
          if (!pybind11::isinstance<pybind11::float_>(item) && !pybind11::isinstance<pybind11::int_>(item))
          {
            return false;
          }
        }
        value.x = seq[0].cast<double>();
        value.y = seq[1].cast<double>();
        value.z = seq[2].cast<double>();
        value.w = seq[3].cast<double>();
        return true;
      }
    }

    return false;
  }
};

} // namespace detail
} // namespace pybind11
