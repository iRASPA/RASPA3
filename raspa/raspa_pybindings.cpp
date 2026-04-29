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
import double2;
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
import loading_data;
import energy_factor;
import energy_status;
import sample_movies;
import property_energy;
import average_energy_type;
import property_energy_histogram;
import property_lambda_probability_histogram;
import property_density_grid;
import property_number_of_molecules_evolution;
import property_volume_evolution;
import property_widom;
import property_gibbs_widom;
import connectivity_table;
import intra_molecular_potentials;
import pressure_data;
import property_pressure;
import widom_data;
import equation_of_states;
import property_conserved_energy_evolution;
import thermostat;
import enthalpy_of_adsorption_data;
import property_enthalpy;
import property_conventional_rdf;
import property_rdf;
import mean_squared_displacement_data;
import property_msd;
import velocity_autocorrelation_function_data;
import property_vacf;


template <typename U, typename T>
std::vector<T> operator*(std::vector<U> v, const T& scalar)
{
  for (auto& element : v)
  {
    element *= scalar;
  }
  return v;
}

template <typename U, typename T>
std::vector<T> operator*(const U& scalar, std::vector<T> v)
{
  for (auto& element : v)
  {
    element *= scalar;
  }
  return v;
}


template<typename T>
std::pair<T, T> operator*(const double& a, const std::pair<T, T>& b)
{
  return std::make_pair(a * b.first, a * b.second);
}

template<typename T>
std::pair<T, T> operator*(const double2& a, const std::pair<T, T>& b)
{
  return std::make_pair(a.x * b.first, a.y * b.second);
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

  pybind11::class_<EquationOfState> eos(m, "EquationOfState");

    pybind11::class_<EquationOfState::FluidInput>(eos, "FluidInput")
        .def(pybind11::init<double, double, double, double, bool>(), 
            pybind11::arg("critical_temperature"), 
            pybind11::arg("critical_pressure"),
            pybind11::arg("acentric_factor"), 
            pybind11::arg("mol_fraction") = 1.0,
            pybind11::arg("swappable") = true)
        .def_readwrite("critical_temperature", &EquationOfState::FluidInput::criticalTemperature)
        .def_readwrite("critical_pressure", &EquationOfState::FluidInput::criticalPressure)
        .def_readwrite("acentric_factor", &EquationOfState::FluidInput::acentricFactor)
        .def_readwrite("mol_fraction", &EquationOfState::FluidInput::molFraction)
        .def_readwrite("swappable", &EquationOfState::FluidInput::swappable);

    pybind11::native_enum<EquationOfState::Type>(eos, "EquationOfStateType", "enum.IntEnum")
        .value("PENG_ROBINSON", EquationOfState::Type::PengRobinson)
        .value("PENG_ROBINSON_GASEM", EquationOfState::Type::PengRobinsonGasem)
        .value("SOAVE_REDLICH_KWONG", EquationOfState::Type::SoaveRedlichKwong)
        .finalize();
    pybind11::native_enum<EquationOfState::MixingRules>(eos, "MixingRules", "enum.IntEnum")
        .value("VAN_DER_WAALS", EquationOfState::MixingRules::VanDerWaals)
        .finalize();

  pybind11::class_<EquationOfState::FluidResult>(eos, "FluidResult")
      .def(pybind11::init<double, std::optional<double>, EquationOfState::FluidState>(), 
          pybind11::arg("compressibility"), 
          pybind11::arg("fugacity_coefficient"),
          pybind11::arg("fluid_state"))
      .def_readwrite("compressibility", &EquationOfState::FluidResult::compressibility)
      .def_readwrite("fugacity_coefficient", &EquationOfState::FluidResult::fugacityCoefficient)
      .def_readwrite("fluid_state", &EquationOfState::FluidResult::fluidState);

    eos.def(pybind11::init<EquationOfState::Type, EquationOfState::MixingRules,
                  double, double, const SimulationBox &, double,
                  std::vector<Component> &>())
       .def_static("compute_fluid_properties", &EquationOfState::computeFluidProperties,
         pybind11::arg("temperature"), 
         pybind11::arg("pressure"),
         pybind11::arg("properties"),
         pybind11::arg("type") = EquationOfState::Type::PengRobinson,
         pybind11::arg("mixing_rules") = EquationOfState::MixingRules::VanDerWaals);

  pybind11::class_<RandomNumber>(m, "RandomNumber")
    .def(pybind11::init<std::size_t>(), 
         pybind11::arg("seed") = 12);

  pybind11::class_<Atom>(m, "Atom")
      .def(pybind11::init<>())
      .def(pybind11::init<double3, double, double, std::uint32_t, std::uint16_t, std::uint8_t, std::uint8_t, std::uint8_t>(),
           pybind11::arg("position"), 
           pybind11::arg("charge") = 0.0,
           pybind11::arg("scaling") = 1.0,
           pybind11::arg("molecule_id") = 0, 
           pybind11::arg("type") = 0, 
           pybind11::arg("component_id") = 0,
           pybind11::arg("group_id") = false, 
           pybind11::arg("is_fractional") = false)
      .def_readwrite("position", &Atom::position)
      .def("__repr__", &Atom::repr);

  pybind11::class_<SimulationBox> simulationBox(m, "SimulationBox");
  pybind11::native_enum<SimulationBox::Type>(simulationBox, "SimulationBoxType", "enum.IntEnum")
      .value("RECTANGULAR", SimulationBox::Type::Rectangular)
      .value("TRICLINIC", SimulationBox::Type::Triclinic)
      .finalize();
  simulationBox
      .def(pybind11::init<double, double, double>(), 
           pybind11::arg("a"), 
           pybind11::arg("b"), 
           pybind11::arg("c"))
      .def(pybind11::init<double, double, double, double, double, double>(), 
           pybind11::arg("a"), 
           pybind11::arg("b"),
           pybind11::arg("c"), 
           pybind11::arg("alpha"), 
           pybind11::arg("beta"), 
           pybind11::arg("gamma"))
      .def(pybind11::init<double, double, double, double, double, double, SimulationBox::Type>(), 
           pybind11::arg("a"),
           pybind11::arg("b"), 
           pybind11::arg("c"), 
           pybind11::arg("alpha"), 
           pybind11::arg("beta"),
           pybind11::arg("gamma"), 
           pybind11::arg("type"))
      .def(pybind11::init<double3x3, SimulationBox::Type>())
      .def_readwrite("type", &SimulationBox::type)
      .def_readonly("length_a", &SimulationBox::lengthA)
      .def_readonly("length_b", &SimulationBox::lengthB)
      .def_readonly("length_c", &SimulationBox::lengthC)
      .def_readonly("volume", &SimulationBox::volume);

  pybind11::class_<VDWParameters>(m, "VDWParameters")
      .def(pybind11::init<double, double>(), pybind11::arg("epsilon"), pybind11::arg("sigma"));

  pybind11::class_<PseudoAtom>(m, "PseudoAtom")
      .def(pybind11::init<std::string, bool, double, double, double, std::size_t, bool, std::string>(),
           pybind11::arg("name"), 
           pybind11::arg("framework_type"), 
           pybind11::arg("mass"), 
           pybind11::arg("charge"),
           pybind11::arg("polarizability") = 0.0, 
           pybind11::arg("atomic_number") = 1, 
           pybind11::arg("print_to_pdb") = true,
           pybind11::arg("source") = "");

  pybind11::class_<ForceField> forceField(m, "ForceField");
  pybind11::native_enum<ForceField::MixingRule>(forceField, "MixingRule", "enum.IntEnum")
      .value("LORENTZ_BERTHELOT", ForceField::MixingRule::Lorentz_Berthelot)
      .value("JORGENSEN", ForceField::MixingRule::Jorgensen)
      .finalize();
  forceField.def(pybind11::init<std::vector<PseudoAtom>, std::vector<VDWParameters>, ForceField::MixingRule, 
                 double, double, double, bool, bool, bool>(),
           pybind11::arg("pseudo_atoms"), 
           pybind11::arg("parameters"), 
           pybind11::arg("mixing_rule") = ForceField::MixingRule::Lorentz_Berthelot,
           pybind11::arg("cutoff_framework_vdw") = 12.0, 
           pybind11::arg("cutoff_molecule_vdw") = 12.0,
           pybind11::arg("cutoff_coulomb") = 12.0, 
           pybind11::arg("shifted") = true,
           pybind11::arg("tail_corrections") = false, 
           pybind11::arg("use_charge") = true)
      .def("find_pseudo_atom", static_cast<std::optional<std::size_t> (ForceField::*)(const std::string &) const>(&ForceField::findPseudoAtom))
      .def("__repr__", &ForceField::repr)
      .def_readonly("pseudo_atoms", &ForceField::pseudoAtoms)
      .def_readonly("vdw_parameters", &ForceField::data)
      .def_readwrite("use_charge", &ForceField::useCharge);

  pybind11::class_<CIFReader> cifReader(m, "CIFReader");
  cifReader
      .def("expand_defined_atoms_to_unit_cell", &CIFReader::expandDefinedAtomsToUnitCell);

  pybind11::class_<Framework> framework(m, "Framework");

  pybind11::native_enum<CIFReader::UseChargesFrom>(framework, "UseChargesFrom", "enum.IntEnum")
      .value("PSEUDO_ATOMS", CIFReader::UseChargesFrom::PseudoAtoms)
      .value("CIF_FILE", CIFReader::UseChargesFrom::CIF_File)
      .value("CHARGE_EQUILIBRATION", CIFReader::UseChargesFrom::ChargeEquilibration)
      .finalize();

  framework
      .def(pybind11::init<const ForceField &, std::string, SimulationBox, std::size_t, 
                          const std::vector<Atom> &, int3>(),
           pybind11::arg("force_field"), 
           pybind11::arg("component_name"),
           pybind11::arg("simulation_box"), 
           pybind11::arg("space_group_hall_number"), 
           pybind11::arg("defined_atoms"), 
           pybind11::arg("number_of_unit_cells"))
      .def_readonly("name", &Framework::name)
      .def("print", &Framework::printStatus)
      .def("__repr__", &Framework::repr);

  pybind11::class_<MCMoveProbabilities>(m, "MCMoveProbabilities")
      .def(pybind11::init<double, double, double, double, double, double, double, double, double, double, double,
                          double, double, double, double, double, double, double, double, double>(),
           pybind11::arg("translation_probability") = 0.0, 
           pybind11::arg("random_translation_probability") = 0.0,
           pybind11::arg("rotation_probability") = 0.0, 
           pybind11::arg("random_rotation_probability") = 0.0,
           pybind11::arg("volume_change_probability") = 0.0, 
           pybind11::arg("reinsertion_cbmc_probability") = 0.0,
           pybind11::arg("partial_reinsertion_cbmc_probability") = 0.0,
           pybind11::arg("identity_change_probability") = 0.0, 
           pybind11::arg("swap_probability") = 0.0,
           pybind11::arg("swap_cbmc_probability") = 0.0, 
           pybind11::arg("swap_cfcmc_probability") = 0.0,
           pybind11::arg("swap_cbcfcmc_probability") = 0.0, 
           pybind11::arg("gibbs_volume_change_probability") = 0.0,
           pybind11::arg("gibbs_swap_cbmc_probability") = 0.0, 
           pybind11::arg("gibbs_swap_cfcmc_probability") = 0.0,
           pybind11::arg("widom_probability") = 0.0, 
           pybind11::arg("widom_cfcmc_probability") = 0.0,
           pybind11::arg("widom_cbcfcmc_probability") = 0.0, 
           pybind11::arg("parallel_tempering_probability") = 0.0,
           pybind11::arg("hybrid_mc_probability") = 0.0)
      .def("join", &MCMoveProbabilities::join)
      .def("__repr__", &MCMoveProbabilities::repr);


  pybind11::class_<PropertyLambdaProbabilityHistogram>(m, "PropertyLambdaProbabilityHistogram")
      .def(pybind11::init<>())
      // dimensionless
      .def_readwrite("bias_factor", &PropertyLambdaProbabilityHistogram::biasFactor)
      // dimensionless
      .def_readwrite("histogram", &PropertyLambdaProbabilityHistogram::histogram)
       // convert result to units of Kelvin
      .def("average_dudlambda", [](PropertyLambdaProbabilityHistogram& p) 
          { return Units::EnergyToKelvin * p.averageDuDlambda();})
      // dimensionless
      .def("result", &PropertyLambdaProbabilityHistogram::result);


  pybind11::class_<ConnectivityTable>(m, "ConnectivityTable")
    .def(pybind11::init<>());

  pybind11::class_<Potentials::IntraMolecularPotentials>(m, "IntraMolecularPotentials")
    .def(pybind11::init<>());

  pybind11::class_<MoveStatistics<double>>(m, "MoveStatisticsDouble")
    .def_readonly("counts", &MoveStatistics<double>::counts)
    .def_readonly("constructed", &MoveStatistics<double>::constructed)
    .def_readonly("accepted", &MoveStatistics<double>::accepted)
    .def_readonly("all_counts", &MoveStatistics<double>::allCounts)
    .def_readonly("total_counts", &MoveStatistics<double>::totalCounts)
    .def_readonly("total_constructed", &MoveStatistics<double>::totalConstructed)
    .def_readonly("total_accepted", &MoveStatistics<double>::totalAccepted)
    .def_readwrite("max_change", &MoveStatistics<double>::maxChange)
    .def_readwrite("target_acceptance", &MoveStatistics<double>::targetAcceptance)
    .def_readwrite("lower_limit", &MoveStatistics<double>::lowerLimit)
    .def_readwrite("upper_limit", &MoveStatistics<double>::upperLimit)
    .def_readwrite("optimize", &MoveStatistics<double>::optimize);

  pybind11::class_<MoveStatistics<double3>>(m, "MoveStatisticsDouble3")
    .def_readonly("counts", &MoveStatistics<double3>::counts)
    .def_readonly("constructed", &MoveStatistics<double3>::constructed)
    .def_readonly("accepted", &MoveStatistics<double3>::accepted)
    .def_readonly("all_counts", &MoveStatistics<double3>::allCounts)
    .def_readonly("total_counts", &MoveStatistics<double3>::totalCounts)
    .def_readonly("total_constructed", &MoveStatistics<double3>::totalConstructed)
    .def_readonly("total_accepted", &MoveStatistics<double3>::totalAccepted)
    .def_readwrite("max_change", &MoveStatistics<double3>::maxChange)
    .def_readwrite("target_acceptance", &MoveStatistics<double3>::targetAcceptance)
    .def_readwrite("lower_limit", &MoveStatistics<double3>::lowerLimit)
    .def_readwrite("upper_limit", &MoveStatistics<double3>::upperLimit)
    .def_readwrite("optimize", &MoveStatistics<double3>::optimize);

  pybind11::class_<Move> move(m, "Move");
    move.def(pybind11::init());
  pybind11::native_enum<Move::Types>(move, "Types", "enum.IntEnum")
      .value("TRANSLATION", Move::Types::Translation)
      .value("RANDOM_TRANSLATION", Move::Types::RandomTranslation)
      .value("ROTATION", Move::Types::Rotation)
      .value("RANDOM_ROTATION", Move::Types::RandomRotation)
      .value("VOLUME_CHANGE", Move::Types::VolumeChange)
      .value("REINSERTION_CBMC", Move::Types::ReinsertionCBMC)
      .value("PARTIAL_REINSERTION_CBMC", Move::Types::PartialReinsertionCBMC)
      .value("IDENTITY_CHANGE_CBMC", Move::Types::IdentityChangeCBMC)
      .value("SWAP", Move::Types::Swap)
      .value("SWAP_CBMC", Move::Types::SwapCBMC)
      .value("SWAP_CFCMC", Move::Types::SwapCFCMC)
      .value("SWAP_CBCFCMC", Move::Types::SwapCBCFCMC)
      .value("GIBBS_VOLUME", Move::Types::GibbsVolume)
      .value("GIBBS_SWAP_CBMC", Move::Types::GibbsSwapCBMC)
      .value("GIBBS_SWAP_CFCMC", Move::Types::GibbsSwapCFCMC)
      .value("WIDOM", Move::Types::Widom)
      .value("WIDOM_CFCMC", Move::Types::WidomCFCMC)
      .value("WIDOM_CBCFCMC", Move::Types::WidomCBCFCMC)
      .value("PARALLEL_TEMPERING", Move::Types::ParallelTempering)
      .value("HYBRID_MC", Move::Types::HybridMC)
      .finalize();

  pybind11::class_<MCMoveStatistics>(m, "MCMoveStatistics")
      .def("__getitem__", [](MCMoveStatistics &self, Move::Types i) { return self[i]; })
      .def("__repr__", &MCMoveStatistics::repr);

  pybind11::class_<WidomData>(m, "WidomData")
      .def_readonly("total", &WidomData::total)
      .def_readonly("excess", &WidomData::excess)
      .def_readonly("ideal_gas", &WidomData::idealGas);

  pybind11::class_<PropertyWidom>(m, "PropertyWidom")
      .def("result", &PropertyWidom::result)
       // convert result to units of Kelvin
      .def("chemical_potential_result", [](PropertyWidom& p, double T) 
          { return Units::EnergyToKelvin * p.chemicalPotentialResult(1.0 / (Units::KB * T));}, pybind11::arg("temperature"))
       // convert result to units of Pascal
      .def("fugacity_result", [](PropertyWidom& p, double T) 
          { return Units::PressureConversionFactor * p.fugacityResult(1.0 / (Units::KB * T));}, pybind11::arg("temperature"));

  pybind11::class_<PropertyGibbsWidom>(m, "PropertyGibbsWidom")
      .def("result", &PropertyGibbsWidom::result)
       // convert result to units of Kelvin
      .def("chemical_potential_result", [](PropertyGibbsWidom& p, double T) 
          { return Units::EnergyToKelvin * p.chemicalPotentialResult(1.0 / (Units::KB * T));}, pybind11::arg("temperature"))
       // convert result to units of Pascal
      .def("fugacity_result", [](PropertyGibbsWidom& p, double T) 
          { return Units::PressureConversionFactor * p.fugacityResult(1.0 / (Units::KB * T));}, pybind11::arg("temperature"));

  pybind11::class_<Component> component(m, "Component");
  pybind11::native_enum<Component::Type>(component, "Type", "enum.IntEnum")
      .value("ADSORBATE", Component::Type::Adsorbate)
      .value("CATION", Component::Type::Cation)
      .finalize();
  component
      .def(pybind11::init<const ForceField &, std::string, double, double, double, std::vector<Atom>,
                          const ConnectivityTable &, const Potentials::IntraMolecularPotentials &, std::size_t,
                          std::size_t, const MCMoveProbabilities &, std::optional<double>, bool, std::vector<double4>>(),
           pybind11::arg("force_field"), 
           pybind11::arg("component_name"),
           pybind11::arg("critical_temperature"), 
           pybind11::arg("critical_pressure"), 
           pybind11::arg("acentric_factor"),
           pybind11::arg("defined_atoms") = std::vector<Atom>(), 
           pybind11::arg("connectivity_table") = ConnectivityTable(),
           pybind11::arg("intra_molecular_potentials") = Potentials::IntraMolecularPotentials(),
           pybind11::arg("number_of_blocks") = 5,
           pybind11::arg("number_of_lambda_bins") = 41,
           pybind11::arg("particle_probabilities") = MCMoveProbabilities(), 
           pybind11::arg("fugacity_coefficient") = std::nullopt,
           pybind11::arg("thermodynamic_integration") = false,
           pybind11::arg("blocking_pockets") = std::vector<double4>())
      .def_readonly("name", &Component::name)
      .def_readwrite("lambda_histogram", &Component::lambdaGC)
      .def_readonly("mc_moves_probabilities", &Component::mc_moves_probabilities)
      .def_readonly("mc_moves_statistics", &Component::mc_moves_statistics)
      .def_readwrite("blocking_pockets", &Component::blockingPockets)
      .def_readwrite("average_rosenbluth_weights", &Component::averageRosenbluthWeights)
      .def_readwrite("average_gibbs_rosenbluth_weights", &Component::averageGibbsRosenbluthWeights)
      .def("print_status", &Component::printStatus)
      .def("__repr__", &Component::repr);


  pybind11::class_<LoadingData>(m, "LoadingData")
      .def(pybind11::init<std::size_t>())
      .def_readonly("number_of_molecules", &LoadingData::numberOfMolecules)
      .def_readonly("number_densities", &LoadingData::numberDensities)
      .def_readonly("inverse_number_densities", &LoadingData::inverseNumberDensities);

  pybind11::class_<SampleMovie>(m, "SampleMovie")
      .def(pybind11::init<std::size_t, std::size_t, bool>(),
           pybind11::arg("system_id"), 
           pybind11::arg("sample_every"), 
           pybind11::arg("restrict_to_box") = true);

  pybind11::class_<PropertyLoading>(m, "PropertyLoading")
      .def(pybind11::init<std::size_t, std::size_t>())
      .def("result", &PropertyLoading::result)
      .def("average_loading_number_of_molecules", &PropertyLoading::averageLoadingNumberOfMolecules)
      .def("write_averages_statistics", &PropertyLoading::writeAveragesStatistics)
      .def("__repr__", &PropertyLoading::repr);

  pybind11::class_<Potentials::EnergyFactor>(m, "EnergyFactor")
      .def_readonly("energy", &Potentials::EnergyFactor::energy)
      .def_readonly("dudlambda", &Potentials::EnergyFactor::dUdlambda);

  pybind11::class_<EnergyStatus>(m, "EnergyStatus")
      .def_readwrite("total_energy", &EnergyStatus::totalEnergy)
      .def("__repr__", &EnergyStatus::repr);

  pybind11::class_<PropertyEnergy>(m, "PropertyEnergy")
       // convert result to units of Kelvin
      .def("result", [](PropertyEnergy& p) { return Units::EnergyToKelvin * p.result();})
      .def("__repr__", &PropertyEnergy::repr);

  pybind11::class_<PressureData>(m, "PressureData")
      .def(pybind11::init<>())
      .def_readonly("total_pressure", &PressureData::totalPressure)
      .def_readonly("excess_pressure", &PressureData::excessPressure)
      .def_readonly("ideal_gas_pressure", &PressureData::idealGasPressure);

  pybind11::class_<PropertyPressure>(m, "PropertyPressure")
      .def("result", [](PropertyPressure& p) { return Units::PressureConversionFactor * p.result();});
      //.def("__repr__", &Pressures::repr);

  pybind11::class_<AverageEnergyType>(m, "AverageEnergyType")
      .def(pybind11::init<double, double, double, double>(), 
           pybind11::arg("total_energy") = 0.0, 
           pybind11::arg("van_der_waals_energy") = 0.0,
           pybind11::arg("coulomb_energy") = 0.0, 
           pybind11::arg("polarization_energy") = 0.0)
      .def_readwrite("total_energy", &AverageEnergyType::totalEnergy)
      .def_readwrite("van_der_waals_energy", &AverageEnergyType::VanDerWaalsEnergy)
      .def_readwrite("coulomb_energy", &AverageEnergyType::CoulombEnergy)
      .def_readwrite("polarization_energy", &AverageEnergyType::polarizationEnergy);


  pybind11::class_<PropertyDensityGrid> property_energy_grid(m, "PropertyDensityGrid");

  pybind11::native_enum<PropertyDensityGrid::Normalization>(property_energy_grid, "Normalization", "enum.IntEnum")
      .value("MAX", PropertyDensityGrid::Normalization::Max)
      .value("NUMBER_DENSITY", PropertyDensityGrid::Normalization::NumberDensity)
      .finalize();

  pybind11::native_enum<PropertyDensityGrid::Binning>(property_energy_grid, "Binning", "enum.IntEnum")
      .value("STANDARD", PropertyDensityGrid::Binning::Standard)
      .value("EQUITABLE", PropertyDensityGrid::Binning::Equitable)
      .finalize();

  property_energy_grid.def(pybind11::init<std::size_t, std::size_t, int3, std::size_t, 
                                          std::size_t, std::vector<std::size_t>,
                                          PropertyDensityGrid::Normalization, PropertyDensityGrid::Binning>(),
           pybind11::arg("number_of_frameworks"), 
           pybind11::arg("number_of_components"),
           pybind11::arg("number_of_grid_points"), 
           pybind11::arg("sample_every"), 
           pybind11::arg("write_every"), 
           pybind11::arg("density_grid_pseudo_atoms_list"),
           pybind11::arg("normalization_type") = PropertyDensityGrid::Normalization::Max,
           pybind11::arg("binningMode") = PropertyDensityGrid::Binning::Standard);


  pybind11::class_<PropertyEnergyHistogram> energy_histogram(m, "PropertyEnergyHistogram");
    energy_histogram
      .def(pybind11::init<std::size_t, std::size_t, std::pair<double, double>, std::size_t, std::size_t>(),
                 pybind11::arg("number_of_blocks"), 
                 pybind11::arg("number_of_bins"), 
                 pybind11::arg("value_range"),
                 pybind11::arg("sample_every"), 
                 pybind11::arg("write_every"))
      // convert result to units of Kelvin
      .def("result", [](PropertyEnergyHistogram& p) { auto [bins, average, error] = p.result(); 
           return std::tuple<std::vector<double>, std::vector<AverageEnergyType>, std::vector<AverageEnergyType>>
                                                                 {Units::EnergyToKelvin * bins, average, error};});

  pybind11::class_<PropertyConventionalRadialDistributionFunction>(m, "PropertyConventionalRadialDistributionFunction")
       .def(pybind11::init<std::size_t, double, std::size_t, std::optional<std::size_t>>(),
                 pybind11::arg("number_of_bins"), 
                 pybind11::arg("range"),
                 pybind11::arg("sample_every"), 
                 pybind11::arg("write_every") = std::nullopt)
       .def("result", &PropertyConventionalRadialDistributionFunction::result);

  pybind11::class_<PropertyRadialDistributionFunction>(m, "PropertyRadialDistributionFunction")
       .def(pybind11::init<std::size_t, std::size_t, std::size_t, double, std::size_t, std::size_t>(),
            pybind11::arg("number_of_blocks"), 
            pybind11::arg("number_of_pseudo_atoms"),
            pybind11::arg("number_of_bins"), 
            pybind11::arg("range"),
            pybind11::arg("sample_every"), 
            pybind11::arg("write_every"))
       .def("result", &PropertyRadialDistributionFunction::result);


  pybind11::class_<MeanSquaredDisplacementData>(m, "MeanSquaredDisplacementData")
      .def_readonly("time", &MeanSquaredDisplacementData::time)
      .def_readonly("xyz", &MeanSquaredDisplacementData::xyz)
      .def_readonly("x", &MeanSquaredDisplacementData::x)
      .def_readonly("y", &MeanSquaredDisplacementData::y)
      .def_readonly("z", &MeanSquaredDisplacementData::z)
      .def_readonly("number_of_samples", &MeanSquaredDisplacementData::numberOfSamples);

  pybind11::class_<PropertyMeanSquaredDisplacement>(m, "PropertyMeanSquaredDisplacement")
       .def(pybind11::init<std::size_t, std::optional<std::size_t>>(),
            pybind11::arg("sample_every"), 
            pybind11::arg("write_every") = std::nullopt)
       .def("result", &PropertyMeanSquaredDisplacement::result);

  pybind11::class_<VelocityAutoCorrelationFunctionData>(m, "VelocityAutoCorrelationFunctionData")
      .def_readonly("time", &VelocityAutoCorrelationFunctionData::time)
      .def_readonly("xyz", &VelocityAutoCorrelationFunctionData::xyz)
      .def_readonly("x", &VelocityAutoCorrelationFunctionData::x)
      .def_readonly("y", &VelocityAutoCorrelationFunctionData::y)
      .def_readonly("z", &VelocityAutoCorrelationFunctionData::z)
      .def_readonly("number_of_samples", &VelocityAutoCorrelationFunctionData::numberOfSamples);

  pybind11::class_<PropertyVelocityAutoCorrelationFunction>(m, "PropertyVelocityAutoCorrelationFunction")
       .def(pybind11::init<std::size_t, std::size_t, std::size_t, std::optional<std::size_t>>(),
            pybind11::arg("number_of_buffers_vacf"), 
            pybind11::arg("buffer_length_vacf"),
            pybind11::arg("sample_every"), 
            pybind11::arg("write_every") = std::nullopt)
       .def("result", &PropertyVelocityAutoCorrelationFunction::result);


  pybind11::class_<PropertyNumberOfMoleculesEvolution>(m, "PropertyNumberOfMoleculesEvolution")
      .def(pybind11::init<std::size_t, std::size_t, std::size_t, std::optional<std::size_t>>(),
            pybind11::arg("number_of_cycles"), 
            pybind11::arg("number_of_components"),
            pybind11::arg("sample_every"), 
            pybind11::arg("write_every") = std::nullopt)
      .def_readonly("result", &PropertyNumberOfMoleculesEvolution::result);

  pybind11::class_<PropertyVolumeEvolution>(m, "PropertyVolumeEvolution")
    .def(pybind11::init<std::size_t, std::size_t, std::optional<std::size_t>>(),
            pybind11::arg("number_of_cycles"), 
            pybind11::arg("sample_every"), 
            pybind11::arg("write_every") = std::nullopt)
      .def_readonly("result", &PropertyVolumeEvolution::result);



  pybind11::class_<System>(m, "System")
      .def(pybind11::init<ForceField, std::optional<SimulationBox>, bool, double, std::optional<double>, double,
                          std::optional<Framework>, std::vector<Component>, std::vector<std::vector<double3>>,
                          std::vector<std::size_t>, std::size_t, MCMoveProbabilities>(),
           pybind11::arg("force_field"), 
           pybind11::arg("simulation_box") = std::nullopt,
           pybind11::arg("has_external_field") = false, 
           pybind11::arg("external_temperature") = 298.0, 
           pybind11::arg("external_pressure") = std::nullopt, 
           pybind11::arg("helium_void_fraction") = 0.0,
           pybind11::arg("framework_components") = std::nullopt,
           pybind11::arg("components") = std::vector<Component>(), 
           pybind11::arg("initial_positions") = std::vector<std::vector<double3>>(),
           pybind11::arg("initial_number_of_molecules") = std::vector<std::size_t>(),
           pybind11::arg("number_of_blocks") = 5,
           pybind11::arg("system_probabilities") = MCMoveProbabilities())
      .def("compute_total_energies", &System::computeTotalEnergies)
      .def("framework_mass", &System::frameworkMass)
      .def_readonly("simulation_box", &System::simulationBox)
      .def_readonly("input_pressure", &System::input_pressure)
      .def_readonly("components", &System::components)
      .def_readonly("number_of_molecules_per_component", &System::numberOfMoleculesPerComponent)
      .def_readonly("mc_moves_statistics", &System::mc_moves_statistics)
      .def_readonly("loadings", &System::loadings)
      .def_readwrite("average_loadings", &System::averageLoadings)
      .def_readonly("average_energies", &System::averageEnergies)
      .def_readonly("average_pressure", &System::averagePressure)
      .def_readonly("thermostat", &System::thermostat)
      .def("set_thermostat", &System::setThermostat,
         pybind11::arg("thermostat") = std::nullopt)
      .def_readonly("property_density_grid", &System::propertyDensityGrid)
      .def("set_property_density_grid", &System::setPropertyDensityGrid,
         pybind11::arg("property:") = std::nullopt)
      .def_readonly("average_energy_histogram", &System::averageEnergyHistogram)
      .def("set_average_energy_histogram", &System::setAverageEnergyHistogram,
         pybind11::arg("property:") = std::nullopt)
      .def_readonly("sample_pdb_movie", &System::samplePDBMovie)
      .def("set_sample_pdb_movie", &System::setSamplePDBMovie,
         pybind11::arg("property:") = std::nullopt)
      .def_readonly("property_conventional_rdf", &System::propertyConventionalRadialDistributionFunction)
      .def("set_property_conventional_rdf", &System::setPropertyConventionalRDF,
         pybind11::arg("property") = std::nullopt)
      .def_readonly("property_msd", &System::propertyMSD)
      .def("set_property_msd", &System::setPropertyMSD,
         pybind11::arg("property") = std::nullopt)
      .def_readonly("property_vacf", &System::propertyVACF)
      .def("set_property_vacf", &System::setPropertyVACF,
         pybind11::arg("property") = std::nullopt)
      .def_readonly("property_rdf", &System::propertyRadialDistributionFunction)
      .def("set_property_rdf", &System::setPropertyRDF,
         pybind11::arg("property") = std::nullopt)
      .def_readonly("property_number_of_molecules_evolution", &System::propertyNumberOfMoleculesEvolution)
      .def("set_property_number_of_molecules_evolution", &System::setPropertyNumberOfMoleculesEvolution,
         pybind11::arg("property") = std::nullopt)
      .def_readonly("property_volume_evolution", &System::propertyVolumeEvolution)
      .def("set_property_volume_evolution", &System::setPropertyVolumeEvolution,
         pybind11::arg("property") = std::nullopt)
      .def_readonly("property_conserved_energy_evolution", &System::propertyConservedEnergyEvolution)
      .def("set_property_conserved_energy_evolution", &System::setPropertyConservedEnergyEvolution,
         pybind11::arg("property") = std::nullopt)
      .def_readwrite("average_enthalpies_of_adsorption", &System::averageEnthalpiesOfAdsorption)
      .def_readwrite("atom_data", &System::atomData)
      .def_readwrite("translational_degrees_of_freedom", &System::translationalDegreesOfFreedom)
      .def_readwrite("rotational_degrees_of_freedom", &System::rotationalDegreesOfFreedom)
      .def("write_mc_move_statistics", &System::writeMCMoveStatistics)
      .def("__repr__", &System::repr);


  pybind11::class_<MonteCarlo>(m, "MonteCarlo")
      .def(pybind11::init<std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t,
                        std::vector<System> &, std::optional<std::size_t>, std::size_t, bool>(),
         pybind11::arg("number_of_cycles"),
         pybind11::arg("number_of_initialization_cycles"),
         pybind11::arg("number_of_equilibration_cycles") = 0, 
         pybind11::arg("print_every") = 5000,
         pybind11::arg("write_binary_restart_every") = 5000, 
         pybind11::arg("rescale_wang_landau_every") = 1000,
         pybind11::arg("optimize_mc_moves_every") = 100, 
         pybind11::arg("systems") = std::vector<System>(),
         pybind11::arg("random_seed") = std::nullopt, 
         pybind11::arg("number_of_blocks") = 5,
         pybind11::arg("output_to_files") = false)
      .def("run", &MonteCarlo::run)
      .def("setup", &MonteCarlo::setup)
      .def("tear_down", &MonteCarlo::tearDown)
      .def("initialize", &MonteCarlo::initialize, pybind11::arg("call_back_function") = pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def("equilibrate", &MonteCarlo::equilibrate, pybind11::arg("call_back_function") = pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def("production", &MonteCarlo::production, pybind11::arg("call_back_function") = pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def("cycle", &MonteCarlo::performCycle)
      .def_readonly("systems", &MonteCarlo::systems);




  pybind11::class_<MolecularDynamics>(m, "MolecularDynamics")
    .def(pybind11::init<std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t,
                        std::vector<System> &, std::optional<std::size_t>, std::size_t, bool>(),
         pybind11::arg("number_of_cycles"), 
         pybind11::arg("number_of_initialization_cycles"),
         pybind11::arg("number_of_equilibration_cycles") = 0, 
         pybind11::arg("print_every") = 5000,
         pybind11::arg("write_binary_restart_every") = 5000, 
         pybind11::arg("rescale_wang_landau_every") = 1000,
         pybind11::arg("optimize_mc_moves_every") = 100, 
         pybind11::arg("systems") = std::vector<System>(),
         pybind11::arg("random_seed") = std::nullopt, 
         pybind11::arg("number_of_blocks") = 5,
         pybind11::arg("output_to_files") = false)
      .def("run", &MolecularDynamics::run)
      .def("setup", &MolecularDynamics::setup)
      .def("tear_down", &MolecularDynamics::tearDown)
      .def("initialize", &MolecularDynamics::initialize, pybind11::arg("call_back_function") = 
                      pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def("equilibrate", &MolecularDynamics::equilibrate, pybind11::arg("call_back_function") = 
                      pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def("production", &MolecularDynamics::production, pybind11::arg("call_back_function") = 
                      pybind11::cpp_function([](void){}), pybind11::arg("call_back_every") = 100)
      .def_readonly("systems", &MolecularDynamics::systems);

  pybind11::class_<RunningEnergy>(m, "RunningEnergy")
      .def("conserved_energy", &RunningEnergy::conservedEnergy)
      .def("potential_energy", &RunningEnergy::potentialEnergy)
      .def("kinetic_energy", &RunningEnergy::kineticEnergy)
      .def_readonly("nose_hoover_energy", &RunningEnergy::NoseHooverEnergy)
      .def_readwrite("molecule_molecule_vdw", &RunningEnergy::moleculeMoleculeVDW)
      .def_readwrite("framework_molecule_vdw", &RunningEnergy::frameworkMoleculeVDW)
      .def("__repr__", &RunningEnergy::repr);

  pybind11::class_<PropertyConservedEnergyEvolution>(m, "PropertyConservedEnergyEvolution")
      .def(pybind11::init<std::size_t, std::size_t, std::optional<std::size_t>>(),
            pybind11::arg("number_of_cycles"), 
            pybind11::arg("sample_every"), 
            pybind11::arg("write_every") = std::nullopt)
       // convert result to units of Kelvin
      .def("result", [](const PropertyConservedEnergyEvolution& p) { return Units::EnergyToKelvin * p.result();});

  pybind11::class_<Thermostat>(m, "Thermostat")
      .def(pybind11::init<std::size_t, std::size_t>(),
            pybind11::arg("thermostat_chain_length") = 5, 
            pybind11::arg("number_of_yoshida_suzuki_steps") = 5);

  pybind11::class_<EnthalpyOfAdsorptionData>(m, "EnthalpyOfAdsorptionData")
      .def(pybind11::init<std::size_t>(),
            pybind11::arg("number_of_components"))
      .def_readonly("values", &EnthalpyOfAdsorptionData::values)
      .def("__getitem__", &EnthalpyOfAdsorptionData::operator[]);

  pybind11::class_<PropertyEnthalpy>(m, "PropertyEnthalpy")
      .def(pybind11::init<std::size_t, std::size_t>(),
            pybind11::arg("number_of_blocks"), 
            pybind11::arg("number_of_components"))
      .def("result", [](const PropertyEnthalpy& p) { return Units::EnergyToKelvin * p.result();});

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
