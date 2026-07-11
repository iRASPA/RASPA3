#include <nanobind/stl/chrono.h>
#include <nanobind/stl/complex.h>
#include <nanobind/ndarray.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/tuple.h>

import std;

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
import simulation_schedule;
import input_reader;
import property_loading;
import loading_data;
import energy_factor;
import energy_status;
import sample_movies;
import property_energy;
import average_energy_type;
import property_energy_histogram;
import property_number_of_molecules_histogram;
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
import partial_molar_properties_data;
import property_partial_molar_properties;
import property_conventional_rdf;
import property_rdf;
import mean_squared_displacement_data;
import property_msd;
import velocity_autocorrelation_function_data;
import property_vacf;

template <typename T>
std::pair<T, T>& operator*=(std::pair<T, T>& p, const double& s)
{
    p.first *= s;
    p.second *= s;
    return p;
}

template <typename T>
std::pair<T, T>& operator*=(const double& s, std::pair<T, T>& p)
{
    p.first *= s;
    p.second *= s;
    return p;
}

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

namespace nanobind {
namespace detail {

template <>
struct type_caster<int3> 
{
  // This macro inserts a lot of boilerplate code and sets the type hint.
  // `io_name` is used to specify different type hints for arguments and return values.
  // The signature of our negate function would then look like:
  // `negate(Sequence[float]) -> tuple[float, float]`
  NB_TYPE_CASTER(int3, io_name("Sequence[int]", "tuple[int, int, int]"));

  // C++ -> Python: convert `Point2D` to `tuple[float, float]`. The second and third arguments
  // are used to indicate the return value policy and parent object (for
  // return_value_policy::reference_internal) and are often ignored by custom casters.
  // The return value should reflect the type hint specified by the second argument of `io_name`.
  static handle from_cpp(const int3 &number, rv_policy /* policy */, cleanup_list * /* cleanup */) noexcept
  {
      return nanobind::make_tuple(number.x, number.y, number.z).release();
  }

  // Python -> C++: convert a `PyObject` into a `int3` and return false upon failure. The
  // second argument indicates whether implicit conversions should be allowed.
  // The accepted types should reflect the type hint specified by the first argument of
  // `io_name`.
  bool from_python(handle src, uint8_t flags, cleanup_list *cleanup) noexcept
  {
    // Check if handle is a Sequence
    if (!nanobind::isinstance<nanobind::sequence>(src)) 
    {
      return false;
    }
    auto seq = nanobind::borrow<nanobind::sequence>(src);
    // Check if exactly two values are in the Sequence
    if (nanobind::len(seq) != 3) 
    {
      return false;
    }
    // Check if each element is either an int
    for (auto item : seq)
    {
      if (!nanobind::isinstance<nanobind::int_>(item))
      {
        return false;
      }
    }
    value.x = nanobind::cast<std::int32_t>(seq[0]);
    value.y = nanobind::cast<std::int32_t>(seq[1]);
    value.z = nanobind::cast<std::int32_t>(seq[2]);
    return true;
  }
};


template <>
struct type_caster<double3>
{
public:
  // Sets positional type hints for Python tooling
  NB_TYPE_CASTER(double3, io_name("Sequence[float] | float", "tuple[float, float, float]"));

  static handle from_cpp(const double3 &number, rv_policy /* policy */, cleanup_list * /* cleanup */) noexcept
  {
    return nanobind::make_tuple(number.x, number.y, number.z).release();
  }

  // Python -> C++
  bool from_python(handle src, uint8_t flags, cleanup_list *cleanup) noexcept
  {
    // Case 1: Single scalar value passed directly (e.g., 5.0 -> [5.0, 5.0, 5.0])
    if (nanobind::isinstance<nanobind::float_>(src) || nanobind::isinstance<nanobind::int_>(src))
    {
      // Fixed: Functional-style nanobind::cast
      double scalar = nanobind::cast<double>(src);
      value.x = scalar;
      value.y = scalar;
      value.z = scalar;
      return true;
    }

    // Case 2: Sequence container passed down
    if (nanobind::isinstance<nanobind::sequence>(src))
    {
      auto seq = nanobind::borrow<nanobind::sequence>(src);

      // Sub-case A: Single-element sequence (e.g., [5.0] -> [5.0, 5.0, 5.0])
      if (nanobind::len(seq) == 1)
      {
        auto item = seq[0];
        if (!nanobind::isinstance<nanobind::float_>(item) && !nanobind::isinstance<nanobind::int_>(item))
        {
          return false;
        }

        double scalar = nanobind::cast<double>(item);
        value.x = scalar;
        value.y = scalar;
        value.z = scalar;
        return true;
      }

      // Sub-case B: Full three-element sequence (e.g., [1.0, 2.0, 3.0])
      if (nanobind::len(seq) == 3)
      {
        for (auto item : seq)
        {
          if (!nanobind::isinstance<nanobind::float_>(item) && !nanobind::isinstance<nanobind::int_>(item))
          {
            return false;
          }
        }

        // Fixed: Functional-style nanobind::cast on indexed sequence elements
        value.x = nanobind::cast<double>(seq[0]);
        value.y = nanobind::cast<double>(seq[1]);
        value.z = nanobind::cast<double>(seq[2]);
        return true;
      }
    }

    return false;
  }
};


template <>
struct type_caster<double4>
{
public:
  // Sets positional type hints for Python tooling
  NB_TYPE_CASTER(double4, io_name("Sequence[float] | float", "tuple[float, float, float, float]"));

  static handle from_cpp(const double4 &number, rv_policy /* policy */, cleanup_list * /* cleanup */) noexcept
  {
    return nanobind::make_tuple(number.x, number.y, number.z, number.w).release();
  }

  // Python -> C++
  bool from_python(handle src, uint8_t flags, cleanup_list *cleanup) noexcept
  {
    // Case 1: Single scalar value passed directly (e.g., 5.0 -> [5.0, 5.0, 5.0])
    if (nanobind::isinstance<nanobind::float_>(src) || nanobind::isinstance<nanobind::int_>(src))
    {
      // Fixed: Functional-style nanobind::cast
      double scalar = nanobind::cast<double>(src);
      value.x = scalar;
      value.y = scalar;
      value.z = scalar;
      value.w = scalar;
      return true;
    }

    // Case 2: Sequence container passed down
    if (nanobind::isinstance<nanobind::sequence>(src))
    {
      auto seq = nanobind::borrow<nanobind::sequence>(src);

      // Sub-case A: Single-element sequence (e.g., [5.0] -> [5.0, 5.0, 5.0])
      if (nanobind::len(seq) == 1)
      {
        auto item = seq[0];
        if (!nanobind::isinstance<nanobind::float_>(item) && !nanobind::isinstance<nanobind::int_>(item))
        {
          return false;
        }

        double scalar = nanobind::cast<double>(item);
        value.x = scalar;
        value.y = scalar;
        value.z = scalar;
        value.w = scalar;
        return true;
      }

      // Sub-case B: Full three-element sequence (e.g., [1.0, 2.0, 3.0])
      if (nanobind::len(seq) == 4)
      {
        for (auto item : seq)
        {
          if (!nanobind::isinstance<nanobind::float_>(item) && !nanobind::isinstance<nanobind::int_>(item))
          {
            return false;
          }
        }

        // Fixed: Functional-style nanobind::cast on indexed sequence elements
        value.x = nanobind::cast<double>(seq[0]);
        value.y = nanobind::cast<double>(seq[1]);
        value.z = nanobind::cast<double>(seq[2]);
        value.w = nanobind::cast<double>(seq[3]);
        return true;
      }
    }

    return false;
  }
};


} // namespace detail
} // namespace nanobind


#define EXPAND_MODULE(name) NB_MODULE(name, m)
EXPAND_MODULE(MODULE_NAME)
{
  nanobind::class_<EquationOfState> eos(m, "EquationOfState");

    nanobind::class_<EquationOfState::FluidInput>(eos, "FluidInput")
        .def(nanobind::init<double, double, double, double, bool>(), 
            nanobind::arg("critical_temperature"), 
            nanobind::arg("critical_pressure"),
            nanobind::arg("acentric_factor"), 
            nanobind::arg("mol_fraction") = 1.0,
            nanobind::arg("swappable") = true)
        .def_rw("critical_temperature", &EquationOfState::FluidInput::criticalTemperature)
        .def_rw("critical_pressure", &EquationOfState::FluidInput::criticalPressure)
        .def_rw("acentric_factor", &EquationOfState::FluidInput::acentricFactor)
        .def_rw("mol_fraction", &EquationOfState::FluidInput::molFraction)
        .def_rw("swappable", &EquationOfState::FluidInput::swappable);

    nanobind::enum_<EquationOfState::Type>(eos, "EquationOfStateType", nanobind::is_arithmetic(), nanobind::is_flag())
        .value("PENG_ROBINSON", EquationOfState::Type::PengRobinson)
        .value("PENG_ROBINSON_GASEM", EquationOfState::Type::PengRobinsonGasem)
        .value("SOAVE_REDLICH_KWONG", EquationOfState::Type::SoaveRedlichKwong);

    nanobind::enum_<EquationOfState::MixingRules>(eos, "MixingRules", nanobind::is_arithmetic(), nanobind::is_flag())
        .value("VAN_DER_WAALS", EquationOfState::MixingRules::VanDerWaals);

    nanobind::enum_<EquationOfState::FluidState>(eos, "FluidState", nanobind::is_arithmetic(), nanobind::is_flag())
        .value("UNKNOWN", EquationOfState::FluidState::Unknown)
        .value("SUPER_CRITICAL_FLUID", EquationOfState::FluidState::SuperCriticalFluid)
        .value("VAPOR", EquationOfState::FluidState::Vapor)
        .value("LIQUID", EquationOfState::FluidState::Liquid)
        .value("VAPOR_LIQUID", EquationOfState::FluidState::VaporLiquid);

  nanobind::class_<EquationOfState::FluidResult>(eos, "FluidResult")
      .def(nanobind::init<double, std::optional<double>, EquationOfState::FluidState>(), 
          nanobind::arg("compressibility"), 
          nanobind::arg("fugacity_coefficient"),
          nanobind::arg("fluid_state"))
      .def_rw("compressibility", &EquationOfState::FluidResult::compressibility)
      .def_rw("fugacity_coefficient", &EquationOfState::FluidResult::fugacityCoefficient)
      .def_rw("fluid_state", &EquationOfState::FluidResult::fluidState);

    eos.def(nanobind::init<EquationOfState::Type, EquationOfState::MixingRules,
                  double, double, const SimulationBox &, double,
                  std::vector<Component> &>())
       .def_static("compute_fluid_properties", &EquationOfState::computeFluidProperties,
         nanobind::arg("temperature"), 
         nanobind::arg("pressure"),
         nanobind::arg("properties"),
         nanobind::arg("type") = EquationOfState::Type::PengRobinson,
         nanobind::arg("mixing_rules") = EquationOfState::MixingRules::VanDerWaals);

  nanobind::class_<RandomNumber>(m, "RandomNumber")
    .def(nanobind::init<std::size_t>(), 
         nanobind::arg("seed") = 12);

  nanobind::class_<Atom>(m, "Atom")
      .def(nanobind::init<>())
      .def(nanobind::init<double3, double, double, std::uint32_t, std::uint16_t, std::uint8_t, std::uint8_t, std::uint8_t>(),
           nanobind::arg("position"), 
           nanobind::arg("charge") = 0.0,
           nanobind::arg("scaling") = 1.0,
           nanobind::arg("molecule_id") = 0, 
           nanobind::arg("type") = 0, 
           nanobind::arg("component_id") = 0,
           nanobind::arg("group_id") = false, 
           nanobind::arg("is_fractional") = false)
      .def_prop_rw("position",
          [](const Atom &a) -> double3 { return a.position; },
          [](Atom &a, const double3 &value) { a.position = value; })
      .def("__repr__", &Atom::repr);

  nanobind::class_<SimulationBox> simulationBox(m, "SimulationBox");
  nanobind::enum_<SimulationBox::Type>(simulationBox, "SimulationBoxType", nanobind::is_arithmetic(), nanobind::is_flag())
      .value("RECTANGULAR", SimulationBox::Type::Rectangular)
      .value("TRICLINIC", SimulationBox::Type::Triclinic);

  simulationBox
      .def(nanobind::init<double, double, double>(), 
           nanobind::arg("a"), 
           nanobind::arg("b"), 
           nanobind::arg("c"))
      .def(nanobind::init<double, double, double, double, double, double>(), 
           nanobind::arg("a"), 
           nanobind::arg("b"),
           nanobind::arg("c"), 
           nanobind::arg("alpha"), 
           nanobind::arg("beta"), 
           nanobind::arg("gamma"))
      .def(nanobind::init<double, double, double, double, double, double, SimulationBox::Type>(), 
           nanobind::arg("a"),
           nanobind::arg("b"), 
           nanobind::arg("c"), 
           nanobind::arg("alpha"), 
           nanobind::arg("beta"),
           nanobind::arg("gamma"), 
           nanobind::arg("type"))
      .def(nanobind::init<double3x3, SimulationBox::Type>())
      .def_rw("type", &SimulationBox::type)
      .def_ro("length_a", &SimulationBox::lengthA)
      .def_ro("length_b", &SimulationBox::lengthB)
      .def_ro("length_c", &SimulationBox::lengthC)
      .def_ro("volume", &SimulationBox::volume);

  nanobind::class_<VDWParameters>(m, "VDWParameters")
      .def(nanobind::init<double, double>(), nanobind::arg("epsilon"), nanobind::arg("sigma"));

  nanobind::class_<PseudoAtom>(m, "PseudoAtom")
      .def(nanobind::init<std::string, bool, double, double, double, std::size_t, bool, std::string>(),
           nanobind::arg("name"), 
           nanobind::arg("framework_type"), 
           nanobind::arg("mass"), 
           nanobind::arg("charge"),
           nanobind::arg("polarizability") = 0.0, 
           nanobind::arg("atomic_number") = 6, 
           nanobind::arg("print_to_pdb") = true,
           nanobind::arg("source") = "");

  nanobind::class_<ForceField> forceField(m, "ForceField");
  nanobind::enum_<ForceField::MixingRule>(forceField, "MixingRule", nanobind::is_arithmetic(), nanobind::is_flag())
      .value("LORENTZ_BERTHELOT", ForceField::MixingRule::Lorentz_Berthelot)
      .value("JORGENSEN", ForceField::MixingRule::Jorgensen);

  forceField.def(nanobind::init<std::vector<PseudoAtom>, std::vector<VDWParameters>, ForceField::MixingRule, 
                 double, double, double, bool, bool, bool>(),
           nanobind::arg("pseudo_atoms"), 
           nanobind::arg("parameters"), 
           nanobind::arg("mixing_rule") = ForceField::MixingRule::Lorentz_Berthelot,
           nanobind::arg("cutoff_framework_vdw") = 12.0, 
           nanobind::arg("cutoff_molecule_vdw") = 12.0,
           nanobind::arg("cutoff_coulomb") = 12.0, 
           nanobind::arg("shifted") = true,
           nanobind::arg("tail_corrections") = false, 
           nanobind::arg("use_charge") = true)
      .def("find_pseudo_atom", static_cast<std::optional<std::size_t> (ForceField::*)(const std::string &) const>(&ForceField::findPseudoAtom))
      .def("__repr__", &ForceField::repr)
      .def_ro("pseudo_atoms", &ForceField::pseudoAtoms)
      .def_ro("vdw_parameters", &ForceField::data)
      .def_rw("use_charge", &ForceField::useCharge);

  nanobind::class_<CIFReader> cifReader(m, "CIFReader");
  cifReader
      .def("expand_defined_atoms_to_unit_cell", &CIFReader::expandDefinedAtomsToUnitCell);

  nanobind::class_<Framework> framework(m, "Framework");

  nanobind::enum_<CIFReader::UseChargesFrom>(framework, "UseChargesFrom", nanobind::is_arithmetic(), nanobind::is_flag())
      .value("PSEUDO_ATOMS", CIFReader::UseChargesFrom::PseudoAtoms)
      .value("CIF_FILE", CIFReader::UseChargesFrom::CIF_File)
      .value("CHARGE_EQUILIBRATION", CIFReader::UseChargesFrom::ChargeEquilibration);

  framework
      .def(nanobind::init<const ForceField &, std::string, SimulationBox, std::size_t, 
                          const std::vector<Atom> &, int3>(),
           nanobind::arg("force_field"), 
           nanobind::arg("component_name"),
           nanobind::arg("simulation_box"), 
           nanobind::arg("space_group_hall_number"), 
           nanobind::arg("defined_atoms"), 
           nanobind::arg("number_of_unit_cells"))
      .def_ro("name", &Framework::name)
      .def("print", &Framework::printStatus)
      .def("__repr__", &Framework::repr);

  nanobind::class_<MCMoveProbabilities>(m, "MCMoveProbabilities")
      .def(nanobind::init<double, double, double, double, double, double, double, double, double, double, double,
                          double, double, double, double, double, double, double, double, double>(),
           nanobind::arg("translation_probability") = 0.0, 
           nanobind::arg("random_translation_probability") = 0.0,
           nanobind::arg("rotation_probability") = 0.0, 
           nanobind::arg("random_rotation_probability") = 0.0,
           nanobind::arg("volume_change_probability") = 0.0, 
           nanobind::arg("reinsertion_cbmc_probability") = 0.0,
           nanobind::arg("partial_reinsertion_cbmc_probability") = 0.0,
           nanobind::arg("identity_change_probability") = 0.0, 
           nanobind::arg("swap_probability") = 0.0,
           nanobind::arg("swap_cbmc_probability") = 0.0, 
           nanobind::arg("swap_cfcmc_probability") = 0.0,
           nanobind::arg("swap_cbcfcmc_probability") = 0.0, 
           nanobind::arg("gibbs_volume_change_probability") = 0.0,
           nanobind::arg("gibbs_swap_cbmc_probability") = 0.0, 
           nanobind::arg("gibbs_swap_cfcmc_probability") = 0.0,
           nanobind::arg("widom_probability") = 0.0, 
           nanobind::arg("widom_cfcmc_probability") = 0.0,
           nanobind::arg("widom_cbcfcmc_probability") = 0.0, 
           nanobind::arg("parallel_tempering_probability") = 0.0,
           nanobind::arg("hybrid_mc_probability") = 0.0)
      .def("join", &MCMoveProbabilities::join)
      .def("__repr__", &MCMoveProbabilities::repr);


  nanobind::class_<PropertyLambdaProbabilityHistogram>(m, "PropertyLambdaProbabilityHistogram")
      .def(nanobind::init<>())
      // dimensionless
      .def_rw("bias_factor", &PropertyLambdaProbabilityHistogram::biasFactor)
      // dimensionless
      .def_rw("histogram", &PropertyLambdaProbabilityHistogram::histogram)
      .def_ro("occupancy_total", &PropertyLambdaProbabilityHistogram::occupancyTotal)
      .def_ro("occupancy_count", &PropertyLambdaProbabilityHistogram::occupancyCount)
       // convert result to units of Kelvin
      .def("average_dudlambda", [](PropertyLambdaProbabilityHistogram& p) 
          { return Units::EnergyToKelvin * p.averageDuDlambda();})
      // dimensionless
      .def("result", &PropertyLambdaProbabilityHistogram::result);


  nanobind::class_<ConnectivityTable>(m, "ConnectivityTable")
    .def(nanobind::init<>());

  nanobind::class_<Potentials::IntraMolecularPotentials>(m, "IntraMolecularPotentials")
    .def(nanobind::init<>());

  nanobind::class_<MoveStatistics<double>>(m, "MoveStatisticsDouble")
    .def_ro("counts", &MoveStatistics<double>::counts)
    .def_ro("constructed", &MoveStatistics<double>::constructed)
    .def_ro("accepted", &MoveStatistics<double>::accepted)
    .def_ro("all_counts", &MoveStatistics<double>::allCounts)
    .def_ro("total_counts", &MoveStatistics<double>::totalCounts)
    .def_ro("total_constructed", &MoveStatistics<double>::totalConstructed)
    .def_ro("total_accepted", &MoveStatistics<double>::totalAccepted)
    .def_rw("max_change", &MoveStatistics<double>::maxChange)
    .def_rw("target_acceptance", &MoveStatistics<double>::targetAcceptance)
    .def_rw("lower_limit", &MoveStatistics<double>::lowerLimit)
    .def_rw("upper_limit", &MoveStatistics<double>::upperLimit)
    .def_rw("optimize", &MoveStatistics<double>::optimize);

  nanobind::class_<MoveStatistics<double3>>(m, "MoveStatisticsDouble3")
    .def_ro("counts", &MoveStatistics<double3>::counts)
    .def_ro("constructed", &MoveStatistics<double3>::constructed)
    .def_ro("accepted", &MoveStatistics<double3>::accepted)
    .def_ro("all_counts", &MoveStatistics<double3>::allCounts)
    .def_ro("total_counts", &MoveStatistics<double3>::totalCounts)
    .def_ro("total_constructed", &MoveStatistics<double3>::totalConstructed)
    .def_ro("total_accepted", &MoveStatistics<double3>::totalAccepted)
    .def_rw("max_change", &MoveStatistics<double3>::maxChange)
    .def_rw("target_acceptance", &MoveStatistics<double3>::targetAcceptance)
    .def_rw("lower_limit", &MoveStatistics<double3>::lowerLimit)
    .def_rw("upper_limit", &MoveStatistics<double3>::upperLimit)
    .def_rw("optimize", &MoveStatistics<double3>::optimize);

  nanobind::class_<Move> move(m, "Move");
    move.def(nanobind::init());
  nanobind::enum_<Move::Types>(move, "Types", nanobind::is_arithmetic(), nanobind::is_flag())
      .value("TRANSLATION", Move::Types::Translation)
      .value("RANDOM_TRANSLATION", Move::Types::RandomTranslation)
      .value("ROTATION", Move::Types::Rotation)
      .value("RANDOM_ROTATION", Move::Types::RandomRotation)
      .value("VOLUME_CHANGE", Move::Types::VolumeChange)
      .value("ANISOTROPIC_VOLUME_CHANGE", Move::Types::AnisotropicVolumeChange)
      .value("REINSERTION_CBMC", Move::Types::ReinsertionCBMC)
      .value("PARTIAL_REINSERTION_CBMC", Move::Types::PartialReinsertionCBMC)
      .value("IDENTITY_CHANGE_CBMC", Move::Types::IdentityChangeCBMC)
      .value("SWAP", Move::Types::Swap)
      .value("SWAP_CBMC", Move::Types::SwapCBMC)
      .value("SWAP_CFCMC", Move::Types::SwapCFCMC)
      .value("SWAP_CBCFCMC", Move::Types::SwapCBCFCMC)
      .value("GIBBS_VOLUME", Move::Types::GibbsVolume)
      .value("GIBBS_SWAP_CBMC", Move::Types::GibbsSwapCBMC)
      .value("GIBBS_CONVENTIONAL_CFCMC", Move::Types::GibbsConventionalCFCMC)
      .value("GIBBS_CONVENTIONAL_CBCFCMC", Move::Types::GibbsConventionalCBCFCMC)
      .value("GIBBS_SWAP_CFCMC", Move::Types::GibbsSwapCFCMC)
      .value("GIBBS_SWAP_CBCFCMC", Move::Types::GibbsSwapCBCFCMC)
      .value("GIBBS_IDENTITY_CHANGE_CBMC", Move::Types::GibbsIdentityChangeCBMC)
      .value("WIDOM", Move::Types::Widom)
      .value("WIDOM_CFCMC", Move::Types::WidomCFCMC)
      .value("WIDOM_CBCFCMC", Move::Types::WidomCBCFCMC)
      .value("PARALLEL_TEMPERING", Move::Types::ParallelTempering)
      .value("HYBRID_MC", Move::Types::HybridMC)
      .value("REACTION_CBMC", Move::Types::ReactionCBMC)
      .value("REACTION_CONVENTIONAL_CFCMC", Move::Types::ReactionConventionalCFCMC)
      .value("REACTION_CONVENTIONAL_CBCFCMC", Move::Types::ReactionConventionalCBCFCMC)
      .value("REACTION_CFCMC", Move::Types::ReactionCFCMC)
      .value("REACTION_CBCFCMC", Move::Types::ReactionCBCFCMC)
      .value("PAIR_SWAP", Move::Types::PairSwap)
      .value("PAIR_SWAP_CBMC", Move::Types::PairSwapCBMC)
      .value("PAIR_SWAP_CFCMC", Move::Types::PairSwapCFCMC)
      .value("PAIR_SWAP_CBCFCMC", Move::Types::PairSwapCBCFCMC);

  nanobind::class_<MCMoveStatistics>(m, "MCMoveStatistics")
      .def("__getitem__", [](MCMoveStatistics &self, Move::Types i) { return self[i]; })
      .def("__repr__", &MCMoveStatistics::repr);

  nanobind::class_<WidomData>(m, "WidomData")
      .def_ro("total", &WidomData::total)
      .def_ro("excess", &WidomData::excess)
      .def_ro("ideal_gas", &WidomData::idealGas);

  nanobind::class_<PropertyWidom>(m, "PropertyWidom")
      .def("result", &PropertyWidom::result)
       // convert result to units of Kelvin
      .def("chemical_potential_result", [](PropertyWidom& p, double T) 
          { return Units::EnergyToKelvin * p.chemicalPotentialResult(1.0 / (Units::KB * T));}, nanobind::arg("temperature"))
       // convert result to units of Pascal
      .def("fugacity_result", [](PropertyWidom& p, double T) 
          { return Units::PressureConversionFactor * p.fugacityResult(1.0 / (Units::KB * T));}, nanobind::arg("temperature"));

  nanobind::class_<PropertyGibbsWidom>(m, "PropertyGibbsWidom")
      .def("result", &PropertyGibbsWidom::result)
       // convert result to units of Kelvin
      .def("chemical_potential_result", [](PropertyGibbsWidom& p, double T) 
          { return Units::EnergyToKelvin * p.chemicalPotentialResult(1.0 / (Units::KB * T));}, nanobind::arg("temperature"))
       // convert result to units of Pascal
      .def("fugacity_result", [](PropertyGibbsWidom& p, double T) 
          { return Units::PressureConversionFactor * p.fugacityResult(1.0 / (Units::KB * T));}, nanobind::arg("temperature"));

  nanobind::class_<Component> component(m, "Component");
  nanobind::enum_<Component::Type>(component, "Type", nanobind::is_arithmetic(), nanobind::is_flag())
      .value("ADSORBATE", Component::Type::Adsorbate)
      .value("CATION", Component::Type::Cation);

  component
      .def(nanobind::init<const ForceField &, std::string, double, double, double, std::vector<Atom>,
                          const ConnectivityTable &, const Potentials::IntraMolecularPotentials &, std::size_t,
                          std::size_t, const MCMoveProbabilities &, std::optional<double>, bool, std::vector<double4>>(),
           nanobind::arg("force_field"), 
           nanobind::arg("component_name"),
           nanobind::arg("critical_temperature"), 
           nanobind::arg("critical_pressure"), 
           nanobind::arg("acentric_factor"),
           nanobind::arg("defined_atoms") = std::vector<Atom>(), 
           nanobind::arg("connectivity_table") = ConnectivityTable(),
           nanobind::arg("intra_molecular_potentials") = Potentials::IntraMolecularPotentials(),
           nanobind::arg("number_of_blocks") = 5,
           nanobind::arg("number_of_lambda_bins") = 41,
           nanobind::arg("particle_probabilities") = MCMoveProbabilities(), 
           nanobind::arg("fugacity_coefficient") = nanobind::none(),
           nanobind::arg("thermodynamic_integration") = false,
           nanobind::arg("blocking_pockets") = std::vector<double4>())
      .def_ro("name", &Component::name)
      .def_rw("lambda_histogram", &Component::lambdaGC)
      .def_ro("mc_moves_probabilities", &Component::mc_moves_probabilities)
      .def_ro("mc_moves_statistics", &Component::mc_moves_statistics)
      .def_rw("blocking_pockets", &Component::blockingPockets)
      .def_rw("average_rosenbluth_weights", &Component::averageRosenbluthWeights)
      .def_rw("average_gibbs_rosenbluth_weights", &Component::averageGibbsRosenbluthWeights)
      .def("print_status", &Component::printStatus)
      .def("__repr__", &Component::repr);

  nanobind::class_<LoadingData>(m, "LoadingData")
      .def(nanobind::init<std::size_t>())
      .def_ro("number_of_molecules", &LoadingData::numberOfMolecules)
      .def_ro("number_densities", &LoadingData::numberDensities)
      .def_ro("inverse_number_densities", &LoadingData::inverseNumberDensities);

  nanobind::class_<SampleMovie>(m, "SampleMovie")
      .def(nanobind::init<std::size_t, std::size_t, bool, std::optional<std::string>>(),
           nanobind::arg("system_id"), 
           nanobind::arg("sample_every"), 
           nanobind::arg("restrict_to_box") = true,
           nanobind::arg("tag") = std::nullopt);

  nanobind::class_<PropertyLoading>(m, "PropertyLoading")
      .def(nanobind::init<std::size_t, std::size_t>())
      .def("result", &PropertyLoading::result)
      .def("average_loading_number_of_molecules", &PropertyLoading::averageLoadingNumberOfMolecules)
      .def("write_averages_statistics", &PropertyLoading::writeAveragesStatistics)
      .def("__repr__", &PropertyLoading::repr);

  nanobind::class_<Potentials::EnergyFactor>(m, "EnergyFactor")
      .def_ro("energy", &Potentials::EnergyFactor::energy)
      .def_ro("dudlambda", &Potentials::EnergyFactor::dUdlambda);

  nanobind::class_<EnergyStatus>(m, "EnergyStatus")
      .def_rw("total_energy", &EnergyStatus::totalEnergy)
      .def("__repr__", &EnergyStatus::repr);

  nanobind::class_<PropertyEnergy>(m, "PropertyEnergy")
       // convert result to units of Kelvin
      .def("result", [](PropertyEnergy& p) { return Units::EnergyToKelvin * p.result();})
      .def("__repr__", &PropertyEnergy::repr);

  nanobind::class_<PressureData>(m, "PressureData")
      .def(nanobind::init<>())
      .def_ro("total_pressure", &PressureData::totalPressure)
      .def_ro("excess_pressure", &PressureData::excessPressure)
      .def_ro("ideal_gas_pressure", &PressureData::idealGasPressure);

  nanobind::class_<PropertyPressure>(m, "PropertyPressure")
      .def("result", [](PropertyPressure& p) { return Units::PressureConversionFactor * p.result();});
      //.def("__repr__", &Pressures::repr);

  nanobind::class_<AverageEnergyType>(m, "AverageEnergyType")
      .def(nanobind::init<double, double, double, double>(), 
           nanobind::arg("total_energy") = 0.0, 
           nanobind::arg("van_der_waals_energy") = 0.0,
           nanobind::arg("coulomb_energy") = 0.0, 
           nanobind::arg("polarization_energy") = 0.0)
      .def_rw("total_energy", &AverageEnergyType::totalEnergy)
      .def_rw("van_der_waals_energy", &AverageEnergyType::VanDerWaalsEnergy)
      .def_rw("coulomb_energy", &AverageEnergyType::CoulombEnergy)
      .def_rw("polarization_energy", &AverageEnergyType::polarizationEnergy);


  nanobind::class_<PropertyDensityGrid> property_energy_grid(m, "PropertyDensityGrid");

  nanobind::enum_<PropertyDensityGrid::Normalization>(property_energy_grid, "Normalization", nanobind::is_arithmetic(), nanobind::is_flag())
      .value("MAX", PropertyDensityGrid::Normalization::Max)
      .value("NUMBER_DENSITY", PropertyDensityGrid::Normalization::NumberDensity);

  nanobind::enum_<PropertyDensityGrid::Binning>(property_energy_grid, "Binning", nanobind::is_arithmetic(), nanobind::is_flag())
      .value("STANDARD", PropertyDensityGrid::Binning::Standard)
      .value("EQUITABLE", PropertyDensityGrid::Binning::Equitable);

  property_energy_grid.def(nanobind::init<std::size_t, std::size_t, int3, std::size_t, 
                                          std::size_t, std::vector<std::size_t>,
                                          PropertyDensityGrid::Normalization, PropertyDensityGrid::Binning>(),
           nanobind::arg("number_of_frameworks"), 
           nanobind::arg("number_of_components"),
           nanobind::arg("number_of_grid_points"), 
           nanobind::arg("sample_every"), 
           nanobind::arg("write_every"), 
           nanobind::arg("density_grid_pseudo_atoms_list"),
           nanobind::arg("normalization_type") = PropertyDensityGrid::Normalization::Max,
           nanobind::arg("binningMode") = PropertyDensityGrid::Binning::Standard);


  nanobind::class_<PropertyEnergyHistogram> energy_histogram(m, "PropertyEnergyHistogram");
    energy_histogram
      .def(nanobind::init<std::size_t, std::size_t, std::pair<double, double>, std::size_t, std::optional<std::size_t>>(),
                 nanobind::arg("number_of_blocks"), 
                 nanobind::arg("number_of_bins"), 
                 nanobind::arg("value_range"),
                 nanobind::arg("sample_every"), 
                 nanobind::arg("write_every") = nanobind::none())
      // convert result to units of Kelvin
      .def("result", [](PropertyEnergyHistogram& p) { auto [bins, average, error] = p.result(); 
           return std::tuple<std::vector<double>, std::vector<AverageEnergyType>, std::vector<AverageEnergyType>>
                                                                 {Units::EnergyToKelvin * bins, average, error};});

  nanobind::class_<PropertyNumberOfMoleculesHistogram>(m, "PropertyNumberOfMoleculesHistogram")
      .def(nanobind::init<std::size_t, std::pair<double, double>, std::size_t, std::optional<std::size_t>>(),
                 nanobind::arg("number_of_blocks"), 
                 nanobind::arg("value_range"),
                 nanobind::arg("sample_every"), 
                 nanobind::arg("write_every") = nanobind::none())
       .def("result", &PropertyNumberOfMoleculesHistogram::result);

  nanobind::class_<PropertyConventionalRadialDistributionFunction>(m, "PropertyConventionalRadialDistributionFunction")
       .def(nanobind::init<std::size_t, double, std::size_t, std::optional<std::size_t>>(),
                 nanobind::arg("number_of_bins"), 
                 nanobind::arg("range"),
                 nanobind::arg("sample_every"), 
                 nanobind::arg("write_every") = nanobind::none())
       .def("result", &PropertyConventionalRadialDistributionFunction::result);

  nanobind::class_<PropertyRadialDistributionFunction>(m, "PropertyRadialDistributionFunction")
       .def(nanobind::init<std::size_t, std::size_t, std::size_t, double, std::size_t, std::size_t>(),
            nanobind::arg("number_of_blocks"), 
            nanobind::arg("number_of_pseudo_atoms"),
            nanobind::arg("number_of_bins"), 
            nanobind::arg("range"),
            nanobind::arg("sample_every"), 
            nanobind::arg("write_every"))
       .def("result", &PropertyRadialDistributionFunction::result);


  nanobind::class_<MeanSquaredDisplacementData>(m, "MeanSquaredDisplacementData")
      .def_ro("time", &MeanSquaredDisplacementData::time)
      .def_ro("xyz", &MeanSquaredDisplacementData::xyz)
      .def_ro("x", &MeanSquaredDisplacementData::x)
      .def_ro("y", &MeanSquaredDisplacementData::y)
      .def_ro("z", &MeanSquaredDisplacementData::z)
      .def_ro("number_of_samples", &MeanSquaredDisplacementData::numberOfSamples);

  nanobind::class_<PropertyMeanSquaredDisplacement>(m, "PropertyMeanSquaredDisplacement")
       .def(nanobind::init<std::size_t, std::optional<std::size_t>>(),
            nanobind::arg("sample_every"), 
            nanobind::arg("write_every") = nanobind::none())
       .def("result", &PropertyMeanSquaredDisplacement::result);

  nanobind::class_<VelocityAutoCorrelationFunctionData>(m, "VelocityAutoCorrelationFunctionData")
      .def_ro("time", &VelocityAutoCorrelationFunctionData::time)
      .def_ro("xyz", &VelocityAutoCorrelationFunctionData::xyz)
      .def_ro("x", &VelocityAutoCorrelationFunctionData::x)
      .def_ro("y", &VelocityAutoCorrelationFunctionData::y)
      .def_ro("z", &VelocityAutoCorrelationFunctionData::z)
      .def_ro("number_of_samples", &VelocityAutoCorrelationFunctionData::numberOfSamples);

  nanobind::class_<PropertyVelocityAutoCorrelationFunction>(m, "PropertyVelocityAutoCorrelationFunction")
       .def(nanobind::init<std::size_t, std::size_t, std::size_t, std::optional<std::size_t>>(),
            nanobind::arg("number_of_buffers_vacf"), 
            nanobind::arg("buffer_length_vacf"),
            nanobind::arg("sample_every"), 
            nanobind::arg("write_every") = nanobind::none())
       .def("result", &PropertyVelocityAutoCorrelationFunction::result);


  nanobind::class_<PropertyNumberOfMoleculesEvolution>(m, "PropertyNumberOfMoleculesEvolution")
      .def(nanobind::init<std::size_t, std::size_t, std::size_t, std::optional<std::size_t>>(),
            nanobind::arg("number_of_production_cycles"), 
            nanobind::arg("number_of_components"),
            nanobind::arg("sample_every"), 
            nanobind::arg("write_every") = nanobind::none())
      .def_ro("result", &PropertyNumberOfMoleculesEvolution::result);

  nanobind::class_<PropertyVolumeEvolution>(m, "PropertyVolumeEvolution")
    .def(nanobind::init<std::size_t, std::size_t, std::optional<std::size_t>>(),
            nanobind::arg("number_of_production_cycles"), 
            nanobind::arg("sample_every"), 
            nanobind::arg("write_every") = nanobind::none())
      .def_ro("result", &PropertyVolumeEvolution::result);

  nanobind::class_<RunningEnergy>(m, "RunningEnergy")
      .def("conserved_energy", &RunningEnergy::conservedEnergy)
      .def("potential_energy", &RunningEnergy::potentialEnergy)
      .def("kinetic_energy", &RunningEnergy::kineticEnergy)
      .def("translational_kinetic_energy", &RunningEnergy::translationalPartKineticEnergy)
      .def("rotational_kinetic_energy", &RunningEnergy::rotationalPartKineticEnergy)
      .def("thermostat_energy", &RunningEnergy::thermostatEnergy)
      .def("coulomb_energy", &RunningEnergy::CoulombEnergy)
      .def("van_der_waals_energy", &RunningEnergy::VanDerWaalsEnergy)
      .def_rw("molecule_molecule_vdw", &RunningEnergy::moleculeMoleculeVDW)
      .def_rw("framework_molecule_vdw", &RunningEnergy::frameworkMoleculeVDW)
      .def("__repr__", &RunningEnergy::repr);

  nanobind::class_<System>(m, "System")
      .def(nanobind::init<ForceField, std::optional<SimulationBox>, bool, double, std::optional<double>, double,
                          std::optional<Framework>, std::vector<Component>, std::vector<std::vector<double3>>,
                          std::vector<std::size_t>, std::size_t, MCMoveProbabilities>(),
           nanobind::arg("force_field"), 
           nanobind::arg("simulation_box") = nanobind::none(),
           nanobind::arg("has_external_field") = false, 
           nanobind::arg("external_temperature") = 298.0, 
           nanobind::arg("external_pressure") = nanobind::none(), 
           nanobind::arg("helium_void_fraction") = 0.0,
           nanobind::arg("framework_components") = nanobind::none(),
           nanobind::arg("components") = std::vector<Component>(), 
           nanobind::arg("initial_positions") = std::vector<std::vector<double3>>(),
           nanobind::arg("initial_number_of_molecules") = std::vector<std::size_t>(),
           nanobind::arg("number_of_blocks") = 5,
           nanobind::arg("system_probabilities") = MCMoveProbabilities())
      .def("compute_total_energies", &System::computeTotalEnergies)
      .def("framework_mass", &System::frameworkMass)
      .def_ro("simulation_box", &System::simulationBox)
      .def_ro("input_pressure", &System::input_pressure)
      .def_ro("components", &System::components)
      .def_ro("number_of_molecules_per_component", &System::numberOfMoleculesPerComponent)
      .def_ro("mc_moves_statistics", &System::mc_moves_statistics)
      .def_ro("loadings", &System::loadings)
      .def_rw("average_loadings", &System::averageLoadings)
      .def_ro("average_energies", &System::averageEnergies)
      .def_ro("average_pressure", &System::averagePressure)
      .def_ro("thermostat", &System::thermostat)
      .def("set_thermostat", &System::setThermostat,
         nanobind::arg("thermostat") = nanobind::none())
      .def_ro("property_density_grid", &System::propertyDensityGrid)
      .def("set_property_density_grid", &System::setPropertyDensityGrid,
         nanobind::arg("property:") = nanobind::none())
      .def_ro("average_energy_histogram", &System::averageEnergyHistogram)
      .def("set_average_energy_histogram", &System::setAverageEnergyHistogram,
         nanobind::arg("property:") = nanobind::none())
      .def_ro("property_number_of_molecules_histogram", &System::averageNumberOfMoleculesHistogram)
      .def("set_number_of_molecules_histogram", &System::setNumberOfMoleculesHistogram,
         nanobind::arg("property:") = nanobind::none())
      .def_ro("sample_pdb_movie", &System::samplePDBMovie)
      .def("set_sample_pdb_movie", &System::setSamplePDBMovie,
         nanobind::arg("property:") = nanobind::none())
      .def_ro("property_conventional_rdf", &System::propertyConventionalRadialDistributionFunction)
      .def("set_property_conventional_rdf", &System::setPropertyConventionalRDF,
         nanobind::arg("property") = nanobind::none())
      .def_ro("property_msd", &System::propertyMSD)
      .def("set_property_msd", &System::setPropertyMSD,
         nanobind::arg("property") = nanobind::none())
      .def_ro("property_vacf", &System::propertyVACF)
      .def("set_property_vacf", &System::setPropertyVACF,
         nanobind::arg("property") = nanobind::none())
      .def_ro("property_rdf", &System::propertyRadialDistributionFunction)
      .def("set_property_rdf", &System::setPropertyRDF,
         nanobind::arg("property") = nanobind::none())
      .def_ro("property_number_of_molecules_evolution", &System::propertyNumberOfMoleculesEvolution)
      .def("set_property_number_of_molecules_evolution", &System::setPropertyNumberOfMoleculesEvolution,
         nanobind::arg("property") = nanobind::none())
      .def_ro("property_volume_evolution", &System::propertyVolumeEvolution)
      .def("set_property_volume_evolution", &System::setPropertyVolumeEvolution,
         nanobind::arg("property") = nanobind::none())
      .def_ro("property_conserved_energy_evolution", &System::propertyConservedEnergyEvolution)
      .def("set_property_conserved_energy_evolution", &System::setPropertyConservedEnergyEvolution,
         nanobind::arg("property") = nanobind::none())
      .def_rw("average_enthalpies_of_adsorption", &System::averageEnthalpiesOfAdsorption)
      .def_rw("average_partial_molar_properties", &System::averagePartialMolarProperties)
      .def_rw("atom_data", &System::atomData)
      .def_rw("translational_degrees_of_freedom", &System::translationalDegreesOfFreedom)
      .def_rw("rotational_degrees_of_freedom", &System::rotationalDegreesOfFreedom)
      .def_rw("time_step", &System::timeStep)
      .def("write_mc_move_statistics", &System::writeMCMoveStatistics)
      .def("__repr__", &System::repr);

    nanobind::class_<MonteCarlo>(m, "MonteCarlo")
      .def("__init__",
         [](MonteCarlo *mc, std::size_t number_of_production_cycles, std::size_t number_of_pre_initialization_cycles,
            std::size_t number_of_initialization_cycles, std::size_t number_of_equilibration_cycles,
            std::size_t print_every, std::size_t write_binary_restart_every, std::size_t rescale_wang_landau_every,
            std::size_t optimize_mc_moves_every, const std::vector<System> &systems,
            std::optional<std::size_t> random_seed, std::size_t number_of_blocks, bool output_to_files)
         {
           new (mc) MonteCarlo(
               SimulationSchedule{number_of_production_cycles, number_of_pre_initialization_cycles,
                                  number_of_initialization_cycles, number_of_equilibration_cycles, print_every,
                                  write_binary_restart_every, rescale_wang_landau_every, optimize_mc_moves_every},
               systems, random_seed, number_of_blocks, output_to_files);
         },
         nanobind::arg("number_of_production_cycles"),
         nanobind::arg("number_of_pre_initialization_cycles") = 0,
         nanobind::arg("number_of_initialization_cycles") = 0,
         nanobind::arg("number_of_equilibration_cycles") = 0,
         nanobind::arg("print_every") = 5000,
         nanobind::arg("write_binary_restart_every") = 5000,
         nanobind::arg("rescale_wang_landau_every") = 1000,
         nanobind::arg("optimize_mc_moves_every") = 100,
         nanobind::arg("systems") = std::vector<System>(),
         nanobind::arg("random_seed") = nanobind::none(),
         nanobind::arg("number_of_blocks") = 5,
         nanobind::arg("output_to_files") = false)
      .def("run", &MonteCarlo::run)
      .def("setup", &MonteCarlo::setup)
      .def("tear_down", &MonteCarlo::tearDown)
      .def("pre_initialize", &MonteCarlo::preInitialize, nanobind::arg("call_back_function") = std::function<void()>([](){}), nanobind::arg("call_back_every") = 100)
      .def("initialize", &MonteCarlo::initialize, nanobind::arg("call_back_function") = std::function<void()>([](){}), nanobind::arg("call_back_every") = 100)
      .def("equilibrate", &MonteCarlo::equilibrate, nanobind::arg("call_back_function") = std::function<void()>([](){}), nanobind::arg("call_back_every") = 100)
      .def("production", &MonteCarlo::production, nanobind::arg("call_back_function") = std::function<void()>([](){}), nanobind::arg("call_back_every") = 100)
      .def("cycle", &MonteCarlo::performCycle)
      .def_ro("systems", &MonteCarlo::systems);


  nanobind::class_<MolecularDynamics>(m, "MolecularDynamics")
    .def("__init__",
         [](MolecularDynamics *md, std::size_t number_of_production_cycles, std::size_t number_of_pre_initialization_cycles,
            std::size_t number_of_initialization_cycles, std::size_t number_of_equilibration_cycles,
            std::size_t print_every, std::size_t write_binary_restart_every, std::size_t rescale_wang_landau_every,
            std::size_t optimize_mc_moves_every, const std::vector<System> &systems,
            std::optional<std::size_t> random_seed, std::size_t number_of_blocks, bool output_to_files)
         {
           new (md) MolecularDynamics(
               SimulationSchedule{number_of_production_cycles, number_of_pre_initialization_cycles,
                                  number_of_initialization_cycles, number_of_equilibration_cycles, print_every,
                                  write_binary_restart_every, rescale_wang_landau_every, optimize_mc_moves_every},
               systems, random_seed, number_of_blocks, output_to_files);
         },
         nanobind::arg("number_of_production_cycles"), 
         nanobind::arg("number_of_pre_initialization_cycles") = 0,
         nanobind::arg("number_of_initialization_cycles") = 0,
         nanobind::arg("number_of_equilibration_cycles") = 0, 
         nanobind::arg("print_every") = 5000,
         nanobind::arg("write_binary_restart_every") = 5000, 
         nanobind::arg("rescale_wang_landau_every") = 1000,
         nanobind::arg("optimize_mc_moves_every") = 100, 
         nanobind::arg("systems") = std::vector<System>(),
         nanobind::arg("random_seed") = nanobind::none(), 
         nanobind::arg("number_of_blocks") = 5,
         nanobind::arg("output_to_files") = false)
      .def("run", &MolecularDynamics::run)
      .def("setup", &MolecularDynamics::setup)
      .def("tear_down", &MolecularDynamics::tearDown)
      .def("pre_initialize", &MolecularDynamics::preInitialize, nanobind::arg("call_back_function") = 
                      nanobind::cpp_function([](void){}), nanobind::arg("call_back_every") = 100)
      .def("initialize", &MolecularDynamics::initialize, nanobind::arg("call_back_function") = 
                      nanobind::cpp_function([](void){}), nanobind::arg("call_back_every") = 100)
      .def("equilibrate", &MolecularDynamics::equilibrate, nanobind::arg("call_back_function") = 
                      nanobind::cpp_function([](void){}), nanobind::arg("call_back_every") = 100)
      .def("production", &MolecularDynamics::production, nanobind::arg("call_back_function") = 
                      nanobind::cpp_function([](void){}), nanobind::arg("call_back_every") = 100)
      .def_ro("systems", &MolecularDynamics::systems);


  nanobind::class_<PropertyConservedEnergyEvolution>(m, "PropertyConservedEnergyEvolution")
      .def(nanobind::init<std::size_t, std::size_t, std::optional<std::size_t>>(),
            nanobind::arg("number_of_production_cycles"), 
            nanobind::arg("sample_every"), 
            nanobind::arg("write_every") = nanobind::none())
       // convert result to units of Kelvin
      .def("result", [](const PropertyConservedEnergyEvolution& p) { return Units::EnergyToKelvin * p.result();});

  nanobind::class_<Thermostat>(m, "Thermostat")
      .def(nanobind::init<std::size_t, std::size_t, double>(),
            nanobind::arg("thermostat_chain_length") = 5, 
            nanobind::arg("number_of_yoshida_suzuki_steps") = 5,
            nanobind::arg("time_scale_parameter") = 0.15);

  nanobind::class_<EnthalpyOfAdsorptionData>(m, "EnthalpyOfAdsorptionData")
      .def(nanobind::init<std::size_t>(),
            nanobind::arg("number_of_components"))
      .def_ro("values", &EnthalpyOfAdsorptionData::values)
      .def("__getitem__", &EnthalpyOfAdsorptionData::operator[]);

  nanobind::class_<PropertyEnthalpy>(m, "PropertyEnthalpy")
      .def(nanobind::init<std::size_t, std::size_t>(),
            nanobind::arg("number_of_blocks"), 
            nanobind::arg("number_of_components"))
      .def("result", [](const PropertyEnthalpy& p) { return Units::EnergyToKelvin * p.result();});

  nanobind::class_<PartialMolarPropertiesData>(m, "PartialMolarPropertiesData")
      .def(nanobind::init<std::size_t>(),
            nanobind::arg("number_of_components"))
      .def_ro("partial_molar_energy", &PartialMolarPropertiesData::partialMolarEnergy)
      .def_ro("partial_molar_volume", &PartialMolarPropertiesData::partialMolarVolume);

  nanobind::class_<PropertyPartialMolarProperties>(m, "PropertyPartialMolarProperties")
      .def(nanobind::init<std::size_t, std::size_t>(),
            nanobind::arg("number_of_blocks"),
            nanobind::arg("number_of_components"))
      .def("average_properties", [](const PropertyPartialMolarProperties& p) { return p.averageProperties(); });
}
