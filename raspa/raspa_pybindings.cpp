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
import framework;
import component;
import system;
import randomnumbers;
import monte_carlo;
import molecular_dynamics;
import simulation_schedule;
import input_reader;
import property_loading;
import energy_dudlambda;
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
import property_pressure;
import equation_of_states;
import property_conserved_energy_evolution;
import thermostat;
import property_enthalpy;
import property_partial_molar_properties;
import property_conventional_rdf;
import property_rdf;
import mean_squared_displacement_data;
import property_msd;
import velocity_autocorrelation_function_data;
import property_vacf;
import property_temperature;
import property_simulationbox;
import thermobarostat;
import minimization_cell_layout;
import archive;
import threadpool;
import spatial_decomposition_settings;
import bond_potential;
import urey_bradley_potential;
import bend_potential;
import inversion_bend_potential;
import torsion_potential;
import bond_bond_potential;
import bond_bend_potential;
import bond_torsion_potential;
import bend_bend_potential;
import bend_torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import cmap_potential;
import chiral_center;
import fragment_graph;
import fragment;
import json;
import stringutils;
import molecule_property_settings;
import property_molecule_properties;
import property_molecule_shape;
import property_molecule_backbone;
import property_end_to_end_acf;
import property_elastic_constants_fluctuation;
import property_radical_waiting_times;
import property_polymer_chains;
import write_lammps_data;
import reaction;
import reactions;
import partition_function;

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

// The conversion of one Python number inside a caster. 'from_python' is noexcept, and nanobind::cast throws for a
// number that does not fit (10**20 for an int, 10**400 for a float): the exception left the caster and ended the
// interpreter in std::terminate. Here it is a failed conversion, a TypeError for the caller.
template <typename Source, typename Number>
bool castNumberNoThrow(const Source &source, Number &number) noexcept
{
  try
  {
    number = nanobind::cast<Number>(source);
    return true;
  }
  catch (...)
  {
    PyErr_Clear();
    return false;
  }
}

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
    if (!castNumberNoThrow(seq[0], value.x)) return false;
    if (!castNumberNoThrow(seq[1], value.y)) return false;
    if (!castNumberNoThrow(seq[2], value.z)) return false;
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
      double scalar{};
      if (!castNumberNoThrow(src, scalar)) return false;
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

        double scalar{};

        if (!castNumberNoThrow(item, scalar)) return false;
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
        if (!castNumberNoThrow(seq[0], value.x)) return false;
        if (!castNumberNoThrow(seq[1], value.y)) return false;
        if (!castNumberNoThrow(seq[2], value.z)) return false;
        return true;
      }
    }

    return false;
  }
};

// A 3x3 matrix crosses as three sequences of three numbers, the box vectors a, b and c (the columns of the cell
// matrix), and comes back as three tuples the same way. Without it the SimulationBox(double3x3) constructor could
// not be called from Python at all.
template <>
struct type_caster<double3x3>
{
 public:
  NB_TYPE_CASTER(double3x3, io_name("Sequence[Sequence[float]]",
                                    "tuple[tuple[float, float, float], "
                                    "tuple[float, float, float], tuple[float, float, float]]"));

  static handle from_cpp(const double3x3& m, rv_policy /* policy */, cleanup_list* /* cleanup */) noexcept
  {
    return nanobind::make_tuple(nanobind::make_tuple(m.ax, m.ay, m.az), nanobind::make_tuple(m.bx, m.by, m.bz),
                                nanobind::make_tuple(m.cx, m.cy, m.cz))
        .release();
  }

  bool from_python(handle src, uint8_t /* flags */, cleanup_list* /* cleanup */) noexcept
  {
    if (!nanobind::isinstance<nanobind::sequence>(src)) return false;
    auto rows = nanobind::borrow<nanobind::sequence>(src);
    if (nanobind::len(rows) != 3) return false;
    double entries[3][3];
    for (std::size_t vector = 0; vector < 3; ++vector)
    {
      auto row = rows[vector];
      if (!nanobind::isinstance<nanobind::sequence>(row)) return false;
      auto items = nanobind::borrow<nanobind::sequence>(row);
      if (nanobind::len(items) != 3) return false;
      for (std::size_t component = 0; component < 3; ++component)
      {
        auto item = items[component];
        if (!nanobind::isinstance<nanobind::float_>(item) && !nanobind::isinstance<nanobind::int_>(item)) return false;
        if (!castNumberNoThrow(item, entries[vector][component])) return false;
      }
    }
    value = double3x3(double3(entries[0][0], entries[0][1], entries[0][2]),
                      double3(entries[1][0], entries[1][1], entries[1][2]),
                      double3(entries[2][0], entries[2][1], entries[2][2]));
    return true;
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
      double scalar{};
      if (!castNumberNoThrow(src, scalar)) return false;
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

        double scalar{};

        if (!castNumberNoThrow(item, scalar)) return false;
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
        if (!castNumberNoThrow(seq[0], value.x)) return false;
        if (!castNumberNoThrow(seq[1], value.y)) return false;
        if (!castNumberNoThrow(seq[2], value.z)) return false;
        if (!castNumberNoThrow(seq[3], value.w)) return false;
        return true;
      }
    }

    return false;
  }
};


} // namespace detail
} // namespace nanobind


// The move types under the names they have in Python. The probability of a move goes by the lower-case name with
// '_probability' appended ('translation_probability' for TRANSLATION).
const std::vector<std::pair<const char*, Move::Types>> moveTypeNames = {
    {"TRANSLATION", Move::Types::Translation},
    {"RANDOM_TRANSLATION", Move::Types::RandomTranslation},
    {"ROTATION", Move::Types::Rotation},
    {"RANDOM_ROTATION", Move::Types::RandomRotation},
    {"TRANSLATION_SMART_MC", Move::Types::TranslationSmartMC},
    {"TRANSLATION_SMART_MC_ALL", Move::Types::TranslationSmartMCAll},
    {"ROTATION_SMART_MC", Move::Types::RotationSmartMC},
    {"ROTATION_SMART_MC_ALL", Move::Types::RotationSmartMCAll},
    {"TRANSLATION_ROTATION_SMART_MC", Move::Types::TranslationRotationSmartMC},
    {"VOLUME_CHANGE", Move::Types::VolumeChange},
    {"ANISOTROPIC_VOLUME_CHANGE", Move::Types::AnisotropicVolumeChange},
    {"REINSERTION_CBMC", Move::Types::ReinsertionCBMC},
    {"PARTIAL_REINSERTION_CBMC", Move::Types::PartialReinsertionCBMC},
    {"IDENTITY_CHANGE_CBMC", Move::Types::IdentityChangeCBMC},
    {"IDENTITY_SWITCH_CBMC", Move::Types::IdentitySwitchCBMC},
    {"SWAP", Move::Types::Swap},
    {"SWAP_CBMC", Move::Types::SwapCBMC},
    {"SWAP_CFCMC", Move::Types::SwapCFCMC},
    {"SWAP_CBCFCMC", Move::Types::SwapCBCFCMC},
    {"GIBBS_VOLUME", Move::Types::GibbsVolume},
    {"GIBBS_SWAP_CBMC", Move::Types::GibbsSwapCBMC},
    {"GIBBS_CONVENTIONAL_CFCMC", Move::Types::GibbsConventionalCFCMC},
    {"GIBBS_CONVENTIONAL_CBCFCMC", Move::Types::GibbsConventionalCBCFCMC},
    {"GIBBS_SWAP_CFCMC", Move::Types::GibbsSwapCFCMC},
    {"GIBBS_SWAP_CBCFCMC", Move::Types::GibbsSwapCBCFCMC},
    {"GIBBS_IDENTITY_CHANGE_CBMC", Move::Types::GibbsIdentityChangeCBMC},
    {"WIDOM", Move::Types::Widom},
    {"WIDOM_CFCMC", Move::Types::WidomCFCMC},
    {"WIDOM_CBCFCMC", Move::Types::WidomCBCFCMC},
    {"PARALLEL_TEMPERING", Move::Types::ParallelTempering},
    {"HYBRID_MC", Move::Types::HybridMC},
    {"REACTION_CBMC", Move::Types::ReactionCBMC},
    {"REACTION_CONVENTIONAL_CFCMC", Move::Types::ReactionConventionalCFCMC},
    {"REACTION_CONVENTIONAL_CBCFCMC", Move::Types::ReactionConventionalCBCFCMC},
    {"REACTION_CFCMC", Move::Types::ReactionCFCMC},
    {"REACTION_CBCFCMC", Move::Types::ReactionCBCFCMC},
    {"PAIR_SWAP", Move::Types::PairSwap},
    {"PAIR_SWAP_CBMC", Move::Types::PairSwapCBMC},
    {"PAIR_SWAP_CFCMC", Move::Types::PairSwapCFCMC},
    {"PAIR_SWAP_CBCFCMC", Move::Types::PairSwapCBCFCMC},
    {"GROUP_SWAP", Move::Types::GroupSwap},
    {"GROUP_SWAP_CBMC", Move::Types::GroupSwapCBMC},
    {"GROUP_SWAP_CFCMC", Move::Types::GroupSwapCFCMC},
    {"GROUP_SWAP_CBCFCMC", Move::Types::GroupSwapCBCFCMC},
    {"TETHERED_PROTON_HOP", Move::Types::TetheredProtonHop},
    {"PIVOT", Move::Types::Pivot},
    {"CRANKSHAFT", Move::Types::Crankshaft},
    {"REPTATION", Move::Types::Reptation},
    {"CONCERTED_ROTATION", Move::Types::ConcertedRotation},
    {"BEAD_DISPLACEMENT", Move::Types::BeadDisplacement},
    {"BEAD_FLIP", Move::Types::BeadFlip},
    {"DOUBLE_BRIDGING", Move::Types::DoubleBridging},
    {"INTRAMOLECULAR_DOUBLE_REBRIDGING", Move::Types::IntramolecularDoubleRebridging},
    {"CROSS_LINK_SWAP", Move::Types::CrossLinkSwap},
    {"CROSS_LINK_FORMATION_SCISSION", Move::Types::CrossLinkFormationScission},
    {"CROSS_LINK_EXCHANGE", Move::Types::CrossLinkExchange},
    {"PIVOT_CBMC", Move::Types::PivotCBMC},
    {"SWAP_MULTI_TRIAL_CBMC", Move::Types::SwapMultiTrialCBMC},
    {"REINSERTION_MULTI_TRIAL_CBMC", Move::Types::ReinsertionMultiTrialCBMC},
    {"GIBBS_SWAP_MULTI_TRIAL_CBMC", Move::Types::GibbsSwapMultiTrialCBMC},
};

// two probabilities had a shorter keyword before every move had one
const std::vector<std::pair<const char*, Move::Types>> moveProbabilityAliases = {
    {"identity_change_probability", Move::Types::IdentityChangeCBMC},
    {"gibbs_volume_change_probability", Move::Types::GibbsVolume},
};

std::string moveProbabilityKeyword(const char* moveName)
{
  std::string keyword(moveName);
  std::ranges::transform(keyword, keyword.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return keyword + "_probability";
}

Move::Types moveTypeOfProbabilityKeyword(const std::string& keyword)
{
  for (const auto& [moveName, moveType] : moveTypeNames)
  {
    if (keyword == moveProbabilityKeyword(moveName)) return moveType;
  }
  for (const auto& [alias, moveType] : moveProbabilityAliases)
  {
    if (keyword == alias) return moveType;
  }
  throw nanobind::type_error(
      std::format("MCMoveProbabilities: unknown keyword '{}'; a probability is the name of a move of Move.Types in "
                  "lower case with '_probability' appended, like 'translation_probability'",
                  keyword)
          .c_str());
}

void setMoveProbability(MCMoveProbabilities& probabilities, Move::Types move, double probability)
{
  if (!(probability >= 0.0))
  {
    throw std::invalid_argument("MCMoveProbabilities: a probability can not be negative");
  }
  probabilities.setProbability(move, probability);
}

// the signature shown for the constructor of MCMoveProbabilities: one keyword for every move
std::string moveProbabilitiesSignature()
{
  std::string signature = "def __init__(self, *";
  for (const auto& [moveName, moveType] : moveTypeNames)
  {
    signature += std::format(", {}: float = 0.0", moveProbabilityKeyword(moveName));
  }
  return signature + ") -> None";
}

// a property of a class that is a fraction, and one that is a count of at least one
template <typename Class>
void defineFraction(nanobind::class_<Class>& bound, const char* name, double Class::* member, const char* doc)
{
  bound.def_prop_rw(
      name, [member](const Class& object) { return object.*member; },
      [member, name](Class& object, double value)
      {
        if (!(value >= 0.0 && value <= 1.0)) throw std::invalid_argument(std::format("{} must be in [0, 1]", name));
        object.*member = value;
      },
      doc);
}

template <typename Class>
void defineCount(nanobind::class_<Class>& bound, const char* name, std::size_t Class::* member, const char* doc)
{
  bound.def_prop_rw(
      name, [member](const Class& object) { return object.*member; },
      [member, name](Class& object, std::size_t value)
      {
        if (value == 0) throw std::invalid_argument(std::format("{} must be at least 1", name));
        object.*member = value;
      },
      doc);
}

// The parameters of every type of a bonded potential, as the potential itself prints them (with their units): a
// term of each type is constructed and its description is kept without the atoms and the values.
template <typename Potential, typename Type, typename NumberOfParameters>
std::string describeParameters(const std::vector<std::pair<const char*, Type>>& types,
                               const NumberOfParameters& numberOfParameters)
{
  using Identifiers = decltype(Potential::identifiers);
  Identifiers identifiers{};
  std::iota(identifiers.begin(), identifiers.end(), std::size_t{0});

  std::string description = "The types and their parameters:\n";
  for (const auto& [name, type] : types)
  {
    const std::size_t count = numberOfParameters[static_cast<std::size_t>(std::to_underlying(type))];
    std::string parameters = count == 0 ? "no parameters" : std::format("{} parameters", count);
    try
    {
      std::string printed = Potential(identifiers, type, std::vector<double>(count, 1.0)).print();
      // "0 - 1 : HARMONIC p_0/k_B=1 [K/A^2], p_1=1 [A]": what follows the name of the type, without the values
      const std::size_t first = printed.find("p_0");
      if (first != std::string::npos)
      {
        printed = printed.substr(first);
        while (!printed.empty() && (printed.back() == '\n' || printed.back() == ' ')) printed.pop_back();
        parameters = std::regex_replace(printed, std::regex("=[-+0-9.eE]+"), "");
      }
    }
    catch (const std::exception&)
    {
    }
    description += std::format("  {}: {}\n", name, parameters);
  }
  return description;
}

// The bonded potentials have one shape: the atoms of the term, a type, and the parameters of that type in the units
// of the molecule files (energies in Kelvin, lengths in Angstrom, angles in degrees). 'numberOfParameters' holds the
// number of parameters of each type, in the order of the enumeration.
template <typename Potential, typename Type, typename NumberOfParameters>
void bindIntraMolecularPotential(nanobind::module_& m, const char* className, const char* typeName,
                                 const std::vector<std::pair<const char*, Type>>& types,
                                 const NumberOfParameters& numberOfParameters)
{
  using Identifiers = decltype(Potential::identifiers);

  // (the documentation of a class is kept by its address)
  static std::deque<std::string> documentation;
  const std::string& parametersOfTypes =
      documentation.emplace_back(describeParameters<Potential, Type>(types, numberOfParameters));

  nanobind::enum_<Type> boundTypes(m, typeName, nanobind::is_arithmetic());
  for (const auto& [name, type] : types) boundTypes.value(name, type);

  const std::string name(className);
  nanobind::class_<Potential>(m, className, parametersOfTypes.c_str())
      .def(
          "__init__",
          [name, numberOfParameters](Potential* potential, Identifiers identifiers, Type type,
                                     std::vector<double> parameters)
          {
            for (std::size_t i = 0; i < identifiers.size(); ++i)
            {
              for (std::size_t j = i + 1; j < identifiers.size(); ++j)
              {
                if (identifiers[i] == identifiers[j])
                {
                  throw std::invalid_argument(std::format("{}: atom {} is named twice", name, identifiers[i]));
                }
              }
            }
            const std::size_t index = static_cast<std::size_t>(std::to_underlying(type));
            if (index >= numberOfParameters.size())
            {
              throw std::invalid_argument(std::format("{}: unknown type", name));
            }
            if (parameters.size() != numberOfParameters[index])
            {
              throw std::invalid_argument(std::format("{}: the type takes {} parameters, {} were given", name,
                                                      numberOfParameters[index], parameters.size()));
            }
            try
            {
              new (potential) Potential(identifiers, type, parameters);
            }
            catch (const std::runtime_error& error)
            {
              // (a type the potential does not support, for instance)
              throw std::invalid_argument(error.what());
            }
          },
          nanobind::arg("identifiers"), nanobind::arg("type"), nanobind::arg("parameters") = std::vector<double>{})
      .def_ro("identifiers", &Potential::identifiers, "The atoms of the term, as indices into the atoms of the component.")
      .def_ro("type", &Potential::type)
      .def_static(
          "describe_parameters", [&parametersOfTypes]() { return parametersOfTypes; },
          "The parameters of every type, with their units.")
      .def("__repr__", &Potential::print);
}

// The connectivity of a component. Without rigid bodies the bonds are the connections, and a table that is given has
// to be that of the bonds. With rigid bodies the table is given: the connections inside a rigid body have no bond
// potential, since the geometry of the body is fixed. As in a molecule file, every bond, bend and torsion of the
// connectivity that is not inside one rigid body needs a potential; a term that names an atom the component does
// not have is refused.
ConnectivityTable connectivityOfComponent(std::size_t numberOfAtoms, const ConnectivityTable& given,
                                          const Potentials::IntraMolecularPotentials& potentials,
                                          const std::vector<std::vector<std::size_t>>& rigidBodies)
{
  auto checkAtoms = [numberOfAtoms](const char* kind, const auto& terms)
  {
    for (const auto& term : terms)
    {
      for (std::size_t atom : term.identifiers)
      {
        if (atom >= numberOfAtoms)
        {
          throw std::invalid_argument(
              std::format("Component: a {} names atom {}, but the component has {} atoms", kind, atom, numberOfAtoms));
        }
      }
    }
  };
  checkAtoms("bond", potentials.bonds);
  checkAtoms("Urey-Bradley term", potentials.ureyBradleys);
  checkAtoms("bend", potentials.bends);
  checkAtoms("inversion bend", potentials.inversionBends);
  checkAtoms("torsion", potentials.torsions);
  checkAtoms("improper torsion", potentials.improperTorsions);
  checkAtoms("bond-bond term", potentials.bondBonds);
  checkAtoms("bond-bend term", potentials.bondBends);
  checkAtoms("bond-torsion term", potentials.bondTorsions);
  checkAtoms("bend-bend term", potentials.bendBends);
  checkAtoms("bend-torsion term", potentials.bendTorsions);
  checkAtoms("CMAP term", potentials.cmaps);
  checkAtoms("Van der Waals pair", potentials.vanDerWaals);
  checkAtoms("Coulomb pair", potentials.coulombs);
  for (const CMAPPotential& term : potentials.cmaps)
  {
    if (term.mapIndex >= potentials.cmapMaps.size())
    {
      throw std::invalid_argument(std::format("Component: a CMAP term names map {}, but there are {} maps",
                                              term.mapIndex, potentials.cmapMaps.size()));
    }
  }

  // the rigid body an atom is part of, if any (a body of one atom is a flexible atom)
  std::vector<std::optional<std::size_t>> bodyOfAtom(numberOfAtoms);
  try
  {
    FragmentGraph::partitionAtoms(numberOfAtoms, rigidBodies);
  }
  catch (const std::exception& error)
  {
    throw std::invalid_argument(std::format("Component: rigid_bodies: {}", error.what()));
  }
  for (std::size_t body = 0; body < rigidBodies.size(); ++body)
  {
    if (rigidBodies[body].size() < 2) continue;
    for (std::size_t atom : rigidBodies[body]) bodyOfAtom[atom] = body;
  }
  auto insideOneRigidBody = [&](const auto& atoms)
  {
    const std::optional<std::size_t> body = bodyOfAtom[atoms[0]];
    return body.has_value() && std::ranges::all_of(atoms, [&](std::size_t atom) { return bodyOfAtom[atom] == body; });
  };
  auto sameAtoms = [](const auto& a, const auto& b)
  { return std::ranges::equal(a, b) || std::ranges::equal(a, b | std::views::reverse); };

  ConnectivityTable bonds(numberOfAtoms);
  for (const BondPotential& bond : potentials.bonds)
  {
    bonds[bond.identifiers[0], bond.identifiers[1]] = true;
    bonds[bond.identifiers[1], bond.identifiers[0]] = true;
  }

  // (an empty table is one that was not given)
  const bool tableGiven = given.numberOfBeads != 0;
  if (tableGiven && given.numberOfBeads != numberOfAtoms)
  {
    throw std::invalid_argument(std::format("Component: connectivity_table is of {} atoms, the component has {}",
                                            given.numberOfBeads, numberOfAtoms));
  }
  const bool hasRigidBodies = std::ranges::any_of(bodyOfAtom, [](const auto& body) { return body.has_value(); });
  if (!hasRigidBodies)
  {
    if (tableGiven && given.table != bonds.table)
    {
      throw std::invalid_argument(
          "Component: connectivity_table differs from the bonds of intra_molecular_potentials; it can be left out");
    }
  }
  else if (!tableGiven)
  {
    throw std::invalid_argument("Component: with rigid_bodies the connectivity_table is given: the connections inside "
                                "a rigid body have no bond potential to take them from");
  }
  const ConnectivityTable table = hasRigidBodies ? given : bonds;

  for (const std::array<std::size_t, 2>& bond : table.findAllBonds())
  {
    if (insideOneRigidBody(bond)) continue;
    if (!bonds[bond[0], bond[1]])
    {
      throw std::invalid_argument(std::format("Component: the bond {}-{} has no bond potential", bond[0], bond[1]));
    }
  }
  for (const BondPotential& bond : potentials.bonds)
  {
    if (!table[bond.identifiers[0], bond.identifiers[1]])
    {
      throw std::invalid_argument(std::format("Component: the bond potential {}-{} is not a connection of "
                                              "connectivity_table",
                                              bond.identifiers[0], bond.identifiers[1]));
    }
  }
  for (const std::array<std::size_t, 3>& bend : table.findAllBends())
  {
    if (insideOneRigidBody(bend)) continue;
    if (std::ranges::none_of(potentials.bends,
                             [&](const BendPotential& potential) { return sameAtoms(potential.identifiers, bend); }))
    {
      throw std::invalid_argument(
          std::format("Component: the bend {}-{}-{} has no bend potential", bend[0], bend[1], bend[2]));
    }
  }
  for (const std::array<std::size_t, 4>& torsion : table.findAllTorsions())
  {
    if (insideOneRigidBody(torsion)) continue;
    if (std::ranges::none_of(potentials.torsions, [&](const TorsionPotential& potential)
                             { return sameAtoms(potential.identifiers, torsion); }))
    {
      throw std::invalid_argument(std::format("Component: the torsion {}-{}-{}-{} has no torsion potential",
                                              torsion[0], torsion[1], torsion[2], torsion[3]));
    }
  }
  return table;
}

// The types of the bonds, bends and torsions under the names they have in Python (the names of the molecule and
// framework-definition files).
const std::vector<std::pair<const char*, BondType>> bondTypeNames = {{"NONE", BondType::None}, {"FIXED", BondType::Fixed}, {"HARMONIC", BondType::Harmonic}, {"CORE_SHELL_SPRING", BondType::CoreShellSpring}, {"MORSE", BondType::Morse}, {"LJ_12_6", BondType::LJ_12_6}, {"LENNARD_JONES", BondType::LennardJones}, {"BUCKINGHAM", BondType::Buckingham}, {"RESTRAINED_HARMONIC", BondType::RestrainedHarmonic}, {"QUARTIC", BondType::Quartic}, {"CFF_QUARTIC", BondType::CFF_Quartic}, {"MM3", BondType::MM3}};
const std::vector<std::pair<const char*, BendType>> bendTypeNames = {{"RIGID", BendType::Rigid}, {"FIXED", BendType::Fixed}, {"HARMONIC", BendType::Harmonic}, {"CORE_SHELL", BendType::CoreShell}, {"QUARTIC", BendType::Quartic}, {"CFF_QUARTIC", BendType::CFF_Quartic}, {"HARMONIC_COSINE", BendType::HarmonicCosine}, {"COSINE", BendType::Cosine}, {"TAFIPOLSKY", BendType::Tafipolsky}, {"MM3", BendType::MM3}};
const std::vector<std::pair<const char*, TorsionType>> torsionTypeNames = {{"FIXED", TorsionType::Fixed}, {"HARMONIC", TorsionType::Harmonic}, {"HARMONIC_COSINE", TorsionType::HarmonicCosine}, {"THREE_COSINE", TorsionType::ThreeCosine}, {"RYCKAERT_BELLEMANS", TorsionType::RyckaertBellemans}, {"TRAPPE", TorsionType::TraPPE}, {"TRAPPE_EXTENDED", TorsionType::TraPPE_Extended}, {"MOD_TRAPPE", TorsionType::ModifiedTraPPE}, {"CVFF", TorsionType::CVFF}, {"CFF", TorsionType::CFF}, {"CFF2", TorsionType::CFF2}, {"OPLS", TorsionType::OPLS}, {"MM3", TorsionType::MM3}, {"FOURIER_SERIES", TorsionType::FourierSeries}, {"FOURIER_SERIES2", TorsionType::FourierSeries2}, {"CVFF_BLOCKED", TorsionType::CVFFBlocked}, {"POLYNOMIAL", TorsionType::Polynomial}};

template <typename Type>
std::string nameOfType(const std::vector<std::pair<const char*, Type>>& names, Type type)
{
  for (const auto& [name, value] : names)
  {
    if (value == type) return name;
  }
  throw std::invalid_argument("FrameworkDefinition: a potential type that is not available");
}

// What makes a framework flexible: the content of a framework-definition file. The potentials are given per
// combination of pseudo-atom types and apply to every bond, bend and torsion of those types the framework has.
struct FrameworkDefinition
{
  nlohmann::basic_json<nlohmann::raspa_map> data{};
};

// A framework from the content of a CIF file, as the input reader makes it: the cell, the space group and the atoms
// are read with the charges from the requested source; without a number of unit cells the smallest supercell is
// taken in which the cutoffs fit (the minimum-image convention); then the definition is applied.
Framework frameworkFromCIF(const ForceField& forceField, const std::string& content, const std::string& name,
                           std::optional<int3> numberOfUnitCells, CIFReader::UseChargesFrom useChargesFrom,
                           const std::optional<FrameworkDefinition>& definition,
                           const std::optional<std::string>& definitionFile)
{
  if (definition.has_value() && definitionFile.has_value())
  {
    throw std::invalid_argument("Framework: give either a definition or a definition_file, not both");
  }
  try
  {
    const auto cif = CIFReader::readCIFString(content, forceField, useChargesFrom);
    if (!cif.has_value())
    {
      throw std::invalid_argument(
          cif.error() == CIFReader::ParseError::invalidForceField
              ? std::format("Framework: '{}' has atoms the force field has no pseudo-atom for", name)
              : std::format("Framework: '{}' could not be read as a CIF file", name));
    }
    const auto& [simulationBox, spaceGroupHallNumber, definedAtoms, fractionalAtomsUnitCell] = cif.value();
    if (definedAtoms.empty())
    {
      throw std::invalid_argument(std::format("Framework: no atoms were read from '{}' as a CIF file", name));
    }

    int3 cells{};
    if (numberOfUnitCells.has_value())
    {
      cells = *numberOfUnitCells;
      if (cells.x < 1 || cells.y < 1 || cells.z < 1)
      {
        throw std::invalid_argument("Framework: number_of_unit_cells must be at least 1 in every direction");
      }
    }
    else
    {
      const double cutOff =
          std::max({forceField.cutOffFrameworkVDWAutomatic ? 0.0 : forceField.cutOffFrameworkVDW,
                      forceField.cutOffMoleculeVDWAutomatic ? 0.0 : forceField.cutOffMoleculeVDW,
                      forceField.useCharge && !forceField.cutOffCoulombAutomatic ? forceField.cutOffCoulomb : 0.0});
      cells = simulationBox.smallestNumberOfUnitCellsForMinimumImagesConvention(cutOff);
      cells.x = std::max(cells.x, std::int32_t{1});
      cells.y = std::max(cells.y, std::int32_t{1});
      cells.z = std::max(cells.z, std::int32_t{1});
    }

    Framework framework(forceField, name, simulationBox, spaceGroupHallNumber, definedAtoms, fractionalAtomsUnitCell,
                        cells);
    if (definition.has_value()) framework.applyFrameworkDefinition(forceField, definition->data);
    if (definitionFile.has_value()) framework.readFrameworkDefinition(forceField, *definitionFile);
    return framework;
  }
  catch (const std::runtime_error& error)
  {
    throw std::invalid_argument(std::format("Framework: {}", error.what()));
  }
}

// What a molecule file says of a flexible molecule besides its potentials: the chiral centres, the fixed atoms of the
// partial reinsertions, the repeat units of a chain, and the two atoms of the end-to-end distance. These are read
// and checked by the reader of the molecule files itself, from a description with the keys of such a file.
void describeFlexibleComponent(Component& component, const std::vector<std::array<std::size_t, 4>>& chiralCenters,
                               const std::vector<std::vector<std::size_t>>& partialReinsertionFixedAtoms,
                               const std::vector<std::vector<std::size_t>>& repeatUnits,
                               const std::optional<std::array<std::size_t, 2>>& endToEndAtoms)
{
  if (component.rigid)
  {
    if (!chiralCenters.empty() || !partialReinsertionFixedAtoms.empty() || !repeatUnits.empty() ||
        endToEndAtoms.has_value())
    {
      throw std::invalid_argument("Component: chiral_centers, partial_reinsertion_fixed_atoms, repeat_units and "
                                  "end_to_end_atoms are of a flexible molecule; this one is a single rigid body");
    }
    return;
  }

  using Description = nlohmann::basic_json<nlohmann::raspa_map>;
  Description description = Description::object();
  if (!chiralCenters.empty()) description["ChiralCenters"] = chiralCenters;
  if (!partialReinsertionFixedAtoms.empty()) description["Partial-reinsertion"] = partialReinsertionFixedAtoms;
  if (!repeatUnits.empty()) description["RepeatUnits"] = repeatUnits;
  if (endToEndAtoms.has_value()) description["EndToEndAtoms"] = *endToEndAtoms;

  try
  {
    component.intraMolecularPotentials.chiralCenters = component.readChiralCenters(description);
    component.partialReinsertionFixedAtoms = component.readPartialReinsertionFixedAtoms(description);
    component.repeatUnits = component.readRepeatUnits(description);
    component.endToEndAtoms = component.determineEndToEndAtoms(description);
    // whether the part that is not fixed can be grown is found when its growth plan is made
    for (const std::vector<std::size_t>& fixedAtoms : component.partialReinsertionFixedAtoms)
    {
      component.growthPlan(fixedAtoms);
    }
  }
  catch (const std::runtime_error& error)
  {
    throw std::invalid_argument(std::format("Component: {}", error.what()));
  }
}

// A reaction that has not been given a move takes the one reaction move of its system that has a probability.
const Move::Types reactionMoveNotGiven = Move::Types::Count;

// Gives a system its reactions, as the input reader does for 'Reactions': a reaction is driven by one reaction move,
// the one it names or the only one of the system that has a probability; its stoichiometry has an entry for every
// component; and the reactions of the fractional-molecule moves have their fractional molecules made.
void setReactions(System& system, std::vector<Reaction> reactions)
{
  if (reactions.empty()) return;

  const std::array reactionMoves{Move::Types::ReactionCBMC, Move::Types::ReactionCFCMC, Move::Types::ReactionCBCFCMC,
                                 Move::Types::ReactionConventionalCFCMC, Move::Types::ReactionConventionalCBCFCMC};
  std::vector<Move::Types> enabled;
  for (const Move::Types move : reactionMoves)
  {
    if (system.mc_moves_probabilities.getProbability(move) > 0.0) enabled.push_back(move);
  }
  if (enabled.empty())
  {
    throw std::invalid_argument("System: reactions are given, but no reaction move has a probability in "
                                "system_probabilities (reaction_cbmc_probability, for instance)");
  }

  const std::size_t numberOfComponents = system.components.size();
  for (std::size_t id = 0; id < reactions.size(); ++id)
  {
    Reaction& reaction = reactions[id];
    reaction.id = id;
    if (reaction.reactantStoichiometry.size() != numberOfComponents ||
        reaction.productStoichiometry.size() != numberOfComponents)
    {
      throw std::invalid_argument(std::format("System: reaction {} has {} reactants and {} products; both lists have "
                                              "an entry for each of the {} components",
                                              id, reaction.reactantStoichiometry.size(),
                                              reaction.productStoichiometry.size(), numberOfComponents));
    }
    if (reaction.reactionMove == reactionMoveNotGiven)
    {
      if (enabled.size() > 1)
      {
        throw std::invalid_argument(std::format("System: several reaction moves have a probability; reaction {} has "
                                                "to name the move that drives it",
                                                id));
      }
      reaction.reactionMove = enabled.front();
    }
    else if (!std::ranges::contains(enabled, reaction.reactionMove))
    {
      throw std::invalid_argument(std::format("System: the move of reaction {} has no probability in "
                                              "system_probabilities",
                                              id));
    }
  }

  // a lambda coordinate that is followed takes a dU/dlambda group: those of the components and those of the reactions
  std::size_t numberOfDUDlambda{0};
  for (const Component& component : system.components)
  {
    for (const PropertyLambdaProbabilityHistogram* histogram :
         {&component.lambdaGC, &component.lambdaGibbs, &component.lambdaPairSwap, &component.lambdaPairSwapCB,
          &component.lambdaGroupSwap, &component.lambdaGroupSwapCB})
    {
      if (histogram->computeDUdlambda) ++numberOfDUDlambda;
    }
  }
  for (const Reaction& reaction : reactions)
  {
    if (reaction.isSerialRxCFC()) numberOfDUDlambda += 1;
    if (reaction.isParallelRxCFC()) numberOfDUDlambda += 2;
  }
  if (numberOfDUDlambda > maximumNumberOfDUDlambdaGroups)
  {
    throw std::invalid_argument(std::format("System: {} lambda coordinates are asked for, at most {} can be followed "
                                            "at the same time",
                                            numberOfDUDlambda, maximumNumberOfDUDlambdaGroups));
  }

  system.reactions.list = std::move(reactions);
  if (system.usesReactionConventionalCFCMC()) system.createReactionFractionalMolecules();
}

// The component of an analysis of molecules: one the system has and that asked for the analysis.
template <typename Property>
void checkSampledComponent(const Property& property, std::size_t component, const char* name)
{
  if (component >= property.numberOfComponents || !property.isSampled(component))
  {
    throw nanobind::index_error(std::format("{}: component {} is not sampled", name, component).c_str());
  }
}

// A system builds its own property from the settings of the one it is given ('set_property_msd' and the like), sized
// for its molecules, components and blocks: the object of the caller is not the one that samples. To let that object
// hold the results nevertheless, the binding remembers which Python object was given to which Python system, and
// assigns the property of the system to it when it is set and whenever a system is copied out of a driver. The
// references are kept here, in a table of the module with the systems as weak keys; the systems themselves know
// nothing of them.
nanobind::handle givenPropertiesTable;

using CopyOfGivenProperty = void (*)(const System&, nanobind::handle);
std::vector<std::pair<const char*, CopyOfGivenProperty>> copiesOfGivenProperties;

template <auto member>
using PropertyOfMember = typename std::remove_cvref_t<decltype(std::declval<System&>().*member)>::value_type;

template <auto member>
void copyPropertyTo(const System& system, nanobind::handle given)
{
  if ((system.*member).has_value()) nanobind::cast<PropertyOfMember<member>&>(given) = *(system.*member);
}

// The properties given to a system, by the name of their setter; none when nothing was given.
nanobind::object givenPropertiesOf(nanobind::handle system) { return givenPropertiesTable.attr("get")(system); }

// Remembers the object given to the setter 'name' of a system (or forgets it, for None), and brings it up to date.
template <auto member>
void attachGivenProperty(nanobind::handle system, const char* name, nanobind::handle property)
{
  nanobind::object given = givenPropertiesOf(system);
  if (property.is_none())
  {
    if (!given.is_none()) given.attr("pop")(name, nanobind::none());
    return;
  }
  if (given.is_none())
  {
    given = nanobind::dict();
    givenPropertiesTable.attr("__setitem__")(system, given);
  }
  given.attr("__setitem__")(name, property);
  copyPropertyTo<member>(nanobind::cast<const System&>(system), property);
}

template <auto member>
std::optional<PropertyOfMember<member>> givenPropertyOrNone(nanobind::handle property, const char* name)
{
  using Property = PropertyOfMember<member>;
  if (property.is_none()) return std::nullopt;
  if (!nanobind::isinstance<Property>(property))
  {
    throw nanobind::type_error(std::format("{}: the property is of another type", name).c_str());
  }
  return nanobind::cast<const Property&>(property);
}

// Binds a setter of a property of the system that remembers the object it is given.
template <auto member, typename Setter>
void bindGivenProperty(nanobind::class_<System>& system, const char* name, Setter setter)
{
  copiesOfGivenProperties.emplace_back(name, &copyPropertyTo<member>);
  system.def(
      name,
      [name, setter](nanobind::handle self, nanobind::handle property)
      {
        std::invoke(setter, nanobind::cast<System&>(self), givenPropertyOrNone<member>(property, name));
        attachGivenProperty<member>(self, name, property);
      },
      nanobind::arg("property") = nanobind::none());
}

// A driver works on its own copy of the systems it is given: the kits are value-based and keep no references to the
// objects of the caller. To let those objects hold the results nevertheless, the systems of the driver are assigned
// back to the Python objects they were copied from whenever a phase returns (copy in, copy out), and their properties
// to the objects that were given to the setters of those systems. Inside a call-back the phase has not returned yet,
// and the systems of the driver itself ('systems') are the ones to read.
template <typename Driver>
struct DriverOfGivenSystems : Driver
{
  using Driver::Driver;

  std::vector<nanobind::object> givenSystems;

  void copyOut()
  {
    if (givenSystems.size() != this->systems.size()) return;

    // an object that was given more than once, a system to the driver or a property to its systems, has no single
    // one to take and is left as it was
    auto once = [](const std::vector<nanobind::object>& objects, nanobind::handle object)
    { return std::ranges::count_if(objects, [&](const nanobind::object& item) { return item.is(object); }) == 1; };

    std::vector<nanobind::object> givenProperties;
    for (std::size_t i = 0; i < givenSystems.size(); ++i)
    {
      if (std::ranges::find_if(givenSystems.begin(), givenSystems.begin() + static_cast<std::ptrdiff_t>(i),
                               [&](const nanobind::object& item) { return item.is(givenSystems[i]); }) !=
          givenSystems.begin() + static_cast<std::ptrdiff_t>(i))
      {
        continue;
      }
      nanobind::object given = givenPropertiesOf(givenSystems[i]);
      if (given.is_none()) continue;
      for (nanobind::handle property : given.attr("values")()) givenProperties.push_back(nanobind::borrow(property));
    }

    for (std::size_t i = 0; i < givenSystems.size(); ++i)
    {
      if (!once(givenSystems, givenSystems[i])) continue;
      nanobind::cast<System&>(givenSystems[i]) = this->systems[i];

      nanobind::object given = givenPropertiesOf(givenSystems[i]);
      if (given.is_none()) continue;
      for (const auto& [name, copy] : copiesOfGivenProperties)
      {
        nanobind::object property = given.attr("get")(name);
        if (!property.is_none() && once(givenProperties, property)) copy(this->systems[i], property);
      }
    }
  }
};

// The constructors of an evolution of one value per cycle: over the cycles given (those of all phases, by either
// keyword), or over those of the driver when none are given.
template <typename Property>
void bindEvolutionConstructors(nanobind::class_<Property>& evolution)
{
  auto construct =
      [](Property* property, std::size_t number_of_cycles, std::size_t sample_every,
         std::optional<std::size_t> write_every)
  {
    if (sample_every == 0) throw std::invalid_argument("sample_every must be positive");
    new (property) Property(number_of_cycles, sample_every, write_every);
  };
  // 'number_of_cycles' is the keyword of the 3.1.0 release
  for (const char* cycles : {"number_of_production_cycles", "number_of_cycles"})
  {
    evolution.def("__init__", construct, nanobind::arg(cycles), nanobind::arg("sample_every") = 1,
                  nanobind::arg("write_every") = nanobind::none());
  }
  evolution
      .def(
          "__init__",
          [construct](Property* property, std::size_t sample_every, std::optional<std::size_t> write_every)
          { construct(property, 0, sample_every, write_every); },
          nanobind::arg("sample_every") = 1, nanobind::arg("write_every") = nanobind::none())
      .def_ro("sample_every", &Property::sampleEvery)
      .def_ro("write_every", &Property::writeEvery);
}

// What only the driver knows is filled in on its copy of a system: the cycles an evolution is recorded over (all
// phases, as the input reader sizes them) when the property was made without them, and the index of the system in
// the names of the LAMMPS files.
const std::size_t systemIdNotGiven = std::numeric_limits<std::size_t>::max();

void completeProperties(System& system, std::size_t systemId, std::size_t numberOfCycles)
{
  if (auto& property = system.propertyNumberOfMoleculesEvolution; property.has_value() && property->totalSamples == 0)
  {
    property = PropertyNumberOfMoleculesEvolution(numberOfCycles, system.components.size(), property->sampleEvery,
                                                  property->writeEvery);
  }
  if (auto& property = system.propertyVolumeEvolution; property.has_value() && property->totalSamples == 0)
  {
    property = PropertyVolumeEvolution(numberOfCycles, property->sampleEvery, property->writeEvery);
  }
  if (auto& property = system.propertyConservedEnergyEvolution; property.has_value() && property->totalSamples == 0)
  {
    property = PropertyConservedEnergyEvolution(numberOfCycles, property->sampleEvery, property->writeEvery);
  }
  if (auto& writer = system.writeLammpsData; writer.has_value() && writer->systemId == systemIdNotGiven)
  {
    writer = WriteLammpsData(systemId, writer->sampleEvery);
  }
}

template <typename Driver>
void bindDriver(nanobind::module_& m, const char* name)
{
  using Bound = DriverOfGivenSystems<Driver>;
  using CallBack = std::function<void()>;

  auto construct =
      [](Bound* driver, std::size_t number_of_production_cycles, std::size_t number_of_pre_initialization_cycles,
         std::size_t number_of_initialization_cycles, std::size_t number_of_equilibration_cycles,
         std::size_t print_every, std::size_t write_binary_restart_every, std::size_t rescale_wang_landau_every,
         std::size_t optimize_mc_moves_every, nanobind::typed<nanobind::sequence, System> systems,
         std::optional<std::size_t> random_seed, std::optional<std::size_t> blocks_given, bool output_to_files)
  {
    std::vector<nanobind::object> given;
    for (nanobind::handle item : systems)
    {
      if (!nanobind::isinstance<System>(item)) throw nanobind::type_error("systems: a sequence of System is expected");
      given.push_back(nanobind::borrow(item));
    }

    // The driver has one number of blocks. Given here, its copies of the systems are sized for it; otherwise it
    // is that of the systems, which then have to agree.
    if (blocks_given == 0uz) throw std::invalid_argument("number_of_blocks must be positive");
    std::size_t number_of_blocks = blocks_given.value_or(5);
    if (!blocks_given.has_value() && !given.empty())
    {
      number_of_blocks = nanobind::cast<const System&>(given.front()).numberOfBlocks();
      for (const nanobind::object& item : given)
      {
        if (nanobind::cast<const System&>(item).numberOfBlocks() != number_of_blocks)
        {
          throw std::invalid_argument("the systems differ in their number of blocks: give 'number_of_blocks' to "
                                      "the driver");
        }
      }
    }

    std::vector<System> copies;
    for (const nanobind::object& item : given)
    {
      copies.push_back(nanobind::cast<const System&>(item));
      if (!copies.back().hasNumberOfBlocks(number_of_blocks)) copies.back().setNumberOfBlocks(number_of_blocks);

      // The analyses of the molecules that the components ask for (their histograms, shape, backbone statistics and
      // end-to-end autocorrelation function) are created here, as the input reader creates them once the driver is
      // known: the autocorrelation function has its lags in picoseconds for molecular dynamics and in cycles for
      // Monte Carlo.
      System& system = copies.back();
      completeProperties(system, copies.size() - 1,
                         number_of_pre_initialization_cycles + number_of_initialization_cycles +
                             number_of_equilibration_cycles + number_of_production_cycles);
      const bool asked =std::ranges::any_of(system.components, [](const Component& component)
                                             {
                                               return component.moleculePropertiesSettings.has_value() ||
                                                      component.moleculeShapeSettings.has_value() ||
                                                      component.moleculeBackboneSettings.has_value() ||
                                                      component.endToEndACFSettings.has_value();
                                             });
      if (asked)
      {
        try
        {
          system.initializeMoleculeProperties(number_of_blocks, std::same_as<Driver, MolecularDynamics>);
        }
        catch (const std::runtime_error& error)
        {
          throw std::invalid_argument(error.what());
        }
      }
    }
    new (driver) Bound(
        SimulationSchedule{number_of_production_cycles, number_of_pre_initialization_cycles,
                           number_of_initialization_cycles, number_of_equilibration_cycles, print_every,
                           write_binary_restart_every, rescale_wang_landau_every, optimize_mc_moves_every},
        copies, random_seed, number_of_blocks, output_to_files);
    driver->givenSystems = std::move(given);
  };

  nanobind::class_<Bound> driver(m, name);

  // the production cycles also go by the keyword 'number_of_cycles' of the 3.1.0 release
  for (const char* cycles : {"number_of_production_cycles", "number_of_cycles"})
  {
    driver.def("__init__", construct, nanobind::arg(cycles), nanobind::arg("number_of_pre_initialization_cycles") = 0,
               nanobind::arg("number_of_initialization_cycles") = 0,
               nanobind::arg("number_of_equilibration_cycles") = 0, nanobind::arg("print_every") = 5000,
               nanobind::arg("write_binary_restart_every") = 5000, nanobind::arg("rescale_wang_landau_every") = 1000,
               nanobind::arg("optimize_mc_moves_every") = 100, nanobind::arg("systems") = nanobind::list(),
               nanobind::arg("random_seed") = nanobind::none(), nanobind::arg("number_of_blocks") = nanobind::none(),
               nanobind::arg("output_to_files") = false);
  }

  auto step = [&driver](const char* pythonName, void (Driver::*function)())
  {
    driver.def(pythonName, [function](Bound& d)
               {
                 (d.*function)();
                 d.copyOut();
               });
  };
  auto phase = [&driver](const char* pythonName, void (Driver::*function)(CallBack, std::size_t))
  {
    driver.def(
        pythonName,
        [function](Bound& d, std::optional<CallBack> call_back_function, std::size_t call_back_every)
        {
          // None, as the type stubs advertise, is no callback (nanobind refused it for a plain std::function)
          (d.*function)(call_back_function.has_value() ? std::move(call_back_function.value()) : CallBack([]() {}),
                        call_back_every);
          d.copyOut();
        },
        nanobind::arg("call_back_function").none() = nanobind::none(), nanobind::arg("call_back_every") = 100);
  };

  step("run", &Driver::run);
  step("setup", &Driver::setup);
  step("tear_down", &Driver::tearDown);
  phase("pre_initialize", &Driver::preInitialize);
  phase("initialize", &Driver::initialize);
  phase("equilibrate", &Driver::equilibrate);
  phase("production", &Driver::production);

  // continues from a binary restart file, as the command-line program does for 'RestartFromBinaryFile': the driver
  // is read back, and what the file does not hold is made again (the output files and the interpolation grids)
  driver.def(
      "restart_from_binary_file",
      [](Bound& d, const std::string& file_name)
      {
        // the device of the multiple-try grower belongs to the machine, not to the state: it is not in the file
        std::vector<PairDevice> devices;
        for (const System& system : d.systems) devices.push_back(system.cbmcDevice);
        const bool outputToFiles = d.outputToFiles;

        // The reader refills the vectors of systems, components and reactions, which frees the objects a Python
        // reference ('driver.systems[0]', 'system.components[0]') still points at: the restored state is moved
        // into the old objects instead, where their counts agree, so that those references stay valid and see it.
        std::vector<System> held = std::move(d.systems);
        readBinaryRestartFile(static_cast<Driver&>(d), file_name);
        if (held.size() == d.systems.size())
        {
          for (std::size_t i = 0; i < held.size(); ++i)
          {
            std::vector<Component> components = std::move(held[i].components);
            std::vector<Reaction> reactions = std::move(held[i].reactions.list);
            held[i] = std::move(d.systems[i]);
            if (components.size() == held[i].components.size())
            {
              for (std::size_t k = 0; k < components.size(); ++k) components[k] = std::move(held[i].components[k]);
              held[i].components = std::move(components);
            }
            if (reactions.size() == held[i].reactions.list.size())
            {
              for (std::size_t k = 0; k < reactions.size(); ++k) reactions[k] = std::move(held[i].reactions.list[k]);
              held[i].reactions.list = std::move(reactions);
            }
          }
          d.systems = std::move(held);
        }

        d.outputToFiles = outputToFiles;
        if (devices.size() == d.systems.size())
        {
          for (std::size_t i = 0; i < devices.size(); ++i) d.systems[i].cbmcDevice = devices[i];
        }
        if (d.outputToFiles) d.createOutputFiles();
        d.createInterpolationGrids();
        d.copyOut();
      },
      nanobind::arg("file_name") = "restart_data.bin",
      "Reads the state of the driver from a binary restart file; 'run' then continues where that run stopped. A "
      "restart file is written every 'write_binary_restart_every' cycles when 'output_to_files' is set.");

  // a single cycle does not copy out: the phases and 'tear_down' do
  if constexpr (requires(Driver d) { d.performCycle(); })
  {
    driver.def("cycle", [](Bound& d) { d.performCycle(); });
  }

  driver.def_prop_ro(
      "number_of_blocks", [](const Bound& d) { return d.estimation.numberOfBins; },
      "The number of blocks the averages and their errors are computed over.");

  driver.def_prop_ro("systems",
                     [](nanobind::handle self)
                     {
                       // the live objects, each held by the parent: a copied vector lost every change made through it
                       nanobind::list list;
                       for (System& item : nanobind::cast<Bound&>(self).systems)
                         list.append(nanobind::cast(&item, nanobind::rv_policy::reference_internal, self));
                       return list;
                     });
}

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

    nanobind::enum_<EquationOfState::Type>(eos, "EquationOfStateType", nanobind::is_arithmetic())
        .value("PENG_ROBINSON", EquationOfState::Type::PengRobinson)
        .value("PENG_ROBINSON_GASEM", EquationOfState::Type::PengRobinsonGasem)
        .value("SOAVE_REDLICH_KWONG", EquationOfState::Type::SoaveRedlichKwong);

    nanobind::enum_<EquationOfState::MixingRules>(eos, "MixingRules", nanobind::is_arithmetic())
        .value("VAN_DER_WAALS", EquationOfState::MixingRules::VanDerWaals);

    nanobind::enum_<EquationOfState::FluidState>(eos, "FluidState", nanobind::is_arithmetic())
        .value("UNKNOWN", EquationOfState::FluidState::Unknown)
        .value("SUPER_CRITICAL_FLUID", EquationOfState::FluidState::SuperCriticalFluid)
        .value("VAPOR", EquationOfState::FluidState::Vapor)
        .value("LIQUID", EquationOfState::FluidState::Liquid)
        .value("VAPOR_LIQUID", EquationOfState::FluidState::VaporLiquid);

    nanobind::class_<EquationOfState::FluidResult>(eos, "FluidResult")
        .def(nanobind::init<double, std::optional<double>, EquationOfState::FluidState>(),
             nanobind::arg("compressibility"), nanobind::arg("fugacity_coefficient"), nanobind::arg("fluid_state"))
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
      .def_rw("charge", &Atom::charge, "The charge of the atom [e].")
      .def_prop_rw(
          "type", [](const Atom& a) { return static_cast<std::size_t>(a.type); },
          [](Atom& a, std::uint16_t type) { a.type = type; }, "The index of the pseudo-atom of the atom.")
      .def_prop_ro("molecule_id", [](const Atom& a) { return static_cast<std::size_t>(a.moleculeId); })
      .def_prop_ro("component_id", [](const Atom& a) { return static_cast<std::size_t>(a.componentId); })
      .def_ro("scaling_vdw", &Atom::scalingVDW, "The scaling of the Van der Waals interactions of the atom.")
      .def_ro("scaling_coulomb", &Atom::scalingCoulomb, "The scaling of the Coulomb interactions of the atom.")
      .def_prop_ro("is_fractional", [](const Atom& a) { return a.isFractional != 0; },
                   "Whether the atom is of a fractional molecule.")
      .def("__repr__", &Atom::repr);

  nanobind::class_<SimulationBox> simulationBox(m, "SimulationBox");
  nanobind::enum_<SimulationBox::Type>(simulationBox, "SimulationBoxType", nanobind::is_arithmetic())
      .value("RECTANGULAR", SimulationBox::Type::Rectangular)
      .value("TRICLINIC", SimulationBox::Type::Triclinic);

  simulationBox
      .def(nanobind::init<double, double, double>(), nanobind::arg("a"), nanobind::arg("b"), nanobind::arg("c"))
      // the angles in degrees (as documented, and as a CIF file gives them); the C++ constructors take radians
      .def(
          "__init__",
          [](SimulationBox* box, double a, double b, double c, double alpha, double beta, double gamma)
          {
            new (box) SimulationBox(a, b, c, alpha * Units::DegreesToRadians, beta * Units::DegreesToRadians,
                                    gamma * Units::DegreesToRadians);
          },
          nanobind::arg("a"), nanobind::arg("b"), nanobind::arg("c"), nanobind::arg("alpha"), nanobind::arg("beta"),
          nanobind::arg("gamma"))
      .def(
          "__init__",
          [](SimulationBox* box, double a, double b, double c, double alpha, double beta, double gamma,
             SimulationBox::Type type)
          {
            new (box) SimulationBox(a, b, c, alpha * Units::DegreesToRadians, beta * Units::DegreesToRadians,
                                    gamma * Units::DegreesToRadians, type);
          },
          nanobind::arg("a"), nanobind::arg("b"), nanobind::arg("c"), nanobind::arg("alpha"), nanobind::arg("beta"),
          nanobind::arg("gamma"), nanobind::arg("type"))
      .def(nanobind::init<double3x3, SimulationBox::Type>())
      .def_rw("type", &SimulationBox::type)
      .def_ro("length_a", &SimulationBox::lengthA)
      .def_ro("length_b", &SimulationBox::lengthB)
      .def_ro("length_c", &SimulationBox::lengthC)
      .def_ro("volume", &SimulationBox::volume);

  nanobind::class_<VDWParameters> vdwParameters(m, "VDWParameters");
  nanobind::enum_<VDWParameters::Type>(vdwParameters, "Type", nanobind::is_arithmetic())
      .value("NONE", VDWParameters::Type::None)
      .value("LENNARD_JONES", VDWParameters::Type::LennardJones)
      .value("BUCKINGHAM", VDWParameters::Type::BuckingHam)
      .value("MORSE", VDWParameters::Type::Morse)
      .value("FEYNMAN_HIBBS", VDWParameters::Type::FeynmannHibbs)
      .value("MM3", VDWParameters::Type::MM3)
      .value("BORN_HUGGINS_MEYER", VDWParameters::Type::BornHugginsMeyer)
      .value("LENNARD_JONES_SHIFTED_FORCE", VDWParameters::Type::LennardJonesShiftedForce)
      .value("POTENTIAL_12_6", VDWParameters::Type::Potential12_6)
      .value("POTENTIAL_12_6_2_0", VDWParameters::Type::Potential12_6_2_0)
      .value("CFF_9_6", VDWParameters::Type::CFF9_6)
      .value("CFF_EPSILON_SIGMA", VDWParameters::Type::CFFEpsilonSigma)
      .value("MATSUOKA_CLEMENTI_YOSHIMINE", VDWParameters::Type::MatsuokaClementiYoshimine)
      .value("GENERIC", VDWParameters::Type::Generic)
      .value("PELLENQ_NICHOLSON", VDWParameters::Type::PellenqNicholson)
      .value("HYDRATED_ION_WATER", VDWParameters::Type::HydratedIonWater)
      .value("MIE", VDWParameters::Type::Mie)
      .value("WEEKS_CHANDLER_ANDERSEN", VDWParameters::Type::WeeksChandlerAndersen)
      .value("LENNARD_JONES_SECOND_ORDER_TAYLOR_SHIFTED", VDWParameters::Type::LennardJonesSecondOrderTaylorShifted)
      .value("LENNARD_JONES_SWITCHED", VDWParameters::Type::LennardJonesSwitched)
      .value("LENNARD_JONES_FORCE_SWITCHED", VDWParameters::Type::LennardJonesForceSwitched)
      .value("REPULSIVE_HARMONIC", VDWParameters::Type::RepulsiveHarmonic);
  vdwParameters
      .def(nanobind::init<double, double>(), nanobind::arg("epsilon"), nanobind::arg("sigma"),
           "A Lennard-Jones potential: epsilon [K] and sigma [Angstrom].")
      .def(
          "__init__",
          [](VDWParameters* parameters, VDWParameters::Type type, const std::vector<double>& values)
          {
            try
            {
              new (parameters) VDWParameters(type, values);
            }
            catch (const std::runtime_error& error)
            {
              throw std::invalid_argument(std::format("VDWParameters: {}", error.what()));
            }
          },
          nanobind::arg("type"), nanobind::arg("parameters"),
          "A potential of the given type, with its parameters as in a force-field file (energies in Kelvin, lengths "
          "in Angstrom).")
      .def_ro("type", &VDWParameters::type)
      .def("__repr__", [](const VDWParameters& parameters)
           { return std::format("VDWParameters({})", VDWParameters::nameOfType(parameters.type)); });

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
  forceField.def_static(
      "from_file",
      [](const std::string& file_name)
      {
        try
        {
          return ForceField(file_name);
        }
        catch (const nlohmann::json::exception& error)
        {
          throw std::invalid_argument(std::format("ForceField: file '{}': {}", file_name, error.what()));
        }
        catch (const std::runtime_error& error)
        {
          throw std::invalid_argument(error.what());
        }
      },
      nanobind::arg("file_name"),
      "Reads a force field from a file (JSON), as the command-line program reads 'force_field.json'. What "
      "'simulation.json' adds to it (the cutoffs, the charge method) is set on the object.");
  nanobind::enum_<ForceField::MixingRule>(forceField, "MixingRule", nanobind::is_arithmetic())
      .value("LORENTZ_BERTHELOT", ForceField::MixingRule::Lorentz_Berthelot)
      .value("JORGENSEN", ForceField::MixingRule::Jorgensen)
      .value("SIXTH_POWER", ForceField::MixingRule::SixthPower);
  nanobind::enum_<ForceField::ChargeMethod>(forceField, "ChargeMethod", nanobind::is_arithmetic())
      .value("EWALD", ForceField::ChargeMethod::Ewald)
      .value("COULOMB", ForceField::ChargeMethod::Coulomb)
      .value("WOLF", ForceField::ChargeMethod::Wolf)
      .value("DAMPED_SHIFTED_FORCE", ForceField::ChargeMethod::DampedShiftedForce)
      .value("MODIFIED_SHIFTED_FORCE", ForceField::ChargeMethod::ModifiedShiftedForce)
      .value("ZERO_DIPOLE", ForceField::ChargeMethod::ZeroDipole);

  using BinaryInteraction = std::tuple<std::string, std::string, VDWParameters>;
  forceField
      .def(
          "__init__",
          [](ForceField* forceField, std::vector<PseudoAtom> pseudo_atoms, std::vector<VDWParameters> parameters,
             ForceField::MixingRule mixing_rule, double cutoff_framework_vdw, double cutoff_molecule_vdw,
             double cutoff_coulomb, bool shifted, bool tail_corrections, bool use_charge,
             const std::vector<BinaryInteraction>& binary_interactions,
             const std::vector<std::optional<VDWParameters>>& parameters_14,
             const std::vector<BinaryInteraction>& binary_interactions_14)
          {
            if (parameters.size() != pseudo_atoms.size())
            {
              throw std::invalid_argument(
                  std::format("ForceField: {} pseudo-atoms, but {} parameters (one for every pseudo-atom, in the same "
                              "order)",
                              pseudo_atoms.size(), parameters.size()));
            }
            if (!parameters_14.empty() && parameters_14.size() != pseudo_atoms.size())
            {
              throw std::invalid_argument(
                  std::format("ForceField: {} pseudo-atoms, but {} entries in parameters_14 (one for every "
                              "pseudo-atom, None where the 1-4 pairs have the regular parameters)",
                              pseudo_atoms.size(), parameters_14.size()));
            }
            // built aside and moved into the Python object at the end: an invalid binary interaction throws out of
            // here, and a force field already constructed in place would never be destroyed
            ForceField built(pseudo_atoms, parameters, mixing_rule, cutoff_framework_vdw, cutoff_molecule_vdw,
                             cutoff_coulomb, shifted, tail_corrections, use_charge);
            ForceField* const self = forceField;
            {
              ForceField* const forceField = &built;

            const bool hasParameters14 =
                std::ranges::any_of(parameters_14, [](const auto& given) { return given.has_value(); });
            if (binary_interactions.empty() && !hasParameters14 && binary_interactions_14.empty())
            {
              new (self) ForceField(std::move(built));
              return;
            }
            // a pair of a zero epsilon or sigma has no interaction, as in the 'BinaryInteractions' of a file (a
            // Lennard-Jones pair of a zero sigma is 1 / sigma^2 in the kernels)
            auto normalized = [](VDWParameters pair)
            {
              if (VDWParameters::isLennardJonesForm(pair.type) && (pair.parameters.x == 0.0 || pair.parameters.y == 0.0))
              {
                pair.type = VDWParameters::Type::None;
              }
              return pair;
            };

            const std::size_t numberOfPseudoAtoms = forceField->numberOfPseudoAtoms;
            auto pairOf = [&](const std::string& nameA, const std::string& nameB)
            {
              const std::optional<std::size_t> indexA = forceField->findPseudoAtom(nameA);
              const std::optional<std::size_t> indexB = forceField->findPseudoAtom(nameB);
              if (!indexA.has_value() || !indexB.has_value())
              {
                throw std::invalid_argument(std::format("ForceField: binary interaction of unknown pseudo-atom '{}'",
                                                        indexA.has_value() ? nameB : nameA));
              }
              return std::pair{*indexA, *indexB};
            };

            // The separate parameters of the 1-4 pairs, as 'parameters14' of a force-field file: the 1-4 table starts
            // as the regular one, takes the parameters that are given, and is mixed like the regular one.
            if (hasParameters14)
            {
              for (std::size_t i = 0; i < parameters_14.size(); ++i)
              {
                if (parameters_14[i].has_value()) forceField->setPair14SelfInteraction(i, *parameters_14[i]);
              }
              forceField->mixTable(forceField->data14);
            }

            // A pair that is given replaces what the mixing rule made of it, as 'BinaryInteractions' of a force-field
            // file does. It is the 1-4 pair as well, unless that pair is given or one of its pseudo-atoms has 1-4
            // parameters of its own (then the mixed 1-4 pair was meant).
            for (const auto& [nameA, nameB, givenPair] : binary_interactions)
            {
              const VDWParameters pair = normalized(givenPair);
              const auto [indexA, indexB] = pairOf(nameA, nameB);
              forceField->data[indexA * numberOfPseudoAtoms + indexB] = pair;
              forceField->data[indexB * numberOfPseudoAtoms + indexA] = pair;
              const bool given14 = std::ranges::any_of(binary_interactions_14, [&](const BinaryInteraction& other)
                                                       {
                                                         const auto [a, b] = pairOf(std::get<0>(other), std::get<1>(other));
                                                         return (a == indexA && b == indexB) || (a == indexB && b == indexA);
                                                       });
              const bool explicit14 = hasParameters14 && (parameters_14[indexA].has_value() || parameters_14[indexB].has_value());
              if (!forceField->data14.empty() && !given14 && !explicit14)
              {
                forceField->data14[indexA * numberOfPseudoAtoms + indexB] = pair;
                forceField->data14[indexB * numberOfPseudoAtoms + indexA] = pair;
              }
            }
            for (const auto& [nameA, nameB, givenPair] : binary_interactions_14)
            {
              const VDWParameters pair = normalized(givenPair);
              const auto [indexA, indexB] = pairOf(nameA, nameB);
              forceField->ensurePair14Table();
              forceField->data14[indexA * numberOfPseudoAtoms + indexB] = pair;
              forceField->data14[indexB * numberOfPseudoAtoms + indexA] = pair;
            }

            // the shifts and the tail corrections are those of the pairs as they are now
            forceField->preComputeDerivedParameters();
            forceField->preComputePotentialShift();
            forceField->preComputeTailCorrection();
            }
            new (self) ForceField(std::move(built));
          },
          nanobind::arg("pseudo_atoms"), nanobind::arg("parameters"),
          nanobind::arg("mixing_rule") = ForceField::MixingRule::Lorentz_Berthelot,
          nanobind::arg("cutoff_framework_vdw") = 12.0, nanobind::arg("cutoff_molecule_vdw") = 12.0,
          nanobind::arg("cutoff_coulomb") = 12.0, nanobind::arg("shifted") = true,
          nanobind::arg("tail_corrections") = false, nanobind::arg("use_charge") = true,
          nanobind::arg("binary_interactions") = std::vector<BinaryInteraction>{},
          nanobind::arg("parameters_14") = std::vector<std::optional<VDWParameters>>{},
          nanobind::arg("binary_interactions_14") = std::vector<BinaryInteraction>{})
      .def("find_pseudo_atom", static_cast<std::optional<std::size_t> (ForceField::*)(const std::string &) const>(&ForceField::findPseudoAtom))
      .def("__repr__", &ForceField::repr)
      .def_ro("pseudo_atoms", &ForceField::pseudoAtoms)
      .def_ro("vdw_parameters", &ForceField::data)
      .def_rw("use_charge", &ForceField::useCharge)
      // The settings below are those of the force field a component or a system is constructed with: both keep
      // their own copy, so a setting is made before the force field is handed to them.
      .def_rw("charge_method", &ForceField::chargeMethod, "How the electrostatic interactions are computed.")
      .def_prop_rw(
          "ewald_precision", [](const ForceField& f) { return f.EwaldPrecision; },
          [](ForceField& f, double precision)
          {
            if (!(precision > 0.0)) throw std::invalid_argument("ewald_precision must be positive");
            f.EwaldPrecision = precision;
            f.automaticEwald = true;
          },
          "The relative precision the Ewald parameters are chosen for [-]; setting it makes the choice automatic.")
      .def_ro("ewald_alpha", &ForceField::EwaldAlpha, "The Ewald damping parameter alpha [1/Angstrom].")
      .def_ro("number_of_wave_vectors", &ForceField::numberOfWaveVectors,
              "The number of wave vectors of the Ewald sum in each direction.")
      .def_ro("automatic_ewald", &ForceField::automaticEwald,
              "Whether alpha and the wave vectors are chosen from 'ewald_precision'.")
      .def(
          "set_ewald_parameters",
          [](ForceField& f, double alpha, int3 numberOfWaveVectors)
          {
            if (!(alpha > 0.0)) throw std::invalid_argument("set_ewald_parameters: alpha must be positive");
            if (numberOfWaveVectors.x < 0 || numberOfWaveVectors.y < 0 || numberOfWaveVectors.z < 0)
            {
              throw std::invalid_argument("set_ewald_parameters: the number of wave vectors can not be negative");
            }
            f.EwaldAlpha = alpha;
            f.numberOfWaveVectors = numberOfWaveVectors;
            f.automaticEwald = false;
          },
          nanobind::arg("alpha"), nanobind::arg("number_of_wave_vectors"),
          "Fixes the Ewald damping parameter alpha [1/Angstrom] and the number of wave vectors, instead of "
          "choosing them from 'ewald_precision'.")
      .def_prop_rw(
          "modified_shifted_force_beta", [](const ForceField& f) { return f.modifiedShiftedForceBeta; },
          [](ForceField& f, double beta)
          {
            if (!(beta > 0.0)) throw std::invalid_argument("modified_shifted_force_beta must be positive");
            f.modifiedShiftedForceBeta = beta;
          },
          "The switching parameter beta of the modified shifted-force charge method [1/Angstrom].")
      .def_rw("omit_ewald_fourier", &ForceField::omitEwaldFourier,
              "Leaves out the Fourier part of the Ewald sum.")
      .def_rw("compute_polarization", &ForceField::computePolarization, "Whether the polarization energy is computed.")
      .def_rw("omit_inter_polarization", &ForceField::omitInterPolarization,
              "Leaves out the polarization between molecules.")
      .def_rw("omit_inter_interactions", &ForceField::omitInterInteractions,
              "Leaves out the interactions between molecules.")
      // a cutoff that is set is a fixed one, and the shifts and tail corrections that depend on it are made again
      .def_prop_rw(
          "cutoff_framework_vdw", [](const ForceField& f) { return f.cutOffFrameworkVDW; },
          [](ForceField& f, double cutoff)
          {
            if (!(cutoff > 0.0)) throw std::invalid_argument("cutoff_framework_vdw must be positive");
            f.cutOffFrameworkVDW = cutoff;
            f.cutOffFrameworkVDWAutomatic = false;
            f.preComputeCutOffDependentParameters();
          },
          "The cutoff of the van der Waals interactions with the framework [Å].")
      .def_prop_rw(
          "cutoff_molecule_vdw", [](const ForceField& f) { return f.cutOffMoleculeVDW; },
          [](ForceField& f, double cutoff)
          {
            if (!(cutoff > 0.0)) throw std::invalid_argument("cutoff_molecule_vdw must be positive");
            f.cutOffMoleculeVDW = cutoff;
            f.cutOffMoleculeVDWAutomatic = false;
            f.preComputeCutOffDependentParameters();
          },
          "The cutoff of the van der Waals interactions between molecules [Å].")
      .def_prop_rw(
          "cutoff_coulomb", [](const ForceField& f) { return f.cutOffCoulomb; },
          [](ForceField& f, double cutoff)
          {
            if (!(cutoff > 0.0)) throw std::invalid_argument("cutoff_coulomb must be positive");
            f.cutOffCoulomb = cutoff;
            f.cutOffCoulombAutomatic = false;
            f.preComputeCutOffDependentParameters();
          },
          "The cutoff of the real-space electrostatics [Å].");

  nanobind::class_<CIFReader> cifReader(m, "CIFReader");
  cifReader
      .def_static("expand_defined_atoms_to_unit_cell", &CIFReader::expandDefinedAtomsToUnitCell,
           nanobind::arg("simulation_box"), nanobind::arg("space_group_hall_number"),
           nanobind::arg("defined_atoms"),
           nanobind::arg("symmetry_operation_strings") = std::vector<std::string>{});

  nanobind::enum_<FrameworkGroupType>(m, "FrameworkGroupType", nanobind::is_arithmetic())
      .value("FIXED", FrameworkGroupType::Fixed)
      .value("RIGID", FrameworkGroupType::Rigid)
      .value("FLEXIBLE", FrameworkGroupType::Flexible);

  using BondDefinition = std::tuple<std::array<std::string, 2>, BondType, std::vector<double>>;
  using BendDefinition = std::tuple<std::array<std::string, 3>, BendType, std::vector<double>>;
  using TorsionDefinition = std::tuple<std::array<std::string, 4>, TorsionType, std::vector<double>>;
  using GroupDefinition = std::pair<FrameworkGroupType, std::vector<std::size_t>>;

  nanobind::class_<FrameworkDefinition>(m, "FrameworkDefinition")
      .def(
          "__init__",
          [](FrameworkDefinition* definition, bool flexible, const std::vector<BondDefinition>& bonds,
             const std::vector<BendDefinition>& bends, const std::vector<TorsionDefinition>& torsions,
             const std::vector<TorsionDefinition>& improper_torsions, const std::vector<GroupDefinition>& groups,
             bool exclude_intra_12_interactions, bool exclude_intra_13_interactions,
             bool exclude_intra_bond_interactions, bool exclude_intra_bend_interactions,
             double intra_14_van_der_waals_scaling, double intra_14_charge_charge_scaling)
          {
            using Description = nlohmann::basic_json<nlohmann::raspa_map>;
            auto potentials = [](const auto& definitions, const auto& names)
            {
              Description list = Description::array();
              for (const auto& [atomTypes, type, parameters] : definitions)
              {
                Description item = Description::array();
                item.push_back(Description(atomTypes));
                item.push_back(Description(nameOfType(names, type)));
                item.push_back(Description(parameters));
                list.push_back(item);
              }
              return list;
            };

            new (definition) FrameworkDefinition();
            Description& data = definition->data;
            data = Description::object();
            data["Type"] = flexible ? "Flexible" : "Rigid";
            if (!bonds.empty()) data["Bonds"] = potentials(bonds, bondTypeNames);
            if (!bends.empty()) data["Bends"] = potentials(bends, bendTypeNames);
            if (!torsions.empty()) data["Torsions"] = potentials(torsions, torsionTypeNames);
            if (!improper_torsions.empty()) data["ImproperTorsions"] = potentials(improper_torsions, torsionTypeNames);
            if (!groups.empty())
            {
              Description list = Description::array();
              for (const auto& [type, atoms] : groups)
              {
                Description group = Description::object();
                group["Type"] = type == FrameworkGroupType::Fixed   ? "Fixed"
                                : type == FrameworkGroupType::Rigid ? "Rigid"
                                                                    : "Flexible";
                group["Atoms"] = atoms;
                list.push_back(group);
              }
              data["Groups"] = list;
            }
            data["ExcludeIntra12Interactions"] = exclude_intra_12_interactions;
            data["ExcludeIntra13Interactions"] = exclude_intra_13_interactions;
            data["ExcludeIntraBondInteractions"] = exclude_intra_bond_interactions;
            data["ExcludeIntraBendInteractions"] = exclude_intra_bend_interactions;
            data["Intra14VanDerWaalsScalingValue"] = intra_14_van_der_waals_scaling;
            data["Intra14ChargeChargeScalingValue"] = intra_14_charge_charge_scaling;
          },
          nanobind::arg("flexible") = true, nanobind::arg("bonds") = std::vector<BondDefinition>{},
          nanobind::arg("bends") = std::vector<BendDefinition>{},
          nanobind::arg("torsions") = std::vector<TorsionDefinition>{},
          nanobind::arg("improper_torsions") = std::vector<TorsionDefinition>{},
          nanobind::arg("groups") = std::vector<GroupDefinition>{},
          nanobind::arg("exclude_intra_12_interactions") = true,
          nanobind::arg("exclude_intra_13_interactions") = true,
          nanobind::arg("exclude_intra_bond_interactions") = false,
          nanobind::arg("exclude_intra_bend_interactions") = false,
          nanobind::arg("intra_14_van_der_waals_scaling") = 0.0,
          nanobind::arg("intra_14_charge_charge_scaling") = 0.0,
          "The definition of a flexible framework. A potential is the names of the pseudo-atoms it is between, its "
          "type and the parameters of that type; it applies to every such bond, bend or torsion of the framework. A "
          "group is a type (fixed, rigid or flexible) and the atoms of the framework it holds.")
      .def("__repr__", [](const FrameworkDefinition& definition) { return definition.data.dump(2); });

  nanobind::class_<Framework> framework(m, "Framework");

  nanobind::enum_<CIFReader::UseChargesFrom>(framework, "UseChargesFrom", nanobind::is_arithmetic())
      .value("PSEUDO_ATOMS", CIFReader::UseChargesFrom::PseudoAtoms)
      .value("CIF_FILE", CIFReader::UseChargesFrom::CIF_File)
      .value("CHARGE_EQUILIBRATION", CIFReader::UseChargesFrom::ChargeEquilibration);

  framework
      .def(
          "__init__",
          [](Framework* framework, const ForceField& force_field, std::string component_name,
             SimulationBox simulation_box, std::size_t space_group_hall_number, const std::vector<Atom>& defined_atoms,
             int3 number_of_unit_cells, const std::optional<FrameworkDefinition>& definition,
             const std::optional<std::string>& definition_file)
          {
            if (definition.has_value() && definition_file.has_value())
            {
              throw std::invalid_argument("Framework: give either a definition or a definition_file, not both");
            }
            new (framework) Framework(force_field, component_name, simulation_box, space_group_hall_number,
                                      defined_atoms, number_of_unit_cells);
            try
            {
              if (definition.has_value()) framework->applyFrameworkDefinition(force_field, definition->data);
              if (definition_file.has_value()) framework->readFrameworkDefinition(force_field, *definition_file);
            }
            catch (const std::runtime_error& error)
            {
              framework->~Framework();
              throw std::invalid_argument(std::format("Framework: {}", error.what()));
            }
          },
          nanobind::arg("force_field"), nanobind::arg("component_name"), nanobind::arg("simulation_box"),
          nanobind::arg("space_group_hall_number"), nanobind::arg("defined_atoms"),
          nanobind::arg("number_of_unit_cells"), nanobind::arg("definition") = nanobind::none(),
          nanobind::arg("definition_file") = nanobind::none())
      .def_static(
          "from_cif",
          [](const ForceField& force_field, const std::string& file_name, std::optional<int3> number_of_unit_cells,
             CIFReader::UseChargesFrom use_charges_from, std::optional<std::string> component_name,
             const std::optional<FrameworkDefinition>& definition, const std::optional<std::string>& definition_file)
          {
            std::string content;
            try
            {
              content = readFileContent(file_name, ".cif");
            }
            catch (const std::runtime_error& error)
            {
              throw std::invalid_argument(std::format("Framework: {}", error.what()));
            }
            const std::string name = component_name.value_or(std::filesystem::path(file_name).stem().string());
            return frameworkFromCIF(force_field, content, name, number_of_unit_cells, use_charges_from, definition,
                                    definition_file);
          },
          nanobind::arg("force_field"), nanobind::arg("file_name"),
          nanobind::arg("number_of_unit_cells") = nanobind::none(),
          nanobind::arg("use_charges_from") = CIFReader::UseChargesFrom::PseudoAtoms,
          nanobind::arg("component_name") = nanobind::none(), nanobind::arg("definition") = nanobind::none(),
          nanobind::arg("definition_file") = nanobind::none(),
          "A framework from a CIF file ('.cif' may be left out; the file is also looked for in RASPA_DIR). Without "
          "a number of unit cells the smallest supercell is taken in which the cutoffs of the force field fit. The "
          "name is that of the file unless one is given.")
      .def_static(
          "from_cif_string",
          [](const ForceField& force_field, const std::string& cif, const std::string& component_name,
             std::optional<int3> number_of_unit_cells, CIFReader::UseChargesFrom use_charges_from,
             const std::optional<FrameworkDefinition>& definition, const std::optional<std::string>& definition_file)
          {
            return frameworkFromCIF(force_field, cif, component_name, number_of_unit_cells, use_charges_from,
                                    definition, definition_file);
          },
          nanobind::arg("force_field"), nanobind::arg("cif"), nanobind::arg("component_name"),
          nanobind::arg("number_of_unit_cells") = nanobind::none(),
          nanobind::arg("use_charges_from") = CIFReader::UseChargesFrom::PseudoAtoms,
          nanobind::arg("definition") = nanobind::none(), nanobind::arg("definition_file") = nanobind::none(),
          "A framework from the content of a CIF file, as 'from_cif' makes it from a file.")
      .def_ro("name", &Framework::name)
      .def_ro("rigid", &Framework::rigid, "Whether no atom of the framework moves.")
      .def_ro("simulation_box", &Framework::simulationBox, "The unit cell.")
      .def_ro("space_group_hall_number", &Framework::spaceGroupHallNumber)
      .def_ro("number_of_unit_cells", &Framework::numberOfUnitCells)
      .def_ro("mass", &Framework::mass, "The mass of the framework [amu].")
      .def_ro("unit_cell_mass", &Framework::unitCellMass, "The mass of one unit cell [amu].")
      .def_ro("net_charge", &Framework::netCharge, "The net charge of the framework [e].")
      .def_ro("defined_atoms", &Framework::definedAtoms, "The atoms of the asymmetric unit, in fractional positions.")
      .def_ro("unit_cell_atoms", &Framework::unitCellAtoms, "The atoms of one unit cell, in Cartesian positions.")
      .def_ro("atoms", &Framework::atoms, "The atoms of the supercell, in Cartesian positions.")
      .def_prop_ro("number_of_atoms", [](const Framework& f) { return f.atoms.size(); })
      .def_prop_ro("number_of_fixed_atoms", &Framework::numberOfFixedAtoms)
      .def_prop_ro("number_of_mobile_atoms", &Framework::numberOfMobileAtoms)
      .def_prop_ro("number_of_rigid_groups", &Framework::numberOfRigidGroups)
      .def_ro("connectivity_table", &Framework::connectivityTable, "The bonds of a flexible framework.")
      .def_ro("intra_molecular_potentials", &Framework::intraMolecularPotentials,
              "The bonded terms of a flexible framework, one for every bond, bend and torsion.")
      .def("print", &Framework::printStatus)
      .def("__repr__", &Framework::repr);

  // the zeolites that are built in: the force field needs the pseudo-atoms 'Si' and 'O'
  for (const auto& [name, maker] :
       std::initializer_list<std::pair<const char*, Framework (*)(const ForceField&, int3)>>{
           {"make_mfi", &Framework::makeMFI},
           {"make_itq_29", &Framework::makeITQ29},
           {"make_fau", &Framework::makeFAU},
           {"make_cha", &Framework::makeCHA}})
  {
    framework.def_static(
        name,
        [maker](const ForceField& force_field, int3 number_of_unit_cells)
        {
          try
          {
            return maker(force_field, number_of_unit_cells);
          }
          catch (const std::runtime_error& error)
          {
            throw std::invalid_argument(std::format("Framework: {}", error.what()));
          }
        },
        nanobind::arg("force_field"), nanobind::arg("number_of_unit_cells") = int3{1, 1, 1});
  }

  // one keyword for every move of Move.Types. The constructor took twenty positional numbers before, of which the
  // last five were handed to the wrong moves ('widom_probability' set the Gibbs identity change, and so on).
  static const std::string moveProbabilitiesSignatureText = moveProbabilitiesSignature();
  nanobind::class_<MCMoveProbabilities>(m, "MCMoveProbabilities")
      .def(
          "__init__",
          [](MCMoveProbabilities* probabilities, nanobind::kwargs keywords)
          {
            // built aside: an unknown keyword or a value that is not a number throws out of here, and an object
            // already constructed in place would never be destroyed
            MCMoveProbabilities built{};
            for (auto [key, value] : keywords)
            {
              const Move::Types move = moveTypeOfProbabilityKeyword(nanobind::cast<std::string>(key));
              if (!nanobind::isinstance<nanobind::float_>(value) && !nanobind::isinstance<nanobind::int_>(value))
              {
                throw nanobind::type_error("MCMoveProbabilities: a probability is a number");
              }
              setMoveProbability(built, move, nanobind::cast<double>(value));
            }
            new (probabilities) MCMoveProbabilities(std::move(built));
          },
          nanobind::sig(moveProbabilitiesSignatureText.c_str()),
          "The relative probabilities of the Monte Carlo moves, one keyword for every move of Move.Types.")
      .def("__getitem__",
           [](const MCMoveProbabilities& probabilities, Move::Types move) { return probabilities.getProbability(move); })
      .def("__setitem__", &setMoveProbability)
      .def("get_probability",
           [](const MCMoveProbabilities& probabilities, Move::Types move) { return probabilities.getProbability(move); },
           nanobind::arg("move"))
      .def("set_probability", &setMoveProbability, nanobind::arg("move"), nanobind::arg("probability"))
      .def("join", &MCMoveProbabilities::join)
      .def("__repr__", &MCMoveProbabilities::repr);

  nanobind::class_<PropertyLambdaProbabilityHistogram>(m, "PropertyLambdaProbabilityHistogram")
      .def(nanobind::init<>())
      .def(nanobind::init<std::size_t, std::size_t>(), nanobind::arg("number_of_blocks"),
           nanobind::arg("number_of_sample_points"))
      // dimensionless
      .def_prop_rw(
          "bias_factor", [](const PropertyLambdaProbabilityHistogram& p) { return p.biasFactor; },
          [](PropertyLambdaProbabilityHistogram& p, const std::vector<double>& value)
          {
            // one bias per lambda bin: the moves index it by bin
            if (value.size() != p.numberOfSamplePoints)
              throw std::invalid_argument(std::format("bias_factor needs {} values (one per lambda bin), got {}",
                                                      p.numberOfSamplePoints, value.size()));
            p.biasFactor = value;
          })
      // dimensionless
      .def_prop_rw(
          "histogram", [](const PropertyLambdaProbabilityHistogram& p) { return p.histogram; },
          [](PropertyLambdaProbabilityHistogram& p, const std::vector<double>& value)
          {
            if (value.size() != p.numberOfSamplePoints)
              throw std::invalid_argument(std::format("histogram needs {} values (one per lambda bin), got {}",
                                                      p.numberOfSamplePoints, value.size()));
            p.histogram = value;
          })
      .def_ro("occupancy_total", &PropertyLambdaProbabilityHistogram::occupancyTotal)
      .def_ro("occupancy_count", &PropertyLambdaProbabilityHistogram::occupancyCount)
      // convert result to units of Kelvin
      .def("average_dudlambda",
           [](PropertyLambdaProbabilityHistogram& p) { return Units::EnergyToKelvin * p.averageDuDlambda(); })
      // dimensionless
      .def("result", &PropertyLambdaProbabilityHistogram::result);

  nanobind::class_<ConnectivityTable>(m, "ConnectivityTable")
      .def(nanobind::init<>())
      .def(
          "__init__",
          [](ConnectivityTable* table, std::size_t number_of_beads,
             const std::vector<std::array<std::size_t, 2>>& bonds)
          {
            for (const std::array<std::size_t, 2>& bond : bonds)
            {
              if (bond[0] >= number_of_beads || bond[1] >= number_of_beads || bond[0] == bond[1])
              {
                throw std::invalid_argument(std::format("ConnectivityTable: {}-{} is not a bond between two of the {} "
                                                        "beads",
                                                        bond[0], bond[1], number_of_beads));
              }
            }
            new (table) ConnectivityTable(number_of_beads);
            for (const std::array<std::size_t, 2>& bond : bonds)
            {
              (*table)[bond[0], bond[1]] = true;
              (*table)[bond[1], bond[0]] = true;
            }
          },
          nanobind::arg("number_of_beads"), nanobind::arg("bonds") = std::vector<std::array<std::size_t, 2>>{})
      .def_ro("number_of_beads", &ConnectivityTable::numberOfBeads)
      .def_prop_ro("bonds", &ConnectivityTable::findAllBonds, "The pairs of connected beads.")
      .def_prop_ro("bends", &ConnectivityTable::findAllBends, "The triplets of beads that form a bend.")
      .def_prop_ro("torsions", &ConnectivityTable::findAllTorsions, "The quadruplets of beads that form a torsion.")
      .def("__repr__", &ConnectivityTable::printConnectivityTable);

  // Two terms of the molecule files are left out: the out-of-plane bend, which has no energy yet (it evaluates to
  // zero), and the in-plane MM3 bend, which the bend potential refuses.
  bindIntraMolecularPotential<BondPotential, BondType>(m, "BondPotential", "BondType",
      bondTypeNames,
      BondPotential::numberOfBondParameters);
  bindIntraMolecularPotential<UreyBradleyPotential, UreyBradleyType>(m, "UreyBradleyPotential", "UreyBradleyType",
      {{"FIXED", UreyBradleyType::Fixed}, {"HARMONIC", UreyBradleyType::Harmonic}, {"CORE_SHELL_SPRING", UreyBradleyType::CoreShellSpring}, {"MORSE", UreyBradleyType::Morse}, {"LJ_12_6", UreyBradleyType::LJ_12_6}, {"LENNARD_JONES", UreyBradleyType::LennardJones}, {"BUCKINGHAM", UreyBradleyType::Buckingham}, {"RESTRAINED_HARMONIC", UreyBradleyType::RestrainedHarmonic}, {"QUARTIC", UreyBradleyType::Quartic}, {"CFF_QUARTIC", UreyBradleyType::CFF_Quartic}, {"MM3", UreyBradleyType::MM3}},
      UreyBradleyPotential::numberOfUreyBradleyParameters);
  bindIntraMolecularPotential<BendPotential, BendType>(m, "BendPotential", "BendType",
      bendTypeNames,
      BendPotential::numberOfBendParameters);
  bindIntraMolecularPotential<InversionBendPotential, InversionBendType>(m, "InversionBendPotential", "InversionBendType",
      {{"HARMONIC", InversionBendType::Harmonic}, {"HARMONIC_COSINE", InversionBendType::HarmonicCosine}, {"PLANAR", InversionBendType::Planar}, {"HARMONIC2", InversionBendType::Harmonic2}, {"HARMONIC_COSINE2", InversionBendType::HarmonicCosine2}, {"PLANAR2", InversionBendType::Planar2}, {"MM3", InversionBendType::MM3}},
      InversionBendPotential::numberOfInversionBendParameters);
  bindIntraMolecularPotential<TorsionPotential, TorsionType>(m, "TorsionPotential", "TorsionType",
      torsionTypeNames,
      TorsionPotential::numberOfTorsionParameters);
  bindIntraMolecularPotential<BondBondPotential, BondBondType>(m, "BondBondPotential", "BondBondType",
      {{"CVFF", BondBondType::CVFF}, {"CFF", BondBondType::CFF}},
      BondBondPotential::numberOfBondBondParameters);
  bindIntraMolecularPotential<BondBendPotential, BondBendType>(m, "BondBendPotential", "BondBendType",
      {{"CVFF", BondBendType::CVFF}, {"CFF", BondBendType::CFF}, {"MM3", BondBendType::MM3}, {"TRUNCATED_HARMONIC", BondBendType::TruncatedHarmonic}, {"SCREENED_HARMONIC", BondBendType::ScreenedHarmonic}, {"SCREENED_VESSAL", BondBendType::ScreenedVessal}, {"TRUNCATED_VESSAL", BondBendType::TruncatedVessal}},
      BondBendPotential::numberOfBondBendParameters);
  bindIntraMolecularPotential<BondTorsionPotential, BondTorsionType>(m, "BondTorsionPotential", "BondTorsionType",
      {{"MM3", BondTorsionType::MM3}},
      BondTorsionPotential::numberOfBondTorsionParameters);
  bindIntraMolecularPotential<BendBendPotential, BendBendType>(m, "BendBendPotential", "BendBendType",
      {{"CVFF", BendBendType::CVFF}, {"CFF", BendBendType::CFF}, {"MM3", BendBendType::MM3}},
      BendBendPotential::numberOfBendBendParameters);
  bindIntraMolecularPotential<BendTorsionPotential, BendTorsionType>(m, "BendTorsionPotential", "BendTorsionType",
      {{"SMOOTHED", BendTorsionType::Smoothed}, {"SMOOTHED_THREE_COSINE", BendTorsionType::SmoothedThreeCosine}, {"NICHOLAS", BendTorsionType::Nicholas}, {"CFF", BendTorsionType::CFF}, {"SMOOTHED_CFF", BendTorsionType::SmoothedCFF}, {"SMOOTHED_CFF2", BendTorsionType::SmoothedCFF2}, {"SMOOTHED_CFF3", BendTorsionType::SmoothedCFF3}, {"CVFF", BendTorsionType::CVFF}},
      BendTorsionPotential::numberOfBendTorsionParameters);

  // the pairs of atoms inside a molecule whose non-bonded interaction is given, instead of following from the bonds
  nanobind::class_<VanDerWaalsPotential>(m, "VanDerWaalsPotential")
      .def(nanobind::init<std::array<std::size_t, 2>, const VDWParameters&, double>(), nanobind::arg("identifiers"),
           nanobind::arg("parameters"), nanobind::arg("scaling") = 1.0)
      .def_ro("identifiers", &VanDerWaalsPotential::identifiers)
      .def_ro("scaling", &VanDerWaalsPotential::scaling);

  nanobind::class_<CoulombPotential>(m, "CoulombPotential")
      .def(
          "__init__",
          [](CoulombPotential* potential, std::array<std::size_t, 2> identifiers, double charge_a, double charge_b,
             double scaling)
          { new (potential) CoulombPotential(identifiers, CoulombType::Coulomb, charge_a, charge_b, scaling); },
          nanobind::arg("identifiers"), nanobind::arg("charge_a"), nanobind::arg("charge_b"),
          nanobind::arg("scaling") = 1.0)
      .def_ro("identifiers", &CoulombPotential::identifiers)
      .def_ro("scaling", &CoulombPotential::scaling);

  // The CMAP correction of protein force fields: an energy map over two consecutive backbone dihedrals, and the
  // terms that use a map.
  nanobind::class_<CMAPMap>(m, "CMAPMap")
      .def(
          "__init__",
          [](CMAPMap* map, std::string name, std::size_t resolution, std::vector<double> energies)
          {
            for (double& energy : energies) energy *= Units::KelvinToEnergy;
            try
            {
              new (map) CMAPMap(std::move(name), resolution, std::move(energies));
            }
            catch (const std::runtime_error& error)
            {
              throw std::invalid_argument(error.what());
            }
          },
          nanobind::arg("name"), nanobind::arg("resolution"), nanobind::arg("energies"),
          "A map of resolution x resolution energies [K] over the two angles; the first angle is the slow index.")
      .def_ro("name", &CMAPMap::name)
      .def_ro("resolution", &CMAPMap::resolution)
      .def("__repr__", &CMAPMap::print);

  nanobind::class_<CMAPPotential>(m, "CMAPPotential")
      .def(nanobind::init<std::array<std::size_t, 5>, std::size_t>(), nanobind::arg("identifiers"),
           nanobind::arg("map_index"),
           "A CMAP term of the five atoms A, B, C, D, E with the map at 'map_index' of the maps of the potentials.")
      .def_ro("identifiers", &CMAPPotential::identifiers)
      .def_ro("map_index", &CMAPPotential::mapIndex)
      .def("__repr__", &CMAPPotential::print);

  nanobind::class_<ChiralCenter> chiralCenter(m, "ChiralCenter");
  nanobind::enum_<ChiralCenter::Chirality>(chiralCenter, "Chirality", nanobind::is_arithmetic())
      .value("S", ChiralCenter::Chirality::S_Chiral)
      .value("R", ChiralCenter::Chirality::R_Chiral);
  chiralCenter
      .def_ro("identifiers", &ChiralCenter::ids, "The center atom and three of its bonded neighbors.")
      .def_ro("chirality", &ChiralCenter::type, "The handedness the atoms have in the positions of the component.");

  // The bonded terms of a flexible molecule. The bonds are its connectivity: every bond, bend and torsion they form
  // needs a potential (Component checks this).
  nanobind::class_<Potentials::IntraMolecularPotentials>(m, "IntraMolecularPotentials")
      .def(
          "__init__",
          [](Potentials::IntraMolecularPotentials* potentials, std::vector<BondPotential> bonds,
             std::vector<BendPotential> bends, std::vector<TorsionPotential> torsions,
             std::vector<TorsionPotential> improper_torsions, std::vector<UreyBradleyPotential> urey_bradleys,
             std::vector<InversionBendPotential> inversion_bends, std::vector<BondBondPotential> bond_bonds,
             std::vector<BondBendPotential> bond_bends, std::vector<BondTorsionPotential> bond_torsions,
             std::vector<BendBendPotential> bend_bends, std::vector<BendTorsionPotential> bend_torsions,
             std::vector<VanDerWaalsPotential> van_der_waals, std::vector<CoulombPotential> coulombs,
             std::vector<CMAPMap> cmap_maps, std::vector<CMAPPotential> cmaps)
          {
            new (potentials) Potentials::IntraMolecularPotentials();
            potentials->bonds = std::move(bonds);
            potentials->bends = std::move(bends);
            potentials->torsions = std::move(torsions);
            potentials->improperTorsions = std::move(improper_torsions);
            potentials->ureyBradleys = std::move(urey_bradleys);
            potentials->inversionBends = std::move(inversion_bends);
            potentials->bondBonds = std::move(bond_bonds);
            potentials->bondBends = std::move(bond_bends);
            potentials->bondTorsions = std::move(bond_torsions);
            potentials->bendBends = std::move(bend_bends);
            potentials->bendTorsions = std::move(bend_torsions);
            potentials->vanDerWaals = std::move(van_der_waals);
            potentials->coulombs = std::move(coulombs);
            potentials->cmapMaps = std::move(cmap_maps);
            potentials->cmaps = std::move(cmaps);
          },
          nanobind::arg("bonds") = std::vector<BondPotential>{}, nanobind::arg("bends") = std::vector<BendPotential>{},
          nanobind::arg("torsions") = std::vector<TorsionPotential>{},
          nanobind::arg("improper_torsions") = std::vector<TorsionPotential>{},
          nanobind::arg("urey_bradleys") = std::vector<UreyBradleyPotential>{},
          nanobind::arg("inversion_bends") = std::vector<InversionBendPotential>{},
          nanobind::arg("bond_bonds") = std::vector<BondBondPotential>{},
          nanobind::arg("bond_bends") = std::vector<BondBendPotential>{},
          nanobind::arg("bond_torsions") = std::vector<BondTorsionPotential>{},
          nanobind::arg("bend_bends") = std::vector<BendBendPotential>{},
          nanobind::arg("bend_torsions") = std::vector<BendTorsionPotential>{},
          nanobind::arg("van_der_waals") = std::vector<VanDerWaalsPotential>{},
          nanobind::arg("coulombs") = std::vector<CoulombPotential>{},
          nanobind::arg("cmap_maps") = std::vector<CMAPMap>{}, nanobind::arg("cmaps") = std::vector<CMAPPotential>{})
      .def_ro("bonds", &Potentials::IntraMolecularPotentials::bonds)
      .def_ro("bends", &Potentials::IntraMolecularPotentials::bends)
      .def_ro("torsions", &Potentials::IntraMolecularPotentials::torsions)
      .def_ro("improper_torsions", &Potentials::IntraMolecularPotentials::improperTorsions)
      .def_ro("urey_bradleys", &Potentials::IntraMolecularPotentials::ureyBradleys)
      .def_ro("inversion_bends", &Potentials::IntraMolecularPotentials::inversionBends)
      .def_ro("bond_bonds", &Potentials::IntraMolecularPotentials::bondBonds)
      .def_ro("bond_bends", &Potentials::IntraMolecularPotentials::bondBends)
      .def_ro("bond_torsions", &Potentials::IntraMolecularPotentials::bondTorsions)
      .def_ro("bend_bends", &Potentials::IntraMolecularPotentials::bendBends)
      .def_ro("bend_torsions", &Potentials::IntraMolecularPotentials::bendTorsions)
      .def_ro("van_der_waals", &Potentials::IntraMolecularPotentials::vanDerWaals)
      .def_ro("coulombs", &Potentials::IntraMolecularPotentials::coulombs)
      .def_ro("cmap_maps", &Potentials::IntraMolecularPotentials::cmapMaps)
      .def_ro("cmaps", &Potentials::IntraMolecularPotentials::cmaps)
      .def_ro("chiral_centers", &Potentials::IntraMolecularPotentials::chiralCenters);

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
    nanobind::enum_<Move::Types> moveTypes(move, "Types", nanobind::is_arithmetic());
    for (const auto& [moveName, moveType] : moveTypeNames) moveTypes.value(moveName, moveType);

    nanobind::class_<MCMoveStatistics>(m, "MCMoveStatistics")
        // the statistics held by the object, not a copy of them: 'stats[MoveTypes.Translation].max_change = 0.5'
        // has to change the move
        .def("__getitem__",
             [](nanobind::handle self, Move::Types i)
             {
               MCMoveStatistics& statistics = nanobind::cast<MCMoveStatistics&>(self);
               return std::visit([&](auto& held)
                                 { return nanobind::cast(&held, nanobind::rv_policy::reference_internal, self); },
                                 statistics.stats[std::to_underlying(i)]);
             })
        .def("__repr__", &MCMoveStatistics::repr);

    nanobind::class_<WidomData>(m, "WidomData")
        .def_ro("total", &WidomData::total)
        .def_ro("excess", &WidomData::excess)
        .def_ro("ideal_gas", &WidomData::idealGas);

    nanobind::class_<PropertyWidom>(m, "PropertyWidom")
        .def("result", &PropertyWidom::result)
        // convert result to units of Kelvin
        .def(
            "chemical_potential_result", [](PropertyWidom& p, double T)
            { return Units::EnergyToKelvin * p.chemicalPotentialResult(1.0 / (Units::KB * T)); },
            nanobind::arg("temperature"))
        // convert result to units of Pascal
        .def(
            "fugacity_result", [](PropertyWidom& p, double T)
            { return Units::PressureConversionFactor * p.fugacityResult(1.0 / (Units::KB * T)); },
            nanobind::arg("temperature"));

    nanobind::class_<PropertyGibbsWidom>(m, "PropertyGibbsWidom")
        .def("result", &PropertyGibbsWidom::result)
        // convert result to units of Kelvin
        .def(
            "chemical_potential_result", [](PropertyGibbsWidom& p, double T)
            { return Units::EnergyToKelvin * p.chemicalPotentialResult(1.0 / (Units::KB * T)); },
            nanobind::arg("temperature"))
        // convert result to units of Pascal
        .def(
            "fugacity_result", [](PropertyGibbsWidom& p, double T)
            { return Units::PressureConversionFactor * p.fugacityResult(1.0 / (Units::KB * T)); },
            nanobind::arg("temperature"));

    // (bound before the component, which has it as a default)
  nanobind::enum_<PartitionFunction::Source>(m, "ThermochemicalDatabase", nanobind::is_arithmetic())
      .value("NASA", PartitionFunction::Source::NASA)
      .value("JANAF", PartitionFunction::Source::JANAF);
  m.def(
      "ln_partition_function",
      [](const std::string& species, double temperature, PartitionFunction::Source database)
      {
        if (!PartitionFunction::contains(species, database))
        {
          throw std::invalid_argument(std::format("no thermochemical data for species '{}' in the database", species));
        }
        if (!(temperature > 0.0)) throw std::invalid_argument("ln_partition_function: the temperature must be positive");
        return PartitionFunction::logPartitionFunction(species, temperature, database);
      },
      nanobind::arg("species"), nanobind::arg("temperature"),
      nanobind::arg("database") = PartitionFunction::Source::NASA,
      "The logarithm of the partition function of an isolated molecule of the species at the temperature [K], from "
      "the NASA polynomials (with the JANAF tables behind them) or from the JANAF tables only.");
  m.def(
      "has_thermochemical_data",
      [](const std::string& species, PartitionFunction::Source database)
      { return PartitionFunction::contains(species, database); },
      nanobind::arg("species"), nanobind::arg("database") = PartitionFunction::Source::NASA);

    nanobind::class_<MoleculePropertiesSettings>(m, "MoleculePropertiesSettings")
        .def(
            "__init__",
            [](MoleculePropertiesSettings* settings, std::size_t sample_every, std::optional<std::size_t> write_every,
               std::size_t number_of_bins, double bond_range, std::optional<double> end_to_end_range,
               double end_to_end_bin_width)
            {
              if (sample_every == 0 || number_of_bins == 0 || !(bond_range > 0.0) || !(end_to_end_bin_width > 0.0))
              {
                throw std::invalid_argument("MoleculePropertiesSettings: sample_every, number_of_bins, bond_range and "
                                            "end_to_end_bin_width must be positive");
              }
              new (settings) MoleculePropertiesSettings();
              settings->sampleEvery = sample_every;
              settings->writeEvery = write_every;
              settings->numberOfBins = number_of_bins;
              settings->bondRange = bond_range;
              settings->endToEndRange = end_to_end_range;
              settings->endToEndBinWidth = end_to_end_bin_width;
            },
            nanobind::arg("sample_every") = 10, nanobind::arg("write_every") = nanobind::none(),
            nanobind::arg("number_of_bins") = 128, nanobind::arg("bond_range") = 4.0,
            nanobind::arg("end_to_end_range") = nanobind::none(), nanobind::arg("end_to_end_bin_width") = 0.25,
            "Histograms of the bond lengths (up to 'bond_range' [Angstrom]), the bend and torsion angles, and the "
            "end-to-end distance (up to 'end_to_end_range', by default the contour length, in bins of "
            "'end_to_end_bin_width' [Angstrom]).")
        .def_ro("sample_every", &MoleculePropertiesSettings::sampleEvery)
        .def_ro("write_every", &MoleculePropertiesSettings::writeEvery)
        .def_ro("number_of_bins", &MoleculePropertiesSettings::numberOfBins)
        .def_ro("bond_range", &MoleculePropertiesSettings::bondRange)
        .def_ro("end_to_end_range", &MoleculePropertiesSettings::endToEndRange)
        .def_ro("end_to_end_bin_width", &MoleculePropertiesSettings::endToEndBinWidth);

    nanobind::class_<MoleculeShapeSettings>(m, "MoleculeShapeSettings")
        .def(
            "__init__",
            [](MoleculeShapeSettings* settings, std::size_t sample_every, std::optional<std::size_t> write_every,
               std::size_t number_of_bins, bool mass_weighted, std::optional<double> radius_of_gyration_range,
               double radius_of_gyration_bin_width)
            {
              if (sample_every == 0 || number_of_bins == 0 || !(radius_of_gyration_bin_width > 0.0))
              {
                throw std::invalid_argument("MoleculeShapeSettings: sample_every, number_of_bins and "
                                            "radius_of_gyration_bin_width must be positive");
              }
              new (settings) MoleculeShapeSettings();
              settings->sampleEvery = sample_every;
              settings->writeEvery = write_every;
              settings->numberOfBins = number_of_bins;
              settings->massWeighted = mass_weighted;
              settings->radiusOfGyrationRange = radius_of_gyration_range;
              settings->radiusOfGyrationBinWidth = radius_of_gyration_bin_width;
            },
            nanobind::arg("sample_every") = 10, nanobind::arg("write_every") = nanobind::none(),
            nanobind::arg("number_of_bins") = 128, nanobind::arg("mass_weighted") = false,
            nanobind::arg("radius_of_gyration_range") = nanobind::none(),
            nanobind::arg("radius_of_gyration_bin_width") = 0.1,
            "The radius of gyration and the shape descriptors, with every atom weighing the same unless "
            "'mass_weighted'.")
        .def_ro("sample_every", &MoleculeShapeSettings::sampleEvery)
        .def_ro("write_every", &MoleculeShapeSettings::writeEvery)
        .def_ro("number_of_bins", &MoleculeShapeSettings::numberOfBins)
        .def_ro("mass_weighted", &MoleculeShapeSettings::massWeighted)
        .def_ro("radius_of_gyration_range", &MoleculeShapeSettings::radiusOfGyrationRange)
        .def_ro("radius_of_gyration_bin_width", &MoleculeShapeSettings::radiusOfGyrationBinWidth);

    nanobind::class_<MoleculeBackboneSettings>(m, "MoleculeBackboneSettings")
        .def(
            "__init__",
            [](MoleculeBackboneSettings* settings, std::size_t sample_every, std::optional<std::size_t> write_every,
               std::size_t number_of_wave_vectors, double wave_vector_lower_limit, double wave_vector_upper_limit)
            {
              if (sample_every == 0 || number_of_wave_vectors == 0)
              {
                throw std::invalid_argument("MoleculeBackboneSettings: sample_every and number_of_wave_vectors must "
                                            "be positive");
              }
              if (!(wave_vector_lower_limit > 0.0) || !(wave_vector_upper_limit > wave_vector_lower_limit))
              {
                throw std::invalid_argument("MoleculeBackboneSettings: wave_vector_lower_limit must be positive and "
                                            "below wave_vector_upper_limit");
              }
              new (settings) MoleculeBackboneSettings();
              settings->sampleEvery = sample_every;
              settings->writeEvery = write_every;
              settings->numberOfWaveVectors = number_of_wave_vectors;
              settings->waveVectorLowerLimit = wave_vector_lower_limit;
              settings->waveVectorUpperLimit = wave_vector_upper_limit;
            },
            nanobind::arg("sample_every") = 10, nanobind::arg("write_every") = nanobind::none(),
            nanobind::arg("number_of_wave_vectors") = 64, nanobind::arg("wave_vector_lower_limit") = 0.01,
            nanobind::arg("wave_vector_upper_limit") = 5.0,
            "The statistics of the backbone; the form factor is on a logarithmic grid of wave vectors [1/Angstrom].")
        .def_ro("sample_every", &MoleculeBackboneSettings::sampleEvery)
        .def_ro("write_every", &MoleculeBackboneSettings::writeEvery)
        .def_ro("number_of_wave_vectors", &MoleculeBackboneSettings::numberOfWaveVectors)
        .def_ro("wave_vector_lower_limit", &MoleculeBackboneSettings::waveVectorLowerLimit)
        .def_ro("wave_vector_upper_limit", &MoleculeBackboneSettings::waveVectorUpperLimit);

    nanobind::class_<EndToEndACFSettings>(m, "EndToEndACFSettings")
        .def(
            "__init__",
            [](EndToEndACFSettings* settings, std::size_t sample_every, std::optional<std::size_t> write_every,
               std::size_t number_of_block_elements)
            {
              if (sample_every == 0 || number_of_block_elements < 2)
              {
                throw std::invalid_argument("EndToEndACFSettings: sample_every must be positive and "
                                            "number_of_block_elements at least 2");
              }
              new (settings) EndToEndACFSettings();
              settings->sampleEvery = sample_every;
              settings->writeEvery = write_every;
              settings->numberOfBlockElements = number_of_block_elements;
            },
            nanobind::arg("sample_every") = 10, nanobind::arg("write_every") = nanobind::none(),
            nanobind::arg("number_of_block_elements") = 25,
            "The autocorrelation function of the end-to-end vector, by the order-N scheme.")
        .def_ro("sample_every", &EndToEndACFSettings::sampleEvery)
        .def_ro("write_every", &EndToEndACFSettings::writeEvery)
        .def_ro("number_of_block_elements", &EndToEndACFSettings::numberOfBlockElements);

    nanobind::class_<Component> component(m, "Component");
    nanobind::enum_<Component::Type>(component, "Type", nanobind::is_arithmetic())
        .value("ADSORBATE", Component::Type::Adsorbate)
        .value("CATION", Component::Type::Cation);

    component
        .def_static(
            "from_file",
            [](const ForceField& force_field, const std::string& file_name, std::optional<std::string> component_name,
               const MCMoveProbabilities& particle_probabilities, Component::Type type,
               std::optional<double> fugacity_coefficient, bool thermodynamic_integration,
               std::size_t number_of_blocks, std::size_t number_of_lambda_bins)
            {
              if (number_of_blocks == 0) throw std::invalid_argument("Component: number_of_blocks must be positive");
              // (a lambda grid of no points has no first bin, and one of a single point no spacing)
              if (number_of_lambda_bins < 2)
                throw std::invalid_argument("Component: number_of_lambda_bins must be at least 2");
              const std::string name = component_name.value_or(std::filesystem::path(file_name).stem().string());
              try
              {
                // the index of the component is that of its place in the system, which the system gives its atoms
                return Component(type, 0, force_field, name, file_name, number_of_blocks, number_of_lambda_bins,
                                 particle_probabilities, fugacity_coefficient, thermodynamic_integration);
              }
              catch (const nlohmann::json::exception& error)
              {
                throw std::invalid_argument(std::format("Component: file '{}': {}", file_name, error.what()));
              }
              catch (const std::runtime_error& error)
              {
                throw std::invalid_argument(error.what());
              }
            },
            nanobind::arg("force_field"), nanobind::arg("file_name"), nanobind::arg("component_name") = nanobind::none(),
            nanobind::arg("particle_probabilities") = MCMoveProbabilities(),
            nanobind::arg("type") = Component::Type::Adsorbate,
            nanobind::arg("fugacity_coefficient") = nanobind::none(),
            nanobind::arg("thermodynamic_integration") = false, nanobind::arg("number_of_blocks") = 5,
            nanobind::arg("number_of_lambda_bins") = 41,
            "Reads a component from a molecule file (JSON), as the command-line program reads the molecule of a "
            "component of 'simulation.json': the extension '.json' may be left out, and a file that is not found is "
            "looked for in the directory of RASPA_DIR. The pseudo-atoms of the molecule have to be those of the "
            "force field. 'component_name' is the name of the file without it.")
        .def("__init__",
             [](Component* component, const ForceField& force_field, std::string component_name,
                double critical_temperature, double critical_pressure, double acentric_factor,
                std::vector<Atom> defined_atoms, const ConnectivityTable& connectivity_table,
                const Potentials::IntraMolecularPotentials& intra_molecular_potentials, std::size_t number_of_blocks,
                std::size_t number_of_lambda_bins, const MCMoveProbabilities& particle_probabilities,
                std::optional<double> fugacity_coefficient, bool thermodynamic_integration,
                std::vector<double4> blocking_pockets, std::vector<std::vector<std::size_t>> rigid_bodies,
                const std::vector<std::array<std::size_t, 4>>& chiral_centers,
                double intra_14_van_der_waals_scaling, double intra_14_charge_charge_scaling,
                const std::vector<std::vector<std::size_t>>& partial_reinsertion_fixed_atoms,
                const std::vector<std::vector<std::size_t>>& repeat_units,
                std::optional<std::array<std::size_t, 2>> end_to_end_atoms)
             {
               for (const Atom& atom : defined_atoms)
               {
                 if (static_cast<std::size_t>(atom.type) >= force_field.pseudoAtoms.size())
                 {
                   throw std::invalid_argument(std::format("Component: an atom has type {}, but the force field has "
                                                           "{} pseudo-atoms",
                                                           static_cast<std::size_t>(atom.type),
                                                           force_field.pseudoAtoms.size()));
                 }
               }
               if (number_of_blocks == 0) throw std::invalid_argument("Component: number_of_blocks must be positive");
               // (a lambda grid of no points has no first bin, and one of a single point no spacing)
               if (number_of_lambda_bins < 2)
                 throw std::invalid_argument("Component: number_of_lambda_bins must be at least 2");
               const ConnectivityTable connectivity = connectivityOfComponent(
                   defined_atoms.size(), connectivity_table, intra_molecular_potentials, rigid_bodies);
               try
               {
                 new (component)
                     Component(force_field, component_name, critical_temperature, critical_pressure, acentric_factor,
                               defined_atoms, connectivity, intra_molecular_potentials, number_of_blocks,
                               number_of_lambda_bins, particle_probabilities, fugacity_coefficient,
                               thermodynamic_integration, blocking_pockets, rigid_bodies,
                               intra_14_van_der_waals_scaling, intra_14_charge_charge_scaling);
               }
               catch (const std::runtime_error& error)
               {
                 throw std::invalid_argument(std::format("Component: {}", error.what()));
               }
               try
               {
                 describeFlexibleComponent(*component, chiral_centers, partial_reinsertion_fixed_atoms, repeat_units,
                                           end_to_end_atoms);
               }
               catch (...)
               {
                 component->~Component();
                 throw;
               }
             },
             nanobind::arg("force_field"), nanobind::arg("component_name"), nanobind::arg("critical_temperature"),
             nanobind::arg("critical_pressure"), nanobind::arg("acentric_factor"),
             nanobind::arg("defined_atoms") = std::vector<Atom>(),
             nanobind::arg("connectivity_table") = ConnectivityTable(),
             nanobind::arg("intra_molecular_potentials") = Potentials::IntraMolecularPotentials(),
             nanobind::arg("number_of_blocks") = 5, nanobind::arg("number_of_lambda_bins") = 41,
             nanobind::arg("particle_probabilities") = MCMoveProbabilities(),
             nanobind::arg("fugacity_coefficient") = nanobind::none(),
             nanobind::arg("thermodynamic_integration") = false,
             nanobind::arg("blocking_pockets") = std::vector<double4>(),
             nanobind::arg("rigid_bodies") = std::vector<std::vector<std::size_t>>(),
             nanobind::arg("chiral_centers") = std::vector<std::array<std::size_t, 4>>(),
             nanobind::arg("intra_14_van_der_waals_scaling") = 0.0,
             nanobind::arg("intra_14_charge_charge_scaling") = 0.0,
             nanobind::arg("partial_reinsertion_fixed_atoms") = std::vector<std::vector<std::size_t>>(),
             nanobind::arg("repeat_units") = std::vector<std::vector<std::size_t>>(),
             nanobind::arg("end_to_end_atoms") = nanobind::none())
        .def_ro("name", &Component::name)
        .def_ro("rigid", &Component::rigid, "Whether the molecule is a single rigid body.")
        .def_prop_ro(
            "number_of_blocks", [](const Component& c) { return c.numberOfBlocks(); },
            "The number of blocks the statistics of the component are kept in. A system sizes its components for "
            "its own.")
        .def_prop_ro(
            "rigid_bodies",
            [](const Component& c)
            {
              std::vector<std::vector<std::size_t>> bodies{};
              for (const Fragment& fragment : c.fragmentGraph.fragments)
              {
                if (fragment.isRigidBody()) bodies.push_back(fragment.atoms);
              }
              return bodies;
            },
            "The groups of atoms that move as one rigid unit.")
        .def_prop_ro(
            "chiral_centers", [](const Component& c) { return c.intraMolecularPotentials.chiralCenters; },
            "The chiral centers, with the handedness they have in the positions of the atoms.")
        .def_ro("partial_reinsertion_fixed_atoms", &Component::partialReinsertionFixedAtoms,
                "The sets of atoms that stay where they are in a partial reinsertion; the others are grown again.")
        .def_ro("repeat_units", &Component::repeatUnits,
                "The repeat units of a chain molecule, for the reptation move.")
        .def_ro("end_to_end_atoms", &Component::endToEndAtoms,
                "The two atoms whose distance is the end-to-end distance of the molecule.")
        .def("backbone_atoms", &Component::backboneAtoms,
             "The atoms on the shortest path of bonds between the end-to-end atoms.")
        .def_ro("intra_14_van_der_waals_scaling", &Component::intra14VanDerWaalsScaling)
        .def_ro("intra_14_charge_charge_scaling", &Component::intra14ChargeChargeScaling)
        .def_ro("connectivity_table", &Component::connectivityTable)
        .def_ro("intra_molecular_potentials", &Component::intraMolecularPotentials)
        .def_prop_rw(
          "lambda_histogram", [](Component& c) -> PropertyLambdaProbabilityHistogram& { return c.lambdaGC; },
          [](Component& c, const PropertyLambdaProbabilityHistogram& value)
          {
            // the moves and the drivers index the histogram by the lambda bins and the blocks of the component
            if (value.numberOfSamplePoints != c.lambdaGC.numberOfSamplePoints ||
                value.numberOfBlocks != c.lambdaGC.numberOfBlocks ||
                value.biasFactor.size() != value.numberOfSamplePoints ||
                value.histogram.size() != value.numberOfSamplePoints)
            {
              throw std::invalid_argument(std::format("lambda_histogram needs {} lambda bins and {} blocks, as the component",
                                                      c.lambdaGC.numberOfSamplePoints, c.lambdaGC.numberOfBlocks));
            }
            c.lambdaGC = value;
          },
          nanobind::rv_policy::reference_internal)
        .def_ro("mc_moves_probabilities", &Component::mc_moves_probabilities)
        .def_ro("mc_moves_statistics", &Component::mc_moves_statistics)
        .def_rw("blocking_pockets", &Component::blockingPockets)
        .def_rw("average_rosenbluth_weights", &Component::averageRosenbluthWeights)
        .def_rw("average_gibbs_rosenbluth_weights", &Component::averageGibbsRosenbluthWeights)
        .def("print_status", &Component::printStatus)
        .def("__repr__", &Component::repr);

    // The settings below are those of the component a system is constructed with: the system keeps its own copy, so
    // a setting is made before the component is handed to it.
    component
        .def_rw("type", &Component::type, "Adsorbate or cation.")
        .def_prop_rw(
            "mol_fraction", [](const Component& c) { return c.molFraction; },
            [](Component& c, double fraction)
            {
              if (!(fraction >= 0.0)) throw std::invalid_argument("mol_fraction can not be negative");
              c.molFraction = fraction;
            },
            "The mole fraction of the component in the reservoir [-].")
        .def_prop_rw(
            "fugacity_coefficient", [](const Component& c) { return c.fugacityCoefficient; },
            [](Component& c, std::optional<double> coefficient)
            {
              if (coefficient.has_value() && !(*coefficient > 0.0))
              {
                throw std::invalid_argument("fugacity_coefficient must be positive (None for the equation of state)");
              }
              c.fugacityCoefficient = coefficient;
              c.fugacityCoefficientGiven = coefficient.has_value();
            },
            nanobind::for_setter(nanobind::arg("value").none()),
            "The fugacity coefficient [-]; None takes it from the equation of state.")
        .def_prop_rw(
            "ln_partition_function", [](const Component& c) { return c.lnPartitionFunction; },
            [](Component& c, double value)
            {
              // a number that is given is not one of a species
              c.lnPartitionFunction = value;
              c.lnPartitionFunctionSpecies.reset();
              c.lnPartitionFunctionFromJANAF = false;
            },
            "The logarithm of the partition function of the isolated molecule, for the reaction moves.")
        .def(
            "set_ln_partition_function_species",
            [](Component& c, const std::string& species, PartitionFunction::Source database)
            {
              if (!PartitionFunction::contains(species, database))
              {
                throw std::invalid_argument(std::format("no thermochemical data for species '{}' in the database",
                                                        species));
              }
              c.lnPartitionFunctionSpecies = species;
              c.lnPartitionFunctionFromJANAF = database == PartitionFunction::Source::JANAF;
            },
            nanobind::arg("species"), nanobind::arg("database") = PartitionFunction::Source::NASA,
            "Takes the partition function from the thermochemical data of a species, at the temperature of the "
            "system the component is given to.")
        .def_ro("ln_partition_function_species", &Component::lnPartitionFunctionSpecies)
        .def_ro("critical_temperature", &Component::criticalTemperature, "[K]")
        .def_ro("critical_pressure", &Component::criticalPressure, "[Pa]")
        .def_ro("acentric_factor", &Component::acentricFactor, "[-]")
        .def_ro("net_charge", &Component::netCharge, "[e]")
        .def_ro("swappable", &Component::swappable, "Whether a move of a system inserts or deletes the component.")
        .def_prop_rw(
            "starting_bead", [](const Component& c) { return c.startingBead; },
            [](Component& c, std::size_t bead)
            {
              if (bead >= c.definedAtoms.size()) throw std::invalid_argument("starting_bead is not an atom of the component");
              c.setStartingBead(bead);
            },
            "The atom the growth of the molecule starts from.")
        .def_rw("ideal_gas_rosenbluth_weight", &Component::idealGasRosenbluthWeight,
                nanobind::for_setter(nanobind::arg("value").none()),
                "The Rosenbluth weight of the isolated molecule [-]; None computes it at the start of the run.")
        .def_prop_rw(
            "ideal_gas_rosenbluth_weight_tolerance",
            [](const Component& c) { return c.idealGasRosenbluthWeightTolerance; },
            [](Component& c, std::optional<double> tolerance)
            {
              if (tolerance.has_value() && !(*tolerance > 0.0))
              {
                throw std::invalid_argument("ideal_gas_rosenbluth_weight_tolerance must be positive");
              }
              c.idealGasRosenbluthWeightTolerance = tolerance;
            },
            nanobind::for_setter(nanobind::arg("value").none()),
            "The relative standard error to which a weight that is not given is computed [-].")
        .def_prop_rw(
            "constraint_tolerance", [](const Component& c) { return c.constraintTolerance; },
            [](Component& c, double tolerance)
            {
              if (!(tolerance > 0.0)) throw std::invalid_argument("constraint_tolerance must be positive");
              c.constraintTolerance = tolerance;
            },
            "The relative error to which molecular dynamics restores a FIXED bond length (SHAKE/RATTLE) [-].")
        .def_prop_rw(
            "constraint_maximum_iterations", [](const Component& c) { return c.constraintMaximumIterations; },
            [](Component& c, std::size_t iterations)
            {
              if (iterations == 0) throw std::invalid_argument("constraint_maximum_iterations must be positive");
              c.constraintMaximumIterations = iterations;
            },
            "The iteration limit of SHAKE and of RATTLE per molecule and step.")
        .def_prop_ro(
            "distance_constraints",
            [](const Component& c)
            {
              std::vector<std::tuple<std::size_t, std::size_t, double>> constraints;
              for (const Component::DistanceConstraint& constraint : c.distanceConstraints)
                constraints.emplace_back(constraint.a, constraint.b, constraint.length);
              return constraints;
            },
            "The distance constraints (atom, atom, length [Å]) that molecular dynamics holds for the FIXED terms.")
        .def_prop_rw(
            "ideal_gas_energy",
            [](const Component& c) -> std::optional<double>
            {
              if (!c.idealGasEnergy.has_value()) return std::nullopt;
              return Units::EnergyToKelvin * *c.idealGasEnergy;
            },
            [](Component& c, std::optional<double> energy)
            {
              c.idealGasEnergy = energy.has_value() ? std::optional<double>(Units::KelvinToEnergy * *energy) : std::nullopt;
            },
            nanobind::for_setter(nanobind::arg("value").none()),
            "The mean intramolecular energy of the isolated molecule [K], the reference of the enthalpy of adsorption.")
        .def_rw("pair_component", &Component::pairComponentId, nanobind::for_setter(nanobind::arg("value").none()),
                "The index of the partner component of the pair-swap moves.")
        .def_rw("maximum_pair_distance", &Component::maximumPairDistance,
                nanobind::for_setter(nanobind::arg("value").none()),
                "The largest separation of the pair of the pair-swap moves [Angstrom].")
        .def_rw("group_components", &Component::groupComponentIds,
                "The indices of the satellite components of the group-swap moves, one entry per molecule.")
        .def_rw("maximum_group_distance", &Component::maximumGroupDistance,
                nanobind::for_setter(nanobind::arg("value").none()),
                "The largest distance of a satellite from the central molecule of the group-swap moves [Angstrom].")
        .def_rw("tethered_proton_hop_sites", &Component::tetheredProtonHopSiteGroups,
                "The sites of the tethered proton-hop move: per tethering site a list of fractional positions.")
        .def_rw("identity_changes", &Component::identityChanges,
                "The indices of the components the identity-change move changes the component into.")
        .def_rw("identity_switches", &Component::identitySwitches,
                "The indices of the partner components of the identity-switch move.")
        .def_rw("gibbs_identity_changes", &Component::gibbsIdentityChanges,
                "The indices of the components the Gibbs identity-change move changes the component into.")
        .def_prop_ro(
            "lambda_gibbs_histogram", [](Component& c) -> PropertyLambdaProbabilityHistogram& { return c.lambdaGibbs; },
            nanobind::rv_policy::reference_internal)
        .def_prop_ro(
            "lambda_pair_swap_histogram",
            [](Component& c) -> PropertyLambdaProbabilityHistogram& { return c.lambdaPairSwap; },
            nanobind::rv_policy::reference_internal)
        .def_prop_ro(
            "lambda_group_swap_histogram",
            [](Component& c) -> PropertyLambdaProbabilityHistogram& { return c.lambdaGroupSwap; },
            nanobind::rv_policy::reference_internal);

    // The analyses of the molecules of the component. A driver creates them for its systems; the results are read
    // from 'property_molecule_properties', 'property_molecule_shape', 'property_molecule_backbone' and
    // 'property_end_to_end_acf' of a system.
    component
        .def_rw("molecule_properties_settings", &Component::moleculePropertiesSettings,
                nanobind::for_setter(nanobind::arg("value").none()),
                "Histograms of the bond lengths, bend angles, torsion angles and the end-to-end distance.")
        .def_rw("molecule_shape_settings", &Component::moleculeShapeSettings,
                nanobind::for_setter(nanobind::arg("value").none()),
                "The radius of gyration and the shape descriptors of the gyration tensor.")
        .def_rw("molecule_backbone_settings", &Component::moleculeBackboneSettings,
                nanobind::for_setter(nanobind::arg("value").none()),
                "The statistics of the backbone of a chain: internal distances, bond correlations, form factor.")
        .def_rw("end_to_end_acf_settings", &Component::endToEndACFSettings,
                nanobind::for_setter(nanobind::arg("value").none()),
                "The autocorrelation function of the end-to-end vector.");

    defineCount(component, "number_of_trial_molecules", &Component::numberOfTrialMolecules,
                "The number of complete trial molecules of the multiple-try moves.");
    defineCount(component, "pivot_cbmc_number_of_trial_angles", &Component::pivotCBMCNumberOfTrialAngles,
                "The number of trial angles of the configurational-bias pivot move.");
    defineCount(component, "crankshaft_max_segment_size", &Component::crankshaftMaxSegmentSize,
                "The largest number of atoms in the rotated segment of the crankshaft move.");
    defineFraction(component, "pivot_randomization_fraction", &Component::pivotRandomizationFraction,
                   "The fraction of pivot attempts that draw the angle from the full circle.");
    defineFraction(component, "pivot_cbmc_randomization_fraction", &Component::pivotCBMCRandomizationFraction,
                   "The fraction of configurational-bias pivot attempts that draw the trial angles from the full circle.");
    defineFraction(component, "crankshaft_randomization_fraction", &Component::crankshaftRandomizationFraction,
                   "The fraction of crankshaft attempts that draw the angle from the full circle.");
    defineFraction(component, "concerted_rotation_randomization_fraction",
                   &Component::concertedRotationRandomizationFraction,
                   "The fraction of concerted-rotation attempts that draw the driver angle from the full circle.");
    defineFraction(component, "bead_flip_randomization_fraction", &Component::beadFlipRandomizationFraction,
                   "The fraction of bead-flip attempts that randomize the rotation completely.");

    nanobind::class_<LoadingData>(m, "LoadingData")
        .def(nanobind::init<std::size_t>(), nanobind::arg("number_of_components"))
        .def_ro("number_of_molecules", &LoadingData::numberOfMolecules)
        .def_ro("number_densities", &LoadingData::numberDensities)
        .def_ro("inverse_number_densities", &LoadingData::inverseNumberDensities);

    nanobind::class_<SampleMovie>(m, "SampleMovie")
        .def(nanobind::init<std::size_t, std::size_t, bool, std::optional<std::string>>(), nanobind::arg("system_id"),
             nanobind::arg("sample_every"), nanobind::arg("restrict_to_box") = true,
             nanobind::arg("tag") = std::nullopt)
        // the movie is written under the index the system has in its driver: 'system_id' is not used
        .def(
            "__init__",
            [](SampleMovie* movie, std::size_t sample_every, bool restrict_to_box, std::optional<std::string> tag)
            {
              if (sample_every == 0) throw std::invalid_argument("SampleMovie: sample_every must be positive");
              new (movie) SampleMovie(0, sample_every, restrict_to_box, tag);
            },
            nanobind::arg("sample_every") = 1, nanobind::arg("restrict_to_box") = true,
            nanobind::arg("tag") = std::nullopt)
        .def_ro("sample_every", &SampleMovie::sampleEvery)
        .def_ro("restrict_to_box", &SampleMovie::restrictToBox)
        .def_ro("tag", &SampleMovie::tag);

    nanobind::class_<PropertyLoading>(m, "PropertyLoading")
        .def(nanobind::init<std::size_t, std::size_t>(), nanobind::arg("number_of_blocks"),
             nanobind::arg("number_of_components"))
        .def("result", &PropertyLoading::average)
        .def("average_loading_number_of_molecules", &PropertyLoading::averageLoadingNumberOfMolecules)
        .def("write_averages_statistics", &PropertyLoading::writeAveragesStatistics)
        .def("__repr__", &PropertyLoading::repr);

    nanobind::class_<EnergyDuDlambda>(m, "EnergyDuDlambda")
        .def_ro("energy", &EnergyDuDlambda::energy)
        .def_ro("dudlambda", &EnergyDuDlambda::dUdlambda);

    nanobind::class_<EnergyStatus>(m, "EnergyStatus")
        .def_rw("total_energy", &EnergyStatus::totalEnergy)
        .def("__repr__", &EnergyStatus::repr);

    nanobind::class_<PropertyEnergy>(m, "PropertyEnergy")
        // convert result to units of Kelvin
        .def("result", [](PropertyEnergy& p) { return Units::EnergyToKelvin * p.average(); })
        .def("__repr__", &PropertyEnergy::repr);

    nanobind::class_<PressureData>(m, "PressureData")
        .def(nanobind::init<>())
        .def_ro("total_pressure", &PressureData::totalPressure)
        .def_ro("excess_pressure", &PressureData::excessPressure)
        .def_ro("ideal_gas_pressure", &PressureData::idealGasPressure);

    nanobind::class_<PropertyPressure>(m, "PropertyPressure")
        .def("result", [](PropertyPressure& p) { return Units::PressureConversionFactor * p.average(); });
    //.def("__repr__", &Pressures::repr);

    nanobind::class_<AverageEnergyType>(m, "AverageEnergyType")
        .def(nanobind::init<double, double, double, double>(), nanobind::arg("total_energy") = 0.0,
             nanobind::arg("van_der_waals_energy") = 0.0, nanobind::arg("coulomb_energy") = 0.0,
             nanobind::arg("polarization_energy") = 0.0)
        .def_rw("total_energy", &AverageEnergyType::totalEnergy)
        .def_rw("van_der_waals_energy", &AverageEnergyType::VanDerWaalsEnergy)
        .def_rw("coulomb_energy", &AverageEnergyType::CoulombEnergy)
        .def_rw("polarization_energy", &AverageEnergyType::polarizationEnergy);

    nanobind::class_<PropertyDensityGrid> property_energy_grid(m, "PropertyDensityGrid");

    nanobind::enum_<PropertyDensityGrid::Normalization>(property_energy_grid, "Normalization",
                                                        nanobind::is_arithmetic())
        .value("MAX", PropertyDensityGrid::Normalization::Max)
        .value("NUMBER_DENSITY", PropertyDensityGrid::Normalization::NumberDensity);

    nanobind::enum_<PropertyDensityGrid::Binning>(property_energy_grid, "Binning", nanobind::is_arithmetic())
        .value("STANDARD", PropertyDensityGrid::Binning::Standard)
        .value("EQUITABLE", PropertyDensityGrid::Binning::Equitable);

    property_energy_grid.def(
        nanobind::init<std::size_t, std::size_t, int3, std::size_t, std::size_t, std::vector<std::size_t>,
                       PropertyDensityGrid::Normalization, PropertyDensityGrid::Binning>(),
        nanobind::arg("number_of_frameworks"), nanobind::arg("number_of_components"),
        nanobind::arg("number_of_grid_points"), nanobind::arg("sample_every"), nanobind::arg("write_every"),
        nanobind::arg("density_grid_pseudo_atoms_list"),
        nanobind::arg("normalization_type") = PropertyDensityGrid::Normalization::Max,
        nanobind::arg("binning_mode") = PropertyDensityGrid::Binning::Standard);
    // the system sizes the grid for its framework and components
    property_energy_grid.def(
        "__init__",
        [](PropertyDensityGrid* grid, int3 number_of_grid_points, std::size_t sample_every, std::size_t write_every,
           std::vector<std::size_t> density_grid_pseudo_atoms_list, PropertyDensityGrid::Normalization normalization,
           PropertyDensityGrid::Binning binning_mode)
        {
          if (sample_every == 0 || write_every == 0)
            throw std::invalid_argument("PropertyDensityGrid: sample_every and write_every must be positive");
          new (grid) PropertyDensityGrid(0, 0, number_of_grid_points, sample_every, write_every,
                                         density_grid_pseudo_atoms_list, normalization, binning_mode);
        },
        nanobind::arg("number_of_grid_points"), nanobind::arg("sample_every") = 10,
        nanobind::arg("write_every") = 5000,
        nanobind::arg("density_grid_pseudo_atoms_list") = std::vector<std::size_t>(),
        nanobind::arg("normalization_type") = PropertyDensityGrid::Normalization::Max,
        nanobind::arg("binning_mode") = PropertyDensityGrid::Binning::Standard);

    nanobind::class_<PropertyEnergyHistogram> energy_histogram(m, "PropertyEnergyHistogram");
    energy_histogram
      .def(nanobind::init<std::size_t, std::size_t, std::pair<double, double>, std::size_t, std::optional<std::size_t>>(),
                 nanobind::arg("number_of_blocks"), 
                 nanobind::arg("number_of_bins"), 
                 nanobind::arg("value_range"),
                 nanobind::arg("sample_every"), 
                 nanobind::arg("write_every") = nanobind::none())
      // the system sizes the histogram for its blocks
      .def("__init__",
           [](PropertyEnergyHistogram* histogram, std::size_t number_of_bins, std::pair<double, double> value_range,
              std::size_t sample_every, std::optional<std::size_t> write_every)
           {
             if (sample_every == 0) throw std::invalid_argument("PropertyEnergyHistogram: sample_every must be positive");
             new (histogram) PropertyEnergyHistogram(1, number_of_bins, value_range, sample_every, write_every);
           },
           nanobind::arg("number_of_bins"), nanobind::arg("value_range"), nanobind::arg("sample_every") = 1,
           nanobind::arg("write_every") = nanobind::none())
      // the bins are in Kelvin already (PropertyEnergyHistogram::bins): converting again scaled them twice
      .def("result", [](PropertyEnergyHistogram& p) { return p.result(); });

  nanobind::class_<PropertyNumberOfMoleculesHistogram>(m, "PropertyNumberOfMoleculesHistogram")
      .def(nanobind::init<std::size_t, std::pair<std::size_t, std::size_t>, std::size_t, std::optional<std::size_t>>(),
                 nanobind::arg("number_of_blocks"), 
                 nanobind::arg("value_range"),
                 nanobind::arg("sample_every"), 
                 nanobind::arg("write_every") = nanobind::none())
      // the system sizes the histogram for its blocks and components
      .def("__init__",
           [](PropertyNumberOfMoleculesHistogram* histogram, std::pair<std::size_t, std::size_t> value_range,
              std::size_t sample_every, std::optional<std::size_t> write_every)
           {
             if (sample_every == 0)
               throw std::invalid_argument("PropertyNumberOfMoleculesHistogram: sample_every must be positive");
             new (histogram) PropertyNumberOfMoleculesHistogram(1, value_range, sample_every, write_every);
           },
           nanobind::arg("value_range"), nanobind::arg("sample_every") = 1,
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
       // the system sizes the histograms for its blocks and the pseudo-atoms of its force field
       .def("__init__",
            [](PropertyRadialDistributionFunction* rdf, std::size_t number_of_bins, double range,
               std::size_t sample_every, std::size_t write_every)
            {
              if (sample_every == 0 || write_every == 0)
                throw std::invalid_argument("PropertyRadialDistributionFunction: sample_every and write_every must "
                                            "be positive");
              new (rdf) PropertyRadialDistributionFunction(1, 0, number_of_bins, range, sample_every, write_every);
            },
            nanobind::arg("number_of_bins"), nanobind::arg("range"), nanobind::arg("sample_every") = 10,
            nanobind::arg("write_every") = 5000)
       .def("result", &PropertyRadialDistributionFunction::result);


  nanobind::class_<MeanSquaredDisplacementData>(m, "MeanSquaredDisplacementData")
      .def_ro("time", &MeanSquaredDisplacementData::time)
      .def_ro("xyz", &MeanSquaredDisplacementData::xyz)
      .def_ro("x", &MeanSquaredDisplacementData::x)
      .def_ro("y", &MeanSquaredDisplacementData::y)
      .def_ro("z", &MeanSquaredDisplacementData::z)
      .def_ro("number_of_samples", &MeanSquaredDisplacementData::numberOfSamples);

  nanobind::class_<PropertyMeanSquaredDisplacement>(m, "PropertyMeanSquaredDisplacement")
       .def(nanobind::init<std::size_t, std::optional<std::size_t>, std::size_t>(),
            nanobind::arg("sample_every"), 
            nanobind::arg("write_every") = nanobind::none(),
            nanobind::arg("number_of_block_elements_msd") = 25)
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


  // The number of molecules per component over the cycles of all phases. The system sizes it for its components
  // (the number given here was never used), and the driver for its cycles when none are given.
  {
    auto construct = [](PropertyNumberOfMoleculesEvolution* property, std::size_t number_of_cycles,
                        std::optional<std::size_t> number_of_components, std::size_t sample_every,
                        std::optional<std::size_t> write_every)
    {
      if (sample_every == 0)
        throw std::invalid_argument("PropertyNumberOfMoleculesEvolution: sample_every must be positive");
      new (property) PropertyNumberOfMoleculesEvolution(number_of_cycles, number_of_components.value_or(0),
                                                        sample_every, write_every);
    };
    nanobind::class_<PropertyNumberOfMoleculesEvolution> evolution(m, "PropertyNumberOfMoleculesEvolution");
    // 'number_of_cycles' is the keyword of the 3.1.0 release
    for (const char* cycles : {"number_of_production_cycles", "number_of_cycles"})
    {
      evolution.def("__init__", construct, nanobind::arg(cycles),
                    nanobind::arg("number_of_components") = nanobind::none(), nanobind::arg("sample_every") = 1,
                    nanobind::arg("write_every") = nanobind::none());
    }
    evolution
        .def(
            "__init__",
            [construct](PropertyNumberOfMoleculesEvolution* property, std::size_t sample_every,
                        std::optional<std::size_t> write_every)
            { construct(property, 0, std::nullopt, sample_every, write_every); },
            nanobind::arg("sample_every") = 1, nanobind::arg("write_every") = nanobind::none())
        .def_ro("sample_every", &PropertyNumberOfMoleculesEvolution::sampleEvery)
        .def_ro("write_every", &PropertyNumberOfMoleculesEvolution::writeEvery)
        .def_ro("result", &PropertyNumberOfMoleculesEvolution::result);
  }

  nanobind::class_<PropertyVolumeEvolution> volumeEvolution(m, "PropertyVolumeEvolution");
  bindEvolutionConstructors<PropertyVolumeEvolution>(volumeEvolution);
  volumeEvolution.def_ro("result", &PropertyVolumeEvolution::result);

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

  nanobind::enum_<MolecularDynamicsEnsemble>(m, "MolecularDynamicsEnsemble", nanobind::is_arithmetic())
      .value("NVE", MolecularDynamicsEnsemble::NVE)
      .value("NVT", MolecularDynamicsEnsemble::NVT)
      .value("NPT", MolecularDynamicsEnsemble::NPT)
      .value("NPTPR", MolecularDynamicsEnsemble::NPTPR)
      .value("MUVT", MolecularDynamicsEnsemble::MuVT)
      .value("MUPT", MolecularDynamicsEnsemble::MuPT)
      .value("MUPTPR", MolecularDynamicsEnsemble::MuPTPR);

  nanobind::enum_<BarostatCoupling>(m, "BarostatCoupling", nanobind::is_arithmetic())
      .value("MOLECULAR", BarostatCoupling::Molecular)
      .value("ATOMIC", BarostatCoupling::Atomic);

  nanobind::enum_<CellMinimizationType>(m, "CellType", nanobind::is_arithmetic())
      .value("FIXED", CellMinimizationType::Fixed)
      .value("ISOTROPIC", CellMinimizationType::Isotropic)
      .value("ANISOTROPIC", CellMinimizationType::Anisotropic)
      .value("MONOCLINIC", CellMinimizationType::Monoclinic)
      .value("REGULAR", CellMinimizationType::Regular)
      .value("REGULAR_UPPER_TRIANGLE", CellMinimizationType::RegularUpperTriangle)
      .value("MONOCLINIC_UPPER_TRIANGLE", CellMinimizationType::MonoclinicUpperTriangle);

  nanobind::enum_<MonoclinicAngleType>(m, "MonoclinicAngleType", nanobind::is_arithmetic())
      .value("ALPHA", MonoclinicAngleType::Alpha)
      .value("BETA", MonoclinicAngleType::Beta)
      .value("GAMMA", MonoclinicAngleType::Gamma);

  // The settings of a barostat, as Thermostat holds those of a thermostat: 'System.set_thermobarostat' builds the
  // barostat of the system from them, with the temperature, the pressure and the time step of the system.
  nanobind::class_<Thermobarostat>(m, "Thermobarostat")
      .def(
          "__init__",
          [](Thermobarostat* barostat, MolecularDynamicsEnsemble ensemble, std::size_t barostat_chain_length,
             std::size_t number_of_yoshida_suzuki_steps, double time_scale_parameter,
             std::size_t number_of_respa_steps, BarostatCoupling coupling, std::optional<CellMinimizationType> cell_type,
             MonoclinicAngleType monoclinic_angle)
          {
            const bool isotropic = molecularDynamicsUsesIsotropicBarostat(ensemble);
            const bool flexible = molecularDynamicsUsesFlexibleBarostat(ensemble);
            if (!isotropic && !flexible)
            {
              throw std::invalid_argument("Thermobarostat: the ensemble has no barostat (NPT, NPTPR, MUPT or MUPTPR)");
            }
            if (barostat_chain_length == 0 || number_of_respa_steps == 0 || number_of_yoshida_suzuki_steps == 0)
            {
              throw std::invalid_argument("Thermobarostat: the chain length and the numbers of steps must be positive");
            }
            if (!(time_scale_parameter > 0.0))
            {
              throw std::invalid_argument("Thermobarostat: time_scale_parameter must be positive");
            }
            // an isotropic barostat scales the volume; a flexible one changes the cell, fully unless a type is given
            CellMinimizationType cellType = CellMinimizationType::Isotropic;
            if (flexible)
            {
              cellType = cell_type.value_or(CellMinimizationType::Regular);
              if (cellType == CellMinimizationType::Fixed) cellType = CellMinimizationType::Regular;
            }
            new (barostat) Thermobarostat();
            barostat->ensemble = ensemble;
            barostat->cellType = cellType;
            barostat->monoclinicAngle = monoclinic_angle;
            barostat->chainLength = barostat_chain_length;
            barostat->numberOfYoshidaSuzukiSteps = number_of_yoshida_suzuki_steps;
            barostat->timeScaleParameterBarostat = time_scale_parameter;
            barostat->numberOfRespaSteps = number_of_respa_steps;
            barostat->coupling = coupling;
          },
          nanobind::arg("ensemble"), nanobind::arg("barostat_chain_length") = 5,
          nanobind::arg("number_of_yoshida_suzuki_steps") = 5, nanobind::arg("time_scale_parameter") = 1.0,
          nanobind::arg("number_of_respa_steps") = 5, nanobind::arg("coupling") = BarostatCoupling::Molecular,
          nanobind::arg("cell_type") = nanobind::none(), nanobind::arg("monoclinic_angle") = MonoclinicAngleType::Beta)
      .def_ro("ensemble", &Thermobarostat::ensemble)
      .def_ro("cell_type", &Thermobarostat::cellType)
      .def_ro("coupling", &Thermobarostat::coupling)
      .def_ro("barostat_chain_length", &Thermobarostat::chainLength)
      .def_ro("number_of_yoshida_suzuki_steps", &Thermobarostat::numberOfYoshidaSuzukiSteps)
      .def_ro("number_of_respa_steps", &Thermobarostat::numberOfRespaSteps)
      .def_ro("time_scale_parameter", &Thermobarostat::timeScaleParameterBarostat);

  // A reaction of the reaction-ensemble moves: how many molecules of each component it takes and gives. It is
  // driven by one reaction move; without one it takes the reaction move of its system that has a probability.
  nanobind::class_<Reaction>(m, "Reaction")
      .def(
          "__init__",
          [](Reaction* reaction, std::vector<std::size_t> reactants, std::vector<std::size_t> products,
             std::optional<Move::Types> move, double maximum_lambda_change, double maximum_lambda_change_products,
             double lambda_switch_point)
          {
            if (move.has_value() && *move != Move::Types::ReactionCBMC && *move != Move::Types::ReactionCFCMC &&
                *move != Move::Types::ReactionCBCFCMC && *move != Move::Types::ReactionConventionalCFCMC &&
                *move != Move::Types::ReactionConventionalCBCFCMC)
            {
              throw std::invalid_argument("Reaction: the move is one of the reaction moves of Move.Types");
            }
            if (reactants.size() != products.size())
            {
              throw std::invalid_argument("Reaction: reactants and products have an entry for every component");
            }
            if (std::ranges::all_of(reactants, [](std::size_t n) { return n == 0; }) ||
                std::ranges::all_of(products, [](std::size_t n) { return n == 0; }))
            {
              throw std::invalid_argument("Reaction: a reaction takes at least one molecule and gives at least one");
            }
            if (!(maximum_lambda_change > 0.0 && maximum_lambda_change <= 1.0) ||
                !(maximum_lambda_change_products > 0.0 && maximum_lambda_change_products <= 1.0) ||
                !(lambda_switch_point >= 0.0 && lambda_switch_point <= 1.0))
            {
              throw std::invalid_argument("Reaction: the maximum lambda changes are in (0, 1] and the switch point "
                                          "in [0, 1]");
            }
            new (reaction) Reaction(0, std::move(reactants), std::move(products));
            reaction->reactionMove = move.value_or(reactionMoveNotGiven);
            reaction->maximumLambdaChange = maximum_lambda_change;
            reaction->maximumLambdaChangeProducts = maximum_lambda_change_products;
            reaction->lambdaSwitchPoint = lambda_switch_point;
          },
          nanobind::arg("reactants"), nanobind::arg("products"), nanobind::arg("move") = nanobind::none(),
          nanobind::arg("maximum_lambda_change") = 0.3, nanobind::arg("maximum_lambda_change_products") = 0.3,
          nanobind::arg("lambda_switch_point") = 0.3,
          "The numbers of molecules of each component the reaction takes ('reactants') and gives ('products'), in "
          "the order of the components of the system.")
      .def_ro("id", &Reaction::id)
      .def_ro("reactants", &Reaction::reactantStoichiometry)
      .def_ro("products", &Reaction::productStoichiometry)
      .def_prop_ro(
          "move",
          [](const Reaction& reaction) -> std::optional<Move::Types>
          {
            if (reaction.reactionMove == reactionMoveNotGiven) return std::nullopt;
            return reaction.reactionMove;
          },
          "The reaction move that drives the reaction.")
      .def_ro("maximum_lambda_change", &Reaction::maximumLambdaChange)
      .def_ro("maximum_lambda_change_products", &Reaction::maximumLambdaChangeProducts)
      .def_ro("lambda_switch_point", &Reaction::lambdaSwitchPoint)
      .def_ro("current_lambda", &Reaction::currentLambda, "The coupling of the fractional molecules of the reaction.")
      .def_prop_ro(
          "lambda_histogram", [](Reaction& reaction) -> PropertyLambdaProbabilityHistogram& { return reaction.lambda; },
          nanobind::rv_policy::reference_internal,
          "The lambda histogram with the fractional molecules on the side of the reactants.")
      .def_prop_ro(
          "lambda_product_side_histogram",
          [](Reaction& reaction) -> PropertyLambdaProbabilityHistogram& { return reaction.lambdaProductSide; },
          nanobind::rv_policy::reference_internal,
          "The lambda histogram with the fractional molecules on the side of the products.")
      .def("__repr__", [](const Reaction& reaction)
           {
             if (reaction.reactionMove == reactionMoveNotGiven)
             {
               return std::format("Reaction(reactants={}, products={})", reaction.reactantStoichiometry,
                                  reaction.productStoichiometry);
             }
             return reaction.printStatus();
           });

  // A histogram of a property comes as the centers of its bins, the probability density, and the error of it.
  nanobind::class_<PropertyMoleculeProperties>(m, "PropertyMoleculeProperties")
      .def_ro("number_of_bonds", &PropertyMoleculeProperties::numberOfBondsPerComponent,
              "The number of bonds of a molecule, per component.")
      .def_ro("number_of_bends", &PropertyMoleculeProperties::numberOfBendsPerComponent)
      .def_ro("number_of_torsions", &PropertyMoleculeProperties::numberOfTorsionsPerComponent)
      .def(
          "bond_histogram",
          [](const PropertyMoleculeProperties& p, std::size_t component, std::size_t bond)
          {
            checkSampledComponent(p, component, "bond_histogram");
            if (bond >= p.numberOfBondsPerComponent[component]) throw nanobind::index_error("bond_histogram: no such bond");
            return p.result(p.bondHistogram, component, bond, p.deltaBondPerComponent[component], 0.0);
          },
          nanobind::arg("component"), nanobind::arg("bond"),
          "The histogram of the length of a bond of the component [Angstrom], by the index of its bond potential.")
      .def(
          "bend_histogram",
          [](const PropertyMoleculeProperties& p, std::size_t component, std::size_t bend)
          {
            checkSampledComponent(p, component, "bend_histogram");
            if (bend >= p.numberOfBendsPerComponent[component]) throw nanobind::index_error("bend_histogram: no such bend");
            return p.result(p.bendHistogram, component, bend, p.deltaBendPerComponent[component], 0.0);
          },
          nanobind::arg("component"), nanobind::arg("bend"), "The histogram of a bend angle [degrees].")
      .def(
          "torsion_histogram",
          [](const PropertyMoleculeProperties& p, std::size_t component, std::size_t torsion)
          {
            checkSampledComponent(p, component, "torsion_histogram");
            if (torsion >= p.numberOfTorsionsPerComponent[component])
            {
              throw nanobind::index_error("torsion_histogram: no such torsion");
            }
            return p.result(p.torsionHistogram, component, torsion, p.deltaTorsionPerComponent[component],
                            -0.5 * p.torsionRange);
          },
          nanobind::arg("component"), nanobind::arg("torsion"),
          "The histogram of a torsion angle, from -180 to 180 [degrees].")
      .def(
          "end_to_end_histogram",
          [](const PropertyMoleculeProperties& p, std::size_t component)
          {
            checkSampledComponent(p, component, "end_to_end_histogram");
            if (!p.endToEndAtomsPerComponent[component].has_value())
            {
              throw nanobind::index_error("end_to_end_histogram: the component has no end-to-end atoms");
            }
            return p.result(p.endToEndHistogram, component, 0, p.deltaEndToEndPerComponent[component], 0.0);
          },
          nanobind::arg("component"), "The histogram of the end-to-end distance [Angstrom].");

  nanobind::class_<PropertyMoleculeShape> moleculeShape(m, "PropertyMoleculeShape");
  nanobind::enum_<PropertyMoleculeShape::Moment>(moleculeShape, "Moment", nanobind::is_arithmetic())
      .value("RADIUS_OF_GYRATION", PropertyMoleculeShape::RadiusOfGyration)
      .value("RADIUS_OF_GYRATION_SQUARED", PropertyMoleculeShape::RadiusOfGyrationSquared)
      .value("RADIUS_OF_GYRATION_FOURTH", PropertyMoleculeShape::RadiusOfGyrationFourth)
      .value("LAMBDA_1", PropertyMoleculeShape::Lambda1)
      .value("LAMBDA_2", PropertyMoleculeShape::Lambda2)
      .value("LAMBDA_3", PropertyMoleculeShape::Lambda3)
      .value("ASPHERICITY", PropertyMoleculeShape::Asphericity)
      .value("ACYLINDRICITY", PropertyMoleculeShape::Acylindricity)
      .value("SHAPE_ANISOTROPY", PropertyMoleculeShape::ShapeAnisotropy)
      .value("SECOND_INVARIANT", PropertyMoleculeShape::SecondInvariant)
      .value("PROLATENESS", PropertyMoleculeShape::Prolateness)
      .value("INVERSE_HYDRODYNAMIC_RADIUS", PropertyMoleculeShape::InverseHydrodynamicRadius)
      .value("END_TO_END_SQUARED", PropertyMoleculeShape::EndToEndSquared)
      .value("TENSOR_XX", PropertyMoleculeShape::TensorXX)
      .value("TENSOR_YY", PropertyMoleculeShape::TensorYY)
      .value("TENSOR_ZZ", PropertyMoleculeShape::TensorZZ)
      .value("TENSOR_XY", PropertyMoleculeShape::TensorXY)
      .value("TENSOR_XZ", PropertyMoleculeShape::TensorXZ)
      .value("TENSOR_YZ", PropertyMoleculeShape::TensorYZ);
  moleculeShape
      .def(
          "result",
          [](const PropertyMoleculeShape& p, std::size_t component, PropertyMoleculeShape::Moment moment)
          {
            checkSampledComponent(p, component, "result");
            return p.momentStatistics(component, moment);
          },
          nanobind::arg("component"), nanobind::arg("moment"),
          "The average of a moment over the molecules of the component and its error: lengths in Angstrom, the "
          "eigenvalues of the gyration tensor, the asphericity and the acylindricity in Angstrom^2.")
      .def(
          "radius_of_gyration_histogram",
          [](const PropertyMoleculeShape& p, std::size_t component)
          {
            checkSampledComponent(p, component, "radius_of_gyration_histogram");
            return p.result(p.radiusOfGyrationHistogram, component, p.deltaRadiusOfGyrationPerComponent[component],
                            0.0);
          },
          nanobind::arg("component"), "The histogram of the radius of gyration [Angstrom].")
      .def(
          "shape_anisotropy_histogram",
          [](const PropertyMoleculeShape& p, std::size_t component)
          {
            checkSampledComponent(p, component, "shape_anisotropy_histogram");
            return p.result(p.shapeAnisotropyHistogram, component, p.deltaShapeAnisotropyPerComponent[component], 0.0);
          },
          nanobind::arg("component"), "The histogram of the relative shape anisotropy, from 0 to 1.")
      .def(
          "prolateness_histogram",
          [](const PropertyMoleculeShape& p, std::size_t component)
          {
            checkSampledComponent(p, component, "prolateness_histogram");
            return p.result(p.prolatenessHistogram, component, p.deltaProlatenessPerComponent[component],
                            p.prolatenessLowerLimit);
          },
          nanobind::arg("component"), "The histogram of the prolateness, from -1/4 to 2.");

  nanobind::class_<PropertyMoleculeBackbone> moleculeBackbone(m, "PropertyMoleculeBackbone");
  nanobind::enum_<PropertyMoleculeBackbone::Moment>(moleculeBackbone, "Moment", nanobind::is_arithmetic())
      .value("BOND_LENGTH", PropertyMoleculeBackbone::BondLength)
      .value("BOND_LENGTH_SQUARED", PropertyMoleculeBackbone::BondLengthSquared)
      .value("END_TO_END_SQUARED", PropertyMoleculeBackbone::EndToEndSquared)
      .value("PROJECTION", PropertyMoleculeBackbone::Projection)
      .value("RADIUS_OF_GYRATION_SQUARED", PropertyMoleculeBackbone::RadiusOfGyrationSquared);
  nanobind::class_<PropertyMoleculeBackbone::Averages>(moleculeBackbone, "Averages")
      .def_ro("internal_distance_squared", &PropertyMoleculeBackbone::Averages::internalDistanceSquared,
              "The mean squared distance of backbone beads k bonds apart [Angstrom^2], for k = 1, 2, ...")
      .def_ro("bond_correlation", &PropertyMoleculeBackbone::Averages::bondCorrelation,
              "The correlation of the directions of backbone bonds k bonds apart, for k = 0, 1, ...")
      .def_ro("form_factor", &PropertyMoleculeBackbone::Averages::formFactor,
              "The single-chain form factor at the wave vectors of 'wave_vectors'.");
  moleculeBackbone
      .def(
          "result",
          [](const PropertyMoleculeBackbone& p, std::size_t component)
          {
            checkSampledComponent(p, component, "result");
            return p.overallAverages(component);
          },
          nanobind::arg("component"), "The averages over the molecules of the component and over the run.")
      .def(
          "moment",
          [](const PropertyMoleculeBackbone& p, std::size_t component, PropertyMoleculeBackbone::Moment moment)
          {
            checkSampledComponent(p, component, "moment");
            return p.statistics(component, [moment](const PropertyMoleculeBackbone::Averages& averages)
                                { return averages.moments[moment]; });
          },
          nanobind::arg("component"), nanobind::arg("moment"), "The average of a moment and its error.")
      .def(
          "persistence_length_from_projection",
          [](const PropertyMoleculeBackbone& p, std::size_t component)
          {
            checkSampledComponent(p, component, "persistence_length_from_projection");
            return p.statistics(component, &PropertyMoleculeBackbone::persistenceLengthFromProjection);
          },
          nanobind::arg("component"),
          "The persistence length from the projection of the chain on its end bonds, and its error [Angstrom].")
      .def(
          "persistence_length_from_fit",
          [](const PropertyMoleculeBackbone& p, std::size_t component)
          {
            checkSampledComponent(p, component, "persistence_length_from_fit");
            return p.statistics(component, &PropertyMoleculeBackbone::persistenceLengthFromFit);
          },
          nanobind::arg("component"),
          "The persistence length from the decay of the bond correlation, and its error [Angstrom].")
      .def(
          "flory_exponent",
          [](const PropertyMoleculeBackbone& p, std::size_t component)
          {
            checkSampledComponent(p, component, "flory_exponent");
            return p.statistics(component, &PropertyMoleculeBackbone::floryExponent);
          },
          nanobind::arg("component"), "The Flory exponent from the internal distances, and its error.")
      .def(
          "wave_vectors",
          [](const PropertyMoleculeBackbone& p, std::size_t component)
          {
            checkSampledComponent(p, component, "wave_vectors");
            return p.waveVectorsPerComponent[component];
          },
          nanobind::arg("component"), "The wave vectors of the form factor [1/Angstrom].")
      .def(
          "backbone",
          [](const PropertyMoleculeBackbone& p, std::size_t component)
          {
            checkSampledComponent(p, component, "backbone");
            return p.backbonePerComponent[component];
          },
          nanobind::arg("component"), "The atoms of the backbone of the component.")
      .def(
          "contour_length",
          [](const PropertyMoleculeBackbone& p, std::size_t component)
          {
            checkSampledComponent(p, component, "contour_length");
            return p.contourLengthPerComponent[component];
          },
          nanobind::arg("component"), "The contour length of the backbone [Angstrom].");

  nanobind::class_<EndToEndAutoCorrelationFunctionData>(m, "EndToEndAutoCorrelationFunctionData")
      .def_ro("time", &EndToEndAutoCorrelationFunctionData::time, "The lag: picoseconds, or cycles for Monte Carlo.")
      .def_ro("acf", &EndToEndAutoCorrelationFunctionData::acf, "<R(0).R(t)> [Angstrom^2].")
      .def_ro("normalized", &EndToEndAutoCorrelationFunctionData::normalized, "<R(0).R(t)> / <R^2>.")
      .def_ro("number_of_samples", &EndToEndAutoCorrelationFunctionData::numberOfSamples)
      .def_ro("end_to_end_squared_acf", &EndToEndAutoCorrelationFunctionData::endToEndSquaredACF)
      .def_ro("radius_of_gyration_squared_acf", &EndToEndAutoCorrelationFunctionData::radiusOfGyrationSquaredACF);

  nanobind::class_<RelaxationTimeEstimates>(m, "RelaxationTimeEstimates")
      .def_ro("integrated", &RelaxationTimeEstimates::integrated)
      .def_ro("integrated_is_lower_bound", &RelaxationTimeEstimates::integratedIsLowerBound)
      .def_ro("one_over_e", &RelaxationTimeEstimates::oneOverE)
      .def_ro("exponential_fit", &RelaxationTimeEstimates::exponentialFit);

  nanobind::class_<ScalarRelaxation>(m, "ScalarRelaxation")
      .def_ro("available", &ScalarRelaxation::available)
      .def_ro("mean", &ScalarRelaxation::mean)
      .def_ro("variance", &ScalarRelaxation::variance)
      .def_ro("times", &ScalarRelaxation::times);

  nanobind::class_<EndToEndRelaxationTimes>(m, "EndToEndRelaxationTimes")
      .def_ro("mean_squared_end_to_end", &EndToEndRelaxationTimes::meanSquaredEndToEnd, "<R^2> [Angstrom^2].")
      .def_ro("longest_lag", &EndToEndRelaxationTimes::longestLag)
      .def_ro("integrated", &EndToEndRelaxationTimes::integrated, "The integral of the normalized function.")
      .def_ro("integrated_is_lower_bound", &EndToEndRelaxationTimes::integratedIsLowerBound)
      .def_ro("one_over_e", &EndToEndRelaxationTimes::oneOverE, "The lag at which the function has decayed to 1/e.")
      .def_ro("exponential_fit", &EndToEndRelaxationTimes::exponentialFit)
      .def_ro("end_to_end_squared", &EndToEndRelaxationTimes::endToEndSquared)
      .def_ro("radius_of_gyration_squared", &EndToEndRelaxationTimes::radiusOfGyrationSquared);

  nanobind::class_<PropertyEndToEndAutoCorrelationFunction>(m, "PropertyEndToEndAutoCorrelationFunction")
      .def("has_data", &PropertyEndToEndAutoCorrelationFunction::hasData, nanobind::arg("component"))
      .def(
          "result",
          [](const PropertyEndToEndAutoCorrelationFunction& p, std::size_t component)
          {
            if (!p.hasData(component)) return std::vector<EndToEndAutoCorrelationFunctionData>{};
            return p.result(component);
          },
          nanobind::arg("component"), "The autocorrelation function of the component, per lag.")
      .def(
          "relaxation_times",
          [](const PropertyEndToEndAutoCorrelationFunction& p, std::size_t component)
          {
            if (!p.hasData(component))
            {
              throw nanobind::index_error("relaxation_times: the component has no data");
            }
            return p.relaxationTimes(component);
          },
          nanobind::arg("component"), "The relaxation times of the end-to-end vector, in the unit of the lags.")
      .def_prop_ro("time_unit", &PropertyEndToEndAutoCorrelationFunction::timeUnit,
                   "The unit of the lags: 'ps' for molecular dynamics, 'cycles' for Monte Carlo.");

  // The elastic constants from the fluctuations of the stress, in molecular dynamics at constant volume and
  // temperature. The result is what the output file holds: the tensors in Voigt order and in GPa, and what follows
  // from them (compliance, moduli, Poisson ratios) when the stiffness can be inverted.
  nanobind::class_<PropertyElasticConstantsFluctuation>(m, "PropertyElasticConstantsFluctuation")
      .def(nanobind::init<std::size_t>(), nanobind::arg("number_of_blocks") = 5)
      .def(
          "result",
          [](const PropertyElasticConstantsFluctuation& p, double eigenvalue_tolerance)
          {
            return nanobind::module_::import_("json").attr("loads")(p.jsonAveragesStatistics(eigenvalue_tolerance).dump());
          },
          nanobind::arg("eigenvalue_tolerance") = 1.0e-8,
          "The averages as a dictionary: the stress, the Born, kinetic and fluctuation terms and the stiffness in "
          "Voigt order [GPa], with the moduli derived from them.")
      .def("write_averages_statistics", &PropertyElasticConstantsFluctuation::writeAveragesStatistics,
           nanobind::arg("eigenvalue_tolerance") = 1.0e-8);

  // the times a radical waits between reactions, of the polymerization driver
  nanobind::class_<PropertyRadicalWaitingTimes> radicalWaitingTimes(m, "PropertyRadicalWaitingTimes");
  nanobind::class_<PropertyRadicalWaitingTimes::Histogram>(radicalWaitingTimes, "Histogram")
      .def_ro("bins_per_decade", &PropertyRadicalWaitingTimes::Histogram::binsPerDecade)
      .def_ro("counts", &PropertyRadicalWaitingTimes::Histogram::counts)
      .def_ro("count", &PropertyRadicalWaitingTimes::Histogram::count)
      .def_ro("maximum", &PropertyRadicalWaitingTimes::Histogram::maximum)
      .def("mean", &PropertyRadicalWaitingTimes::Histogram::mean)
      .def("geometric_mean", &PropertyRadicalWaitingTimes::Histogram::geometricMean)
      .def("median", &PropertyRadicalWaitingTimes::Histogram::median)
      .def(
          "bin_edges",
          [](const PropertyRadicalWaitingTimes::Histogram& histogram)
          {
            std::vector<double> edges;
            for (std::size_t bin = 0; bin <= histogram.numberOfBins(); ++bin) edges.push_back(histogram.lowerEdge(bin));
            return edges;
          },
          "The edges of the bins [ps], logarithmically spaced.");
  radicalWaitingTimes
      .def(
          "__init__",
          [](PropertyRadicalWaitingTimes* property, std::size_t bins_per_decade, std::size_t write_every)
          {
            if (bins_per_decade == 0 || write_every == 0)
            {
              throw std::invalid_argument("PropertyRadicalWaitingTimes: bins_per_decade and write_every must be "
                                          "positive");
            }
            new (property) PropertyRadicalWaitingTimes(bins_per_decade, write_every);
          },
          nanobind::arg("bins_per_decade") = 10, nanobind::arg("write_every") = 5000)
      .def_ro("bins_per_decade", &PropertyRadicalWaitingTimes::binsPerDecade)
      .def_ro("write_every", &PropertyRadicalWaitingTimes::writeEvery)
      .def("closed", &PropertyRadicalWaitingTimes::closed, "The histogram of the intervals that a reaction ended.")
      .def_ro("by_reaction", &PropertyRadicalWaitingTimes::byReaction,
              "The histograms of the intervals, per reaction that ended them.");

  // the size and shape of the chains of the polymerization driver, per chain length
  nanobind::class_<PropertyPolymerChains>(m, "PropertyPolymerChains")
      .def(
          "__init__",
          [](PropertyPolymerChains* property, std::size_t sample_every, std::size_t write_every,
             const std::string& sampling, std::size_t lag_steps, std::size_t conversion_bins,
             std::size_t minimum_length, bool mass_weighted)
          {
            PropertyPolymerChains::Sampling mode{};
            if (sampling == "periodic")
            {
              mode = PropertyPolymerChains::Sampling::Periodic;
            }
            else if (sampling == "before_reaction")
            {
              mode = PropertyPolymerChains::Sampling::BeforeReaction;
            }
            else
            {
              throw std::invalid_argument("PropertyPolymerChains: sampling must be 'periodic' or 'before_reaction'");
            }
            new (property) PropertyPolymerChains(sample_every, write_every, mode, lag_steps, conversion_bins,
                                                 minimum_length, mass_weighted);
          },
          nanobind::arg("sample_every") = 1000, nanobind::arg("write_every") = 5000,
          nanobind::arg("sampling") = "periodic", nanobind::arg("lag_steps") = 0,
          nanobind::arg("conversion_bins") = 10, nanobind::arg("minimum_length") = 2,
          nanobind::arg("mass_weighted") = false)
      .def_ro("sample_every", &PropertyPolymerChains::sampleEvery)
      .def_ro("write_every", &PropertyPolymerChains::writeEvery)
      .def_prop_ro("sampling",
                   [](const PropertyPolymerChains& p)
                   {
                     return p.sampling == PropertyPolymerChains::Sampling::Periodic ? "periodic" : "before_reaction";
                   })
      .def_ro("lag_steps", &PropertyPolymerChains::lagSteps)
      .def_ro("conversion_bins", &PropertyPolymerChains::numberOfConversionBins)
      .def_ro("minimum_length", &PropertyPolymerChains::minimumLength)
      .def_ro("mass_weighted", &PropertyPolymerChains::massWeighted)
      .def_ro("sampling_events", &PropertyPolymerChains::numberOfSamplingEvents)
      .def_ro("skipped_ramping", &PropertyPolymerChains::skippedRamping)
      .def_ro("skipped_spanning", &PropertyPolymerChains::skippedSpanning)
      .def_prop_ro(
          "moments",
          [](const PropertyPolymerChains&)
          {
            return std::vector<std::string>(PropertyPolymerChains::momentNames.begin(),
                                            PropertyPolymerChains::momentNames.end());
          },
          "The names of the accumulated per-chain values, in the order of 'sum' and 'sumOfSquares'.")
      .def(
          "entries",
          [](const PropertyPolymerChains& p)
          {
            nanobind::list entries;
            for (const auto& [key, accumulator] : p.accumulators)
            {
              nanobind::dict entry;
              entry["length"] = key.length;
              entry["topology"] = std::string(PropertyPolymerChains::topologyName(key.topology));
              entry["living"] = key.living;
              entry["conversionBin"] = key.conversionBin;
              entry["count"] = accumulator.count;
              entry["sum"] = std::vector<double>(accumulator.sum.begin(), accumulator.sum.end());
              entry["sumOfSquares"] = std::vector<double>(accumulator.sumOfSquares.begin(), accumulator.sumOfSquares.end());
              entries.append(entry);
            }
            return entries;
          },
          "The raw sums per key (length, topology, living, conversion bin), as the 'entries' of the JSON output; "
          "combine runs with raspalib.polymer_chains.combine.");

  // writes the configuration as LAMMPS data files ('s<system>.data' and its tables) every 'sample_every' cycles
  nanobind::class_<WriteLammpsData>(m, "WriteLammpsData")
      .def(
          "__init__",
          [](WriteLammpsData* writer, std::optional<std::size_t> system_id, std::size_t sample_every)
          {
            if (sample_every == 0) throw std::invalid_argument("WriteLammpsData: sample_every must be positive");
            if (system_id.has_value())
            {
              new (writer) WriteLammpsData(*system_id, sample_every);
              return;
            }
            // the files are named when the driver gives the index of the system
            new (writer) WriteLammpsData();
            writer->sampleEvery = sample_every;
            writer->systemId = systemIdNotGiven;
          },
          nanobind::arg("system_id") = nanobind::none(), nanobind::arg("sample_every") = 10)
      .def_prop_ro(
          "system_id",
          [](const WriteLammpsData& writer) -> std::optional<std::size_t>
          {
            if (writer.systemId == systemIdNotGiven) return std::nullopt;
            return writer.systemId;
          },
          "The index in the names of the files; None until a driver gives the index of the system.")
      .def_ro("sample_every", &WriteLammpsData::sampleEvery);

  nanobind::class_<PropertyTemperature>(m, "PropertyTemperature")
      .def("result", [](const PropertyTemperature& p) { return p.average(); },
           "The average temperature and its error [K].");

  nanobind::class_<PropertySimulationBox>(m, "PropertySimulationBox")
      .def("result", [](const PropertySimulationBox& p) { return p.average(); },
           "The average simulation box and its error.");

  // the objects given to the setters of the properties of a system, per system (see 'attachGivenProperty')
  m.attr("_given_properties") = nanobind::module_::import_("weakref").attr("WeakKeyDictionary")();
  givenPropertiesTable = m.attr("_given_properties");
  copiesOfGivenProperties.clear();

  nanobind::class_<System> system(m, "System", nanobind::is_weak_referenceable());
  system
      .def("__init__",
           [](System* system, ForceField force_field, std::optional<SimulationBox> simulation_box,
              bool has_external_field, double external_temperature, std::optional<double> external_pressure,
              double helium_void_fraction, std::optional<Framework> framework_components,
              std::vector<Component> components, std::vector<std::vector<double3>> initial_positions,
              std::vector<std::size_t> initial_number_of_molecules, std::optional<std::size_t> blocks_given,
              MCMoveProbabilities system_probabilities, std::vector<Reaction> reactions)
           {
             // The number of blocks is that of the components when it is not given (they then have to agree), and
             // the components are sized for it when it is.
             if (blocks_given == 0uz) throw std::invalid_argument("System: number_of_blocks must be positive");
             std::size_t number_of_blocks = blocks_given.value_or(5);
             if (!blocks_given.has_value() && !components.empty())
             {
               number_of_blocks = components.front().numberOfBlocks();
               for (const Component& component : components)
               {
                 if (component.numberOfBlocks() != number_of_blocks)
                 {
                   throw std::invalid_argument("System: the components differ in their number of blocks: give "
                                               "'number_of_blocks' to the system");
                 }
               }
             }
             for (Component& component : components) component.setNumberOfBlocks(number_of_blocks);

             // the partition function of a species is that of the temperature of the system
             for (Component& component : components) component.updateLnPartitionFunction(external_temperature);

             new (system) System(force_field, simulation_box, has_external_field, external_temperature,
                                 external_pressure, helium_void_fraction, framework_components, components,
                                 initial_positions, initial_number_of_molecules, number_of_blocks,
                                 system_probabilities);
             try
             {
               setReactions(*system, std::move(reactions));
             }
             catch (const std::runtime_error& error)
             {
               system->~System();
               throw std::invalid_argument(std::format("System: {}", error.what()));
             }
             catch (...)
             {
               system->~System();
               throw;
             }
           },
           nanobind::arg("force_field"), nanobind::arg("simulation_box") = nanobind::none(),
           nanobind::arg("has_external_field") = false, nanobind::arg("external_temperature") = 298.0,
           nanobind::arg("external_pressure") = nanobind::none(), nanobind::arg("helium_void_fraction") = 0.0,
           nanobind::arg("framework_components") = nanobind::none(),
           nanobind::arg("components") = std::vector<Component>(),
           nanobind::arg("initial_positions") = std::vector<std::vector<double3>>(),
           nanobind::arg("initial_number_of_molecules") = std::vector<std::size_t>(),
           nanobind::arg("number_of_blocks") = nanobind::none(),
           nanobind::arg("system_probabilities") = MCMoveProbabilities(),
           nanobind::arg("reactions") = std::vector<Reaction>())
      .def_prop_ro("reactions",
                   [](nanobind::handle self)
                   {
                     // the live objects, each held by the parent
                     nanobind::list list;
                     for (Reaction& item : nanobind::cast<System&>(self).reactions.list)
                       list.append(nanobind::cast(&item, nanobind::rv_policy::reference_internal, self));
                     return list;
                   })
      .def_prop_ro(
          "number_of_blocks", [](const System& s) { return s.numberOfBlocks(); },
          "The number of blocks the averages of the system are kept in. A driver sizes its systems for its own.")
      .def_ro("mc_moves_probabilities", &System::mc_moves_probabilities,
              "The probabilities of the moves of the system (volume, Gibbs volume, reactions).")
      // (the member takes a pointer to an overlap flag, which was bound as a boolean that had to be given)
      .def("compute_total_energies", [](System& system) { return system.computeTotalEnergies(); },
           "The energies of the configuration as it is, in internal units of 10 J/mol.")
      .def("framework_mass", &System::frameworkMass)
      .def_ro("simulation_box", &System::simulationBox)
      .def_ro("input_pressure", &System::input_pressure)
      .def_prop_ro("components",
                   [](nanobind::handle self)
                   {
                     // the live objects, each held by the parent: a copied vector lost every change made through it
                     nanobind::list list;
                     for (Component& item : nanobind::cast<System&>(self).components)
                       list.append(nanobind::cast(&item, nanobind::rv_policy::reference_internal, self));
                     return list;
                   })
      .def_ro("number_of_molecules_per_component", &System::numberOfMoleculesPerComponent)
      .def_ro("mc_moves_statistics", &System::mc_moves_statistics)
      .def_ro("loadings", &System::loadings)
      .def_rw("average_loadings", &System::averageLoadings)
      .def_ro("average_energies", &System::averageEnergies)
      .def_ro("average_pressure", &System::averagePressure)
      .def_ro("thermostat", &System::thermostat)
      .def(
          "set_thermostat",
          [](System& system, const std::optional<Thermostat>& thermostat)
          {
            system.setThermostat(thermostat);
            // the ensemble says what the thermostat makes of it, as it does for an input file: constant energy
            // becomes constant temperature, and back when the thermostat is removed
            if (thermostat.has_value() && system.molecularDynamicsEnsemble == MolecularDynamicsEnsemble::NVE)
            {
              system.molecularDynamicsEnsemble = MolecularDynamicsEnsemble::NVT;
            }
            if (!thermostat.has_value() && system.molecularDynamicsEnsemble == MolecularDynamicsEnsemble::NVT)
            {
              system.molecularDynamicsEnsemble = MolecularDynamicsEnsemble::NVE;
            }
          },
          nanobind::arg("thermostat") = nanobind::none())
      .def_ro("property_density_grid", &System::propertyDensityGrid)
      .def_ro("average_energy_histogram", &System::averageEnergyHistogram)
      .def_ro("property_number_of_molecules_histogram", &System::averageNumberOfMoleculesHistogram)
      .def_ro("sample_pdb_movie", &System::samplePDBMovie)
      .def_ro("property_conventional_rdf", &System::propertyConventionalRadialDistributionFunction)
      .def_ro("property_msd", &System::propertyMSD)
      .def_ro("property_vacf", &System::propertyVACF)
      .def_ro("property_rdf", &System::propertyRadialDistributionFunction)
      .def_ro("property_number_of_molecules_evolution", &System::propertyNumberOfMoleculesEvolution)
      .def_ro("property_molecule_properties", &System::propertyMoleculeProperties,
              "The histograms of the molecules of the components that have 'molecule_properties_settings'.")
      .def_ro("property_molecule_shape", &System::propertyMoleculeShape)
      .def_ro("property_molecule_backbone", &System::propertyMoleculeBackbone)
      .def_ro("property_end_to_end_acf", &System::propertyEndToEndACF)
      .def_ro("property_elastic_constants_fluctuation", &System::propertyElasticConstantsFluctuation);
  // bound by hand for its 'sample_every', and registered like the setters of bindGivenProperty so that a driver
  // copies the sampled property back to the object that was given
  copiesOfGivenProperties.emplace_back("set_property_elastic_constants_fluctuation",
                                       &copyPropertyTo<&System::propertyElasticConstantsFluctuation>);
  system
      .def(
          "set_property_elastic_constants_fluctuation",
          [](nanobind::handle self, nanobind::handle given, std::size_t sample_every)
          {
            const char* name = "set_property_elastic_constants_fluctuation";
            System& system = nanobind::cast<System&>(self);
            const auto property = givenPropertyOrNone<&System::propertyElasticConstantsFluctuation>(given, name);
            if (!property.has_value())
            {
              system.propertyElasticConstantsFluctuation.reset();
              attachGivenProperty<&System::propertyElasticConstantsFluctuation>(self, name, given);
              return;
            }
            if (sample_every == 0) throw std::invalid_argument("sample_every must be positive");
            if (system.hasExternalField || system.forceField.computePolarization)
            {
              throw std::invalid_argument("the elastic constants from fluctuations do not support an external field "
                                          "or polarization");
            }
            // the blocks are those of the system, whatever the given property was made with (the sampling indexes
            // the blocks of the system)
            system.propertyElasticConstantsFluctuation =
                PropertyElasticConstantsFluctuation(system.averageEnergies.numberOfBlocks);
            system.elasticConstantsSampleEvery = sample_every;
            attachGivenProperty<&System::propertyElasticConstantsFluctuation>(self, name, given);
          },
          nanobind::arg("property") = nanobind::none(), nanobind::arg("sample_every") = 100,
          "Samples the elastic constants from the fluctuations of the stress, every 'sample_every' steps of "
          "molecular dynamics at constant volume and temperature.")
      .def_ro("property_radical_waiting_times", &System::propertyRadicalWaitingTimes)
      .def_ro("property_polymer_chains", &System::propertyPolymerChains)
      .def_ro("write_lammps_data", &System::writeLammpsData)
      .def_ro("property_volume_evolution", &System::propertyVolumeEvolution)
      .def_ro("property_conserved_energy_evolution", &System::propertyConservedEnergyEvolution)
      .def_rw("average_enthalpies_of_adsorption", &System::averageEnthalpiesOfAdsorption)
      .def_rw("average_partial_molar_properties", &System::averagePartialMolarProperties)
      // read only: the atoms are kept in step with their dynamics, molecules and groups
      .def_ro("atom_data", &System::atomData)
      .def_rw("translational_degrees_of_freedom", &System::translationalDegreesOfFreedom)
      .def_rw("rotational_degrees_of_freedom", &System::rotationalDegreesOfFreedom)
      .def_prop_rw(
          "time_step", [](const System& system) { return system.timeStep; },
          [](System& system, double dt) { system.setTimeStep(dt); },
          "The MD time step [ps]; set, it is also that of the thermostat, the barostat and the correlation functions.")
      .def("write_mc_move_statistics", &System::writeMCMoveStatistics)
      .def_ro("temperature", &System::temperature, "The external temperature [K].")
      .def_prop_rw(
          "external_pressure_diagonal", [](const System& s) { return s.input_pressureTensorDiagonal; },
          [](System& s, double3 pressure)
          {
            s.input_pressureTensorDiagonal = pressure;
            s.pressureTensorDiagonal = pressure / Units::PressureConversionFactor;
          },
          "The external pressure along x, y and z [Pa], for the anisotropic volume move.")
      .def_rw("molecular_dynamics_ensemble", &System::molecularDynamicsEnsemble,
              "The ensemble of molecular dynamics. 'set_thermobarostat' sets it to that of the barostat.")
      .def_ro("thermobarostat", &System::thermobarostat)
      .def("set_thermobarostat", &System::setThermobarostat, nanobind::arg("thermobarostat") = nanobind::none(),
           "Builds the barostat of the system from the given settings, with the temperature, the pressure and the "
           "time step the system has at this moment; None removes it. The ensembles with a barostat are "
           "thermostatted as well: see 'set_thermostat'.")
      .def_prop_rw(
          "number_of_hybrid_mc_steps", [](const System& s) { return s.numberOfHybridMCSteps; },
          [](System& s, std::size_t steps)
          {
            if (steps == 0) throw std::invalid_argument("number_of_hybrid_mc_steps must be at least 1");
            s.numberOfHybridMCSteps = steps;
          },
          "The number of molecular-dynamics steps of a hybrid Monte Carlo move.")
      .def_ro("average_temperature", &System::averageTemperature)
      .def_ro("average_translational_temperature", &System::averageTranslationalTemperature)
      .def_ro("average_rotational_temperature", &System::averageRotationalTemperature)
      .def_ro("average_simulation_box", &System::averageSimulationBox)
      .def("__repr__", &System::repr);

  // the setters of the properties: the system builds its own from the settings of the one given, and the given one
  // holds the results when a driver returns
  bindGivenProperty<&System::propertyDensityGrid>(system, "set_property_density_grid",
                                                  &System::setPropertyDensityGrid);
  bindGivenProperty<&System::averageEnergyHistogram>(system, "set_average_energy_histogram",
                                                     &System::setAverageEnergyHistogram);
  bindGivenProperty<&System::averageNumberOfMoleculesHistogram>(system, "set_number_of_molecules_histogram",
                                                                &System::setNumberOfMoleculesHistogram);
  bindGivenProperty<&System::samplePDBMovie>(system, "set_sample_pdb_movie", &System::setSamplePDBMovie);
  bindGivenProperty<&System::propertyConventionalRadialDistributionFunction>(
      system, "set_property_conventional_rdf", &System::setPropertyConventionalRDF);
  bindGivenProperty<&System::propertyMSD>(system, "set_property_msd", &System::setPropertyMSD);
  bindGivenProperty<&System::propertyVACF>(system, "set_property_vacf", &System::setPropertyVACF);
  bindGivenProperty<&System::propertyRadialDistributionFunction>(system, "set_property_rdf",
                                                                 &System::setPropertyRDF);
  bindGivenProperty<&System::propertyNumberOfMoleculesEvolution>(system, "set_property_number_of_molecules_evolution",
                                                                 &System::setPropertyNumberOfMoleculesEvolution);
  bindGivenProperty<&System::propertyVolumeEvolution>(system, "set_property_volume_evolution",
                                                      &System::setPropertyVolumeEvolution);
  bindGivenProperty<&System::propertyConservedEnergyEvolution>(system, "set_property_conserved_energy_evolution",
                                                               &System::setPropertyConservedEnergyEvolution);
  bindGivenProperty<&System::propertyRadicalWaitingTimes>(system, "set_property_radical_waiting_times",
                                                          &System::setPropertyRadicalWaitingTimes);
  bindGivenProperty<&System::propertyPolymerChains>(system, "set_property_polymer_chains",
                                                    &System::setPropertyPolymerChains);
  bindGivenProperty<&System::writeLammpsData>(
      system, "set_write_lammps_data",
      [](System& s, const std::optional<WriteLammpsData>& writer) { s.writeLammpsData = writer; });

  // The thread pool is one for the process. It is started once: a later call with another number of threads keeps
  // the threads that are running, and one thread switches to serial work.
  m.def(
      "set_number_of_threads",
      [](std::size_t number_of_threads)
      {
        if (number_of_threads == 0) throw std::invalid_argument("set_number_of_threads: at least one thread");
        auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type>::instance();
        pool.init(number_of_threads, number_of_threads > 1 ? ThreadPool::ThreadingType::ThreadPool
                                                           : ThreadPool::ThreadingType::Serial);
      },
      nanobind::arg("number_of_threads"),
      "Sets the number of threads of the thread pool of the process, as 'NumberOfThreads' does for the "
      "command-line program.");
  m.def(
      "number_of_threads",
      []() -> std::size_t
      {
        auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type>::instance();
        // the calling thread works too
        return pool.getThreadingType() == ThreadPool::ThreadingType::ThreadPool ? pool.getThreadCount() + 1 : 1;
      },
      "The number of threads the thread pool works with.");

  bindDriver<MonteCarlo>(m, "MonteCarlo");
  bindDriver<MolecularDynamics>(m, "MolecularDynamics");

  nanobind::class_<PropertyConservedEnergyEvolution> conservedEnergyEvolution(m, "PropertyConservedEnergyEvolution");
  bindEvolutionConstructors<PropertyConservedEnergyEvolution>(conservedEnergyEvolution);
  conservedEnergyEvolution
      // the running energies as they are, in internal units like every RunningEnergy (System.compute_total_energies);
      // their __repr__ prints Kelvin (multiplying here as well printed EnergyToKelvin squared)
      .def("result", &PropertyConservedEnergyEvolution::result);

  nanobind::class_<Thermostat>(m, "Thermostat")
      .def(nanobind::init<std::size_t, std::size_t, double>(),
            nanobind::arg("thermostat_chain_length") = 5, 
            nanobind::arg("number_of_yoshida_suzuki_steps") = 5,
            nanobind::arg("time_scale_parameter") = 0.15)
      .def_prop_rw(
          "langevin", [](const Thermostat& t) { return t.isLangevin(); },
          [](Thermostat& t, bool langevin)
          { t.type = langevin ? Thermostat::Type::Langevin : Thermostat::Type::NoseHoover; },
          "Whether the temperature is held by a Langevin thermostat (friction and random kicks) instead of the "
          "Nose-Hoover chain.")
      .def_prop_rw(
          "langevin_friction", [](const Thermostat& t) { return t.langevinFriction; },
          [](Thermostat& t, double friction)
          {
            if (!(friction > 0.0)) throw std::invalid_argument("langevin_friction must be positive");
            t.langevinFriction = friction;
          },
          "The friction coefficient of the Langevin thermostat [1/ps].");

  nanobind::class_<EnthalpyOfAdsorptionData>(m, "EnthalpyOfAdsorptionData")
      .def(nanobind::init<std::size_t>(), nanobind::arg("number_of_components"))
      .def_ro("values", &EnthalpyOfAdsorptionData::values)
      .def("__getitem__",
           [](const EnthalpyOfAdsorptionData& data, std::ptrdiff_t index)
           {
             // a Python sequence ends with an IndexError (iteration stops on it), not past the end of the values, and
             // a negative index counts from the end
             const std::ptrdiff_t size = static_cast<std::ptrdiff_t>(data.values.size());
             if (index < 0) index += size;
             if (index < 0 || index >= size) throw nanobind::index_error();
             return data.values[static_cast<std::size_t>(index)];
           })
      .def("__len__", [](const EnthalpyOfAdsorptionData& data) { return data.values.size(); });

  nanobind::class_<PropertyEnthalpy>(m, "PropertyEnthalpy")
      .def(nanobind::init<std::size_t, std::size_t>(),
            nanobind::arg("number_of_blocks"), 
            nanobind::arg("number_of_components"))
      // 'ideal_gas_energies': the average intramolecular energy of an isolated molecule of each swappable component
      // [K], subtracted as in the output file (needed for flexible molecules); the list is empty where the
      // enthalpy is undefined
      .def("result",
           [](const PropertyEnthalpy& p, const std::vector<double>& idealGasEnergies)
           {
             std::vector<double> internal(idealGasEnergies);
             for (double& energy : internal) energy *= Units::KelvinToEnergy;
             return Units::EnergyToKelvin * p.result(internal);
           },
           nanobind::arg("ideal_gas_energies") = std::vector<double>());

  nanobind::class_<PartialMolarPropertiesData>(m, "PartialMolarPropertiesData")
      .def(nanobind::init<std::size_t>(),
            nanobind::arg("number_of_components"))
      .def_ro("partial_molar_energy", &PartialMolarPropertiesData::partialMolarEnergy)
      .def_ro("partial_molar_volume", &PartialMolarPropertiesData::partialMolarVolume);

  nanobind::class_<PropertyPartialMolarProperties>(m, "PropertyPartialMolarProperties")
      .def(nanobind::init<std::size_t, std::size_t>(),
            nanobind::arg("number_of_blocks"),
            nanobind::arg("number_of_components"))
      // convert the partial molar energies to units of Kelvin
      .def("average_properties",
           [](const PropertyPartialMolarProperties& p)
           {
             auto [average, error] = p.average();
             for (double& value : average.partialMolarEnergy) value *= Units::EnergyToKelvin;
             for (double& value : error.partialMolarEnergy) value *= Units::EnergyToKelvin;
             return std::pair{average, error};
           });
}
