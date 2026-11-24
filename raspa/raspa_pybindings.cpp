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
        entries[vector][component] = nanobind::cast<double>(item);
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

// A driver works on its own copy of the systems it is given: the kits are value-based and keep no references to the
// objects of the caller. To let those objects hold the results nevertheless, the systems of the driver are assigned
// back to the Python objects they were copied from whenever a phase returns (copy in, copy out). Inside a call-back
// the phase has not returned yet, and the systems of the driver itself ('systems') are the ones to read.
template <typename Driver>
struct DriverOfGivenSystems : Driver
{
  using Driver::Driver;

  std::vector<nanobind::object> givenSystems;

  void copyOut()
  {
    if (givenSystems.size() != this->systems.size()) return;
    for (std::size_t i = 0; i < givenSystems.size(); ++i)
    {
      // an object that was given more than once has no single system to take
      if (std::ranges::count_if(givenSystems,
                                [&](const nanobind::object& given) { return given.is(givenSystems[i]); }) != 1)
      {
        continue;
      }
      nanobind::cast<System&>(givenSystems[i]) = this->systems[i];
    }
  }
};

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
         std::optional<std::size_t> random_seed, std::size_t number_of_blocks, bool output_to_files)
  {
    std::vector<nanobind::object> given;
    std::vector<System> copies;
    for (nanobind::handle item : systems)
    {
      if (!nanobind::isinstance<System>(item)) throw nanobind::type_error("systems: a sequence of System is expected");
      given.push_back(nanobind::borrow(item));
      copies.push_back(nanobind::cast<const System&>(item));
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
               nanobind::arg("random_seed") = nanobind::none(), nanobind::arg("number_of_blocks") = 5,
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
        pythonName, [function](Bound& d, CallBack call_back_function, std::size_t call_back_every)
        {
          (d.*function)(std::move(call_back_function), call_back_every);
          d.copyOut();
        },
        nanobind::arg("call_back_function") = CallBack([]() {}), nanobind::arg("call_back_every") = 100);
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

        readBinaryRestartFile(static_cast<Driver&>(d), file_name);

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
            new (forceField) ForceField(pseudo_atoms, parameters, mixing_rule, cutoff_framework_vdw,
                                        cutoff_molecule_vdw, cutoff_coulomb, shifted, tail_corrections, use_charge);

            const bool hasParameters14 =
                std::ranges::any_of(parameters_14, [](const auto& given) { return given.has_value(); });
            if (binary_interactions.empty() && !hasParameters14 && binary_interactions_14.empty()) return;

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
            for (const auto& [nameA, nameB, pair] : binary_interactions)
            {
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
            for (const auto& [nameA, nameB, pair] : binary_interactions_14)
            {
              const auto [indexA, indexB] = pairOf(nameA, nameB);
              forceField->ensurePair14Table();
              forceField->data14[indexA * numberOfPseudoAtoms + indexB] = pair;
              forceField->data14[indexB * numberOfPseudoAtoms + indexA] = pair;
            }

            // the shifts and the tail corrections are those of the pairs as they are now
            forceField->preComputeDerivedParameters();
            forceField->preComputePotentialShift();
            forceField->preComputeTailCorrection();
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
      .def_ro("cutoff_framework_vdw", &ForceField::cutOffFrameworkVDW)
      .def_ro("cutoff_molecule_vdw", &ForceField::cutOffMoleculeVDW)
      .def_ro("cutoff_coulomb", &ForceField::cutOffCoulomb);

  nanobind::class_<CIFReader> cifReader(m, "CIFReader");
  cifReader
      .def_static("expand_defined_atoms_to_unit_cell", &CIFReader::expandDefinedAtomsToUnitCell,
           nanobind::arg("simulation_box"), nanobind::arg("space_group_hall_number"),
           nanobind::arg("defined_atoms"),
           nanobind::arg("symmetry_operation_strings") = std::vector<std::string>{});

  nanobind::class_<Framework> framework(m, "Framework");

  nanobind::enum_<CIFReader::UseChargesFrom>(framework, "UseChargesFrom", nanobind::is_arithmetic())
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

  // one keyword for every move of Move.Types. The constructor took twenty positional numbers before, of which the
  // last five were handed to the wrong moves ('widom_probability' set the Gibbs identity change, and so on).
  static const std::string moveProbabilitiesSignatureText = moveProbabilitiesSignature();
  nanobind::class_<MCMoveProbabilities>(m, "MCMoveProbabilities")
      .def(
          "__init__",
          [](MCMoveProbabilities* probabilities, nanobind::kwargs keywords)
          {
            new (probabilities) MCMoveProbabilities();
            for (auto [key, value] : keywords)
            {
              const Move::Types move = moveTypeOfProbabilityKeyword(nanobind::cast<std::string>(key));
              if (!nanobind::isinstance<nanobind::float_>(value) && !nanobind::isinstance<nanobind::int_>(value))
              {
                throw nanobind::type_error("MCMoveProbabilities: a probability is a number");
              }
              setMoveProbability(*probabilities, move, nanobind::cast<double>(value));
            }
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
      {{"NONE", BondType::None}, {"FIXED", BondType::Fixed}, {"HARMONIC", BondType::Harmonic}, {"CORE_SHELL_SPRING", BondType::CoreShellSpring}, {"MORSE", BondType::Morse}, {"LJ_12_6", BondType::LJ_12_6}, {"LENNARD_JONES", BondType::LennardJones}, {"BUCKINGHAM", BondType::Buckingham}, {"RESTRAINED_HARMONIC", BondType::RestrainedHarmonic}, {"QUARTIC", BondType::Quartic}, {"CFF_QUARTIC", BondType::CFF_Quartic}, {"MM3", BondType::MM3}},
      BondPotential::numberOfBondParameters);
  bindIntraMolecularPotential<UreyBradleyPotential, UreyBradleyType>(m, "UreyBradleyPotential", "UreyBradleyType",
      {{"FIXED", UreyBradleyType::Fixed}, {"HARMONIC", UreyBradleyType::Harmonic}, {"CORE_SHELL_SPRING", UreyBradleyType::CoreShellSpring}, {"MORSE", UreyBradleyType::Morse}, {"LJ_12_6", UreyBradleyType::LJ_12_6}, {"LENNARD_JONES", UreyBradleyType::LennardJones}, {"BUCKINGHAM", UreyBradleyType::Buckingham}, {"RESTRAINED_HARMONIC", UreyBradleyType::RestrainedHarmonic}, {"QUARTIC", UreyBradleyType::Quartic}, {"CFF_QUARTIC", UreyBradleyType::CFF_Quartic}, {"MM3", UreyBradleyType::MM3}},
      UreyBradleyPotential::numberOfUreyBradleyParameters);
  bindIntraMolecularPotential<BendPotential, BendType>(m, "BendPotential", "BendType",
      {{"RIGID", BendType::Rigid}, {"FIXED", BendType::Fixed}, {"HARMONIC", BendType::Harmonic}, {"CORE_SHELL", BendType::CoreShell}, {"QUARTIC", BendType::Quartic}, {"CFF_QUARTIC", BendType::CFF_Quartic}, {"HARMONIC_COSINE", BendType::HarmonicCosine}, {"COSINE", BendType::Cosine}, {"TAFIPOLSKY", BendType::Tafipolsky}, {"MM3", BendType::MM3}},
      BendPotential::numberOfBendParameters);
  bindIntraMolecularPotential<InversionBendPotential, InversionBendType>(m, "InversionBendPotential", "InversionBendType",
      {{"HARMONIC", InversionBendType::Harmonic}, {"HARMONIC_COSINE", InversionBendType::HarmonicCosine}, {"PLANAR", InversionBendType::Planar}, {"HARMONIC2", InversionBendType::Harmonic2}, {"HARMONIC_COSINE2", InversionBendType::HarmonicCosine2}, {"PLANAR2", InversionBendType::Planar2}, {"MM3", InversionBendType::MM3}},
      InversionBendPotential::numberOfInversionBendParameters);
  bindIntraMolecularPotential<TorsionPotential, TorsionType>(m, "TorsionPotential", "TorsionType",
      {{"FIXED", TorsionType::Fixed}, {"HARMONIC", TorsionType::Harmonic}, {"HARMONIC_COSINE", TorsionType::HarmonicCosine}, {"THREE_COSINE", TorsionType::ThreeCosine}, {"RYCKAERT_BELLEMANS", TorsionType::RyckaertBellemans}, {"TRAPPE", TorsionType::TraPPE}, {"TRAPPE_EXTENDED", TorsionType::TraPPE_Extended}, {"MOD_TRAPPE", TorsionType::ModifiedTraPPE}, {"CVFF", TorsionType::CVFF}, {"CFF", TorsionType::CFF}, {"CFF2", TorsionType::CFF2}, {"OPLS", TorsionType::OPLS}, {"MM3", TorsionType::MM3}, {"FOURIER_SERIES", TorsionType::FourierSeries}, {"FOURIER_SERIES2", TorsionType::FourierSeries2}, {"CVFF_BLOCKED", TorsionType::CVFFBlocked}, {"POLYNOMIAL", TorsionType::Polynomial}},
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

    nanobind::class_<Component> component(m, "Component");
    nanobind::enum_<Component::Type>(component, "Type", nanobind::is_arithmetic())
        .value("ADSORBATE", Component::Type::Adsorbate)
        .value("CATION", Component::Type::Cation);

    component
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
             nanobind::arg("tag") = std::nullopt);

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
      .def(nanobind::init<std::size_t, std::pair<std::size_t, std::size_t>, std::size_t, std::optional<std::size_t>>(),
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


  nanobind::class_<PropertyNumberOfMoleculesEvolution>(m, "PropertyNumberOfMoleculesEvolution")
      .def(nanobind::init<std::size_t, std::size_t, std::size_t, std::optional<std::size_t>>(),
            nanobind::arg("number_of_production_cycles"), 
            nanobind::arg("number_of_components"),
            nanobind::arg("sample_every"), 
            nanobind::arg("write_every") = nanobind::none())
      // the keyword 'number_of_cycles' of the 3.1.0 release
      .def(nanobind::init<std::size_t, std::size_t, std::size_t, std::optional<std::size_t>>(),
            nanobind::arg("number_of_cycles"), 
            nanobind::arg("number_of_components"),
            nanobind::arg("sample_every"), 
            nanobind::arg("write_every") = nanobind::none())
      .def_ro("result", &PropertyNumberOfMoleculesEvolution::result);

  nanobind::class_<PropertyVolumeEvolution>(m, "PropertyVolumeEvolution")
    .def(nanobind::init<std::size_t, std::size_t, std::optional<std::size_t>>(),
            nanobind::arg("number_of_production_cycles"), 
            nanobind::arg("sample_every"), 
            nanobind::arg("write_every") = nanobind::none())
    // the keyword 'number_of_cycles' of the 3.1.0 release
    .def(nanobind::init<std::size_t, std::size_t, std::optional<std::size_t>>(),
            nanobind::arg("number_of_cycles"), 
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

  nanobind::class_<PropertyTemperature>(m, "PropertyTemperature")
      .def("result", [](const PropertyTemperature& p) { return p.average(); },
           "The average temperature and its error [K].");

  nanobind::class_<PropertySimulationBox>(m, "PropertySimulationBox")
      .def("result", [](const PropertySimulationBox& p) { return p.average(); },
           "The average simulation box and its error.");

  nanobind::class_<System>(m, "System")
      .def(nanobind::init<ForceField, std::optional<SimulationBox>, bool, double, std::optional<double>, double,
                          std::optional<Framework>, std::vector<Component>, std::vector<std::vector<double3>>,
                          std::vector<std::size_t>, std::size_t, MCMoveProbabilities>(),
           nanobind::arg("force_field"), nanobind::arg("simulation_box") = nanobind::none(),
           nanobind::arg("has_external_field") = false, nanobind::arg("external_temperature") = 298.0,
           nanobind::arg("external_pressure") = nanobind::none(), nanobind::arg("helium_void_fraction") = 0.0,
           nanobind::arg("framework_components") = nanobind::none(),
           nanobind::arg("components") = std::vector<Component>(),
           nanobind::arg("initial_positions") = std::vector<std::vector<double3>>(),
           nanobind::arg("initial_number_of_molecules") = std::vector<std::size_t>(),
           nanobind::arg("number_of_blocks") = 5, nanobind::arg("system_probabilities") = MCMoveProbabilities())
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
      .def("set_thermostat", &System::setThermostat, nanobind::arg("thermostat") = nanobind::none())
      .def_ro("property_density_grid", &System::propertyDensityGrid)
      .def("set_property_density_grid", &System::setPropertyDensityGrid, nanobind::arg("property") = nanobind::none())
      .def_ro("average_energy_histogram", &System::averageEnergyHistogram)
      .def("set_average_energy_histogram", &System::setAverageEnergyHistogram,
           nanobind::arg("property") = nanobind::none())
      .def_ro("property_number_of_molecules_histogram", &System::averageNumberOfMoleculesHistogram)
      .def("set_number_of_molecules_histogram", &System::setNumberOfMoleculesHistogram,
           nanobind::arg("property") = nanobind::none())
      .def_ro("sample_pdb_movie", &System::samplePDBMovie)
      .def("set_sample_pdb_movie", &System::setSamplePDBMovie, nanobind::arg("property") = nanobind::none())
      .def_ro("property_conventional_rdf", &System::propertyConventionalRadialDistributionFunction)
      .def("set_property_conventional_rdf", &System::setPropertyConventionalRDF,
           nanobind::arg("property") = nanobind::none())
      .def_ro("property_msd", &System::propertyMSD)
      .def("set_property_msd", &System::setPropertyMSD, nanobind::arg("property") = nanobind::none())
      .def_ro("property_vacf", &System::propertyVACF)
      .def("set_property_vacf", &System::setPropertyVACF, nanobind::arg("property") = nanobind::none())
      .def_ro("property_rdf", &System::propertyRadialDistributionFunction)
      .def("set_property_rdf", &System::setPropertyRDF, nanobind::arg("property") = nanobind::none())
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
      // read only: the atoms are kept in step with their dynamics, molecules and groups
      .def_ro("atom_data", &System::atomData)
      .def_rw("translational_degrees_of_freedom", &System::translationalDegreesOfFreedom)
      .def_rw("rotational_degrees_of_freedom", &System::rotationalDegreesOfFreedom)
      .def_rw("time_step", &System::timeStep)
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

  nanobind::class_<PropertyConservedEnergyEvolution>(m, "PropertyConservedEnergyEvolution")
      .def(nanobind::init<std::size_t, std::size_t, std::optional<std::size_t>>(),
           nanobind::arg("number_of_production_cycles"), nanobind::arg("sample_every"),
           nanobind::arg("write_every") = nanobind::none())
      // the keyword 'number_of_cycles' of the 3.1.0 release
      .def(nanobind::init<std::size_t, std::size_t, std::optional<std::size_t>>(),
           nanobind::arg("number_of_cycles"), nanobind::arg("sample_every"),
           nanobind::arg("write_every") = nanobind::none())
      // the running energies as they are, in internal units like every RunningEnergy (System.compute_total_energies);
      // their __repr__ prints Kelvin (multiplying here as well printed EnergyToKelvin squared)
      .def("result", &PropertyConservedEnergyEvolution::result);

  nanobind::class_<Thermostat>(m, "Thermostat")
      .def(nanobind::init<std::size_t, std::size_t, double>(),
            nanobind::arg("thermostat_chain_length") = 5, 
            nanobind::arg("number_of_yoshida_suzuki_steps") = 5,
            nanobind::arg("time_scale_parameter") = 0.15);

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
