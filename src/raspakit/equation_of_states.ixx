module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <fstream>
#include <vector>
#endif

export module equation_of_states;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import component;
import simulationbox;

/**
 * \brief Represents the equation of state for fluid mixtures in the simulation.
 *
 * The EquationOfState struct encapsulates the properties and behaviors associated with
 * computing the thermodynamic properties of fluid mixtures based on different equations of state.
 * It supports various types of equations of state, mixing rules for multi-component systems, and
 * fluid states. The struct provides methods to compute component fluid properties and supports
 * serialization through friend operators.
 */
export struct EquationOfState
{
  /**
   * \brief Enumeration of supported equations of state.
   */
  enum class Type : int
  {
    PengRobinson = 0,       ///< Peng-Robinson equation of state.
    PengRobinsonGasem = 1,  ///< Peng-Robinson-Gasem equation of state.
    SoaveRedlichKwong = 2   ///< Soave-Redlich-Kwong equation of state.
  };

  /**
   * \brief Enumeration of mixing rules for multi-component systems.
   */
  enum class MultiComponentMixingRules : int
  {
    VanDerWaals = 0  ///< van der Waals mixing rules.
  };

  /**
   * \brief Enumeration of possible fluid states.
   */
  enum class FluidState : int
  {
    Unknown = 0,             ///< Fluid state is unknown.
    SuperCriticalFluid = 1,  ///< Supercritical fluid state.
    Vapor = 2,               ///< Vapor state.
    Liquid = 3,              ///< Liquid state.
    VaporLiquid = 4          ///< Vapor-liquid coexistence.
  };

  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  EquationOfState::FluidState fluidState{EquationOfState::FluidState::Unknown};  ///< Current fluid state.
  EquationOfState::Type equationOfState{EquationOfState::Type::PengRobinson};    ///< Type of equation of state used.
  EquationOfState::MultiComponentMixingRules multiComponentMixingRules{
      EquationOfState::MultiComponentMixingRules::VanDerWaals};  ///< Mixing rules for multi-component systems.

  /**
   * \brief Default constructor for the EquationOfState struct.
   *
   * Initializes an EquationOfState object with default values.
   */
  EquationOfState() = default;

  /**
   * \brief Constructs an EquationOfState with specified parameters.
   *
   * Initializes an EquationOfState object and computes the fluid properties
   * for the given components based on the specified equation of state and mixing rules.
   *
   * \param type The type of equation of state to use.
   * \param multiComponentMixingRules The mixing rules for multi-component systems.
   * \param temperature The temperature of the system in Kelvin.
   * \param pressure The pressure of the system in Pascal.
   * \param simulationBox The simulation box containing the system dimensions.
   * \param heliumVoidFraction The fraction of void space occupied by helium.
   * \param components The components present in the fluid mixture.
   */
  EquationOfState(EquationOfState::Type type, EquationOfState::MultiComponentMixingRules multiComponentMixingRules,
                  double temperature, double pressure, const SimulationBox &simulationBox, double heliumVoidFraction,
                  std::vector<Component> &components);

  /**
   * \brief Computes the fluid properties for each component in the mixture.
   *
   * Calculates various thermodynamic properties such as fugacity coefficients,
   * amount of excess molecules, bulk fluid density, and compressibility factors
   * for the components based on the specified equation of state and mixing rules.
   *
   * \param equationOfState The type of equation of state to use.
   * \param multiComponentMixingRules The mixing rules for multi-component systems.
   * \param temperature The temperature of the system in Kelvin.
   * \param pressure The pressure of the system in Pascal.
   * \param simulationBox The simulation box containing the system dimensions.
   * \param heliumVoidFraction The fraction of void space occupied by helium.
   * \param components The components present in the fluid mixture.
   */
  void computeComponentFluidProperties(EquationOfState::Type equationOfState,
                                       EquationOfState::MultiComponentMixingRules multiComponentMixingRules,
                                       double temperature, double pressure, const SimulationBox &simulationBox,
                                       double heliumVoidFraction, std::vector<Component> &components);

  /**
   * \brief Serializes the EquationOfState object to an output archive.
   *
   * \param archive The output archive to serialize to.
   * \param s The EquationOfState object to serialize.
   * \return Reference to the output archive.
   */
  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EquationOfState &s);

  /**
   * \brief Deserializes the EquationOfState object from an input archive.
   *
   * \param archive The input archive to deserialize from.
   * \param s The EquationOfState object to deserialize into.
   * \return Reference to the input archive.
   */
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EquationOfState &s);
};
