module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>
#endif

export module reaction;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import property_lambda_probability_histogram;
import json;

/**
 * \brief Represents a chemical reaction within the simulation.
 *
 * The Reaction struct encapsulates the properties and behaviors of a chemical reaction,
 * including its unique identifier, stoichiometry of reactants and products, a histogram
 * for lambda probabilities, and methods for printing and serialization.
 */
export struct Reaction
{
  /**
   * \brief Default constructor for the Reaction struct.
   *
   * Initializes a Reaction object with default values.
   */
  Reaction() {};

  /**
   * \brief Constructs a Reaction with specified parameters.
   *
   * Initializes a Reaction with the provided identifier and stoichiometry for reactants and products.
   *
   * \param id The unique identifier for the reaction.
   * \param reactantStoichiometry A vector containing the stoichiometry of reactants.
   * \param productStoichiometry A vector containing the stoichiometry of products.
   */
  Reaction(std::size_t id, std::vector<std::size_t> reactantStoichiometry,
           std::vector<std::size_t> productStoichiometry)
      : id(id), reactantStoichiometry(reactantStoichiometry), productStoichiometry(productStoichiometry), lambda(5, 21)
  {
  }

  bool operator==(Reaction const &) const = default;

  std::uint64_t versionNumber{1};  ///< Version number of the Reaction struct.

  std::size_t id;                                  ///< Unique identifier for the reaction.
  std::vector<std::size_t> reactantStoichiometry;  ///< Stoichiometry of reactants.
  std::vector<std::size_t> productStoichiometry;   ///< Stoichiometry of products.

  PropertyLambdaProbabilityHistogram lambda;  ///< Histogram for lambda probabilities.

  /**
   * \brief Returns a string representation of the Reaction status.
   *
   * Generates a string containing information about the Reaction's current state.
   *
   * \return A string representing the Reaction's status.
   */
  std::string printStatus() const;

  /**
   * \brief Returns a JSON representation of the Reaction status.
   *
   * Generates a JSON object containing information about the Reaction's current state.
   *
   * \return A JSON object representing the Reaction's status.
   */
  nlohmann::json jsonStatus() const;

  // Friend functions for serialization
  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Reaction &r);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Reaction &r);
};
