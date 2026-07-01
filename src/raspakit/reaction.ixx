module;

export module reaction;

import std;

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
      : id(id),
        reactantStoichiometry(std::move(reactantStoichiometry)),
        productStoichiometry(std::move(productStoichiometry)),
        lambda(5, 21),
        lambdaProductSide(5, 21)
  {
  }

  bool operator==(Reaction const &) const = default;

  std::uint64_t versionNumber{3};  ///< Version number of the Reaction struct.

  std::size_t id;                                  ///< Unique identifier for the reaction.
  /// Stoichiometry of reactants (one entry per adsorbate component; values may exceed 1).
  std::vector<std::size_t> reactantStoichiometry;
  /// Stoichiometry of products (one entry per adsorbate component; values may exceed 1).
  std::vector<std::size_t> productStoichiometry;

  PropertyLambdaProbabilityHistogram lambda;             ///< Bias histogram (reactant-side fractionals, serial Rx/CFC).
  PropertyLambdaProbabilityHistogram lambdaProductSide;  ///< Bias histogram (product-side fractionals, serial Rx/CFC).

  double currentLambda{0.0};           ///< Current coupling parameter λ ∈ [0, 1].
  double maximumLambdaChange{0.3};     ///< Maximum random change in λ (reactant-side, serial or parallel).
  double maximumLambdaChangeProducts{0.3};  ///< Maximum λ change when product-side fractionals are present (serial).
  double lambdaSwitchPoint{0.5};       ///< λ threshold for fractional vs whole-molecule reaction moves (serial).

  bool serialRxCFC{false};              ///< Use serial Rx/CFC (single fractional side) instead of parallel Rx/CFC.
  bool fractionalSideIsReactants{true};  ///< δ: true = reactant fractionals present, false = product fractionals.

  /// Molecule indices (per component) for reactant fractional molecules used in CFC-RXMC.
  std::vector<std::vector<std::size_t>> reactantFractionalMoleculeIds;
  /// Molecule indices (per component) for product fractional molecules used in CFC-RXMC.
  std::vector<std::vector<std::size_t>> productFractionalMoleculeIds;

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
