module;

export module van_der_waals_potential;

import std;

import archive;
import double3;
import vdwparameters;

/**
 * \brief An explicit intramolecular van der Waals pair term.
 *
 * Pairs the two atoms with a copy of the force field's pair potential between their pseudo-atom types
 * ('VDWParameters': the functional form, its parameters, the shift at the cut-off and the soft-core
 * constants) times the pair 'scaling' (one for an ordinary pair, the 1-4 scaling for a 1-4 pair).
 *
 * The energy, gradient and Hessian routines evaluate these pairs through the force field itself
 * (Potentials::intraMolecularVDW with the atom types), so for them only 'identifiers', 'scaling' and
 * 'pair14' matter. The self-contained 'calculateEnergy' is the force-field-free evaluation used where no force
 * field is at hand (the CBMC lookahead guide): the untruncated, unshifted potential at full coupling,
 * dispatched on 'parameters.type' like every other evaluation.
 */
export struct VanDerWaalsPotential
{
  std::uint64_t versionNumber{3};  ///< Version number for serialization.

  std::array<std::size_t, 2> identifiers{0, 0};  ///< Identifiers of the two atoms forming the pair.
  double scaling{1.0};                           ///< The pair scaling (1-4 scaling, 1 for ordinary pairs).
  VDWParameters parameters{};                    ///< The force-field pair potential between the two atom types.
  bool pair14{false};                            ///< Whether the pair uses the 1-4 pair parameters (ForceField::pair14).

  VanDerWaalsPotential() = default;

  /**
   * \brief Constructs a pair term from the force field's pair potential.
   *
   * \param identifiers The two atoms.
   * \param parameters The force field's pair potential between their types (ForceField::operator()).
   * \param scaling The pair scaling.
   */
  VanDerWaalsPotential(std::array<std::size_t, 2> identifiers, const VDWParameters &parameters, double scaling)
      : identifiers(identifiers), scaling(scaling), parameters(parameters)
  {
  }

  /**
   * \brief Constructs a pair term from raw parameters in force-field file units.
   *
   * \param identifiers The two atoms.
   * \param type The functional form.
   * \param vector_parameters The parameters in the force-field file conventions of 'type' (energies in Kelvin).
   * \param scaling The pair scaling.
   */
  VanDerWaalsPotential(std::array<std::size_t, 2> identifiers, VDWParameters::Type type,
                       const std::vector<double> &vector_parameters, double scaling)
      : identifiers(identifiers), scaling(scaling), parameters(type, vector_parameters)
  {
  }

  bool operator==(VanDerWaalsPotential const &) const = default;

  /// Multiplies the energy parameters of the pair potential by 'factor' (solute tempering).
  void scaleEnergy(double factor) { parameters.scaleEnergy(factor); }

  /**
   * \brief Generates a string representation of the pair term.
   *
   * \return "A - B : <form> p0: ..., scaling: ...".
   */
  std::string print() const;

  /**
   * \brief The scaled, untruncated and unshifted pair energy at full coupling.
   *
   * Dispatches on the functional form of 'parameters' (VDWParameters::potentialEnergyAtFullCoupling).
   */
  double calculateEnergy(const double3 &posA, const double3 &posB) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VanDerWaalsPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VanDerWaalsPotential &b);
};
