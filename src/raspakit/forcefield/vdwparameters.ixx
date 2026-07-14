module;

export module vdwparameters;

import std;

import archive;
import double4;
import units;

/**
 * \brief Represents the van der Waals parameters for particle interactions.
 *
 * The VDWParameters struct encapsulates parameters for different types of van der Waals potentials,
 * including Lennard-Jones, Buckingham, Morse, Feynmann-Hibbs, MM3, Born-Huggins-Meyer, and many
 * potentials ported from RASPA2 (shifted-force Lennard-Jones, 12-6, 12-6-2-0, CFF, Matsuoka-Clementi-
 * Yoshimine, generic, Pellenq-Nicholson, hydrated-ion-water, Mie, and Weeks-Chandler-Andersen).
 * It includes the potential parameters, a shift value calculated at the cutoff distance,
 * tail correction energy, and the type of potential.
 *
 * All potentials support the continuous-fractional (CFCMC) lambda-scaling used for the
 * Lennard-Jones potential: the potential is evaluated at the softened distance
 *   r_soft^6 = r^6 + 0.5 (1 - lambda)^2 w^6
 * and multiplied by lambda. Here 'w' is the soft-core reference diameter of the potential
 * (the distance where the fully-coupled potential crosses zero, i.e. sigma for Lennard-Jones).
 * At lambda = 1 the original potential is recovered exactly; for lambda < 1 the energy remains
 * finite (and repulsive) in the limit r -> 0.
 */
export struct VDWParameters
{
  /**
   * \brief Enumeration of van der Waals potential types.
   */
  enum class Type : std::size_t
  {
    None = 0,
    LennardJones = 1,
    BuckingHam = 2,
    Morse = 3,
    FeynmannHibbs = 4,
    MM3 = 5,
    BornHugginsMeyer = 6,
    LennardJonesShiftedForce = 7,
    Potential12_6 = 8,
    Potential12_6_2_0 = 9,
    CFF9_6 = 10,
    CFFEpsilonSigma = 11,
    MatsuokaClementiYoshimine = 12,
    Generic = 13,
    PellenqNicholson = 14,
    HydratedIonWater = 15,
    Mie = 16,
    WeeksChandlerAndersen = 17,
    LennardJonesSecondOrderTaylorShifted = 18,
    RepulsiveHarmonic = 100
  };

  // Hot members first: the Lennard-Jones code path only touches 'parameters', 'shift', and 'type',
  // which together occupy the first 48 bytes of the struct.
  double4 parameters;  ///< The primary potential parameters (p0..p3, RASPA2 ordering, without shift).
  double shift;        ///< The potential energy shift calculated at the cutoff distance.
  Type type{0};        ///< The type of van der Waals potential.

  /// Additional parameters (p4, p5) and derived constants.
  /// Layout: x,y hold overflow input parameters or precomputed constants
  ///   FeynmannHibbs:            x = FeynmannHibbs pre-factor  hbar^2/(24 mu kB T)  [A^2]
  ///   LennardJonesShiftedForce and LennardJonesSecondOrderTaylorShifted:
  ///                                x = (sigma/rc)^6, y = rc
  ///   Generic:                  x = p4, y = p5
  ///   PellenqNicholson:         x = p4 (C10)
  ///   HydratedIonWater:         x = p4 (C12)
  ///   BornHugginsMeyer:         x = p4 (C8)
  /// w always holds the soft-core reference diameter to the sixth power (w^6) used in the
  /// continuous-fractional lambda-scaling.
  double4 parameters2;

  double tailCorrectionEnergy;    ///< The tail correction energy for the potential.
  double tailCorrectionPressure;  ///< The tail correction pressure for the potential.

  /**
   * \brief Default constructor for VDWParameters.
   *
   * Initializes the parameters with zero values.
   */
  VDWParameters()
      : parameters(0.0, 0.0, 0.0, 0.0),
        shift(0.0),
        type(Type::LennardJones),
        parameters2(0.0, 0.0, 0.0, 0.0),
        tailCorrectionEnergy(0.0),
        tailCorrectionPressure(0.0)
  {
  }

  VDWParameters(double4 parameters, double shift, double tailCorrectionEnergy, double tailCorrectionPressure, Type type)
      : parameters(parameters),
        shift(shift),
        type(type),
        parameters2(0.0, 0.0, 0.0, 0.0),
        tailCorrectionEnergy(tailCorrectionEnergy),
        tailCorrectionPressure(tailCorrectionPressure)
  {
  }

  /**
   * \brief Constructs a Lennard-Jones VDWParameter structure.
   *
   * Initializes the VDWParameters with the given epsilon and sigma values for a Lennard-Jones potential.
   *
   * \param epsilon The strength parameter of the potential in units of Kelvin.
   * \param sigma The size parameter of the potential in units of Angstrom.
   */
  VDWParameters(double epsilon, double sigma)
      : parameters(double4(epsilon * Units::KelvinToEnergy, sigma, 0.0, 0.0)),
        shift(0.0),
        type(Type::LennardJones),
        parameters2(0.0, 0.0, 0.0, 0.0),
        tailCorrectionEnergy(0.0),
        tailCorrectionPressure(0.0)
  {
  }

  VDWParameters(double epsilon, double sigma, Type type)
      : parameters(double4(epsilon * Units::KelvinToEnergy, sigma, 0.0, 0.0)),
        shift(0.0),
        type(type),
        parameters2(0.0, 0.0, 0.0, 0.0),
        tailCorrectionEnergy(0.0),
        tailCorrectionPressure(0.0)
  {
  }

  /**
   * \brief Constructs a VDWParameters structure from a list of raw input parameters.
   *
   * The parameters follow the RASPA2 conventions per potential type (without the trailing shift
   * parameter, which is handled by the 'shift' member). Parameters carrying units of energy are
   * converted from Kelvin to internal energy units; lengths, inverse lengths, exponents, and
   * masses are taken as-is.
   *
   * \param type The type of van der Waals potential.
   * \param values The raw parameter values as given in the force field file.
   */
  VDWParameters(Type type, const std::vector<double> &values);

  /**
   * \brief Returns interaction type from given string
   *
   * \param interactionType string to convert to interaction type
   */
  static Type stringToEnum(std::string interactionType);

  /**
   * \brief Returns a human-readable name for the given potential type.
   */
  static std::string nameOfType(Type type);

  /**
   * \brief Metadata for the input parameters of a potential type.
   */
  struct ParameterMetadata
  {
    std::size_t count;             ///< The number of input parameters.
    std::array<bool, 6> isEnergy;  ///< Whether each parameter carries units of energy (Kelvin input).
  };

  /**
   * \brief Returns the parameter metadata (count and unit conversion mask) for a potential type.
   */
  static ParameterMetadata parameterMetadata(Type type);

  /**
   * \brief The unshifted potential energy at full coupling (lambda = 1).
   *
   * Evaluates the plain potential (no shift applied) at squared distance rr. Derived
   * constants (Feynmann-Hibbs pre-factor, shifted-force cutoff constants) must have been
   * computed via computeDerivedParameters() before calling this function.
   *
   * \param rr The squared distance.
   * \return The unshifted potential energy.
   */
  double potentialEnergyAtFullCoupling(double rr) const;

  /**
   * \brief Computes derived constants for the potential.
   *
   * Fills in the derived slots of 'parameters2': the Feynmann-Hibbs temperature pre-factor,
   * the shifted-force cutoff constants, and the soft-core reference diameter w^6 used in the
   * continuous-fractional lambda-scaling. The soft-core diameter is sigma for sigma-based
   * potentials, derived analytically for inverse-power potentials, and located numerically
   * (as the zero-crossing of the fully-coupled potential) for the remaining forms.
   *
   * \param cutOff The cutoff distance.
   * \param temperature The external temperature (used for the Feynmann-Hibbs potential).
   */
  void computeDerivedParameters(double cutOff, double temperature);

  /**
   * \brief Computes the potential energy shift at the cutoff distance.
   *
   * Calculates the shift in potential energy at the specified cutoff distance to ensure continuity
   * of the potential. Requires computeDerivedParameters() to have been called first.
   *
   * \param cutOff The cutoff distance at which to compute the shift.
   */
  void computeShiftAtCutOff(double cutOff)
  {
    shift = 0.0;
    shift = potentialEnergyAtFullCoupling(cutOff * cutOff);
  }

  double strengthParameter() const
  {
    switch (type)
    {
      case Type::LennardJones:
      case Type::FeynmannHibbs:
      case Type::LennardJonesShiftedForce:
      case Type::LennardJonesSecondOrderTaylorShifted:
      case Type::CFFEpsilonSigma:
      case Type::MM3:
      case Type::WeeksChandlerAndersen:
        return parameters.x;
      default:
        return 0.0;
    }
  };

  double sizeParameter() const
  {
    switch (type)
    {
      case Type::LennardJones:
      case Type::FeynmannHibbs:
      case Type::LennardJonesShiftedForce:
      case Type::LennardJonesSecondOrderTaylorShifted:
      case Type::CFFEpsilonSigma:
      case Type::MM3:
      case Type::WeeksChandlerAndersen:
        return parameters.y;
      default:
        return 0.0;
    }
  };

  bool operator==(const VDWParameters &other) const;
  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VDWParameters &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VDWParameters &p);
};

/**
 * \brief Computes the Tang-Toennies damping coefficients for the Pellenq-Nicholson potential.
 *
 * \param r The distance.
 * \param b The exponential decay parameter of the Born repulsion.
 * \param f6 The damping factor of the r^-6 term (output).
 * \param f8 The damping factor of the r^-8 term (output).
 * \param f10 The damping factor of the r^-10 term (output).
 */
export inline void computeTangToenniesDampingCoefficients(double r, double b, double &f6, double &f8, double &f10)
{
  double val = std::exp(-b * r);
  double br = b * r;

  double power = 1.0;
  double factorial = 1.0;
  double sum = val;
  for (std::size_t k = 1; k <= 6; ++k)
  {
    power *= br;
    factorial *= static_cast<double>(k);
    sum += power * val / factorial;
  }
  f6 = 1.0 - sum;

  for (std::size_t k = 7; k <= 8; ++k)
  {
    power *= br;
    factorial *= static_cast<double>(k);
    sum += power * val / factorial;
  }
  f8 = 1.0 - sum;

  for (std::size_t k = 9; k <= 10; ++k)
  {
    power *= br;
    factorial *= static_cast<double>(k);
    sum += power * val / factorial;
  }
  f10 = 1.0 - sum;
}
