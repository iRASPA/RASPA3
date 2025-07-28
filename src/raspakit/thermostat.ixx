module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <print>
#include <tuple>
#include <vector>
#endif

export module thermostat;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import units;
import randomnumbers;
import molecule;

/**
 * \brief Represents a Nose-Hoover thermostat for molecular dynamics simulations.
 *
 * The Thermostat struct encapsulates the properties and methods required to implement a Nose-Hoover
 * thermostat for temperature control in molecular dynamics simulations. It includes parameters for
 * temperature, thermostat chain length, Yoshida-Suzuki steps, and degrees of freedom. The struct
 * provides methods to initialize the thermostat, perform Nose-Hoover NVT integration steps, and
 * compute the thermostat's energy contribution.
 */
export struct Thermostat
{
  /**
   * \brief Default constructor for the Thermostat struct.
   *
   * Initializes a Thermostat object with default values.
   */
  Thermostat() {}

  /**
   * \brief Constructs a Thermostat with specified parameters.
   *
   * Initializes a Thermostat with the provided temperature, thermostat chain length, number of
   * Yoshida-Suzuki steps, time step deltaT, and degrees of freedom for translation and rotation.
   *
   * \param temperature The target temperature for the thermostat.
   * \param thermostatChainLength The length of the thermostat chain.
   * \param numberOfYoshidaSuzukiSteps The number of Yoshida-Suzuki steps for integration.
   * \param deltaT The simulation time step.
   * \param translationalDegreesOfFreedom The number of translational degrees of freedom.
   * \param rotationalDgreesOfFreedom The number of rotational degrees of freedom.
   */
  Thermostat(double temperature, std::size_t thermostatChainLength, std::size_t numberOfYoshidaSuzukiSteps,
             double deltaT, std::size_t translationalDegreesOfFreedom, std::size_t rotationalDgreesOfFreedom);

  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  double temperature;                         ///< The target temperature for the thermostat.
  std::size_t thermostatChainLength;          ///< The length of the thermostat chain.
  double timeScaleParameterThermostat{0.15};  ///< Time scale parameter for the thermostat.
  std::size_t numberOfRespaSteps{5};          ///< Number of RESPA steps.
  std::size_t numberOfYoshidaSuzukiSteps{5};  ///< Number of Yoshida-Suzuki steps for integration.
  double deltaT{};                            ///< The simulation time step.

  std::size_t
      translationalCenterOfMassConstraint{};  ///< Constraint on translational center of mass degrees of freedom.
  std::size_t translationalDegreesOfFreedom;  ///< Number of translational degrees of freedom.
  std::size_t rotationalDgreesOfFreedom;      ///< Number of rotational degrees of freedom.

  std::vector<double>
      thermostatDegreesOfFreedomTranslation;          ///< Degrees of freedom for the translational thermostat chain.
  std::vector<double> thermostatForceTranslation;     ///< Forces for the translational thermostat chain.
  std::vector<double> thermostatVelocityTranslation;  ///< Velocities for the translational thermostat chain.
  std::vector<double> thermostatPositionTranslation;  ///< Positions for the translational thermostat chain.
  std::vector<double> thermostatMassTranslation;      ///< Masses for the translational thermostat chain.

  std::vector<double> thermostatDegreesOfFreedomRotation;  ///< Degrees of freedom for the rotational thermostat chain.
  std::vector<double> thermostatForceRotation;             ///< Forces for the rotational thermostat chain.
  std::vector<double> thermostatVelocityRotation;          ///< Velocities for the rotational thermostat chain.
  std::vector<double> thermostatPositionRotation;          ///< Positions for the rotational thermostat chain.
  std::vector<double> thermostatMassRotation;              ///< Masses for the rotational thermostat chain.

  std::vector<double> w;  ///< Weights for Yoshida-Suzuki integration steps.

  /**
   * \brief Initializes the thermostat with random velocities.
   *
   * Sets up the thermostat chain by initializing velocities using a random number generator.
   *
   * \param random A reference to a RandomNumber generator.
   */
  void initialize(RandomNumber &random);

  /**
   * \brief Performs a Nose-Hoover NVT integration step.
   *
   * Updates the thermostat variables and computes scaling factors for translational and rotational
   * degrees of freedom based on the provided kinetic energies.
   *
   * \param UKineticTranslation The kinetic energy of translation.
   * \param UKineticRotation The kinetic energy of rotation.
   * \return A pair of scaling factors for translation and rotation.
   */
  std::pair<double, double> NoseHooverNVT(double UKineticTranslation, double UKineticRotation);

  /**
   * \brief Computes the total energy of the thermostat.
   *
   * Calculates the energy contribution from both translational and rotational thermostat chains.
   *
   * \return The total energy of the thermostat.
   */
  double getEnergy();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Thermostat &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Thermostat &s);
};
