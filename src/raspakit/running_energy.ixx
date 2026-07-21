module;

export module running_energy;

import std;

import archive;
import scaling;
import json;

/**
 * \brief Maximum number of simultaneously tracked thermodynamic-integration lambdas.
 *
 * Each fractional molecule that participates in thermodynamic integration is tagged with a
 * group id (Atom::groupId, 1-based; 0 means untracked). The dU/dlambda contributions are
 * accumulated separately per group so that up to this many independent CFCMC lambdas can be
 * followed at the same time.
 */
export inline constexpr std::size_t maximumNumberOfDUDlambdaGroups = 4;

/**
 * \brief Accumulates energy components during simulation.
 *
 * The RunningEnergy struct stores various components of energy calculated during a simulation run, including potential
 * and kinetic energies, as well as contributions from various interactions such as van der Waals and Coulomb forces. It
 * provides methods to calculate total energies, perform arithmetic operations, and output energy information.
 */
export struct RunningEnergy
{
  /**
   * \brief Default constructor initializing all energy components to zero.
   *
   * Initializes all the energy components and related terms to zero, preparing the RunningEnergy object for
   * accumulation during simulation.
   */
  RunningEnergy()
      : externalFieldVDW(0.0),
        frameworkMoleculeVDW(0.0),
        moleculeMoleculeVDW(0.0),
        externalFieldCharge(0.0),
        frameworkMoleculeCharge(0.0),
        moleculeMoleculeCharge(0.0),
        ewald_fourier(0.0),
        ewald_self(0.0),
        ewald_exclusion(0.0),
        bond(0.0),
        ureyBradley(0.0),
        bend(0.0),
        inversionBend(0.0),
        outOfPlaneBend(0.0),
        torsion(0.0),
        improperTorsion(0.0),
        bondBond(0.0),
        bondBend(0.0),
        bondTorsion(0.0),
        bendBend(0.0),
        bendTorsion(0.0),
        intraVDW(0.0),
        intraCoul(0.0),
        tail(0.0),
        polarization(0.0),
        dudlambdaVDW{},
        dudlambdaCharge{},
        dudlambdaEwald{},
        translationalKineticEnergy(0.0),
        rotationalKineticEnergy(0.0),
        NoseHooverEnergy(0.0),
        thermobarostatEnergy(0.0)
  {
  }

  /**
   * \brief Generates a formatted string of Monte Carlo energy components.
   *
   * Returns a string containing the total potential energy and individual energy components relevant for Monte Carlo
   * simulations, formatted for display.
   *
   * \return A string representing the energy components in Monte Carlo simulations.
   */
  std::string printMC() const;

  /**
   * \brief Generates a formatted string of molecular dynamics energy components.
   *
   * Returns a string containing the total potential and kinetic energies, as well as individual energy components
   * relevant for molecular dynamics simulations, formatted for display.
   *
   * \return A string representing the energy components in molecular dynamics simulations.
   */
  std::string printMD() const;

  /**
   * \brief Generates a formatted string comparing energy components between two Monte Carlo states.
   *
   * Computes the difference between the current RunningEnergy and another instance, and returns a string that shows the
   * comparison of energy components, highlighting any drifts or discrepancies.
   *
   * \param other Another RunningEnergy instance to compare with.
   * \return A string representing the differences in energy components between two Monte Carlo states.
   */
  std::string printMCDiff(RunningEnergy& other) const;

  /**
   * \brief Generates a formatted string comparing energy components between two molecular dynamics states.
   *
   * Computes the difference between the current RunningEnergy and another instance, and returns a string that shows the
   * comparison of energy components, highlighting any drifts or discrepancies in molecular dynamics simulations.
   *
   * \param other Another RunningEnergy instance to compare with.
   * \return A string representing the differences in energy components between two molecular dynamics states.
   */
  std::string printMDDiff(RunningEnergy& other) const;

  /**
   * \brief Generates a labeled formatted string of Monte Carlo energy components.
   *
   * Returns a string containing the total potential energy and individual energy components relevant for Monte Carlo
   * simulations, with an additional label for identification.
   *
   * \param label A string label to include in the output.
   * \return A labeled string representing the energy components in Monte Carlo simulations.
   */
  std::string printMC(const std::string& label) const;

  /**
   * \brief Generates a labeled formatted string of molecular dynamics energy components, including conserved energy and
   * drift.
   *
   * Returns a string containing the total potential and kinetic energies, conserved energy, drift from a reference
   * energy, and individual energy components relevant for molecular dynamics simulations, with an additional label for
   * identification.
   *
   * \param label A string label to include in the output.
   * \param referenceEnergy The reference energy to compute the drift.
   * \return A labeled string representing the energy components in molecular dynamics simulations.
   */
  std::string printMD(const std::string& label, double referenceEnergy) const;

  /**
   * \brief Generates a JSON object of Monte Carlo energy components.
   *
   * Returns a JSON object containing the total potential energy and individual energy components relevant for Monte
   * Carlo simulations.
   *
   * \return A JSON object representing the energy components in Monte Carlo simulations.
   */
  nlohmann::json jsonMC() const;

  /**
   * \brief Generates a JSON object of molecular dynamics energy components.
   *
   * Returns a JSON object containing the total potential and kinetic energies, as well as individual energy components
   * relevant for molecular dynamics simulations.
   *
   * \return A JSON object representing the energy components in molecular dynamics simulations.
   */
  nlohmann::json jsonMD() const;

  /**
   * \brief Computes the total potential energy.
   *
   * Sums up all the potential energy components, including interactions with external fields, framework-molecule
   * interactions, molecule-molecule interactions, Ewald sums, intramolecular energies, tail corrections, and
   * polarization.
   *
   * \return The total potential energy.
   */
  inline double potentialEnergy() const
  {
    return externalFieldVDW + frameworkMoleculeVDW + moleculeMoleculeVDW + externalFieldCharge +
           frameworkMoleculeCharge + moleculeMoleculeCharge + ewald_fourier + ewald_self + ewald_exclusion + bond +
           ureyBradley + bend + inversionBend + outOfPlaneBend + torsion + improperTorsion + bondBond + bondBend +
           bondTorsion + bendBend + bendTorsion + intraVDW + intraCoul + tail + polarization;
  }

  /**
   * \brief Computes the total Coulombic energy.
   *
   * Sums up the Coulombic energy components, including framework-molecule Coulomb interactions, molecule-molecule
   * Coulomb interactions, and Ewald sums.
   *
   * \return The total Coulombic energy.
   */
  inline double CoulombEnergy() const
  {
    return frameworkMoleculeCharge + moleculeMoleculeCharge + ewald_fourier + ewald_self + ewald_exclusion;
  }

  inline double VanDerWaalsEnergy() const
  {
    return frameworkMoleculeVDW + moleculeMoleculeVDW;
  }

  /**
   * \brief Computes the total kinetic energy.
   *
   * Sums up the translational and rotational kinetic energies.
   *
   * \return The total kinetic energy.
   */
  inline double kineticEnergy() const { return translationalKineticEnergy + rotationalKineticEnergy; }
  inline double translationalPartKineticEnergy() const { return translationalKineticEnergy; }
  inline double rotationalPartKineticEnergy() const { return rotationalKineticEnergy; }
  inline double thermostatEnergy() const { return NoseHooverEnergy; }
  inline double extendedThermobarostatEnergy() const { return thermobarostatEnergy; }

  /**
   * \brief Computes the total Ewald Fourier energy.
   *
   * Sums up the Fourier components of the Ewald sum, including the Fourier term, self-interaction correction, and
   * exclusion term.
   *
   * \return The total Ewald Fourier energy.
   */
  inline double ewaldFourier() const { return ewald_fourier + ewald_self + ewald_exclusion; }

  /**
   * \brief Computes the conserved energy in molecular dynamics simulations.
   *
   * Sums up the total potential energy, kinetic energy, and energy associated with the Nose-Hoover thermostat/barostat
   * to compute the conserved energy in molecular dynamics simulations.
   *
   * \return The conserved energy.
   */
  inline double conservedEnergy() const
  {
    return externalFieldVDW + frameworkMoleculeVDW + moleculeMoleculeVDW + externalFieldCharge +
           frameworkMoleculeCharge + moleculeMoleculeCharge + ewald_fourier + ewald_self + ewald_exclusion + bond +
           ureyBradley + bend + inversionBend + outOfPlaneBend + torsion + improperTorsion + bondBond + bondBend +
           bondTorsion + bendBend + bendTorsion + intraVDW + intraCoul + tail + polarization +
           translationalKineticEnergy + rotationalKineticEnergy + NoseHooverEnergy + thermobarostatEnergy;
  }

  /**
   * \brief Computes the derivative of the potential energy with respect to lambda for a TI group.
   *
   * Calculates the derivative of the potential energy with respect to the scaling parameter lambda of the
   * thermodynamic-integration group, considering both van der Waals and Coulombic interactions.
   *
   * \param lambda The scaling parameter of the group.
   * \param group The 1-based thermodynamic-integration group id (0 means untracked and returns zero).
   * \return The derivative of the potential energy with respect to lambda.
   */
  inline double dudlambda(double lambda, std::size_t group = 1) const
  {
    if (group == 0 || group > maximumNumberOfDUDlambdaGroups) return 0.0;
    return Scaling::scalingVDWDerivative(lambda) * dudlambdaVDW[group - 1] +
           Scaling::scalingCoulombDerivative(lambda) * (dudlambdaCharge[group - 1] + dudlambdaEwald[group - 1]);
  }

  /**
   * \brief Accumulates a pair van der Waals dU/dlambda contribution into the per-group accumulators.
   *
   * The potential functions return the symmetric derivative factor X with dU/d(scalingA) = scalingB * X
   * and dU/d(scalingB) = scalingA * X. Each side is routed into the accumulator of the group its atom
   * belongs to (group id 0 means the atom is not tracked).
   */
  inline void addDudlambdaVDW(std::uint8_t groupIdA, std::uint8_t groupIdB, double scalingA, double scalingB,
                              double factor)
  {
    if (groupIdA != 0) dudlambdaVDW[groupIdA - 1] += scalingB * factor;
    if (groupIdB != 0) dudlambdaVDW[groupIdB - 1] += scalingA * factor;
  }

  /// Accumulates a pair real-space Coulomb dU/dlambda contribution (see addDudlambdaVDW).
  inline void addDudlambdaCharge(std::uint8_t groupIdA, std::uint8_t groupIdB, double scalingA, double scalingB,
                                 double factor)
  {
    if (groupIdA != 0) dudlambdaCharge[groupIdA - 1] += scalingB * factor;
    if (groupIdB != 0) dudlambdaCharge[groupIdB - 1] += scalingA * factor;
  }

  /// Accumulates a pair Ewald dU/dlambda contribution (see addDudlambdaVDW).
  inline void addDudlambdaEwald(std::uint8_t groupIdA, std::uint8_t groupIdB, double scalingA, double scalingB,
                                double factor)
  {
    if (groupIdA != 0) dudlambdaEwald[groupIdA - 1] += scalingB * factor;
    if (groupIdB != 0) dudlambdaEwald[groupIdB - 1] += scalingA * factor;
  }

  /// Total van der Waals dU/dlambda summed over all groups (diagnostics only).
  inline double totalDudlambdaVDW() const
  {
    double sum = 0.0;
    for (double value : dudlambdaVDW) sum += value;
    return sum;
  }

  /// Total real-space Coulomb dU/dlambda summed over all groups (diagnostics only).
  inline double totalDudlambdaCharge() const
  {
    double sum = 0.0;
    for (double value : dudlambdaCharge) sum += value;
    return sum;
  }

  /// Total Ewald dU/dlambda summed over all groups (diagnostics only).
  inline double totalDudlambdaEwald() const
  {
    double sum = 0.0;
    for (double value : dudlambdaEwald) sum += value;
    return sum;
  }

  /**
   * \brief Resets all energy components to zero.
   *
   * Sets all the energy components and related terms to zero, effectively clearing the accumulated energies.
   */
  inline void zero()
  {
    externalFieldVDW = 0.0;
    frameworkMoleculeVDW = 0.0;
    moleculeMoleculeVDW = 0.0;
    externalFieldCharge = 0.0;
    frameworkMoleculeCharge = 0.0;
    moleculeMoleculeCharge = 0.0;
    ewald_fourier = 0.0;
    ewald_self = 0.0;
    ewald_exclusion = 0.0;
    bond = 0.0, ureyBradley = 0.0, bend = 0.0, inversionBend = 0.0, outOfPlaneBend = 0.0, torsion = 0.0,
    improperTorsion = 0.0, bondBond = 0.0, bondBend = 0.0, bondTorsion = 0.0, bendBend = 0.0, bendTorsion = 0.0,
    intraVDW = 0.0;
    intraCoul = 0.0;
    tail = 0.0;
    polarization = 0.0;
    dudlambdaVDW.fill(0.0);
    dudlambdaCharge.fill(0.0);
    dudlambdaEwald.fill(0.0);
    translationalKineticEnergy = 0.0;
    rotationalKineticEnergy = 0.0;
    NoseHooverEnergy = 0.0;
    thermobarostatEnergy = 0.0;
  }

  inline RunningEnergy& operator+=(const RunningEnergy& b)
  {
    externalFieldVDW += b.externalFieldVDW;
    frameworkMoleculeVDW += b.frameworkMoleculeVDW;
    moleculeMoleculeVDW += b.moleculeMoleculeVDW;
    externalFieldCharge += b.externalFieldCharge;
    frameworkMoleculeCharge += b.frameworkMoleculeCharge;
    moleculeMoleculeCharge += b.moleculeMoleculeCharge;
    ewald_fourier += b.ewald_fourier;
    ewald_self += b.ewald_self;
    ewald_exclusion += b.ewald_exclusion;
    bond += b.bond;
    ureyBradley += b.ureyBradley;
    bend += b.bend;
    inversionBend += b.inversionBend;
    outOfPlaneBend += b.outOfPlaneBend;
    torsion += b.torsion;
    improperTorsion += b.improperTorsion;
    bondBond += b.bondBond;
    bondBend += b.bondBend;
    bondTorsion += b.bondTorsion;
    bendBend += b.bendBend;
    bendTorsion += b.bendTorsion;
    intraVDW += b.intraVDW;
    intraCoul += b.intraCoul;
    tail += b.tail;
    polarization += b.polarization;
    for (std::size_t i = 0; i < maximumNumberOfDUDlambdaGroups; ++i)
    {
      dudlambdaVDW[i] += b.dudlambdaVDW[i];
      dudlambdaCharge[i] += b.dudlambdaCharge[i];
      dudlambdaEwald[i] += b.dudlambdaEwald[i];
    }
    translationalKineticEnergy += b.translationalKineticEnergy;
    rotationalKineticEnergy += b.rotationalKineticEnergy;
    NoseHooverEnergy += b.NoseHooverEnergy;
    thermobarostatEnergy += b.thermobarostatEnergy;

    return *this;
  }

  inline RunningEnergy& operator-=(const RunningEnergy& b)
  {
    externalFieldVDW -= b.externalFieldVDW;
    frameworkMoleculeVDW -= b.frameworkMoleculeVDW;
    moleculeMoleculeVDW -= b.moleculeMoleculeVDW;
    externalFieldCharge -= b.externalFieldCharge;
    frameworkMoleculeCharge -= b.frameworkMoleculeCharge;
    moleculeMoleculeCharge -= b.moleculeMoleculeCharge;
    ewald_fourier -= b.ewald_fourier;
    ewald_self -= b.ewald_self;
    ewald_exclusion -= b.ewald_exclusion;
    bond -= b.bond;
    ureyBradley -= b.ureyBradley;
    bend -= b.bend;
    inversionBend -= b.inversionBend;
    outOfPlaneBend -= b.outOfPlaneBend;
    torsion -= b.torsion;
    improperTorsion -= b.improperTorsion;
    bondBond -= b.bondBond;
    bondBend -= b.bondBend;
    bondTorsion -= b.bondTorsion;
    bendBend -= b.bendBend;
    bendTorsion -= b.bendTorsion;
    intraVDW -= b.intraVDW;
    intraCoul -= b.intraCoul;
    tail -= b.tail;
    polarization -= b.polarization;
    for (std::size_t i = 0; i < maximumNumberOfDUDlambdaGroups; ++i)
    {
      dudlambdaVDW[i] -= b.dudlambdaVDW[i];
      dudlambdaCharge[i] -= b.dudlambdaCharge[i];
      dudlambdaEwald[i] -= b.dudlambdaEwald[i];
    }
    translationalKineticEnergy -= b.translationalKineticEnergy;
    rotationalKineticEnergy -= b.rotationalKineticEnergy;
    NoseHooverEnergy -= b.NoseHooverEnergy;
    thermobarostatEnergy -= b.thermobarostatEnergy;

    return *this;
  }

  inline RunningEnergy operator-() const
  {
    RunningEnergy v{};
    v.externalFieldVDW = -externalFieldVDW;
    v.frameworkMoleculeVDW = -frameworkMoleculeVDW;
    v.moleculeMoleculeVDW = -moleculeMoleculeVDW;
    v.externalFieldCharge = -externalFieldCharge;
    v.frameworkMoleculeCharge = -frameworkMoleculeCharge;
    v.moleculeMoleculeCharge = -moleculeMoleculeCharge;
    v.ewald_fourier = -ewald_fourier;
    v.ewald_self = -ewald_self;
    v.ewald_exclusion = -ewald_exclusion;
    v.bond = -bond;
    v.ureyBradley = -ureyBradley;
    v.bend = -bend;
    v.inversionBend = -inversionBend;
    v.outOfPlaneBend = -outOfPlaneBend;
    v.torsion = -torsion;
    v.improperTorsion = -improperTorsion;
    v.bondBond = -bondBond;
    v.bondBend = -bondBend;
    v.bondTorsion = -bondTorsion;
    v.bendBend = -bendBend;
    v.bendTorsion = -bendTorsion;
    v.intraVDW = -intraVDW;
    v.intraCoul = -intraCoul;
    v.tail = -tail;
    v.polarization = -polarization;
    for (std::size_t i = 0; i < maximumNumberOfDUDlambdaGroups; ++i)
    {
      v.dudlambdaVDW[i] = -dudlambdaVDW[i];
      v.dudlambdaCharge[i] = -dudlambdaCharge[i];
      v.dudlambdaEwald[i] = -dudlambdaEwald[i];
    }
    v.translationalKineticEnergy = -translationalKineticEnergy;
    v.rotationalKineticEnergy = -rotationalKineticEnergy;
    v.NoseHooverEnergy = -NoseHooverEnergy;
    v.thermobarostatEnergy = -thermobarostatEnergy;

    return v;
  }

  inline RunningEnergy& operator*=(const double& b)
  {
    externalFieldVDW *= b;
    frameworkMoleculeVDW *= b;
    moleculeMoleculeVDW *= b;
    externalFieldCharge *= b;
    frameworkMoleculeCharge *= b;
    moleculeMoleculeCharge *= b;
    ewald_fourier *= b;
    ewald_self *= b;
    ewald_exclusion *= b;
    bond *= b;
    ureyBradley *= b;
    bend *= b;
    inversionBend *= b;
    outOfPlaneBend *= b;
    torsion *= b;
    improperTorsion *= b;
    bondBond *= b;
    bondBend *= b;
    bondTorsion *= b;
    bendBend *= b;
    bendTorsion *= b;
    intraVDW *= b;
    intraCoul *= b;
    tail *= b;
    polarization *= b;
    for (std::size_t i = 0; i < maximumNumberOfDUDlambdaGroups; ++i)
    {
      dudlambdaVDW[i] *= b;
      dudlambdaCharge[i] *= b;
      dudlambdaEwald[i] *= b;
    }
    translationalKineticEnergy *= b;
    rotationalKineticEnergy *= b;
    NoseHooverEnergy *= b;
    thermobarostatEnergy *= b;

    return *this;
  }


  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  double externalFieldVDW;         ///< Energy from van der Waals interactions with external field.
  double frameworkMoleculeVDW;     ///< Energy from van der Waals interactions between framework and molecules.
  double moleculeMoleculeVDW;      ///< Energy from van der Waals interactions between molecules.
  double externalFieldCharge;      ///< Energy from Coulomb interactions with external field.
  double frameworkMoleculeCharge;  ///< Energy from Coulomb interactions between framework and molecules.
  double moleculeMoleculeCharge;   ///< Energy from Coulomb interactions between molecules.
  double ewald_fourier;            ///< Fourier component of Ewald sum for Coulomb interactions.
  double ewald_self;               ///< Self-interaction correction in Ewald sum.
  double ewald_exclusion;          ///< Exclusion term in Ewald sum for Coulomb interactions.
  double bond;
  double ureyBradley;
  double bend;
  double inversionBend;
  double outOfPlaneBend;
  double torsion;
  double improperTorsion;
  double bondBond;
  double bondBend;
  double bondTorsion;
  double bendBend;
  double bendTorsion;
  double intraVDW;                    ///< Intramolecular van der Waals energy.
  double intraCoul;                   ///< Intramolecular Coulomb energy.
  double tail;                        ///< Tail correction energy for van der Waals interactions.
  double polarization;                ///< Energy contribution from polarization effects.
  std::array<double, maximumNumberOfDUDlambdaGroups>
      dudlambdaVDW;  ///< Per-group derivative of van der Waals energy with respect to lambda.
  std::array<double, maximumNumberOfDUDlambdaGroups>
      dudlambdaCharge;  ///< Per-group derivative of Coulomb energy with respect to lambda (real space).
  std::array<double, maximumNumberOfDUDlambdaGroups>
      dudlambdaEwald;  ///< Per-group derivative of Coulomb energy with respect to lambda (Ewald sum).
  double translationalKineticEnergy;  ///< Translational kinetic energy.
  double rotationalKineticEnergy;     ///< Rotational kinetic energy.
  double NoseHooverEnergy;            ///< Energy associated with Nose-Hoover thermostat/barostat.
  double thermobarostatEnergy;        ///< Pressure-volume, cell kinetic, and barostat-chain energy.

  /**
   * \brief Returns a string representation of the RunningEnergy object.
   *
   * \return A string describing the RunningEnergy object.
   */
  std::string repr() const;


  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const RunningEnergy& c);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, RunningEnergy& c);
};

export inline RunningEnergy operator+(const RunningEnergy& a, const RunningEnergy& b)
{
  RunningEnergy m{};
  m.externalFieldVDW = a.externalFieldVDW + b.externalFieldVDW;
  m.frameworkMoleculeVDW = a.frameworkMoleculeVDW + b.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = a.moleculeMoleculeVDW + b.moleculeMoleculeVDW;
  m.externalFieldCharge = a.externalFieldCharge + b.externalFieldCharge;
  m.frameworkMoleculeCharge = a.frameworkMoleculeCharge + b.frameworkMoleculeCharge;
  m.moleculeMoleculeCharge = a.moleculeMoleculeCharge + b.moleculeMoleculeCharge;
  m.ewald_fourier = a.ewald_fourier + b.ewald_fourier;
  m.ewald_self = a.ewald_self + b.ewald_self;
  m.ewald_exclusion = a.ewald_exclusion + b.ewald_exclusion;
  m.bond = a.bond + b.bond;
  m.ureyBradley = a.ureyBradley + b.ureyBradley;
  m.bend = a.bend + b.bend;
  m.inversionBend = a.inversionBend + b.inversionBend;
  m.outOfPlaneBend = a.outOfPlaneBend + b.outOfPlaneBend;
  m.torsion = a.torsion + b.torsion;
  m.improperTorsion = a.improperTorsion + b.improperTorsion;
  m.bondBond = a.bondBond + b.bondBond;
  m.bondBend = a.bondBend + b.bondBend;
  m.bondTorsion = a.bondTorsion + b.bondTorsion;
  m.bendBend = a.bendBend + b.bendBend;
  m.bendTorsion = a.bendTorsion + b.bendTorsion;
  m.intraVDW = a.intraVDW + b.intraVDW;
  m.intraCoul = a.intraCoul + b.intraCoul;
  m.tail = a.tail + b.tail;
  m.polarization = a.polarization + b.polarization;
  for (std::size_t i = 0; i < maximumNumberOfDUDlambdaGroups; ++i)
  {
    m.dudlambdaVDW[i] = a.dudlambdaVDW[i] + b.dudlambdaVDW[i];
    m.dudlambdaCharge[i] = a.dudlambdaCharge[i] + b.dudlambdaCharge[i];
    m.dudlambdaEwald[i] = a.dudlambdaEwald[i] + b.dudlambdaEwald[i];
  }
  m.translationalKineticEnergy = a.translationalKineticEnergy + b.translationalKineticEnergy;
  m.rotationalKineticEnergy = a.rotationalKineticEnergy + b.rotationalKineticEnergy;
  m.NoseHooverEnergy = a.NoseHooverEnergy + b.NoseHooverEnergy;
  m.thermobarostatEnergy = a.thermobarostatEnergy + b.thermobarostatEnergy;

  return m;
}

export inline RunningEnergy operator-(const RunningEnergy& a, const RunningEnergy& b)
{
  RunningEnergy m{};
  m.externalFieldVDW = a.externalFieldVDW - b.externalFieldVDW;
  m.frameworkMoleculeVDW = a.frameworkMoleculeVDW - b.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = a.moleculeMoleculeVDW - b.moleculeMoleculeVDW;
  m.externalFieldCharge = a.externalFieldCharge - b.externalFieldCharge;
  m.frameworkMoleculeCharge = a.frameworkMoleculeCharge - b.frameworkMoleculeCharge;
  m.moleculeMoleculeCharge = a.moleculeMoleculeCharge - b.moleculeMoleculeCharge;
  m.ewald_fourier = a.ewald_fourier - b.ewald_fourier;
  m.ewald_self = a.ewald_self - b.ewald_self;
  m.ewald_exclusion = a.ewald_exclusion - b.ewald_exclusion;
  m.bond = a.bond - b.bond;
  m.ureyBradley = a.ureyBradley - b.ureyBradley;
  m.bend = a.bend - b.bend;
  m.inversionBend = a.inversionBend - b.inversionBend;
  m.outOfPlaneBend = a.outOfPlaneBend - b.outOfPlaneBend;
  m.torsion = a.torsion - b.torsion;
  m.improperTorsion = a.improperTorsion - b.improperTorsion;
  m.bondBond = a.bondBond - b.bondBond;
  m.bondBend = a.bondBend - b.bondBend;
  m.bondTorsion = a.bondTorsion - b.bondTorsion;
  m.bendBend = a.bendBend - b.bendBend;
  m.bendTorsion = a.bendTorsion - b.bendTorsion;
  m.intraVDW = a.intraVDW - b.intraVDW;
  m.intraCoul = a.intraCoul - b.intraCoul;
  m.tail = a.tail - b.tail;
  m.polarization = a.polarization - b.polarization;
  for (std::size_t i = 0; i < maximumNumberOfDUDlambdaGroups; ++i)
  {
    m.dudlambdaVDW[i] = a.dudlambdaVDW[i] - b.dudlambdaVDW[i];
    m.dudlambdaCharge[i] = a.dudlambdaCharge[i] - b.dudlambdaCharge[i];
    m.dudlambdaEwald[i] = a.dudlambdaEwald[i] - b.dudlambdaEwald[i];
  }
  m.translationalKineticEnergy = a.translationalKineticEnergy - b.translationalKineticEnergy;
  m.rotationalKineticEnergy = a.rotationalKineticEnergy - b.rotationalKineticEnergy;
  m.NoseHooverEnergy = a.NoseHooverEnergy - b.NoseHooverEnergy;
  m.thermobarostatEnergy = a.thermobarostatEnergy - b.thermobarostatEnergy;
  return m;
}

export inline RunningEnergy operator*(double a, const RunningEnergy b)
{
  RunningEnergy m{};
  m.externalFieldVDW = a * b.externalFieldVDW;
  m.frameworkMoleculeVDW = a * b.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = a * b.moleculeMoleculeVDW;
  m.externalFieldCharge = a * b.externalFieldCharge;
  m.frameworkMoleculeCharge = a * b.frameworkMoleculeCharge;
  m.moleculeMoleculeCharge = a * b.moleculeMoleculeCharge;
  m.ewald_fourier = a * b.ewald_fourier;
  m.ewald_self = a * b.ewald_self;
  m.ewald_exclusion = a * b.ewald_exclusion;
  m.bond = a * b.bond;
  m.ureyBradley = a * b.ureyBradley;
  m.bend = a * b.bend;
  m.inversionBend = a * b.inversionBend;
  m.outOfPlaneBend = a * b.outOfPlaneBend;
  m.torsion = a * b.torsion;
  m.improperTorsion = a * b.improperTorsion;
  m.bondBond = a * b.bondBond;
  m.bondBend = a * b.bondBend;
  m.bondTorsion = a * b.bondTorsion;
  m.bendBend = a * b.bendBend;
  m.bendTorsion = a * b.bendTorsion;
  m.intraVDW = a * b.intraVDW;
  m.intraCoul = a * b.intraCoul;
  m.tail = a * b.tail;
  m.polarization = a * b.polarization;
  for (std::size_t i = 0; i < maximumNumberOfDUDlambdaGroups; ++i)
  {
    m.dudlambdaVDW[i] = a * b.dudlambdaVDW[i];
    m.dudlambdaCharge[i] = a * b.dudlambdaCharge[i];
    m.dudlambdaEwald[i] = a * b.dudlambdaEwald[i];
  }
  m.translationalKineticEnergy = a * b.translationalKineticEnergy;
  m.rotationalKineticEnergy = a * b.rotationalKineticEnergy;
  m.NoseHooverEnergy = a * b.NoseHooverEnergy;
  m.thermobarostatEnergy = a * b.thermobarostatEnergy;

  return m;
}

export inline RunningEnergy operator*(const RunningEnergy a, double b)
{
  RunningEnergy m{};
  m.externalFieldVDW = b * a.externalFieldVDW;
  m.frameworkMoleculeVDW = b * a.frameworkMoleculeVDW;
  m.moleculeMoleculeVDW = b * a.moleculeMoleculeVDW;
  m.externalFieldCharge = b * a.externalFieldCharge;
  m.frameworkMoleculeCharge = b * a.frameworkMoleculeCharge;
  m.moleculeMoleculeCharge = b * a.moleculeMoleculeCharge;
  m.ewald_fourier = b * a.ewald_fourier;
  m.ewald_self = b * a.ewald_self;
  m.ewald_exclusion = b * a.ewald_exclusion;
  m.bond = b * a.bond;
  m.ureyBradley = b * a.ureyBradley;
  m.bend = b * a.bend;
  m.inversionBend = b * a.inversionBend;
  m.outOfPlaneBend = b * a.outOfPlaneBend;
  m.torsion = b * a.torsion;
  m.improperTorsion = b * a.improperTorsion;
  m.bondBond = b * a.bondBond;
  m.bondBend = b * a.bondBend;
  m.bondTorsion = b * a.bondTorsion;
  m.bendBend = b * a.bendBend;
  m.bendTorsion = b * a.bendTorsion;
  m.intraVDW = b * a.intraVDW;
  m.intraCoul = b * a.intraCoul;
  m.tail = b * a.tail;
  m.polarization = b * a.polarization;
  for (std::size_t i = 0; i < maximumNumberOfDUDlambdaGroups; ++i)
  {
    m.dudlambdaVDW[i] = b * a.dudlambdaVDW[i];
    m.dudlambdaCharge[i] = b * a.dudlambdaCharge[i];
    m.dudlambdaEwald[i] = b * a.dudlambdaEwald[i];
  }
  m.translationalKineticEnergy = b * a.translationalKineticEnergy;
  m.rotationalKineticEnergy = b * a.rotationalKineticEnergy;
  m.NoseHooverEnergy = b * a.NoseHooverEnergy;
  m.thermobarostatEnergy = b * a.thermobarostatEnergy;

  return m;
}
