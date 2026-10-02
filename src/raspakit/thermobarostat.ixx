module;

export module thermobarostat;

import std;

import archive;
import double3;
import double3x3;
import randomnumbers;
import units;
import minimization_cell_layout;

export enum class MolecularDynamicsEnsemble : std::uint8_t { NVE, NVT, NPT, NPTPR, MuVT, MuPT, MuPTPR };

export std::optional<MolecularDynamicsEnsemble> molecularDynamicsEnsembleFromString(std::string_view value);
export std::string molecularDynamicsEnsembleName(MolecularDynamicsEnsemble ensemble);
export bool molecularDynamicsUsesThermostat(MolecularDynamicsEnsemble ensemble);
export bool molecularDynamicsUsesIsotropicBarostat(MolecularDynamicsEnsemble ensemble);
export bool molecularDynamicsUsesFlexibleBarostat(MolecularDynamicsEnsemble ensemble);
export bool molecularDynamicsHasParticleExchange(MolecularDynamicsEnsemble ensemble);

// How the barostat couples the cell to the particles. 'Molecular' drives every molecule through its centre of
// mass (the molecular virial and the centre-of-mass kinetic energy enter the cell equation of motion), so the
// pressure the barostat drives to the set point is the molecular pressure that the code reports. 'Atomic' couples
// every flexible atom individually (atomic virial and atomic kinetic energy); the corresponding estimator is then
// the atomic pressure. Rigid molecules and rigid groups are always coupled through their centre of mass.
export enum class BarostatCoupling : std::uint8_t { Molecular, Atomic };

export std::optional<BarostatCoupling> barostatCouplingFromString(std::string_view value);
export std::string barostatCouplingName(BarostatCoupling coupling);

export struct Thermobarostat
{
  std::uint64_t versionNumber{2};
  MolecularDynamicsEnsemble ensemble{MolecularDynamicsEnsemble::NVE};
  CellMinimizationType cellType{CellMinimizationType::Isotropic};
  MonoclinicAngleType monoclinicAngle{MonoclinicAngleType::Beta};
  BarostatCoupling coupling{BarostatCoupling::Molecular};
  double temperature{300.0};
  double pressure{};
  double timeStep{0.0005};
  double timeScaleParameterBarostat{1.0};
  std::size_t translationalDegreesOfFreedom{};
  std::size_t cellDegreesOfFreedom{1};
  std::size_t chainLength{5};
  std::size_t numberOfRespaSteps{5};
  std::size_t numberOfYoshidaSuzukiSteps{5};

  // NPT uses x=ln(V) with velocity xdot=d ln(V)/dt; NPTPR uses the logarithmic cell-rate matrix.
  // 'logVolumeMass' is the Martyna-Tobias-Klein mass W=(N_f+3) k_B T tau_b^2 of the strain
  // epsilon=ln(V)/3, so the variable x=3 epsilon has mass W/9: its kinetic energy is W xdot^2/18 and its
  // equation of motion is xddot = 3 G_epsilon / W (see 'logVolumeKineticEnergy' and the NPT drivers).
  double logVolumePosition{};
  double logVolumeVelocity{};
  double logVolumeMass{1.0};
  double3x3 cellVelocity{};
  double cellMass{1.0};

  std::vector<double> chainDegreesOfFreedom{};
  std::vector<double> chainForce{};
  std::vector<double> chainVelocity{};
  std::vector<double> chainPosition{};
  std::vector<double> chainMass{};
  std::vector<double> yoshidaSuzukiWeights{};

  Thermobarostat() = default;
  Thermobarostat(MolecularDynamicsEnsemble ensemble, CellMinimizationType cellType,
                 MonoclinicAngleType monoclinicAngle, double temperature, double pressure, double timeStep,
                 std::size_t translationalDegreesOfFreedom, std::size_t chainLength = 5,
                 std::size_t numberOfYoshidaSuzukiSteps = 5, double timeScaleParameterBarostat = 1.0);

  void initialize(RandomNumber& random);
  void refreshDegreesOfFreedom(RandomNumber& random, std::size_t translationalDegreesOfFreedom, double volume);
  double chainStep(double kineticEnergy);
  // kinetic energy of the barostat (isotropic: of x=ln(V) with mass W/9; flexible: of the cell-rate matrix)
  double barostatKineticEnergy() const;
  double energy(double volume) const;
  void reverseMomenta();

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Thermobarostat& value);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Thermobarostat& value);
};

export std::size_t thermobarostatCellDegreesOfFreedom(CellMinimizationType type);
export double3x3 projectCellTensor(const double3x3& tensor, CellMinimizationType type,
                                  MonoclinicAngleType angle = MonoclinicAngleType::Beta);
export double3x3 cellForce(const double3x3& configurationalVirial, const double3x3& kineticStress,
                          double volume, double externalPressure, double mass, CellMinimizationType type,
                          MonoclinicAngleType angle = MonoclinicAngleType::Beta);
export double sinhc(double value);
export double3x3 matrixExponential(const double3x3& matrix);
export double3x3 matrixPhi1(const double3x3& matrix);
export void propagateCellAndPosition(double3x3& cell, std::span<double3> positions,
                                     std::span<const double3> velocities, const double3x3& cellVelocity, double dt,
                                     bool upperTriangular);
export double3x3 velocityPropagator(const double3x3& cellVelocity, double dt,
                                   std::size_t translationalDegreesOfFreedom);
