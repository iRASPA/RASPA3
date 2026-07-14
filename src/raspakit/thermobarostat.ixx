module;

export module thermobarostat;

import std;

import archive;
import double3;
import double3x3;
import randomnumbers;
import units;
import minimization_cell_layout;

export enum class MolecularDynamicsEnsemble : std::uint8_t { NVE, NVT, NPT, NPTPR };

export std::optional<MolecularDynamicsEnsemble> molecularDynamicsEnsembleFromString(std::string_view value);
export std::string molecularDynamicsEnsembleName(MolecularDynamicsEnsemble ensemble);

export struct Thermobarostat
{
  std::uint64_t versionNumber{1};
  MolecularDynamicsEnsemble ensemble{MolecularDynamicsEnsemble::NVE};
  CellMinimizationType cellType{CellMinimizationType::Isotropic};
  MonoclinicAngleType monoclinicAngle{MonoclinicAngleType::Beta};
  double temperature{300.0};
  double pressure{};
  double timeStep{0.0005};
  double timeScaleParameterBarostat{1.0};
  std::size_t translationalDegreesOfFreedom{};
  std::size_t cellDegreesOfFreedom{1};
  std::size_t chainLength{5};
  std::size_t numberOfRespaSteps{5};
  std::size_t numberOfYoshidaSuzukiSteps{5};

  // NPT uses eta=ln(V); NPTPR uses the logarithmic cell-rate matrix.
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
  double chainStep(double kineticEnergy);
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
