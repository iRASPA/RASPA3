module;

module thermobarostat;

import std;

import archive;
import double3;
import double3x3;
import randomnumbers;
import units;
import minimization_cell_layout;

namespace
{
bool equalCaseInsensitive(std::string_view lhs, std::string_view rhs)
{
  return lhs.size() == rhs.size() &&
         std::ranges::equal(lhs, rhs, [](char a, char b) {
           return std::tolower(static_cast<unsigned char>(a)) == std::tolower(static_cast<unsigned char>(b));
         });
}

double matrixNorm(const double3x3& matrix)
{
  double result{};
  for (std::size_t column = 0; column != 3; ++column)
  {
    double sum{};
    for (std::size_t row = 0; row != 3; ++row) sum += std::abs(matrix.mm[column][row]);
    result = std::max(result, sum);
  }
  return result;
}

std::vector<double> yoshidaSuzukiWeights(std::size_t count)
{
  switch (count)
  {
    case 1:
      return {1.0};
    case 3:
    {
      const double a = 1.0 / (2.0 - std::cbrt(2.0));
      return {a, 1.0 - 2.0 * a, a};
    }
    case 5:
    {
      const double a = 1.0 / (4.0 - std::cbrt(4.0));
      return {a, a, 1.0 - 4.0 * a, a, a};
    }
    case 7:
      return {0.784513610477560, 0.235573213359357, -1.17767998417887,
              1.0 - 2.0 * (0.784513610477560 + 0.235573213359357 - 1.17767998417887),
              -1.17767998417887, 0.235573213359357, 0.784513610477560};
    case 9:
      return {0.192, 0.5549108184097836197, 0.1246596199418886442, -0.8431820635969335053,
              1.0 - 2.0 * (0.192 + 0.5549108184097836197 + 0.1246596199418886442 -
                           0.8431820635969335053),
              -0.8431820635969335053, 0.1246596199418886442, 0.5549108184097836197, 0.192};
    default:
      throw std::runtime_error("Yoshida-Suzuki steps must be 1, 3, 5, 7, or 9");
  }
}

void zeroLowerTriangle(double3x3& matrix)
{
  matrix.ay = 0.0;
  matrix.az = 0.0;
  matrix.bz = 0.0;
}
}  // namespace

std::optional<MolecularDynamicsEnsemble> molecularDynamicsEnsembleFromString(std::string_view value)
{
  for (MolecularDynamicsEnsemble ensemble : {MolecularDynamicsEnsemble::NVE, MolecularDynamicsEnsemble::NVT,
                                             MolecularDynamicsEnsemble::NPT, MolecularDynamicsEnsemble::NPTPR,
                                             MolecularDynamicsEnsemble::MuVT, MolecularDynamicsEnsemble::MuPT,
                                             MolecularDynamicsEnsemble::MuPTPR})
  {
    if (equalCaseInsensitive(value, molecularDynamicsEnsembleName(ensemble))) return ensemble;
  }
  return std::nullopt;
}

std::string molecularDynamicsEnsembleName(MolecularDynamicsEnsemble ensemble)
{
  switch (ensemble)
  {
    case MolecularDynamicsEnsemble::NVE:
      return "NVE";
    case MolecularDynamicsEnsemble::NVT:
      return "NVT";
    case MolecularDynamicsEnsemble::NPT:
      return "NPT";
    case MolecularDynamicsEnsemble::NPTPR:
      return "NPTPR";
    case MolecularDynamicsEnsemble::MuVT:
      return "MuVT";
    case MolecularDynamicsEnsemble::MuPT:
      return "MuPT";
    case MolecularDynamicsEnsemble::MuPTPR:
      return "MuPTPR";
  }
  std::unreachable();
}

bool molecularDynamicsUsesThermostat(MolecularDynamicsEnsemble ensemble)
{
  return ensemble != MolecularDynamicsEnsemble::NVE;
}

bool molecularDynamicsUsesIsotropicBarostat(MolecularDynamicsEnsemble ensemble)
{
  return ensemble == MolecularDynamicsEnsemble::NPT || ensemble == MolecularDynamicsEnsemble::MuPT;
}

bool molecularDynamicsUsesFlexibleBarostat(MolecularDynamicsEnsemble ensemble)
{
  return ensemble == MolecularDynamicsEnsemble::NPTPR || ensemble == MolecularDynamicsEnsemble::MuPTPR;
}

bool molecularDynamicsHasParticleExchange(MolecularDynamicsEnsemble ensemble)
{
  return ensemble == MolecularDynamicsEnsemble::MuVT || ensemble == MolecularDynamicsEnsemble::MuPT ||
         ensemble == MolecularDynamicsEnsemble::MuPTPR;
}

std::size_t thermobarostatCellDegreesOfFreedom(CellMinimizationType type)
{
  switch (type)
  {
    case CellMinimizationType::Isotropic:
      return 1;
    case CellMinimizationType::Anisotropic:
      return 3;
    case CellMinimizationType::Monoclinic:
    case CellMinimizationType::MonoclinicUpperTriangle:
      return 4;
    case CellMinimizationType::Regular:
    case CellMinimizationType::RegularUpperTriangle:
      return 6;
    case CellMinimizationType::Fixed:
      return 0;
  }
  std::unreachable();
}

Thermobarostat::Thermobarostat(MolecularDynamicsEnsemble ensemble, CellMinimizationType cellType,
                               MonoclinicAngleType monoclinicAngle, double temperature, double pressure,
                               double timeStep, std::size_t translationalDegreesOfFreedom, std::size_t chainLength,
                               std::size_t numberOfYoshidaSuzukiSteps, double timeScaleParameterBarostat)
    : ensemble(ensemble),
      cellType(molecularDynamicsUsesIsotropicBarostat(ensemble) ? CellMinimizationType::Isotropic : cellType),
      monoclinicAngle(monoclinicAngle),
      temperature(temperature),
      pressure(pressure),
      timeStep(timeStep),
      timeScaleParameterBarostat(timeScaleParameterBarostat),
      translationalDegreesOfFreedom(translationalDegreesOfFreedom),
      cellDegreesOfFreedom(thermobarostatCellDegreesOfFreedom(this->cellType)),
      chainLength(chainLength),
      numberOfYoshidaSuzukiSteps(numberOfYoshidaSuzukiSteps),
      chainDegreesOfFreedom(chainLength),
      chainForce(chainLength),
      chainVelocity(chainLength),
      chainPosition(chainLength),
      chainMass(chainLength),
      yoshidaSuzukiWeights(::yoshidaSuzukiWeights(numberOfYoshidaSuzukiSteps))
{
  if (!molecularDynamicsUsesIsotropicBarostat(ensemble) && !molecularDynamicsUsesFlexibleBarostat(ensemble))
    throw std::runtime_error("Thermobarostat requires an NPT, NPTPR, MuPT, or MuPTPR ensemble");
  if (molecularDynamicsUsesFlexibleBarostat(ensemble) && cellDegreesOfFreedom == 0)
    throw std::runtime_error("NPTPR/MuPTPR requires a variable CellType");
  if (!(temperature > 0.0) || !std::isfinite(temperature))
    throw std::runtime_error("Thermobarostat temperature must be positive and finite");
  if (!(timeStep > 0.0) || !std::isfinite(timeStep))
    throw std::runtime_error("Thermobarostat timestep must be positive and finite");
  if (!(timeScaleParameterBarostat > 0.0) || !std::isfinite(timeScaleParameterBarostat))
    throw std::runtime_error("Thermobarostat time scale must be positive and finite");
  if (!std::isfinite(pressure)) throw std::runtime_error("Thermobarostat pressure must be finite");
  if (chainLength == 0) throw std::runtime_error("BarostatChainLength must be positive");
  const double particleScale = static_cast<double>(translationalDegreesOfFreedom + 3) * Units::KB * temperature;
  logVolumeMass = particleScale * timeScaleParameterBarostat * timeScaleParameterBarostat;
  cellMass = particleScale * timeScaleParameterBarostat * timeScaleParameterBarostat / 3.0;
}

void Thermobarostat::refreshDegreesOfFreedom(RandomNumber& random, std::size_t newTranslationalDegreesOfFreedom,
                                             double volume)
{
  if (!(volume > 0.0) || !std::isfinite(volume)) throw std::runtime_error("Thermobarostat requires positive volume");
  translationalDegreesOfFreedom = newTranslationalDegreesOfFreedom;

  const double particleScale = static_cast<double>(translationalDegreesOfFreedom + 3) * Units::KB * temperature;
  logVolumeMass = particleScale * timeScaleParameterBarostat * timeScaleParameterBarostat;
  cellMass = particleScale * timeScaleParameterBarostat * timeScaleParameterBarostat / 3.0;
  logVolumePosition = std::log(volume);

  if (molecularDynamicsUsesFlexibleBarostat(ensemble))
  {
    for (std::size_t column = 0; column != 3; ++column)
      for (std::size_t row = 0; row != 3; ++row)
        cellVelocity.mm[column][row] = random.Gaussian() * std::sqrt(Units::KB * temperature / cellMass);
    cellVelocity = projectCellTensor(cellVelocity, cellType, monoclinicAngle);
  }
}

void Thermobarostat::initialize(RandomNumber& random)
{
  if (chainLength == 0) throw std::runtime_error("BarostatChainLength must be positive");
  if (numberOfRespaSteps == 0) throw std::runtime_error("NumberOfRespaSteps must be positive");
  refreshDegreesOfFreedom(random, translationalDegreesOfFreedom, std::exp(logVolumePosition));
  if (molecularDynamicsUsesIsotropicBarostat(ensemble))
    logVolumeVelocity = random.Gaussian() * std::sqrt(Units::KB * temperature / logVolumeMass);

  chainDegreesOfFreedom[0] = static_cast<double>(cellDegreesOfFreedom) * Units::KB * temperature;
  for (std::size_t i = 1; i != chainLength; ++i) chainDegreesOfFreedom[i] = Units::KB * temperature;
  for (std::size_t i = 0; i != chainLength; ++i)
  {
    chainMass[i] = chainDegreesOfFreedom[i] * timeScaleParameterBarostat * timeScaleParameterBarostat;
    chainVelocity[i] = random.Gaussian() * std::sqrt(chainDegreesOfFreedom[i] / chainMass[i]);
  }
}

double Thermobarostat::chainStep(double kineticEnergy)
{
  if (chainLength == 0) return 1.0;
  chainForce[0] = (2.0 * kineticEnergy - chainDegreesOfFreedom[0]) / chainMass[0];
  double scale = 1.0;
  for (std::size_t respa = 0; respa != numberOfRespaSteps; ++respa)
  {
    for (double weight : yoshidaSuzukiWeights)
    {
      const double quarter = weight * timeStep / (4.0 * static_cast<double>(numberOfRespaSteps));
      chainVelocity.back() += chainForce.back() * quarter;
      for (std::size_t k = chainLength - 1; k != 0; --k)
      {
        const double attenuation = std::exp(-0.5 * quarter * chainVelocity[k]);
        chainVelocity[k - 1] = chainVelocity[k - 1] * attenuation * attenuation +
                               chainForce[k - 1] * attenuation * quarter;
      }
      const double attenuation = std::exp(-2.0 * quarter * chainVelocity[0]);
      scale *= attenuation;
      chainForce[0] = (2.0 * kineticEnergy * scale * scale - chainDegreesOfFreedom[0]) / chainMass[0];
      for (std::size_t k = 0; k != chainLength; ++k)
        chainPosition[k] += 2.0 * quarter * chainVelocity[k];
      for (std::size_t k = 0; k + 1 != chainLength; ++k)
      {
        const double coupling = std::exp(-0.5 * quarter * chainVelocity[k + 1]);
        chainVelocity[k] = chainVelocity[k] * coupling * coupling + chainForce[k] * coupling * quarter;
        chainForce[k + 1] =
            (chainMass[k] * chainVelocity[k] * chainVelocity[k] - chainDegreesOfFreedom[k + 1]) /
            chainMass[k + 1];
      }
      chainVelocity.back() += chainForce.back() * quarter;
    }
  }
  return scale;
}

double Thermobarostat::energy(double volume) const
{
  double result = pressure * volume;
  if (molecularDynamicsUsesIsotropicBarostat(ensemble))
    result += 0.5 * logVolumeMass * logVolumeVelocity * logVolumeVelocity;
  else
  {
    double normSquared{};
    for (std::size_t column = 0; column != 3; ++column)
      for (std::size_t row = 0; row != 3; ++row)
        normSquared += cellVelocity.mm[column][row] * cellVelocity.mm[column][row];
    result += 0.5 * cellMass * normSquared;
  }
  for (std::size_t i = 0; i != chainLength; ++i)
    result += 0.5 * chainMass[i] * chainVelocity[i] * chainVelocity[i] +
              chainDegreesOfFreedom[i] * chainPosition[i];
  return result;
}

void Thermobarostat::reverseMomenta()
{
  logVolumeVelocity = -logVolumeVelocity;
  cellVelocity = -cellVelocity;
  for (double& velocity : chainVelocity) velocity = -velocity;
}

double3x3 projectCellTensor(const double3x3& tensor, CellMinimizationType type, MonoclinicAngleType angle)
{
  double3x3 result{};
  switch (type)
  {
    case CellMinimizationType::Fixed:
      return result;
    case CellMinimizationType::Isotropic:
    {
      const double average = tensor.trace() / 3.0;
      return double3x3(average, average, average);
    }
    case CellMinimizationType::Anisotropic:
      return double3x3(tensor.ax, tensor.by, tensor.cz);
    case CellMinimizationType::Regular:
      result = 0.5 * (tensor + tensor.transpose());
      return result;
    case CellMinimizationType::RegularUpperTriangle:
      result = tensor;
      zeroLowerTriangle(result);
      return result;
    case CellMinimizationType::Monoclinic:
      result = double3x3(tensor.ax, tensor.by, tensor.cz);
      if (angle == MonoclinicAngleType::Alpha) result.bz = result.cy = 0.5 * (tensor.bz + tensor.cy);
      if (angle == MonoclinicAngleType::Beta) result.az = result.cx = 0.5 * (tensor.az + tensor.cx);
      if (angle == MonoclinicAngleType::Gamma) result.ay = result.bx = 0.5 * (tensor.ay + tensor.bx);
      return result;
    case CellMinimizationType::MonoclinicUpperTriangle:
      result = double3x3(tensor.ax, tensor.by, tensor.cz);
      if (angle == MonoclinicAngleType::Alpha) result.cy = tensor.cy;
      if (angle == MonoclinicAngleType::Beta) result.cx = tensor.cx;
      if (angle == MonoclinicAngleType::Gamma) result.bx = tensor.bx;
      return result;
  }
  std::unreachable();
}

double3x3 cellForce(const double3x3& configurationalVirial, const double3x3& kineticStress, double volume,
                    double externalPressure, double mass, CellMinimizationType type, MonoclinicAngleType angle)
{
  double3x3 force = configurationalVirial + kineticStress -
                    double3x3(externalPressure * volume, externalPressure * volume, externalPressure * volume);
  return (1.0 / mass) * projectCellTensor(force, type, angle);
}

double sinhc(double value)
{
  if (std::abs(value) < 1.0e-6)
  {
    const double square = value * value;
    return 1.0 + square / 6.0 + square * square / 120.0;
  }
  return std::sinh(value) / value;
}

double3x3 matrixExponential(const double3x3& matrix)
{
  const double norm = matrixNorm(matrix);
  const int scaling = norm > 0.5 ? std::max(0, static_cast<int>(std::ceil(std::log2(norm / 0.5)))) : 0;
  const double3x3 scaled = matrix * std::ldexp(1.0, -scaling);
  double3x3 result = double3x3::identity();
  double3x3 term = result;
  for (std::size_t order = 1; order != 30; ++order)
  {
    term = term * scaled * (1.0 / static_cast<double>(order));
    result += term;
    if (matrixNorm(term) < 1.0e-16) break;
  }
  for (int i = 0; i != scaling; ++i) result = result * result;
  return result;
}

double3x3 matrixPhi1(const double3x3& matrix)
{
  double3x3 result = double3x3::identity();
  double3x3 power = double3x3::identity();
  double factorial = 1.0;
  for (std::size_t order = 1; order != 40; ++order)
  {
    power = power * matrix;
    factorial *= static_cast<double>(order + 1);
    const double3x3 term = power * (1.0 / factorial);
    result += term;
    if (matrixNorm(term) < 1.0e-16) break;
  }
  return result;
}

void propagateCellAndPosition(double3x3& cell, std::span<double3> positions, std::span<const double3> velocities,
                              const double3x3& cellVelocity, double dt, bool upperTriangular)
{
  if (positions.size() != velocities.size()) throw std::runtime_error("Position/velocity span size mismatch");
  const double3x3 argument = dt * cellVelocity;
  double3x3 exponential = matrixExponential(argument);
  double3x3 drift = dt * matrixPhi1(argument);
  for (std::size_t i = 0; i != positions.size(); ++i) positions[i] = exponential * positions[i] + drift * velocities[i];
  cell = exponential * cell;
  if (upperTriangular)
  {
    zeroLowerTriangle(cell);
    zeroLowerTriangle(exponential);
  }
}

double3x3 velocityPropagator(const double3x3& cellVelocity, double dt,
                             std::size_t translationalDegreesOfFreedom)
{
  double3x3 coupling = cellVelocity;
  const double correction = cellVelocity.trace() /
                            static_cast<double>(std::max<std::size_t>(1, translationalDegreesOfFreedom));
  coupling.ax += correction;
  coupling.by += correction;
  coupling.cz += correction;
  return matrixExponential(-dt * coupling);
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Thermobarostat& value)
{
  archive << value.versionNumber;
  archive << static_cast<std::uint8_t>(value.ensemble);
  archive << static_cast<std::uint8_t>(value.cellType);
  archive << static_cast<std::uint8_t>(value.monoclinicAngle);
  archive << value.temperature << value.pressure << value.timeStep << value.timeScaleParameterBarostat;
  archive << value.translationalDegreesOfFreedom << value.cellDegreesOfFreedom << value.chainLength;
  archive << value.numberOfRespaSteps << value.numberOfYoshidaSuzukiSteps;
  archive << value.logVolumePosition << value.logVolumeVelocity << value.logVolumeMass;
  archive << value.cellVelocity << value.cellMass;
  archive << value.chainDegreesOfFreedom << value.chainForce << value.chainVelocity << value.chainPosition;
  archive << value.chainMass << value.yoshidaSuzukiWeights;
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Thermobarostat& value)
{
  std::uint64_t version{};
  std::uint8_t ensemble{}, cellType{}, angle{};
  archive >> version;
  if (version > value.versionNumber) throw std::runtime_error("Unsupported Thermobarostat restart version");
  archive >> ensemble >> cellType >> angle;
  value.ensemble = static_cast<MolecularDynamicsEnsemble>(ensemble);
  value.cellType = static_cast<CellMinimizationType>(cellType);
  value.monoclinicAngle = static_cast<MonoclinicAngleType>(angle);
  archive >> value.temperature >> value.pressure >> value.timeStep >> value.timeScaleParameterBarostat;
  archive >> value.translationalDegreesOfFreedom >> value.cellDegreesOfFreedom >> value.chainLength;
  archive >> value.numberOfRespaSteps >> value.numberOfYoshidaSuzukiSteps;
  archive >> value.logVolumePosition >> value.logVolumeVelocity >> value.logVolumeMass;
  archive >> value.cellVelocity >> value.cellMass;
  archive >> value.chainDegreesOfFreedom >> value.chainForce >> value.chainVelocity >> value.chainPosition;
  archive >> value.chainMass >> value.yoshidaSuzukiWeights;
  if (!molecularDynamicsUsesIsotropicBarostat(value.ensemble) &&
      !molecularDynamicsUsesFlexibleBarostat(value.ensemble))
    throw std::runtime_error("Invalid Thermobarostat ensemble in restart");
  if (value.cellDegreesOfFreedom != thermobarostatCellDegreesOfFreedom(value.cellType))
    throw std::runtime_error("Inconsistent Thermobarostat cell degrees of freedom in restart");
  if (value.chainLength == 0 || value.numberOfRespaSteps == 0)
    throw std::runtime_error("Invalid Thermobarostat chain configuration in restart");
  if (value.chainDegreesOfFreedom.size() != value.chainLength || value.chainForce.size() != value.chainLength ||
      value.chainVelocity.size() != value.chainLength || value.chainPosition.size() != value.chainLength ||
      value.chainMass.size() != value.chainLength)
    throw std::runtime_error("Inconsistent Thermobarostat chain lengths in restart");
  if (value.yoshidaSuzukiWeights.size() != value.numberOfYoshidaSuzukiSteps)
    throw std::runtime_error("Inconsistent Thermobarostat Yoshida-Suzuki weights in restart");
  return archive;
}
