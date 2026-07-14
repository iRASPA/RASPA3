module;

module property_elastic_constants_fluctuation;

import std;
import archive;
import units;
import elastic_tensor;

namespace
{
template <std::size_t N, typename Operation>
std::array<double, N> transformArray(const std::array<double, N>& a, const std::array<double, N>& b,
                                     Operation operation)
{
  std::array<double, N> result{};
  std::transform(a.begin(), a.end(), b.begin(), result.begin(), operation);
  return result;
}

template <std::size_t N, typename Operation>
std::array<double, N> transformArray(const std::array<double, N>& a, Operation operation)
{
  std::array<double, N> result{};
  std::transform(a.begin(), a.end(), result.begin(), operation);
  return result;
}

template <std::size_t N>
void addArray(std::array<double, N>& a, const std::array<double, N>& b)
{
  std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<>());
}

template <std::size_t N>
nlohmann::json convertedArray(const std::array<double, N>& values, double conversion)
{
  nlohmann::json result = nlohmann::json::array();
  for (double value : values) result.push_back(conversion * value);
  return result;
}
}  // namespace

ElasticFluctuationData operator-(const ElasticFluctuationData& a, const ElasticFluctuationData& b)
{
  ElasticFluctuationData result{};
  result.configurationalStress = transformArray(a.configurationalStress, b.configurationalStress, std::minus<>());
  result.kineticStress = transformArray(a.kineticStress, b.kineticStress, std::minus<>());
  result.born = transformArray(a.born, b.born, std::minus<>());
  result.kinetic = transformArray(a.kinetic, b.kinetic, std::minus<>());
  result.fluctuation = transformArray(a.fluctuation, b.fluctuation, std::minus<>());
  result.stiffness = transformArray(a.stiffness, b.stiffness, std::minus<>());
  result.pressureCorrection = transformArray(a.pressureCorrection, b.pressureCorrection, std::minus<>());
  result.tangentStiffness = transformArray(a.tangentStiffness, b.tangentStiffness, std::minus<>());
  return result;
}

ElasticFluctuationData operator*(const ElasticFluctuationData& a, const ElasticFluctuationData& b)
{
  ElasticFluctuationData result{};
  result.configurationalStress = transformArray(a.configurationalStress, b.configurationalStress, std::multiplies<>());
  result.kineticStress = transformArray(a.kineticStress, b.kineticStress, std::multiplies<>());
  result.born = transformArray(a.born, b.born, std::multiplies<>());
  result.kinetic = transformArray(a.kinetic, b.kinetic, std::multiplies<>());
  result.fluctuation = transformArray(a.fluctuation, b.fluctuation, std::multiplies<>());
  result.stiffness = transformArray(a.stiffness, b.stiffness, std::multiplies<>());
  result.pressureCorrection = transformArray(a.pressureCorrection, b.pressureCorrection, std::multiplies<>());
  result.tangentStiffness = transformArray(a.tangentStiffness, b.tangentStiffness, std::multiplies<>());
  return result;
}

ElasticFluctuationData operator*(double a, const ElasticFluctuationData& b)
{
  ElasticFluctuationData result{};
  result.configurationalStress = transformArray(b.configurationalStress, [a](double x) { return a * x; });
  result.kineticStress = transformArray(b.kineticStress, [a](double x) { return a * x; });
  result.born = transformArray(b.born, [a](double x) { return a * x; });
  result.kinetic = transformArray(b.kinetic, [a](double x) { return a * x; });
  result.fluctuation = transformArray(b.fluctuation, [a](double x) { return a * x; });
  result.stiffness = transformArray(b.stiffness, [a](double x) { return a * x; });
  result.pressureCorrection = transformArray(b.pressureCorrection, [a](double x) { return a * x; });
  result.tangentStiffness = transformArray(b.tangentStiffness, [a](double x) { return a * x; });
  return result;
}

ElasticFluctuationData sqrt(const ElasticFluctuationData& a)
{
  return ElasticFluctuationData{
      .configurationalStress = transformArray(a.configurationalStress, [](double x) { return std::sqrt(x); }),
      .kineticStress = transformArray(a.kineticStress, [](double x) { return std::sqrt(x); }),
      .born = transformArray(a.born, [](double x) { return std::sqrt(x); }),
      .kinetic = transformArray(a.kinetic, [](double x) { return std::sqrt(x); }),
      .fluctuation = transformArray(a.fluctuation, [](double x) { return std::sqrt(x); }),
      .stiffness = transformArray(a.stiffness, [](double x) { return std::sqrt(x); }),
      .pressureCorrection = transformArray(a.pressureCorrection, [](double x) { return std::sqrt(x); }),
      .tangentStiffness = transformArray(a.tangentStiffness, [](double x) { return std::sqrt(x); })};
}

ElasticFluctuationData& operator+=(ElasticFluctuationData& a, const ElasticFluctuationData& b)
{
  addArray(a.configurationalStress, b.configurationalStress);
  addArray(a.kineticStress, b.kineticStress);
  addArray(a.born, b.born);
  addArray(a.kinetic, b.kinetic);
  addArray(a.fluctuation, b.fluctuation);
  addArray(a.stiffness, b.stiffness);
  addArray(a.pressureCorrection, b.pressureCorrection);
  addArray(a.tangentStiffness, b.tangentStiffness);
  return a;
}

ElasticFluctuationTerms::ElasticFluctuationTerms(const std::array<double, 6>& configurationalStressValue,
                                                 const std::array<double, 6>& kineticStressValue,
                                                 const std::array<double, 36>& bornValue, double betaValue,
                                                 double volumeValue, double kineticEntitiesValue)
    : configurationalStress(configurationalStressValue),
      kineticStress(kineticStressValue),
      born(bornValue),
      beta(betaValue),
      volume(volumeValue),
      kineticEntities(kineticEntitiesValue)
{
  for (std::size_t i = 0; i < elasticVoigtSize; ++i)
    for (std::size_t j = 0; j < elasticVoigtSize; ++j)
      elasticAt(stressProducts, i, j) = configurationalStress[i] * configurationalStress[j];
}

ElasticFluctuationTerms& ElasticFluctuationTerms::operator+=(const ElasticFluctuationTerms& other)
{
  addArray(configurationalStress, other.configurationalStress);
  addArray(kineticStress, other.kineticStress);
  addArray(stressProducts, other.stressProducts);
  addArray(born, other.born);
  beta += other.beta;
  volume += other.volume;
  kineticEntities += other.kineticEntities;
  return *this;
}

ElasticFluctuationData ElasticFluctuationTerms::compositeProperty() const
{
  ElasticFluctuationData result{};
  result.configurationalStress = configurationalStress;
  result.kineticStress = kineticStress;
  result.born = born;
  if (!(beta > 0.0) || !(volume > 0.0)) return result;

  for (std::size_t i = 0; i < elasticVoigtSize; ++i)
  {
    for (std::size_t j = 0; j < elasticVoigtSize; ++j)
    {
      const double covariance = elasticAt(stressProducts, i, j) - configurationalStress[i] * configurationalStress[j];
      elasticAt(result.fluctuation, i, j) = -beta * volume * covariance;
    }
  }
  const double idealPressure = kineticEntities / (beta * volume);
  for (std::size_t normal = 0; normal < 3; ++normal) elasticAt(result.kinetic, normal, normal) = 2.0 * idealPressure;
  for (std::size_t shear = 3; shear < elasticVoigtSize; ++shear)
    elasticAt(result.kinetic, shear, shear) = idealPressure;

  for (std::size_t index = 0; index < result.stiffness.size(); ++index)
    result.stiffness[index] = result.born[index] + result.kinetic[index] + result.fluctuation[index];
  symmetrizeElasticTensor(result.stiffness);

  const double meanPressure = (configurationalStress[0] + kineticStress[0] + configurationalStress[1] +
                               kineticStress[1] + configurationalStress[2] + kineticStress[2]) /
                              3.0;
  for (std::size_t normal = 0; normal < 3; ++normal)
    elasticAt(result.pressureCorrection, normal, normal) = -meanPressure;
  for (std::size_t shear = 3; shear < elasticVoigtSize; ++shear)
    elasticAt(result.pressureCorrection, shear, shear) = -0.5 * meanPressure;
  for (std::size_t index = 0; index < result.tangentStiffness.size(); ++index)
    result.tangentStiffness[index] = result.stiffness[index] + result.pressureCorrection[index];
  return result;
}

ElasticFluctuationTerms operator*(double a, const ElasticFluctuationTerms& b)
{
  ElasticFluctuationTerms result{};
  result.configurationalStress = transformArray(b.configurationalStress, [a](double x) { return a * x; });
  result.kineticStress = transformArray(b.kineticStress, [a](double x) { return a * x; });
  result.stressProducts = transformArray(b.stressProducts, [a](double x) { return a * x; });
  result.born = transformArray(b.born, [a](double x) { return a * x; });
  result.beta = a * b.beta;
  result.volume = a * b.volume;
  result.kineticEntities = a * b.kineticEntities;
  return result;
}

ElasticFluctuationTerms operator/(const ElasticFluctuationTerms& a, double b) { return (1.0 / b) * a; }

std::string PropertyElasticConstantsFluctuation::writeAveragesStatistics(double eigenvalueTolerance) const
{
  const auto [mean, confidence] = average();
  const bool reduced = Units::unitSystem == Units::System::ReducedUnits;
  const double conversion = reduced ? 1.0 : 1.0e-9 * Units::PressureConversionFactor;
  const std::string_view unit = reduced ? "reduced pressure" : "GPa";
  const ElasticDerivedProperties derived = deriveElasticProperties(mean.tangentStiffness, eigenvalueTolerance);
  const double sampleCount = std::accumulate(bookKeeping.begin(), bookKeeping.end(), 0.0,
                                             [](double count, const auto& block) { return count + block.second; });
  std::string output = std::format(
      "\nIsothermal elastic constants from NVT stress fluctuations\n"
      "===============================================================================\n"
      "Voigt order: xx yy zz yz xz xy; shear strains use engineering convention.\n"
      "Samples: {:.0f}; blocks: {}\n\n"
      "Mean configurational stress [{}]: {: .7e} {: .7e} {: .7e} {: .7e} {: .7e} {: .7e}\n"
      "Mean kinetic stress [{}]:       {: .7e} {: .7e} {: .7e} {: .7e} {: .7e} {: .7e}\n\n"
      "Affine Born term [{}]\n{}\n"
      "Kinetic term [{}]\n{}\n"
      "Stress-covariance term [{}]\n{}\n"
      "Helmholtz isothermal stiffness [{}]\n{}\n"
      "Hydrostatic prestress correction [{}]\n{}\n"
      "Prestress-corrected tangent stiffness [{}]\n{}\n"
      "Tangent-stiffness confidence interval [{}]\n{}\n",
      sampleCount, numberOfBlocks, unit, conversion * mean.configurationalStress[0],
      conversion * mean.configurationalStress[1], conversion * mean.configurationalStress[2],
      conversion * mean.configurationalStress[3], conversion * mean.configurationalStress[4],
      conversion * mean.configurationalStress[5], unit, conversion * mean.kineticStress[0],
      conversion * mean.kineticStress[1], conversion * mean.kineticStress[2], conversion * mean.kineticStress[3],
      conversion * mean.kineticStress[4], conversion * mean.kineticStress[5], unit,
      elasticMatrixString(mean.born, conversion), unit, elasticMatrixString(mean.kinetic, conversion), unit,
      elasticMatrixString(mean.fluctuation, conversion), unit, elasticMatrixString(mean.stiffness, conversion), unit,
      elasticMatrixString(mean.pressureCorrection, conversion), unit,
      elasticMatrixString(mean.tangentStiffness, conversion), unit,
      elasticMatrixString(confidence.tangentStiffness, conversion));
  if (derived.complianceAvailable)
  {
    const double complianceConversion = reduced ? 1.0 : 1.0 / conversion;
    output += std::format(
        "Tangent compliance [{}]\n{}\n"
        "Young moduli: {: .7e} {: .7e} {: .7e} [{}]\n"
        "Bulk modulus (Voigt/Reuss/Hill): {: .7e} {: .7e} {: .7e} [{}]\n"
        "Shear modulus (Voigt/Reuss/Hill): {: .7e} {: .7e} {: .7e} [{}]\n",
        reduced ? "inverse reduced pressure" : "GPa^-1", elasticMatrixString(derived.compliance, complianceConversion),
        conversion * derived.youngModuli[0], conversion * derived.youngModuli[1], conversion * derived.youngModuli[2],
        unit, conversion * derived.bulkModulusVoigt, conversion * derived.bulkModulusReuss,
        conversion * derived.bulkModulusHill, unit, conversion * derived.shearModulusVoigt,
        conversion * derived.shearModulusReuss, conversion * derived.shearModulusHill, unit);
  }
  else
    output += "Compliance and derived Reuss/Hill moduli are unavailable because the tensor is singular.\n";
  if (sampleCount < 3.0 * static_cast<double>(numberOfBlocks))
    output += "Warning: fewer than three elastic samples per block; confidence intervals are not reliable.\n";
  return output;
}

nlohmann::json PropertyElasticConstantsFluctuation::jsonAveragesStatistics(double eigenvalueTolerance) const
{
  const auto [mean, confidence] = average();
  const double conversion =
      Units::unitSystem == Units::System::ReducedUnits ? 1.0 : 1.0e-9 * Units::PressureConversionFactor;
  const ElasticDerivedProperties derived = deriveElasticProperties(mean.tangentStiffness, eigenvalueTolerance);
  nlohmann::json status;
  status["voigtOrder"] = {"xx", "yy", "zz", "yz", "xz", "xy"};
  status["configurationalStress"] = convertedArray(mean.configurationalStress, conversion);
  status["kineticStress"] = convertedArray(mean.kineticStress, conversion);
  status["born"] = convertedArray(mean.born, conversion);
  status["kinetic"] = convertedArray(mean.kinetic, conversion);
  status["fluctuation"] = convertedArray(mean.fluctuation, conversion);
  status["stiffness"] = convertedArray(mean.stiffness, conversion);
  status["pressureCorrection"] = convertedArray(mean.pressureCorrection, conversion);
  status["tangentStiffness"] = convertedArray(mean.tangentStiffness, conversion);
  status["tangentStiffnessConfidence"] = convertedArray(confidence.tangentStiffness, conversion);
  status["complianceAvailable"] = derived.complianceAvailable;
  status["stabilityEigenvalues"] = convertedArray(derived.stabilityEigenvalues, conversion);
  if (derived.complianceAvailable)
  {
    const double complianceConversion = Units::unitSystem == Units::System::ReducedUnits ? 1.0 : 1.0 / conversion;
    status["compliance"] = convertedArray(derived.compliance, complianceConversion);
    status["youngModuli"] = convertedArray(derived.youngModuli, conversion);
    status["poissonRatios"] = derived.poissonRatios;
    status["bulkModulus"] = {{"voigt", conversion * derived.bulkModulusVoigt},
                             {"reuss", conversion * derived.bulkModulusReuss},
                             {"hill", conversion * derived.bulkModulusHill}};
    status["shearModulus"] = {{"voigt", conversion * derived.shearModulusVoigt},
                              {"reuss", conversion * derived.shearModulusReuss},
                              {"hill", conversion * derived.shearModulusHill}};
  }
  return status;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ElasticFluctuationTerms& terms)
{
  archive << terms.versionNumber << terms.configurationalStress << terms.kineticStress << terms.stressProducts
          << terms.born << terms.beta << terms.volume << terms.kineticEntities;
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ElasticFluctuationTerms& terms)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > terms.versionNumber)
    throw std::runtime_error("Invalid PropertyElasticConstantsFluctuation archive version");
  archive >> terms.configurationalStress >> terms.kineticStress >> terms.stressProducts >> terms.born >> terms.beta >>
      terms.volume >> terms.kineticEntities;
  return archive;
}
