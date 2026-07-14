module;

export module potential_vdw_rare_derivatives;

import std;

import forcefield;
import vdwparameters;

namespace Potentials::Detail
{
export struct RareVDWDerivatives
{
  double energy;
  double dUdlambda;
  double radialFirstDerivative;
  double radialSecondDerivative;
};

namespace Internal
{
struct RadialJet
{
  double value;
  double first;
  double second;

  constexpr RadialJet(double value_, double first_ = 0.0, double second_ = 0.0)
      : value(value_), first(first_), second(second_)
  {
  }

  static constexpr RadialJet variable(double value) { return RadialJet(value, 1.0, 0.0); }
};

constexpr RadialJet operator+(const RadialJet& a, const RadialJet& b)
{
  return RadialJet(a.value + b.value, a.first + b.first, a.second + b.second);
}

constexpr RadialJet operator-(const RadialJet& a, const RadialJet& b)
{
  return RadialJet(a.value - b.value, a.first - b.first, a.second - b.second);
}

constexpr RadialJet operator-(const RadialJet& a)
{
  return RadialJet(-a.value, -a.first, -a.second);
}

constexpr RadialJet operator*(const RadialJet& a, const RadialJet& b)
{
  return RadialJet(a.value * b.value, a.first * b.value + a.value * b.first,
                   a.second * b.value + 2.0 * a.first * b.first + a.value * b.second);
}

inline RadialJet inverse(const RadialJet& a)
{
  double inverseValue = 1.0 / a.value;
  double inverseSquared = inverseValue * inverseValue;
  return RadialJet(inverseValue, -a.first * inverseSquared,
                   2.0 * a.first * a.first * inverseSquared * inverseValue - a.second * inverseSquared);
}

inline RadialJet operator/(const RadialJet& a, const RadialJet& b)
{
  return a * inverse(b);
}

inline RadialJet exponential(const RadialJet& a)
{
  double value = std::exp(a.value);
  return RadialJet(value, value * a.first, value * (a.second + a.first * a.first));
}

inline RadialJet power(const RadialJet& a, double exponent)
{
  double value = std::pow(a.value, exponent);
  double factor = exponent * std::pow(a.value, exponent - 1.0);
  double curvature = exponent * (exponent - 1.0) * std::pow(a.value, exponent - 2.0);
  return RadialJet(value, factor * a.first, curvature * a.first * a.first + factor * a.second);
}

constexpr RadialJet compose(const RadialJet& outer, const RadialJet& inner)
{
  return RadialJet(outer.value, outer.first * inner.first,
                   outer.second * inner.first * inner.first + outer.first * inner.second);
}

inline RadialJet softenedRadius(const RadialJet& radius, double inverseScaling, double w6)
{
  return power(power(radius, 6.0) + 0.5 * inverseScaling * inverseScaling * w6, 1.0 / 6.0);
}

inline RareVDWDerivatives softCoreResult(const RadialJet& radialTerm, const RadialJet& softRadius,
                                         double scaling, double inverseScaling, double w6)
{
  RadialJet spatialTerm = compose(radialTerm, softRadius);
  double softRadius5 = std::pow(softRadius.value, 5);
  double dUdlambda =
      radialTerm.value - scaling * inverseScaling * w6 * radialTerm.first / (6.0 * softRadius5);
  return RareVDWDerivatives{scaling * radialTerm.value, dUdlambda, scaling * spatialTerm.first,
                            scaling * spatialTerm.second};
}

inline RareVDWDerivatives linearScalingResult(const RadialJet& radialTerm, double scaling)
{
  return RareVDWDerivatives{scaling * radialTerm.value, radialTerm.value, scaling * radialTerm.first,
                            scaling * radialTerm.second};
}

inline RadialJet tangToennies(std::size_t order, const RadialJet& br)
{
  RadialJet sum(1.0);
  RadialJet term(1.0);
  for (std::size_t k = 1; k <= order; ++k)
  {
    term = term * br / static_cast<double>(k);
    sum = sum + term;
  }
  return 1.0 - exponential(-br) * sum;
}
}  // namespace Internal

export [[clang::noinline]] [[clang::preserve_most]] inline RareVDWDerivatives evaluateRareVDWDerivatives(
    const ForceField& forcefield, double scalingA, double scalingB, double rr, std::size_t typeA,
    std::size_t typeB, VDWParameters::Type potentialType)
{
  using namespace Internal;

  const VDWParameters& p = forcefield(typeA, typeB);
  double scaling = scalingA * scalingB;
  double inverseScaling = 1.0 - scaling;
  RadialJet radius = RadialJet::variable(std::sqrt(rr));

  if (potentialType == VDWParameters::Type::Morse)
  {
    RadialJet exponentialTerm = exponential(-p.parameters.y * (radius - p.parameters.z));
    RadialJet term = p.parameters.x * ((1.0 - exponentialTerm) * (1.0 - exponentialTerm) - 1.0) - p.shift;
    return linearScalingResult(term, scaling);
  }

  if (potentialType == VDWParameters::Type::MatsuokaClementiYoshimine)
  {
    RadialJet term = p.parameters.x * exponential(-p.parameters.y * radius) +
                     p.parameters.z * exponential(-p.parameters.w * radius) - p.shift;
    return linearScalingResult(term, scaling);
  }

  if (potentialType == VDWParameters::Type::RepulsiveHarmonic)
  {
    double range = p.parameters.y;
    if (radius.value >= range) return RareVDWDerivatives{0.0, 0.0, 0.0, 0.0};
    RadialJet displacement = 1.0 - radius / range;
    RadialJet term = 0.5 * p.parameters.x * displacement * displacement;
    return RareVDWDerivatives{term.value, 0.0, term.first, term.second};
  }

  if (potentialType == VDWParameters::Type::WeeksChandlerAndersen)
  {
    double sigma2 = p.parameters.y * p.parameters.y;
    if (rr > std::cbrt(2.0) * sigma2) return RareVDWDerivatives{0.0, 0.0, 0.0, 0.0};

    double sigma6 = sigma2 * sigma2 * sigma2;
    double eps4 = 4.0 * p.parameters.x;
    RadialJet softRadius = softenedRadius(radius, inverseScaling, sigma6);
    RadialJet radial = RadialJet::variable(softRadius.value);
    RadialJet u6 = sigma6 / power(radial, 6.0);
    double gcut = 1.0 / (2.0 + 0.5 * inverseScaling * inverseScaling);
    RadialJet term = eps4 * (u6 * u6 - u6) - eps4 * (gcut * gcut - gcut);
    RareVDWDerivatives result = softCoreResult(term, softRadius, scaling, inverseScaling, sigma6);
    result.dUdlambda -= scaling * eps4 * inverseScaling * (2.0 * gcut - 1.0) * gcut * gcut;
    return result;
  }

  double w6 = p.parameters2.w;
  RadialJet softRadius = softenedRadius(radius, inverseScaling, w6);
  RadialJet r = RadialJet::variable(softRadius.value);
  RadialJet term(0.0);

  switch (potentialType)
  {
    case VDWParameters::Type::BuckingHam:
      term = p.parameters.x * exponential(-p.parameters.y * r) - p.parameters.z / power(r, 6.0) - p.shift;
      break;
    case VDWParameters::Type::FeynmannHibbs:
    {
      double eps4 = 4.0 * p.parameters.x;
      RadialJet u6 = std::pow(p.parameters.y, 6) / power(r, 6.0);
      term = eps4 * (u6 * u6 - u6) +
             eps4 * p.parameters2.x * (132.0 * u6 * u6 - 30.0 * u6) / power(r, 2.0) - p.shift;
      break;
    }
    case VDWParameters::Type::MM3:
    {
      RadialJet ratio = p.parameters.y / r;
      if (ratio.value > 3.02)
      {
        term = p.parameters.x * 192.270 * ratio * ratio - p.shift;
      }
      else
      {
        term = p.parameters.x *
                   (1.84e5 * exponential(-12.0 * r / p.parameters.y) - 2.25 * power(ratio, 6.0)) -
               p.shift;
      }
      break;
    }
    case VDWParameters::Type::BornHugginsMeyer:
      term = p.parameters.x * exponential(p.parameters.y * (p.parameters.z - r)) -
             p.parameters.w / power(r, 6.0) - p.parameters2.x / power(r, 8.0) - p.shift;
      break;
    case VDWParameters::Type::LennardJonesShiftedForce:
    {
      double eps4 = 4.0 * p.parameters.x;
      double c6 = p.parameters2.x;
      double cutoff = p.parameters2.y;
      double linearCoefficient = 12.0 * c6 * c6 - 6.0 * c6;
      RadialJet u6 = std::pow(p.parameters.y, 6) / power(r, 6.0);
      term = eps4 * (u6 * u6 - u6 - c6 * (c6 - 1.0) + linearCoefficient * (r - cutoff) / cutoff);
      break;
    }
    case VDWParameters::Type::LennardJonesSecondOrderTaylorShifted:
    {
      double eps4 = 4.0 * p.parameters.x;
      double c6 = p.parameters2.x;
      double cutoff = p.parameters2.y;
      double linearCoefficient = 12.0 * c6 * c6 - 6.0 * c6;
      double quadraticCoefficient = 156.0 * c6 * c6 - 42.0 * c6;
      RadialJet displacement = (r - cutoff) / cutoff;
      RadialJet u6 = std::pow(p.parameters.y, 6) / power(r, 6.0);
      term = eps4 * (u6 * u6 - u6 - c6 * (c6 - 1.0) + linearCoefficient * displacement -
                     0.5 * quadraticCoefficient * displacement * displacement);
      break;
    }
    case VDWParameters::Type::Potential12_6:
      term = p.parameters.x / power(r, 12.0) - p.parameters.y / power(r, 6.0) - p.shift;
      break;
    case VDWParameters::Type::Potential12_6_2_0:
      term = p.parameters.x / power(r, 12.0) + p.parameters.y / power(r, 6.0) +
             p.parameters.z / power(r, 2.0) + p.parameters.w - p.shift;
      break;
    case VDWParameters::Type::CFF9_6:
      term = p.parameters.x / power(r, 9.0) - p.parameters.y / power(r, 6.0) - p.shift;
      break;
    case VDWParameters::Type::CFFEpsilonSigma:
    {
      RadialJet ratio = p.parameters.y / r;
      term = p.parameters.x * (2.0 * power(ratio, 9.0) - 3.0 * power(ratio, 6.0)) - p.shift;
      break;
    }
    case VDWParameters::Type::Generic:
      term = p.parameters.x * exponential(-p.parameters.y * r) - p.parameters.z / power(r, 4.0) -
             p.parameters.w / power(r, 6.0) - p.parameters2.x / power(r, 8.0) -
             p.parameters2.y / power(r, 10.0) - p.shift;
      break;
    case VDWParameters::Type::PellenqNicholson:
    {
      RadialJet br = p.parameters.y * r;
      term = p.parameters.x * exponential(-br) -
             tangToennies(6, br) * p.parameters.z / power(r, 6.0) -
             tangToennies(8, br) * p.parameters.w / power(r, 8.0) -
             tangToennies(10, br) * p.parameters2.x / power(r, 10.0) - p.shift;
      break;
    }
    case VDWParameters::Type::HydratedIonWater:
      term = p.parameters.x * exponential(-p.parameters.y * r) - p.parameters.z / power(r, 4.0) -
             p.parameters.w / power(r, 6.0) - p.parameters2.x / power(r, 12.0) - p.shift;
      break;
    case VDWParameters::Type::Mie:
      term = p.parameters.x / power(r, p.parameters.y) - p.parameters.z / power(r, p.parameters.w) - p.shift;
      break;
    default:
      return RareVDWDerivatives{0.0, 0.0, 0.0, 0.0};
  }

  return softCoreResult(term, softRadius, scaling, inverseScaling, w6);
}
}  // namespace Potentials::Detail
