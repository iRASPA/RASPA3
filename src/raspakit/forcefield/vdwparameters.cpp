module;

module vdwparameters;

import std;

import archive;
import stringutils;
import units;
import double4;

VDWParameters::Type VDWParameters::stringToEnum(const std::string interactionType)
{
  static const std::unordered_map<std::string, Type> strToEnumMap = {
      {"none", Type::None},
      {"lennardjones", Type::LennardJones},
      {"lennard-jones", Type::LennardJones},
      {"buckingham", Type::BuckingHam},
      {"morse", Type::Morse},
      {"feynmannhibbs", Type::FeynmannHibbs},
      {"feynman-hibbs-lennard-jones", Type::FeynmannHibbs},
      {"mm3", Type::MM3},
      {"mm3-vdw", Type::MM3},
      {"bornhugginsmeyer", Type::BornHugginsMeyer},
      {"born-huggins-meyer", Type::BornHugginsMeyer},
      {"lennard-jones-shifted-force", Type::LennardJonesShiftedForce},
      {"lennardjonesshiftedforce", Type::LennardJonesShiftedForce},
      {"12-6", Type::Potential12_6},
      {"potential-12-6", Type::Potential12_6},
      {"12-6-2-0", Type::Potential12_6_2_0},
      {"potential-12-6-2-0", Type::Potential12_6_2_0},
      {"cff-9-6", Type::CFF9_6},
      {"cff9-6", Type::CFF9_6},
      {"cff-eps-sigma", Type::CFFEpsilonSigma},
      {"cffepssigma", Type::CFFEpsilonSigma},
      {"matsuoka-clementi-yoshimine", Type::MatsuokaClementiYoshimine},
      {"mcy", Type::MatsuokaClementiYoshimine},
      {"generic", Type::Generic},
      {"pellenq-nicholson", Type::PellenqNicholson},
      {"pellenqnicholson", Type::PellenqNicholson},
      {"hydrated-ion-water", Type::HydratedIonWater},
      {"hydratedionwater", Type::HydratedIonWater},
      {"mie", Type::Mie},
      {"wca", Type::WeeksChandlerAndersen},
      {"weeks-chandler-andersen", Type::WeeksChandlerAndersen},
      {"repulsiveharmonic", Type::RepulsiveHarmonic},
      {"repulsive-harmonic", Type::RepulsiveHarmonic}};

  std::string lowerInteractionType = toLower(interactionType);
  auto it = strToEnumMap.find(lowerInteractionType);
  if (it != strToEnumMap.end())
  {
    return it->second;
  }
  else
  {
    throw std::invalid_argument("Invalid string for interaction type conversion");
  }
}

std::string VDWParameters::nameOfType(VDWParameters::Type type)
{
  switch (type)
  {
    case Type::None:
      return "None";
    case Type::LennardJones:
      return "Lennard-Jones";
    case Type::BuckingHam:
      return "Buckingham";
    case Type::Morse:
      return "Morse";
    case Type::FeynmannHibbs:
      return "Feynman-Hibbs Lennard-Jones";
    case Type::MM3:
      return "MM3";
    case Type::BornHugginsMeyer:
      return "Born-Huggins-Meyer";
    case Type::LennardJonesShiftedForce:
      return "Lennard-Jones shifted-force";
    case Type::Potential12_6:
      return "12-6";
    case Type::Potential12_6_2_0:
      return "12-6-2-0";
    case Type::CFF9_6:
      return "CFF 9-6";
    case Type::CFFEpsilonSigma:
      return "CFF eps-sigma";
    case Type::MatsuokaClementiYoshimine:
      return "Matsuoka-Clementi-Yoshimine";
    case Type::Generic:
      return "Generic";
    case Type::PellenqNicholson:
      return "Pellenq-Nicholson";
    case Type::HydratedIonWater:
      return "Hydrated-ion-water";
    case Type::Mie:
      return "Mie";
    case Type::WeeksChandlerAndersen:
      return "Weeks-Chandler-Andersen";
    case Type::RepulsiveHarmonic:
      return "Repulsive harmonic";
  }
  return "Unknown";
}

VDWParameters::ParameterMetadata VDWParameters::parameterMetadata(VDWParameters::Type type)
{
  switch (type)
  {
    case Type::None:
      return {0, {}};
    case Type::LennardJones:
      // p0: epsilon [K], p1: sigma [A]
      return {2, {true, false}};
    case Type::BuckingHam:
      // p0: A [K], p1: b [A^-1], p2: C6 [K A^6]
      return {3, {true, false, true}};
    case Type::Morse:
      // p0: D [K], p1: a [A^-1], p2: r0 [A]
      return {3, {true, false, false}};
    case Type::FeynmannHibbs:
      // p0: epsilon [K], p1: sigma [A], p2: reduced mass [amu]
      return {3, {true, false, false}};
    case Type::MM3:
      // p0: epsilon [K], p1: sigma [A]
      return {2, {true, false}};
    case Type::BornHugginsMeyer:
      // p0: A [K], p1: b [A^-1], p2: sigma [A], p3: C6 [K A^6], p4: C8 [K A^8]
      return {5, {true, false, false, true, true}};
    case Type::LennardJonesShiftedForce:
      // p0: epsilon [K], p1: sigma [A]
      return {2, {true, false}};
    case Type::Potential12_6:
      // p0: [K A^12], p1: [K A^6]
      return {2, {true, true}};
    case Type::Potential12_6_2_0:
      // p0: [K A^12], p1: [K A^6], p2: [K A^2], p3: [K]
      return {4, {true, true, true, true}};
    case Type::CFF9_6:
      // p0: [K A^9], p1: [K A^6]
      return {2, {true, true}};
    case Type::CFFEpsilonSigma:
      // p0: epsilon [K], p1: sigma [A]
      return {2, {true, false}};
    case Type::MatsuokaClementiYoshimine:
      // p0: A [K], p1: a [A^-1], p2: B [K], p3: b [A^-1]
      return {4, {true, false, true, false}};
    case Type::Generic:
      // p0: A [K], p1: b [A^-1], p2: C4 [K A^4], p3: C6 [K A^6], p4: C8 [K A^8], p5: C10 [K A^10]
      return {6, {true, false, true, true, true, true}};
    case Type::PellenqNicholson:
      // p0: A [K], p1: b [A^-1], p2: C6 [K A^6], p3: C8 [K A^8], p4: C10 [K A^10]
      return {5, {true, false, true, true, true}};
    case Type::HydratedIonWater:
      // p0: A [K], p1: b [A^-1], p2: C4 [K A^4], p3: C6 [K A^6], p4: C12 [K A^12]
      return {5, {true, false, true, true, true}};
    case Type::Mie:
      // p0: [K A^p1], p1: repulsive exponent [-], p2: [K A^p3], p3: attractive exponent [-]
      return {4, {true, false, true, false}};
    case Type::WeeksChandlerAndersen:
      // p0: epsilon [K], p1: sigma [A]
      return {2, {true, false}};
    case Type::RepulsiveHarmonic:
      // p0: k [K], p1: r0 [A]
      return {2, {true, false}};
  }
  return {0, {}};
}

VDWParameters::VDWParameters(VDWParameters::Type type, const std::vector<double> &values)
    : parameters(0.0, 0.0, 0.0, 0.0),
      shift(0.0),
      type(type),
      parameters2(0.0, 0.0, 0.0, 0.0),
      tailCorrectionEnergy(0.0),
      tailCorrectionPressure(0.0)
{
  ParameterMetadata metadata = parameterMetadata(type);
  if (values.size() < metadata.count)
  {
    throw std::runtime_error(std::format("Potential '{}' requires {} parameters, but only {} were given\n",
                                         nameOfType(type), metadata.count, values.size()));
  }

  for (std::size_t i = 0; i < std::min(values.size(), std::size_t{6}); ++i)
  {
    double value = values[i];
    if (i < metadata.count && metadata.isEnergy[i])
    {
      value *= Units::KelvinToEnergy;
    }
    if (i < 4)
    {
      parameters[i] = value;
    }
    else
    {
      parameters2[i - 4] = value;
    }
  }
}

double VDWParameters::potentialEnergyAtFullCoupling(double rr) const
{
  switch (type)
  {
    case Type::LennardJones:
    {
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      double arg1 = parameters.x;
      double arg2 = parameters.y * parameters.y;
      double temp = (rr / arg2);
      double rri3 = 1.0 / (temp * temp * temp);
      return 4.0 * arg1 * (rri3 * (rri3 - 1.0));
    }
    case Type::BuckingHam:
    {
      // p_0*exp(-p_1*r)-p_2/r^6
      double r = std::sqrt(rr);
      double rri6 = 1.0 / (rr * rr * rr);
      return parameters.x * std::exp(-parameters.y * r) - parameters.z * rri6;
    }
    case Type::Morse:
    {
      // p_0*{(1.0-exp[-p_1*(r-p_2)])^2-1.0}
      double r = std::sqrt(rr);
      double expTerm = std::exp(-parameters.y * (r - parameters.z));
      return parameters.x * ((1.0 - expTerm) * (1.0 - expTerm) - 1.0);
    }
    case Type::FeynmannHibbs:
    {
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 mu K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      double arg1 = parameters.x;
      double arg2 = parameters.y * parameters.y;
      double fh = parameters2.x;
      double temp = (rr / arg2);
      double rri3 = 1.0 / (temp * temp * temp);
      return 4.0 * arg1 * (rri3 * (rri3 - 1.0) + fh * (rri3 * (132.0 * rri3 - 30.0) / rr));
    }
    case Type::MM3:
    {
      // eps*[1.84e5*exp(-12/P)-2.25*P^6]  if P<=3.02   with P=sigma/r
      // eps*192.27*P^2                    if P>3.02
      double r = std::sqrt(rr);
      double P = parameters.y / r;
      if (P > 3.02)
      {
        return parameters.x * 192.270 * P * P;
      }
      double P2 = P * P;
      double P6 = P2 * P2 * P2;
      return parameters.x * (1.84e5 * std::exp(-12.0 / P) - 2.25 * P6);
    }
    case Type::BornHugginsMeyer:
    {
      // p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8
      double r = std::sqrt(rr);
      double rri6 = 1.0 / (rr * rr * rr);
      double rri8 = rri6 / rr;
      return parameters.x * std::exp(parameters.y * (parameters.z - r)) - parameters.w * rri6 - parameters2.x * rri8;
    }
    case Type::LennardJonesShiftedForce:
    {
      // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]-[(p_1/rc)^12-(p_1/rc)^6]+[12*(p_1/rc)^12-6*(p_1/rc)^6]*(rc-r)/rc}
      double r = std::sqrt(rr);
      double arg1 = parameters.x;
      double arg2 = parameters.y * parameters.y;
      double temp = (rr / arg2);
      double rri3 = 1.0 / (temp * temp * temp);
      double rri3_c = parameters2.x;
      double cutOff = parameters2.y;
      return 4.0 * arg1 *
             (rri3 * (rri3 - 1.0) - rri3_c * (rri3_c - 1.0) +
              (12.0 * rri3_c * rri3_c - 6.0 * rri3_c) * (r - cutOff) / cutOff);
    }
    case Type::Potential12_6:
    {
      // p_0/r^12-p_1/r^6
      double rri6 = 1.0 / (rr * rr * rr);
      return parameters.x * rri6 * rri6 - parameters.y * rri6;
    }
    case Type::Potential12_6_2_0:
    {
      // p_0/r^12+p_1/r^6+p_2/r^2+p_3
      double rri6 = 1.0 / (rr * rr * rr);
      return parameters.x * rri6 * rri6 + parameters.y * rri6 + parameters.z / rr + parameters.w;
    }
    case Type::CFF9_6:
    {
      // p_0/r^9-p_1/r^6
      double rri6 = 1.0 / (rr * rr * rr);
      double rri9 = rri6 * std::sqrt(rri6);
      return parameters.x * rri9 - parameters.y * rri6;
    }
    case Type::CFFEpsilonSigma:
    {
      // p_0*[2*(p_1/r)^9-3*(p_1/r)^6]
      double arg2 = parameters.y * parameters.y;
      double temp = arg2 / rr;
      double si6 = temp * temp * temp;
      double si9 = si6 * std::sqrt(si6);
      return parameters.x * (2.0 * si9 - 3.0 * si6);
    }
    case Type::MatsuokaClementiYoshimine:
    {
      // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
      double r = std::sqrt(rr);
      return parameters.x * std::exp(-parameters.y * r) + parameters.z * std::exp(-parameters.w * r);
    }
    case Type::Generic:
    {
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10
      double r = std::sqrt(rr);
      double rri2 = 1.0 / rr;
      double rri4 = rri2 * rri2;
      double rri6 = rri4 * rri2;
      double rri8 = rri6 * rri2;
      double rri10 = rri8 * rri2;
      return parameters.x * std::exp(-parameters.y * r) - parameters.z * rri4 - parameters.w * rri6 -
             parameters2.x * rri8 - parameters2.y * rri10;
    }
    case Type::PellenqNicholson:
    {
      // p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10
      double r = std::sqrt(rr);
      double f6, f8, f10;
      computeTangToenniesDampingCoefficients(r, parameters.y, f6, f8, f10);
      double rri2 = 1.0 / rr;
      double rri6 = rri2 * rri2 * rri2;
      double rri8 = rri6 * rri2;
      double rri10 = rri8 * rri2;
      return parameters.x * std::exp(-parameters.y * r) - f6 * parameters.z * rri6 - f8 * parameters.w * rri8 -
             f10 * parameters2.x * rri10;
    }
    case Type::HydratedIonWater:
    {
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12
      double r = std::sqrt(rr);
      double rri2 = 1.0 / rr;
      double rri4 = rri2 * rri2;
      double rri6 = rri4 * rri2;
      double rri12 = rri6 * rri6;
      return parameters.x * std::exp(-parameters.y * r) - parameters.z * rri4 - parameters.w * rri6 -
             parameters2.x * rri12;
    }
    case Type::Mie:
    {
      // p_0/r^p_1-p_2/r^p_3
      double r = std::sqrt(rr);
      return parameters.x / std::pow(r, parameters.y) - parameters.z / std::pow(r, parameters.w);
    }
    case Type::WeeksChandlerAndersen:
    {
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+p_0 if r<2^(1/6)*p_1, otherwise 0
      double arg2 = parameters.y * parameters.y;
      if (rr > std::cbrt(2.0) * arg2)
      {
        return 0.0;
      }
      double temp = (rr / arg2);
      double rri3 = 1.0 / (temp * temp * temp);
      return 4.0 * parameters.x * (rri3 * (rri3 - 1.0)) + parameters.x;
    }
    case Type::RepulsiveHarmonic:
    {
      double r = std::sqrt(rr);
      double temp = (1.0 - r / parameters.y);
      return r >= parameters.y ? 0.0 : 0.5 * parameters.x * temp * temp;
    }
    case Type::None:
      return 0.0;
  }
  return 0.0;
}

void VDWParameters::computeDerivedParameters(double cutOff, double temperature)
{
  // type-specific derived constants (must be set before evaluating the potential)
  switch (type)
  {
    case Type::FeynmannHibbs:
    {
      // hbar^2/(24 mu kB T), with mu the reduced mass in atomic mass units
      parameters2.x =
          parameters.z > 0.0 ? Units::FeymannHibbsConversionFactor / (parameters.z * temperature) : 0.0;
      break;
    }
    case Type::LennardJonesShiftedForce:
    {
      double temp = (parameters.y * parameters.y) / (cutOff * cutOff);
      parameters2.x = temp * temp * temp;  // (sigma/rc)^6
      parameters2.y = cutOff;
      break;
    }
    default:
      break;
  }

  // soft-core reference diameter w^6 for the continuous-fractional lambda-scaling
  switch (type)
  {
    case Type::None:
    case Type::RepulsiveHarmonic:
    case Type::Morse:
    case Type::MatsuokaClementiYoshimine:
      // finite at r -> 0: no soft-core needed
      parameters2.w = 0.0;
      break;
    case Type::LennardJones:
    case Type::FeynmannHibbs:
    case Type::MM3:
    case Type::LennardJonesShiftedForce:
    case Type::CFFEpsilonSigma:
    case Type::WeeksChandlerAndersen:
    {
      double sigma2 = parameters.y * parameters.y;
      parameters2.w = sigma2 * sigma2 * sigma2;
      break;
    }
    case Type::Potential12_6:
    {
      parameters2.w = (parameters.x > 0.0 && parameters.y > 0.0) ? parameters.x / parameters.y : 1.0;
      break;
    }
    case Type::CFF9_6:
    {
      double w3 = (parameters.x > 0.0 && parameters.y > 0.0) ? parameters.x / parameters.y : 1.0;
      parameters2.w = w3 * w3;
      break;
    }
    case Type::Mie:
    {
      if (parameters.x > 0.0 && parameters.z > 0.0 && parameters.y > parameters.w)
      {
        parameters2.w = std::pow(parameters.x / parameters.z, 6.0 / (parameters.y - parameters.w));
      }
      else
      {
        parameters2.w = 1.0;
      }
      break;
    }
    case Type::BuckingHam:
    case Type::BornHugginsMeyer:
    case Type::Generic:
    case Type::PellenqNicholson:
    case Type::HydratedIonWater:
    case Type::Potential12_6_2_0:
    {
      // locate the outermost zero-crossing of the repulsive wall numerically
      double wall = 0.0;
      double rOuter = cutOff;
      double energyOuter = potentialEnergyAtFullCoupling(rOuter * rOuter);
      const double step = 0.01;
      for (double r = rOuter - step; r > 0.2; r -= step)
      {
        double energy = potentialEnergyAtFullCoupling(r * r);
        if (energyOuter <= 0.0 && energy > 0.0)
        {
          // bisect between r (positive) and r + step (non-positive)
          double lo = r, hi = r + step;
          for (std::size_t iteration = 0; iteration < 50; ++iteration)
          {
            double mid = 0.5 * (lo + hi);
            if (potentialEnergyAtFullCoupling(mid * mid) > 0.0)
            {
              lo = mid;
            }
            else
            {
              hi = mid;
            }
          }
          wall = 0.5 * (lo + hi);
          break;
        }
        energyOuter = energy;
        rOuter = r;
      }
      if (wall <= 0.0)
      {
        wall = 1.0;
      }
      double wall2 = wall * wall;
      parameters2.w = wall2 * wall2 * wall2;
      break;
    }
  }
}

bool VDWParameters::operator==(const VDWParameters &other) const
{
  return (parameters == other.parameters && parameters2 == other.parameters2 && shift == other.shift &&
          tailCorrectionEnergy == other.tailCorrectionEnergy &&
          tailCorrectionPressure == other.tailCorrectionPressure && type == other.type);
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VDWParameters &p)
{
  archive << p.parameters;
  archive << p.parameters2;
  archive << p.shift;
  archive << p.tailCorrectionEnergy;
  archive << p.tailCorrectionPressure;
  archive << p.type;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VDWParameters &p)
{
  archive >> p.parameters;
  archive >> p.parameters2;
  archive >> p.shift;
  archive >> p.tailCorrectionEnergy;
  archive >> p.tailCorrectionPressure;
  archive >> p.type;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("VDWParameters: Error in binary restart\n"));
  }
#endif

  return archive;
}
