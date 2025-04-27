module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <numbers>
#include <optional>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <vector>
#endif

module vdwparameters;

#ifndef USE_LEGACY_HEADERS
import <filesystem>;
import <fstream>;
import <cstdlib>;
import <iostream>;
import <sstream>;
import <vector>;
import <array>;
import <map>;
import <cmath>;
import <string>;
import <string_view>;
import <optional>;
import <numbers>;
import <algorithm>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <type_traits>;
import <iterator>;
import <functional>;
import <print>;
import <unordered_map>;
#endif

import archive;

std::string toLower(const std::string &str)
{
  std::string result = str;
  std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c) { return std::tolower(c); });
  return result;
}

VDWParameters::Type VDWParameters::stringToEnum(const std::string interactionType)
{
  static const std::unordered_map<std::string, Type> strToEnumMap = {{"none", Type::None},
                                                                     {"lennardjones", Type::LennardJones},
                                                                     {"lennard-jones", Type::LennardJones},
                                                                     {"buckingham", Type::BuckingHam},
                                                                     {"morse", Type::Morse},
                                                                     {"feynmannhibbs", Type::FeynmannHibbs},
                                                                     {"mm3", Type::MM3},
                                                                     {"bornhugginsmeyer", Type::BornHugginsMeyer}};

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

bool VDWParameters::operator==(const VDWParameters &other) const
{
  return (parameters == other.parameters && shift == other.shift &&
          tailCorrectionEnergy == other.tailCorrectionEnergy &&
          tailCorrectionPressure == other.tailCorrectionPressure && type == other.type);
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VDWParameters &p)
{
  archive << p.parameters;
  archive << p.shift;
  archive << p.tailCorrectionEnergy;
  archive << p.tailCorrectionPressure;
  archive << p.type;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VDWParameters &p)
{
  archive >> p.parameters;
  archive >> p.shift;
  archive >> p.tailCorrectionEnergy;
  archive >> p.tailCorrectionPressure;
  archive >> p.type;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("VDWParameters: Error in binary restart\n"));
  }
#endif

  return archive;
}
