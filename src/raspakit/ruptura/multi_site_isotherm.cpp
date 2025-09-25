module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <ostream>
#include <print>
#include <source_location>
#include <sstream>
#include <vector>
#endif

module multi_site_isotherm;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import special_functions;
import isotherm;
import stringutils;

std::string MultiSiteIsotherm::print() const
{
  std::ostringstream stream;
  std::print(stream, "    number of isotherm sites:  {}\n", numberOfSites);
  for (std::size_t i = 0; i < numberOfSites; ++i)
  {
    std::print(stream, "{}", sites[i].print());
  }
  return stream.str();
}

std::string MultiSiteIsotherm::printAsInputFormat() const
{
  std::ostringstream stream;
  for (std::size_t i = 0; i < numberOfSites; ++i)
  {
    std::print(stream, "{}", sites[i].printAsInputFormat());
  }
  return stream.str();
}

void MultiSiteIsotherm::add(const Isotherm &isotherm)
{
  siteParameterIndex.push_back(numberOfParameters);
  sites.push_back(isotherm);

  numberOfParameters += isotherm.numberOfParameters;
  for (std::size_t i = 0; i < isotherm.numberOfParameters; ++i)
  {
    parameterIndices.emplace_back(sites.size() - 1, i);
  }
}

// returns the inverse-pressure (1/P) that corresponds to the given reduced_grand_potential psi
// advantage: for isotherms with zero equilibrium constant the result would be infinite, but the inverse is zero
double MultiSiteIsotherm::inversePressureForPsi(double reduced_grand_potential, double &cachedP0) const
{
  const double tiny = 1.0e-15;

  double left_bracket;
  double right_bracket;

  // For a single Langmuir or Langmuir-Freundlich site, the inverse can be handled analytically
  if (numberOfSites == 1)
  {
    return sites[0].inversePressureForPsi(reduced_grand_potential, cachedP0);
  }

  // from here on, work with pressure, and return 1.0 / pressure at the end of the routine
  double p_start;
  if (cachedP0 <= 0.0)
  {
    p_start = 5.0;
  }
  else
  {
    // use the last value of Pi0
    p_start = cachedP0;
  }

  // use bisection algorithm
  double s = psiForPressure(p_start);

  std::size_t nr_steps = 0;
  left_bracket = p_start;
  right_bracket = p_start;

  if (s < reduced_grand_potential)
  {
    // find the bracket on the right
    do
    {
      right_bracket *= 2.0;
      s = psiForPressure(right_bracket);

      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::cout << "reduced_grand_potential: " << reduced_grand_potential << std::endl;
        std::cout << "psi: " << s << std::endl;
        std::cout << "p_start: " << p_start << std::endl;
        std::cout << "Left bracket: " << left_bracket << std::endl;
        std::cout << "Right bracket: " << right_bracket << std::endl;
        throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum < 1) does NOT converge\n");
      }
    } while (s < reduced_grand_potential);
  }
  else
  {
    // find the bracket on the left
    do
    {
      left_bracket *= 0.5;
      s = psiForPressure(left_bracket);

      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::cout << "reduced_grand_potential: " << reduced_grand_potential << std::endl;
        std::cout << "psi: " << s << std::endl;
        std::cout << "p_start: " << p_start << std::endl;
        std::cout << "Left bracket: " << left_bracket << std::endl;
        std::cout << "Right bracket: " << right_bracket << std::endl;
        throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum > 1) does NOT converge\n");
      }
    } while (s > reduced_grand_potential);
  }

  do
  {
    double middle = 0.5 * (left_bracket + right_bracket);
    s = psiForPressure(middle);

    if (s > reduced_grand_potential)
      right_bracket = middle;
    else
      left_bracket = middle;

    ++nr_steps;
    if (nr_steps > 100000)
    {
      std::cout << "Left bracket: " << left_bracket << std::endl;
      std::cout << "Right bracket: " << right_bracket << std::endl;
      throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum < 1) does NOT converge\n");
    }
  } while (std::abs(left_bracket - right_bracket) / std::abs(left_bracket + right_bracket) > tiny);

  double middle = 0.5 * (left_bracket + right_bracket);

  //  Store the last value of Pi0
  cachedP0 = middle;

  return 1.0 / middle;
}

double MultiSiteIsotherm::fitness() const
{
  const double penaltyCost = 50.0;
  for (std::size_t i = 0; i < numberOfSites; ++i)
  {
    if (sites[i].isUnphysical()) return penaltyCost;
  }
  return 0.0;
}

std::string MultiSiteIsotherm::gnuplotFunctionString([[maybe_unused]] char s) const
{
  std::ostringstream stream;
  for (std::size_t i = 0; i < numberOfSites; ++i)
  {
    // +1 because gnuplot start counting from 1
    stream << sites[i].gnuplotFunctionString(s, siteParameterIndex[i] + 1);
    if (i < numberOfSites - 1)
    {
      stream << "+";
    }
  }
  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MultiSiteIsotherm &c)
{
  archive << c.versionNumber;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MultiSiteIsotherm &c)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > c.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Component' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Component: Error in binary restart\n"));
  }
#endif

  return archive;
}
