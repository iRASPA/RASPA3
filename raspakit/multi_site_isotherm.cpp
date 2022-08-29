module;

module multi_site_isotherm;


import <sstream>;
import <cmath>;
import <vector>;
import <iostream>;

import special_functions;
import isotherm;
import print;

void MultiSiteIsotherm::add(const Isotherm &isotherm, const std::vector<double> &p)
{
  siteParameterIndex.push_back(parameters.size());
  parameters.insert(parameters.end(), p.begin(), p.end());
  sites.push_back(isotherm);
}

std::string MultiSiteIsotherm::print() const
{
  std::ostringstream stream;
  std::print(stream, std::print("    number of isotherm sites:  {}\n", numberOfSites));
  for(size_t i = 0; i < numberOfSites; ++i)
  {
    std::print(stream, sites[i].print(&parameters[siteParameterIndex[i]]));
  }
  return stream.str();
}

std::string MultiSiteIsotherm::printAsInputFormat() const
{
  std::ostringstream stream;
  for(size_t i = 0; i < numberOfSites; ++i)
  {
    std::print(stream, sites[i].printAsInputFormat(&parameters[siteParameterIndex[i]]));
  }
  return stream.str();
}

// returns the inverse-pressure (1/P) that corresponds to the given reduced_grand_potential psi
// advantage: for isotherms with zero equilibrium constant the result would be infinite, but the inverse is zero
double MultiSiteIsotherm::inversePressureForPsi(double reduced_grand_potential, double &cachedP0) const
{
  const double tiny = 1.0e-15;

  double left_bracket;
  double right_bracket;

  // For a single Langmuir or Langmuir-Freundlich site, the inverse can be handled analytically
  if(numberOfSites == 1)
  {
    switch(sites[0].type)
    {
      case Isotherm::Type::Langmuir:
      {
        double denominator = std::exp(reduced_grand_potential / parameters[0]) - 1.0;
        return parameters[1] / denominator;
      }
      case Isotherm::Type::Henry:
      {
        return parameters[0] / reduced_grand_potential;
      }
      case Isotherm::Type::Freundlich:
      {
        return std::pow((parameters[0] * parameters[1])/reduced_grand_potential, parameters[1]);
      }
      case Isotherm::Type::Sips:
      {
        return parameters[1] / std::pow((std::exp(reduced_grand_potential/
                (parameters[2] * parameters[0])) - 1.0), parameters[2]);
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        double denominator = std::exp(reduced_grand_potential * parameters[2] / parameters[0]) - 1.0;
        return std::pow(parameters[1] / denominator, 1.0 / parameters[2]);
      }
      default:
        break;
    }
  }

  // from here on, work with pressure, and return 1.0 / pressure at the end of the routine
  double p_start;
  if(cachedP0 <= 0.0)
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

  size_t nr_steps = 0;
  left_bracket = p_start;
  right_bracket = p_start;

  if(s < reduced_grand_potential)
  {
    // find the bracket on the right
    do 
    {
      right_bracket *= 2.0;
      s = psiForPressure(right_bracket);

      ++nr_steps;
      if(nr_steps>100000)
      {
        std::cout << "reduced_grand_potential: " << reduced_grand_potential << std::endl;
        std::cout << "psi: " << s << std::endl;
        std::cout << "p_start: " << p_start << std::endl;
        std::cout << "Left bracket: " << left_bracket << std::endl;
        std::cout << "Right bracket: " << right_bracket << std::endl;
        throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum < 1) does NOT converge\n");
      }
    }
    while(s < reduced_grand_potential);
  }
  else
  {
    // find the bracket on the left
    do 
    {
      left_bracket *= 0.5;
      s = psiForPressure(left_bracket);

      ++nr_steps;
      if(nr_steps>100000)
      {
        std::cout << "reduced_grand_potential: " << reduced_grand_potential << std::endl;
        std::cout << "psi: " << s << std::endl;
        std::cout << "p_start: " << p_start << std::endl;
        std::cout << "Left bracket: " << left_bracket << std::endl;
        std::cout << "Right bracket: " << right_bracket << std::endl;
        throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum > 1) does NOT converge\n");
      }
    }
    while(s > reduced_grand_potential);
  }

  do
  {
    double middle = 0.5 * (left_bracket + right_bracket);
    s = psiForPressure(middle);

    if(s > reduced_grand_potential)
       right_bracket = middle;
    else
       left_bracket = middle;

    ++nr_steps;
    if(nr_steps>100000)
    {
      std::cout << "Left bracket: " << left_bracket << std::endl;
      std::cout << "Right bracket: " << right_bracket << std::endl;
      throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum < 1) does NOT converge\n");
    }
  }
  while(std::abs(left_bracket - right_bracket) / std::abs(left_bracket + right_bracket) > tiny);

  double middle = 0.5 * (left_bracket + right_bracket);

  //  Store the last value of Pi0
  cachedP0 = middle;

  return 1.0 / middle;
}

double MultiSiteIsotherm::fitness() const
{
  const double penaltyCost = 50.0;
  for(size_t i = 0; i < numberOfSites; ++i)
  {
    if(sites[i].isUnphysical(&parameters[siteParameterIndex[i]])) return penaltyCost;
  }
  return 0.0;
}

std::string MultiSiteIsotherm::gnuplotFunctionString(char s) const
{
  std::ostringstream stream;
  for(size_t i = 0; i < numberOfSites; ++i)
  {
    // +1 because gnuplot start counting from 1
    stream << sites[i].gnuplotFunctionString(s, siteParameterIndex[i] + 1);
    if(i<numberOfSites-1)
    {
      stream << "+";
    }
  }
  return stream.str();
}
