module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>
#endif

export module isotherm;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import special_functions;
import randomnumbers;

constexpr std::size_t maxTerms = 5;

// Langmuir:
// parameter 0: K
// parameter 1: N
//
// Langmuir-Freundlich
// parameter 0: K
// parameter 1: N
// parameter 2: power

export struct Isotherm
{
  enum class Type
  {
    Langmuir = 0,
    Anti_Langmuir = 1,
    BET = 2,
    Henry = 3,
    Freundlich = 4,
    Sips = 5,
    Langmuir_Freundlich = 6,
    Redlich_Peterson = 7,
    Toth = 8,
    Unilan = 9,
    OBrien_Myers = 10,
    Quadratic = 11,
    Temkin = 12
  };

  Isotherm(Isotherm::Type type, const std::vector<double> &values, std::size_t numberOfValues);

  bool operator==(Isotherm const &) const = default;

  Isotherm::Type type;
  std::vector<double> parameters;
  std::size_t numberOfParameters;

  std::string print() const;
  std::string printAsInputFormat() const;

  inline double value(double pressure) const
  {
    switch (type)
    {
      case Isotherm::Type::Langmuir:
      {
        double temp = parameters[1] * pressure;
        return parameters[0] * temp / (1.0 + temp);
      }
      case Isotherm::Type::Anti_Langmuir:
      {
        return parameters[0] * pressure / (1.0 - parameters[1] * pressure);
      }
      case Isotherm::Type::BET:
      {
        return parameters[0] * parameters[1] * pressure /
               ((1.0 - parameters[2] * pressure) * (1.0 - parameters[2] + parameters[1] * pressure));
      }
      case Isotherm::Type::Henry:
      {
        return parameters[0] * pressure;
      }
      case Isotherm::Type::Freundlich:
      {
        return parameters[0] * std::pow(pressure, 1.0 / parameters[1]);
      }
      case Isotherm::Type::Sips:
      {
        double temp = std::pow(parameters[1] * pressure, 1.0 / parameters[2]);
        return parameters[0] * temp / (1.0 + temp);
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        double temp = parameters[1] * std::pow(pressure, parameters[2]);
        return parameters[0] * temp / (1.0 + temp);
      }
      case Isotherm::Type::Redlich_Peterson:
      {
        return parameters[0] * pressure / (1.0 + parameters[1] * std::pow(pressure, parameters[2]));
      }
      case Isotherm::Type::Toth:
      {
        double temp = parameters[1] * pressure;
        return parameters[0] * temp / std::pow(1.0 + std::pow(temp, parameters[2]), 1.0 / parameters[2]);
      }
      case Isotherm::Type::Unilan:
      {
        double temp1 = 1.0 + parameters[1] * std::exp(parameters[2]) * pressure;
        double temp2 = 1.0 + parameters[1] * std::exp(-parameters[2]) * pressure;
        return parameters[0] * (0.5 / parameters[2]) * std::log(temp1 / temp2);
      }
      case Isotherm::Type::OBrien_Myers:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = 1.0 + temp1;
        return parameters[0] *
               (temp1 / temp2 + parameters[2] * parameters[2] * temp1 * (1.0 - temp1) / (temp2 * temp2 * temp2));
      }
      case Isotherm::Type::Quadratic:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = parameters[2] * pressure * pressure;
        return parameters[0] * (temp1 + 2.0 * temp2) / (1.0 + temp1 + temp2);
      }
      case Isotherm::Type::Temkin:
      {
        double temp = parameters[1] * pressure;
        double temp1 = temp / (1.0 + temp);
        return parameters[0] * (temp1 + parameters[2] * temp1 * temp1 * (temp1 - 1.0));
      }
      default:
        throw std::runtime_error("Error: unkown isotherm type");
    }
  }

  // the reduced grand potential psi (spreading pressure) for this pressure
  inline double psiForPressure(double pressure) const
  {
    switch (type)
    {
      case Isotherm::Type::Langmuir:
      {
        return parameters[0] * std::log(1.0 + parameters[1] * pressure);
      }
      case Isotherm::Type::Anti_Langmuir:
      {
        return -(parameters[0] / parameters[1]) * std::log(1.0 - parameters[1] * pressure);
      }
      case Isotherm::Type::BET:
      {
        return (parameters[0] * parameters[1]) *
               std::log((1.0 - parameters[2] + parameters[1] * pressure) /
                        ((1.0 - parameters[2]) * (1.0 - parameters[2] * pressure))) /
               (parameters[1] + parameters[2] - parameters[2] * parameters[2]);
      }
      case Isotherm::Type::Henry:
      {
        return parameters[0] * pressure;
      }
      case Isotherm::Type::Freundlich:
      {
        return parameters[0] * parameters[1] * std::pow(pressure, 1.0 / parameters[1]);
      }
      case Isotherm::Type::Sips:
      {
        return parameters[2] * parameters[0] * std::log(1.0 + std::pow(parameters[1] * pressure, 1.0 / parameters[2]));
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        return (parameters[0] / parameters[2]) * std::log(1.0 + parameters[1] * std::pow(pressure, parameters[2]));
      }
      case Isotherm::Type::Redlich_Peterson:
      {
        if (parameters[1] * std::pow(pressure, parameters[2]) < 1.0)
        {
          return parameters[0] * pressure *
                 hypergeometric2F1(1.0, 1.0 / parameters[2], 1.0 + 1.0 / parameters[2],
                                   -parameters[1] * std::pow(pressure, parameters[2]));
        }
        else
        {
          double prefactor = parameters[0] / parameters[2];
          double temp = std::numbers::pi / (std::pow(parameters[1], 1.0 / parameters[2]) *
                                            std::sin(std::numbers::pi * 1.0 / parameters[2]));

          double term1 = -1.0 / (parameters[1] * std::pow(pressure, parameters[2]));
          double numerator = 1.0;
          double sum = 0.0;
          // quickly converging sum
          for (std::size_t k = 1; k <= 15; k++)
          {
            numerator *= term1;
            sum += numerator / (static_cast<double>(k) * parameters[2] - 1.0);
          }
          return prefactor * (temp + pressure * parameters[2] * sum);
        }
      }
      case Isotherm::Type::Toth:
      {
        double temp = parameters[1] * pressure;
        double theta = temp / std::pow(1.0 + std::pow(temp, parameters[2]), 1.0 / parameters[2]);
        double theta_pow = std::pow(theta, parameters[2]);
        double psi = parameters[0] * (theta - (theta / parameters[2]) * std::log(1.0 - theta_pow));

        // use the first 100 terms of the sum
        double temp1 = parameters[0] * theta;
        double temp2 = 0.0;
        for (std::size_t k = 1; k <= 100; ++k)
        {
          temp1 *= theta_pow;
          temp2 += parameters[2];
          psi -= temp1 / (temp2 * (temp2 + 1.0));
        }

        return psi;
      }
      case Isotherm::Type::Unilan:
      {
        return (0.5 * parameters[0] / parameters[2]) * (li2(-parameters[1] * std::exp(-parameters[2]) * pressure) -
                                                        li2(-parameters[1] * std::exp(parameters[2]) * pressure));
      }
      case Isotherm::Type::OBrien_Myers:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = 1.0 + temp1;
        return parameters[0] * (std::log(temp2) + 0.5 * parameters[2] * parameters[2] * temp1 / (temp2 * temp2));
      }
      case Isotherm::Type::Quadratic:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = parameters[2] * pressure * pressure;
        return parameters[0] * std::log(1.0 + temp1 + temp2);
      }
      case Isotherm::Type::Temkin:
      {
        double temp = parameters[1] * pressure;
        double temp1 = temp / (1.0 + temp);
        return parameters[0] * (std::log(1.0 + temp) - 0.5 * parameters[2] * temp1 * temp1);
      }
      default:
        throw std::runtime_error("Error: unkown isotherm type");
    }
  }

  inline double inversePressureForPsi(double reduced_grand_potential, double &cachedP0) const
  {
    switch (type)
    {
      case Isotherm::Type::Langmuir:
      {
        double denominator = std::exp(reduced_grand_potential / parameters[0]) - 1.0;
        return parameters[1] / denominator;
      }
      case Isotherm::Type::Anti_Langmuir:
      {
        double denominator = 1.0 - std::exp(-parameters[1] * reduced_grand_potential / parameters[0]);
        return parameters[1] / denominator;
      }
      case Isotherm::Type::Henry:
      {
        return parameters[0] / reduced_grand_potential;
      }
      case Isotherm::Type::Freundlich:
      {
        return std::pow((parameters[0] * parameters[1]) / reduced_grand_potential, parameters[1]);
      }
      case Isotherm::Type::Sips:
      {
        return parameters[1] /
               std::pow((std::exp(reduced_grand_potential / (parameters[2] * parameters[0])) - 1.0), parameters[2]);
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        double denominator = std::exp(reduced_grand_potential * parameters[2] / parameters[0]) - 1.0;
        return std::pow(parameters[1] / denominator, 1.0 / parameters[2]);
      }
      default:
      {
        const double tiny = 1.0e-15;

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
        double left_bracket = p_start;
        double right_bracket = p_start;

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
              throw std::runtime_error(
                  "Error (Inverse bisection): initial bracketing (for sum < 1) "
                  "does NOT converge\n");
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
              throw std::runtime_error(
                  "Error (Inverse bisection): initial bracketing (for sum > 1) does "
                  "NOT converge\n");
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
            throw std::runtime_error(
                "Error (Inverse bisection): initial bracketing (for sum < 1) does "
                "NOT converge\n");
          }
        } while (std::abs(left_bracket - right_bracket) / std::abs(left_bracket + right_bracket) > tiny);

        double middle = 0.5 * (left_bracket + right_bracket);

        //  Store the last value of Pi0
        cachedP0 = middle;

        return 1.0 / middle;
      }
    }
  }

  void randomize(RandomNumber &random, double maximumLoading);

  bool isUnphysical() const;

  std::string gnuplotFunctionString(char s, std::size_t i) const;
};
