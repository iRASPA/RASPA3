export module isotherm;

import <cstddef>;
import <array>;
import <vector>;
import <cmath>;
import <iostream>;
import <string>;
import <exception>;
import <numbers>;

import special_functions;
import randomnumbers;

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
    BET = 1,
    Henry = 2,
    Freundlich = 3,
    Sips = 4,
    Langmuir_Freundlich = 5,
    Redlich_Peterson = 6,
    Toth = 7,
    Unilan = 8,
    OBrien_Myers = 9,
    Quadratic = 10,
    Temkin = 11
  };

  enum class MixturePredictionMethod
  {
    IAST = 0,
    ExplicitLangmuir = 1
  };

  struct PressureRange 
  {
    double startPressure;
    double endPressure;
    size_t numberOfPoints;
  };

  Isotherm(Isotherm::Type type);

  Isotherm() noexcept = default;
  Isotherm(const Isotherm &a) noexcept = default;
  Isotherm& operator=(const Isotherm& a) noexcept = default;
  Isotherm(Isotherm&& a) noexcept = default;
  Isotherm& operator=(Isotherm&& a) noexcept = default;
  ~Isotherm() noexcept = default;

  Isotherm::Type type;

  std::string print(const double terms[]) const;
  std::string printAsInputFormat(const double terms[]) const;

  inline double value(double pressure, const double terms[]) const
  {
    switch(type)
    {
      case Isotherm::Type::Langmuir:
      {
        double temp = terms[1] * pressure;
        return terms[0] * temp / (1.0 + temp);
      }
      case Isotherm::Type::BET:
      {
        return terms[0] * terms[1] * pressure /
               ((1.0 - terms[2] * pressure) * (1.0 - terms[2] + terms[1] * pressure));
      }
      case Isotherm::Type::Henry:
      {
        return terms[0] * pressure;
      }
      case Isotherm::Type::Freundlich:
      {
        return terms[0] * std::pow(pressure, 1.0/terms[1]);
      }
      case Isotherm::Type::Sips:
      {
        double temp = std::pow(terms[1] * pressure, 1.0/terms[2]);
        return terms[0] * temp / (1.0 + temp);
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        double temp = terms[1] * std::pow(pressure, terms[2]);
        return terms[0] * temp / (1.0 + temp);
      }
      case Isotherm::Type::Redlich_Peterson:
      {
        return terms[0] * pressure / (1.0 + terms[1] * std::pow(pressure, terms[2]));
      }
      case Isotherm::Type::Toth:
      {
        double temp = terms[1] * pressure;
        return terms[0] * temp / std::pow(1.0 + std::pow(temp, terms[2]), 1.0/terms[2]);
      }
      case Isotherm::Type::Unilan:
      {
        double temp1 =  1.0 + terms[1] * std::exp(terms[2]) * pressure;
        double temp2 =  1.0 + terms[1] * std::exp(-terms[2]) * pressure;
        return terms[0] * (0.5 / terms[2]) * std::log(temp1 / temp2);
      }
      case Isotherm::Type::OBrien_Myers:
      {
        double temp1 = terms[1] * pressure;
        double temp2 = 1.0 + temp1;
        return terms[0] * (temp1 / temp2 + terms[2] * terms[2] * temp1 * (1.0 - temp1) / (temp2 * temp2 * temp2));
      }
      case Isotherm::Type::Quadratic:
      {
        double temp1 = terms[1] * pressure;
        double temp2 = terms[2] * pressure * pressure;
        return terms[0] * (temp1 + 2.0 * temp2) / (1.0 + temp1 + temp2);
      }
      case Isotherm::Type::Temkin:
      {
        double temp = terms[1] * pressure;
        double temp1 = temp / (1.0 + temp);
        return terms[0] * (temp1 + terms[2] * temp1 * temp1 * (temp1 - 1.0));
      }
      default:
        throw std::runtime_error("Error: unkown isotherm type");
    }
  }

  // the reduced grand potential psi (spreading pressure) for this pressure
  inline double psiForPressure(double pressure, const double terms[]) const
  {
    switch(type)
    {
      case Isotherm::Type::Langmuir:
      {
        return terms[0] * std::log(1.0 + terms[1] * pressure);
      }
      case Isotherm::Type::BET:
      {
        return (terms[0] * terms[1]) * std::log((1.0 - terms[2] + terms[1] * pressure) /
                                               ((1.0 - terms[2]) * (1.0 - terms[2] * pressure))) / 
               (terms[1] + terms[2] - terms[2] * terms[2]);
      }
      case Isotherm::Type::Henry:
      {
        return terms[0] * pressure;
      }
      case Isotherm::Type::Freundlich:
      {
        return terms[0] * terms[1] * std::pow(pressure, 1.0/terms[1]);
      }
      case Isotherm::Type::Sips:
      {
        return terms[2] * terms[0] * std::log(1.0 + std::pow(terms[1]*pressure, 1.0/terms[2]));
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        return terms[0] * terms[3] * std::log(1.0 + terms[1] * std::pow(pressure, terms[2]));
      }
      case Isotherm::Type::Redlich_Peterson:
      {
        if(terms[1]  * std::pow(pressure,terms[2]) < 1.0)
        {
          return terms[0] * pressure * hypergeometric2F1(1.0, 1.0/terms[2], 1.0 + 1.0/terms[2],
                       -terms[1] * std::pow(pressure,terms[2]));
        }
        else 
        {
          double prefactor = terms[0] / terms[2];
          double temp = std::numbers::pi / (std::pow(terms[1], 1.0/terms[2]) * std::sin(std::numbers::pi * 1.0/terms[2]));

          double term1 = -1.0/(terms[1] * std::pow(pressure, terms[2]));
          double numerator = 1.0;
          double sum=0.0;
          // quickly converging sum
          for(size_t k = 1; k <= 15; k++)
          {
            numerator *= term1;
            sum += numerator / (static_cast<double>(k) * terms[2] - 1.0);
          }
          return prefactor * (temp + pressure * terms[2] * sum);
        }
      }
      case Isotherm::Type::Toth:
      {
        double temp = terms[1] * pressure;
        double theta = temp / std::pow(1.0 + std::pow(temp, terms[2]), 1.0/terms[2]);
        double theta_pow = std::pow(theta, terms[2]);
        double psi = terms[0] * (theta - (theta / terms[2]) * std::log(1.0-theta_pow));

        // use the first 100 terms of the sum
        double temp1 = terms[0] * theta;
        double temp2 = 0.0;
        for(size_t k = 1; k <= 100; ++k)
        {
          temp1 *= theta_pow;
          temp2 += terms[2];
          psi -= temp1 / (temp2 * (temp2 + 1.0));
        }

        return psi;
      }
      case Isotherm::Type::Unilan:
      {
        return (0.5 * terms[0] / terms[2]) * (li2(-terms[1] * std::exp(-terms[2]) * pressure) - 
                                              li2(-terms[1] * std::exp(terms[2]) * pressure));
      }
      case Isotherm::Type::OBrien_Myers:
      {
        double temp1 = terms[1] * pressure;
        double temp2 = 1.0 + temp1;
        return terms[0] * (std::log(temp2) + 0.5 * terms[2] * terms[2] * temp1 / (temp2 * temp2));
      }
      case Isotherm::Type::Quadratic:
      {
        double temp1 = terms[1] * pressure;
        double temp2 = terms[2] * pressure * pressure;
        return terms[0] * std::log(1.0 + temp1 + temp2);
      }
      case Isotherm::Type::Temkin:
      {
        double temp = terms[1] * pressure;
        double temp1 = temp / (1.0 + temp);
        return terms[0] * (std::log(1.0 + temp) - 0.5 * terms[2] * temp1 * temp1);
      }
      default:
        throw std::runtime_error("Error: unkown isotherm type");
    }
  }

  void randomize(double params[], double maximumLoading);

  bool isUnphysical(const double terms[]) const;

  std::string gnuplotFunctionString(char s, size_t i) const;
};

