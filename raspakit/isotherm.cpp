module;

module isotherm;

import <cstdlib>;
import <iostream>;
import <sstream>;
import <cmath>;

import randomnumbers;
import print;

Isotherm::Isotherm(Isotherm::Type type):
    type(type)
{
}

std::string Isotherm::print(const double terms[]) const
{
  std::ostringstream stream;

  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      std::print(stream, "    Langmuir isotherm\n");
      std::print(stream, "        q_sat: {}\n", terms[0]);
      std::print(stream, "        b:     {}\n", terms[1]);
      break;
    }
    case Isotherm::Type::BET:
    {
      std::print(stream, "    BET isotherm\n");
      std::print(stream, "        q_sat: {}\n", terms[0]);
      std::print(stream, "        b:     {}\n", terms[1]);
      std::print(stream, "        c:     {}\n", terms[2]);
      break;
    }
    case Isotherm::Type::Henry:
    {
      std::print(stream, "    Henry isotherm\n");
      std::print(stream, "        a:     {}\n", terms[0]);
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      std::print(stream, "    Freundlich isotherm\n");
      std::print(stream, "        a:     {}\n", terms[0]);
      std::print(stream, "        nu:    {}\n", terms[1]);
      break;
    }
    case Isotherm::Type::Sips:
    {
      std::print(stream, "    Sips isotherm\n");
      std::print(stream, "        q_sat: {}\n", terms[0]);
      std::print(stream, "        b:     {}\n", terms[1]);
      std::print(stream, "        nu:    {}\n", terms[2]);
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      std::print(stream, "    Langmuir-Freundlich isotherm\n");
      std::print(stream, "        q_sat: {}\n", terms[0]);
      std::print(stream, "        b:     {}\n", terms[1]);
      std::print(stream, "        nu:    {}\n", terms[2]);
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      std::print(stream, "    Redlich-Peterson isotherm\n");
      std::print(stream, "        a:     {}\n", terms[0]);
      std::print(stream, "        b:     {}\n", terms[1]);
      std::print(stream, "        nu:    {}\n", terms[2]);
      break;
    }
    case Isotherm::Type::Toth:
    {
      std::print(stream, "    Toth isotherm\n");
      std::print(stream, "        q_sat: {}\n", terms[0]);
      std::print(stream, "        b:     {}\n", terms[1]);
      std::print(stream, "        nu:    {}\n", terms[2]);
      break;
    }
    case Isotherm::Type::Unilan:
    {
      std::print(stream, "    Unilan isotherm\n");
      std::print(stream, "        q_sat: {}\n", terms[0]);
      std::print(stream, "        b:     {}\n", terms[1]);
      std::print(stream, "        eta:   {}\n", terms[2]);
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      std::print(stream, "    O'Brian & Myers isotherm\n");
      std::print(stream, "        q_sat: {}\n", terms[0]);
      std::print(stream, "        b:     {}\n", terms[1]);
      std::print(stream, "        sigma: {}\n", terms[2]);
      break;
    }
    case Isotherm::Type::Quadratic:
    {
      std::print(stream, "    Quadratic isotherm\n");
      std::print(stream, "        q_sat: {}\n", terms[0]);
      std::print(stream, "        b:     {}\n", terms[1]);
      std::print(stream, "        c:     {}\n", terms[2]);
      break;
    }
    case Isotherm::Type::Temkin:
    {
      std::print(stream, "    Temkin isotherm\n");
      std::print(stream, "        q_sat: {}\n", terms[0]);
      std::print(stream, "        b:     {}\n", terms[1]);
      std::print(stream, "        c:     {}\n", terms[2]);
      break;
    }
    default:
      break;
  }
  return stream.str();
}

std::string Isotherm::printAsInputFormat(const double terms[]) const
{
  std::ostringstream stream;

  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      std::print(stream, "                Langmuir              {:<11g}  {:<11g}\n", terms[0], terms[1]);
      break;
    }
    case Isotherm::Type::BET:
    {
      std::print(stream, "                BET                   {:<11g}  {:<11g}  {:<11g}\n", terms[0], terms[1], terms[2]);
      break;
    }
    case Isotherm::Type::Henry:
    {
      std::print(stream, "                Henry                 {:<11g}\n", terms[0]);
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      std::print(stream, "                Freundlich            {:<11g}  {:<11g}\n", terms[0], terms[1]);
      break;
    }
    case Isotherm::Type::Sips:
    {
      std::print(stream, "                Sips                  {:<11g}  {:<11g}  {:<11g}\n", terms[0], terms[1], terms[2]);
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      std::print(stream, "                Langmuir-Freundlich   {:<11g}  {:<11g}  {:<11g}\n", terms[0], terms[1], terms[2]);
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      std::print(stream, "                Redlich-Peterson      {:<11g}  {:<11g}  {:<11g}\n", terms[0], terms[1], terms[2]);
      break;
    }
    case Isotherm::Type::Toth:
    {
      std::print(stream, "                Toth                  {:<11g}  {:<11g}  {:<11g}\n", terms[0], terms[1], terms[2]);
      break;
    }
    case Isotherm::Type::Unilan:
    {
      std::print(stream, "                Unilan                {:<11g}  {:<11g}  {:<11g}\n", terms[0], terms[1], terms[2]);
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      std::print(stream, "                O'Brien&Myers         {:<11g}  {:<11g}  {:<11g}\n", terms[0], terms[1], terms[2]);
      break;
    }
    case Isotherm::Type::Quadratic:
    {
      std::print(stream, "                Quadratic             {:<11g}  {:<11g}  {:<11g}\n", terms[0], terms[1], terms[2]);
      break;
    }
    case Isotherm::Type::Temkin:
    {
      std::print(stream, "                Temkin                {:<11g}  {:<11g}  {:<11g}\n", terms[0], terms[1], terms[2]);
      break;
    }
    default:
      break;
  }
  return stream.str();
}

bool Isotherm::isUnphysical(const double terms[]) const
{
  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      if(terms[0] < 0 || terms[0] > 1.0e20 || terms[1] < 0.0 || terms[1] > 1.0e10 ) return true;
      return false;
    }
    case Isotherm::Type::BET:
    {
      return false;
    }
    case Isotherm::Type::Henry:
    {
      if(terms[0] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Freundlich:
    {
      if(terms[0] < 0.0 || terms[1] < 0.0 || terms[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Sips:
    {
      if(terms[0] < 0 || terms[0] > 1.0e20 || terms[1] < 0.0 || terms[1] > 1.0e10 || terms[2] < 0.0 || terms[2] > 100.0 ) return true;
      return false;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      if(terms[0] < 1.0E-20 || terms[0] > 1.0e20 || terms[1] < 0.0 || terms[1] > 1.0e10 || terms[2] < 0.0 || terms[2] > 100.0 ) return true;
      return false;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      if(terms[0] < 0.0 || terms[1] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Toth:
    {
      if(terms[0] < 0 || terms[1] < 0.0 || terms[2] < 0.0 || terms[2] > 100.0 ) return true;
      return false;
    }
    case Isotherm::Type::Unilan:
    {
      if(terms[0] < 0.0 || terms[1] < 0.0 || terms[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      if(terms[0] < 0.0 || terms[1] < 0.0 || terms[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Quadratic:
    {
      if(terms[0] < 0.0 || terms[1] < 0.0 || terms[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Temkin:
    {
      if(terms[0] < 0.0 || terms[1] < 0.0 || terms[2] < 0.0) return true;
      return false;
    }
    default:
      throw std::runtime_error("Error: unkown isotherm type");
  }
}

void Isotherm::randomize(double params[], double maximumLoading)
{
  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      params[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      params[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      break;
    }
    case Isotherm::Type::BET:
    {
      params[0] = 10.0 * RandomNumber::Uniform();
      params[1] = 10.0 * RandomNumber::Uniform();
      params[2] = 10.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Henry:
    {
      params[0] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      params[0] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      params[1] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Sips:
    {
      params[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      params[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      params[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      params[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      params[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      params[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      params[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      params[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      params[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Toth:
    {
      params[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      params[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      params[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Unilan:
    {
      params[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      params[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      params[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      params[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      params[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      params[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Quadratic:
    {
      params[0] = 2.1 * maximumLoading * RandomNumber::Uniform();
      params[1] = 10.0 * RandomNumber::Uniform();
      params[2] = 10.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Temkin:
    {
      params[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      params[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      params[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    default:
      throw std::runtime_error("Error: unkown isotherm type");
  }
}

std::string Isotherm::gnuplotFunctionString(char c, size_t i) const
{
  char stringBuffer[1024];

  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      sprintf(stringBuffer, "%c[%ld]*%c[%ld]*x/(1.0+%c[%ld]*x)", c, i, c, i+1, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::BET:
    {
      sprintf(stringBuffer, "%c[%ld]*%c[%ld]*x/((1.0-%c[%ld]*x)*(1.0-%c[%ld]+%c[%ld]*x))", 
              c, i, c, i+1, c, i+2, c, i+2, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::Henry:
    {
      sprintf(stringBuffer, "%c[%ld]*x", c, i);
      return stringBuffer;
    }
    case Isotherm::Type::Freundlich:
    {
      sprintf(stringBuffer, "%c[%ld]*x**[%ld]", c, i, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::Sips:
    {
      sprintf(stringBuffer, "%c[%ld]*((%c[%ld]*x)**(1.0/%c[%ld]))/(1.0+(%c[%ld]*x)**(1.0/%c[%ld]))",
              c, i, c, i+1, c, i+2, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      sprintf(stringBuffer, "%c[%ld]*%c[%ld]*x**%c[%ld]/(1.0+%c[%ld]*x**%c[%ld])", 
              c, i, c, i+1, c, i+2, c, i +1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      sprintf(stringBuffer, "%c[%ld]*x/(1.0+%c[%ld]*x**%c[%ld])", c, i, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Toth:
    {
      sprintf(stringBuffer, "%c[%ld]*%c[%ld]*x/((1.0+(%c[%ld]*x)**%c[%ld])**(1.0/%c[%ld]))",
              c, i, c, i+1, c, i+1, c, i+2, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Unilan:
    {
      sprintf(stringBuffer, "(%c[%ld]/(2.0*%c[%ld]))*log((1.0+%c[%ld]*exp(%c[%ld])*x)/(1.0+%c[%ld]*exp(-%c[%ld])*x))",
              c, i, c, i+2, c, i+1, c, i+2, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      sprintf(stringBuffer, "%c[%ld]*(%c[%ld]*x/(1.0+%c[%ld]*x) + (%c[%ld]**2)*%c[%ld]*x*(1.0-%c[%ld]*x)/(2.0*(1.0+%c[%ld]*x)**3))",
          c, i, c, i+1, c, i+1, c, i+2, c, i+1, c, i+1, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::Quadratic:
    {
      sprintf(stringBuffer, "%c[%ld]*(%c[%ld]*x+2.0*%c[%ld]*x**2)/(1.0+%c[%ld]*x+%c[%ld]*x**2)",
              c, i, c, i+1, c, i+2, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Temkin:
    {
      sprintf(stringBuffer, "%c[%ld]*(%c[%ld]*x/(1.0+%c[%ld]*x))+%c[%ld]*%c[%ld]*((%c[%ld]*x/(1.0+%c[%ld]*x))**2)*(%c[%ld]*x/(1.0+%c[%ld]*x)-1.0)", 
          c, i, c, i+1, c, i+1, c, i, c, i+2, c, i+1, c, i+1, c, i+1, c, i+1);
      return stringBuffer;
    }
    default:
      throw std::runtime_error("Error: unkown isotherm type");
  }
}
