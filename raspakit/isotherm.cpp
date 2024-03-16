module;

module isotherm;

import <cstdlib>;
import <iostream>;
import <sstream>;
import <cmath>;
import <vector>;
import <print>;

import randomnumbers;
import stringutils;


Isotherm::Isotherm(Isotherm::Type t, const std::vector<double> &values, size_t numberOfValues):
    type(t),
    parameters(values),
    numberOfParameters(numberOfValues)
{
}

std::string Isotherm::print() const
{
  std::ostringstream stream;

  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      std::print(stream, "    Langmuir isotherm\n");
      std::print(stream, "        q_sat: {}\n", parameters[0]);
      std::print(stream, "        b:     {}\n", parameters[1]);
      break;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      std::print(stream, "    Anti-Langmuir isotherm\n");
      std::print(stream, "        a:     {}\n", parameters[0]);
      std::print(stream, "        b:     {}\n", parameters[1]);
      break;
    }
    case Isotherm::Type::BET:
    {
      std::print(stream, "    BET isotherm\n");
      std::print(stream, "        q_sat: {}\n", parameters[0]);
      std::print(stream, "        b:     {}\n", parameters[1]);
      std::print(stream, "        c:     {}\n", parameters[2]);
      break;
    }
    case Isotherm::Type::Henry:
    {
      std::print(stream, "    Henry isotherm\n");
      std::print(stream, "        a:     {}\n", parameters[0]);
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      std::print(stream, "    Freundlich isotherm\n");
      std::print(stream, "        a:     {}\n", parameters[0]);
      std::print(stream, "        nu:    {}\n", parameters[1]);
      break;
    }
    case Isotherm::Type::Sips:
    {
      std::print(stream, "    Sips isotherm\n");
      std::print(stream, "        q_sat: {}\n", parameters[0]);
      std::print(stream, "        b:     {}\n", parameters[1]);
      std::print(stream, "        nu:    {}\n", parameters[2]);
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      std::print(stream, "    Langmuir-Freundlich isotherm\n");
      std::print(stream, "        q_sat: {}\n", parameters[0]);
      std::print(stream, "        b:     {}\n", parameters[1]);
      std::print(stream, "        nu:    {}\n", parameters[2]);
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      std::print(stream, "    Redlich-Peterson isotherm\n");
      std::print(stream, "        a:     {}\n", parameters[0]);
      std::print(stream, "        b:     {}\n", parameters[1]);
      std::print(stream, "        nu:    {}\n", parameters[2]);
      break;
    }
    case Isotherm::Type::Toth:
    {
      std::print(stream, "    Toth isotherm\n");
      std::print(stream, "        q_sat: {}\n", parameters[0]);
      std::print(stream, "        b:     {}\n", parameters[1]);
      std::print(stream, "        nu:    {}\n", parameters[2]);
      break;
    }
    case Isotherm::Type::Unilan:
    {
      std::print(stream, "    Unilan isotherm\n");
      std::print(stream, "        q_sat: {}\n", parameters[0]);
      std::print(stream, "        b:     {}\n", parameters[1]);
      std::print(stream, "        eta:   {}\n", parameters[2]);
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      std::print(stream, "    O'Brian & Myers isotherm\n");
      std::print(stream, "        q_sat: {}\n", parameters[0]);
      std::print(stream, "        b:     {}\n", parameters[1]);
      std::print(stream, "        sigma: {}\n", parameters[2]);
      break;
    }
    case Isotherm::Type::Quadratic:
    {
      std::print(stream, "    Quadratic isotherm\n");
      std::print(stream, "        q_sat: {}\n", parameters[0]);
      std::print(stream, "        b:     {}\n", parameters[1]);
      std::print(stream, "        c:     {}\n", parameters[2]);
      break;
    }
    case Isotherm::Type::Temkin:
    {
      std::print(stream, "    Temkin isotherm\n");
      std::print(stream, "        q_sat: {}\n", parameters[0]);
      std::print(stream, "        b:     {}\n", parameters[1]);
      std::print(stream, "        c:     {}\n", parameters[2]);
      break;
    }
    default:
      break;
  }

  return stream.str();
}

std::string Isotherm::printAsInputFormat() const
{
  std::ostringstream stream;

  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      std::print(stream, "                Langmuir              {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1]);
      break;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      std::print(stream, "                Anti-Langmuir         {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1]);
      break;
    }
    case Isotherm::Type::BET:
    {
      std::print(stream, "                BET                   {:<11g}  {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1], parameters[2]);
      break;
    }
    case Isotherm::Type::Henry:
    {
      std::print(stream, "                Henry                 {:<11g}\n", 
                         parameters[0]);
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      std::print(stream, "                Freundlich            {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1]);
      break;
    }
    case Isotherm::Type::Sips:
    {
      std::print(stream, "                Sips                  {:<11g}  {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1], parameters[2]);
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      std::print(stream, "                Langmuir-Freundlich   {:<11g}  {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1], parameters[2]);
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      std::print(stream, "                Redlich-Peterson      {:<11g}  {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1], parameters[2]);
      break;
    }
    case Isotherm::Type::Toth:
    {
      std::print(stream, "                Toth                  {:<11g}  {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1], parameters[2]);
      break;
    }
    case Isotherm::Type::Unilan:
    {
      std::print(stream, "                Unilan                {:<11g}  {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1], parameters[2]);
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      std::print(stream, "                O'Brien&Myers         {:<11g}  {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1], parameters[2]);
      break;
    }
    case Isotherm::Type::Quadratic:
    {
      std::print(stream, "                Quadratic             {:<11g}  {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1], parameters[2]);
      break;
    }
    case Isotherm::Type::Temkin:
    {
      std::print(stream, "                Temkin                {:<11g}  {:<11g}  {:<11g}\n", 
                         parameters[0], parameters[1], parameters[2]);
      break;
    }
    default:
      break;
  }
  return stream.str();
}

bool Isotherm::isUnphysical() const
{
  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      if(parameters[0] < 0 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10 ) return true;
      return false;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      if(parameters[0] < 0 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10 ) return true;
      return false;
    }
    case Isotherm::Type::BET:
    {
      return false;
    }
    case Isotherm::Type::Henry:
    {
      if(parameters[0] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Freundlich:
    {
      if(parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Sips:
    {
      if(parameters[0] < 0 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10 || 
         parameters[2] < 0.0 || parameters[2] > 100.0 ) return true;
      return false;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      if(parameters[0] < 1.0e-20 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10 || 
         parameters[2] < 0.0 || parameters[2] > 100.0 ) return true;
      return false;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      if(parameters[0] < 0.0 || parameters[1] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Toth:
    {
      if(parameters[0] < 0 || parameters[1] < 0.0 || parameters[2] < 0.0 || parameters[2] > 100.0 ) return true;
      return false;
    }
    case Isotherm::Type::Unilan:
    {
      if(parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      if(parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Quadratic:
    {
      if(parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Temkin:
    {
      if(parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    default:
      throw std::runtime_error("Error: unkown isotherm type\n");
  }
}

void Isotherm::randomize(RandomNumber &random, double maximumLoading)
{
  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      parameters[0] = 1.1 * maximumLoading * random.uniform();
      parameters[1] = std::pow(random.uniform(), 10.0 * 2.0 * (random.uniform() - 1.0));
      break;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      parameters[0] = 1.1 * maximumLoading * random.uniform();
      parameters[1] = std::pow(random.uniform(), 10.0 * 2.0 * (random.uniform() - 1.0));
      break;
    }
    case Isotherm::Type::BET:
    {
      parameters[0] = 10.0 * random.uniform();
      parameters[1] = 10.0 * random.uniform();
      parameters[2] = 10.0 * random.uniform();
      break;
    }
    case Isotherm::Type::Henry:
    {
      parameters[0] = std::pow(random.uniform(), 10.0 * 2.0 * (random.uniform() - 1.0));
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      parameters[0] = std::pow(random.uniform(), 10.0 * 2.0 * (random.uniform() - 1.0));
      parameters[1] = 0.1 + 2.0 * random.uniform();
      break;
    }
    case Isotherm::Type::Sips:
    {
      parameters[0] = 1.1 * maximumLoading * random.uniform();
      parameters[1] = std::pow(random.uniform(), 10.0 * 2.0 * (random.uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * random.uniform();
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      parameters[0] = 1.1 * maximumLoading * random.uniform();
      parameters[1] = std::pow(random.uniform(), 10.0 * 2.0 * (random.uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * random.uniform();
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      parameters[0] = 1.1 * maximumLoading * random.uniform();
      parameters[1] = std::pow(random.uniform(), 10.0 * 2.0 * (random.uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * random.uniform();
      break;
    }
    case Isotherm::Type::Toth:
    {
      parameters[0] = 1.1 * maximumLoading * random.uniform();
      parameters[1] = std::pow(random.uniform(), 10.0 * 2.0 * (random.uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * random.uniform();
      break;
    }
    case Isotherm::Type::Unilan:
    {
      parameters[0] = 1.1 * maximumLoading * random.uniform();
      parameters[1] = std::pow(random.uniform(), 10.0 * 2.0 * (random.uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * random.uniform();
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      parameters[0] = 1.1 * maximumLoading * random.uniform();
      parameters[1] = std::pow(random.uniform(), 10.0 * 2.0 * (random.uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * random.uniform();
      break;
    }
    case Isotherm::Type::Quadratic:
    {
      parameters[0] = 2.1 * maximumLoading * random.uniform();
      parameters[1] = 10.0 * random.uniform();
      parameters[2] = 10.0 * random.uniform();
      break;
    }
    case Isotherm::Type::Temkin:
    {
      parameters[0] = 1.1 * maximumLoading * random.uniform();
      parameters[1] = std::pow(random.uniform(), 10.0 * 2.0 * (random.uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * random.uniform();
      break;
    }
    default:
      throw std::runtime_error("Error: unkown isotherm type\n");
  }
}

std::string Isotherm::gnuplotFunctionString(char c, size_t i) const
{
  char stringBuffer[1024];

  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*%c[%zd]*x/(1.0+%c[%zd]*x)", c, i, c, i+1, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*x/(1.0-%c[%zd]*x)", c, i, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::BET:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*%c[%zd]*x/((1.0-%c[%zd]*x)*(1.0-%c[%zd]+%c[%zd]*x))", 
              c, i, c, i+1, c, i+2, c, i+2, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::Henry:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*x", c, i);
      return stringBuffer;
    }
    case Isotherm::Type::Freundlich:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*x**[%zd]", c, i, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::Sips:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*((%c[%zd]*x)**(1.0/%c[%zd]))/(1.0+(%c[%zd]*x)**(1.0/%c[%zd]))",
              c, i, c, i+1, c, i+2, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*%c[%zd]*x**%c[%zd]/(1.0+%c[%zd]*x**%c[%zd])", 
              c, i, c, i+1, c, i+2, c, i +1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*x/(1.0+%c[%zd]*x**%c[%zd])", c, i, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Toth:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*%c[%zd]*x/((1.0+(%c[%zd]*x)**%c[%zd])**(1.0/%c[%zd]))",
              c, i, c, i+1, c, i+1, c, i+2, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Unilan:
    {
      snprintf(stringBuffer, 1024, "(%c[%zd]/(2.0*%c[%zd]))*log((1.0+%c[%zd]*exp(%c[%zd])*x)/"
                                   "(1.0+%c[%zd]*exp(-%c[%zd])*x))", c, i, c, i+2, c, i+1, c, i+2, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*(%c[%zd]*x/(1.0+%c[%zd]*x) + (%c[%zd]**2)*%c[%zd]*x*(1.0-%c[%zd]*x)/"
                                   "(2.0*(1.0+%c[%zd]*x)**3))", c, i, c, i+1, c, i+1, c, i+2, c, i+1, c, i+1, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::Quadratic:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*(%c[%zd]*x+2.0*%c[%zd]*x**2)/(1.0+%c[%zd]*x+%c[%zd]*x**2)",
              c, i, c, i+1, c, i+2, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Temkin:
    {
      snprintf(stringBuffer, 1024, "%c[%zd]*(%c[%zd]*x/(1.0+%c[%zd]*x))+%c[%zd]*%c[%zd]*((%c[%zd]*x/(1.0+%c[%zd]*x))"
                                   "**2)*(%c[%zd]*x/(1.0+%c[%zd]*x)-1.0)", 
                                   c, i, c, i+1, c, i+1, c, i, c, i+2, c, i+1, c, i+1, c, i+1, c, i+1);
      return stringBuffer;
    }
    default:
      throw std::runtime_error("Error: unkown isotherm type\n");
  }
}
