module;

module partition_function;

import std;

import janaf;
import nasa_polynomials;

namespace PartitionFunction
{

bool contains(std::string_view speciesName, Source source)
{
  switch (source)
  {
    case Source::NASA:
      return NASAPolynomials::contains(speciesName) || JANAF::contains(speciesName);
    case Source::JANAF:
      return JANAF::contains(speciesName);
  }
  std::unreachable();
}

double logPartitionFunction(std::string_view speciesName, double temperature, Source source)
{
  switch (source)
  {
    case Source::JANAF:
      return JANAF::logPartitionFunction(speciesName, temperature);
    case Source::NASA:
      break;
  }

  if (NASAPolynomials::contains(speciesName))
  {
    try
    {
      return NASAPolynomials::logPartitionFunction(speciesName, temperature);
    }
    catch (const std::runtime_error &)
    {
      // typically the temperature is below the 200 K lower bound of the NASA polynomials;
      // the JANAF tables extend down to 100 K
      if (!JANAF::contains(speciesName))
      {
        throw;
      }
    }
    return JANAF::logPartitionFunction(speciesName, temperature);
  }

  if (JANAF::contains(speciesName))
  {
    return JANAF::logPartitionFunction(speciesName, temperature);
  }

  throw std::runtime_error(
      std::format("[PartitionFunction]: species '{}' not found in the NASA-polynomial (Burcat) database "
                  "or the JANAF tables; specify 'LnPartitionFunction' as a number instead\n",
                  speciesName));
}

}  // namespace PartitionFunction
