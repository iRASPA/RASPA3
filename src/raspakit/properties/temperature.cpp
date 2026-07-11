module;

module property_temperature;

import std;

import archive;
import units;
import stringutils;

std::string PropertyTemperature::writeAveragesStatistics(const std::string tag)
{
  std::ostringstream stream;
  double conv = Units::TemperatureConversionFactor;

  std::pair<double, double> temperatureAverages = average();
  for (std::size_t i = 0; i < numberOfBlocks; ++i)
  {
    double blockAverage = averaged(i);
    std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    {} temperature  {: .6e} +/- {: .6e} [K]\n", tag, conv * temperatureAverages.first,
             temperatureAverages.second);
  std::print(stream, "\n");

  return stream.str();
}
