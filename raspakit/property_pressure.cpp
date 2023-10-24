module;

module property_pressure;

import <string>;
import <iostream>;
import <sstream>;
import <tuple>;
import <vector>;
import <print>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;

import archive;
import double3;
import double3x3;
import stringutils;
import units;

std::string PropertyPressure::writeAveragesStatistics() const
{
  std::ostringstream stream;

  double conv = Units::PressureConversionFactor;

  std::pair<double, double> pressureAverage = averageExcessPressure();
  std::print(stream, "Pressure averages and statistics:\n");
  std::print(stream, "===============================================================================\n\n");

  std::pair<double3x3, double3x3> currentPressureTensor = averagePressureTensor();
  double3x3 pressureTensor = 1e-5 * Units::PressureConversionFactor * currentPressureTensor.first;
  double3x3 pressureTensorError = 1e-5 * Units::PressureConversionFactor * currentPressureTensor.second;
  std::print(stream, "Average pressure tensor: \n");
  std::print(stream, "-------------------------------------------------------------------------------\n");
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n",
    pressureTensor.ax, pressureTensor.bx, pressureTensor.cx, pressureTensorError.ax, pressureTensorError.bx, pressureTensorError.cx);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n",
    pressureTensor.ay, pressureTensor.by, pressureTensor.cy, pressureTensorError.ay, pressureTensorError.by, pressureTensorError.cy);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n\n",
    pressureTensor.az, pressureTensor.bz, pressureTensor.cz, pressureTensorError.az, pressureTensorError.bz, pressureTensorError.cz);

  std::pair<double, double> pressureIdealGasAverage = averageIdealGasPressure();
  for (size_t i = 0; i < bookKeepingIdealGasPressure.size(); ++i)
  {
    double blockAverage = averagedIdealGasPressure(i);
    std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Ideal gas pressure  {: .6e} +/- {: .6e} [Pa]\n",
    conv * pressureIdealGasAverage.first, pressureIdealGasAverage.second);
  std::print(stream, "                        {: .6e} +/- {: .6e} [bar]\n",
    1e-5 * conv * pressureIdealGasAverage.first, 1e-5 * conv * pressureIdealGasAverage.second);
  std::print(stream, "\n\n");

  for (size_t i = 0; i < bookKeepingExcessPressure.size(); ++i)
  {
    double blockAverage = averagedExcessPressure(i);
    std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Excess pressure  {: .6e} +/- {: .6e} [Pa]\n",
    conv * pressureAverage.first, conv * pressureAverage.second);
  std::print(stream, "                     {: .6e} +/- {: .6e} [bar]\n",
    1e-5 * conv * pressureAverage.first, 1e-5 * conv * pressureAverage.second);
  std::print(stream, "\n\n");

  std::pair<double, double> pressureTotalAverage = averagePressure();
  for (size_t i = 0; i < bookKeepingExcessPressure.size(); ++i)
  {
    double blockAverage = averagedPressure(i);
    std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Pressure average  {: .6e} +/- {: .6e} [Pa]\n",
    conv * pressureTotalAverage.first, conv * pressureTotalAverage.second);
  std::print(stream, "                      {: .6e} +/- {: .6e} [bar]\n",
    1e-5 * conv * pressureTotalAverage.first, 1e-5 * conv * pressureTotalAverage.second);
  std::print(stream, "\n\n");

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPressure &e)
{
  archive << e.versionNumber;

  archive << e.numberOfBlocks;
  archive << e.bookKeepingExcessPressure;
  archive << e.bookKeepingIdealGasPressure;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPressure &e)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > e.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyPressure' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> e.numberOfBlocks;
  archive >> e.bookKeepingExcessPressure;
  archive >> e.bookKeepingIdealGasPressure;

  return archive;
}
