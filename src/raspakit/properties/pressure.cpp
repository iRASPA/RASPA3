module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#endif

module property_pressure;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <iostream>;
import <sstream>;
import <fstream>;
import <tuple>;
import <vector>;
import <array>;
import <map>;
import <algorithm>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <print>;
#endif

import archive;
import double3;
import double3x3;
import stringutils;
import units;
import json;

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
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.ax, pressureTensor.bx,
             pressureTensor.cx, pressureTensorError.ax, pressureTensorError.bx, pressureTensorError.cx);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.ay, pressureTensor.by,
             pressureTensor.cy, pressureTensorError.ay, pressureTensorError.by, pressureTensorError.cy);
  std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n\n", pressureTensor.az, pressureTensor.bz,
             pressureTensor.cz, pressureTensorError.az, pressureTensorError.bz, pressureTensorError.cz);

  std::pair<double, double> pressureIdealGasAverage = averageIdealGasPressure();
  for (size_t i = 0; i < bookKeepingIdealGasPressure.size(); ++i)
  {
    double blockAverage = averagedIdealGasPressure(i);
    std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Ideal gas pressure  {: .6e} +/- {: .6e} [Pa]\n", conv * pressureIdealGasAverage.first,
             pressureIdealGasAverage.second);
  std::print(stream, "                        {: .6e} +/- {: .6e} [bar]\n", 1e-5 * conv * pressureIdealGasAverage.first,
             1e-5 * conv * pressureIdealGasAverage.second);
  std::print(stream, "\n\n");

  for (size_t i = 0; i < bookKeepingExcessPressure.size(); ++i)
  {
    double blockAverage = averagedExcessPressure(i);
    std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Excess pressure  {: .6e} +/- {: .6e} [Pa]\n", conv * pressureAverage.first,
             conv * pressureAverage.second);
  std::print(stream, "                     {: .6e} +/- {: .6e} [bar]\n", 1e-5 * conv * pressureAverage.first,
             1e-5 * conv * pressureAverage.second);
  std::print(stream, "\n\n");

  std::pair<double, double> pressureTotalAverage = averagePressure();
  for (size_t i = 0; i < bookKeepingExcessPressure.size(); ++i)
  {
    double blockAverage = averagedPressure(i);
    std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Pressure average  {: .6e} +/- {: .6e} [Pa]\n", conv * pressureTotalAverage.first,
             conv * pressureTotalAverage.second);
  std::print(stream, "                      {: .6e} +/- {: .6e} [bar]\n", 1e-5 * conv * pressureTotalAverage.first,
             1e-5 * conv * pressureTotalAverage.second);
  std::print(stream, "\n\n");

  return stream.str();
}

nlohmann::json PropertyPressure::jsonAveragesStatistics() const
{
  nlohmann::json status;
  double conv = Units::PressureConversionFactor;

  std::pair<double3x3, double3x3> currentPressureTensor = averagePressureTensor();
  status["averagePressureTensor"]["mean"] = conv * currentPressureTensor.first;
  status["averagePressureTensor"]["confidence"] = conv * currentPressureTensor.second;

  std::pair<double, double> pressureIdealGasAverage = averageIdealGasPressure();

  std::vector<double> blocksIdealGasPressure(numberOfBlocks);
  std::transform(bookKeepingIdealGasPressure.begin(), bookKeepingIdealGasPressure.end(), blocksIdealGasPressure.begin(),
                 [conv](const std::pair<double, double> &book) { return conv * book.first / book.second; });
  status["idealGasPressure"]["block"] = blocksIdealGasPressure;
  status["idealGasPressure"]["mean"] = conv * pressureIdealGasAverage.first;
  status["idealGasPressure"]["confidence"] = conv * pressureIdealGasAverage.second;

  std::pair<double, double> pressureExcessAverage = averageExcessPressure();
  std::vector<double> blocksExcessPressure(numberOfBlocks);
  std::transform(bookKeepingExcessPressure.begin(), bookKeepingExcessPressure.end(), blocksExcessPressure.begin(),
                 [conv](const std::pair<double3x3, double> &book)
                 { return conv * book.first.trace() / (3.0 * book.second); });
  status["excessPressure"]["block"] = blocksExcessPressure;
  status["excessPressure"]["mean"] = conv * pressureExcessAverage.first;
  status["excessPressure"]["confidence"] = conv * pressureExcessAverage.second;

  std::pair<double, double> pressureTotalAverage = averagePressure();
  std::vector<double> blocksPressure(numberOfBlocks);
  std::transform(bookKeepingIdealGasPressure.begin(), bookKeepingIdealGasPressure.end(),
                 bookKeepingExcessPressure.begin(), blocksPressure.begin(),
                 [conv](const std::pair<double, double> ideal, const std::pair<double3x3, double> excess)
                 { return conv * ((ideal.first / ideal.second) + (excess.first.trace() / (3.0 * excess.second))); });
  status["pressure"]["block"] = blocksPressure;
  status["pressure"]["mean"] = conv * pressureTotalAverage.first;
  status["pressure"]["confidence"] = conv * pressureTotalAverage.second;

  return status;
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
  if (versionNumber > e.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyPressure' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> e.numberOfBlocks;
  archive >> e.bookKeepingExcessPressure;
  archive >> e.bookKeepingIdealGasPressure;

  return archive;
}
