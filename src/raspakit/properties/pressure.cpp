module;

module property_pressure;

import std;

import archive;
import double3;
import double3x3;
import pressure_data;
import stringutils;
import units;
import json;

std::string PropertyPressure::writeAveragesStatistics() const
{
  std::ostringstream stream;

  double conv = Units::PressureConversionFactor;

  std::print(stream, "Pressure averages and statistics:\n");
  std::print(stream, "===============================================================================\n\n");

  std::pair<PressureData, PressureData> average_pressure = result();

  switch (Units::unitSystem)
  {
    case Units::System::RASPA:
    {
      double3x3 pressureTensor = 1e-5 * Units::PressureConversionFactor * average_pressure.first.totalPressureTensor;
      double3x3 pressureTensorError = 1e-5 * Units::PressureConversionFactor * average_pressure.second.totalPressureTensor;
      std::print(stream, "Average pressure tensor: \n");
      std::print(stream, "-------------------------------------------------------------------------------\n");
      std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.ax,
                 pressureTensor.bx, pressureTensor.cx, pressureTensorError.ax, pressureTensorError.bx,
                 pressureTensorError.cx);
      std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n", pressureTensor.ay,
                 pressureTensor.by, pressureTensor.cy, pressureTensorError.ay, pressureTensorError.by,
                 pressureTensorError.cy);
      std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [bar]\n\n", pressureTensor.az,
                 pressureTensor.bz, pressureTensor.cz, pressureTensorError.az, pressureTensorError.bz,
                 pressureTensorError.cz);

      for (std::size_t i = 0; i < bookKeepingPressure.size(); ++i)
      {
        double blockAverage = averagedPressures(i).idealGasPressure;
        std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Ideal gas pressure  {: .6e} +/- {: .6e} [Pa]\n", conv * average_pressure.first.idealGasPressure,
                 average_pressure.second.idealGasPressure);
      std::print(stream, "                        {: .6e} +/- {: .6e} [bar]\n",
                 1e-5 * conv * average_pressure.first.idealGasPressure, 1e-5 * conv * average_pressure.second.idealGasPressure);
      std::print(stream, "\n\n");

      for (std::size_t i = 0; i < bookKeepingPressure.size(); ++i)
      {
        double blockAverage = averagedPressures(i).excessPressure;
        std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Excess pressure  {: .6e} +/- {: .6e} [Pa]\n", conv * average_pressure.first.excessPressure,
                 conv * average_pressure.second.excessPressure);
      std::print(stream, "                     {: .6e} +/- {: .6e} [bar]\n", 1e-5 * conv * average_pressure.first.excessPressure,
                 1e-5 * conv * average_pressure.second.excessPressure);
      std::print(stream, "\n\n");

      for (std::size_t i = 0; i < bookKeepingPressure.size(); ++i)
      {
        double blockAverage = averagedPressures(i).totalPressure;
        std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Pressure average  {: .6e} +/- {: .6e} [Pa]\n", conv * average_pressure.first.totalPressure,
                 conv * average_pressure.second.totalPressure);
      std::print(stream, "                      {: .6e} +/- {: .6e} [bar]\n", 1e-5 * conv * average_pressure.first.totalPressure,
                 1e-5 * conv * average_pressure.second.totalPressure);
      std::print(stream, "\n\n");
    }
    break;
    case Units::System::ReducedUnits:
    {
      double3x3 pressureTensor = Units::PressureConversionFactor * average_pressure.first.totalPressureTensor;
      double3x3 pressureTensorError = Units::PressureConversionFactor * average_pressure.second.totalPressureTensor;
      std::print(stream, "Average pressure tensor: \n");
      std::print(stream, "-------------------------------------------------------------------------------\n");
      std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [{}]\n", pressureTensor.ax,
                 pressureTensor.bx, pressureTensor.cx, pressureTensorError.ax, pressureTensorError.bx,
                 pressureTensorError.cx, Units::unitOfPressureString);
      std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [{}]\n", pressureTensor.ay,
                 pressureTensor.by, pressureTensor.cy, pressureTensorError.ay, pressureTensorError.by,
                 pressureTensorError.cy, Units::unitOfPressureString);
      std::print(stream, "{: .4e} {: .4e} {: .4e} +/- {:.4e} {:.4e} {:.4e} [{}]\n\n", pressureTensor.az,
                 pressureTensor.bz, pressureTensor.cz, pressureTensorError.az, pressureTensorError.bz,
                 pressureTensorError.cz, Units::unitOfPressureString);

      for (std::size_t i = 0; i < bookKeepingPressure.size(); ++i)
      {
        double blockAverage = averagedPressures(i).idealGasPressure;
        std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Ideal gas pressure  {: .6e} +/- {: .6e} [{}]\n", conv * average_pressure.first.idealGasPressure,
                 average_pressure.second.idealGasPressure, Units::unitOfPressureString);
      std::print(stream, "\n\n");

      for (std::size_t i = 0; i < bookKeepingPressure.size(); ++i)
      {
        double blockAverage = averagedPressures(i).excessPressure;
        std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Excess pressure  {: .6e} +/- {: .6e} [{}]\n", conv * average_pressure.first.excessPressure,
                 conv * average_pressure.second.excessPressure, Units::unitOfPressureString);
      std::print(stream, "\n\n");

      for (std::size_t i = 0; i < bookKeepingPressure.size(); ++i)
      {
        double blockAverage = averagedPressures(i).totalPressure;
        std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Pressure average  {: .6e} +/- {: .6e} [{}]\n", conv * average_pressure.first.totalPressure,
                 conv * average_pressure.second.totalPressure, Units::unitOfPressureString);
      std::print(stream, "\n\n");
    }
    break;
  }

  return stream.str();
}

nlohmann::json PropertyPressure::jsonAveragesStatistics() const
{
  nlohmann::json status;
  double conv = Units::PressureConversionFactor;

  std::pair<PressureData, PressureData> average_pressure = result();

  status["averagePressureTensor"]["mean"] = conv * average_pressure.first.totalPressureTensor;
  status["averagePressureTensor"]["confidence"] = conv * average_pressure.second.totalPressureTensor;

  /* TODO
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
  */

  return status;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPressure &e)
{
  archive << e.versionNumber;

  archive << e.numberOfBlocks;
  archive << e.bookKeepingPressure;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPressure &e)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > e.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyPressure' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> e.numberOfBlocks;
  archive >> e.bookKeepingPressure;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyPressure: Error in binary restart\n"));
  }
#endif

  return archive;
}
