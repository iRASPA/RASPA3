module;

module system;

import double3x3;
import component;
import energy_factor;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import averages;
import print;
import units;
import property_pressure;

import <iostream>;
import <random>;
import <sstream>;

void System::writePressureAveragesStatistics(std::ostream &stream) const
{
  double conv = Units::PressureConversionFactor;

  std::pair<double, double> pressureAverage = averagePressure.averageExcessPressure();
  std::print(stream, "Pressure averages and statistics:\n");
  std::print(stream, "===============================================================================\n\n");

  std::pair<double3x3, double3x3> currentPressureTensor = averagePressure.averagePressureTensor();
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

  std::pair<double, double> pressureIdealGasAverage = averagePressure.averageIdealGasPressure();
  for (size_t i = 0; i < averagePressure.bookKeepingIdealGasPressure.size(); ++i)
  {
      double blockAverage = averagePressure.averagedIdealGasPressure(i);
      std::print(stream, "    Block[ {:2d}] {}\n", i, conv * blockAverage);
  }
  std::print(stream, "    -----------------------------------------------------------------------\n");
  std::print(stream, "    Ideal gas pressure  {} +/- {} [Pa]\n",
          conv * pressureIdealGasAverage.first, pressureIdealGasAverage.second);
  std::print(stream, "                        {} +/- {} [bar]\n",
          1e-5 * conv * pressureIdealGasAverage.first, 1e-5 * conv * pressureIdealGasAverage.second);
  std::print(stream, "\n\n");

  for (size_t i = 0; i < averagePressure.bookKeepingExcessPressure.size(); ++i)
  {
      double blockAverage = averagePressure.averagedExcessPressure(i);
      std::print(stream, "    Block[ {:2d}] {}\n", i, conv * blockAverage);
  }
  std::print(stream, "    -----------------------------------------------------------------------\n");
  std::print(stream, "    Excess pressure  {} +/- {} [Pa]\n",
          conv * pressureAverage.first, conv * pressureAverage.second);
  std::print(stream, "                     {} +/- {} [bar]\n",
          1e-5 * conv * pressureAverage.first, 1e-5 * conv * pressureAverage.second);
  std::print(stream, "\n\n");

  std::pair<double, double> pressureTotalAverage = averagePressure.averagePressure();
  for (size_t i = 0; i < averagePressure.bookKeepingExcessPressure.size(); ++i)
  {
      double blockAverage = averagePressure.averagedPressure(i);
      std::print(stream, "    Block[ {:2d}] {}\n", i, conv * blockAverage);
  }
  std::print(stream, "    -----------------------------------------------------------------------\n");
  std::print(stream, "    Pressure average  {} +/- {} [Pa]\n",
          conv * pressureTotalAverage.first, pressureTotalAverage.second);
  std::print(stream, "                      {} +/- {} [bar]\n",
          1e-5 * conv * pressureTotalAverage.first, 1e-5 * conv * pressureTotalAverage.second);
  std::print(stream, "\n\n");
}

