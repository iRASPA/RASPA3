module;

module dudlambda;

import <vector>;
import <iostream>;
import <cmath>;
import <string>;
import <sstream>;
import <format>;
import <algorithm>;
import <numeric>;
import <cmath>;
import <numbers>;

import units;
import print;
import averages;

void dUdLambda::WangLandauIteration(dUdLambda::WangLandauPhase phase)
{
  switch (phase)
  {
  case dUdLambda::WangLandauPhase::Initialize:
    WangLandauScalingFactor = 1.0;
    std::fill(histogram.begin(), histogram.end(), 0.0);
    std::fill(biasFactor.begin(), biasFactor.end(), 0.0);
    break;
  case dUdLambda::WangLandauPhase::Sample:
    biasFactor[currentBin] -= WangLandauScalingFactor;
    histogram[currentBin] += 1.0;
    break;
  case dUdLambda::WangLandauPhase::AdjustBiasingFactors:
  {
    std::vector<double>::iterator minValueIterator = std::min_element(histogram.begin(), histogram.end());
    double minimumValue = *minValueIterator;
    if (minimumValue > 0.01)
    {
      double sumOfHistogram = std::accumulate(histogram.begin(), histogram.end(), 0.0);
      if (minimumValue / sumOfHistogram > 0.01)
      {
        WangLandauScalingFactor *= 0.5;
      }
    }
    std::fill(histogram.begin(), histogram.end(), 0.0);
  }
  break;
  case dUdLambda::WangLandauPhase::Finalize:
    std::fill(histogram.begin(), histogram.end(), 0.0);

    double normalize = biasFactor[0];
    for (double& bias : biasFactor)
    {
      bias -= normalize;
    }
    break;
  }
}


std::string dUdLambda::writeAveragesStatistics(double beta) const
{
  std::ostringstream stream;

  std::print(stream, "Thermodynamic integration (dU/dlambda)\n");
  std::print(stream, "===============================================================================\n\n");


  double conv = Units::EnergyToKelvin;
  std::pair<std::vector<double3>, std::vector<double3>> dudlambda = dUdlambdaBookKeeping.averageDuDlambda();
  for (size_t binIndex = 0; binIndex < dUdlambdaBookKeeping.numberOfBins; ++binIndex)
  {
    std::print(stream, "{}{:2d}-{:5f} (lambda) <dU/dlambda>: {: .6e} +/- {:.6e} [-]\n",
      "    ", binIndex, static_cast<double>(binIndex) * delta, conv * dudlambda.first[binIndex].x, conv * dudlambda.second[binIndex].x);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Excess chemical potential: integral du/dlambda over lambda (trapezoidal rule)\n");
  for (size_t blockIndex = 0; blockIndex < dUdlambdaBookKeeping.numberOfBlocks; ++blockIndex)
  {
    double blockAverage = dUdlambdaBookKeeping.averagedExcessChemicalPotential(blockIndex);
    std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, Units::EnergyToKelvin * blockAverage);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::pair<double, double> averageExcessChemicalPotentialDUDlambda = dUdlambdaBookKeeping.averageExcessChemicalPotential();
  std::pair<double, double> averageIdealGasChemicalPotentialDUDlambda = dUdlambdaBookKeeping.averageIdealGasChemicalPotential(beta);
  std::pair<double, double> averageTotalChemicalPotentialDUDlambda = dUdlambdaBookKeeping.averageTotalChemicalPotential(beta);
  std::pair<double, double> averageFugacityDUDlambda = dUdlambdaBookKeeping.averageFugacity(beta);
  std::print(stream, "    Excess chemical potential:   {: .6e} +/- {: .6e} [K]\n",
    Units::EnergyToKelvin * averageExcessChemicalPotentialDUDlambda.first,
    Units::EnergyToKelvin * averageExcessChemicalPotentialDUDlambda.second);
  std::print(stream, "    Ideal chemical potential:    {: .6e} +/- {: .6e} [K]\n",
    Units::EnergyToKelvin * averageIdealGasChemicalPotentialDUDlambda.first,
    Units::EnergyToKelvin * averageIdealGasChemicalPotentialDUDlambda.second);
  std::print(stream, "    Total chemical potential:    {: .6e} +/- {: .6e} [K]\n",
    Units::EnergyToKelvin * averageTotalChemicalPotentialDUDlambda.first,
    Units::EnergyToKelvin * averageTotalChemicalPotentialDUDlambda.second);
  
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Excess chemical potential:   {: .6e} +/- {: .6e} [kJ/mol]\n",
    Units::EnergyToKJPerMol * averageExcessChemicalPotentialDUDlambda.first,
    Units::EnergyToKJPerMol * averageExcessChemicalPotentialDUDlambda.second);
  std::print(stream, "    Ideal chemical potential:    {: .6e} +/- {: .6e} [kJ/mol]\n",
    Units::EnergyToKJPerMol * averageIdealGasChemicalPotentialDUDlambda.first,
    Units::EnergyToKJPerMol * averageIdealGasChemicalPotentialDUDlambda.second);
  std::print(stream, "    Total chemical potential:    {: .6e} +/- {: .6e} [kJ/mol]\n",
    Units::EnergyToKJPerMol * averageTotalChemicalPotentialDUDlambda.first,
    Units::EnergyToKJPerMol * averageTotalChemicalPotentialDUDlambda.second);

  return stream.str();
}