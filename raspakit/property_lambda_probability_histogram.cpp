module;

#if defined(_WIN32)
  import <cassert>;
#else
  #include <assert.h>
#endif

module property_lambda_probability_histogram;

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
import <optional>;
import <print>;

import double3;

import units;
import stringutils;
import averages;

void PropertyLambdaProbabilityHistogram::WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase phase, bool containsTheFractionalMolecule, [[maybe_unused]] double value)
{
  switch (phase)
  {
  case PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize:
    WangLandauScalingFactor = 0.01;
    std::fill(histogram.begin(), histogram.end(), 0.0);
    std::fill(biasFactor.begin(), biasFactor.end(), 0.0);
    break;
  case PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample:
    assert(currentBin >= 0 && currentBin < numberOfBins);
    if(containsTheFractionalMolecule)
    {
      biasFactor[currentBin] -= WangLandauScalingFactor;
      histogram[currentBin] += 1.0;
    }
    break;
  case PropertyLambdaProbabilityHistogram::WangLandauPhase::AdjustBiasingFactors:
  {
    WangLandauScalingFactor *= 0.5;
    std::fill(histogram.begin(), histogram.end(), 0.0);
    occupancyCount = 0;
    occupancyTotal = 0;
    //std::vector<double>::iterator maxValueIterator = std::max_element(histogram.begin(), histogram.end());
    //double maximumValue = *maxValueIterator;
    //
    //bool allBinsAtLeast90PercentOfThisValue = true;
    //for(size_t i = 0; i < histogram.size(); ++i)
    //{ 
    //  if (histogram[i] < 0.9 * maximumValue)
    //  {
    //    allBinsAtLeast90PercentOfThisValue = false;
    //    break;
    //  }
    //}
    //if (allBinsAtLeast90PercentOfThisValue)
    //{
    //  std::cout << "allBinsAtLeast90PercentOfThisValue" << std::endl;
    //  WangLandauScalingFactor *= 0.5;
    //  std::fill(histogram.begin(), histogram.end(), 0.0);
    //  occupancyCount = 0;
    //  occupancyTotal = 0;
    //}
  }
  break;
  case PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize:
    std::fill(histogram.begin(), histogram.end(), 0.0);
    break;
  }
}


std::string PropertyLambdaProbabilityHistogram::writeAveragesStatistics(double beta, std::optional<double> imposedChemicalPotential, std::optional<double> imposedFugacity) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  size_t lastBin = numberOfBins - 1;

  std::print(stream, "    Lambda histogram and bias:\n");
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::pair<std::vector<double>, std::vector<double>> histogram_avg = averageProbabilityHistogram();
  for (size_t i = 0; i < numberOfBins; ++i)
  {
    std::print(stream, "{}{:2d}-{:4f} (lambda) P: {: .5e} +/- {:.5e} bias: {: .5e} [-]\n",
      "    ", i, static_cast<double>(i) * delta,
      histogram_avg.first[i],
      histogram_avg.second[i],
      biasFactor[i]);
  }
  std::print(stream, "\n\n");

  std::pair<double, double> idealGasChemicalPotential = averageIdealGasChemicalPotential(beta);

  std::pair<std::vector<double>, std::vector<double>> freeEnergy = averageLandauFreeEnergyHistogram(beta);
  std::pair<double, double> excessChemicalPotential = averageExcessChemicalPotential(beta);
  double excessChemicalPotentialBias = (biasFactor[lastBin] - biasFactor[0]) / beta;
  std::pair<double, double> totalChemicalPotential = averageTotalChemicalPotential(beta, excessChemicalPotentialBias);

  std::pair<double, double> measuredFugacity = averageFugacity(beta, excessChemicalPotentialBias);

  std::print(stream, "    Lambda statistics:\n");
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  for (size_t i = 0; i < numberOfBins; ++i)
  {
    std::print(stream, "{}{:2d}-{:4f} (lambda) Free energy: {:.6e} +/- {:.6e} [K]\n",
      "    ", i, static_cast<double>(i) * delta,
      conv * freeEnergy.first[i], conv * freeEnergy.second[i]);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Excess chemical potential: (ln(P(lambda=1))-ln(P(lambda=0)))/Beta\n");
  for (size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
  {
    double blockAverage = averagedExcessChemicalPotential(blockIndex, beta);
    std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, conv * (blockAverage + excessChemicalPotentialBias));
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Excess chemical potential:    {: .6e} +/- {: .6e} [K]\n",
    Units::EnergyToKelvin * (excessChemicalPotential.first + excessChemicalPotentialBias),
    Units::EnergyToKelvin * excessChemicalPotential.second);
  std::print(stream, "    Ideal gas chemical potential: {: .6e} +/- {: .6e} [K]\n",
    Units::EnergyToKelvin * idealGasChemicalPotential.first,
    Units::EnergyToKelvin * idealGasChemicalPotential.second);
  std::print(stream, "    Total chemical potential:     {: .6e} +/- {: .6e} [K]\n",
    Units::EnergyToKelvin * totalChemicalPotential.first,
    Units::EnergyToKelvin * totalChemicalPotential.second);
  if (imposedChemicalPotential)
  {
    std::print(stream, "    Imposed chemical potential:   {: .6e} [K]\n", Units::EnergyToKelvin * imposedChemicalPotential.value());
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    Excess chemical potential:    {: .6e} +/- {: .6e} [kJ/mol]\n",
    Units::EnergyToKJPerMol * (excessChemicalPotential.first + excessChemicalPotentialBias),
    Units::EnergyToKJPerMol * excessChemicalPotential.second);
  std::print(stream, "    Ideal gas chemical potential: {: .6e} +/- {: .6e} [kJ/mol]\n",
    Units::EnergyToKJPerMol * idealGasChemicalPotential.first,
    Units::EnergyToKJPerMol * idealGasChemicalPotential.second);
  std::print(stream, "    Total chemical potential:     {: .6e} +/- {: .6e} [kJ/mol]\n",
    Units::EnergyToKJPerMol * totalChemicalPotential.first,
    Units::EnergyToKJPerMol * totalChemicalPotential.second);
  if (imposedChemicalPotential)
  {
    std::print(stream, "    Imposed chemical potential:   {: .6e} [K]\n", Units::EnergyToKJPerMol * imposedChemicalPotential.value());
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  if (imposedFugacity)
  {
    std::print(stream, "    Imposed fugacity:             {: .6e} [Pa]\n", Units::PressureConversionFactor * imposedFugacity.value());
  }
  std::print(stream, "    Measured fugacity:            {: .6e} +/- {: .6e} [Pa]\n",
    Units::PressureConversionFactor * measuredFugacity.first,
    Units::PressureConversionFactor * measuredFugacity.second);
 
  std::print(stream, "\n\n");
  return stream.str();
}

std::string PropertyLambdaProbabilityHistogram::writeDUdLambdaStatistics(double beta, std::optional<double> imposedChemicalPotential, std::optional<double> imposedFugacity) const
{
  std::ostringstream stream;

  if (computeDUdlambda)
  {
    std::print(stream, "    Thermodynamic integration (dU/dlambda)\n");
    std::print(stream, "    ===========================================================================\n\n");
    
    double conv = Units::EnergyToKelvin;
    std::pair<std::vector<double3>, std::vector<double3>> dudlambda = averageDuDlambda();
    for (size_t binIndex = 0; binIndex < numberOfBins; ++binIndex)
    {
      std::print(stream, "{}{:2d}-{:5f} (lambda) <dU/dlambda>: {: .6e} +/- {:.6e} [-]\n",
        "    ", binIndex, static_cast<double>(binIndex) * delta, conv * dudlambda.first[binIndex].x, conv * dudlambda.second[binIndex].x);
    }
    std::print(stream, "    ---------------------------------------------------------------------------\n");
    std::print(stream, "    Excess chemical potential: integral du/dlambda over lambda (Simpson's rule)\n");
    for (size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
    {
      double blockAverage = averagedExcessChemicalPotentialDUdlambda(blockIndex);
      std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, Units::EnergyToKelvin * blockAverage);
    }
    std::print(stream, "    ---------------------------------------------------------------------------\n");
    std::pair<double, double> averageExcessChemicalPotentialDUDlambda = averageExcessChemicalPotentialDUdlambda();
    std::pair<double, double> averageIdealGasChemicalPotentialDUDlambda = averageIdealGasChemicalPotential(beta);
    std::pair<double, double> averageTotalChemicalPotentialDUDlambda = averageTotalChemicalPotential(beta);
    std::pair<double, double> averageFugacityDUDlambda = averageFugacityDUdlambda(beta);
    std::print(stream, "    Excess chemical potential:   {: .6e} +/- {: .6e} [K]\n",
      Units::EnergyToKelvin * averageExcessChemicalPotentialDUDlambda.first,
      Units::EnergyToKelvin * averageExcessChemicalPotentialDUDlambda.second);
    std::print(stream, "    Ideal chemical potential:    {: .6e} +/- {: .6e} [K]\n",
      Units::EnergyToKelvin * averageIdealGasChemicalPotentialDUDlambda.first,
      Units::EnergyToKelvin * averageIdealGasChemicalPotentialDUDlambda.second);
    std::print(stream, "    Total chemical potential:    {: .6e} +/- {: .6e} [K]\n",
      Units::EnergyToKelvin * averageTotalChemicalPotentialDUDlambda.first,
      Units::EnergyToKelvin * averageTotalChemicalPotentialDUDlambda.second);
    
    if (imposedChemicalPotential)
    {
      std::print(stream, "    Imposed chemical potential:  {: .6e} [K]\n", Units::EnergyToKelvin * imposedChemicalPotential.value());
    }
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
    if (imposedChemicalPotential)
    {
      std::print(stream, "    Imposed chemical potential:  {: .6e} [K]\n", Units::EnergyToKJPerMol * imposedChemicalPotential.value());
    }
    std::print(stream, "    ---------------------------------------------------------------------------\n");
    if (imposedFugacity)
    {
      std::print(stream, "    Imposed fugacity:            {: .6e} [Pa]\n", Units::PressureConversionFactor * imposedFugacity.value());
    }
    std::print(stream, "    Measured fugacity:           {: .6e} +/- {: .6e} [Pa]\n",
      Units::PressureConversionFactor * averageFugacityDUDlambda.first,
      Units::PressureConversionFactor * averageFugacityDUDlambda.second);
    
    std::print(stream, "\n\n");
  }

  return stream.str();
}
