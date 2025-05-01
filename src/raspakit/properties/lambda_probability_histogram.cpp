module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <numeric>
#include <optional>
#include <print>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>
#endif

module property_lambda_probability_histogram;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <map>;
import <array>;
import <iostream>;
import <cmath>;
import <string>;
import <sstream>;
import <fstream>;
import <format>;
import <algorithm>;
import <numeric>;
import <cmath>;
import <numbers>;
import <optional>;
import <complex>;
import <print>;
import <ranges>;
import <numeric>;
import <filesystem>;
#endif

import double3;
import archive;
import units;
import stringutils;
import averages;
import json;

void PropertyLambdaProbabilityHistogram::WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase phase,
                                                             bool containsTheFractionalMolecule,
                                                             [[maybe_unused]] double value)
{
  switch (phase)
  {
    case PropertyLambdaProbabilityHistogram::WangLandauPhase::Initialize:
      WangLandauScalingFactor = 0.01;
      std::fill(histogram.begin(), histogram.end(), 0.0);
      break;
    case PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample:
      // assert(currentBin >= 0 && currentBin < numberOfSamplePoints);
      if (containsTheFractionalMolecule)
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
      // std::vector<double>::iterator maxValueIterator = std::max_element(histogram.begin(), histogram.end());
      // double maximumValue = *maxValueIterator;
      //
      // bool allBinsAtLeast90PercentOfThisValue = true;
      // for(size_t i = 0; i < histogram.size(); ++i)
      //{
      //   if (histogram[i] < 0.9 * maximumValue)
      //   {
      //     allBinsAtLeast90PercentOfThisValue = false;
      //     break;
      //   }
      // }
      // if (allBinsAtLeast90PercentOfThisValue)
      //{
      //   std::cout << "allBinsAtLeast90PercentOfThisValue" << std::endl;
      //   WangLandauScalingFactor *= 0.5;
      //   std::fill(histogram.begin(), histogram.end(), 0.0);
      //   occupancyCount = 0;
      //   occupancyTotal = 0;
      // }
    }
    break;
    case PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize:
      std::fill(histogram.begin(), histogram.end(), 0.0);
      break;
  }
}

std::pair<std::vector<double>, std::vector<double>>
PropertyLambdaProbabilityHistogram::normalizedAverageProbabilityHistogram()
{
  std::pair<std::vector<double>, std::vector<double>> histogram_avg = averageProbabilityHistogram();

  double totalHistogram = std::accumulate(histogram_avg.first.begin(), histogram_avg.first.end(), 0.0);
  double normalize = static_cast<double>(numberOfSamplePoints) / totalHistogram;

  for (double &value : histogram_avg.first)
  {
    value *= normalize;
  }
  for (double &value : histogram_avg.second)
  {
    value *= normalize;
  }

  return histogram_avg;
}

std::string PropertyLambdaProbabilityHistogram::writeAveragesStatistics(double beta,
                                                                        std::optional<double> imposedChemicalPotential,
                                                                        std::optional<double> imposedFugacity) const
{
  std::ostringstream stream;

  double conv = Units::EnergyToKelvin;
  size_t lastBin = numberOfSamplePoints - 1;

  std::print(stream, "    Lambda histogram and bias:\n");
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::pair<std::vector<double>, std::vector<double>> histogram_avg = averageProbabilityHistogram();

  double totalHistogram = std::accumulate(histogram_avg.first.begin(), histogram_avg.first.end(), 0.0);
  double normalize = static_cast<double>(numberOfSamplePoints) / totalHistogram;
  for (size_t i = 0; i < numberOfSamplePoints; ++i)
  {
    std::print(stream, "{}{:2d}-{:4f} (lambda) P: {: .5e} +/- {:.5e} bias: {: .5e} [-]\n", "    ", i,
               static_cast<double>(i) * delta, normalize * histogram_avg.first[i], normalize * histogram_avg.second[i],
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

  double minimum_free_energy = 0.0;
  auto minimum_iterator = std::ranges::min_element(freeEnergy.first);
  if (minimum_iterator != freeEnergy.first.end())
  {
    minimum_free_energy = conv * (*minimum_iterator);
  }

  switch (Units::unitSystem)
  {
    case Units::System::RASPA:
    {
      for (size_t i = 0; i < numberOfSamplePoints; ++i)
      {
        std::print(stream, "{}{:2d}-{:4f} (lambda) Free energy: {:.6e} +/- {:.6e} [K]\n", "    ", i,
                   static_cast<double>(i) * delta, conv * freeEnergy.first[i] - minimum_free_energy,
                   conv * freeEnergy.second[i]);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential: (ln(P(lambda=1))-ln(P(lambda=0)))/Beta\n");
      for (size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
      {
        double blockAverage = averagedExcessChemicalPotential(blockIndex, beta);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex,
                   conv * (blockAverage + excessChemicalPotentialBias));
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
        std::print(stream, "    Imposed chemical potential:   {: .6e} [K]\n",
                   Units::EnergyToKelvin * imposedChemicalPotential.value());
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
        std::print(stream, "    Imposed chemical potential:   {: .6e} [kJ/mol]\n",
                   Units::EnergyToKJPerMol * imposedChemicalPotential.value());
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      if (imposedFugacity)
      {
        std::print(stream, "    Imposed fugacity:             {: .6e} [Pa]\n",
                   Units::PressureConversionFactor * imposedFugacity.value());
      }
      std::print(stream, "    Measured fugacity:            {: .6e} +/- {: .6e} [Pa]\n",
                 Units::PressureConversionFactor * measuredFugacity.first,
                 Units::PressureConversionFactor * measuredFugacity.second);
      break;
    }
    case Units::System::ReducedUnits:
    {
      for (size_t i = 0; i < numberOfSamplePoints; ++i)
      {
        std::print(stream, "{}{:2d}-{:4f} (lambda) Free energy: {:.6e} +/- {:.6e} [-]\n", "    ", i,
                   static_cast<double>(i) * delta, beta * (freeEnergy.first[i] - minimum_free_energy),
                   beta * freeEnergy.second[i]);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential: (ln(P(lambda=1))-ln(P(lambda=0)))\n");
      for (size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
      {
        double blockAverage = averagedExcessChemicalPotential(blockIndex, beta);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex,
                   beta * (blockAverage + excessChemicalPotentialBias));
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Beta * Excess chemical potential:    {: .6e} +/- {: .6e} [-]\n",
                 beta * (excessChemicalPotential.first + excessChemicalPotentialBias),
                 beta * excessChemicalPotential.second);
      std::print(stream, "    Beta * Ideal gas chemical potential: {: .6e} +/- {: .6e} [-]\n",
                 beta * idealGasChemicalPotential.first, beta * idealGasChemicalPotential.second);
      std::print(stream, "    Beta * Total chemical potential:     {: .6e} +/- {: .6e} [-]\n",
                 beta * totalChemicalPotential.first, beta * totalChemicalPotential.second);
      if (imposedChemicalPotential)
      {
        std::print(stream, "    Beta * Imposed chemical potential:   {: .6e} [-]\n",
                   beta * imposedChemicalPotential.value());
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      if (imposedFugacity)
      {
        std::print(stream, "    Imposed fugacity:             {: .6e} [{}]\n",
                   Units::PressureConversionFactor * imposedFugacity.value(), Units::unitOfPressureString);
      }
      std::print(stream, "    Measured fugacity:            {: .6e} +/- {: .6e} [{}]\n",
                 Units::PressureConversionFactor * measuredFugacity.first,
                 Units::PressureConversionFactor * measuredFugacity.second, Units::unitOfPressureString);
      break;
    }
  }

  std::print(stream, "\n\n");
  return stream.str();
}

std::string PropertyLambdaProbabilityHistogram::writeDUdLambdaStatistics(double beta,
                                                                         std::optional<double> imposedChemicalPotential,
                                                                         std::optional<double> imposedFugacity) const
{
  std::ostringstream stream;

  if (computeDUdlambda)
  {
    std::print(stream, "    Thermodynamic integration (dU/dlambda)\n");
    std::print(stream, "    ===========================================================================\n\n");

    double conv = Units::EnergyToKelvin;
    std::pair<std::vector<double3>, std::vector<double3>> dudlambda = averageDuDlambda();
    for (size_t binIndex = 0; binIndex < numberOfSamplePoints; ++binIndex)
    {
      std::print(stream, "{}{:2d}-{:5f} (lambda) <dU/dlambda>: {: .6e} +/- {:.6e} [K/-]\n", "    ", binIndex,
                 static_cast<double>(binIndex) * delta, conv * dudlambda.first[binIndex].x,
                 conv * dudlambda.second[binIndex].x);
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
      std::print(stream, "    Imposed chemical potential:  {: .6e} [K]\n",
                 Units::EnergyToKelvin * imposedChemicalPotential.value());
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
      std::print(stream, "    Imposed chemical potential:  {: .6e} [kJ/mol]\n",
                 Units::EnergyToKJPerMol * imposedChemicalPotential.value());
    }
    std::print(stream, "    ---------------------------------------------------------------------------\n");
    if (imposedFugacity)
    {
      std::print(stream, "    Imposed fugacity:            {: .6e} [Pa]\n",
                 Units::PressureConversionFactor * imposedFugacity.value());
    }
    std::print(stream, "    Measured fugacity:           {: .6e} +/- {: .6e} [Pa]\n",
               Units::PressureConversionFactor * averageFugacityDUDlambda.first,
               Units::PressureConversionFactor * averageFugacityDUDlambda.second);

    std::print(stream, "\n\n");
  }

  return stream.str();
}

void PropertyLambdaProbabilityHistogram::writeBiasingFile(std::filesystem::path path)
{
  std::ofstream stream(path, std::ios::out);

  nlohmann::json json;
  std::vector<double> shifted{biasFactor};
  for (size_t i = 0; i < shifted.size(); ++i)
  {
    shifted[i] -= shifted.back();
  }
  json["bias"] = shifted;

  stream << json.dump(2);
}

void PropertyLambdaProbabilityHistogram::readBiasingFile(std::filesystem::path path)
{
  std::ifstream stream(path, std::ios::in);

  if (!std::filesystem::exists(path))
  {
    throw std::runtime_error(std::format("[Input reader]: File '{}' not found\n", path.string()));
  }

  std::ifstream input(path);

  nlohmann::basic_json<nlohmann::raspa_map> parsed_lambda_bias_data{};

  parsed_lambda_bias_data = nlohmann::json::parse(input);

  if (parsed_lambda_bias_data.contains("bias") && parsed_lambda_bias_data["bias"].is_array())
  {
    std::vector<double> lambda_bias_data = parsed_lambda_bias_data["bias"].get<std::vector<double>>();

    // resample bias-factor to take different vector size into account
    for (size_t i = 0; i < biasFactor.size(); ++i)
    {
      double lambda = static_cast<double>(i) * delta;

      size_t index = static_cast<size_t>(lambda * static_cast<double>(lambda_bias_data.size() - 1));
      biasFactor[i] = lambda_bias_data[index];
    }
  }
}

nlohmann::json PropertyLambdaProbabilityHistogram::jsonAveragesStatistics(
    double beta, std::optional<double> imposedChemicalPotential, std::optional<double> imposedFugacity) const
{
  nlohmann::json status;

  size_t lastBin = numberOfSamplePoints - 1;

  /*
  std::pair<std::vector<double>, std::vector<double>> histogram_avg = averageProbabilityHistogram();
  for (size_t i = 0; i < numberOfSamplePoints; ++i)
  {
    std::print(stream, "{}{:2d}-{:4f} (lambda) P: {: .5e} +/- {:.5e} bias: {: .5e} [-]\n", "    ", i,
               static_cast<double>(i) * delta, histogram_avg.first[i], histogram_avg.second[i], biasFactor[i]);
  }
  std::print(stream, "\n\n");
  */

  std::pair<double, double> idealGasChemicalPotential = averageIdealGasChemicalPotential(beta);

  std::pair<std::vector<double>, std::vector<double>> freeEnergy = averageLandauFreeEnergyHistogram(beta);
  std::pair<double, double> excessChemicalPotential = averageExcessChemicalPotential(beta);
  double excessChemicalPotentialBias = (biasFactor[lastBin] - biasFactor[0]) / beta;
  std::pair<double, double> totalChemicalPotential = averageTotalChemicalPotential(beta, excessChemicalPotentialBias);

  std::pair<double, double> measuredFugacity = averageFugacity(beta, excessChemicalPotentialBias);

  /*
  std::print(stream, "    Lambda statistics:\n");
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  for (size_t i = 0; i < numberOfSamplePoints; ++i)
  {
    std::print(stream, "{}{:2d}-{:4f} (lambda) Free energy: {:.6e} +/- {:.6e} [K]\n", "    ", i,
               static_cast<double>(i) * delta, conv * freeEnergy.first[i], conv * freeEnergy.second[i]);
  }
  */
  std::vector<double> excessChemicalPotentialBlocks(numberOfBlocks);
  for (size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
  {
    excessChemicalPotentialBlocks[blockIndex] =
        Units::EnergyToKelvin * averagedExcessChemicalPotential(blockIndex, beta);
  }

  status["excessChemicalPotential"]["blocks"]["[K]"] = excessChemicalPotentialBlocks;
  status["excessChemicalPotential"]["mean"]["[K]"] =
      Units::EnergyToKelvin * (excessChemicalPotential.first + excessChemicalPotentialBias);
  status["excessChemicalPotential"]["confidence"]["[K]"] = Units::EnergyToKelvin * excessChemicalPotential.second;
  status["idealGasChemicalPotential"]["mean"]["[K]"] = Units::EnergyToKelvin * idealGasChemicalPotential.first;
  status["idealGasChemicalPotential"]["confidence"]["[K]"] = Units::EnergyToKelvin * idealGasChemicalPotential.second;
  status["totalGasChemicalPotential"]["mean"]["[K]"] = Units::EnergyToKelvin * totalChemicalPotential.first;
  status["totalGasChemicalPotential"]["confidence"]["[K]"] = Units::EnergyToKelvin * totalChemicalPotential.second;

  status["excessChemicalPotential"]["mean"]["[kJ/mol]"] =
      Units::EnergyToKJPerMol * (excessChemicalPotential.first + excessChemicalPotentialBias);
  status["excessChemicalPotential"]["confidence"]["[kJ/mol]"] =
      Units::EnergyToKJPerMol * excessChemicalPotential.second;
  status["idealGasChemicalPotential"]["mean"]["[kJ/mol]"] = Units::EnergyToKJPerMol * idealGasChemicalPotential.first;
  status["idealGasChemicalPotential"]["confidence"]["[kJ/mol]"] =
      Units::EnergyToKJPerMol * idealGasChemicalPotential.second;
  status["totalGasChemicalPotential"]["mean"]["[kJ/mol]"] = Units::EnergyToKJPerMol * totalChemicalPotential.first;
  status["totalGasChemicalPotential"]["confidence"]["[kJ/mol]"] =
      Units::EnergyToKJPerMol * totalChemicalPotential.second;

  if (imposedChemicalPotential)
  {
    status["imposedChemicalPotential"]["[K]"] = Units::EnergyToKelvin * imposedChemicalPotential.value();
    status["imposedChemicalPotential"]["[kJ/mol]"] = Units::EnergyToKJPerMol * imposedChemicalPotential.value();
  }
  if (imposedFugacity)
  {
    status["imposedFugacity"] = Units::PressureConversionFactor * imposedFugacity.value();
  }
  status["measuredFugacity"]["mean"] = Units::PressureConversionFactor * measuredFugacity.first;
  status["measuredFugacity"]["confidence"] = Units::PressureConversionFactor * measuredFugacity.second;

  return status;
}

nlohmann::json PropertyLambdaProbabilityHistogram::jsonDUdLambdaStatistics(
    double beta, std::optional<double> imposedChemicalPotential, std::optional<double> imposedFugacity) const
{
  nlohmann::json status;

  if (computeDUdlambda)
  {
    std::pair<std::vector<double3>, std::vector<double3>> dudlambda = averageDuDlambda();

    /*
    for (size_t binIndex = 0; binIndex < numberOfSamplePoints; ++binIndex)
    {
      std::print(stream, "{}{:2d}-{:5f} (lambda) <dU/dlambda>: {: .6e} +/- {:.6e} [-]\n", "    ", binIndex,
                 static_cast<double>(binIndex) * delta, conv * dudlambda.first[binIndex].x,
                 conv * dudlambda.second[binIndex].x);
    }
    */
    std::vector<double> excessChemicalPotentialBlocks;
    for (size_t blockIndex = 0; blockIndex < numberOfBlocks; ++blockIndex)
    {
      double blockAverage = averagedExcessChemicalPotentialDUdlambda(blockIndex);
      excessChemicalPotentialBlocks[blockIndex] = Units::EnergyToKelvin * blockAverage;
    }
    std::pair<double, double> averageExcessChemicalPotentialDUDlambda = averageExcessChemicalPotentialDUdlambda();
    std::pair<double, double> averageIdealGasChemicalPotentialDUDlambda = averageIdealGasChemicalPotential(beta);
    std::pair<double, double> averageTotalChemicalPotentialDUDlambda = averageTotalChemicalPotential(beta);
    std::pair<double, double> averageFugacityDUDlambda = averageFugacityDUdlambda(beta);

    status["excessChemicalPotential"]["blocks"] = excessChemicalPotentialBlocks;
    status["excessChemicalPotential"]["mean"]["[K]"] =
        Units::EnergyToKelvin * averageExcessChemicalPotentialDUDlambda.first;
    status["excessChemicalPotential"]["confidence"]["[K]"] =
        Units::EnergyToKelvin * averageExcessChemicalPotentialDUDlambda.second;
    status["idealChemicalPotential"]["mean"]["[K]"] =
        Units::EnergyToKelvin * averageIdealGasChemicalPotentialDUDlambda.first;
    status["idealChemicalPotential"]["confidence"]["[K]"] =
        Units::EnergyToKelvin * averageIdealGasChemicalPotentialDUDlambda.second;
    status["totalChemicalPotential"]["mean"]["[K]"] =
        Units::EnergyToKelvin * averageTotalChemicalPotentialDUDlambda.first;
    status["totalChemicalPotential"]["confidence"]["[K]"] =
        Units::EnergyToKelvin * averageTotalChemicalPotentialDUDlambda.second;

    status["excessChemicalPotential"]["mean"]["[kJ/mol]"] =
        Units::EnergyToKJPerMol * averageExcessChemicalPotentialDUDlambda.first;
    status["excessChemicalPotential"]["confidence"]["[kJ/mol]"] =
        Units::EnergyToKJPerMol * averageExcessChemicalPotentialDUDlambda.second;
    status["idealChemicalPotential"]["mean"]["[kJ/mol]"] =
        Units::EnergyToKJPerMol * averageIdealGasChemicalPotentialDUDlambda.first;
    status["idealChemicalPotential"]["confidence"]["[kJ/mol]"] =
        Units::EnergyToKJPerMol * averageIdealGasChemicalPotentialDUDlambda.second;
    status["totalChemicalPotential"]["mean"]["[kJ/mol]"] =
        Units::EnergyToKJPerMol * averageTotalChemicalPotentialDUDlambda.first;
    status["totalChemicalPotential"]["confidence"]["[kJ/mol]"] =
        Units::EnergyToKJPerMol * averageTotalChemicalPotentialDUDlambda.second;

    if (imposedChemicalPotential)
    {
      status["imposedChemicalPotential"]["[K]"] = Units::EnergyToKelvin * imposedChemicalPotential.value();
      status["imposedChemicalPotential"]["[kJ/mol]"] = Units::EnergyToKJPerMol * imposedChemicalPotential.value();
    }

    if (imposedFugacity)
    {
      status["imposedFugacity"] = Units::PressureConversionFactor * imposedFugacity.value();
    }
    status["measuredFugacity"]["mean"] = Units::PressureConversionFactor * averageFugacityDUDlambda.first;
    status["measuredFugacity"]["confidence"] = Units::PressureConversionFactor * averageFugacityDUDlambda.second;
  }

  return status;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyLambdaProbabilityHistogram &p)
{
  archive << p.versionNumber;
  archive << p.numberOfBlocks;

  archive << p.numberOfSamplePoints;
  archive << p.jump_bins;
  archive << p.currentBin;
  archive << p.delta;

  archive << p.WangLandauScalingFactor;

  archive << p.histogram;
  archive << p.biasFactor;

  archive << p.bookKeepingLambda;

  archive << p.bookKeepingDensity;

  archive << p.computeDUdlambda;
  archive << p.bookKeepingDUdlambda;

  archive << p.occupancyCount;
  archive << p.occupancyTotal;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyLambdaProbabilityHistogram &p)
{
  archive >> p.versionNumber;
  archive >> p.numberOfBlocks;

  archive >> p.numberOfSamplePoints;
  archive >> p.jump_bins;
  archive >> p.currentBin;
  archive >> p.delta;

  archive >> p.WangLandauScalingFactor;

  archive >> p.histogram;
  archive >> p.biasFactor;

  archive >> p.bookKeepingLambda;

  archive >> p.bookKeepingDensity;

  archive >> p.computeDUdlambda;
  archive >> p.bookKeepingDUdlambda;

  archive >> p.occupancyCount;
  archive >> p.occupancyTotal;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyLambdaProbabilityHistogram: Error in binary restart\n"));
  }
#endif

  return archive;
}
