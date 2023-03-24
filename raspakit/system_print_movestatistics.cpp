module;

module system;


import double3;
import component;
import print;
import lambda;
import simulationbox;
import units;
import loadings;
import averages;
import property_lambda_probability_histogram;
import property_dudlambda;
import property_widom;
import property_loading;
import move_statistics;
import mc_moves_probabilities_system;
import mc_moves_probabilities_particles;

import <iostream>;
import <random>;
import <sstream>;

inline void formatStatistics(std::ostream &stream, const std::string name, const MoveStatistics<double>& move)
{
  std::print(stream, "{} total:        {:10}\n", name, move.counts);
  std::print(stream, "{} constructed:  {:10}\n", name, move.constructed);
  std::print(stream, "{} accepted:     {:10}\n", name, move.accepted);
  std::print(stream, "{} fraction:     {:10f}\n", name, move.accepted / std::max(1.0, double(move.counts)));
  std::print(stream, "{} max-change:   {:10f}\n\n", name, move.maxChange);
}

inline void formatStatistics(std::ostream &stream, const std::string name, const MoveStatistics<double3> &move)
{
  std::print(stream, "{} total:        {:10} {:10} {:10}\n", name, move.counts.x, move.counts.y, move.counts.z);
  std::print(stream, "{} constructed:  {:10} {:10} {:10}\n", name, move.constructed.x, move.constructed.y, move.constructed.z);
  std::print(stream, "{} accepted:     {:10} {:10} {:10}\n", name, move.accepted.x, move.accepted.y, move.accepted.z);
  std::print(stream, "{} fraction:     {:10f} {:10f} {:10f}\n", name, move.accepted.x / std::max(1.0, double(move.counts.x)),
                 move.accepted.y / std::max(1.0, double(move.counts.y)), move.accepted.z / std::max(1.0, double(move.counts.z)));
  std::print(stream, "{} max-change:   {:10f} {:10f} {:10f}\n\n", name, move.maxChange.x, move.maxChange.y, move.maxChange.z);
}

void System::writeMCMoveStatistics(std::ostream &stream) const
{
  if (mc_moves_probabilities.probabilityVolumeMove > 0.0) formatStatistics(stream, "Volume", mc_moves_probabilities.statistics_VolumeMove);
  if (mc_moves_probabilities.probabilityGibbsVolumeMove > 0.0) formatStatistics(stream, "Gibbs Volume", mc_moves_probabilities.statistics_GibbsVolumeMove);

  for (size_t componentId = 0; const Component& component: components)
  {
    std::print(stream,"Component {} [{}]\n", componentId, component.name);

    std::print(stream, component.mc_moves_probabilities.writeMCMoveStatistics());

    double conv = Units::EnergyToKelvin;
    double imposedChemicalPotential = std::log(beta * component.molFraction * pressure) / beta;
    double imposedFugacity = component.molFraction * pressure;

    if(component.hasFractionalMolecule)
    {
      size_t numberOfBins = component.lambda.numberOfBins;
      size_t lastBin = numberOfBins-1;

      std::print(stream, "    Lambda histogram and bias:\n");
      std::print(stream, "    -------------------------------------------------------------------------------\n");
      std::pair<std::vector<double>, std::vector<double>> histogram = component.lambda.newHistogram.averageProbabilityHistogram();
      for (size_t i = 0; i < component.lambda.numberOfBins; ++i)
      {
        std::print(stream, "{}{:2d}-{:5f} (lambda)  P: {: .6e} +/- {:.6e} bias: {: .6e} [-]\n", 
                "    ",i, static_cast<double>(i) *  component.lambda.delta,
                histogram.first[i],
                histogram.second[i],
                component.lambda.biasFactor[i]);
      }
      std::print(stream, "\n\n");

      std::pair<Loadings, Loadings> loadingData = averageLoadings.averageLoading();
      std::pair<double, double> idealGasChemicalPotential = component.lambda.newHistogram.averageIdealGasChemicalPotential(beta);

      std::pair<std::vector<double>, std::vector<double>> freeEnergy = component.lambda.newHistogram.averageLandauFreeEnergyHistogram(beta);
      std::pair<double, double> excessChemicalPotential = component.lambda.newHistogram.averageExcessChemicalPotential(beta);
      double excessChemicalPotentialBias = (component.lambda.biasFactor[lastBin] - component.lambda.biasFactor[0]) / beta;
      std::pair<double, double> totalChemicalPotential = component.lambda.newHistogram.averageTotalChemicalPotential(beta, excessChemicalPotentialBias);

      std::pair<double, double> measuredFugacity = component.lambda.newHistogram.averageFugacity(beta, excessChemicalPotentialBias);

      std::print(stream, "    Lambda statistics:\n");
      std::print(stream, "    -------------------------------------------------------------------------------\n");
      for (size_t i = 0; i < component.lambda.numberOfBins; ++i)
      {
        std::print(stream, "{}{:2d}-{:5f} (lambda) Free energy: {: .6e} +/- {:.6e} [K]\n",
            "    ",i, static_cast<double>(i) *  component.lambda.delta, 
            conv * freeEnergy.first[i], conv * freeEnergy.second[i]);
      }
      std::print(stream, "    -------------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential: (ln(P(lambda=1))-ln(P(lambda=0)))/Beta\n");
      for (size_t blockIndex = 0; blockIndex < component.lambda.newHistogram.numberOfBlocks; ++blockIndex)
      {
        double blockAverage = component.lambda.newHistogram.averagedExcessChemicalPotential(blockIndex, beta);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, conv * (blockAverage + excessChemicalPotentialBias));
      }
      std::print(stream, "    -------------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential:    {} +/- {} [K]\n", 
              Units::EnergyToKelvin * (excessChemicalPotential.first + excessChemicalPotentialBias),
              Units::EnergyToKelvin * excessChemicalPotential.second);
      std::print(stream, "    Ideal gas chemical potential: {} +/- {} [K]\n", 
              Units::EnergyToKelvin * idealGasChemicalPotential.first, 
              Units::EnergyToKelvin * idealGasChemicalPotential.second);
      std::print(stream, "    Total chemical potential:     {} +/- {} [K]\n", 
              Units::EnergyToKelvin * totalChemicalPotential.first,
              Units::EnergyToKelvin * totalChemicalPotential.second);
      if(component.swapable)
      {
        std::print(stream, "    Imposed chemical potential:   {} [K]\n", Units::EnergyToKelvin * imposedChemicalPotential);
      }
      std::print(stream, "    -------------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential:    {} +/- {} [kJ/mol]\n", 
               Units::EnergyToKJPerMol * (excessChemicalPotential.first + excessChemicalPotentialBias),
               Units::EnergyToKJPerMol * excessChemicalPotential.second);
      std::print(stream, "    Ideal gas chemical potential: {} +/- {} [kJ/mol]\n", 
              Units::EnergyToKJPerMol * idealGasChemicalPotential.first, 
              Units::EnergyToKJPerMol * idealGasChemicalPotential.second);
      std::print(stream, "    Total chemical potential:     {} +/- {} [kJ/mol]\n", 
              Units::EnergyToKJPerMol * totalChemicalPotential.first,
              Units::EnergyToKJPerMol * totalChemicalPotential.second);
      if(component.swapable)
      {
        std::print(stream, "    Imposed chemical potential:   {} [K]\n", Units::EnergyToKJPerMol * imposedChemicalPotential);
        std::print(stream, "    -------------------------------------------------------------------------------\n");
        std::print(stream, "    Imposed fugacity:             {} [Pa]\n", Units::PressureConversionFactor * imposedFugacity);
        std::print(stream, "    Measured fugacity:            {} +/- {} [Pa]\n", 
                Units::PressureConversionFactor * measuredFugacity.first,
                Units::PressureConversionFactor * measuredFugacity.second);
      }
      std::print(stream, "\n\n");

      std::print(stream, "    dU/dlambda statistics:\n");
      std::print(stream, "    -------------------------------------------------------------------------------\n");
      std::pair<std::vector<double3>, std::vector<double3>> dudlambda = component.lambda.dUdlambdaBookKeeping.averageDuDlambda();
      for (size_t binIndex = 0; binIndex < component.lambda.dUdlambdaBookKeeping.numberOfBins; ++binIndex)
      {
        std::print(stream, "{}{:2d}-{:5f} (lambda) <dU/dlambda>: {: .6e} +/- {:.6e} [-]\n",
            "    ",binIndex, static_cast<double>(binIndex) *  component.lambda.delta, conv * dudlambda.first[binIndex].x, conv * dudlambda.second[binIndex].x);
      }
      std::print(stream, "    -----------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential: integral du/dlambda over lambda (trapezoidal rule)\n");
      for (size_t blockIndex = 0; blockIndex < component.lambda.dUdlambdaBookKeeping.numberOfBlocks; ++blockIndex)
      {
        double blockAverage = component.lambda.dUdlambdaBookKeeping.averagedExcessChemicalPotential(blockIndex);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, Units::EnergyToKelvin * blockAverage);
      }
      std::print(stream, "    -------------------------------------------------------------------------------\n");
      std::pair<double, double> averageExcessChemicalPotentialDUDlambda = component.lambda.dUdlambdaBookKeeping.averageExcessChemicalPotential();
      std::pair<double, double> averageIdealGasChemicalPotentialDUDlambda = component.lambda.dUdlambdaBookKeeping.averageIdealGasChemicalPotential(beta);
      std::pair<double, double> averageTotalChemicalPotentialDUDlambda = component.lambda.dUdlambdaBookKeeping.averageTotalChemicalPotential(beta);
      std::pair<double, double> averageFugacityDUDlambda = component.lambda.dUdlambdaBookKeeping.averageFugacity(beta);
      std::print(stream, "    Excess chemical potential:   {} +/- {} [K]\n", 
              Units::EnergyToKelvin * averageExcessChemicalPotentialDUDlambda.first, 
              Units::EnergyToKelvin * averageExcessChemicalPotentialDUDlambda.second);
      std::print(stream, "    Ideal chemical potential:    {} +/- {} [K]\n", 
              Units::EnergyToKelvin *  averageIdealGasChemicalPotentialDUDlambda.first,
              Units::EnergyToKelvin *  averageIdealGasChemicalPotentialDUDlambda.second);
      std::print(stream, "    Total chemical potential:    {} +/- {} [K]\n", 
              Units::EnergyToKelvin *  averageTotalChemicalPotentialDUDlambda.first,
              Units::EnergyToKelvin *  averageTotalChemicalPotentialDUDlambda.second);
      if(component.swapable)
      {
        std::print(stream, "    Imposed chemical potential:   {} [K]\n", Units::EnergyToKelvin * imposedChemicalPotential);
      }
      std::print(stream, "    -----------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential:   {} +/- {} [kJ/mol]\n", 
              Units::EnergyToKJPerMol *  averageExcessChemicalPotentialDUDlambda.first,
              Units::EnergyToKJPerMol *  averageExcessChemicalPotentialDUDlambda.second);
      std::print(stream, "    Ideal chemical potential:    {} +/- {} [kJ/mol]\n", 
              Units::EnergyToKJPerMol *  averageIdealGasChemicalPotentialDUDlambda.first,
              Units::EnergyToKJPerMol *  averageIdealGasChemicalPotentialDUDlambda.second);
      std::print(stream, "    Total chemical potential:    {} +/- {} [kJ/mol]\n", 
              Units::EnergyToKJPerMol *  averageTotalChemicalPotentialDUDlambda.first,
              Units::EnergyToKJPerMol *  averageTotalChemicalPotentialDUDlambda.second);
      if(component.swapable)
      {
        std::print(stream, "    Imposed chemical potential:   {} [kJ/mol]\n", Units::EnergyToKJPerMol * imposedChemicalPotential);
        std::print(stream, "    -------------------------------------------------------------------------------\n");
        std::print(stream, "    Imposed fugacity:             {} [Pa]\n", Units::PressureConversionFactor * imposedFugacity);
        std::print(stream, "    Measured fugacity:            {} +/- {} [Pa]\n", 
                Units::PressureConversionFactor * averageFugacityDUDlambda.first,
                Units::PressureConversionFactor * averageFugacityDUDlambda.second);
      }
      std::print(stream, "\n\n");
    }

    if(component.mc_moves_probabilities.probabilityWidomMove > 0.0)
    {
      std::print(stream, "    Widom insertion Rosenbluth weight  statistics:\n");
      std::print(stream, "    -------------------------------------------------------------------------------\n");
      for (size_t blockIndex = 0; blockIndex < component.averageRosenbluthWeights.numberOfBlocks; ++blockIndex)
      {
        double blockAverage = component.averageRosenbluthWeights.averagedRosenbluthWeight(blockIndex);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, blockAverage);
      }
      std::print(stream, "    -----------------------------------------------------------------------\n");
      std::pair<double, double> averageRosenbluthWeight = component.averageRosenbluthWeights.averageRosenbluthWeight();
      std::print(stream, "    Average Rosenbluth weight:   {} +/- {} [-]\n", 
              averageRosenbluthWeight.first, averageRosenbluthWeight.second);
      std::print(stream, "\n\n");


      std::print(stream, "    Widom insertion chemical potential  statistics:\n");
      std::print(stream, "    -------------------------------------------------------------------------------\n");
      for (size_t blockIndex = 0; blockIndex < component.averageRosenbluthWeights.numberOfBlocks; ++blockIndex)
      {
        double blockAverage = component.averageRosenbluthWeights.averagedExcessChemicalPotential(blockIndex, beta);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, conv * blockAverage);
      }
      std::print(stream, "    -----------------------------------------------------------------------\n");
      std::pair<double, double> averageExcessWidomChemicalPotential = component.averageRosenbluthWeights.averageExcessChemicalPotential(beta);
      std::pair<double, double> averageIdealGasWidomChemicalPotential = component.averageRosenbluthWeights.averageIdealGasChemicalPotential(beta);
      std::pair<double, double> averageTotalWidomChemicalPotential = component.averageRosenbluthWeights.averageTotalChemicalPotential(beta);
      std::pair<double, double> averageWidomFugacity = component.averageRosenbluthWeights.averageFugacity(beta);
      std::print(stream, "    Excess chemical potential:   {} +/- {} [K]\n", 
              Units::EnergyToKelvin *  averageExcessWidomChemicalPotential.first,
              Units::EnergyToKelvin *  averageExcessWidomChemicalPotential.second);
      std::print(stream, "    Ideal chemical potential:    {} +/- {} [K]\n", 
              Units::EnergyToKelvin *  averageIdealGasWidomChemicalPotential.first,
              Units::EnergyToKelvin *  averageIdealGasWidomChemicalPotential.second);
      std::print(stream, "    Total chemical potential:    {} +/- {} [K]\n", 
              Units::EnergyToKelvin *  averageTotalWidomChemicalPotential.first,
              Units::EnergyToKelvin *  averageTotalWidomChemicalPotential.second);
      if(component.swapable)
      {
        std::print(stream, "    Imposed chemical potential:   {} [K]\n", Units::EnergyToKelvin * imposedChemicalPotential);
      }
      std::print(stream, "    -----------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential:   {} +/- {} [kJ/mol]\n", 
              Units::EnergyToKJPerMol *  averageExcessWidomChemicalPotential.first,
              Units::EnergyToKJPerMol *  averageExcessWidomChemicalPotential.second);
      std::print(stream, "    Ideal chemical potential:    {} +/- {} [kJ/mol]\n", 
              Units::EnergyToKJPerMol *  averageIdealGasWidomChemicalPotential.first,
              Units::EnergyToKJPerMol *  averageIdealGasWidomChemicalPotential.second);
      std::print(stream, "    Total chemical potential:    {} +/- {} [kJ/mol]\n", 
              Units::EnergyToKJPerMol *  averageTotalWidomChemicalPotential.first,
              Units::EnergyToKJPerMol *  averageTotalWidomChemicalPotential.second);
      if(component.swapable)
      {
        std::print(stream, "    Imposed chemical potential:   {} [kJ/mol]\n", Units::EnergyToKJPerMol * imposedChemicalPotential);
        std::print(stream, "    -------------------------------------------------------------------------------\n");
        std::print(stream, "    Imposed fugacity:             {} [Pa]\n", Units::PressureConversionFactor * imposedFugacity);
        std::print(stream, "    Measured fugacity:            {} +/- {} [Pa]\n", 
                Units::PressureConversionFactor * averageWidomFugacity.first,
                Units::PressureConversionFactor * averageWidomFugacity.second);
      }
    }
    std::print(stream, "\n\n");

    ++componentId;
  }

 
  std::print(stream, "\n\n");
}
