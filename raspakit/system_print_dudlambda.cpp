module;

module system;

import <string>;
import <sstream>;

import print;

std::string System::writedudLambdaStatistics() const
{
  std::ostringstream stream;

  std::print(stream, "dU/dlambda statistics:\n");
  std::print(stream, "---------------------------------------------------------------------------\n");


  return stream.str();
}




/*
      std::print(stream, "    dU/dlambda statistics:\n");
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::pair<std::vector<double3>, std::vector<double3>> dudlambda = component.lambda.dUdlambdaBookKeeping.averageDuDlambda();
      for (size_t binIndex = 0; binIndex < component.lambda.dUdlambdaBookKeeping.numberOfBins; ++binIndex)
      {
        std::print(stream, "{}{:2d}-{:5f} (lambda) <dU/dlambda>: {: .6e} +/- {:.6e} [-]\n",
            "    ",binIndex, static_cast<double>(binIndex) *  component.lambda.delta, conv * dudlambda.first[binIndex].x, conv * dudlambda.second[binIndex].x);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential: integral du/dlambda over lambda (trapezoidal rule)\n");
      for (size_t blockIndex = 0; blockIndex < component.lambda.dUdlambdaBookKeeping.numberOfBlocks; ++blockIndex)
      {
        double blockAverage = component.lambda.dUdlambdaBookKeeping.averagedExcessChemicalPotential(blockIndex);
        std::print(stream, "        Block[ {:2d}] {}\n", blockIndex, Units::EnergyToKelvin * blockAverage);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::pair<double, double> averageExcessChemicalPotentialDUDlambda = component.lambda.dUdlambdaBookKeeping.averageExcessChemicalPotential();
      std::pair<double, double> averageIdealGasChemicalPotentialDUDlambda = component.lambda.dUdlambdaBookKeeping.averageIdealGasChemicalPotential(beta);
      std::pair<double, double> averageTotalChemicalPotentialDUDlambda = component.lambda.dUdlambdaBookKeeping.averageTotalChemicalPotential(beta);
      std::pair<double, double> averageFugacityDUDlambda = component.lambda.dUdlambdaBookKeeping.averageFugacity(beta);
      std::print(stream, "    Excess chemical potential:   {: .6e} +/- {: .6e} [K]\n",
              Units::EnergyToKelvin * averageExcessChemicalPotentialDUDlambda.first,
              Units::EnergyToKelvin * averageExcessChemicalPotentialDUDlambda.second);
      std::print(stream, "    Ideal chemical potential:    {: .6e} +/- {: .6e} [K]\n",
              Units::EnergyToKelvin *  averageIdealGasChemicalPotentialDUDlambda.first,
              Units::EnergyToKelvin *  averageIdealGasChemicalPotentialDUDlambda.second);
      std::print(stream, "    Total chemical potential:    {: .6e} +/- {: .6e} [K]\n",
              Units::EnergyToKelvin *  averageTotalChemicalPotentialDUDlambda.first,
              Units::EnergyToKelvin *  averageTotalChemicalPotentialDUDlambda.second);
      if(component.swapable)
      {
        std::print(stream, "    Imposed chemical potential:   {: .6e} [K]\n", Units::EnergyToKelvin * imposedChemicalPotential);
      }
      std::print(stream, "    ---------------------------------------------------------------------------\n");
      std::print(stream, "    Excess chemical potential:   {: .6e} +/- {: .6e} [kJ/mol]\n",
              Units::EnergyToKJPerMol *  averageExcessChemicalPotentialDUDlambda.first,
              Units::EnergyToKJPerMol *  averageExcessChemicalPotentialDUDlambda.second);
      std::print(stream, "    Ideal chemical potential:    {: .6e} +/- {: .6e} [kJ/mol]\n",
              Units::EnergyToKJPerMol *  averageIdealGasChemicalPotentialDUDlambda.first,
              Units::EnergyToKJPerMol *  averageIdealGasChemicalPotentialDUDlambda.second);
      std::print(stream, "    Total chemical potential:    {: .6e} +/- {: .6e} [kJ/mol]\n",
              Units::EnergyToKJPerMol *  averageTotalChemicalPotentialDUDlambda.first,
              Units::EnergyToKJPerMol *  averageTotalChemicalPotentialDUDlambda.second);
      if(component.swapable)
      {
        std::print(stream, "    Imposed chemical potential:   {: .6e} [kJ/mol]\n", Units::EnergyToKJPerMol * imposedChemicalPotential);
        std::print(stream, "    ---------------------------------------------------------------------------\n");
        std::print(stream, "    Imposed fugacity:             {: .6e} [Pa]\n", Units::PressureConversionFactor * imposedFugacity);
        std::print(stream, "    Measured fugacity:            {: .6e} +/- {: .6e} [Pa]\n",
                Units::PressureConversionFactor * averageFugacityDUDlambda.first,
                Units::PressureConversionFactor * averageFugacityDUDlambda.second);
      }*/