module;

module equation_of_states;

import std;

import archive;
import units;
import cubic;
import component;
import simulationbox;

EquationOfState::EquationOfState(EquationOfState::Type type, EquationOfState::MultiComponentMixingRules rules,
                                 double temperature, double pressure, const SimulationBox &simulationBox,
                                 double heliumVoidFraction, std::vector<Component> &components)
    : equationOfState(type), multiComponentMixingRules(rules)
{
  computeComponentFluidProperties(equationOfState, multiComponentMixingRules, temperature, pressure, simulationBox,
                                  heliumVoidFraction, components);
}

// T in Kelvin
// p in Pascal
void EquationOfState::computeComponentFluidProperties(EquationOfState::Type type,
                                                      EquationOfState::MultiComponentMixingRules rules,
                                                      double temperature, double pressure,
                                                      const SimulationBox &simulationBox, double heliumVoidFraction,
                                                      std::vector<Component> &components)
{
  std::vector<double> a(components.size());
  std::vector<double> b(components.size());
  std::vector<double> A(components.size());
  std::vector<double> B(components.size());

  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (components[i].swappable)
    {
      double Tc = components[i].criticalTemperature;
      double Pc = components[i].criticalPressure;
      double w = components[i].acentricFactor;
      double Tr = temperature / Tc;
      double temp = 0.0;
      double kappa, alpha;
      switch (equationOfState)
      {
        case EquationOfState::Type::SoaveRedlichKwong:
          kappa = 0.480 + 1.574 * w - 0.176 * w * w;
          temp = (1.0 + kappa * (1.0 - std::sqrt(Tr)));
          alpha = temp * temp;
          temp = Units::MolarGasConstant * Tc;
          a[i] = 0.42748 * alpha * temp * temp / Pc;
          b[i] = 0.08664 * temp / Pc;
          temp = Units::MolarGasConstant * temperature;
          A[i] = a[i] * pressure / (temp * temp);
          B[i] = b[i] * pressure / temp;
          break;
        case EquationOfState::Type::PengRobinsonGasem:
          kappa = 0.134 + 0.508 * w - 0.0467 * w * w;
          alpha = std::exp((2.0 + 0.836 * Tr) * (1.0 - std::pow(Tr, kappa)));
          temp = Units::MolarGasConstant * Tc;
          a[i] = 0.45724 * alpha * temp * temp / Pc;
          b[i] = 0.07780 * temp / Pc;
          temp = Units::MolarGasConstant * temperature;
          A[i] = a[i] * pressure / (temp * temp);
          B[i] = b[i] * pressure / temp;
          break;
        case EquationOfState::Type::PengRobinson:
          kappa = 0.37464 + 1.54226 * w - 0.26992 * w * w;
          temp = 1.0 + kappa * (1.0 - std::sqrt(Tr));
          alpha = temp * temp;
          temp = Units::MolarGasConstant * Tc;
          a[i] = 0.45724 * alpha * temp * temp / Pc;
          b[i] = 0.07780 * temp / Pc;
          temp = Units::MolarGasConstant * temperature;
          A[i] = a[i] * pressure / (temp * temp);
          B[i] = b[i] * pressure / temp;
          break;
      }
    }
  }

  double Amix = 0.0;
  double Bmix = 0.0;
  std::vector<std::vector<double>> binaryInteractionParameter(components.size(),
                                                              std::vector<double>(components.size()));
  std::vector<std::vector<double>> aij(components.size(), std::vector<double>(components.size()));
  std::vector<std::vector<double>> Aij(components.size(), std::vector<double>(components.size()));
  switch (rules)
  {
    default:
    case MultiComponentMixingRules::VanDerWaals:
      for (std::size_t i = 0; i < components.size(); i++)
        for (std::size_t j = 0; j < components.size(); j++)
        {
          if (components[i].swappable && components[j].swappable)
          {
            aij[i][j] = (1.0 - binaryInteractionParameter[i][j]) * std::sqrt(a[i] * a[j]);
            Aij[i][j] = (1.0 - binaryInteractionParameter[i][j]) * std::sqrt(A[i] * A[j]);
          }
        }
      Amix = 0.0;
      Bmix = 0.0;
      for (std::size_t i = 0; i < components.size(); ++i)
      {
        if (components[i].swappable)
        {
          Bmix += components[i].molFraction * b[i];
          for (std::size_t j = 0; j < components.size(); ++j)
          {
            if (components[j].swappable)
            {
              Amix += components[i].molFraction * components[j].molFraction * aij[i][j];
            }
          }
        }
      }
      double temp = Units::MolarGasConstant * temperature;
      Amix *= pressure / (temp * temp);
      Bmix *= pressure / temp;
      break;
  }

  std::vector<double> coefficients(4);
  switch (type)
  {
    case EquationOfState::Type::SoaveRedlichKwong:
      coefficients[3] = 1.0;
      coefficients[2] = -1.0;
      coefficients[1] = Amix - Bmix - Bmix * Bmix;
      coefficients[0] = -Amix * Bmix;
      break;
    case EquationOfState::Type::PengRobinsonGasem:
    case EquationOfState::Type::PengRobinson:
    default:
      coefficients[3] = 1.0;
      coefficients[2] = Bmix - 1.0;
      coefficients[1] = Amix - 3.0 * Bmix * Bmix - 2.0 * Bmix;
      coefficients[0] = -(Amix * Bmix - Bmix * Bmix - Bmix * Bmix * Bmix);
      break;
  }

  // solve equation-of-state
  std::vector<double> compressibility(3);
  int numberOfSolutions = 0;
  cubic(coefficients.data(), compressibility.data(), &numberOfSolutions);

  // sort the compressibilities, Compressibility[0] the highest, Compressibility[2] the lowest
  if (compressibility[0] < compressibility[1])
  {
    double temp = compressibility[0];
    compressibility[0] = compressibility[1];
    compressibility[1] = temp;
  }
  if (compressibility[1] < compressibility[2])
  {
    double temp = compressibility[1];
    compressibility[1] = compressibility[2];
    compressibility[2] = temp;
  }
  if (compressibility[0] < compressibility[1])
  {
    double temp = compressibility[0];
    compressibility[0] = compressibility[1];
    compressibility[1] = temp;
  }

  // In cases where three roots are obtained, the highest value corresponds to the solution for vapour phase,
  // and the lowest value corresponds to the solution for liquid phase. (The intermediate solution has no
  // physical meaning). Note that the solutions in this case don't imply that the liquid and vapour phases
  // are in equilibrium with each other. Rather, one of the phases will be stable and the other will be a
  // metastable state. Only at the special condition that T = Tsat and P = Psat are vapour and liquid in
  // equilibrium with each other.

  // In other cases, only a single root is obtained. This may correspond to liquid, vapour or supercritical fluid,
  // depending on the relative magnitudes of P and Pc, and T and Tc.

  std::vector<std::vector<double>> fugacityCoefficients(
      components.size(), std::vector<double>(static_cast<std::size_t>(numberOfSolutions)));
  switch (equationOfState)
  {
    case EquationOfState::Type::SoaveRedlichKwong:
      for (std::size_t i = 0; i < components.size(); i++)
      {
        if (components[i].swappable)
        {
          for (std::size_t j = 0; j < static_cast<std::size_t>(numberOfSolutions); ++j)
          {
            double temp = 0.0;
            for (std::size_t k = 0; k < components.size(); k++) temp += 2.0 * components[k].molFraction * Aij[i][k];

            fugacityCoefficients[i][j] =
                std::exp((B[i] / Bmix) * (compressibility[j] - 1.0) - std::log(compressibility[j] - Bmix) -
                         Amix * (temp / Amix - B[i] / Bmix) * std::log(1.0 + Bmix / compressibility[j]) / Bmix);
          }
        }
      }
      break;
    case EquationOfState::Type::PengRobinsonGasem:
    case EquationOfState::Type::PengRobinson:
    default:
      for (std::size_t i = 0; i < components.size(); ++i)
      {
        if (components[i].swappable)
        {
          for (std::size_t j = 0; j < static_cast<std::size_t>(numberOfSolutions); ++j)
          {
            double temp = 0.0;
            for (std::size_t k = 0; k < components.size(); ++k) temp += 2.0 * components[k].molFraction * Aij[i][k];

            fugacityCoefficients[i][j] =
                std::exp((B[i] / Bmix) * (compressibility[j] - 1.0) - std::log(compressibility[j] - Bmix) -
                         (Amix / (2.0 * std::sqrt(2.0) * Bmix)) * (temp / Amix - B[i] / Bmix) *
                             std::log((compressibility[j] + (1.0 + std::sqrt(2.0)) * Bmix) /
                                      (compressibility[j] + (1.0 - std::sqrt(2)) * Bmix)));
          }
        }
      }
      break;
  }

  /*
  if(ExcessVolume[CurrentSystem]>0.0)
    excess_volume=ExcessVolume[CurrentSystem]*heliumVoidFraction[CurrentSystem];
  else
    excess_volume=Volume[CurrentSystem]*heliumVoidFraction[CurrentSystem];
    */
  double excess_volume = simulationBox.volume * heliumVoidFraction;

  // get the gas-phase fugacity coefficient
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (components[i].swappable)
    {
      if (numberOfSolutions == 1)
      {
        if (!components[i].fugacityCoefficient.has_value())
        {
          components[i].fugacityCoefficient = fugacityCoefficients[i][0];
        }
        components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant * excess_volume *
                                                Units::AngstromCubed * pressure /
                                                (compressibility[0] * Units::MolarGasConstant * temperature);
        components[i].bulkFluidDensity =
            (components[i].molFraction * pressure / (compressibility[0] * Units::MolarGasConstant * temperature)) *
            components[i].totalMass / 1000.0;
        components[i].compressibility = compressibility[0];
        if ((temperature > components[i].criticalTemperature) && (pressure > components[i].criticalPressure))
          fluidState = FluidState::SuperCriticalFluid;
        else if ((temperature < components[i].criticalTemperature) && (pressure < components[i].criticalPressure))
          fluidState = FluidState::Vapor;
        else if ((temperature < components[i].criticalTemperature) && (pressure > components[i].criticalPressure))
          fluidState = FluidState::Liquid;
      }
      else
      {
        if (compressibility[2] > 0.0)
        {
          if (fugacityCoefficients[i][0] < fugacityCoefficients[i][2])
          {
            // vapour (stable) and liquid (metastable)
            fluidState = FluidState::Vapor;
            if (!components[i].fugacityCoefficient.has_value())
            {
              components[i].fugacityCoefficient = fugacityCoefficients[i][0];
            }
            components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant *
                                                    excess_volume * Units::AngstromCubed * pressure /
                                                    (compressibility[0] * Units::MolarGasConstant * temperature);
            components[i].bulkFluidDensity =
                (components[i].molFraction * pressure / (compressibility[0] * Units::MolarGasConstant * temperature)) *
                components[i].totalMass / 1000.0;
            components[i].compressibility = compressibility[0];
          }
          else if (fugacityCoefficients[i][0] > fugacityCoefficients[i][2])
          {
            // vapour (metastable) and liquid (stable)
            fluidState = FluidState::Liquid;
            if (!components[i].fugacityCoefficient.has_value())
            {
              components[i].fugacityCoefficient = fugacityCoefficients[i][2];
            }
            components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant *
                                                    excess_volume * Units::AngstromCubed * pressure /
                                                    (compressibility[2] * Units::MolarGasConstant * temperature);
            components[i].bulkFluidDensity =
                (components[i].molFraction * pressure / (compressibility[2] * Units::MolarGasConstant * temperature)) *
                components[i].totalMass / 1000.0;
            components[i].compressibility = compressibility[2];
          }
          else
          {
            // vapour (stable) and liquid (stable)
            fluidState = FluidState::VaporLiquid;
            if (!components[i].fugacityCoefficient.has_value())
            {
              components[i].fugacityCoefficient = fugacityCoefficients[i][0];
            }
            components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant *
                                                    excess_volume * Units::AngstromCubed * pressure /
                                                    (compressibility[0] * Units::MolarGasConstant * temperature);
            components[i].bulkFluidDensity =
                (components[i].molFraction * pressure / (compressibility[0] * Units::MolarGasConstant * temperature)) *
                components[i].totalMass / 1000.0;
            components[i].compressibility = compressibility[0];
          }
        }
        else
        {
          if (!components[i].fugacityCoefficient.has_value())
          {
            components[i].fugacityCoefficient = fugacityCoefficients[i][0];
          }
          components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant * excess_volume *
                                                  Units::AngstromCubed * pressure /
                                                  (compressibility[0] * Units::MolarGasConstant * temperature);
          components[i].bulkFluidDensity =
              (components[i].molFraction * pressure / (compressibility[0] * Units::MolarGasConstant * temperature)) *
              components[i].totalMass / 1000.0;
          components[i].compressibility = compressibility[0];
          if ((temperature > components[i].criticalTemperature) && (pressure > components[i].criticalPressure))
            fluidState = FluidState::SuperCriticalFluid;
          else if ((temperature < components[i].criticalTemperature) && (pressure < components[i].criticalPressure))
            fluidState = FluidState::Vapor;
          else if ((temperature < components[i].criticalTemperature) && (pressure > components[i].criticalPressure))
            fluidState = FluidState::Liquid;
        }
      }
    }
    else
    {
      components[i].fugacityCoefficient = std::nullopt;
      components[i].partialPressure = -1.0;
    }
  }
}



std::vector<EquationOfState::FluidResult> EquationOfState::computeFluidProperties(
                                                      double temperature, double pressure,
                                                      const std::vector<EquationOfState::FluidInput> &equationOfStateProperties,
                                                      EquationOfState::Type type = EquationOfState::Type::PengRobinson,
                                                      EquationOfState::MultiComponentMixingRules rules = EquationOfState::MultiComponentMixingRules::VanDerWaals)
{
  std::vector<double> a(equationOfStateProperties.size());
  std::vector<double> b(equationOfStateProperties.size());
  std::vector<double> A(equationOfStateProperties.size());
  std::vector<double> B(equationOfStateProperties.size());

  std::vector<FluidResult> result(equationOfStateProperties.size());

  for (std::size_t i = 0; i < equationOfStateProperties.size(); ++i)
  {
    double Tc = equationOfStateProperties[i].criticalTemperature;
    double Pc = equationOfStateProperties[i].criticalPressure;
    double w = equationOfStateProperties[i].acentricFactor;

    if (equationOfStateProperties[i].swappable)
    {
      double Tr = temperature / Tc;
      double temp = 0.0;
      double kappa, alpha;
      switch (type)
      {
        case EquationOfState::Type::SoaveRedlichKwong:
          kappa = 0.480 + 1.574 * w - 0.176 * w * w;
          temp = (1.0 + kappa * (1.0 - std::sqrt(Tr)));
          alpha = temp * temp;
          temp = Units::MolarGasConstant * Tc;
          a[i] = 0.42748 * alpha * temp * temp / Pc;
          b[i] = 0.08664 * temp / Pc;
          temp = Units::MolarGasConstant * temperature;
          A[i] = a[i] * pressure / (temp * temp);
          B[i] = b[i] * pressure / temp;
          break;
        case EquationOfState::Type::PengRobinsonGasem:
          kappa = 0.134 + 0.508 * w - 0.0467 * w * w;
          alpha = std::exp((2.0 + 0.836 * Tr) * (1.0 - std::pow(Tr, kappa)));
          temp = Units::MolarGasConstant * Tc;
          a[i] = 0.45724 * alpha * temp * temp / Pc;
          b[i] = 0.07780 * temp / Pc;
          temp = Units::MolarGasConstant * temperature;
          A[i] = a[i] * pressure / (temp * temp);
          B[i] = b[i] * pressure / temp;
          break;
        case EquationOfState::Type::PengRobinson:
          kappa = 0.37464 + 1.54226 * w - 0.26992 * w * w;
          temp = 1.0 + kappa * (1.0 - std::sqrt(Tr));
          alpha = temp * temp;
          temp = Units::MolarGasConstant * Tc;
          a[i] = 0.45724 * alpha * temp * temp / Pc;
          b[i] = 0.07780 * temp / Pc;
          temp = Units::MolarGasConstant * temperature;
          A[i] = a[i] * pressure / (temp * temp);
          B[i] = b[i] * pressure / temp;
          break;
      }
    }
  }

  double Amix = 0.0;
  double Bmix = 0.0;
  std::vector<std::vector<double>> binaryInteractionParameter(equationOfStateProperties.size(),
                                                              std::vector<double>(equationOfStateProperties.size()));
  std::vector<std::vector<double>> aij(equationOfStateProperties.size(), std::vector<double>(equationOfStateProperties.size()));
  std::vector<std::vector<double>> Aij(equationOfStateProperties.size(), std::vector<double>(equationOfStateProperties.size()));
  switch (rules)
  {
    default:
    case MultiComponentMixingRules::VanDerWaals:
      for (std::size_t i = 0; i < equationOfStateProperties.size(); i++)
      {
        for (std::size_t j = 0; j < equationOfStateProperties.size(); j++)
        {
          if (equationOfStateProperties[i].swappable && equationOfStateProperties[j].swappable)
          {
            aij[i][j] = (1.0 - binaryInteractionParameter[i][j]) * std::sqrt(a[i] * a[j]);
            Aij[i][j] = (1.0 - binaryInteractionParameter[i][j]) * std::sqrt(A[i] * A[j]);
          }
        }
      }
      Amix = 0.0;
      Bmix = 0.0;
      for (std::size_t i = 0; i < equationOfStateProperties.size(); ++i)
      {
        if (equationOfStateProperties[i].swappable)
        {
          Bmix += equationOfStateProperties[i].molFraction * b[i];
          for (std::size_t j = 0; j < equationOfStateProperties.size(); ++j)
          {
            if (equationOfStateProperties[j].swappable)
            {
              Amix += equationOfStateProperties[i].molFraction * equationOfStateProperties[j].molFraction * aij[i][j];
            }
          }
        }
      }
      double temp = Units::MolarGasConstant * temperature;
      Amix *= pressure / (temp * temp);
      Bmix *= pressure / temp;
      break;
  }

  std::vector<double> coefficients(4);
  switch (type)
  {
    case EquationOfState::Type::SoaveRedlichKwong:
      coefficients[3] = 1.0;
      coefficients[2] = -1.0;
      coefficients[1] = Amix - Bmix - Bmix * Bmix;
      coefficients[0] = -Amix * Bmix;
      break;
    case EquationOfState::Type::PengRobinsonGasem:
    case EquationOfState::Type::PengRobinson:
    default:
      coefficients[3] = 1.0;
      coefficients[2] = Bmix - 1.0;
      coefficients[1] = Amix - 3.0 * Bmix * Bmix - 2.0 * Bmix;
      coefficients[0] = -(Amix * Bmix - Bmix * Bmix - Bmix * Bmix * Bmix);
      break;
  }

  // solve equation-of-state
  std::vector<double> compressibility(3);
  int numberOfSolutions = 0;
  cubic(coefficients.data(), compressibility.data(), &numberOfSolutions);

  // sort the compressibilities, Compressibility[0] the highest, Compressibility[2] the lowest
  if (compressibility[0] < compressibility[1])
  {
    double temp = compressibility[0];
    compressibility[0] = compressibility[1];
    compressibility[1] = temp;
  }
  if (compressibility[1] < compressibility[2])
  {
    double temp = compressibility[1];
    compressibility[1] = compressibility[2];
    compressibility[2] = temp;
  }
  if (compressibility[0] < compressibility[1])
  {
    double temp = compressibility[0];
    compressibility[0] = compressibility[1];
    compressibility[1] = temp;
  }

  // In cases where three roots are obtained, the highest value corresponds to the solution for vapour phase,
  // and the lowest value corresponds to the solution for liquid phase. (The intermediate solution has no
  // physical meaning). Note that the solutions in this case don't imply that the liquid and vapour phases
  // are in equilibrium with each other. Rather, one of the phases will be stable and the other will be a
  // metastable state. Only at the special condition that T = Tsat and P = Psat are vapour and liquid in
  // equilibrium with each other.

  // In other cases, only a single root is obtained. This may correspond to liquid, vapour or supercritical fluid,
  // depending on the relative magnitudes of P and Pc, and T and Tc.

  std::vector<std::vector<double>> fugacityCoefficients(
      equationOfStateProperties.size(), std::vector<double>(static_cast<std::size_t>(numberOfSolutions)));
  switch (type)
  {
    case EquationOfState::Type::SoaveRedlichKwong:
      for (std::size_t i = 0; i < equationOfStateProperties.size(); i++)
      {
        if (equationOfStateProperties[i].swappable)
        {
          for (std::size_t j = 0; j < static_cast<std::size_t>(numberOfSolutions); ++j)
          {
            double temp = 0.0;
            for (std::size_t k = 0; k < equationOfStateProperties.size(); k++) 
            {
              temp += 2.0 * equationOfStateProperties[k].molFraction * Aij[i][k];
            }

            fugacityCoefficients[i][j] =
                std::exp((B[i] / Bmix) * (compressibility[j] - 1.0) - std::log(compressibility[j] - Bmix) -
                         Amix * (temp / Amix - B[i] / Bmix) * std::log(1.0 + Bmix / compressibility[j]) / Bmix);
          }
        }
      }
      break;
    case EquationOfState::Type::PengRobinsonGasem:
    case EquationOfState::Type::PengRobinson:
    default:
      for (std::size_t i = 0; i < equationOfStateProperties.size(); ++i)
      {
        if (equationOfStateProperties[i].swappable)
        {
          for (std::size_t j = 0; j < static_cast<std::size_t>(numberOfSolutions); ++j)
          {
            double temp = 0.0;
            for (std::size_t k = 0; k < equationOfStateProperties.size(); ++k) 
            {
              temp += 2.0 * equationOfStateProperties[k].molFraction * Aij[i][k];
            }

            fugacityCoefficients[i][j] =
                std::exp((B[i] / Bmix) * (compressibility[j] - 1.0) - std::log(compressibility[j] - Bmix) -
                         (Amix / (2.0 * std::sqrt(2.0) * Bmix)) * (temp / Amix - B[i] / Bmix) *
                             std::log((compressibility[j] + (1.0 + std::sqrt(2.0)) * Bmix) /
                                      (compressibility[j] + (1.0 - std::sqrt(2)) * Bmix)));
          }
        }
      }
      break;
  }

  /*
  if(ExcessVolume[CurrentSystem]>0.0)
    excess_volume=ExcessVolume[CurrentSystem]*heliumVoidFraction[CurrentSystem];
  else
    excess_volume=Volume[CurrentSystem]*heliumVoidFraction[CurrentSystem];
    */
  //double excess_volume = simulationBox.volume * heliumVoidFraction;

  // get the gas-phase fugacity coefficient
  for (std::size_t i = 0; i < equationOfStateProperties.size(); ++i)
  {
    if (equationOfStateProperties[i].swappable)
    {
      if (numberOfSolutions == 1)
      {
        if (!result[i].fugacityCoefficient.has_value())
        {
          result[i].fugacityCoefficient = fugacityCoefficients[i][0];
        }
        //result[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant * excess_volume *
        //                                        Units::AngstromCubed * pressure /
        //                                        (compressibility[0] * Units::MolarGasConstant * temperature);
        //result[i].bulkFluidDensity =
        //    (mol_fraction_i * pressure / (compressibility[0] * Units::MolarGasConstant * temperature)) *
        //    components[i].totalMass / 1000.0;
        result[i].compressibility = compressibility[0];
        if ((temperature > equationOfStateProperties[i].criticalTemperature) && (pressure > equationOfStateProperties[i].criticalPressure))
          result[i].fluidState = FluidState::SuperCriticalFluid;
        else if ((temperature < equationOfStateProperties[i].criticalTemperature) && (pressure < equationOfStateProperties[i].criticalPressure))
          result[i].fluidState = FluidState::Vapor;
        else if ((temperature < equationOfStateProperties[i].criticalTemperature) && (pressure > equationOfStateProperties[i].criticalPressure))
          result[i].fluidState = FluidState::Liquid;
      }
      else
      {
        if (compressibility[2] > 0.0)
        {
          if (fugacityCoefficients[i][0] < fugacityCoefficients[i][2])
          {
            // vapour (stable) and liquid (metastable)
            result[i].fluidState = FluidState::Vapor;
            if (!result[i].fugacityCoefficient.has_value())
            {
              result[i].fugacityCoefficient = fugacityCoefficients[i][0];
            }
            //components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant *
            //                                        excess_volume * Units::AngstromCubed * pressure /
            //                                        (compressibility[0] * Units::MolarGasConstant * temperature);
            //components[i].bulkFluidDensity =
            //    (components[i].molFraction * pressure / (compressibility[0] * Units::MolarGasConstant * temperature)) *
            //    components[i].totalMass / 1000.0;
            result[i].compressibility = compressibility[0];
          }
          else if (fugacityCoefficients[i][0] > fugacityCoefficients[i][2])
          {
            // vapour (metastable) and liquid (stable)
            result[i].fluidState = FluidState::Liquid;
            if (!result[i].fugacityCoefficient.has_value())
            {
              result[i].fugacityCoefficient = fugacityCoefficients[i][2];
            }
            //components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant *
            //                                        excess_volume * Units::AngstromCubed * pressure /
            //                                        (compressibility[2] * Units::MolarGasConstant * temperature);
            //components[i].bulkFluidDensity =
            //    (components[i].molFraction * pressure / (compressibility[2] * Units::MolarGasConstant * temperature)) *
            //    components[i].totalMass / 1000.0;
            result[i].compressibility = compressibility[2];
          }
          else
          {
            // vapour (stable) and liquid (stable)
            result[i].fluidState = FluidState::VaporLiquid;
            if (!result[i].fugacityCoefficient.has_value())
            {
              result[i].fugacityCoefficient = fugacityCoefficients[i][0];
            }
            //components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant *
            //                                        excess_volume * Units::AngstromCubed * pressure /
            //                                        (compressibility[0] * Units::MolarGasConstant * temperature);
            //components[i].bulkFluidDensity =
            //    (components[i].molFraction * pressure / (compressibility[0] * Units::MolarGasConstant * temperature)) *
            //    components[i].totalMass / 1000.0;
            result[i].compressibility = compressibility[0];
          }
        }
        else
        {
          if (!result[i].fugacityCoefficient.has_value())
          {
            result[i].fugacityCoefficient = fugacityCoefficients[i][0];
          }
          //components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant * excess_volume *
          //                                        Units::AngstromCubed * pressure /
          //                                        (compressibility[0] * Units::MolarGasConstant * temperature);
          //components[i].bulkFluidDensity =
          //    (components[i].molFraction * pressure / (compressibility[0] * Units::MolarGasConstant * temperature)) *
          //    components[i].totalMass / 1000.0;
          result[i].compressibility = compressibility[0];
          if ((temperature > equationOfStateProperties[i].criticalTemperature) && (pressure > equationOfStateProperties[i].criticalPressure))
            result[i].fluidState = FluidState::SuperCriticalFluid;
          else if ((temperature < equationOfStateProperties[i].criticalTemperature) && (pressure < equationOfStateProperties[i].criticalPressure))
            result[i].fluidState = FluidState::Vapor;
          else if ((temperature < equationOfStateProperties[i].criticalTemperature) && (pressure > equationOfStateProperties[i].criticalPressure))
            result[i].fluidState = FluidState::Liquid;
        }
      }
    }
    else
    {
      result[i].fugacityCoefficient = std::nullopt;
    }
  }

  return result;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EquationOfState &s)
{
  archive << s.versionNumber;

  archive << s.fluidState;
  archive << s.equationOfState;
  archive << s.multiComponentMixingRules;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EquationOfState &s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > s.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'EquationOfState' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> s.fluidState;
  archive >> s.equationOfState;
  archive >> s.multiComponentMixingRules;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("EquationOfState: Error in binary restart\n"));
  }
#endif

  return archive;
}
