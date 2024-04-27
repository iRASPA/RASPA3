module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <cmath>
#include <optional>
#include <iostream>
#endif

module equation_of_states;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <cmath>;
import <optional>;
import <iostream>;
#endif

import units;
import cubic;
import component;
import simulationbox;

EquationOfState::EquationOfState(EquationOfState::Type type,
                                 EquationOfState::MultiComponentMixingRules rules,
                                 double temperature,
                                 double pressure,
                                 const SimulationBox &simulationBox,
                                 double HeliumVoidFraction,
                                 std::vector<Component> &components):
  equationOfState(type),
  multiComponentMixingRules(rules)
{
  computeComponentFluidProperties(equationOfState, multiComponentMixingRules, temperature, pressure,
                                  simulationBox, HeliumVoidFraction, components);
}

// T in Kelvin
// p in Pascal
void EquationOfState::computeComponentFluidProperties(EquationOfState::Type type,
                                          EquationOfState::MultiComponentMixingRules rules, 
                                          double temperature,
                                          double pressure,
                                          const SimulationBox &simulationBox,
                                          double HeliumVoidFraction,
                                          std::vector<Component> &components)
{
  std::vector<double> a(components.size());
  std::vector<double> b(components.size());
  std::vector<double> A(components.size());
  std::vector<double> B(components.size());

  for (size_t i = 0; i < components.size(); ++i)
  {
    if (components[i].swapable)
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
        alpha = exp((2.0 + 0.836 * Tr) * (1.0 - std::pow(Tr, kappa)));
        temp = Units::MolarGasConstant * Tc;
        a[i] = 0.45724 * alpha * temp * temp / Pc;
        b[i] = 0.07780 * temp / Pc;
        temp = Units::MolarGasConstant * temperature;
        A[i] = a[i] * pressure / (temp * temp);
        B[i] = b[i] * pressure  / temp;
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
  std::vector<std::vector<double>> binaryInteractionParameter(components.size(), std::vector<double>(components.size()));
  std::vector<std::vector<double>> aij(components.size(), std::vector<double>(components.size()));
  std::vector<std::vector<double>> Aij(components.size(), std::vector<double>(components.size()));
  switch (rules)
  {
  default:
  case MultiComponentMixingRules::VanDerWaals:
    for (size_t i = 0; i < components.size(); i++)
      for (size_t j = 0; j < components.size(); j++)
      {
        if (components[i].swapable && components[j].swapable)
        {
          aij[i][j] = (1.0 - binaryInteractionParameter[i][j]) * std::sqrt(a[i] * a[j]);
          Aij[i][j] = (1.0 - binaryInteractionParameter[i][j]) * std::sqrt(A[i] * A[j]);
        }
      }
    Amix = 0.0;
    Bmix = 0.0;
    for (size_t i = 0; i < components.size(); ++i)
    {
      if (components[i].swapable)
      {
        Bmix += components[i].molFraction * b[i];
        for (size_t j = 0; j < components.size(); ++j)
        {
          if (components[j].swapable)
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

  std::vector<std::vector<double>> fugacityCoefficients(components.size(), std::vector<double>(static_cast<size_t>(numberOfSolutions)));
  switch (equationOfState)
  {
  case EquationOfState::Type::SoaveRedlichKwong:
    for (size_t i = 0; i < components.size(); i++)
    {
      if (components[i].swapable)
      {
        for (size_t j = 0; j < static_cast<size_t>(numberOfSolutions); ++j)
        {
          double temp = 0.0;
          for (size_t k = 0; k < components.size(); k++)
            temp += 2.0 * components[k].molFraction * Aij[i][k];

          fugacityCoefficients[i][j] = std::exp((B[i] / Bmix) * (compressibility[j] - 1.0) - std::log(compressibility[j] - Bmix) -
              Amix * (temp / Amix - B[i] / Bmix) * std::log(1.0 + Bmix / compressibility[j]) / Bmix);
        }
      }
    }
    break;
  case EquationOfState::Type::PengRobinsonGasem:
  case EquationOfState::Type::PengRobinson:
  default:
    for (size_t i = 0; i < components.size(); ++i)
    {
      if (components[i].swapable)
      {
        for (size_t j = 0; j < static_cast<size_t>(numberOfSolutions); ++j)
        {
          double temp = 0.0;
          for (size_t k = 0; k < components.size(); ++k)
            temp += 2.0 * components[k].molFraction * Aij[i][k];

          fugacityCoefficients[i][j] = 
            std::exp((B[i] / Bmix) * (compressibility[j] - 1.0) - std::log(compressibility[j] - Bmix)
            - (Amix / (2.0 * std::sqrt(2.0) * Bmix)) * (temp / Amix - B[i] / Bmix) *
            log((compressibility[j] + (1.0 + std::sqrt(2.0)) * Bmix) / 
            (compressibility[j] + (1.0 - std::sqrt(2)) * Bmix)));
        }
      }
    }
    break;
  }


  /*
  if(ExcessVolume[CurrentSystem]>0.0)
    excess_volume=ExcessVolume[CurrentSystem]*HeliumVoidFraction[CurrentSystem];
  else
    excess_volume=Volume[CurrentSystem]*HeliumVoidFraction[CurrentSystem];
    */
  double excess_volume = simulationBox.volume * HeliumVoidFraction;


  // get the gas-phase fugacity coefficient
  for(size_t i=0;i<components.size(); ++i)
  {
    if(components[i].swapable)
    {
      if(numberOfSolutions==1)
      {
        if (!components[i].fugacityCoefficient.has_value())
        {
          components[i].fugacityCoefficient = fugacityCoefficients[i][0];
        }
        components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant *
           excess_volume * Units::AngstromCubed * pressure / (compressibility[0] * Units::MolarGasConstant * temperature);
        components[i].bulkFluidDensity = (components[i].molFraction * pressure / (compressibility[0] * Units::MolarGasConstant * temperature))*
                                                       components[i].totalMass / 1000.0;
        components[i].compressibility = compressibility[0];
        if((temperature > components[i].criticalTemperature) && (pressure > components[i].criticalPressure))
          fluidState = FluidState::SuperCriticalFluid;
        else if((temperature < components[i].criticalTemperature) && (pressure < components[i].criticalPressure))
          fluidState = FluidState::Vapor;
        else if((temperature < components[i].criticalTemperature) && (pressure > components[i].criticalPressure))
          fluidState = FluidState::Liquid;
      }
      else
      {
        if(compressibility[2]>0.0)
        {
          if(fugacityCoefficients[i][0] < fugacityCoefficients[i][2])
          {
            // vapour (stable) and liquid (metastable)
            fluidState = FluidState::Vapor;
            if (!components[i].fugacityCoefficient.has_value())
            {
              components[i].fugacityCoefficient = fugacityCoefficients[i][0];
            }
            components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant *
              excess_volume * Units::AngstromCubed * pressure / (compressibility[0] * Units::MolarGasConstant * temperature);
            components[i].bulkFluidDensity = (components[i].molFraction * pressure / (compressibility[0] * Units::MolarGasConstant * temperature))*
                                                          components[i].totalMass/1000.0;
            components[i].compressibility = compressibility[0];
          }
          else if(fugacityCoefficients[i][0] > fugacityCoefficients[i][2])
          {
            // vapour (metastable) and liquid (stable)
            fluidState = FluidState::Liquid;
            if (!components[i].fugacityCoefficient.has_value())
            {
              components[i].fugacityCoefficient = fugacityCoefficients[i][2];
            }
            components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant *
              excess_volume * Units::AngstromCubed * pressure / (compressibility[2] * Units::MolarGasConstant * temperature);
            components[i].bulkFluidDensity = (components[i].molFraction * pressure / (compressibility[2] * Units::MolarGasConstant * temperature))*
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
              excess_volume * Units::AngstromCubed * pressure / (compressibility[0] * Units::MolarGasConstant * temperature);
            components[i].bulkFluidDensity = (components[i].molFraction * pressure / (compressibility[0] * Units::MolarGasConstant * temperature))*
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
          components[i].amountOfExcessMolecules = components[i].molFraction * Units::AvogadroConstant*
            excess_volume * Units::AngstromCubed * pressure / (compressibility[0] * Units::MolarGasConstant * temperature);
          components[i].bulkFluidDensity = (components[i].molFraction * pressure / (compressibility[0] * Units::MolarGasConstant * temperature))*
                                                        components[i].totalMass/1000.0;
          components[i].compressibility = compressibility[0];
          if((temperature > components[i].criticalTemperature) && (pressure > components[i].criticalPressure))
            fluidState = FluidState::SuperCriticalFluid;
          else if((temperature < components[i].criticalTemperature) && (pressure < components[i].criticalPressure))
            fluidState = FluidState::Vapor;
          else if((temperature < components[i].criticalTemperature) && (pressure > components[i].criticalPressure))
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

  /*
  // The term cm^3 (STP) is not a unit of volume but a unit to express the number of gas molecules
  // It is essentially a term in units of pV, e.g. cm^3 at 1 atm and 273 K.
  // The volume of 1 mole of a gas is 22.4 liters

  for(i=0;i<NumberOfComponents;i++)
  {
    // molec/uc -> cm^3 (STP)/g
    Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]=1e6*(MOLAR_GAS_CONSTANT*273.15/ATM_TO_PA)*number_of_unit_cells/FrameworkMass;

    // molec/uc -> cm^3 (STP)/cm^3
    Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]=number_of_unit_cells*(MOLAR_GAS_CONSTANT*273.15/ATM_TO_PA)
                                                           /(Volume[CurrentSystem]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT);

    // mol/kg -> cm^3 (STP)/g
    Components[i].MOL_PER_KG_TO_CC_STP_G[CurrentSystem]=1e3*(MOLAR_GAS_CONSTANT*273.15/ATM_TO_PA);

    // mol/kg -> cm^3 (STP)/cm^3
    Components[i].MOL_PER_KG_TO_CC_STP_CC[CurrentSystem]=(MOLAR_GAS_CONSTANT*273.15/ATM_TO_PA)*FrameworkMass
                                                          /(1e3*Volume[CurrentSystem]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT);
  }

  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].Swapable)
    {
      if(Components[i].PartialPressure[CurrentSystem]<0.0)
      {
        fprintf(stderr, "Error: Partial pressure of component %d [%s] is NOT set !!\n",
          i,Components[i].Name);
        exit(-1);
      }
    }
  }
}

  */

}
