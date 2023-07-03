module;

module loadings;


import <string>;
import <optional>;
import <iostream>;
import <sstream>;
import <vector>;
import <ostream>;

import print;
import component;
import units;

std::string Loadings::printStatus(const Component& comp, std::optional<double> frameworkMass) const
{
  std::ostringstream stream;

  if (frameworkMass.has_value())
  {
    std::print(stream, "Component {} ({})\n", comp.componentId, comp.name);

    const double toMolePerKg = 1000.0 / frameworkMass.value();
    const double toMgPerG = 1000.0 * comp.mass / frameworkMass.value();
    const double densityConversionFactor = 1.0 / (1000.0 * Units::Angstrom  * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant);

    switch(comp.type)
    {
      case Component::Type::Framework:
        std::print(stream, "    framework density:   {: .6e} kg/m^3\n", densityConversionFactor * comp.mass * numberDensities[comp.componentId]);
        break;
      default:
        double loading = numberOfMolecules[comp.componentId];
        std::print(stream, "    absolute adsorption: {: .6e} molec/uc\n", loading);
        std::print(stream, "                         {: .6e} mol/kg-framework\n", loading * toMolePerKg);
        std::print(stream, "                         {: .6e} mg/g-framework\n", loading * toMgPerG);

        double excess_loading = numberOfMolecules[comp.componentId] - comp.amountOfExcessMolecules;
        std::print(stream, "    excess adsorption:   {: .6e} molec/uc\n", excess_loading);
        std::print(stream, "                         {: .6e} mol/kg-framework\n", excess_loading * toMolePerKg);
        std::print(stream, "                         {: .6e} mg/g-framework\n", excess_loading * toMgPerG);
        break;
    }
  }
  else
  {
    const double densityConversionFactor = 1.0 / (1000.0 * Units::Angstrom  * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant);
    std::print(stream, "Component {} ({})\n", comp.componentId, comp.name);
    std::print(stream, "    molecules:        {: .6e} molec/uc\n", numberOfMolecules[comp.componentId]);
    std::print(stream, "    number density:   {: .6e} molec/A^3\n", numberDensities[comp.componentId]);
    std::print(stream, "    density:          {: .6e} kg/m^3\n", densityConversionFactor * comp.mass * numberDensities[comp.componentId]);
  }

  return stream.str();
}

std::string Loadings::printStatus(const Component& comp, const Loadings& average, const Loadings& error, std::optional<double> frameworkMass) const
{
  std::ostringstream stream;

  if (frameworkMass)
  {
    const double toMolePerKg = 1000.0 / frameworkMass.value();
    const double toMgPerKg = 1000.0 * comp.mass / frameworkMass.value();
    const double densityConversionFactor = 1.0 / (1000.0 * Units::Angstrom  * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant);

    std::print(stream, "Component {} ({})\n", comp.componentId, comp.name);

    switch(comp.type)
    {
      case Component::Type::Framework:
        std::print(stream, "    density:   {: .6e} kg/m^3\n", densityConversionFactor * comp.mass * numberDensities[comp.componentId]);
        break;
      default:
        double loading = numberOfMolecules[comp.componentId];
        double loading_avg = average.numberOfMolecules[comp.componentId];
        double loading_error = error.numberOfMolecules[comp.componentId];
        std::print(stream, "    absolute adsorption: {:.6e} molec/uc ({:.6e} +/- {:.6e})]\n",
            loading, loading_avg, loading_error);
        std::print(stream, "                         {:.6e} mol/kg   ({:.6e} +/- {:.6e})]\n",
            loading * toMolePerKg, loading_avg * toMolePerKg, loading_error * toMolePerKg);
        std::print(stream, "                         {:.6e} mg/g     ({:.6e} +/- {:.6e})]\n",
            loading * toMgPerKg, loading_avg * toMgPerKg, loading_error * toMgPerKg);

        double excess_loading = numberOfMolecules[comp.componentId] - comp.amountOfExcessMolecules;
        double excess_loading_avg = average.numberOfMolecules[comp.componentId] - comp.amountOfExcessMolecules;
        double excess_loading_error = error.numberOfMolecules[comp.componentId];
        std::print(stream, "    excess adsorption:   {:.6e} molec/uc ({:.6e} +/- {:.6e})]\n", 
                excess_loading, excess_loading_avg, excess_loading_error);
        std::print(stream, "                         {:.6e} mol/kg   ({:.6e} +/- {:.6e})]\n",
                excess_loading * toMolePerKg, excess_loading_avg * toMolePerKg, excess_loading_error * toMolePerKg);
        std::print(stream, "                         {:.6e} mg/g     ({:.6e} +/- {:.6e})]\n",
                excess_loading * toMgPerKg, excess_loading_avg * toMgPerKg, excess_loading_error * toMgPerKg);
      break;
    }
  }
  else
  {
    const double densityConversionFactor = 1.0 / (1000.0 * Units::Angstrom  * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant);

    std::print(stream, "Component {} ({})\n", comp.componentId, comp.name);
    std::print(stream, "    molecules:      {:.6e} molec/uc   ({:.6e} +/- {:.6e})\n",
        numberOfMolecules[comp.componentId], average.numberOfMolecules[comp.componentId], error.numberOfMolecules[comp.componentId]);
    std::print(stream, "    number density: {:.6e} molec./A^3 ({:.6e} +/- {:.6e})\n",
        numberDensities[comp.componentId], average.numberDensities[comp.componentId], error.numberDensities[comp.componentId]);
    std::print(stream, "    density:        {:.6e} kg/m^3     ({:.6e} +/- {:.6e})\n",
        densityConversionFactor * comp.mass * numberDensities[comp.componentId], densityConversionFactor * comp.mass * average.numberDensities[comp.componentId],
        densityConversionFactor * comp.mass * error.numberDensities[comp.componentId]);
  }

  return stream.str();
}

