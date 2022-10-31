module;

module mixture_prediction;

import <vector>;
import <span>;
import <cmath>;
import <string>;
import <iostream>;
import <iomanip>;
import <fstream>;
import <limits>;
import <algorithm>;
import <numeric>;
import <ostream>;
import <filesystem>;

import print;
import atom;
import component;
import isotherm;
import multi_site_isotherm;
import system;
import simulationbox;
import pressure_range;

bool LangmuirLoadingSorter(const Component &lhs, const Component &rhs)
{
  if(lhs.isCarrierGas) return false;
  if(rhs.isCarrierGas) return true;
  return lhs.isotherm.parameters[0] < rhs.isotherm.parameters[0];
}

std::vector<size_t> sortedIndices(const std::span<const Component> &v)
{
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) 
       {
         if(v[i1].isCarrierGas) return false;
         if(v[i2].isCarrierGas) return true;
         return v[i1].isotherm.parameters[0] < v[i2].isotherm.parameters[0];
       });

  return idx;
}


MixturePrediction::MixturePrediction(const System &system, double pressureStart, double pressureEnd, 
                   size_t numberOfPressurePoints, MixturePrediction::PressureScale pressureScale):
   system(system),
   displayName(system.components.front().name),
   components(system.spanOfAdsorbateComponents()),
   sortedComponentIndices(sortedIndices(components)),
   sortedComponents(components.begin(),components.end()),
   Ncomp(components.size()),
   predictionMethod(system.mixturePredictionMethod),
   alpha1(Ncomp),
   alpha2(Ncomp),
   alpha_prod(Ncomp),
   x(Ncomp),
   temperature(system.simulationBox.temperature),
   pressureStart(pressureStart),
   pressureEnd(pressureEnd),
   numberOfPressurePoints(numberOfPressurePoints),
   pressureScale(pressureScale)
{
  if(predictionMethod == Isotherm::MixturePredictionMethod::ExplicitLangmuir)
  {
    std::sort(sortedComponents.begin(), sortedComponents.end(), &LangmuirLoadingSorter);
  }
}

MixturePrediction::MixturePrediction(const System &system) :
                                     system(system),
                                     displayName(system.components.front().name),
                                     components(system.spanOfAdsorbateComponents()),
                                     sortedComponentIndices(sortedIndices(components)),
                                     sortedComponents(components.begin(),components.end()),
                                     Ncomp(components.size()),
                                     predictionMethod(system.mixturePredictionMethod),
                                     alpha1(Ncomp),
                                     alpha2(Ncomp),
                                     alpha_prod(Ncomp),
                                     x(Ncomp),
                                     temperature(system.simulationBox.temperature)
                                     //pressureStart(inputreader.pressureStart),
                                     //pressureEnd(inputreader.pressureEnd),
                                     //numberOfPressurePoints(inputreader.numberOfPressurePoints),
                                     //pressureScale(PressureScale(inputreader.pressureScale))
{
  if(predictionMethod == Isotherm::MixturePredictionMethod::ExplicitLangmuir)
  {
    std::sort(sortedComponents.begin(), sortedComponents.end(), &LangmuirLoadingSorter);
  }
}

std::pair<size_t, size_t> MixturePrediction::predictMixture(const std::vector<double> &Yi,
                                                            const double &P,
                                                            std::vector<double> &Xi,
                                                            std::vector<double> &Ni,
                                                            double *cachedP0,
                                                            double &cachedPsi)
{
  switch(predictionMethod)
  {
    case Isotherm::MixturePredictionMethod::IAST:
    default:
      return computeIAST(Yi, P, Xi, Ni, cachedP0, cachedPsi);
    case Isotherm::MixturePredictionMethod::ExplicitLangmuir:
      return computeExplicitLangmuir(Yi, P, Xi, Ni);
  }
}

void MixturePrediction::writeHeader(std::ostream &stream) const
{
  switch(predictionMethod)
  {
    case Isotherm::MixturePredictionMethod::IAST:
    default:
      std::print(stream, "Component data\n");
      std::print(stream, "=======================================================\n");
      for(size_t i = 0; i < Ncomp; ++i)
      {
        components[i].printBreakthroughStatus(stream);
        std::print(stream, "\n");
      }
      break;
    case Isotherm::MixturePredictionMethod::ExplicitLangmuir:
      std::print(stream, "Sorted component data\n");
      std::print(stream, "=======================================================\n");
      for(size_t i = 0; i < Ncomp; ++i)
      {
        sortedComponents[i].get().printBreakthroughStatus(stream);
        std::print(stream, "\n");
      }
      break;
  }
}

// Yi  = gas phase molefraction
// P   = total pressure
// Xi  = adsorbed phase molefraction
// Ni  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeIAST(const std::vector<double> &Yi,
                                                         const double &P,
                                                         std::vector<double> &Xi, 
                                                         std::vector<double> &Ni,
                                                         double *cachedP0,
                                                         double &cachedPsi)
{
  const double tiny = 1.0e-15;

  if(P < 0.0)
  {
    printErrorStatus(0.0, 0.0, P, Yi, cachedP0);
    throw std::runtime_error("Error (IAST): negative total pressure\n");
  }

  double sumYi = 0.0;
  for(size_t i = 0; i < Ncomp; ++i)
  {
    sumYi += Yi[i];
  }
  if(std::abs(sumYi-1.0) > 1e-15)
  {
    printErrorStatus(0.0, sumYi, P, Yi, cachedP0);
    throw std::runtime_error("Error (IAST): sum Yi at IAST start not unity\n");
  }


  // if only an inert component present
  // this happens at the beginning of the simulation when the whole column is filled with the carrier gas
  for(size_t i = 0; i < Ncomp; ++i)
  {
    if(components[i].isCarrierGas)
    {
      if(std::abs(Yi[i] - 1.0) < tiny)
      {
        for(size_t j = 0; j < Ncomp; ++j)
        {
          Xi[i]  = 0.0;
          Ni[i]  = 0.0;
        }

        // do not count it for the IAST statistics
        return std::make_pair(0, 0);
      }
    }
  }

  double initial_psi = 0.0;
  for(size_t i = 0; i < Ncomp; ++i)
  {
    initial_psi += Yi[i] * components[i].isotherm.psiForPressure(P);
  }

  if(initial_psi < tiny)
  {
    // nothing is adsorbing
    for(size_t i = 0; i < Ncomp; ++i)
    {
      Xi[i]  = 0.0;
      Ni[i]  = 0.0;
    }

    // do not count it for the IAST statistics
    return std::make_pair(0, 0);
  }

  // condition 1: same reduced grand potential for all components (done by using a single variable)
  // condition 2: mol-fractions add up to unity

  double psi = 0.0;
  size_t nr_steps=0;
  if(cachedPsi > tiny)
  {
    initial_psi = cachedPsi;
  }
  // for this initial estimate 'initial_psi' compute the sum of mol-fractions
  double sumXi = 0.0;
  for(size_t i = 0; i < Ncomp; ++i)
  {
    sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(initial_psi, cachedP0[i]);
  }

  // initialize the bisection algorithm
  double left_bracket = initial_psi;
  double right_bracket = initial_psi;
  if(sumXi > 1.0)
  {
    do
    {
      right_bracket *= 2.0;

      sumXi = 0.0;
      for(size_t i = 0; i < Ncomp; ++i)
      {
        sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(right_bracket, cachedP0[i]);
      }
      ++nr_steps;
      if(nr_steps>100000)
      {
        std::cout << "Left bracket: " << left_bracket << std::endl;
        std::cout << "Right bracket: " << right_bracket << std::endl;
        printErrorStatus(0.0, sumXi, P, Yi, cachedP0);
        throw std::runtime_error("Error (IAST bisection): initial bracketing (for sum > 1) does NOT converge\n");
      }
    } while(sumXi > 1.0);
  }
  else
  {

    // Make an initial estimate for the reduced grandpotential when the
    // sum of the molefractions is larger than 1
    do 
    {
      left_bracket *= 0.5;

      sumXi = 0.0;
      for(size_t i = 0; i < Ncomp; ++i)
      {
        sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(left_bracket, cachedP0[i]);
      }
      ++nr_steps;
      if(nr_steps>100000)
      {
        std::cout << "Left bracket: " << left_bracket << std::endl;
        std::cout << "Right bracket: " << right_bracket << std::endl;
        printErrorStatus(0.0, sumXi, P, Yi, cachedP0);
        throw std::runtime_error("Error (IAST bisection): initial bracketing (for sum < 1) does NOT converge\n");
      }
    }
    while(sumXi < 1.0);
  }

  // bisection algorithm
  size_t numberOfIASTSteps = 0;
  do 
  {
    psi  = 0.5 * (left_bracket + right_bracket);

    sumXi = 0.0;
    for(size_t i = 0; i < Ncomp; ++i)
    {
      sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(psi, cachedP0[i]);
    }

    if(sumXi > 1.0)
    {
      left_bracket = psi;
    }
    else
    {
      right_bracket = psi;
    }

    ++numberOfIASTSteps;
    if(numberOfIASTSteps>100000)
    {
      throw std::runtime_error("Error (IAST bisection): NO convergence\n");
    }
  }
  while(std::abs(left_bracket - right_bracket) / std::abs(left_bracket + right_bracket) > tiny); // convergence test
 
  psi = 0.5 * (left_bracket + right_bracket);


  sumXi = 0.0;
  for(size_t i = 0; i < Ncomp; ++i)
  {
    sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(psi, cachedP0[i]);
  }

  //printErrorStatus(psi, sumXi, P, Yi, cachedP0);

  // cache the value of psi for subsequent use
  cachedPsi = psi;

  // calculate mol-fractions in adsorbed phase and total loading
  double inverse_q_total = 0.0;
  for(size_t i = 0; i < Ncomp; ++i)
  {
    double ip = components[i].isotherm.inversePressureForPsi(psi, cachedP0[i]);
    Xi[i] = Yi[i] * P * ip / sumXi;

    if(Xi[i] > tiny)
    {
      inverse_q_total += Xi[i] / components[i].isotherm.value(1.0 / ip);
    }
    else
    {
      Xi[i] = 0.0;
    }
  }

  // calculate loading for all of the components
  if (inverse_q_total == 0.0)
  {
    for(size_t i = 0; i < Ncomp; ++i)
    {
      Ni[i] = 0.0;
    }
  }
  else
  {
    for(size_t i = 0; i < Ncomp; ++i)
    {
      Ni[i] = Xi[i] / inverse_q_total;
    }
  }
  return std::make_pair(numberOfIASTSteps, 1);
}

// solve the mixed-langmuir equations derived by Assche et al.
// T. R. Van Assche, G.V. Baron, and J. F. Denayer
// An explicit multicomponent adsorption isotherm model:
// Accounting for the size-effect for components with Langmuir adsorption behavior. 
// Adsorption, 24(6), 517-530 (2018)

// An explicit multicomponent adsorption isotherm model: accounting for the 
// size-effect for components with Langmuir adsorption behavior

// In the input file molecules must be added in the following order:
// Largest molecule should be the first component or the component with 
// smallest saturation(Nimax) loading should be the first component
// Last component is the carrier gas

// At present, only single site isotherms are considered for pure components

std::pair<size_t, size_t> MixturePrediction::computeExplicitLangmuir(const std::vector<double> &Yi,
                                                                     const double &P,
                                                                     std::vector<double> &Xi, 
                                                                     std::vector<double> &Ni)
{
  x[0] = 1.0;
  for(size_t i = 1; i < Ncomp; ++i)
  {
    x[i] = sortedComponents[i].get().isotherm.parameters[0] / sortedComponents[i - 1].get().isotherm.parameters[0];
  }

  size_t index = sortedComponentIndices[Ncomp - 1];
  alpha1[Ncomp - 1] = std::pow((1.0 + sortedComponents[Ncomp - 1].get().isotherm.parameters[1] * Yi[index] * P), x[Ncomp - 1]);
  alpha2[Ncomp - 1] = 1.0 + sortedComponents[Ncomp - 1].get().isotherm.parameters[1] * Yi[index] * P;
  for(size_t i = Ncomp - 2; i > 0; i--)
  {
    index = sortedComponents[i].get().componentId - system.numberOfFrameworks;
    alpha1[i] = std::pow((alpha1[i + 1] + sortedComponents[i].get().isotherm.parameters[1] * Yi[index] * P), x[i]);
    alpha2[i] = alpha1[i + 1] + sortedComponents[i].get().isotherm.parameters[1] * Yi[index] * P;
  }
  index = sortedComponentIndices[0];
  alpha1[0] = alpha1[1] + sortedComponents[0].get().isotherm.parameters[1] * Yi[index] * P;
  alpha2[0] = alpha1[1] + sortedComponents[0].get().isotherm.parameters[1] * Yi[index] * P;

  double beta = alpha2[0];

  alpha_prod[0] = 1.0;
  for(size_t i = 1; i < Ncomp; ++i)
  {
    alpha_prod[i] = (alpha1[i] / alpha2[i]) * alpha_prod[i - 1];
  }

  for(size_t i = 0; i < Ncomp; ++i)
  {
    index = sortedComponentIndices[i];
    Ni[index] = sortedComponents[i].get().isotherm.parameters[0] * sortedComponents[i].get().isotherm.parameters[1] * 
                 Yi[index] * P * alpha_prod[i] / beta;
  }
  double N = 0.0;
  for(size_t i = 0; i < Ncomp; ++i)
  {
    N += Ni[i];
  }
  for(size_t i = 0; i < Ncomp; ++i)
  {
    Xi[i] = Ni[i] / N;
  }

  return std::make_pair(1,1);
}


void MixturePrediction::run(std::ostream &stream)
{
  std::vector<double> Yi(Ncomp);
  std::vector<double> Xi(Ncomp);
  std::vector<double> Ni(Ncomp);
  std::vector<double> cachedP0(Ncomp);
  double cachedPsi;

  for(size_t i = 0; i < Ncomp; ++i)
  {
    Yi[i] = components[i].molFraction;
  }
  
  std::vector<double> pressures(system.pressure_range.pressures());

  std::string directoryNameString = std::print("MixturePrediction/System_{}/", system.systemId);;
  std::filesystem::path directoryName{ directoryNameString };
  std::filesystem::create_directories(directoryName);

  createPureComponentsPlotScript(directoryNameString);
  createMixturePlotScript(directoryNameString);
  createMixtureAdsorbedMolFractionPlotScript(directoryNameString);

  // create the output files
  std::vector<std::ofstream> streams;
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = directoryNameString + "component_" + std::to_string(components[i].componentId) + "_" + components[i].name + ".data";
    streams.emplace_back(std::ofstream{ fileName });
  }

  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(streams[i], "# column 1: total pressure [Pa]\n");
    std::print(streams[i], "# column 2: pure component isotherm value\n");
    std::print(streams[i], "# column 3: mixture component isotherm value\n");
    std::print(streams[i], "# column 4: gas-phase mol-fraction\n");
    std::print(streams[i], "# column 5: adsorbed phase mol-fraction\n");
  }

  for(size_t i = 0; i < system.pressure_range.numberOfPoints; ++i)
  {
    std::pair<double, double> performance = predictMixture(Yi, pressures[i], Xi, Ni, &cachedP0[0], cachedPsi);
    std::print(stream, "Pressure: {:10e}  iterations: {}\n", pressures[i], performance.first);

    for (size_t j = 0; j < Ncomp; j++)
    {
      std::print(streams[j], "{} {} {} {} {}\n", pressures[i], components[j].isotherm.value(pressures[i]), Ni[j], Yi[j], Xi[j]);
    }
  }
}


void MixturePrediction::createPureComponentsPlotScript(std::string directoryNameString)
{
  std::ofstream stream(directoryNameString + "plot_pure_components");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set xlabel 'Total bulk fluid phase fugacity, {{/Helvetica-Italic f}} / Pa' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Absolute loading, {{/Helvetica-Italic q}}_i' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set bmargin 4\n");
  if(pressureScale == PressureScale::Log) 
  {
    std::print(stream, "set key top left width 2 samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10' maxcolumns 2\n");
    std::print(stream, "set log x\n");
    std::print(stream, "set format x '10^{{%T}}'\n");
  }
  else 
  {
    std::print(stream, "set key outside right width 2 samplen 2.5 height 0.5 spacing 1.5 font \"Helvetica, 10\" maxcolumns 2\n");
  }
  std::print(stream, "set key title '{} {{/:Italic T}}={} K'\n", displayName, temperature);

  std::print(stream, "set output 'pure_component_isotherms.pdf'\n");
  std::print(stream, "set term pdf color solid\n");

  std::print(stream, "plot \\\n");
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(components[i].componentId) + "_" + components[i].name + ".data";
    std::print(stream, "     '{}' us ($1):($2) title '{}' with po{}\n", 
               fileName, components[i].name, (i < Ncomp - 1 ? ",\\" : ""));
  }
}

void MixturePrediction::createMixturePlotScript(std::string directoryNameString)
{
  std::ofstream stream(directoryNameString + "plot_mixture");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set xlabel 'Total bulk fluid phase fugacity, {{/Helvetica-Italic f}} / Pa' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Absolute loading, {{/Helvetica-Italic q}}_i' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set bmargin 4\n");
  if(pressureScale == PressureScale::Log) 
  {
    std::print(stream, "set key top left samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10' maxcolumns 2\n");
    std::print(stream, "set log x\n");
    std::print(stream, "set format x '10^{{%T}}'\n");
  }
  else 
  {
    std::print(stream, "set key outside right samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10' maxcolumns 2\n");
  }
  std::print(stream, "set key title '{} {{/:Italic T}}={} K\n", displayName, temperature);

  std::print(stream, "set output 'iast_mixture.pdf'\n");
  std::print(stream, "set term pdf color solid\n");

  std::print(stream, "plot \\\n");
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(components[i].componentId) + "_" + components[i].name + ".data";
    std::print(stream, "    '{}' us ($1):($3) title '{} (y_i={})' with po{}\n", 
               fileName, components[i].name,  components[i].molFraction, (i < Ncomp - 1 ? ",\\" : ""));
  }
}

void MixturePrediction::createMixtureAdsorbedMolFractionPlotScript(std::string directoryNameString)
{
  std::ofstream stream(directoryNameString + "plot_mixture_mol_fractions");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set xlabel 'Total bulk fluid phase fugacity, {{/Helvetica-Italic f}} / Pa' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Adsorbed mol-fraction, {{/Helvetica-Italic Y}}_i / [-]' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set bmargin 4\n");
  if(pressureScale == PressureScale::Log) 
  {
    std::print(stream, "set key outside right samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10' maxcolumns 2\n");
    std::print(stream, "set log x\n");
    std::print(stream, "set format x '10^{{%T}}'\n");
  }
  else 
  {
    std::print(stream, "set key outside right samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10' maxcolumns 2\n");
  }
  std::print(stream, "set key title '{} {{/:Italic T}}={} K'\n", displayName, temperature);

  std::print(stream, "set output 'iast_mixture_mol_fractions.pdf'\n");
  std::print(stream, "set term pdf color solid\n");

  std::print(stream, "plot \\\n");
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(components[i].componentId) + "_" + components[i].name + ".data";
    std::print(stream, "    '{}' us ($1):($5) title '{} (y_i={})' with po{}\n", 
               fileName, components[i].name, components[i].molFraction, (i < Ncomp - 1 ? ",\\" : ""));
  }
}

void MixturePrediction::printErrorStatus(double psi, double sum, double P, const std::vector<double> Yi, double cachedP0[]) const
{
  std::cout << "psi: " << psi << std::endl;
  std::cout << "sum: " << sum << std::endl;
  for(size_t i = 0; i < Ncomp; ++i)
    std::cout << "cachedP0: " << cachedP0[i] << std::endl;
  for(size_t i = 0; i < Ncomp; ++i)
  {
    double value = components[i].isotherm.inversePressureForPsi(psi, cachedP0[i]);
    std::cout << "inversePressure: " << value << std::endl;
  }
  std::cout << "P: " << P << std::endl;
  for(size_t i = 0; i < Ncomp; ++i)
  {
    std::cout << "Yi[i] "<< i << " " << Yi[i] << std::endl;
  }
}
