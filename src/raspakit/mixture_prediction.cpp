module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <ostream>
#include <print>
#include <span>
#include <sstream>
#include <string>
#include <vector>
#endif

module mixture_prediction;

#ifndef USE_LEGACY_HEADERS
import <functional>;
import <vector>;
import <span>;
import <cmath>;
import <string>;
import <iostream>;
import <iomanip>;
import <fstream>;
import <sstream>;
import <limits>;
import <algorithm>;
import <numeric>;
import <ostream>;
import <filesystem>;
import <print>;
#endif

import stringutils;
import atom;
import component;
import isotherm;
import multi_site_isotherm;
import system;
import simulationbox;
import pressure_range;
import bond_potential;

bool LangmuirLoadingSorter(Component const &lhs, Component const &rhs)
{
  if (lhs.isCarrierGas) return false;
  if (rhs.isCarrierGas) return true;
  return lhs.isotherm.sites[0].parameters[0] < rhs.isotherm.sites[0].parameters[0];
}

// allow std::pairs to be added
template <typename T, typename U>
std::pair<T, U> operator+(const std::pair<T, U> &l, const std::pair<T, U> &r)
{
  return {l.first + r.first, l.second + r.second};
}
template <typename T, typename U>
std::pair<T, U> &operator+=(std::pair<T, U> &l, const std::pair<T, U> &r)
{
  l.first += r.first;
  l.second += r.second;
  return l;
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
                if (v[i1].isCarrierGas) return false;
                if (v[i2].isCarrierGas) return true;
                return true;
                return v[i1].isotherm.parameters(0) < v[i2].isotherm.parameters(0);
              });

  return idx;
}

MixturePrediction::MixturePrediction(const System &system)
    : system(system),
      displayName(system.components.front().name),
      components(system.components),
      sortedComponentIndices(sortedIndices(components)),
      sortedComponents(components.begin(), components.end()),
      Ncomp(components.size()),
      Nsorted(components.size() - system.numberOfCarrierGases),
      predictionMethod(system.mixturePredictionMethod),
      maxIsothermTerms(system.maxIsothermTerms),
      alpha1(Ncomp),
      alpha2(Ncomp),
      alpha_prod(Ncomp),
      x(Ncomp),
      pstar(Nsorted),
      psi(Nsorted),
      G(Nsorted),
      delta(Nsorted),
      Phi(Nsorted * Nsorted),
      temperature(system.temperature)
{
  if (predictionMethod == MultiSiteIsotherm::PredictionMethod::EI)
  {
    std::sort(sortedComponents.begin(), sortedComponents.end(), &LangmuirLoadingSorter);
  }
  else if (predictionMethod == MultiSiteIsotherm::PredictionMethod::SEI)
  {
    for (size_t i = 0; i < maxIsothermTerms; ++i)
    {
      for (size_t j = 0; j < Ncomp; ++j)
      {
        if (j != system.carrierGasComponent)
        {
          segregatedSortedComponents[i][j].isotherm.sites[0] = components[j].isotherm.sites[i];
          segregatedSortedComponents[i][j].isotherm.numberOfSites = 1;
        }
      }
    }
    for (size_t i = 0; i < maxIsothermTerms; ++i)
    {
      std::sort(segregatedSortedComponents[i].begin(), segregatedSortedComponents[i].end(), &LangmuirLoadingSorter);
    }
  }
  else
  {
    auto it = sortedComponents.begin() + static_cast<std::ptrdiff_t>(system.carrierGasComponent);
    std::rotate(it, it + 1, sortedComponents.end());
  }
}

std::string MixturePrediction::writeHeader() const
{
  std::ostringstream stream;

  /*

  switch(predictionMethod)
  {
    case MultiSiteIsotherm::PredictionMethod::IAST:
    default:
      std::print(stream, "Component data\n");
      std::print(stream, "=======================================================\n");
      for(size_t i = 0; i < Ncomp; ++i)
      {
        std::print(stream, components[i].printBreakthroughStatus());
        std::print(stream, "\n");
      }
      break;
    case MultiSiteIsotherm::PredictionMethod::EI:
      std::print(stream, "Sorted component data\n");
      std::print(stream, "=======================================================\n");
      for(size_t i = 0; i < Ncomp; ++i)
      {
        std::print(stream, sortedComponents[i].get().printBreakthroughStatus());
        std::print(stream, "\n");
      }
      break;
  }
  */

  return stream.str();
}

std::pair<size_t, size_t> MixturePrediction::predictMixture(const std::vector<double> &Yi, const double &P,
                                                            std::vector<double> &Xi, std::vector<double> &Ni,
                                                            double *cachedP0, double *cachedPsi)
{
  const double tiny = 1.0e-10;

  if (P < 0.0)
  {
    printErrorStatus(0.0, 0.0, P, Yi, cachedP0);
    throw std::runtime_error("Error (IAST): negative total pressure\n");
  }

  double sumYi = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    sumYi += Yi[i];
  }
  if (std::abs(sumYi - 1.0) > 1e-15)
  {
    printErrorStatus(0.0, sumYi, P, Yi, cachedP0);
    throw std::runtime_error("Error (IAST): sum Yi at IAST start not unity\n");
  }

  // if only an inert component present
  // this happens at the beginning of the simulation when the whole column is filled with the carrier gas
  if (std::abs(Yi[system.carrierGasComponent] - 1.0) < tiny)
  {
    for (size_t i = 0; i < Ncomp; ++i)
    {
      Xi[i] = 0.0;
      Ni[i] = 0.0;
    }

    // do not count it for the IAST statistics
    return std::make_pair(0, 0);
  }

  switch (predictionMethod)
  {
    case MultiSiteIsotherm::PredictionMethod::IAST:
    default:
      switch (iastMethod)
      {
        case IASTMethod::FastIAST:
        default:
          return computeFastIAST(Yi, P, Xi, Ni, cachedP0, cachedPsi);
        case IASTMethod::NestedLoopBisection:
          return computeIASTNestedLoopBisection(Yi, P, Xi, Ni, cachedP0, cachedPsi);
      }
    case MultiSiteIsotherm::PredictionMethod::SIAST:
      switch (iastMethod)
      {
        case IASTMethod::FastIAST:
        default:
          return computeFastSIAST(Yi, P, Xi, Ni, cachedP0, cachedPsi);
        case IASTMethod::NestedLoopBisection:
          return computeSIASTNestedLoopBisection(Yi, P, Xi, Ni, cachedP0, cachedPsi);
      }
    case MultiSiteIsotherm::PredictionMethod::EI:
      return computeExplicitIsotherm(Yi, P, Xi, Ni);
    case MultiSiteIsotherm::PredictionMethod::SEI:
      return computeSegratedExplicitIsotherm(Yi, P, Xi, Ni);
  }
}

// Yi  = gas phase molefraction
// P   = total pressure
// Xi  = adsorbed phase molefraction
// Ni  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeFastIAST(const std::vector<double> &Yi, const double &P,
                                                             std::vector<double> &Xi, std::vector<double> &Ni,
                                                             double *cachedP0, double *cachedPsi)
{
  const double tiny = 1.0e-13;

  size_t numberOfIASTSteps = 0;

  std::fill(pstar.begin(), pstar.end(), 0.0);
  std::fill(G.begin(), G.end(), 0.0);
  std::fill(delta.begin(), delta.end(), 0.0);
  std::fill(Phi.begin(), Phi.end(), 0.0);

  if (cachedPsi[0] > 0.0)
  {
    for (size_t i = 0; i < Nsorted; ++i)
    {
      pstar[i] = cachedP0[sortedComponents[i].get().componentId];
    }
  }
  else
  {
    double initial_psi = 0.0;
    for (size_t i = 0; i < Nsorted; ++i)
    {
      double temp_psi =
          Yi[sortedComponents[i].get().componentId] * sortedComponents[i].get().isotherm.psiForPressure(P);
      initial_psi += temp_psi;
    }
    cachedPsi[0] = initial_psi;

    double cachevalue = 0.0;
    for (size_t i = 0; i < Nsorted; ++i)
    {
      pstar[i] = 1.0 / sortedComponents[i].get().isotherm.inversePressureForPsi(initial_psi, cachevalue);
    }
  }

  double error = 1.0;
  double sum_xi = 0.0;
  do
  {
    // compute G
    for (size_t i = 0; i < Nsorted - 1; ++i)
    {
      G[i] = sortedComponents[i].get().isotherm.psiForPressure(pstar[i]) -
             sortedComponents[Nsorted - 1].get().isotherm.psiForPressure(pstar[Nsorted - 1]);
    }

    G[Nsorted - 1] = 0.0;
    for (size_t i = 0; i < Nsorted; i++)
    {
      G[Nsorted - 1] += Yi[sortedComponents[i].get().componentId] * P / pstar[i];
    }
    G[Nsorted - 1] -= 1.0;

    // compute Jacobian matrix Phi
    for (size_t i = 0; i < Nsorted - 1; i++)
    {
      Phi[i + i * Nsorted] = sortedComponents[i].get().isotherm.value(pstar[i]) / pstar[i];
    }
    for (size_t i = 0; i < Nsorted - 1; i++)
    {
      Phi[i + (Nsorted - 1) * Nsorted] =
          -sortedComponents[Nsorted - 1].get().isotherm.value(pstar[Nsorted - 1]) / pstar[Nsorted - 1];
    }
    for (size_t i = 0; i < Nsorted; i++)
    {
      Phi[(Nsorted - 1) + i * Nsorted] = -Yi[sortedComponents[i].get().componentId] * P / (pstar[i] * pstar[i]);
    }

    // corrections
    for (size_t i = 0; i < Nsorted - 1; i++)
    {
      Phi[(Nsorted - 1) + (Nsorted - 1) * Nsorted] -=
          Phi[(Nsorted - 1) + i * Nsorted] * Phi[i + (Nsorted - 1) * Nsorted] / Phi[i + i * Nsorted];
      G[Nsorted - 1] -= Phi[(Nsorted - 1) + i * Nsorted] * G[i] / Phi[i + i * Nsorted];
    }

    // compute delta
    delta[Nsorted - 1] = G[Nsorted - 1] / Phi[(Nsorted - 1) + (Nsorted - 1) * Nsorted];

    // trick to loop downward from Nsorted - 2 to and including zero (still using size_t as index)
    for (size_t i = Nsorted - 1; i-- != 0;)
    {
      delta[i] = (G[i] - delta[Nsorted - 1] * Phi[i + (Nsorted - 1) * Nsorted]) / Phi[i + i * Nsorted];
    }

    // update pstar
    for (size_t i = 0; i < Nsorted; i++)
    {
      double newvalue = pstar[i] - delta[i];
      if (newvalue > 0.0)
        pstar[i] = newvalue;
      else
      {
        pstar[i] = 0.5 * pstar[i];
      }
    }

    // compute error in psi's
    for (size_t i = 0; i < Nsorted; i++)
    {
      psi[i] = sortedComponents[i].get().isotherm.psiForPressure(pstar[i]);
    }

    sum_xi = 0.0;
    for (size_t i = 0; i < Nsorted; ++i)
    {
      sum_xi += Yi[sortedComponents[i].get().componentId] * P / std::max(pstar[i], 1e-15);
    }

    double avg = std::accumulate(std::begin(psi), std::end(psi), 0.0) / static_cast<double>(psi.size());

    double accum = 0.0;
    std::for_each(std::begin(psi), std::end(psi), [&](const double d) { accum += (d - avg) * (d - avg); });

    error = std::sqrt(accum / static_cast<double>(psi.size() - 1));

    numberOfIASTSteps++;
  } while (!(((error < tiny) && (std::fabs(sum_xi - 1.0) < 1e-10)) || (numberOfIASTSteps >= 50)));

  for (size_t i = 0; i < Nsorted; ++i)
  {
    cachedP0[sortedComponents[i].get().componentId] = pstar[i];
  }

  for (size_t i = 0; i < Nsorted; ++i)
  {
    Xi[sortedComponents[i].get().componentId] =
        Yi[sortedComponents[i].get().componentId] * P / std::max(pstar[i], 1e-15);
  }
  if (system.numberOfCarrierGases > 0)
  {
    Xi[system.carrierGasComponent] = 0.0;
  }

  double sum = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    sum += Xi[i];
  }
  for (size_t i = 0; i < Ncomp; ++i)
  {
    Xi[i] /= sum;
  }

  double inverse_q_total = 0.0;
  for (size_t i = 0; i < Nsorted; ++i)
  {
    inverse_q_total += Xi[sortedComponents[i].get().componentId] / sortedComponents[i].get().isotherm.value(pstar[i]);
  }
  for (size_t i = 0; i < Ncomp; ++i)
  {
    Ni[i] = Xi[i] / inverse_q_total;
  }
  if (system.numberOfCarrierGases > 0)
  {
    Ni[system.carrierGasComponent] = 0.0;
  }

  return std::make_pair(numberOfIASTSteps, 1);
}

// Yi  = gas phase molefraction
// P   = total pressure
// Xi  = adsorbed phase molefraction
// Ni  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeFastSIAST(const std::vector<double> &Yi, const double &P,
                                                              std::vector<double> &Xi, std::vector<double> &Ni,
                                                              double *cachedP0, double *cachedPsi)
{
  std::fill(Xi.begin(), Xi.end(), 0.0);
  std::fill(Ni.begin(), Ni.end(), 0.0);

  std::pair<size_t, size_t> acc;
  for (size_t i = 0; i < maxIsothermTerms; ++i)
  {
    acc += computeFastSIAST(i, Yi, P, Xi, Ni, cachedP0, cachedPsi);
  }

  double N = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    N += Ni[i];
  }
  for (size_t i = 0; i < Ncomp; ++i)
  {
    Xi[i] = Ni[i] / N;
  }

  return acc;
}

// computes IAST per term
// Yi  = gas phase molefraction
// P   = total pressure
// Xi  = adsorbed phase molefraction
// Ni  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeFastSIAST(size_t site, const std::vector<double> &Yi,
                                                              const double &P, std::vector<double> &Xi,
                                                              std::vector<double> &Ni, double *cachedP0,
                                                              double *cachedPsi)
{
  const double tiny = 1.0e-13;

  size_t numberOfIASTSteps = 0;

  std::fill(pstar.begin(), pstar.end(), 0.0);
  std::fill(G.begin(), G.end(), 0.0);
  std::fill(delta.begin(), delta.end(), 0.0);
  std::fill(Phi.begin(), Phi.end(), 0.0);

  if (cachedPsi[site] > tiny)
  {
    for (size_t i = 0; i < Nsorted; ++i)
    {
      pstar[i] = cachedP0[sortedComponents[i].get().componentId + site * Ncomp];
    }
  }
  else
  {
    double initial_psi = 0.0;
    for (size_t i = 0; i < Nsorted; ++i)
    {
      double temp_psi =
          Yi[sortedComponents[i].get().componentId] * sortedComponents[i].get().isotherm.psiForPressure(site, P);
      initial_psi += temp_psi;
    }
    cachedPsi[site] = initial_psi;

    double cachevalue = 0.0;
    for (size_t i = 0; i < Nsorted; ++i)
    {
      pstar[i] = 1.0 / sortedComponents[i].get().isotherm.inversePressureForPsi(site, initial_psi, cachevalue);
    }
  }

  double error = 1.0;
  double sum_xi = 1.0;
  do
  {
    // compute G
    for (size_t i = 0; i < Nsorted - 1; ++i)
    {
      G[i] = sortedComponents[i].get().isotherm.psiForPressure(site, pstar[i]) -
             sortedComponents[Nsorted - 1].get().isotherm.psiForPressure(site, pstar[Nsorted - 1]);
    }

    G[Nsorted - 1] = 0.0;
    for (size_t i = 0; i < Nsorted; i++)
    {
      G[Nsorted - 1] += Yi[sortedComponents[i].get().componentId] * P / pstar[i];
    }
    G[Nsorted - 1] -= 1.0;

    // compute Jacobian matrix Phi
    for (size_t i = 0; i < Nsorted - 1; i++)
    {
      Phi[i + i * Nsorted] = sortedComponents[i].get().isotherm.value(site, pstar[i]) / pstar[i];
    }
    for (size_t i = 0; i < Nsorted - 1; i++)
    {
      Phi[i + (Nsorted - 1) * Nsorted] =
          -sortedComponents[Nsorted - 1].get().isotherm.value(site, pstar[Nsorted - 1]) / pstar[Nsorted - 1];
    }
    for (size_t i = 0; i < Nsorted; i++)
    {
      Phi[(Nsorted - 1) + i * Nsorted] = -Yi[sortedComponents[i].get().componentId] * P / (pstar[i] * pstar[i]);
    }

    // corrections
    for (size_t i = 0; i < Nsorted - 1; i++)
    {
      Phi[(Nsorted - 1) + (Nsorted - 1) * Nsorted] -=
          Phi[(Nsorted - 1) + i * Nsorted] * Phi[i + (Nsorted - 1) * Nsorted] / Phi[i + i * Nsorted];
      G[Nsorted - 1] -= Phi[(Nsorted - 1) + i * Nsorted] * G[i] / Phi[i + i * Nsorted];
    }

    // compute delta
    delta[Nsorted - 1] = G[Nsorted - 1] / Phi[(Nsorted - 1) + (Nsorted - 1) * Nsorted];

    // trick to loop downward from Nsorted - 2 to and including zero (still using size_t as index)
    for (size_t i = Nsorted - 1; i-- != 0;)
    {
      delta[i] = (G[i] - delta[Nsorted - 1] * Phi[i + (Nsorted - 1) * Nsorted]) / Phi[i + i * Nsorted];
    }

    // update pstar
    for (size_t i = 0; i < Nsorted; i++)
    {
      double newvalue = pstar[i] - delta[i];
      if (newvalue > 0.0)
        pstar[i] = newvalue;
      else
      {
        pstar[i] = 0.5 * pstar[i];
      }
    }

    // compute error in psi's
    for (size_t i = 0; i < Nsorted; i++)
    {
      psi[i] = sortedComponents[i].get().isotherm.psiForPressure(site, pstar[i]);
    }

    sum_xi = 0.0;
    for (size_t i = 0; i < Nsorted; ++i)
    {
      sum_xi += Yi[sortedComponents[i].get().componentId] * P / std::max(pstar[i], 1e-15);
    }

    double avg = std::accumulate(std::begin(psi), std::end(psi), 0.0) / static_cast<double>(psi.size());

    double accum = 0.0;
    std::for_each(std::begin(psi), std::end(psi), [&](const double d) { accum += (d - avg) * (d - avg); });

    error = std::sqrt(accum / static_cast<double>(psi.size() - 1));

    numberOfIASTSteps++;
  } while (!(((error < tiny) && (std::fabs(sum_xi - 1.0) < 1e-10)) || (numberOfIASTSteps >= 50)));

  for (size_t i = 0; i < Nsorted; ++i)
  {
    cachedP0[sortedComponents[i].get().componentId + site * Ncomp] = pstar[i];
  }

  for (size_t i = 0; i < Nsorted; ++i)
  {
    Xi[sortedComponents[i].get().componentId] =
        Yi[sortedComponents[i].get().componentId] * P / std::max(pstar[i], 1e-15);
  }
  if (system.numberOfCarrierGases > 0)
  {
    Xi[system.carrierGasComponent] = 0.0;
  }

  double sum = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    sum += Xi[i];
  }
  for (size_t i = 0; i < Ncomp; ++i)
  {
    Xi[i] /= sum;
  }

  double inverse_q_total = 0.0;
  for (size_t i = 0; i < Nsorted; ++i)
  {
    inverse_q_total +=
        Xi[sortedComponents[i].get().componentId] / sortedComponents[i].get().isotherm.value(site, pstar[i]);
  }
  for (size_t i = 0; i < Ncomp; ++i)
  {
    Ni[i] += Xi[i] / inverse_q_total;
  }
  if (system.numberOfCarrierGases > 0)
  {
    Ni[system.carrierGasComponent] = 0.0;
  }

  return std::make_pair(numberOfIASTSteps, 1);
}

// Yi  = gas phase molefraction
// P   = total pressure
// Xi  = adsorbed phase molefraction
// Ni  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeIASTNestedLoopBisection(const std::vector<double> &Yi,
                                                                            const double &P, std::vector<double> &Xi,
                                                                            std::vector<double> &Ni, double *cachedP0,
                                                                            double *cachedPsi)
{
  const double tiny = 1.0e-15;

  double initial_psi = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    initial_psi += Yi[i] * components[i].isotherm.psiForPressure(P);
  }

  if (initial_psi < tiny)
  {
    // nothing is adsorbing
    for (size_t i = 0; i < Ncomp; ++i)
    {
      Xi[i] = 0.0;
      Ni[i] = 0.0;
    }

    // do not count it for the IAST statistics
    return std::make_pair(0, 0);
  }

  // condition 1: same reduced grand potential for all components (done by using a single variable)
  // condition 2: mol-fractions add up to unity

  double psi_value = 0.0;
  size_t nr_steps = 0;
  if (cachedPsi[0] > tiny)
  {
    initial_psi = cachedPsi[0];
  }
  // for this initial estimate 'initial_psi' compute the sum of mol-fractions
  double sumXi = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(initial_psi, cachedP0[i]);
  }

  // initialize the bisection algorithm
  double left_bracket = initial_psi;
  double right_bracket = initial_psi;
  if (sumXi > 1.0)
  {
    do
    {
      right_bracket *= 2.0;

      sumXi = 0.0;
      for (size_t i = 0; i < Ncomp; ++i)
      {
        sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(right_bracket, cachedP0[i]);
      }
      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::cout << "Left bracket: " << left_bracket << std::endl;
        std::cout << "Right bracket: " << right_bracket << std::endl;
        printErrorStatus(0.0, sumXi, P, Yi, cachedP0);
        throw std::runtime_error("Error (IAST bisection): initial bracketing (for sum > 1) does NOT converge\n");
      }
    } while (sumXi > 1.0);
  }
  else
  {
    // Make an initial estimate for the reduced grandpotential when the
    // sum of the molefractions is larger than 1
    do
    {
      left_bracket *= 0.5;

      sumXi = 0.0;
      for (size_t i = 0; i < Ncomp; ++i)
      {
        sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(left_bracket, cachedP0[i]);
      }
      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::cout << "Left bracket: " << left_bracket << std::endl;
        std::cout << "Right bracket: " << right_bracket << std::endl;
        printErrorStatus(0.0, sumXi, P, Yi, cachedP0);
        throw std::runtime_error("Error (IAST bisection): initial bracketing (for sum < 1) does NOT converge\n");
      }
    } while (sumXi < 1.0);
  }

  // bisection algorithm
  size_t numberOfIASTSteps = 0;
  do
  {
    psi_value = 0.5 * (left_bracket + right_bracket);

    sumXi = 0.0;
    for (size_t i = 0; i < Ncomp; ++i)
    {
      sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(psi_value, cachedP0[i]);
    }

    if (sumXi > 1.0)
    {
      left_bracket = psi_value;
    }
    else
    {
      right_bracket = psi_value;
    }

    ++numberOfIASTSteps;
    if (numberOfIASTSteps > 100000)
    {
      throw std::runtime_error("Error (IAST bisection): NO convergence\n");
    }
  } while (std::abs(left_bracket - right_bracket) / std::abs(left_bracket + right_bracket) > tiny);  // convergence test

  psi_value = 0.5 * (left_bracket + right_bracket);

  sumXi = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(psi_value, cachedP0[i]);
  }

  // cache the value of psi for subsequent use
  cachedPsi[0] = psi_value;

  // calculate mol-fractions in adsorbed phase and total loading
  double inverse_q_total = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    double ip = components[i].isotherm.inversePressureForPsi(psi_value, cachedP0[i]);
    Xi[i] = Yi[i] * P * ip / sumXi;

    if (Xi[i] > tiny)
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
    for (size_t i = 0; i < Ncomp; ++i)
    {
      Ni[i] = 0.0;
    }
  }
  else
  {
    for (size_t i = 0; i < Ncomp; ++i)
    {
      Ni[i] = Xi[i] / inverse_q_total;
    }
  }

  return std::make_pair(numberOfIASTSteps, 1);
}

// Yi  = gas phase molefraction
// P   = total pressure
// Xi  = adsorbed phase molefraction
// Ni  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeSIASTNestedLoopBisection(const std::vector<double> &Yi,
                                                                             const double &P, std::vector<double> &Xi,
                                                                             std::vector<double> &Ni, double *cachedP0,
                                                                             double *cachedPsi)
{
  std::fill(Xi.begin(), Xi.end(), 0.0);
  std::fill(Ni.begin(), Ni.end(), 0.0);

  std::pair<size_t, size_t> acc;
  for (size_t i = 0; i < maxIsothermTerms; ++i)
  {
    acc += computeSIASTNestedLoopBisection(i, Yi, P, Xi, Ni, cachedP0, cachedPsi);
  }

  double N = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    N += Ni[i];
  }
  for (size_t i = 0; i < Ncomp; ++i)
  {
    Xi[i] = Ni[i] / N;
  }

  return acc;
}

// computes IAST per term
// Yi  = gas phase molefraction
// P   = total pressure
// Xi  = adsorbed phase molefraction
// Ni  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeSIASTNestedLoopBisection(size_t site, const std::vector<double> &Yi,
                                                                             const double &P, std::vector<double> &Xi,
                                                                             std::vector<double> &Ni, double *cachedP0,
                                                                             double *cachedPsi)
{
  const double tiny = 1.0e-15;

  double initial_psi = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    initial_psi += Yi[i] * components[i].isotherm.psiForPressure(site, P);
  }

  if (initial_psi < tiny)
  {
    // nothing is adsorbing
    // do not count it for the IAST statistics
    return std::make_pair(0, 0);
  }

  // condition 1: same reduced grand potential for all components (done by using a single variable)
  // condition 2: mol-fractions add up to unity

  double psi_value = 0.0;
  size_t nr_steps = 0;
  if (cachedPsi[site] > tiny)
  {
    initial_psi = cachedPsi[site];
  }
  // for this initial estimate 'initial_psi' compute the sum of mol-fractions
  double sumXi = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(site, initial_psi, cachedP0[i + Ncomp * site]);
  }

  // initialize the bisection algorithm
  double left_bracket = initial_psi;
  double right_bracket = initial_psi;
  if (sumXi > 1.0)
  {
    do
    {
      right_bracket *= 2.0;

      sumXi = 0.0;
      for (size_t i = 0; i < Ncomp; ++i)
      {
        sumXi +=
            Yi[i] * P * components[i].isotherm.inversePressureForPsi(site, right_bracket, cachedP0[i + Ncomp * site]);
      }
      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::cout << "Left bracket: " << left_bracket << std::endl;
        std::cout << "Right bracket: " << right_bracket << std::endl;
        printErrorStatus(0.0, sumXi, P, Yi, cachedP0);
        throw std::runtime_error("Error (IAST bisection): initial bracketing (for sum > 1) does NOT converge\n");
      }
    } while (sumXi > 1.0);
  }
  else
  {
    // Make an initial estimate for the reduced grandpotential when the
    // sum of the molefractions is larger than 1
    do
    {
      left_bracket *= 0.5;

      sumXi = 0.0;
      for (size_t i = 0; i < Ncomp; ++i)
      {
        sumXi +=
            Yi[i] * P * components[i].isotherm.inversePressureForPsi(site, left_bracket, cachedP0[i + Ncomp * site]);
      }
      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::cout << "Left bracket: " << left_bracket << std::endl;
        std::cout << "Right bracket: " << right_bracket << std::endl;
        printErrorStatus(0.0, sumXi, P, Yi, cachedP0);
        throw std::runtime_error("Error (IAST bisection): initial bracketing (for sum < 1) does NOT converge\n");
      }
    } while (sumXi < 1.0);
  }

  // bisection algorithm
  size_t numberOfIASTSteps = 0;
  do
  {
    psi_value = 0.5 * (left_bracket + right_bracket);

    sumXi = 0.0;
    for (size_t i = 0; i < Ncomp; ++i)
    {
      sumXi += Yi[i] * P * components[i].isotherm.inversePressureForPsi(site, psi_value, cachedP0[i + Ncomp * site]);
    }

    if (sumXi > 1.0)
    {
      left_bracket = psi_value;
    }
    else
    {
      right_bracket = psi_value;
    }

    ++numberOfIASTSteps;
    if (numberOfIASTSteps > 100000)
    {
      throw std::runtime_error("Error (IAST bisection): NO convergence\n");
    }
  } while (std::abs(left_bracket - right_bracket) / std::abs(left_bracket + right_bracket) > tiny);  // convergence test

  psi_value = 0.5 * (left_bracket + right_bracket);

  // cache the value of psi for subsequent use
  cachedPsi[site] = psi_value;

  // calculate mol-fractions in adsorbed phase and total loading
  double inverse_q_total = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    double ip = components[i].isotherm.inversePressureForPsi(site, psi_value, cachedP0[i + Ncomp * site]);
    Xi[i] = Yi[i] * P * ip;

    if (Xi[i] > tiny)
    {
      inverse_q_total += Xi[i] / components[i].isotherm.value(site, 1.0 / ip);
    }
  }

  // calculate loading for all of the components
  if (inverse_q_total > 0.0)
  {
    for (size_t i = 0; i < Ncomp; ++i)
    {
      Ni[i] += Xi[i] / inverse_q_total;
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

std::pair<size_t, size_t> MixturePrediction::computeExplicitIsotherm(const std::vector<double> &Yi, const double &P,
                                                                     std::vector<double> &Xi, std::vector<double> &Ni)
{
  x[0] = 1.0;
  for (size_t i = 1; i < Ncomp; ++i)
  {
    x[i] = sortedComponents[i].get().isotherm.sites[0].parameters[0] /
           sortedComponents[i - 1].get().isotherm.sites[0].parameters[0];
  }

  alpha1[Ncomp - 1] = std::pow((1.0 + sortedComponents[Ncomp - 1].get().isotherm.sites[0].parameters[1] *
                                          Yi[sortedComponents[Ncomp - 1].get().componentId] * P),
                               x[Ncomp - 1]);
  alpha2[Ncomp - 1] = 1.0 + sortedComponents[Ncomp - 1].get().isotherm.sites[0].parameters[1] *
                                Yi[sortedComponents[Ncomp - 1].get().componentId] * P;
  for (size_t i = Ncomp - 2; i > 0; i--)
  {
    alpha1[i] = std::pow((alpha1[i + 1] + sortedComponents[i].get().isotherm.sites[0].parameters[1] *
                                              Yi[sortedComponents[i].get().componentId] * P),
                         x[i]);
    alpha2[i] = alpha1[i + 1] + sortedComponents[i].get().isotherm.sites[0].parameters[1] *
                                    Yi[sortedComponents[i].get().componentId] * P;
  }
  alpha1[0] = alpha1[1] +
              sortedComponents[0].get().isotherm.sites[0].parameters[1] * Yi[sortedComponents[0].get().componentId] * P;
  alpha2[0] = alpha1[1] +
              sortedComponents[0].get().isotherm.sites[0].parameters[1] * Yi[sortedComponents[0].get().componentId] * P;

  double beta = alpha2[0];

  alpha_prod[0] = 1.0;
  for (size_t i = 1; i < Ncomp; ++i)
  {
    alpha_prod[i] = (alpha1[i] / alpha2[i]) * alpha_prod[i - 1];
  }

  for (size_t i = 0; i < Ncomp; ++i)
  {
    size_t index = sortedComponents[i].get().componentId;
    Ni[index] = sortedComponents[i].get().isotherm.sites[0].parameters[0] *
                sortedComponents[i].get().isotherm.sites[0].parameters[1] * Yi[index] * P * alpha_prod[i] / beta;
  }
  double N = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    N += Ni[i];
  }
  for (size_t i = 0; i < Ncomp; ++i)
  {
    Xi[i] = Ni[i] / N;
  }

  return std::make_pair(1, 1);
}

std::pair<size_t, size_t> MixturePrediction::computeSegratedExplicitIsotherm(const std::vector<double> &Yi,
                                                                             const double &P, std::vector<double> &Xi,
                                                                             std::vector<double> &Ni)
{
  std::fill(Xi.begin(), Xi.end(), 0.0);
  std::fill(Ni.begin(), Ni.end(), 0.0);

  std::pair<size_t, size_t> acc;
  for (size_t i = 0; i < maxIsothermTerms; ++i)
  {
    acc += computeSegratedExplicitIsotherm(i, Yi, P, Xi, Ni);
  }

  double N = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    N += Ni[i];
  }
  for (size_t i = 0; i < Ncomp; ++i)
  {
    Xi[i] = Ni[i] / N;
  }

  return acc;
}

std::pair<size_t, size_t> MixturePrediction::computeSegratedExplicitIsotherm(size_t site, const std::vector<double> &Yi,
                                                                             const double &P, std::vector<double> &Xi,
                                                                             std::vector<double> &Ni)
{
  x[0] = 1.0;
  for (size_t i = 1; i < Ncomp; ++i)
  {
    x[i] = segregatedSortedComponents[site][i].isotherm.sites[0].parameters[0] /
           segregatedSortedComponents[site][i - 1].isotherm.sites[0].parameters[0];
  }

  alpha1[Ncomp - 1] = std::pow((1.0 + segregatedSortedComponents[site][Ncomp - 1].isotherm.sites[0].parameters[1] *
                                          Yi[segregatedSortedComponents[site][Ncomp - 1].componentId] * P),
                               x[Ncomp - 1]);
  alpha2[Ncomp - 1] = 1.0 + segregatedSortedComponents[site][Ncomp - 1].isotherm.sites[0].parameters[1] *
                                Yi[segregatedSortedComponents[site][Ncomp - 1].componentId] * P;
  for (size_t i = Ncomp - 2; i > 0; i--)
  {
    alpha1[i] = std::pow((alpha1[i + 1] + segregatedSortedComponents[site][i].isotherm.sites[0].parameters[1] *
                                              Yi[segregatedSortedComponents[site][i].componentId] * P),
                         x[i]);
    alpha2[i] = alpha1[i + 1] + segregatedSortedComponents[site][i].isotherm.sites[0].parameters[1] *
                                    Yi[segregatedSortedComponents[site][i].componentId] * P;
  }
  alpha1[0] = alpha1[1] + segregatedSortedComponents[site][0].isotherm.sites[0].parameters[1] *
                              Yi[segregatedSortedComponents[site][0].componentId] * P;
  alpha2[0] = alpha1[1] + segregatedSortedComponents[site][0].isotherm.sites[0].parameters[1] *
                              Yi[segregatedSortedComponents[site][0].componentId] * P;

  double beta = alpha2[0];

  alpha_prod[0] = 1.0;
  for (size_t i = 1; i < Ncomp; ++i)
  {
    alpha_prod[i] = (alpha1[i] / alpha2[i]) * alpha_prod[i - 1];
  }

  for (size_t i = 0; i < Ncomp; ++i)
  {
    size_t index = segregatedSortedComponents[site][i].componentId;
    Ni[index] += segregatedSortedComponents[site][i].isotherm.sites[0].parameters[0] *
                 segregatedSortedComponents[site][i].isotherm.sites[0].parameters[1] * Yi[index] * P * alpha_prod[i] /
                 beta;
  }
  double N = 0.0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    N += Ni[i];
  }
  for (size_t i = 0; i < Ncomp; ++i)
  {
    Xi[i] = Ni[i] / N;
  }

  return std::make_pair(1, 1);
}

void MixturePrediction::print() const
{
  std::cout << "Component data\n";
  std::cout << "=======================================================\n";
  std::cout << "maximum isotherm terms:        " << maxIsothermTerms << "\n";
  for (size_t i = 0; i < Ncomp; ++i)
  {
    // std::cout << sortedComponents[i].get().print(i);
    std::cout << "\n";
  }
}

void MixturePrediction::run(std::ostream &stream)
{
  std::vector<double> Yi(Ncomp);
  std::vector<double> Xi(Ncomp);
  std::vector<double> Ni(Ncomp);
  std::vector<double> cachedP0(Ncomp * maxIsothermTerms);
  std::vector<double> cachedPsi(maxIsothermTerms);

  for (size_t i = 0; i < Ncomp; ++i)
  {
    Yi[i] = components[i].molFraction;
  }

  // create the output files
  std::filesystem::create_directory("MixturePrediction");
  std::filesystem::create_directory(std::format("MixturePrediction/System_{}", system.systemId));
  std::vector<std::ofstream> streams;
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = std::format("MixturePrediction/System_{}/component_{}_{}.txt", system.systemId,
                                       std::to_string(i), components[i].name);
    // FIX TODO
    // streams.emplace_back(std::ofstream{ fileName });
  }

  for (size_t i = 0; i < Ncomp; i++)
  {
    streams[i] << "# column 1: total pressure [Pa]\n";
    streams[i] << "# column 2: pure component isotherm value\n";
    streams[i] << "# column 3: mixture component isotherm value\n";
    streams[i] << "# column 4: gas-phase mol-fraction y_i\n";
    streams[i] << "# column 5: adsorbed phase mol-fraction x_i\n";
    streams[i] << "# column 6: hypothetical pressure p_i^*\n";
    streams[i] << "# column 7: reduced grand potential psi_i\n";
    streams[i] << std::setprecision(14);
  }

  for (const double &pressure : system.pressure_range.pressures())
  {
    std::pair<double, double> performance = predictMixture(Yi, pressure, Xi, Ni, &cachedP0[0], &cachedPsi[0]);
    std::print(stream, "Pressure: {:10e}  iterations: {}\n", pressure, performance.first);

    for (size_t j = 0; j < Ncomp; j++)
    {
      double p_star = Yi[j] * pressure / Xi[j];
      std::print(streams[j], "{} {} {} {} {} {}\n", pressure, components[j].isotherm.value(pressure), Ni[j], Yi[j],
                 Xi[j], components[j].isotherm.psiForPressure(p_star));
    }
  }
}

void MixturePrediction::createPureComponentsPlotScript()
{
  std::filesystem::create_directory("MixturePrediction");
  std::filesystem::create_directory(std::format("MixturePrediction/System_{}", system.systemId));
  std::ofstream stream(std::format("MixturePrediction/System_{}/plot_pure_components", system.systemId));

  stream << "set encoding utf8\n";
  stream << "set xlabel 'Total bulk fluid phase fugacity, {/Helvetica-Italic f} / Pa' font \"Helvetica,18\"\n";
  stream << "set ylabel 'Absolute loading, {/Helvetica-Italic q}_i' offset 0.0,0 font \"Helvetica,18\"\n";
  stream << "set bmargin 4\n";
  if (system.pressure_range.scale == PressureRange::Scale::Log)
  {
    stream << "set key top left width 2 samplen 2.5 height 0.5 spacing 1.5 font \"Helvetica, 10\" maxcolumns 2\n";
    stream << "set log x\n";
    stream << "set format x \"10^{%T}\"\n";
    stream << "set xrange[" << system.pressure_range.pressureStart << ":]\n";
  }
  else
  {
    stream << "set key outside right width 2 samplen 2.5 height 0.5 spacing 1.5 font \"Helvetica, 10\" maxcolumns 2\n";
  }
  stream << "set key title '" << displayName << " {/:Italic T}=" << temperature << " K'\n";

  stream << "set output 'pure_component_isotherms.pdf'\n";
  stream << "set term pdf color solid\n";

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 0.5 lw 2 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 0.5 lw 2 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 0.5 lw 2 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 0.5 lw 2 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 0.5 lw 2 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 0.5 lw 2 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 0.5 lw 2 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 0.5 lw 2 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 0.5 lw 2 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 0.5 lw 2 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 0.5 lw 2 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 0.5 lw 2 lc rgb '0x000000'\n";

  stream << "plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + components[i].name + ".txt";
    stream << "    "
           << "\"" << fileName << "\""
           << " us ($1):($2)"
           << " title \"" << components[i].name << "\""
           << " with po" << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
}

void MixturePrediction::createMixturePlotScript()
{
  std::filesystem::create_directory("MixturePrediction");
  std::filesystem::create_directory(std::format("MixturePrediction/System_{}", system.systemId));
  std::ofstream stream(std::format("MixturePrediction/System_{}/plot_mixture", system.systemId));

  stream << "set encoding utf8\n";
  stream << "set xlabel 'Total bulk fluid phase fugacity, {/Helvetica-Italic f} / Pa' font \"Helvetica,18\"\n";
  stream << "set ylabel 'Absolute loading, {/Helvetica-Italic q}_i' offset 0.0,0 font \"Helvetica,18\"\n";
  stream << "set bmargin 4\n";
  if (system.pressure_range.scale == PressureRange::Scale::Log)
  {
    stream << "set key top left samplen 2.5 height 0.5 spacing 1.5 font \"Helvetica, 10\" maxcolumns 2\n";
    stream << "set log x\n";
    stream << "set format x \"10^{%T}\"\n";
    stream << "set xrange[" << system.pressure_range.pressureStart << ":]\n";
  }
  else
  {
    stream << "set key outside right samplen 2.5 height 0.5 spacing 1.5 font \"Helvetica, 10\" maxcolumns 2\n";
  }
  stream << "set key title '" << displayName << " {/:Italic T}=" << temperature << " K'\n";

  stream << "set output 'mixture_prediction.pdf'\n";
  stream << "set term pdf color solid\n";

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 0.5 lw 2 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 0.5 lw 2 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 0.5 lw 2 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 0.5 lw 2 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 0.5 lw 2 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 0.5 lw 2 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 0.5 lw 2 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 0.5 lw 2 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 0.5 lw 2 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 0.5 lw 2 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 0.5 lw 2 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 0.5 lw 2 lc rgb '0x000000'\n";

  stream << "plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + components[i].name + ".txt";
    stream << "    "
           << "\"" << fileName << "\""
           << " us ($1):($3)"
           << " title \"" << components[i].name << " (y_i=" << components[i].molFraction << ")\""
           << " with po" << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
}

void MixturePrediction::createMixtureAdsorbedMolFractionPlotScript()
{
  std::filesystem::create_directory("MixturePrediction");
  std::filesystem::create_directory(std::format("MixturePrediction/System_{}", system.systemId));
  std::ofstream stream(std::format("MixturePrediction/System_{}/plot_mixture_mol_fractions", system.systemId));

  stream << "set encoding utf8\n";
  stream << "set xlabel 'Total bulk fluid phase fugacity, {/Helvetica-Italic f} / Pa' font \"Helvetica,18\"\n";
  stream << "set ylabel 'Adsorbed mol-fraction, {/Helvetica-Italic Y}_i / [-]' offset 0.0,0 font \"Helvetica,18\"\n";
  stream << "set bmargin 4\n";
  if (system.pressure_range.scale == PressureRange::Scale::Log)
  {
    stream << "set key outside right samplen 2.5 height 0.5 spacing 1.5 font \"Helvetica, 10\" maxcolumns 2\n";
    stream << "set log x\n";
    stream << "set format x \"10^{%T}\"\n";
    stream << "set xrange[" << system.pressure_range.pressureStart << ":]\n";
  }
  else
  {
    stream << "set key outside right samplen 2.5 height 0.5 spacing 1.5 font \"Helvetica, 10\" maxcolumns 2\n";
  }
  stream << "set key title '" << displayName << " {/:Italic T}=" << temperature << " K'\n";

  stream << "set output 'mixture_prediction_mol_fractions.pdf'\n";
  stream << "set term pdf color solid\n";

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 0.5 lw 2 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 0.5 lw 2 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 0.5 lw 2 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 0.5 lw 2 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 0.5 lw 2 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 0.5 lw 2 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 0.5 lw 2 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 0.5 lw 2 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 0.5 lw 2 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 0.5 lw 2 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 0.5 lw 2 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 0.5 lw 2 lc rgb '0x000000'\n";

  stream << "plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + components[i].name + ".txt";
    stream << "    "
           << "\"" << fileName << "\""
           << " us ($1):($5)"
           << " title \"" << components[i].name << " (y_i=" << components[i].molFraction << ")\""
           << " with po" << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
}

void MixturePrediction::createPlotScript()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::filesystem::create_directory("MixturePrediction");
  std::filesystem::create_directory(std::format("MixturePrediction/System_{}", system.systemId));
  std::ofstream stream_graphs(std::format("MixturePrediction/System_{}/make_graphs.bat", system.systemId));
  stream_graphs << "set PATH=%PATH%;C:\\Program Files\\gnuplot\\bin;C:\\Program Files\\ffmpeg-master-"
                   "latest-win64-gpl\\bin;C:\\Program Files\\ffmpeg\\bin\n";
  stream_graphs << "gnuplot.exe plot_pure_components\n";
  stream_graphs << "gnuplot.exe plot_mixture\n";
  stream_graphs << "gnuplot.exe plot_mixture_mol_fractions\n";
  std::filesystem::path path{"make_graphs.bat"};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::filesystem::create_directory("MixturePrediction");
  std::filesystem::create_directory(std::format("MixturePrediction/System_{}", system.systemId));
  std::ofstream stream_graphs(std::format("MixturePrediction/System_{}/make_graphs", system.systemId));
  stream_graphs << "#!/bin/sh\n";
  stream_graphs << "cd -- \"$(dirname \"$0\")\"\n";
  stream_graphs << "gnuplot plot_pure_components\n";
  stream_graphs << "gnuplot plot_mixture\n";
  stream_graphs << "gnuplot plot_mixture_mol_fractions\n";
  std::filesystem::path path{std::format("MixturePrediction/System_{}/make_graphs", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif
}

void MixturePrediction::printErrorStatus(double psi_value, double sum, double P, const std::vector<double> Yi,
                                         double cachedP0[])
{
  std::cout << "psi: " << psi_value << std::endl;
  std::cout << "sum: " << sum << std::endl;
  for (size_t i = 0; i < Ncomp; ++i) std::cout << "cachedP0: " << cachedP0[i] << std::endl;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    double value = components[i].isotherm.inversePressureForPsi(psi_value, cachedP0[i]);
    std::cout << "inversePressure: " << value << std::endl;
  }
  std::cout << "P: " << P << std::endl;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    std::cout << "Yi[i] " << i << " " << Yi[i] << std::endl;
  }
}
