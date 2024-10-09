module;

#ifdef USE_LEGACY_HEADERS
#include <functional>
#include <ostream>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#endif

export module mixture_prediction;

#ifndef USE_LEGACY_HEADERS
import <functional>;
import <vector>;
import <span>;
import <tuple>;
import <string>;
import <ostream>;
#endif

import atom;
import isotherm;
import multi_site_isotherm;
import component;
import system;
import bond_potential;

export class MixturePrediction
{
 public:
  enum class IASTMethod
  {
    FastIAST = 0,
    NestedLoopBisection = 1
  };

  MixturePrediction(const System &system);

  std::string writeHeader() const;

  void print() const;
  void run(std::ostream &stream);
  void createPureComponentsPlotScript();
  void createMixturePlotScript();
  void createMixtureAdsorbedMolFractionPlotScript();
  void createPlotScript();

  // Yi  = gas phase mol-fraction
  // P   = total pressure
  // Xi  = adsorbed phase mol-fraction
  // Ni  = number of adsorbed molecules of component i
  std::pair<size_t, size_t> predictMixture(const std::vector<double> &Yi, const double &P, std::vector<double> &Xi,
                                           std::vector<double> &Ni, double *cachedP0, double *cachedPsi);

 private:
  const System &system;
  std::string displayName;
  std::vector<Component> components;
  std::vector<size_t> sortedComponentIndices;
  std::vector<std::reference_wrapper<const Component>> sortedComponents;
  const size_t Ncomp;
  const size_t Nsorted;
  // size_t numberOfCarrierGases{ 0 };
  // size_t carrierGasComponent{ 0 };
  MultiSiteIsotherm::PredictionMethod predictionMethod;
  IASTMethod iastMethod;
  size_t maxIsothermTerms{2};
  std::vector<std::vector<Component>> segregatedSortedComponents;

  std::vector<double> alpha1;
  std::vector<double> alpha2;
  std::vector<double> alpha_prod;
  std::vector<double> x;

  std::vector<double> pstar;
  std::vector<double> psi;
  std::vector<double> G;
  std::vector<double> delta;
  std::vector<double> Phi;

  double temperature{300.0};

  std::pair<size_t, size_t> computeFastIAST(const std::vector<double> &Yi, const double &P, std::vector<double> &Xi,
                                            std::vector<double> &Ni, double *cachedP0, double *cachedPsi);
  std::pair<size_t, size_t> computeFastSIAST(const std::vector<double> &Yi, const double &P, std::vector<double> &Xi,
                                             std::vector<double> &Ni, double *cachedP0, double *cachedPsi);
  std::pair<size_t, size_t> computeFastSIAST(size_t term, const std::vector<double> &Yi, const double &P,
                                             std::vector<double> &Xi, std::vector<double> &Ni, double *cachedP0,
                                             double *cachedPsi);

  std::pair<size_t, size_t> computeIASTNestedLoopBisection(const std::vector<double> &Yi, const double &P,
                                                           std::vector<double> &Xi, std::vector<double> &Ni,
                                                           double *cachedP0, double *cachedPsi);
  std::pair<size_t, size_t> computeSIASTNestedLoopBisection(const std::vector<double> &Yi, const double &P,
                                                            std::vector<double> &Xi, std::vector<double> &Ni,
                                                            double *cachedP0, double *cachedPsi);
  std::pair<size_t, size_t> computeSIASTNestedLoopBisection(size_t term, const std::vector<double> &Yi, const double &P,
                                                            std::vector<double> &Xi, std::vector<double> &Ni,
                                                            double *cachedP0, double *cachedPsi);
  std::pair<size_t, size_t> computeExplicitIsotherm(const std::vector<double> &Yi, const double &P,
                                                    std::vector<double> &Xi, std::vector<double> &Ni);
  std::pair<size_t, size_t> computeSegratedExplicitIsotherm(const std::vector<double> &Yi, const double &P,
                                                            std::vector<double> &Xi, std::vector<double> &Ni);
  std::pair<size_t, size_t> computeSegratedExplicitIsotherm(size_t site, const std::vector<double> &Yi, const double &P,
                                                            std::vector<double> &Xi, std::vector<double> &Ni);

  void printErrorStatus(double psi, double sum, double P, const std::vector<double> Yi, double cachedP0[]);
};
