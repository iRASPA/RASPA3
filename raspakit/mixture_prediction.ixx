export module mixture_prediction;

import <vector>;
import <span>;
import <tuple>;
import <string>;
import <ostream>;

import atom;
import isotherm;
import multi_site_isotherm;
import component;
import system;

export struct MixturePrediction
{
  enum class PressureScale
  {
    Log = 0,
    Normal = 1
  };

  MixturePrediction(const System &system, double pressureStart, double pressureEnd,
                    size_t numberOfPressurePoints, MixturePrediction::PressureScale pressureScale);

  MixturePrediction(const System &system);

  void writeHeader(std::ostream &stream) const;

  // Yi  = gas phase mol-fraction
  // P   = total pressure
  // Xi  = adsorbed phase mol-fraction
  // Ni  = number of adsorbed molecules of component i
  std::pair<size_t, size_t> predictMixture(const std::vector<double> &Yi,
                                           const double &P,
                                           std::vector<double> &Xi,
                                           std::vector<double> &Ni,
                                           double *cachedP0,
                                           double &cachedPsi);
  void run(std::ostream &stream);
  void createPureComponentsPlotScript(std::string directoryNameString);
  void createMixturePlotScript(std::string directoryNameString);
  void createMixtureAdsorbedMolFractionPlotScript(std::string directoryNameString);

  const System &system;
  std::string displayName;
  std::span<const Component> components;
  std::vector<size_t> sortedComponentIndices;
  std::vector<std::reference_wrapper<const Component>> sortedComponents;
  const size_t Ncomp;
  Isotherm::MixturePredictionMethod predictionMethod;

  std::vector<double> alpha1;
  std::vector<double> alpha2;
  std::vector<double> alpha_prod;
  std::vector<double> x;

  double temperature{ 300.0 };
  double pressureStart{ 1e3 };
  double pressureEnd{ 1e7 };
  size_t numberOfPressurePoints{ 100 };
  PressureScale pressureScale{ PressureScale::Log };

  std::pair<size_t, size_t> computeIAST(const std::vector<double> &Yi,
                                    const double &P,
                                    std::vector<double> &Xi,
                                    std::vector<double> &Ni,
                                    double *cachedP0,
                                    double &cachedPsi);
  std::pair<size_t, size_t> computeExplicitLangmuir(const std::vector<double> &Yi,
                                                    const double &P,
                                                    std::vector<double> &Xi,
                                                    std::vector<double> &Ni);

  void printErrorStatus(double psi, double sum, double P, const std::vector<double> Yi, double cachedP0[]) const;
};

