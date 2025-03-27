module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#endif

export module isotherm_fitting;

#ifndef USE_LEGACY_HEADERS
import <tuple>;
import <array>;
import <vector>;
import <string>;
import <unordered_map>;
import <iostream>;
#endif

import randomnumbers;
import input_reader;
import component;
import multi_site_isotherm;
import system;

export struct IsothermFitting
{
  struct DNA
  {
    DNA(std::string genotype, MultiSiteIsotherm phenotype, double fitness)
        : genotype(genotype), phenotype(phenotype), fitness(fitness), hash(std::hash<std::string>{}(genotype))
    {
    }
    DNA() noexcept = default;
    DNA(const DNA &a) noexcept = default;
    DNA &operator=(const DNA &a) noexcept = default;
    DNA(DNA &&a) noexcept = default;
    DNA &operator=(DNA &&a) noexcept = default;
    ~DNA() noexcept = default;

    std::string genotype;
    MultiSiteIsotherm phenotype;
    double fitness;
    size_t hash;

    bool operator<(const DNA &other) const { return fitness < other.fitness; }
  };

  IsothermFitting(System &system) noexcept;

  System &system;
  RandomNumber random;
  std::vector<std::pair<double, double>> readData(size_t componentId);
  void printSolution(size_t Id);
  void run(std::ostream &stream);
  void createPlotScript(size_t componentId, const DNA &citizen);
  void createPlotScript();

  std::string writeHeader();
  void writeComponentIsothermFittingStatus(std::ostream &stream,
                                           const std::vector<std::pair<double, double>> &rawData) const;

  DNA newCitizen(size_t Id, const std::vector<std::pair<double, double>> &rawData);
  void updateCitizen(DNA &citizen, const std::vector<std::pair<double, double>> &rawData);
  double fitness(const MultiSiteIsotherm &phenotype, const std::vector<std::pair<double, double>> &rawData);
  double RCorrelation(const MultiSiteIsotherm &phenotype, const std::vector<std::pair<double, double>> &rawData);
  size_t biodiversity(size_t Id, const std::vector<DNA> &citizens);
  void nuclearDisaster(size_t Id, const std::vector<std::pair<double, double>> &rawData);
  void elitism();
  void mutate(DNA &Mutant, size_t Id);
  void crossover(size_t Id, size_t s1, size_t s2, size_t i1, size_t i2, size_t j1, size_t j2);
  void chooseRandomly(size_t kk1, size_t kk2, size_t jj1, size_t jj2, size_t &ii1, size_t &ii2);
  void mate(size_t Id, const std::vector<std::pair<double, double>> &rawData);
  void sortByFitness();
  void writeCitizen(std::ostream &stream, size_t componentId, size_t citizen, size_t step, size_t variety,
                    size_t fullfilledCondition, const std::vector<std::pair<double, double>> &rawData);
  DNA fit(std::ostream &stream, size_t Id, const std::vector<std::pair<double, double>> &rawData);
  const DNA simplex(std::ostream &stream, DNA citizen, double scale,
                    const std::vector<std::pair<double, double>> &rawData);
  void printSolution(std::ostream &stream, size_t componentId, const DNA &citizen);

  std::vector<MultiSiteIsotherm> isotherms;
  double maximumLoading{0.0};

  bool fittingFlag{false};
  bool physicalConstrainsFlag{false};
  bool seedFlag{false};
  bool pressureRangeFlag{false};
  bool refittingFlag{false};

  size_t GA_Size;             // population size
  double GA_MutationRate;     // mutation rate
  double GA_EliteRate;        // elitists population rate
  double GA_MotleyCrowdRate;  // pirates population rate
  double GA_DisasterRate;
  size_t GA_Elitists;   // number of elitists
  size_t GA_Motleists;  // number of pirates

  std::vector<DNA> popAlpha;
  std::vector<DNA> popBeta;
  std::vector<DNA> &parents;
  std::vector<DNA> &children;
};
