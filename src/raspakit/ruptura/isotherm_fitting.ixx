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
import std;
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
    std::size_t hash;

    bool operator<(const DNA &other) const { return fitness < other.fitness; }
  };

  IsothermFitting(System &system) noexcept;

  System &system;
  RandomNumber random;
  std::vector<std::pair<double, double>> readData(std::size_t componentId);
  void printSolution(std::size_t Id);
  void run(std::ostream &stream);
  void createPlotScript(std::size_t componentId, const DNA &citizen);
  void createPlotScript();

  std::string writeHeader();
  void writeComponentIsothermFittingStatus(std::ostream &stream,
                                           const std::vector<std::pair<double, double>> &rawData) const;

  DNA newCitizen(std::size_t Id, const std::vector<std::pair<double, double>> &rawData);
  void updateCitizen(DNA &citizen, const std::vector<std::pair<double, double>> &rawData);
  double fitness(const MultiSiteIsotherm &phenotype, const std::vector<std::pair<double, double>> &rawData);
  double RCorrelation(const MultiSiteIsotherm &phenotype, const std::vector<std::pair<double, double>> &rawData);
  std::size_t biodiversity(std::size_t Id, const std::vector<DNA> &citizens);
  void nuclearDisaster(std::size_t Id, const std::vector<std::pair<double, double>> &rawData);
  void elitism();
  void mutate(DNA &Mutant, std::size_t Id);
  void crossover(std::size_t Id, std::size_t s1, std::size_t s2, std::size_t i1, std::size_t i2, std::size_t j1,
                 std::size_t j2);
  void chooseRandomly(std::size_t kk1, std::size_t kk2, std::size_t jj1, std::size_t jj2, std::size_t &ii1,
                      std::size_t &ii2);
  void mate(std::size_t Id, const std::vector<std::pair<double, double>> &rawData);
  void sortByFitness();
  void writeCitizen(std::ostream &stream, std::size_t componentId, std::size_t citizen, std::size_t step,
                    std::size_t variety, std::size_t fullfilledCondition,
                    const std::vector<std::pair<double, double>> &rawData);
  DNA fit(std::ostream &stream, std::size_t Id, const std::vector<std::pair<double, double>> &rawData);
  const DNA simplex(std::ostream &stream, DNA citizen, double scale,
                    const std::vector<std::pair<double, double>> &rawData);
  void printSolution(std::ostream &stream, std::size_t componentId, const DNA &citizen);

  std::vector<MultiSiteIsotherm> isotherms;
  double maximumLoading{0.0};

  bool fittingFlag{false};
  bool physicalConstrainsFlag{false};
  bool seedFlag{false};
  bool pressureRangeFlag{false};
  bool refittingFlag{false};

  std::size_t GA_Size;        // population size
  double GA_MutationRate;     // mutation rate
  double GA_EliteRate;        // elitists population rate
  double GA_MotleyCrowdRate;  // pirates population rate
  double GA_DisasterRate;
  std::size_t GA_Elitists;   // number of elitists
  std::size_t GA_Motleists;  // number of pirates

  std::vector<DNA> popAlpha;
  std::vector<DNA> popBeta;
  std::vector<DNA> &parents;
  std::vector<DNA> &children;
};
