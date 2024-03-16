module;

module isotherm_fitting;

import <string>;
import <vector>;
import <iostream>;
import <fstream>;
import <filesystem>;
import <sstream>;
import <iterator>;
import <map>;
import <algorithm>;
import <exception>;
import <cmath>;
import <cstdlib>;
import <bitset>;
import <cstring>;
import <climits>;
import <unordered_set>;
import <chrono>;
import <print>;
import <optional>;

import randomnumbers;
import stringutils;
import special_functions;
import component;
import system;
import input_reader;
import isotherm;
import multi_site_isotherm;
import system;
import simulationbox;


IsothermFitting::IsothermFitting(System &system) noexcept:
  system(system),
  random(std::nullopt),
  isotherms(system.components.size()),
  GA_Size(static_cast<size_t>(std::pow(2.0, 12.0))),
  GA_MutationRate( 1.0/3.0 ),
  GA_EliteRate( 0.15 ),
  GA_MotleyCrowdRate( 0.25 ),
  GA_DisasterRate( 0.001 ),
  GA_Elitists( static_cast<size_t>(static_cast<double>(GA_Size) * GA_EliteRate) ),
  GA_Motleists( static_cast<size_t>(static_cast<double>(GA_Size) * (1.0 - GA_MotleyCrowdRate)) ),
  popAlpha(GA_Size),
  popBeta(GA_Size),
  parents(popAlpha),
  children(popBeta)
{
  for(size_t i = 0 ; i < system.components.size(); ++i)
  {
    isotherms[i] = system.components[i].isotherm;
  }
}

std::string IsothermFitting::writeHeader()
{
  std::ostringstream stream;

  for(size_t i = 0; i < system.components.size(); ++i)
  {
    std::print(stream, "Number of isotherm parameters: {}\n", system.components[i].isotherm.numberOfParameters);
    std::print(stream, "{}", system.components[i].isotherm.print());
    std::print(stream, "\n");
  }
  std::print(stream, "\n");

  return stream.str();
}

void IsothermFitting::writeComponentIsothermFittingStatus(std::ostream &stream, 
                                                          const std::vector<std::pair<double, double>> &rawData) const
{
  std::print(stream, "Found {} data points\n", rawData.size());
  for(const std::pair<double, double> &data : rawData)
  {
    std::print(stream, "pressure: {:.8e}  loading: {}\n", data.first, data.second);
  }
  std::print(stream, "\n");

  if(!rawData.empty())
  {
    std::pair<double, double> pressureRange = std::make_pair(rawData.front().first, rawData.back().first);
    std::print(stream,"Lowest pressure:     {:.8e}\n", pressureRange.first);
    std::print(stream,"Highest pressure:    {:.8e}\n", pressureRange.second);
  }
  std::print(stream, "\n\n");
}


void IsothermFitting::run(std::ostream &stream)
{
  for(size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if((system.components[componentId].type != Component::Type::Framework) &&
       (!system.components[componentId].filename.empty()))
    {
      const std::vector<std::pair<double, double>> rawData = readData(componentId);
      writeComponentIsothermFittingStatus(stream, rawData);

      const DNA bestCitizen = fit(stream, componentId, rawData);
      const DNA optimizedBestCitizen = simplex(stream, bestCitizen, 1.0, rawData);
      optimizedBestCitizen.phenotype.print();

      printSolution(stream, componentId, optimizedBestCitizen);

      createPlotScript(componentId, optimizedBestCitizen);
    }
  }
}

std::vector<std::pair<double, double>> IsothermFitting::readData(size_t componentId)
{
  std::string filename = system.components[componentId].filename;
  std::ifstream fileInput{ filename };
  std::string errorOpeningFile = "[IsothermFitting] File '" + filename + "' exists, but error opening file\n";
  if (!fileInput) throw std::runtime_error(errorOpeningFile);

  size_t columnPressure = system.components[componentId].columnPressure;
  size_t columnLoading = system.components[componentId].columnLoading;

  std::string line{};

  maximumLoading = 0.0;
  std::vector<std::pair<double, double>> rawData{};
  while (std::getline(fileInput, line))
  {
    std::string trimmedLine = trim(line);
    if(!startsWith(trimmedLine, "#"))
    {
      if (!line.empty())
      {
        std::istringstream iss(line);

        std::vector<std::string> results((std::istream_iterator<std::string>(iss)),
                                 std::istream_iterator<std::string>());
        if(columnPressure < results.size() &&
           columnLoading < results.size())
        {
          double pressure;
          double loading;
          std::istringstream s(results[columnPressure]);
          s >> pressure;
          std::istringstream t(results[columnLoading]);
          t >> loading;
          if(loading > maximumLoading)
          {
            maximumLoading = loading;
          }
          rawData.push_back(std::make_pair(pressure, loading));
        }
      }
    }
  }

  if(rawData.empty())
  {
    throw std::runtime_error("Error: no pressure points found\n");
  }

  // sort the pressures
  std::sort(rawData.begin(), rawData.end());

  return rawData;
}

// create a new citizen in the Ensemble
IsothermFitting::DNA IsothermFitting::newCitizen(size_t Id, const std::vector<std::pair<double, double>> &rawData)
{
    DNA citizen;

  citizen.phenotype = isotherms[Id].randomized(random, maximumLoading);

  citizen.genotype.clear();
  citizen.genotype.reserve((sizeof(double) * CHAR_BIT) *
                  citizen.phenotype.numberOfParameters);
  for(size_t i = 0; i < citizen.phenotype.numberOfParameters; ++i)
  {
    // convert from double to bitset
    uint64_t p;
    std::memcpy(&p, &citizen.phenotype.parameters(i), sizeof(double));
    std::bitset<sizeof(double) * CHAR_BIT> bitset(p);

    // add the bit-string to the genotype representation
    citizen.genotype += bitset.to_string();
  }

  citizen.hash = std::hash<MultiSiteIsotherm>{}(citizen.phenotype);
  citizen.fitness = fitness(citizen.phenotype, rawData);

  return citizen;
}

void IsothermFitting::updateCitizen(DNA &citizen, const std::vector<std::pair<double, double>> &rawData)
{
  citizen.fitness = fitness(citizen.phenotype, rawData);
}

// in latest clang, using -ffastmath optimizes isnan away
inline bool my_isnan(double val) 
{
  union { double f; uint64_t x; } u = { val };
  return (u.x << 1) > (0x7ff0000000000000u << 1);
}

double IsothermFitting::fitness(const MultiSiteIsotherm &phenotype, 
                                const std::vector<std::pair<double, double>> &rawData)
{
  double fitnessValue = phenotype.fitness();
  size_t m = rawData.size();              // number of observations
  size_t p = phenotype.numberOfParameters; // number of adjustable parameters
  for(std::pair<double, double> dataPoint: rawData)
  {
    double pressure = dataPoint.first;
    double loading = dataPoint.second;
    double difference = loading - phenotype.value(pressure);
    double weight = 1.0;
    fitnessValue += weight * difference * difference;
  }
  fitnessValue = sqrt(fitnessValue / static_cast<double>(m - p));

  if(my_isnan(fitnessValue)) fitnessValue = 99999999.999999;
  if(fitnessValue==0.0000000000) fitnessValue = 99999999.999999;

  return fitnessValue;
}

double IsothermFitting::RCorrelation(const MultiSiteIsotherm &phenotype, 
                                     const std::vector<std::pair<double, double>> &rawData)
{
  double RCorrelationValue = phenotype.fitness();
  size_t m = rawData.size();
  double loading_avg_o = 0.0;
  double loading_avg_e = 0.0;
  double tmp1 = 0.0;
  double tmp2 = 0.0;
  double tmp3 = 0.0;

  for(std::pair<double, double> dataPoint: rawData)
  {
    double pressure = dataPoint.first;
    double loading = dataPoint.second;
    loading_avg_o += loading / static_cast<double>(m);
    loading_avg_e += phenotype.value(pressure) / static_cast<double>(m);
  }

  for(std::pair<double, double> dataPoint: rawData)
  {
    double pressure = dataPoint.first;
    double loading = dataPoint.second;
    tmp1 += (loading-loading_avg_o)*(phenotype.value(pressure)-loading_avg_e);
    tmp2 += (loading-loading_avg_o)*(loading-loading_avg_o);
    tmp3 += (phenotype.value(pressure)-loading_avg_e)*(phenotype.value(pressure)-loading_avg_e);
  }
  RCorrelationValue = tmp1/sqrt(tmp2*tmp3);

  return RCorrelationValue;
}


size_t IsothermFitting::biodiversity([[maybe_unused]] size_t Id, const std::vector<DNA> &citizens)
{
  std::map<size_t, size_t> counts;
  for(const DNA &dna: citizens)
  {
    if(counts.find(dna.hash) != counts.end())
    {
      ++counts[dna.hash];
    }
    else
    {
      counts[dna.hash] = 1;
    }
  }
  size_t biodiversity = 0;
  for(const std::pair<size_t, size_t> value: counts)
  {
    if(value.second > 1)
    {
      biodiversity += value.second;
    }
  }

  return biodiversity;
}

void IsothermFitting::nuclearDisaster(size_t Id, const std::vector<std::pair<double, double>> &rawData)
{
  for(size_t i = 1; i < children.size(); ++i)
  {
    children[i] = newCitizen(Id, rawData);
  }
}

void IsothermFitting::elitism()
{
  std::copy(parents.begin(), parents.begin() + static_cast<std::vector<DNA>::difference_type>(GA_Elitists), 
            children.begin());
}

void IsothermFitting::mutate(DNA &mutant, [[maybe_unused]] size_t Id)
{
  mutant.genotype.clear();
  mutant.genotype.reserve((sizeof(double) * CHAR_BIT) * mutant.phenotype.numberOfParameters);
  for(size_t i = 0; i < mutant.phenotype.numberOfParameters; ++i)
  {
    // convert from double to bitset
    uint64_t p;
    std::memcpy(&p, &mutant.phenotype.parameters(i), sizeof(double));
    std::bitset<sizeof(double) * CHAR_BIT> bitset(p);

    // mutation: randomly flip bit
    bitset.flip(std::size_t((sizeof(double) * CHAR_BIT) * random.uniform()));

    // convert from bitset to double
    p = bitset.to_ullong();
    std::memcpy(&mutant.phenotype.parameters(i), &p, sizeof(double));

    // add the bit-string to the genotype representation
    mutant.genotype += bitset.to_string();
  }

  // calculate the hah-value from the entire bit-string
  mutant.hash = std::hash<MultiSiteIsotherm>{}(mutant.phenotype);
}

// [s1:s2] range of children
// [i1:i2] range of parent1
// [j1:j2] range of parent2
//----------------------------------------------
//  parent1    parent2                children
//    *          *                      *
//  00|000000  11|111111        ->    00|111111
//----------------------------------------------
void IsothermFitting::crossover(size_t Id, size_t s1,size_t s2, size_t i1, size_t i2, size_t j1, size_t j2)
{
    size_t k1,k2;
  double  tmp1;
  for(size_t i = s1; i < s2; ++i)
  {
    chooseRandomly(i1, i2, j1, j2, k1, k2);
    tmp1 = random.uniform();
    // choose between single cross-over using bit-strings or random parameter-swap
    if(tmp1 < 0.490)
      // One-point crossover:
      // --------------------
    {
      // remove the extreme values 0 and 32*Npar - 1 (they are not valid for crossover)
      size_t bitStringSize = (sizeof(double) * CHAR_BIT) * isotherms[Id].numberOfParameters;
      size_t spos = random.integer(1, bitStringSize - 2);
      children[i].genotype = parents[k1].genotype.substr(0, spos) +
                             parents[k2].genotype.substr(spos, bitStringSize - spos);

      // convert the bit-strings to doubles
      for(size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
      {
        size_t pos = j * (sizeof(double) * CHAR_BIT);
        size_t size = sizeof(double) * CHAR_BIT;
        std::bitset<sizeof(double) * CHAR_BIT> bitset(children[i].genotype, pos, size);
        uint64_t p = bitset.to_ullong();
        std::memcpy(&children[i].phenotype.parameters(j), &p, sizeof(double));
      }
    }
    else if ( tmp1 < 0.499)
      // Two-point crossover:
      // --------------------
      {
      size_t bitStringSize = (sizeof(double) * CHAR_BIT) * isotherms[Id].numberOfParameters;
      size_t spos1 = random.integer(1, bitStringSize - 3);
      size_t spos2 = random.integer(spos1, bitStringSize - 2);
      children[i].genotype = parents[k1].genotype.substr(0, spos1) +
                             parents[k2].genotype.substr(spos1, spos2 - spos1) +
                             parents[k1].genotype.substr(spos2, bitStringSize - spos2);
      // convert the bit-strings to doubles
      for(size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
      {
        size_t pos = j * (sizeof(double) * CHAR_BIT);
        size_t size = sizeof(double) * CHAR_BIT;
        std::bitset<sizeof(double) * CHAR_BIT> bitset(children[i].genotype, pos, size);
        uint64_t p = bitset.to_ullong();
        std::memcpy(&children[i].phenotype.parameters(j), &p, sizeof(double));
      }
    }
    else if (tmp1 < 0.500)
    {
      // Uniform crossover:
      // ------------------
      size_t bitStringSize = (sizeof(double) * CHAR_BIT) * isotherms[Id].numberOfParameters;
      size_t rolling_k = k1;
      for (size_t j = 0; j < bitStringSize; j++)
      {
        if(random.uniform() < 0.25 )
        {
           if(rolling_k == k1)
           {
             rolling_k = k2;
           }
           else
           {
             rolling_k = k1;
           }
        }
        children[i].genotype.substr(j,1) = parents[rolling_k].genotype.substr(j,1);
      }
      // convert the bit-strings to doubles
      for(size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
      {
        size_t pos = j * (sizeof(double) * CHAR_BIT);
        size_t size = sizeof(double) * CHAR_BIT;
        std::bitset<sizeof(double) * CHAR_BIT> bitset(children[i].genotype, pos, size);
        uint64_t p = bitset.to_ullong();
        std::memcpy(&children[i].phenotype.parameters(j), &p, sizeof(double));
      }
    }
    else
    {
      children[i].genotype.clear();
      for(size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
      {
        // randomly choose whether the parameter comes from parent k1 or k2
        if(random.uniform() < 0.5)
        {
          children[i].phenotype.parameters(j) = parents[k1].phenotype.parameters(j);
        }
        else
        {
          children[i].phenotype.parameters(j) = parents[k2].phenotype.parameters(j);
        }

        // convert from double to bitString
        uint64_t p;
        std::memcpy(&p, &children[i].phenotype.parameters(j), sizeof(double));
        std::bitset<sizeof(double) * CHAR_BIT> bitset(p);
        children[i].genotype += bitset.to_string();
      }
    }

    children[i].hash = std::hash<MultiSiteIsotherm>{}(children[i].phenotype);
  }
}

void IsothermFitting::chooseRandomly(size_t kk1,size_t kk2,size_t jj1,size_t jj2, size_t &ii1, size_t &ii2)
{
  ii1  = random.integer(kk1, kk2);
  ii2  = random.integer(jj1, jj2);
  while ( ii1 == ii2 )
  {
    ii2 = random.integer(jj1,jj2);
  };
}

void IsothermFitting::mate(size_t Id, const std::vector<std::pair<double, double>> &rawData)
{
    // retain the first 25% of the children
  elitism();

  // mates from GA_Elitists to (GA_Size - GA_Elitists)
  crossover(Id, GA_Elitists, GA_Elitists + static_cast<size_t>(static_cast<double>(GA_Motleists) * 0.5), 0, 
            GA_Elitists, 0, GA_Elitists);
  crossover(Id, GA_Elitists + static_cast<size_t>(static_cast<double>(GA_Motleists) * 0.5) + 1, 
            GA_Size - GA_Elitists, 0, GA_Elitists, GA_Elitists, GA_Size - 1);

  // mutation from GA_Elitists to (GA_Size - GA_Elitists) with "GA_MutationRate" probability
  for(size_t i = GA_Elitists; i < GA_Size - GA_Elitists; ++i)
  {
    if(random.uniform() < GA_MutationRate)
    {
      mutate(children[i], Id);
    }
    updateCitizen(children[i], rawData);
  }

  // replace the last GA_Elitists (the worst) of the children by new children
  for(size_t i = GA_Size - GA_Elitists; i < GA_Size; ++i)
  {
    children[i] = newCitizen(Id, rawData);
  }

  // replace the last (GA_Size - 1) children by new children
  if(random.uniform() < GA_DisasterRate)
  {
    nuclearDisaster(Id, rawData);
  }
}

void IsothermFitting::sortByFitness()
{
  std::sort(parents.begin(), parents.end());
}


void IsothermFitting::writeCitizen(std::ostream &stream, size_t componentId, size_t citizen, size_t step, 
                           size_t variety, size_t fullfilledCondition, 
                           const std::vector<std::pair<double, double>> &rawData)
{
  if(fullfilledCondition > 0)
  {
    std::print(stream, 
        "mol: {:2d} step: {:5d} Fitness: {:8.5f} R^2: {:8.5f} Similarity: {:5d}/{:<5d} Finishing: {:3d}/{:<3d}\n",
        componentId, step, parents[citizen].fitness, std::pow(RCorrelation(parents[citizen].phenotype, rawData),2), 
        variety, GA_Size, fullfilledCondition, 100);
  }
  else 
  {
    std::print(stream,
         "mol: {:2d} step: {:5d} Fitness: {:8.5f} R^2: {:8.5f} Similarity: {:2d}/{:<5d}\n", 
         componentId, step, parents[citizen].fitness, std::pow(RCorrelation(parents[citizen].phenotype, rawData),2), 
         variety, GA_Size);
  }
  std::print(stream, "number of parameters: {}\n", parents[citizen].phenotype.numberOfParameters);
  for(size_t i = 0; i < parents[citizen].phenotype.numberOfParameters; ++i)
  { 
    std::print(stream, "      genotype: {}  parameter: {}\n",
               parents[citizen].genotype.substr(64*i,64),
               parents[citizen].phenotype.parameters(i));
  }
  std::print(stream, "\n");
}


IsothermFitting::DNA IsothermFitting::fit(std::ostream &stream, size_t componentId, 
                                          const std::vector<std::pair<double, double>> &rawData)
{
  size_t optimisationStep{ 0 };
  const size_t maxOptimisationStep{ 1000 };
  size_t fullFilledConditionStep{ 0 };
  const size_t maxFullfilledConditionStep{ 100 };
  const double minimumFitness{ 5.0e-1 };
  double tempFitnessValue{ 999.0 };
  size_t tempVarietyValue{ 0 };
  const double toleranceEqualFitness{ 1e-3 };
  const size_t minstep{ 10 };

  fullFilledConditionStep = 0;
  optimisationStep = 0;

  for(size_t i = 0; i < popAlpha.size(); ++i)
  {
    popAlpha[i] = newCitizen(componentId, rawData);
    popBeta[i] = newCitizen(componentId, rawData);
  }

  parents = popAlpha;
  children = popBeta;

  sortByFitness();

  // print Initial (and unsorted) population
  writeCitizen(stream, componentId, 0, 0, 0, 0, rawData);

  if(refittingFlag)
  {
    std::cout << "Refitting activated\n";
    for(size_t citizen = 0; citizen < 2; ++citizen)
    {
      parents[citizen].genotype.clear();
      parents[citizen].genotype.reserve((sizeof(double) * CHAR_BIT) *
                      parents[citizen].phenotype.numberOfParameters);
      for(size_t i = 0; i < parents[citizen].phenotype.numberOfParameters; ++i)
      {
        // convert from double to bitset
        uint64_t p;
        std::memcpy(&p, &parents[citizen].phenotype.parameters(i), sizeof(double));
        std::bitset<sizeof(double) * CHAR_BIT> bitset(p);

        // add the bit-string to the genotype representation
        parents[citizen].genotype += bitset.to_string();
      }
      parents[citizen].hash = std::hash<MultiSiteIsotherm>{}(parents[citizen].phenotype);
      updateCitizen(parents[citizen], rawData);
    }
    std::copy(parents.begin(), parents.end(), children.begin());

    mate(componentId, rawData);

    std::swap(parents, children);

    sortByFitness();
  }

  isotherms[componentId] = parents[0].phenotype;


  std::cout << "Starting Genetic Algorithm optimization\n";

  bool continueCondition = true;
  do
  {
    sortByFitness();

    tempVarietyValue = biodiversity(componentId, children);

    writeCitizen(stream, componentId, 0, optimisationStep,tempVarietyValue, fullFilledConditionStep, rawData);

    if( optimisationStep >= minstep &&
        parents[0].fitness <= minimumFitness &&
        std::abs(parents[0].fitness - tempFitnessValue) <= toleranceEqualFitness )
    {
      fullFilledConditionStep += 1;
    }
    else
    {
      fullFilledConditionStep = 0;
    }

    if (optimisationStep >= maxOptimisationStep ||
       fullFilledConditionStep >= maxFullfilledConditionStep )
    {
      continueCondition = false;
    }

    // pairing
    mate(componentId, rawData);

    std::swap(parents, children);

    // take the best fitness value:
    tempFitnessValue = parents[0].fitness;

    optimisationStep += 1;

  }while(continueCondition);

  writeCitizen(stream, componentId, 0, optimisationStep,tempVarietyValue, fullFilledConditionStep, rawData);

  return parents[0];
}

// The Nelder-Mead method uses a simplex (a hyper-tetrahedron of n+1 vertices in n dimensions)
// Advantage: it does not use derivates, works well, can tolerate some noise
// Disadvantage: it is not garanteed to converge
// The steps of the method are:
// 1) Sort: according to the fitness
// 2) Reflect: get rid of the worst point, replace it by something better
// 3) Extend: if better then extend it even further
// 4) Contract: if not better then contract it 
// 5) Shrink: if still not better we shrink towards the best performing point
// 6) Check convergence
const IsothermFitting::DNA IsothermFitting::simplex(std::ostream &stream, DNA citizen, double scale, 
                                                    const std::vector<std::pair<double, double>> &rawData)
{
  size_t n = citizen.phenotype.numberOfParameters;
  std::vector<std::vector<double>> v(n + 1, std::vector<double>(n)); // holds vertices of simplex
  std::vector<double> f(n + 1);  // value of function at each vertex
  std::vector<double> vr(n);     // reflection - coordinates
  std::vector<double> ve(n);     // expansion - coordinates
  std::vector<double> vc(n);     // contraction - coordinates
  std::vector<double> vm(n);     // centroid - coordinates
  std::vector<double> vtmp(n);   // temporary array passed to FUNC
  size_t vs;                     // vertex with the smallest value
  size_t vh;                     // vertex with next smallest value
  size_t vg;                     // vertex with largest value
  double fr;                     // value of function at reflection point
  double fe;                     // value of function at expansion point
  double fc;                     // value of function at contraction point
  size_t iprint{ 0 };
  size_t MAX_IT = 1000000;
  double EPSILON = 1.0e-4;
  double ALPHA = 1.0;
  double BETA = 0.5;
  double GAMMA = 2.0;

  std::print(stream, "\nMinimising the cost function using the Nelder-Mead SIMPLEX method:\n\n");

  for(size_t i = 0 ; i < n; ++i)
  {
    v[0][i] = citizen.phenotype.parameters(i);
  }

  // values used to create initial simplex
  double pn = scale * (std::sqrt(static_cast<double>(n) + 1.0) - 1.0 + static_cast<double>(n)) / 
              (static_cast<double>(n) * sqrt(2.0));
  double qn = scale * (std::sqrt(static_cast<double>(n) + 1.0) - 1.0) / (static_cast<double>(n) * sqrt(2.0));
  for(size_t i = 1; i <= n; ++i)
  {
    for(size_t j = 0; j < n; ++j)
    {
      if (i - 1 == j)
      {
        v[i][j] = pn + citizen.phenotype.parameters(j);
      }
      else 
      {
        v[i][j] = qn + citizen.phenotype.parameters(j);
      }
    }
  }

  for(size_t i = 0; i <= n; ++i)
  {
    for(size_t j = 0; j < n; ++j)
    {
      citizen.phenotype.parameters(j) = v[i][j];
    }
    f[i] = fitness(citizen.phenotype, rawData);
  }

  // print out the initial simplex
  // print out the initial function values
  // find the index of the smallest value for printing
  if (iprint == 0)
  {
    vs = 0;
    for(size_t j = 0; j <= n; ++j)
    {
      if (f[j] < f[vs])
      {
        vs = j;
      }
    }

    std::print(stream, "Initial Values from genetic algorithm:\n");

    for(size_t j = 0; j < n; ++j)
    {
      std::print(stream, "    {}\n",v[vs][j]);
    }
    std::print(stream, "Fit: {}\n\n", f[vs]);
  }

  for(size_t itr = 1; itr <= MAX_IT; ++itr)
  {
    // Step 1: Sort
    // ====================================================================
    std::vector<size_t> sortIndexes = sort_indexes(f);
    vs = sortIndexes[0];   // index of smallest
    vg = sortIndexes[n];   // index of largest
    vh = sortIndexes[n-1]; // index of second largest
  
    // calculate the center point of every point except for the worst one
    for(size_t j = 0; j < n; ++j)
    {
      double cent = 0.0;
      for(size_t i = 0; i <= n; ++i)
      {
        if (i != vg)
        {
          cent = cent + v[i][j];
        }
      }
      vm[j] = cent / static_cast<double>(n);
    }

    // Step 2: Reflect vg to new vertex vr
    // ====================================================================
    for(size_t j = 0; j < n; ++j)
    {
      vr[j] = (1.0 + ALPHA) * vm[j] - ALPHA * v[vg][j];
      citizen.phenotype.parameters(j) = vr[j];
    }
    fr = fitness(citizen.phenotype, rawData);

    if ((fr <= f[vh]) && (fr > f[vs]))
    {
      for(size_t j = 0; j < n; ++j)
      {
        v[vg][j] = vr[j];
      }
      f[vg] = fr;
    }

    // Step 3: Extend a step further in this direction
    // ====================================================================
    if (fr <= f[vs])
    {
      for(size_t j = 0; j < n; ++j)
      {
        ve[j] = GAMMA * vr[j] + (1.0 - GAMMA) * vm[j];
        citizen.phenotype.parameters(j) = ve[j];
      }
      fe = fitness(citizen.phenotype, rawData);

      // by making fe < fr as opposed to fe < f(vs), Rosenbrocks function
      // takes 62 iterations as opposed to 64.

      if (fe < fr)
      {
        for(size_t j = 0; j < n; ++j)
        {
          v[vg][j] = ve[j];
        }
        f[vg] = fe;
      }
      else
      {
        for(size_t j = 0; j < n; ++j)
        {
          v[vg][j] = vr[j];
        }
        f[vg] = fr;
      }
    }

    if (fr > f[vh])
    {
      // Step 4: Contraction
      // ====================================================================
      for(size_t j = 0; j < n; ++j)
      {
        vc[j] = BETA * v[vg][j] + (1.0 - BETA) * vm[j];
        citizen.phenotype.parameters(j) = vc[j];
      }
      fc = fitness(citizen.phenotype, rawData);
      if (fc < f[vg])
      {
        for(size_t j = 0; j < n; ++j)
        {
          v[vg][j] = vc[j];
        }
        f[vg] = fc;
      }
      else
      {
        // Step 4: Shrink
        // ====================================================================
        // at this point the contraction is not successful,
        // we must halve the distance from vs to all the
        // vertices of the simplex and then continue.
        for(size_t row=0; row <= n; ++row)
        {
          if (row != vs)
          {
            for(size_t j=0; j < n; ++j)
            {
              v[row][j] = v[vs][j] + 0.5 * (v[row][j] - v[vs][j]);
            }
          }
        }
        for(size_t m = 0; m < n; ++m)
        {
          vtmp[m] = v[vg][m];
          citizen.phenotype.parameters(m) = vtmp[m];
        }
        f[vg] = fitness(citizen.phenotype, rawData);

        for(size_t m = 0; m < n; ++m)
        {
          vtmp[m] = v[vh][m];
          citizen.phenotype.parameters(m) = vtmp[m];
        }
        f[vh] = fitness(citizen.phenotype, rawData);
      }
    }

    // Step 6: Test for convergence
    // ====================================================================
    double fsum = 0.0;
    for(size_t j = 0; j <= n; ++j)
    {
      fsum = fsum + f[j];
    }
    double favg = fsum / (static_cast<double>(n) + 1.0);

    if (favg < EPSILON || itr == MAX_IT)
    {
      // print out the value at each iteration

      if(itr != MAX_IT)
      {
        std::print(stream, "Nelder-Mead has converged: {} < {} \n\n", favg, EPSILON);
      }
      else
      {
        std::print(stream, "Reached maximum number of steps: {} = {}\n\n", itr, MAX_IT);
      }

      // find the index of the smallest value
      sortIndexes = sort_indexes(f);
      vs = sortIndexes[0];   // index of smallest
      for(size_t m = 0; m < n; ++m)
      {
        citizen.phenotype.parameters(m) = v[vs][m];
      }
      double min = fitness(citizen.phenotype, rawData);

      std::print(stream, "Final Values: \n");
      for(size_t j = 0; j < n; ++j)
      {
        std::print(stream, "    {}\n", v[vs][j]);
      }
      std::print(stream, "Fit: {}\n\n", min);

      return citizen;
    }
  }

  return citizen;
}

void IsothermFitting::printSolution(std::ostream &stream,  size_t componentId, const DNA &citizen)
{
  const Component &component = system.components[componentId];
  std::print(stream, "Component {:2d} MoleculeName             {}\n", componentId, component.name);
  std::print(stream, "             FileName                 {}\n", component.filename);
  std::print(stream, "             ColumnPressure           {:<3d}\n", component.columnPressure);
  std::print(stream, "             ColumnLoading            {:<3d}\n", component.columnLoading);
  std::print(stream, "             NumberOfIsothermSites    {:<3d}\n", component.isotherm.sites.size());
  std::print(stream, "{}\n\n", citizen.phenotype.printAsInputFormat());
}

void IsothermFitting::createPlotScript(size_t componentId, const DNA &citizen)
{
  std::string rawDataFileName = system.components[componentId].filename;
  std::string componentName = system.components[componentId].name;
  double T = system.temperature;

  std::filesystem::path directoryName = std::format("IsothermFitting/System_{}/", system.systemId);
  std::filesystem::path plotFileName = std::format("IsothermFitting/System_{}/{}", system.systemId,
      "plot_fit_component_" + std::to_string(componentId) + "_" + componentName);
  std::filesystem::create_directories(directoryName);
  std::filesystem::copy(rawDataFileName, std::format("IsothermFitting/System_{}/{}", system.systemId, rawDataFileName),
                        std::filesystem::copy_options::update_existing);

  std::ofstream stream(plotFileName);

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set xlabel 'Pressure, {{/Helvetica-Italic p}} / [Pa]' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Absolute loading, {{/Helvetica-Italic q}} / [mol/kg]' offset 0.0,0 "
                     "font 'Helvetica,18'\n");
  std::print(stream, "set bmargin 4\n");
  std::print(stream, "set yrange[0:]\n");
  if(system.components[componentId].pressureScale == Component::PressureScale::Log)
  {
    std::print(stream, "set log x\n");
  }

  std::print(stream, "set key  right bottom vertical samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n");
  std::print(stream, "set key title '{} ({}) {{/:Italic T}}={} K'\n", 
             componentName, system.components.front().name, T);

  std::print(stream, "set output 'isotherms_fit_{}.pdf'\n", componentName);
  std::print(stream, "set term pdf color solid\n");

  std::print(stream, "array s[{}]\n", system.components[componentId].isotherm.numberOfParameters);
  for(size_t i = 0; i < system.components[componentId].isotherm.numberOfParameters; ++i)
  {
    std::print(stream, "s[{}]={}\n", i+1, system.components[componentId].isotherm.parameters(i));
  }
  std::print(stream, "array p[{}]\n", citizen.phenotype.numberOfParameters);
  for(size_t i = 0; i < citizen.phenotype.numberOfParameters; ++i)
  {
    std::print(stream, "p[{}]={}\n", i + 1, citizen.phenotype.parameters(i));
  }
  std::print(stream, "plot \\\n");
  std::print(stream, "{} title 'start f(x)' with li dt 2 lw 2,\\\n",
             system.components[componentId].isotherm.gnuplotFunctionString('s'));
  std::print(stream, "{} title 'fit f(x)' with li lw 2,\\\n",
             citizen.phenotype.gnuplotFunctionString('p'));
  size_t columnPressure = system.components[componentId].columnPressure + 1;
  size_t columnLoading = system.components[componentId].columnLoading + 1;
  std::print(stream, "'{}' us {}:{} title 'raw data' with po pt 5 ps 0.5\n", rawDataFileName, columnPressure, columnLoading);
}

void IsothermFitting::createPlotScript()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::filesystem::create_directory("IsothermFitting");
  std::filesystem::create_directory(std::format("IsothermFitting/System_{}", system.systemId));
  std::ofstream stream_graphs(std::format("IsothermFitting/System_{}/make_graphs.bat", system.systemId));
  stream_graphs << "set PATH=%PATH%;C:\\Program Files\\gnuplot\\bin;C:\\Program Files\\ffmpeg-master-latest-"
                   "win64-gpl\\bin;C:\\Program Files\\ffmpeg\\bin\n";
  for(size_t i = 0; i < system.components.size(); ++i)
  { 
    if(system.components[i].isotherm.numberOfParameters > 0)
    {
      std::print(stream_graphs, "gnuplot.exe plot_fit_component_{}_{}\n", 
                                std::to_string(i), system.components[i].name);
    }
  }
  std::filesystem::path path{ "make_graphs.bat" };
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::filesystem::create_directory("IsothermFitting");
  std::filesystem::create_directory(std::format("IsothermFitting/System_{}", system.systemId));
  std::ofstream stream_graphs(std::format("IsothermFitting/System_{}/make_graphs", system.systemId));
  stream_graphs << "#!/bin/sh\n";
  stream_graphs << "cd -- \"$(dirname \"$0\")\"\n";
  for (size_t i = 0; i < system.components.size(); ++i)
  {
    if (system.components[i].isotherm.numberOfParameters > 0)
    {
      std::print(stream_graphs, "gnuplot.exe plot_fit_component_{}_{}\n", 
                                std::to_string(i), system.components[i].name);
    }
  }
  std::filesystem::path path{ "make_graphs" };
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif
}
