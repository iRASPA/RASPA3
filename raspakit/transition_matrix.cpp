module;

module transition_matrix;

import <cmath>;
import <cstddef>;
import <vector>;
import <algorithm>;
import <utility>;
import <filesystem>;
import <fstream>;

import double3;

void TransitionMatrix::initialize()
{
  if(!doTMMC) return;

  cmatrix.resize(maxMacrostate - minMacrostate + 1);
  bias.resize(maxMacrostate - minMacrostate + 1); 
  std::fill(bias.begin(), bias.end(), 1.0);
  lnpi.resize(maxMacrostate - minMacrostate + 1);
  forward_lnpi.resize(maxMacrostate - minMacrostate + 1);
  reverse_lnpi.resize(maxMacrostate - minMacrostate + 1);
  histogram.resize(maxMacrostate - minMacrostate + 1);
}

void TransitionMatrix::update(double3 Pacc, size_t newN, size_t oldN)
{
  if(!doTMMC) return;

  Pacc.clamp(0.0, 1.0);

  if(rejectOutofBound && ((newN > maxMacrostate) || (newN < minMacrostate))) return;

  cmatrix[oldN - minMacrostate] += Pacc;
  histogram[newN - minMacrostate]++;
};

// return the biasing Factor
double TransitionMatrix::biasFactor(size_t newN, size_t oldN)
{
  if(!doTMMC || !useBias || !useTMBias) return 1.0;

  double TMMCBias = bias[newN - minMacrostate] - bias[oldN - minMacrostate];
  return std::exp(TMMCBias); 
};

//From Vince Shen's pseudo code//
void TransitionMatrix::adjustBias()
{
  if(!doTMMC || !useBias || !useTMBias) return;

  if((numberOfSteps % updateTMEvery != 0) || numberOfSteps == 0) return;

  printf("Adjusting bias\n");
  numberOfUpdates++;

  //First step is to get the lowest and highest visited states in terms of loading//
  size_t minVisited = 0; 
  size_t maxVisited = 0;
  size_t nonzeroCount=0;
  //Zhao's special note: length of the vectors for TMMC = maxMacrostate - minMacrostate + 1;
  //The a, minVisited, and maxVisited here do not go out of bound of the vector//
  for(size_t a = 0; a < maxMacrostate - minMacrostate + 1; a++)
  {
    if(histogram[a] != 0)
    {
      if(nonzeroCount==0) minVisited = a;
      maxVisited = a;
      nonzeroCount++;
    }
  }
  //printf("minVisited: %zu, maxVisited: %zu\n", minVisited, maxVisited);
  lnpi[minVisited] = 0.0;
  double maxlnpi = lnpi[minVisited];
  //Update the lnpi for the sampled region//
  //x: -1; y: 0; z: +1//
  for(size_t a = minVisited; a < maxVisited; a++)
  {
    //Zhao's note: add protection to avoid numerical issues//
    if(cmatrix[a].z   != 0) lnpi[a+1] = lnpi[a]   + std::log(cmatrix[a].z)   - std::log(cmatrix[a].x   + cmatrix[a].y   + cmatrix[a].z);   //Forward//
    forward_lnpi[a+1] = lnpi[a+1];
    if(cmatrix[a+1].x != 0) lnpi[a+1] = lnpi[a+1] - std::log(cmatrix[a+1].x) + std::log(cmatrix[a+1].x + cmatrix[a+1].y + cmatrix[a+1].z); //Reverse//
    reverse_lnpi[a+1] = lnpi[a+1];
    //printf("Loading: %zu, a+1, %zu, lnpi[a]: %.5f, lnpi[a+1]: %.5f\n", a, a+1, lnpi[a], lnpi[a+1]);
    if(lnpi[a+1] > maxlnpi) maxlnpi = lnpi[a+1];
  }
  //For the unsampled states, fill them with the minVisited/maxVisited stats//
  for(size_t a = 0; a < minVisited; a++) lnpi[a] = lnpi[minVisited];
  for(size_t a = maxVisited; a < maxMacrostate - minMacrostate + 1; a++) lnpi[a] = lnpi[maxVisited];
  //Normalize//
  double normalFactor = 0.0;
  for(size_t a = 0; a < maxMacrostate - minMacrostate + 1; a++) lnpi[a] -= maxlnpi;
  for(size_t a = 0; a < maxMacrostate - minMacrostate + 1; a++) normalFactor += std::exp(lnpi[a]); //sum of exp(lnpi)//
  //printf("Normalize Factor (Before): %.5f\n", normalFactor);
  normalFactor = -std::log(normalFactor); //Take log of normalFactor//
  //printf("Normalize Factor (After):  %.5f\n", normalFactor);
  for(size_t a = 0; a < maxMacrostate - minMacrostate + 1; a++)
  {
    lnpi[a] += normalFactor; //Zhao's note: mind the sign//
    bias[a]= -lnpi[a];
  }
};

//Clear Collection matrix stats (used after initialization cycles)
void TransitionMatrix::clearCMatrix()
{
  if(!doTMMC || !rezeroAfterInitialization) return;
  numberOfSteps = 0;
  double3 temp = {0.0, 0.0, 0.0};
  std::fill(cmatrix.begin(),  cmatrix.end(),  temp);
  std::fill(histogram.begin(),  histogram.end(),  0.0);
  std::fill(lnpi.begin(), lnpi.end(), 0.0);
  std::fill(bias.begin(), bias.end(), 1.0);
};

void TransitionMatrix::writeReport()
{
  std::ofstream textTMMCFile{};
  std::filesystem::path cwd = std::filesystem::current_path();

  std::string dirname="TMMC/";
  //std::string dirname=std::print("TMMC/System_{}/", systemId);
  std::string fname  = dirname + "/" + "TMMC_Statistics.data";

  std::filesystem::path directoryName = cwd /dirname;
  std::filesystem::path fileName = cwd /fname;
  std::filesystem::create_directories(directoryName);
  textTMMCFile = std::ofstream(fileName, std::ios::out);

  if(doTMMC)
  {
    textTMMCFile << "Performed " << numberOfSteps << " Steps\n";
    textTMMCFile << "Collection Matrix Updated " << numberOfUpdates << " Times\n";
    textTMMCFile << "Min Macrostate : " << minMacrostate << '\n';
    textTMMCFile << "Max Macrostate : " << maxMacrostate << '\n';
    textTMMCFile << "N CM[-1] CM[0] CM[1] bias lnpi Forward_lnpi Reverse_lnpi histogram" << '\n';
    for(size_t j = minMacrostate; j < maxMacrostate + 1; j++)
    {
      size_t newj = j - minMacrostate;
      textTMMCFile << j << " " << cmatrix[newj].x << " " << cmatrix[newj].y << " " << cmatrix[newj].z << " " << bias[newj] << " " 
                   << lnpi[newj] << " " << forward_lnpi[newj] << " " << reverse_lnpi[newj] << " " << histogram[newj] << '\n';
    }
  }
};
