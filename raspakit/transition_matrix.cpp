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

void TransitionMatrix::Initialize_TMMC_Vectors()
{
  if(!DoTMMC) return;
  CMatrix.resize(MaxMacrostate - MinMacrostate + 1);
  TMBias.resize(MaxMacrostate - MinMacrostate + 1); std::fill(TMBias.begin(), TMBias.end(), 1.0);
  lnpi.resize(MaxMacrostate - MinMacrostate + 1);
  forward_lnpi.resize(MaxMacrostate - MinMacrostate + 1);
  reverse_lnpi.resize(MaxMacrostate - MinMacrostate + 1);
  Histogram.resize(MaxMacrostate - MinMacrostate + 1);
}

void TransitionMatrix::update(double Pacc, size_t N, TransitionMatrix::MoveType moveType)
{
  if(!DoTMMC) return;
  if(Pacc > 1.0) Pacc = 1.0;
  switch(moveType)
  {
    case MoveType::Translation: 
    case MoveType::Rotation: 
    case MoveType::Reinsertion: 
    case MoveType::CF_LambdaChange:
    {
      if(RejectOutofBound && ((N > MaxMacrostate) || (N < MinMacrostate))) return;
      N -= MinMacrostate;
      CMatrix[N].y += Pacc;
      Histogram[N] ++;
      break;
    }
    case MoveType::Insertion: 
    case MoveType::CF_Insertion:
    {
      size_t OldN = N;
      N -= MinMacrostate;
      size_t NewN   = N + 1;
      CMatrix[N].z += Pacc;     //Insertion is the third value//
      CMatrix[N].y += 1.0 - Pacc;
      if(RejectOutofBound && ((OldN + 1) > MaxMacrostate)) return;
      Histogram[NewN] ++;
      break;
    }
    case MoveType::Deletion: 
    case MoveType::CF_Deletion:
    {
      size_t OldN = N;
      N -= MinMacrostate;
      size_t NewN   = N - 1;
      CMatrix[N].x += Pacc;  //Deletion is the first value//
      CMatrix[N].y += 1-Pacc;
      if(RejectOutofBound && ((OldN - 1) < MinMacrostate)) return;
      Histogram[NewN] ++;
      break;
    }
  }
};

//The bias is added to the preFactor//
double TransitionMatrix::biasFactor(size_t oldN, TransitionMatrix::MoveType moveType)
{
  if(!DoTMMC || !DoUseBias || !UseTMBias) return 1.0;
  if(oldN < MinMacrostate || oldN > MaxMacrostate) return 1.0; //No bias for macrostate out of the bound
  switch(moveType)
  {
    case MoveType::Translation: 
    case MoveType::Rotation: 
    case MoveType::Reinsertion: 
    case MoveType::CF_LambdaChange:
    {
      //Do not need the bias for moves that does not change the macrostate//
      return 1.0;
    }
    case MoveType::Insertion: 
    case MoveType::CF_Insertion:
    {
      if(RejectOutofBound && ((oldN + 1) > MaxMacrostate)) return 1.0;
      oldN -= MinMacrostate;
      double TMMCBias = TMBias[oldN + 1] - TMBias[oldN];
      return std::exp(TMMCBias); //See if Minus sign works//
    }
    case MoveType::Deletion: 
    case MoveType::CF_Deletion:
    {
      if(RejectOutofBound && ((oldN - 1) < MinMacrostate)) return 1.0;
      oldN -= MinMacrostate;
      double TMMCBias = TMBias[oldN - 1] - TMBias[oldN];
      return std::exp(TMMCBias); //See if Minus sign works//
    }
  }
  return 1.0;
};

//From Vince Shen's pseudo code//
void TransitionMatrix::AdjustTMBias()
{
  if(!DoTMMC || !DoUseBias || !UseTMBias) return;
  if((TMMC_Steps % UpdateTMEvery != 0) || TMMC_Steps == 0) return;
  printf("Adjusting TMBias\n");
  TMMC_n_Update ++;
  //First step is to get the lowest and highest visited states in terms of loading//
  size_t MinVisited = 0; size_t MaxVisited = 0;
  size_t nonzeroCount=0;
  //Zhao's special note: length of the vectors for TMMC = MaxMacrostate - MinMacrostate + 1;
  //The a, MinVisited, and MaxVisited here do not go out of bound of the vector//
  for(size_t a = 0; a < MaxMacrostate - MinMacrostate + 1; a++)
  {
    if(Histogram[a] != 0)
    {
      if(nonzeroCount==0) MinVisited = a;
      MaxVisited = a;
      nonzeroCount++;
    }
  }
  //printf("MinVisited: %zu, MaxVisited: %zu\n", MinVisited, MaxVisited);
  lnpi[MinVisited] = 0.0;
  double Maxlnpi = lnpi[MinVisited];
  //Update the lnpi for the sampled region//
  //x: -1; y: 0; z: +1//
  for(size_t a = MinVisited; a < MaxVisited; a++)
  {
    //Zhao's note: add protection to avoid numerical issues//
    if(CMatrix[a].z   != 0) lnpi[a+1] = lnpi[a]   + std::log(CMatrix[a].z)   - std::log(CMatrix[a].x   + CMatrix[a].y   + CMatrix[a].z);   //Forward//
    forward_lnpi[a+1] = lnpi[a+1];
    if(CMatrix[a+1].x != 0) lnpi[a+1] = lnpi[a+1] - std::log(CMatrix[a+1].x) + std::log(CMatrix[a+1].x + CMatrix[a+1].y + CMatrix[a+1].z); //Reverse//
    reverse_lnpi[a+1] = lnpi[a+1];
    //printf("Loading: %zu, a+1, %zu, lnpi[a]: %.5f, lnpi[a+1]: %.5f\n", a, a+1, lnpi[a], lnpi[a+1]);
    if(lnpi[a+1] > Maxlnpi) Maxlnpi = lnpi[a+1];
  }
  //For the unsampled states, fill them with the MinVisited/MaxVisited stats//
  for(size_t a = 0; a < MinVisited; a++) lnpi[a] = lnpi[MinVisited];
  for(size_t a = MaxVisited; a < MaxMacrostate - MinMacrostate + 1; a++) lnpi[a] = lnpi[MaxVisited];
  //Normalize//
  double NormalFactor = 0.0;
  for(size_t a = 0; a < MaxMacrostate - MinMacrostate + 1; a++) lnpi[a] -= Maxlnpi;
  for(size_t a = 0; a < MaxMacrostate - MinMacrostate + 1; a++) NormalFactor += std::exp(lnpi[a]); //sum of exp(lnpi)//
  //printf("Normalize Factor (Before): %.5f\n", NormalFactor);
  NormalFactor = -std::log(NormalFactor); //Take log of NormalFactor//
  //printf("Normalize Factor (After):  %.5f\n", NormalFactor);
  for(size_t a = 0; a < MaxMacrostate - MinMacrostate + 1; a++)
  {
    lnpi[a] += NormalFactor; //Zhao's note: mind the sign//
    TMBias[a]= -lnpi[a];
  }
};

//Clear Collection matrix stats (used after initialization cycles)
void TransitionMatrix::ClearCMatrix()
{
  if(!DoTMMC || !RezeroAfterInitialization) return;
  TMMC_Steps = 0;
  double3 temp = {0.0, 0.0, 0.0};
  std::fill(CMatrix.begin(),  CMatrix.end(),  temp);
  std::fill(Histogram.begin(),  Histogram.end(),  0.0);
  std::fill(lnpi.begin(), lnpi.end(), 0.0);
  std::fill(TMBias.begin(), TMBias.end(), 1.0);
};

void TransitionMatrix::writeTMMCReport()
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

  if(DoTMMC)
  {
    textTMMCFile << "Performed " << TMMC_Steps << " Steps\n";
    textTMMCFile << "Collection Matrix Updated " << TMMC_n_Update << " Times\n";
    textTMMCFile << "Min Macrostate : " << MinMacrostate << '\n';
    textTMMCFile << "Max Macrostate : " << MaxMacrostate << '\n';
    textTMMCFile << "N CM[-1] CM[0] CM[1] TMBias lnpi Forward_lnpi Reverse_lnpi Histogram" << '\n';
    for(size_t j = MinMacrostate; j < MaxMacrostate + 1; j++)
    {
      size_t newj = j - MinMacrostate;
      textTMMCFile << j << " " << CMatrix[newj].x << " " << CMatrix[newj].y << " " << CMatrix[newj].z << " " << TMBias[newj] << " " << lnpi[newj] << " " << forward_lnpi[newj] << " " << reverse_lnpi[newj] << " " << Histogram[newj] << '\n';
    }
  }
};
