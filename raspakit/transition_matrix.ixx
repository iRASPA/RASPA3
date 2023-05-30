export module transition_matrix;

import double3;

import <cmath>;
import <cstddef>;
import <vector>;

export struct TransitionMatrix
{
  std::vector<double3> cmatrix; //x = deletion, y = other move that doesnot change macrostate, z = insertion
  std::vector<double>  bias;
  std::vector<double>  lnpi;         //For TM//
  std::vector<double>  forward_lnpi; //Debugging//
  std::vector<double>  reverse_lnpi;
  std::vector<size_t>  histogram;

  size_t numberOfSteps   = { 0 }; //Since the code is cycle-based, use an additional counter for TMMC
  size_t minMacrostate   = { 0 };
  size_t maxMacrostate   = { 100 };
  size_t updateTMEvery   = { 1000000 };
  size_t numberOfUpdates = { 0 };

  bool doTMMC = { false };
  bool useBias = { false }; //Whether or not to use TM Bias for changing macrostates
  bool useTMBias = { true };  //Whether to use TM for the Bias
  bool rejectOutofBound = { true }; //Whether to reject the move out of the bound of macrostate
  bool rezeroAfterInitialization = { false };

  void initialize();
  void updateMatrix(double3 Pacc, size_t oldN);
  void updateHistogram(size_t N);
  double biasFactor(size_t newN, size_t oldN);
  void adjustBias();
  void clearCMatrix();
  void writeStatistics();
};

