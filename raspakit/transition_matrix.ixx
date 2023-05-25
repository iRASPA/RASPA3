export module transition_matrix;

import double3;

import <cmath>;
import <cstddef>;
import <vector>;

export struct TransitionMatrix
{
  enum class MoveType : size_t
  {
    Translation = 0,
    Rotation = 1,
    Insertion = 2,
    Deletion = 3,
    Reinsertion = 4,
    CF_LambdaChange = 5,
    CF_Insertion = 6,
    CF_Deletion = 7
  };

  std::vector<double3> CMatrix; //x = deletion, y = other move that doesnot change macrostate, z = insertion//
  std::vector<double>  TMBias;
  std::vector<double>  lnpi;    //For TM//
  std::vector<double>  forward_lnpi; //Debugging//
  std::vector<double>  reverse_lnpi;
  std::vector<double>  Histogram;

  size_t TMMC_Steps    = { 0 }; //Since the code is cycle-based, use an additional counter for TMMC//
  size_t MinMacrostate = { 0 };
  size_t MaxMacrostate = { 100 };
  size_t UpdateTMEvery = { 1000000 };
  size_t TMMC_n_Update = { 0 };

  bool DoTMMC = { false };
  bool DoUseBias = { false }; //Whether or not to use TM Bias for changing macrostates//
  bool UseTMBias = { true };  //Whether to use TM for the Bias//
  bool RejectOutofBound = { true }; //Whether to reject the move out of the bound of macrostate//
  bool RezeroAfterInitialization = { false };

  void Initialize_TMMC_Vectors();
  void update(double Pacc, size_t oldN, TransitionMatrix::MoveType moveType);
  double biasFactor(size_t oldN, TransitionMatrix::MoveType moveType);
  void AdjustTMBias();
  void ClearCMatrix();
  void writeTMMCReport();
};

