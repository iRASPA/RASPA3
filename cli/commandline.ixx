module;

export module commandline;

import std;

import archive;
import threadpool;
import input_reader;
import forcefield;
import monte_carlo;
import monte_carlo_transition_matrix;
import molecular_dynamics;
//import mixture_prediction_simulation;
//import isotherm_fitting_simulation;
//import multi_site_isotherm;
import opencl;
import getopt;
#ifdef BUILD_LIBTORCH
import libtorch_test;
#endif

export namespace CommandLine
{
enum State : std::uint8_t
{
  None = 0,
  Help = 1,
  OpenCL = 2,
  Input = 3,
  SurfaceArea = 4,
  VoidFraction = 5,
  TessellationComputation = 6,
  PSD = 7,
  EnergyGrid = 8,
  PSD_BV = 9, // Pore Size Distribution using Ban, Vlugt method
  Last = 10
};

void run(int argc, char* argv[]);
}  // namespace CommandLine
